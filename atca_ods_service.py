# ATCA Operational Data Sharing (ODS) service runner
# (C) Jamie Stevens 2024
#
# This code queries ATCA monitoring systems and tries to predict,
# as far ahead as is reasonable, the likely pointing positions
# of the ATCA, and creates a JSON file with that information.
# The JSON file follows the NRAO ODS format, which can be queried
# by satellite operators.

from astropy.time import Time
import datetime
import ephem
from html.parser import HTMLParser
import json
import math
import requests
import statistics
import time
import cabb_scheduler as cabb


def bat2mjd(bat=None, dUT=0):
    fbat = float(bat) / 1e6
    fbat = fbat - dUT
    fbat = fbat / 60 / 60 / 24
    return fbat

def stringToDegrees(t="00:00:00", hours=False):
    # Convert a sexagesimal string into degrees.
    e = t.split(":")
    if len(e) != 3:
        # Not enough information.
        return 0
    d = float(e[0])
    m = float(e[1])
    s = float(e[2])
    if t.startswith("-"):
        w = -1
    else:
        w = 1
    d = d * w
    r = d + (m / 60) + (s / 3600)
    if (hours):
        r = r * 15
    return r * w

def unwrapTime(tnow=None, tstr=None):
    # We create a time that has to be before the current time, from
    # a MoniCA value.
    if (tnow is None or tstr is None or tstr == "" or tstr == "null"):
        return 0
    els = tstr.split(":")
    tscan = datetime.datetime(tnow.year, tnow.month, tnow.day,
                              int(els[0]), int(els[1]), int(els[2]))
    if tnow < tscan:
        # This will happen when (eg) the scan time is 23:58:00 and the time
        # now is 00:02:00 the next day. We have to use the previous day for
        # the scan time.
        tscan = tscan - datetime.timedelta(days=1)
    return tscan

def scanLengthToDelta(scan):
    scanLength = scan.getScanLength()
    scanTimes = scanLength.split(":")
    scanDelta = datetime.timedelta(hours=int(scanTimes[0]), minutes=int(scanTimes[1]),
                                   seconds=int(scanTimes[2]))
    return scanDelta

class ScheduleParser(HTMLParser):
    def __init__(self):
        HTMLParser.__init__(self)
        self.schedData = None
        self.foundBlock = False
        
    def handle_starttag(self, tag, attrs):
        if tag == "pre":
            self.foundBlock = True

    def handle_endtag(self, tag):
        if self.foundBlock == True and tag == "pre":
            self.foundBlock = False
            
    def handle_data(self, data):
        if self.foundBlock == True:
            self.schedData = data

# A class for dealing with ATCA information.
class AtcaOds():
    def __init__(self):
        self.monicaSitePoints = {
            'frequency1': "site.misc.obs.freq1", # IF1 central frequency
            'frequency2': "site.misc.obs.freq2", # IF2 central frequency
            'sourceName': "site.misc.obs.source", # current source name
            'rightAscension': "site.misc.obs.target1", # current source RA
            'declination': "site.misc.obs.target2", # current source Dec
            'scheduleFile': "site.misc.obs.schedFile", # the current running schedule
            'weatherStowed': "site.misc.WeatherStow", # whether the array is weather stowed
            'arrayName': "site.misc.array", # the name of the array
            'scanNumber': "site.misc.obs.scan", # the current scan number
            'loopNumber': "site.misc.obs.loop",
            'scheduleCommand': "site.misc.obs.command", # what command is running
            'scanStartTime': "site.misc.obs.scanStart", # the start time of the current scan
            'scanEndTime': "site.misc.obs.scanEnd", # the end time of the current scan
            'dUTC': "caclock.misc.clock.dUTC" # the current dUTC between TAI and UTC
        }
        self.monicaAntennaPoints = {
            'caobsState': "misc.obs.caobsAntState",
            'azError': "servo.AzError",
            'elError': "servo.ElError",
            'wrap': "servo.AzWrap",
            'trackTime': "servo.TrackTime",
            'cycleTime': "cycle.Period",
            'cycleReference': "cycle.RefEpoch"
        }
        self.antennaNames = [ "ca01", "ca02", "ca03", "ca04", "ca05", "ca06" ]
        self.monica = self.prepareMonica()
        self.loadedSchedule = None
        self.scheduleFile = ""

    def prepareMonica(self, *args):
        # Set up a MoniCA connection for the parameters that we'll need later.
        monica = cabb.monica_information.monicaServer({
            'webserverName': 'www-int.narrabri.atnf.csiro.au'
        })
        
        pointsToAdd = [ self.monicaSitePoints[a] for a in self.monicaSitePoints ]
        for a in self.monicaAntennaPoints:
            for i in range(0, len(self.antennaNames)):
                pointsToAdd.append(self.antennaNames[i] + "." + self.monicaAntennaPoints[a])
        monica.addPoints(pointsToAdd)
                
                
        # Do an initial update.
        monica.updatePoints()
        return monica

    def getMonicaValue(self, pointname="", antenna=""):
        if ((pointname not in self.monicaSitePoints) and
            (pointname not in self.monicaAntennaPoints)):
            return None
        if pointname in self.monicaSitePoints:
            return self.monica.getPointByName(self.monicaSitePoints[pointname]).getValue()
        if pointname in self.monicaAntennaPoints and antenna in self.antennaNames:
            p = self.monica.getPointByName(antenna + "." + self.monicaAntennaPoints[pointname])
            if p is not None:
                return p.getValue()
        return None
            
    def getSchedule(self, schedFile=None):
        # Get the schedule from the server in the same way that the scheduler does.
        if schedFile is None:
            return None
        response = requests.get(url="https://www.narrabri.atnf.csiro.au/cgi-bin/sched/cabbsched",
                                params={ "cmd": "load",
                                         "file": schedFile })
        schedParser = ScheduleParser()
        schedParser.feed(response.text)
        schedule = cabb.schedule()
        schedule.parse(schedParser.schedData)
        schedParser.close()
        response.close()
        self.loadedSchedule = schedule
        self.scheduleFile = schedFile

    def getCurrentScan(self):
        if (self.loadedSchedule is None):
            return None
        try:
            scanNum = int(self.getMonicaValue('scanNumber'))
        except ValueError:
            scanNum = 1
        print(" current scan is %d" % scanNum)
        try:
            scan = self.loadedSchedule.getScan(scanNum - 1)
        except IndexError:
            scan = None
        return scan

    def getScanLength(self):
        scan = self.getCurrentScan()
        if scan is None:
            return 0
        return scanLengthToDelta(scan)
    
    def evaluateScan(self, scan=None, evTime=None):
        if scan is None or evTime is None:
            return None
        # Convert the scan into a fixed body.
        scanObject = ephem.FixedBody()
        scanObject._epoch = '2000'
        scanObject._ra = scan.getRightAscension()
        scanObject._dec = scan.getDeclination()

    def getArrayInformation(self):
        # We return the JSON keys needed for the array configuration.
        ods_obj = { 'site_lat_deg': -30.3128846, 'site_lon_deg': 149.5501388, 'site_el_m': 236.87,
                    'site_id': "atca_" + self.getMonicaValue('arrayName') }
        
        return ods_obj

    def getCurrentSourceDetails(self):
        robj = {
            'src_id': self.getMonicaValue('sourceName'),
            'src_ra_j2000_deg': stringToDegrees(self.getMonicaValue('rightAscension'),
                                                hours=True),
            'src_dec_j2000_deg': stringToDegrees(self.getMonicaValue('declination'))
        }
        return robj

    def getActiveAntennas(self):
        # We work out which antennas are currently active, and return that list.
        activeAnts = []
        for ant in self.antennaNames:
            status = self.getMonicaValue('caobsState', ant)
            if ((status != "DISABLED") and (status != "DETACHED") and (status != "DRIVE_ERR")):
                activeAnts.append(ant)
        return activeAnts

    def getIntegrationTime(self, activeAnts=[]):
        # Work out the cycle time for the active antennas.
        cycTime = 0
        for ant in activeAnts:
            c = float(self.getMonicaValue('cycleTime', ant))
            if (c > cycTime):
                cycTime = c
        return { 'corr_integ_time_sec': cycTime / 1e6, 'src_is_pulsar_bool': False }

    def getFrequencyConfiguration(self):
        # Return parameters associated with the frequency configuration.
        # Determine the frequency range.
        try:
            f1 = float(self.getMonicaValue('frequency1'))
            f2 = float(self.getMonicaValue('frequency2'))
        except ValueError:
            freq_low = 0
            freq_high = 0
            beam_width = 0
        else:
            freq_low = f1 - 1024
            freq_high = f1 + 1024
            f2_low = f2 - 1024
            f2_high = f2 + 1024
            if f2_low < freq_low:
                freq_low = f2_low
            if f2_high > freq_high:
                freq_high = f2_high
            freq_low = freq_low * 1e6
            freq_high = freq_high * 1e6
        
            # We use the lowest frequency to determine how big the beam radius is.
            wlength = 299792458 / freq_low
            beam_width = wlength * 180 / (22 * 2 * math.pi)
        
        return { 'freq_lower_hz': freq_low, 'freq_upper_hz': freq_high,
                 'src_radius': beam_width }

    def checkStowed(self, activeAnts=[]):
        stowedAnts = []
        for ant in activeAnts:
            state = monica.getPointByName(ant + "." + monicaAntennaPoints['caobsState']).getValue()
            if ((state == "STOWED") or (state == "STOWING") or (state == "PARKED") or
                (state == "PARKING")):
                stowedAnts.append(ant)
        if len(stowedAnts) > 0:
            # Figure out a reason why the antennas are stowed.
            weatherStow = False
            return "stowed"

    def getCurrentStatus(self):

        # Update MoniCA.
        self.monica.updatePoints()

        schedFile = self.getMonicaValue('scheduleFile')
        if (self.scheduleFile != schedFile):
            # The schedule file has changed.
            self.getSchedule(schedFile)

        dutc = int(self.getMonicaValue('dUTC'))
        time1 = Time(bat2mjd(int(self.getMonicaValue('trackTime', 'ca03'), 0),
                             dUT=dutc), format='mjd')
        
        robj = self.getArrayInformation()

        # The current source name, position and frequency configuration are always
        # just stored as is.
        robj.update(self.getCurrentSourceDetails())
        robj.update(self.getFrequencyConfiguration())
    
        # We consider ourselves in use if one antenna is not stowed/parked and caobs
        # is doing something with it.
        activeAnts = self.getActiveAntennas()
        somethingActive = (len(activeAnts) > 0)

        # The integration time will default to 0 if nothing is active.
        robj.update(self.getIntegrationTime(activeAnts))
    
        if somethingActive:
        
            # What is the time now?
            tnow = datetime.datetime.utcnow()
            # We get the start and (likely) end times of the current scan
            tscanstart = unwrapTime(tnow, self.getMonicaValue('scanStartTime'))
            tscanend = unwrapTime(tnow, self.getMonicaValue('scanEndTime'))
            print("full scan start time is %s" % tscanstart)
            print("full scan end time is %s" % tscanend)
            scanlength = self.getScanLength()
            print("scan is %s long" % scanlength)
        else:
            # Put in a note that nothing is active.
            robj['notes'] = "No active antennas at this time"

        return robj

# A class for dealing with ASKAP information.
class AskapOds():
    def __init__(self):
        self.endPoint = "https://prod-api.vlbi.atnf.tools/get_status/askap"
        self.status = None

    def getArrayInformation(self):
        # We return the JSON keys needed for the array configuration.
        ods_obj = { 'site_lat_deg': -26.697, 'site_lon_deg': 116.631425, 'site_el_m': 361.0,
                    'site_id': "askap", 'corr_integ_time_sec': 2.0 }
        return ods_obj

    def getCurrentSourceDetails(self):
        if self.status is None:
            return {}

        # We have to work out which antennas are working.
        antState = self.status['state']
        goodAntennas = [ i for i,x in enumerate(antState) if x is not None ]

        # Now get the RA and Dec for those good antennas.
        goodRa = [ stringToDegrees(self.status['rightAscensionICRF'][i]) for i in goodAntennas ]
        goodDec = [ stringToDegrees(self.status['declinationICRF'][i]) for i in goodAntennas ]
        # Find the average RA and Dec
        avgRa = statistics.mean(goodRa)
        avgDec = statistics.mean(goodDec)

        # Work out the observation timing.
        start_dt = datetime.datetime.fromisoformat(self.status['schedblock']['startTime'])
        delta_dt = datetime.timedelta(seconds=self.status['schedblock']['duration'])
        end_dt = start_dt + delta_dt
        timeFormat = "%Y-%m-%dT%H:%M:%S.%f"
        
        robj = {
            'src_id': self.status['schedblock']['alias'],
            'src_ra_j2000_deg': avgRa, 'src_dec_j2000_deg': avgDec,
            'src_is_pulsar_bool': False, 'slew_sec': 0,
            'trk_rate_ra_deg_per_sec': 0, 'trk_rate_dec_deg_per_sec': 0,
            'src_start_utc': start_dt.strftime(timeFormat),
            'src_end_utc': end_dt.strftime(timeFormat),
            'notes': "current observation"
        }
        return robj

    def getFrequencyConfiguration(self):
        if self.status is None:
            return {}

        freqs = self.status['configuration']['frequencies']
        lowfreqs = [ f - 180 for f in freqs ]
        highfreqs = [ f + 180 for f in freqs ]
        freq_low = min(lowfreqs) * 1e6
        freq_high = max(highfreqs) * 1e6

        return { 'freq_lower_hz': freq_low, 'freq_upper_hz': freq_high,
                 'src_radius': 5 }
    
    def getCurrentStatus(self):
        # Get the JSON from our endpoint.
        try:
            response = requests.get(url=self.endPoint)
            self.status = json.loads(response.text)
        except requests.exceptions.ConnectionError:
            print("Connection not made, using old data.")
        #print(self.status)

        robj = self.getArrayInformation()
        robj.update(self.getCurrentSourceDetails())
        robj.update(self.getFrequencyConfiguration())

        return robj
        
if __name__ == "__main__":

    # Do our prep work.
    # Set up the MoniCA server interface.
    #monica = prepareMonica()
    atca = AtcaOds()
    askap = AskapOds()

    # Directory containing the ODS objects.
    outDir = "/var/www/vhosts/www.narrabri.atnf.csiro.au/htdocs/ods"
    
    # Begin the main loop.
    while True:

        # The output JSON.
        askapOdsOutput = []

                     
        #print("CA01 tracking at %s" % time1.iso)

        # The first thing we do is get what the telescope is currently doing,
        # and this will be the first entry in the ODS array.
        #odsOutput.append(atca.getCurrentStatus())
        askapOdsOutput.append(askap.getCurrentStatus())

        askapOutObj = json.dumps({ "ods_data": askapOdsOutput})
        print(askapOutObj)

        # Write this out to a file.
        askapFile = outDir + "/ods.json"
        with open(askapFile, "w") as fp:
            print(askapOutObj, file=fp)
        
        time.sleep(30)
