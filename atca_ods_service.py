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
import time
import cabb_scheduler as cabb

monicaSitePoints = {
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

monicaAntennaPoints = {
    'caobsState': "misc.obs.caobsAntState",
    'azError': "servo.AzError",
    'elError': "servo.ElError",
    'wrap': "servo.AzWrap",
    'trackTime': "servo.TrackTime",
    'cycleTime': "cycle.Period",
    'cycleReference': "cycle.RefEpoch"
}
antennaNames = [ "ca01", "ca02", "ca03", "ca04", "ca05", "ca06" ]

def bat2mjd(bat=None, dUT=0):
    fbat = float(bat) / 1e6
    fbat = fbat - dUT
    fbat = fbat / 60 / 60 / 24
    return fbat

def prepareMonica(*args):
    # Set up a MoniCA connection for the parameters that we'll need later.
    monica = cabb.monica_information.monicaServer({
        'webserverName': 'www-int.narrabri.atnf.csiro.au'
    })

    global monicaSitePoints
    global monicaAntennaPoints
    global antennaNames
    pointsToAdd = [ monicaSitePoints[a] for a in monicaSitePoints ]
    for a in monicaAntennaPoints:
        for i in range(0, len(antennaNames)):
            pointsToAdd.append(antennaNames[i] + "." + monicaAntennaPoints[a])
    monica.addPoints(pointsToAdd)
    
    
    # Do an initial update.
    monica.updatePoints()
    return monica

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


def getSchedule(schedFile=None):
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
    #print("schedule has %d scans" % schedule.getNumberOfScans())
    return schedule

def evaluateScan(scan=None, evTime=None):
    if scan is None or evTime is None:
        return None
    # Convert the scan into a fixed body.
    scanObject = ephem.FixedBody()
    scanObject._epoch = '2000'
    scanObject._ra = scan.getRightAscension()
    scanObject._dec = scan.getDeclination()

def getArrayInformation():
    # We return the JSON keys needed for the array configuration.
    ods_obj = { 'site_lat_deg': -30.3128846, 'site_lon_deg': 149.5501388, 'site_el_m': 236.87,
                'site_id': "atca_" + monica.getPointByName(monicaSitePoints['arrayName']).getValue() }
    
    return ods_obj

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

def getCurrentSourceDetails():
    robj = {
        'src_id': monica.getPointByName(monicaSitePoints['sourceName']).getValue(),
        'src_ra_j2000_deg': stringToDegrees(monica.getPointByName(monicaSitePoints['rightAscension']).getValue(),
                                            hours=True),
        'src_dec_j2000_deg': stringToDegrees(monica.getPointByName(monicaSitePoints['declination']).getValue())
    }
    return robj
    

def getActiveAntennas():
    # We work out which antennas are currently active, and return that list.
    activeAnts = []
    for ant in antennaNames:
        status = monica.getPointByName(ant + "." + monicaAntennaPoints['caobsState']).getValue()
        if ((status != "DISABLED") and (status != "DETACHED") and (status != "DRIVE_ERR")):
            activeAnts.append(ant)
    return activeAnts

def getIntegrationTime(activeAnts=[]):
    # Work out the cycle time for the active antennas.
    cycTime = 0
    for ant in activeAnts:
        c = float(monica.getPointByName(ant + "." + monicaAntennaPoints['cycleTime']).getValue())
        if (c > cycTime):
            cycTime = c
    return { 'corr_integ_time_sec': cycTime / 1e6, 'src_is_pulsar_bool': False }

def getFrequencyConfiguration():
    # Return parameters associated with the frequency configuration.
    # Determine the frequency range.
    f1 = float(monica.getPointByName(monicaSitePoints['frequency1']).getValue())
    f2 = float(monica.getPointByName(monicaSitePoints['frequency2']).getValue())
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
        
    return { 'freq_lower_hz': freq_low * 1e6, 'freq_upper_hz': freq_high * 1e6,
             'src_radius': beam_width }

def checkStowed(activeAnts=[]):
    stowedAnts = []
    for ant in activeAnts:
        state = monica.getPointByName(ant + "." + monicaAntennaPoints['caobsState']).getValue()
        if ((state == "STOWED") or (state == "STOWING") or (state == "PARKED") or
            (state == "PARKING")):
            stowedAnts.append(ant)
    if len(stowedAnts) > 0:
        # Figure out a reason why the antennas are stowed.
        weatherStow = False

def getCurrentStatus():

    robj = getArrayInformation()

    # The current source name, position and frequency configuration are always
    # just stored as is.
    robj.update(getCurrentSourceDetails())
    robj.update(getFrequencyConfiguration())
    
    # We consider ourselves in use if one antenna is not stowed/parked and caobs
    # is doing something with it.
    activeAnts = getActiveAntennas()
    somethingActive = (len(activeAnts) > 0)

    # The integration time will default to 0 if nothing is active.
    robj.update(getIntegrationTime(activeAnts))
    
    if somethingActive:
        
        # What is the time now?
        tnow = datetime.datetime.utcnow()
        print(tnow)
        # We get the start and (likely) end times of the current scan
        scanStartTime = monica.getPointByName(monicaSitePoints['scanStartTime']).getValue()
        scanEndTime = monica.getPointByName(monicaSitePoints['scanEndTime']).getValue()
        startels = scanStartTime.split(":")
        endels = scanEndTime.split(":")
        # Create a datetime from the base of the current time and the time that
        # the scan started.
        tscanstart = datetime.datetime(tnow.year, tnow.month, tnow.day,
                                       int(startels[0]), int(startels[1]), int(startels[2]))
        if tnow < tscanstart:
            # This will happen when (eg) the scan start time is 23:58:00 and the time
            # now is 00:02:00 the next day. We have to use the previous day for the
            # scan start time.
            tscanstart = tscanstart - datetime.timedelta(days=1)
        print("full scan start time is %s" % tscanstart)
        # Create a datetime from the base of the current time and the time that
        # the scan will (likely) end.
        tscanend = datetime.datetime(tnow.year, tnow.month, tnow.day,
                                     int(endels[0]), int(endels[1]), int(endels[2]))
        if tnow < tscanend:
            # This will happen when (eg) the scan start time is 23:58:00 and the time
            # now is 00:02:00 the next day. We have to use the previous day for the
            # scan end time.
            tscanend = tscanend - datetime.timedelta(days=1)
        print("full scan end time is %s" % tscanend)
    else:
        # Put in a note that nothing is active.
        robj['notes'] = "No active antennas at this time"

    return robj
    
if __name__ == "__main__":

    # Do our prep work.
    # Set up the MoniCA server interface.
    monica = prepareMonica()
    
    # Begin the main loop.
    while True:

        # The output JSON.
        odsOutput = []

        monica.updatePoints()
        schedFile = monica.getPointByName(monicaSitePoints['scheduleFile']).getValue()
        #print("schedule file is %s" % schedFile)
              
        currentSchedule = getSchedule(schedFile)
        dutc = int(monica.getPointByName(monicaSitePoints['dUTC']).getValue())
        time1 = Time(bat2mjd(int(monica.getPointByName("ca03." +
                                                       monicaAntennaPoints['trackTime']).getValue(), 0),
                             dUT=dutc), format='mjd')
                     
        #print("CA01 tracking at %s" % time1.iso)

        # The first thing we do is get what the telescope is currently doing,
        # and this will be the first entry in the ODS array.
        odsOutput.append(getCurrentStatus())

        outObj = json.dumps(odsOutput)
        print(outObj)
        
        time.sleep(30)
