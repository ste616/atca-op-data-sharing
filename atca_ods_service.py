# ATCA Operational Data Sharing (ODS) service runner
# (C) Jamie Stevens 2024
#
# This code queries ATCA monitoring systems and tries to predict,
# as far ahead as is reasonable, the likely pointing positions
# of the ATCA, and creates a JSON file with that information.
# The JSON file follows the NRAO ODS format, which can be queried
# by satellite operators.

import ephem
from html.parser import HTMLParser
import json
import requests
import time
import cabb_scheduler as cabb

monicaSitePoints = {
    'frequency1': "site.misc.obs.freq1", # IF1 central frequency
    'frequency2': "site.misc.obs.freq2", # IF2 central frequency
    'scheduleFile': "site.misc.obs.schedFile", # the current running schedule
    'weatherStowed': "site.misc.WeatherStow", # whether the array is weather stowed
    'arrayName': "site.misc.array", # the name of the array
    'scanNumber': "site.misc.obs.scan", # the current scan number
    'loopNumber': "site.misc.obs.loop",
    'scheduleCommand': "site.misc.obs.command" # what command is running
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
    print("schedule has %d scans" % schedule.getNumberOfScans())
    return schedule

def evaluateScan(scan=None, evTime=None):
    if scan is None or evTime is None:
        return None
    # Convert the scan into a fixed body.
    scanObject = ephem.FixedBody()
    scanObject._epoch = '2000'
    scanObject._ra = scan.getRightAscension()
    scanObject._dec = scan.getDeclination()

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
        print("schedule file is %s" % schedFile)
              
        currentSchedule = getSchedule(schedFile)
        time1 = monica.getPointByName("ca01." +
                                      monicaAntennaPoints['trackTime']).getValue()
        print("CA01 tracking at " + time1)
        
        
        time.sleep(30)
