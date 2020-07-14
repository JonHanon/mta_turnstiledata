import csv
from datetime import datetime, timedelta
from pytz import timezone
import numpy
import matplotlib
from matplotlib import pyplot, dates
from matplotlib.dates import date2num, num2date
import dateutil.relativedelta as REL
from scipy.interpolate import PchipInterpolator as pchip

def overlap(List1 , List2):
    overlap = []
    for key in List1: 
        if key in List2: 
            overlap += [key]
    return overlap

class recordCA:
    def __init__(self, CA, name, line):
        self.turnstiles = {}
        self.entries = {}
        self.exits = {}
        self.name = name
        self.line = line
        self.padded_str = (CA+': ').ljust(7)+name.rjust(20)+','+line

class recordSCP:
    def __init__(self):
        self.db_entries = {}
        self.db_exits = {}
    def insert(self,timestamp,entries,exits):
        self.db_entries[timestamp] = entries
        self.db_exits[timestamp] = exits

class recordStation:
    def __init__(self, name, line):
        self.CAs = []
        self.entries = {}
        self.exits = {}
        self.name = name
        self.line = line
        self.title = name + ' (' + line + ')'

def read_files(files):
    record = {}
    for file in files:
        read_file(record, file)
    return record

def read_file(record, file):
    with open(file, 'r') as csvfile:
        header = [h.strip() for h in csvfile.readline().split(',')]
        reader = csv.DictReader(csvfile, fieldnames=header)
        ny_tz = timezone('US/Eastern')
        prevSCP = '' # initialized for first iteration
        prevCA = ''
        for row in reader:
            if row['DIVISION'] == 'PTH':
                continue
            CA = row['C/A']
            if CA[0:3] == 'JFK':
                continue
            SCP = row['SCP']
            if prevSCP != SCP: # new turnstile
                if prevCA != CA: # new station identifier
                    line = row['LINENAME']
                    if CA not in record:
                        record[CA] = recordCA(CA, row['STATION'], ''.join(sorted(line)))
##                    print(record[CA].padded_str)
                if SCP not in record[CA].turnstiles:
                    record[CA].turnstiles[SCP] = recordSCP()
            datetime_str = str(row['DATE'])+' '+str(row['TIME'])
            ts = ny_tz.localize(datetime.strptime(datetime_str, "%m/%d/%Y %H:%M:%S"))
            record[CA].turnstiles[SCP].insert(ts,int(row['ENTRIES']), int(row['EXITS']))
            prevSCP = SCP
            prevCA = CA

def generate_CA_splines(record):
    for CA, rCA in record.items():
        SCPentrySpline = {}
        SCPexitSpline = {}
        timestamps = {}
        for SCP, rSCP in rCA.turnstiles.items():
            time_from_last = {}
            cum_entries = {}
            cum_exits = {}
            firstIter = True # too complicated not to use a boolean
            for ts in rSCP.db_entries:
                pop_prev = False
                entries = rSCP.db_entries[ts]
                exits = rSCP.db_exits[ts]
                if firstIter: # first iteration has no previous
                    firstIter = False
                    ts_prev = ts
                    timestamps[ts] = True
                    time_from_last = timedelta(days=25000)
                    cum_entries[ts] = 0
                    cum_exits[ts] = 0
                    entries_prev = entries
                    exits_prev = exits
                    if len(rSCP.db_entries) == 1: # can't spline one point
                        ts = ts+timedelta(seconds=1)
                        timestamps[ts] = True
                        rSCP.db_entries[ts] = rSCP.db_entries[ts_prev]
                        rSCP.db_exits[ts] = rSCP.db_exits[ts_prev]
                        cum_entries[ts] = 0
                        cum_exits[ts] = 0
                        time_from_last = timedelta(seconds=1)
                        break
                    continue
                timediff = ts-ts_prev
                if timediff < timedelta(hours=1): # dirty data
                    if time_from_last_prev < timedelta(hours=4):
                        timestamps.pop(ts_prev) # add on to previous
                        timediff += time_from_last_prev
                        pop_prev = True
                    else: # add on to next
                        continue
                else: # reset count to use same code for all cases
                    diff_entries = 0
                    diff_exits = 0
                timestamps[ts] = True
                time_from_last = timediff
                diff_entries += min(abs(entries - entries_prev), entries)
                cum_entries[ts] = cum_entries[ts_prev]+diff_entries
                diff_exits += min(abs(exits - exits_prev), exits)
                cum_exits[ts] = cum_exits[ts_prev]+diff_exits
                if diff_entries > 50000:
                    print(rCA.name + ' ' + CA + ' ' + SCP + ' ' + ts.strftime('%m/%d/%Y %H:%M:%S') + ' entries ' + str(diff_entries))
                if diff_exits > 50000:
                    print(rCA.name + ' ' + CA + ' ' + SCP + ' ' + ts.strftime('%m/%d/%Y %H:%M:%S') + ' exits ' + str(diff_exits))
                if pop_prev:
                    cum_entries.pop(ts_prev)
                    cum_exits.pop(ts_prev)
                ts_prev = ts
                entries_prev = entries
                exits_prev = exits
                time_from_last_prev = time_from_last
##            if rCA.name == 'TIMES SQ-42 ST' or rCA.name == '42ST-PORT AUTH':
##                print(SCP + ' ' + CA)
##                for entry in 
##                print(cum_exits)
##                print()
            times = date2num([*cum_entries.keys()])
            SCPentrySpline[SCP] = pchip(times, [*cum_entries.values()])
            SCPexitSpline[SCP] = pchip(times, [*cum_exits.values()])
            del rSCP.db_entries, rSCP.db_exits # clear RAM as we go
        for ts in sorted(timestamps.keys()):
            rCA.entries[ts] = 0
            rCA.exits[ts] = 0
            ts_num = date2num(ts)
            for SCP, rSCP in rCA.turnstiles.items():
                entAdd = SCPentrySpline[SCP](ts_num)
                exAdd = SCPexitSpline[SCP](ts_num)
                if entAdd != entAdd: # splines are NaN outside bounds
                    continue
                rCA.entries[ts] += entAdd
                rCA.exits[ts] += exAdd
        del rCA.turnstiles # clear RAM as we go
        times = date2num([*rCA.entries.keys()])
        rCA.entrySpline = pchip(times,[*rCA.entries.values()])
        rCA.exitSpline = pchip(times,[*rCA.exits.values()])
        del rCA.exits # we still need rCA.entries' keys one more time

def findStations(record, junctionDict = {}, junctionNames = {}):
    stations = {}
    for CA, rCA in record.items():
        (name, line) = (rCA.name,rCA.line)
        if (name, line) in junctionDict: # Merge junctions
            (name, line) = junctionDict[(name, line)]
        if (name, line) not in stations:
            stations[(name, line)] = recordStation(name, line)
        stations[(name, line)].CAs.append(CA)
    for station in junctionNames: # Set a merged title for junctions
        if station in stations:
            stations[station].title = junctionNames[station]
    for (name, line), rStation in stations.items():
        if (name, line) in junctionDict:
            (name, line) = junctionDict[(name, line)]
        timestamps = {}
        for CA in rStation.CAs:
            for ts in record[CA].entries:
                timestamps[ts] = True
            del record[CA].entries
        for ts in sorted(timestamps.keys()):
            rStation.entries[ts] = 0
            rStation.exits[ts] = 0
    ny_tz = timezone('US/Eastern')
    for station, rStation in stations.items():
##        print(rStation.title + ' ' + str(rStation.CAs))
##        print('Station: ' + rStation.title)
        for ts in rStation.entries:
##            print('Timestamp: ' + ts.strftime('%m/%d/%Y %H:%M:%S'))
            ts_num = date2num(ts)
            for CA in rStation.CAs:
                rCA = record[CA]
                entAdd = rCA.entrySpline(ts_num)
                exAdd = rCA.exitSpline(ts_num)
                if entAdd != entAdd: # dirty data outside of the bounds
                    continue
                rStation.entries[ts] += entAdd
                rStation.exits[ts] += exAdd
        times = date2num([*rStation.entries.keys()])
        rStation.entrySpline = pchip(times,[*rStation.entries.values()])
        rStation.exitSpline = pchip(times,[*rStation.exits.values()])
    return stations
            
def countStations(stations, boolPrint=False):
    numLines = 0
    for station in stations:
        if boolPrint:
            print(station.name + ' ' + station.line)
        numLines += 1
    return numLines

########## Hardcoded Information
junctionArray = [['FRANKLIN AV / BOTANIC GARDEN (2345S)',                  ('FRANKLIN AV', '2345S'), ('BOTANIC GARDEN', '2345S')],\
                 ['JAY ST-METROTEC (ACFR)',                                ('JAY ST-METROTEC', 'R'), ('JAY ST-METROTEC', 'ACF')],\
                 ['LORIMER ST / METROPOLITAN AV (GL)',                     ('LORIMER ST', 'GL'), ('METROPOLITAN AV', 'GL')],\
                 ['14 ST / 8 AV (ACEL)',                                   ('14 ST', 'ACEL'), ('8 AV', 'ACEL')],\
                 ['14 ST / 6 AV (123FLM)',                                 ('14 ST', '123FLM'), ('6 AV', '123FLM')],\
                 ['34 ST-PENN STA (123ACE)',                               ('34 ST-PENN STA', '123ACE'), ('34 ST-PENN STA', 'ACE'),\
                                                                           ('34 ST-PENN STA', '123')],\
                 ['42 ST-BRYANT PK / 5 AVE (7BDFM)',                       ('42 ST-BRYANT PK', '7BDFM'), ('5 AVE', '7BDFM')],\
                 ['B\'WAY-LAFAYETTE / BLEECKER ST (6BDFQ)',                ('B\'WAY-LAFAYETTE', '6BDFQ'), ('BLEECKER ST', '6DF')],\
                 ['BROOKLYN BRIDGE / CHAMBERS ST (456JZ)',                 ('BROOKLYN BRIDGE', '456JZ'), ('CHAMBERS ST', '456JZ')],\
                 ['CHAMBERS ST / WTC / PARK PLACE / CORTLANDT (23ACENRW)', ('CHAMBERS ST', '23ACE'), ('WORLD TRADE CTR', '23ACE'),\
                                                                           ('PARK PLACE', '23ACE'), ('CORTLANDT ST', 'NRW')],\
                 ['LEXINGTON AV/53 / 51 ST (6EM)',                         ('LEXINGTON AV/53', '6EM'), ('51 ST', '6')],\
                 ['SOUTH FERRY / WHITEHALL (1RW)',                         ('SOUTH FERRY', '1RW'), ('WHITEHALL S-FRY', '1RW')],\
                 ['TIMES SQ-42 ST-PORT AUTH (1237ACEGNQRSW)',              ('TIMES SQ-42 ST', '1237ACENQRSW'), ('TIMES SQ-42 ST', '1237ACENQRS'),\
                                                                           ('42 ST-PORT AUTH', '1237ACEGNRSW'), ('42 ST-PORT AUTH', '1237ACENQRSW')],\
                 ['74-BROADWAY / JKSN HT-ROOSVLT (7EFMR)',                 ('74 ST-BROADWAY', '7EFMR'), ('JKSN HT-ROOSVLT', '7EFMR')],\
                 ['COURT SQ-23 ST (7EGM)',                                 ('COURT SQ', 'EGM'), ('COURT SQ', '7'), ('COURT SQ-23 ST', 'EGM')]]

########## Preprocessing
junctionDict = {}
junctionNames = {}
for j in junctionArray:
##    print(j[0])
##    print(j[1])
    junctionNames[j[1]] = j[0]
    for s in j[2:]:
        junctionDict[s] = j[1]
##        print(s)
##    print()

########## Main Code
matplotlib.rcParams['timezone'] = 'US/Eastern'
start_date = datetime(2019,1,26)
rd = REL.relativedelta(days=1, weekday=REL.SA) # get next Saturday
start_date += rd
end_date = datetime(2019,2,22)
end_date += rd
    
datestr_list = []
while start_date <= end_date:
    datestr_list += [start_date.strftime('%y%m%d')]
    start_date += timedelta(7)

files = []
for s in datestr_list:
    files += ['turnstile_'+s+'.txt']
print(files)
record = read_files(files)
generate_CA_splines(record)
stations = findStations(record, junctionDict, junctionNames)
#stations = findStations(record)

##for st, name in junctionNames.items():
##    if st in stations:
##        print(stations[st].title + ' found: ' + stations[st].name)
##    else:
##        print(str(st) + ' not found')

##numLines = countStations(stations)
##print('there are ' + str(numLines) + ' stations\n')

plot_stations = [('14 ST-UNION SQ','456LNQRW'), ('34 ST-PENN STA', '123ACE'), ('TIMES SQ-42 ST', '1237ACENQRSW')]
for station in plot_stations:
    data = stations[station]
    ts_obj = [*data.entries.keys()]
    ts_val = date2num(ts_obj)
    tplot = numpy.linspace(ts_val[0],ts_val[-1], 10000)

    fig = pyplot.figure()
    ax = pyplot.gca()
    pyplot.plot_date(tplot, data.entrySpline.derivative()(tplot)/24, fmt='b-') # derivative is in date2num, which is per day
    #pyplot.plot_date(ts_val, [*data.entries.values()])
    fig.suptitle(data.title + ' Entries ' + ts_obj[0].strftime('%m/%d/%Y') + ' to ' + ts_obj[-1].strftime('%m/%d/%Y'))
    pyplot.ylabel('Entries per hour')
    xfmt = dates.DateFormatter('%Y-%m-%d %H:%M:%S')
    ax.xaxis.set_major_formatter(xfmt)
    ##pyplot.xticks(ticks)
    fig.autofmt_xdate()
    fig.set_size_inches(16, 8)
    pyplot.show()

    fig = pyplot.figure()
    ax = pyplot.gca()
    pyplot.plot_date(tplot, data.exitSpline.derivative()(tplot)/24, fmt='b-') # derivative is in date2num, which is per day
    fig.suptitle(data.title + ' Exits ' + ts_obj[0].strftime('%m/%d/%Y') + ' to ' + ts_obj[-1].strftime('%m/%d/%Y'))
    pyplot.ylabel('Exits per hour')
    xfmt = dates.DateFormatter('%Y-%m-%d %H:%M:%S')
    ax.xaxis.set_major_formatter(xfmt)
    ##pyplot.xticks(ticks)
    fig.autofmt_xdate()
    fig.set_size_inches(16, 8)
    pyplot.show()
