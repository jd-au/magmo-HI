# Prepares an extended list of days with the patterns to match the RPFITS files.
#
# This takes in the file magmo-days.csv and produces the file magmo-days-full.csv
# The new file also contains a set of three patterns for each days which will
# match the day's RPFITS files.

# Author James Dempsey
# Date 9 Jul 2016

import csv
import datetime
from dateutil.relativedelta import relativedelta

# metadata
daysStartPostUtcDay = [1, 2, 3, 4, 5, 6, 42]
daysStartAfterEightPmUtc = [22, 23, 24]
daysPrevOnly = [31]
daysCurrOnlyShort = [43]

# functions
def getDayData():
    """
    Read the magmo-days.csv file and produce a list of day rows.
    :return: The list of day rows
    """

    dayData = []
    with open('magmo-days.csv', 'rb') as magmodays:
        reader = csv.reader(magmodays)
        next(reader)
        for row in reader:
            dayData.append(row)
    return dayData

def getFileDateFragments(dayRow):
    """
    Calculate three identifying date fragments which will match all data
    gathered for the day and no other data. This handles conversion from the
    Narrabri local dates in the original csv to the UTC date times used in
    RPFITS file names
    :param dayRow: The row of data for a day (will be day number and local date)
    :return: The date fragments for finding files named using ISO format UTC date times.
    """
    dayNum = int(dayRow[0])
    dayDateTime = datetime.datetime.strptime(dayRow[1], "%d/%m/%Y")
    #print dayDateTime
    dayDate = dayDateTime.date()
    prevDayDate = dayDate - relativedelta(days=1)
    #print prevDayDate.isoformat()
    days = []
    if dayNum in daysStartPostUtcDay:
        for i in range(0,3):
            days.append(dayDate.isoformat())
    elif dayNum in daysStartAfterEightPmUtc:
        days.append(prevDayDate.isoformat() + "_2")
        days.append(dayDate.isoformat() + "_0")
        days.append(dayDate.isoformat() + "_1")
    elif dayNum in daysPrevOnly:
        days.append(prevDayDate.isoformat() + "_1")
        days.append(prevDayDate.isoformat() + "_2")
        days.append(prevDayDate.isoformat() + "_2")
    elif dayNum in daysCurrOnlyShort:
        days.append(dayDate.isoformat() + "_0")
        days.append(dayDate.isoformat() + "_1")
        days.append(dayDate.isoformat() + "_1")
    else:
        days.append(prevDayDate.isoformat()+"_1")
        days.append(prevDayDate.isoformat()+"_2")
        days.append(dayDate.isoformat()+"_0")
    return days


# ### Script starts here ###

dayData = getDayData()
for row in dayData:
    dateFrag = getFileDateFragments(row)
    for frag in dateFrag:
        row.append(frag)

with open('magmo-days-full.csv', 'wb') as magmodays:
    magmodays.write('day,date,pattern1,pattern2,pattern3\n')
    for row in dayData:
        entry = ','.join(row) + '\n'
        magmodays.write(entry)

print "Processed %s days." % len(dayData)