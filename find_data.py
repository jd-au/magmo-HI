# Script to query ATOA to find the RPFITS files representing a day's observations.
# Each file is then downloaded.
# Note in the future, an interfsce will be provided in ATOA to allow a list of file
# urls which do not require login to be downloaded and then used in a tool such as
# wget,

# Author: James Dempsey
# Date: 6 August, 2016

import magmo
import sys
import time
# REST requests
import urllib, urllib2, base64
# VO Table parsing
from astropy.io import votable
# For hidden entry of password
import getpass
import re
import requests

# Constants
atoa_tap_url = 'http://atoavo.atnf.csiro.au/tap/sync'
atoa_login_url = 'http://atoa.atnf.csiro.au/login'
atoa_download_service = 'http://atoa.atnf.csiro.au/listDownload.jsp'
obs_prog_id = 'C2291'

chunk_size = 4*1024 # bytes


# http://horus.atnf.csiro.au:8081/atoavo/tap/sync?request=doQuery&lang=ADQL&format=votable&query=select+*+from+ivoa.obscore


def adql_query(url, query_string, filename, username=None, password=None, file_write_mode='w'):
    """ Do an adql query, and write the resulting VO Table to a file """
    req = urllib2.Request(url)
    # Uses basic auth to securely access the data access information for the image cube
    if username is not None:
        base64string = base64.encodestring('%s:%s' % (username, password)).replace('\n', '')
        req.add_header("Authorization", "Basic %s" % base64string)
    data = urllib.urlencode({'query':query_string, 'request':'doQuery', 'lang':'ADQL', 'format':'votable'})
    u = urllib2.urlopen(req, data)
    queryResult = u.read()
    queryResult = re.sub('character varying\([0-9]+\)', 'char" arraysize="*', queryResult)
    with open(filename,file_write_mode) as f:
        f.write(queryResult)


def query_atoa(day_row):
    obs_ids = []
    base_query = "SELECT distinct access_url " \
                 + "FROM ivoa.obscore where obs_collection = 'C2291' " +\
                   "and frequency in (1421.0, 1420.5) and data_flag < 999 "

    temp_dir = 'temp'
    magmo.ensure_dir_exists(temp_dir)
    temp_file = temp_dir + "/query-result.xml"
    day_select = ' and ('
    for i in range(2,len(day_row)):
        if i > 2:
            day_select += ' or '
        day_select += "obs_id like '" + day_row[i] + "%." + obs_prog_id + "'"
    day_select += ')'
    query = base_query + day_select + " order by 1"

    adql_query(atoa_tap_url, query, temp_file)
    result_votable = votable.parse(temp_file, pedantic=False)
    results = result_votable.get_first_table().array
    for row in results:
        obs_id = row['access_url']
        # print obs_id
        if obs_id is not None:
            obs_ids.append(obs_id)

    return obs_ids


def login_to_atoa(userid, password):
    session = requests.session()

    # This is the form data that the page sends when logging in
    login_data = {
        'j_username': userid,
        'j_password': password,
        'submit': 'Login',
        '_action': 'login',
    }

    # Authenticate
    r = session.post(atoa_login_url, data=login_data)
    print r.headers
    return session


def get_download_urls(obs_ids, opener):
    data = ''
    for id in obs_ids:
        data += id + '\n'
    form_data = urllib.urlencode({'filelist': data, 'bundle':'textlist'})
    resp = opener.open(atoa_download_service, form_data)
    urls = resp.read()
    return urls


def download_files(urls, session):
    for url in urls:
        filename = url[url.find('fname')+6:]
        print 'Downloading file ', filename
        continue
        r = session.get(url, stream=True)

        with open(filename, 'wb') as fd:
            for chunk in r.iter_content(chunk_size):
                fd.write(chunk)
                print "."


def main():
    # Read day parameter
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("Incorrect number of parameters.")
        print("Usage: python find_data.py day userid [password]")
        exit(1)
    day = sys.argv[1]
    userid = sys.argv[2]
    password = None
    if len(sys.argv) > 3:
        password = sys.argv[3]
    else:
        password = getpass.getpass("Enter your OPAL password: ")

    start = time.time()

    day_row = magmo.get_day_file_data(day)
    if day_row is None:
        print "Day %s is not defined." % (day)
        exit(1)
    print "#### Started finding data for MAGMO day %s at %s ####" % \
          (day, time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start)))
    print day_row

    # Query ATOA for desired ids
    obs_ids = query_atoa(day_row)

    session = login_to_atoa(userid, password)
    download_files(obs_ids, session)
    #urls = get_download_urls(obs_ids, opener)
    #print urls


    # Report
    end = time.time()
    print '#### Finding completed at %s ####' \
          % (time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)))
    print 'Processed in %.02f s' % (end - start)
    exit(0)


# Run the script if it is called from the command line
if __name__ == "__main__":
    main()