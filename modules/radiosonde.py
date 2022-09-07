"""
Filename:    radiosonde.py
Author:      Deanna Nash, dlnash@ucsb.edu
Description: Functions used to download and process University of Wyoming Radiosonde Data
"""
import pandas as pd
from bs4 import BeautifulSoup
import urllib3

def download_university_wyoming_radiosonde(dates, stn, download=True):
    # Download, process, and Plot Sounding Data from University of Wyoming
    # https://kbkb-wx-python.blogspot.com/2015/07/plotting-sounding-data-from-university.html
    path_to_data = '../data/'
    filenames = []
    for i, date in enumerate(dates):
        year = date.strftime('%Y')
        month = date.strftime('%m')
        day = date.strftime('%d')
        hour = date.strftime('%H')
        if download == True: 
            # 1)
            # Wyoming URL to download Sounding from
            url = 'http://weather.uwyo.edu/cgi-bin/sounding?region=naconf&TYPE=TEXT%3ALIST&YEAR='+year+'&MONTH='+month+'&FROM='+day+hour+'&TO='+day+hour+'&STNM='+stn
            http = urllib3.PoolManager()
            response = http.request('GET', url)

            # 2)
            # Remove the html tags
            soup = BeautifulSoup(response.data.decode('utf-8'))
            data_text = soup.get_text()

            # 3)
            # Split the content by new line.
            splitted = data_text.split("\n",data_text.count("\n"))

            # 4)
            # Write this splitted text to a .txt document
            Sounding_filename = path_to_data + 'radiosonde/' + str(stn)+'.'+str(year)+str(month)+str(day)+str(hour)+'.txt'
            filenames.append(Sounding_filename)
            f = open(Sounding_filename,'w')
            for line in splitted[4:]:
                f.write(line+'\n')
            f.close()
        else:
            # get sounding filenames
            Sounding_filename = path_to_data + 'radiosonde/' + str(stn)+'.'+str(year)+str(month)+str(day)+str(hour)+'.txt'
            filenames.append(Sounding_filename)
            
    return filenames

def lines_that_contain(string, fp):
    return [line for line in fp if string in line]

def get_lat_lon_wyoming(f):
    '''Find the lines that contain the station latitude and longitude and return the values as float'''
    with open(f, "r") as fp:
        for line in lines_that_contain("Station latitude:", fp):
            line_lat = line
    with open(f, "r") as fp:        
        for line in lines_that_contain("Station longitude:", fp):
            line_lon = line

    slat = line_lat[-6:-1]
    slon = line_lon[-6:-1]
    
    return float(slat), float(slon)

def get_skipfooter_wyoming(f):
    '''Determine the number of lines that need to be skipped at the end when reading the sounding data'''
    ## find where the table data ends
    lookup = 'Station information and sounding indices'
    with open(f) as fp:
        for num, line in enumerate(fp, 1):
            if lookup in line:
                y = num
    # count the total number of lines
    with open(f) as fp:
        x = len(fp.readlines())
    # this is the number of lines at the end of the file that need to be skipped
    skipfoot = (x-y)+1 # add one bc python indexes at 0
    
    return skipfoot
