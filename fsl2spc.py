# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 11:01:20 2017

@author: eric.lenning
"""

"""
FSL Rawinsonde data format - https://ruc.noaa.gov/raobs/fsl_format-new.html

The official FSL data format is similar to the format used by the National
Severe Storms Forecast Center (NSSFC) in Kansas City.  The first 4 lines of 
the sounding are identification and information lines. All additional lines 
are data lines.  An entry of 32767 (original format) or 99999 (new format) 
indicates that the information is either missing, not reported, or not 
applicable.

                            ---COLUMN NUMBER---

      1          2          3          4          5          6           7
 LINTYP
                                header lines
    254        HOUR        DAY      MONTH       YEAR    (blank)     (blank)
      1       WBAN#       WMO#        LAT D      LON D     ELEV       RTIME
      2       HYDRO       MXWD      TROPL      LINES     TINDEX      SOURCE
      3     (blank)      STAID    (blank)    (blank)      SONDE     WSUNITS
                        
                                data lines
      9    PRESSURE     HEIGHT       TEMP      DEWPT   WIND DIR    WIND SPD
      4
      5
      6
      7
      8    
      



      5   9490    518     96     16  99999  99999
      6   9387    609  99999  99999    320     21
      
And SPC data format:

%TITLE%
 OAX   140616/1900 

   LEVEL       HGHT       TEMP       DWPT       WDIR       WSPD
-------------------------------------------------------------------
%RAW%
 1000.00,    34.00,  -9999.00,  -9999.00,  -9999.00,  -9999.00
 965.00,    350.00,     27.80,     23.80,    150.00,     23.00
 962.00,    377.51,     27.40,     22.80,  -9999.00,  -9999.00
 936.87,    610.00,     25.51,     21.72,    145.00,     29.99
 925.00,    722.00,     24.60,     21.20,    150.00,     33.99
...
%END%

"""
import sys
import re
import numpy as np
from SPCSounding import SPCSounding
import os.path
import glob
import calendar
abbr_to_num = {name.upper(): num for num, name in enumerate(calendar.month_abbr) if num}
                
def addToSPCSounding(sounding, groups):
    # Pressure and height (m) should always be present
    mb = int(groups[1])/10.0
    
    hght = int(groups[2])
    
    tmpc = int(groups[3])/10.0        
            
    dptc = int(groups[4])/10.0

    # Don't add level if pressure or height are missing
    # OR if tmpc or dptc are missing
    
    if ((mb == 9999.9) or (hght == 99999) or (tmpc == 9999.9) or (dptc == 9999.9)):
        return False
            
    wdir = int(groups[5])
                
    wspd = int(groups[6])
                      
    if wspd == 99999:
        wspd = -9999        
            
    if wdir == 99999:
        wdir = -9999
        
    sounding.addLevel(mb, hght, tmpc, dptc, wdir, wspd)
    
    return True

def parseFSL(fsl_filename,out_filename):
    fsl_file = open(fsl_filename,'r').read()
    
    out_file = open(out_filename,'w')
    
    data = np.array([l.strip() for l in fsl_file.split('\n')])    
    
    hour = 99
    day = 99
    month = 'MMM'
    year = 9999
    site = 'MMM'
    lat = ""
    lon = ""
    elev = ""
    
    sounding = None
    levelsAdded = False
    
    for d in data:
        #254     12      1      DEC    2016
        m = re.match("254\s+(\d+)\s+(\d+)\s+(\w+)\s+(\d{4})",d)
        if (m):
            
            # If a sounding already exists this '254' match is a new sounding.
            # Create a profile and analyze the sounding just completed.
            
            if (sounding and levelsAdded):
                sounding.printSPCSounding(out_file)    
                
            #This is redundant for the very first sounding but resets this object
            #when subsequent new soundings start.
            
            sounding = None
            levelsAdded = False
            
            hour, day, month, year = m.groups()
            
            hour = "{:02d}00".format(int(hour))
            
            month_num = abbr_to_num[month]
            yr = int(year) % 100
            
            date = "{:02d}{:02d}{:02d}".format(int(month_num),int(day),int(yr))
            
        #1  13723  72317  36.08N 79.95W   277   1103
        m = re.match("1\s+\d+\s+\d+\s+(\d+(\.\d*)?|\.\d+)[NS]\s*(\d+(\.\d*)?|\.\d+)[EW]\s+(\d+)\s+",d)
        if (m):              
            #print m.groups()
            lat, latDec, lon, lonDec, elev = m.groups() # not used?
            #print lat, lon, elev

        #2     70     93   1110    146  99999      3
        m = re.match("2\s+\d+\s+\d+\s+\d+\s+(\d+)\s+",d)
        if (m):                
            lines = m.groups()[0] # not used

        #3           GSO                99999     kt
        m = re.match("3\s+(\w+)\s+\d+\s+(\w{2})",d)
        if (m):                
           site, wsunits = m.groups() # not used  
           
           # Now we should have all we need for a new sounding    
           sounding = SPCSounding(site, date, hour)

        #9   9770    277    106     46    310     10
        #4  10000     79  99999  99999  99999  99999
        #5   9680    354    110     40  99999  99999
        m = re.match("([4-9])\s+(\d+)\s+(\d+)\s+(\-*\d+)\s+(\-*\d+)\s+(\d+)\s+(\d+)",d)
        if (m):      
            
            if (addToSPCSounding(sounding, m.groups())):
                levelsAdded = True
            
                levelCode = int(m.groups()[0])
  
                if (levelCode == 9):
                    sfcHght = int(m.groups()[2])
                        
                    if ((sfcHght != int(elev))):
                        print "For {} Stn and Sfc hght are: {} {}".format(site, elev, sfcHght)
                        # sys.exit()
    
    out_file.close()
"""
pathname = os.path.normpath("C:/Users/eric.lenning/Python/Soundings/")
filename = "spc*"
spcfiles = os.path.normpath(pathname + "/2 - spcFormat\\" + filename)

for spcFilename in glob.glob(spcfiles): # pathname + filename):
    print "Reading from", spcFilename
    spc_path, spc_name = os.path.split(os.path.abspath(spcFilename))
    outFilename1 = os.path.normpath(pathname + "/4 - pklFormat\\pkl_" + spc_name)
    outFilename2 = os.path.normpath(pathname + "/3 - expanded spcFormat\\exp_" + spc_name)
    print "Writing to", outFilename1, outFilename2
"""

pathname = os.path.normpath("C:/Users/eric.lenning/Python/Soundings/")
filename = "/raob*"
fslFiles = os.path.normpath(pathname + "/1 - fslFormat\\" + filename)

for fslFilename in glob.glob(fslFiles):
    print "Reading from", fslFilename
    fsl_path, fsl_name = os.path.split(os.path.abspath(fslFilename))
    outFilename = os.path.normpath(pathname + "/2 - spcFormat\\spc" + fsl_name)
    print "Writing to", outFilename
    parseFSL(fslFilename,outFilename)

"""    
try: 
    fslFilename = sys.argv[1]
    print fslFilename
except:
    print "Usage: {} <fslFilename>".format(sys.argv[0])
"""