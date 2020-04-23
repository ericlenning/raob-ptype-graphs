# -*- coding: utf-8 -*-
"""
Created on Sun Sep 24 01:50:17 2017

@author: eric.lenning
"""

class FSLSounding:
    """Represents data from a FSL-formatted RAOB."""
    
    def __init__(self, site, date, hour):
        self.site = site  # three-letter station ID
        self.date = date  # YYMMDD
        self.hour = hour  # HHMM in UTC
        
        self.levels = []
    
    def addLevel(self, mb, hght, temp, dwpt, wdir, wspd):
        level = {
            "mb": mb,
            "hght" : hght,
            "temp" : temp,
            "dwpt" : dwpt,
            "wdir" : wdir,
            "wspd" : wspd,
        }
        
        self.levels.append(level)

    """        
    %TITLE%
    OAX   140616/1900 

       LEVEL       HGHT       TEMP       DWPT       WDIR       WSPD
       -------------------------------------------------------------------
       %RAW%
       1000.00,    34.00,  -9999.00,  -9999.00,  -9999.00,  -9999.00
       965.00,    350.00,     27.80,     23.80,    150.00,     23.00
       """
 
    def printSPCSounding(self):
        print "%TITLE%"
        print(" {}   {}/{:02d}".format(self.site, self.date, int(self.hour)))
        print ""
        print "   LEVEL       HGHT       TEMP       DWPT       WDIR       WSPD"
        print "-------------------------------------------------------------------"
        print "%RAW%"
        
        for l in self.levels:
            
            print(" {:<7.2f}, {:8.2f}, {:8.2f}, {:8.2f}, {:8.2f}, {:8.2f}".format(
                  l['mb'], l['hght'], l['temp'], l['dwpt'], l['wdir'], l['wspd']))
        
        print "%END%"
        
        
    