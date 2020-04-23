# -*- coding: utf-8 -*-
"""
Created on Sun Sep 24 01:50:17 2017

@author: eric.lenning
"""
from operator import itemgetter
#from numpy import mean

class SPCSounding(object):
    """Represents data from a FSL-formatted RAOB."""
    
    def __init__(self, site, date, hour):
        self.site = site  # three-letter station ID
        self.date = date  # YYMMDD
        self.hour = hour  # HHMM in UTC
        
        self.levels = []
        
        self.cleaned = False
    
    # Assume checks for bad data have already been performed, as have 
    # conversions to proper units.
    
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
 
    def sortAndClean(self, unsortedLevels):
        # Sort to get ascending heights and descending pressures
        sortedLevels = sorted(unsortedLevels, key=itemgetter('hght'))
        sortedLevels = sorted(sortedLevels, key=itemgetter('mb'),reverse=True)
        
        # now seek to merge duplicate lines (same pressure levels)
        
        prevLevel = {}
        dups = []
        
        for thisLevel in sortedLevels:
            
            if (bool(prevLevel) and (thisLevel['mb'] == prevLevel['mb'])):
                
                thisLevel['hght'] = .5 * (thisLevel['hght'] + prevLevel['hght']) # hght should always be valid
                
                for k in thisLevel.keys():
                    if (thisLevel[k] > -9999):
                        if (prevLevel[k] > -9999):
                            #if (thisLevel[k] != prevLevel[k]):
                            #    print "Averaging {}: {} {}".format(k, thisLevel[k], prevLevel[k])                                
                            thisLevel[k] = 0.5 * (thisLevel[k] + prevLevel[k])
                    elif (prevLevel[k] > -9999):
                        thisLevel[k] = prevLevel[k]
                
                dups.append(prevLevel)
                
            prevLevel = thisLevel
        
        for d in dups:
            sortedLevels.remove(d)
        
        # Now do the same to remove repeated heights
        prevLevel = {}
        dups = []
        
        for thisLevel in sortedLevels:
            
            if (bool(prevLevel) and (thisLevel['hght'] <= prevLevel['hght'])):

                 # in case next level is higher than this but lower than previous
                thisLevel['hght'] = prevLevel['hght']
                
                dups.append(thisLevel)
                
            prevLevel = thisLevel
        
        for d in dups:
            sortedLevels.remove(d)
            
        self.cleaned = True
        
        return sortedLevels
            
    def printSPCSounding(self,outfh):
                
        outfh.write("%TITLE%\n")
        outfh.write(" {}   {}/{}\n".format(self.site, self.date, self.hour))
        outfh.write("\n")
        outfh.write("   LEVEL       HGHT       TEMP       DWPT       WDIR       WSPD\n")
        outfh.write("-------------------------------------------------------------------\n")
        outfh.write("%RAW%\n")        
        
        if (not self.cleaned):
            sortedLevels = self.sortAndClean(self.levels)
        
        for thisLevel in sortedLevels:
            
            levelString = "{:.2f},".format(thisLevel['mb'])
            
            outfh.write(" {:<10s}{:9.2f},  {:8.2f},  {:8.2f},  {:8.2f},  {:8.2f}\n".format(
                    levelString, thisLevel['hght'], thisLevel['temp'], 
                    thisLevel['dwpt'], thisLevel['wdir'], thisLevel['wspd']))            
        
        outfh.write("%END%\n")
    
    def getSPCSoundingText(self):
        soundingText = "%TITLE%\n"
        soundingText += " {}   {}/{}\n".format(self.site, self.date, self.hour)
        soundingText += "\n"
        soundingText += "   LEVEL       HGHT       TEMP       DWPT       WDIR       WSPD\n"
        soundingText += "-------------------------------------------------------------------\n"
        soundingText += "%RAW%\n"
        
        if (not self.cleaned):
            sortedLevels = self.sortAndClean(self.levels)
        
        for thisLevel in sortedLevels:
            
            levelString = "{:.2f},".format(thisLevel['mb'])
            
            soundingText += " {:<10s}{:9.2f},  {:8.2f},  {:8.2f},  {:8.2f},  {:8.2f}\n".format(
                    levelString, thisLevel['hght'], thisLevel['temp'], 
                    thisLevel['dwpt'], thisLevel['wdir'], thisLevel['wspd'])    
        
        soundingText += "%END%\n"
        
        return soundingText
    