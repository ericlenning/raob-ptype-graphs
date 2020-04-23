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
#import sys
import re
import datetime
import numpy as np
#from SPCSounding import SPCSounding
import os.path
from StringIO import StringIO
import glob
import sharppy.sharptab.profile as profile
import sharppy.sharptab.thermo as thermo
import sharppy.sharptab.interp as interp
import sharppy.sharptab.utils as utils
import sharppy.sharptab.params as params
import cPickle as pickle
from spcProfile import spcExpProfile
     
# Return pressure levels at bottom and top of DGZ as defined by the portion of
# the sounding between -12 and -18 C.  For this we will use the highest -12 and
# last occurrence of -18 above that level.
    
def dgz(prof):
    
    baseT = -12
    topT = -18
    
    above_m12 = [i for i,v in enumerate(prof.tmpc) if v > baseT]
    above_m18 = [i for i,v in enumerate(prof.tmpc) if v > topT]
    
    lastMinus12 = above_m12[-1] # next(x[0] for x in enumerate(prof.tmpc) if x[1] <= -12)
    
    #firstMinus18 = next(x[0] for x in enumerate(prof.tmpc[firstMinus12:]) if x[1] <= -17)
    
    ind = lastMinus12

    # We know the index of the last temperature warmer than -12.  Interpolate
    # between this and the first temperature colder than that (at ind + 1) to
    # find the pressure at -12 itself.
    
    pbot = np.power(10, np.interp(baseT, [prof.tmpc[ind+1], prof.tmpc[ind]],
                            [prof.logp[ind+1], prof.logp[ind]]))    

    lastMinus18 = above_m18[-1]
        
    ind = lastMinus18 # + lastMinus12

    # It is possible that the last temperature warmer than -12 is also the last
    # warmer than -18.  This should be okay, because we are interpolating to get
    # the actual pressure at -12 and -18.
    
    #if ind == lastMinus12:
    #    ind += 1

    # We know the index of the last temperature warmer than -18.  Interpolate
    # between this and the first temperature colder than that (at ind + 1) to
    # find the pressure at -18 itself.
        
    ptop = np.power(10, np.interp(topT, [prof.tmpc[ind+1], prof.tmpc[ind]],
                            [prof.logp[ind+1], prof.logp[ind]]))
    
    if pbot < ptop:
        print "DGZ error"
    
    return pbot, ptop
           
# parse text for one SPC sounding and return tuple containing p, h, T, Td, wdir, wspd
# lists for the sounding

def parseSPC(spc_text):

    ## split the text based on newlines
    data = np.array([l.strip() for l in spc_text.split('\n')])

    loc = ""
    datestring = ""
    
    m = re.match("(\w+)\s+(\d+/\d+)",data[1])

    if (m):
        print m.groups()
        loc, datestring = m.groups()

    ## necessary index points
    title_idx = np.where( data == '%TITLE%')[0][0]
    start_idx = np.where( data == '%RAW%' )[0] + 1
    finish_idx = np.where( data == '%END%')[0]
    
    ## put it all together for StringIO
    level_data = '\n'.join(data[start_idx[0] : finish_idx[0]][:])
    sound_data = StringIO( level_data )

    ## read the data into arrays
    p, h, T, Td, wdir, wspd = np.genfromtxt( sound_data, delimiter=',', comments="%", unpack=True )

    return p, h, T, Td, wdir, wspd, loc, datestring

def insertLevels(prof, zeroHt, zeroPres, level):   
    
    prof.dwpc = np.insert(prof.dwpc,level,
                          [interp.dwpt(prof,zeroPres)])
    prof.vtmp = np.insert(prof.vtmp,level,
                          [interp.vtmp(prof,zeroPres)])
    prof.thetae = np.insert(prof.thetae,level,
                          [interp.thetae(prof,zeroPres)])
    #prof.wetbulb = np.insert(prof.wetbulb,level,
    #                      [interp.generic_interp_pres(np.log10(zeroPres), prof.logp[::-1], prof.wetbulb[::-1])])
    try:
        dir,mag = interp.vec(prof,zeroPres)
        prof.wdir = np.insert(prof.wdir,level,[dir])
        prof.wspd = np.insert(prof.wspd,level,[mag])
        prof.u, prof.v = utils.vec2comp(prof.wdir, prof.wspd)
    except:
        prof.wdir = np.insert(prof.wdir,level,[0])
        prof.wspd = np.insert(prof.wspd,level,[0])
        prof.u, prof.v = utils.vec2comp(prof.wdir, prof.wspd)        

    prof.hght = np.insert(prof.hght,level,[zeroHt])
    prof.pres = np.insert(prof.pres,level,[zeroPres])
    prof.logp = np.log10(prof.pres.copy())    
    
    return prof
    
def interpZeroTemps(prof):
    
    tLength = len(prof.tmpc)
    
    level = 1    
    
    while (level < tLength):
        prevLevel = level - 1
        nextLevel = level + 1
        
        prevT = prof.tmpc[prevLevel]
        thisT = prof.tmpc[level]
        
        if (thisT * prevT < 0):
            # print "Crossed zero:", prevT, thisT
            
            zeroHt = -9999
            
            if (thisT < prevT):
                zeroHt = np.interp(0,np.flipud(prof.tmpc[prevLevel:nextLevel]),np.flipud(prof.hght[prevLevel:nextLevel]))
            else:
                zeroHt = np.interp(0,prof.tmpc[prevLevel:nextLevel],prof.hght[prevLevel:nextLevel])
            
            zeroPres = interp.pres(prof,zeroHt)
            
            prof.tmpc = np.insert(prof.tmpc,level,[0])        

            prof.wetbulb = np.insert(prof.wetbulb,level,
                          [interp.generic_interp_pres(np.log10(zeroPres), prof.logp[::-1], prof.wetbulb[::-1])])
            
            prof = insertLevels(prof,zeroHt,zeroPres,level)
            
            tLength += 1
            
            # Double increment j since you just added a new element
            level = level + 2
        else:
            level = level + 1

    return prof
    
def interpZeroWetbulbs(prof):
    
    tLength = len(prof.wetbulb)
    
    level = 1    
    
    while (level < tLength):
        prevLevel = level - 1
        nextLevel = level + 1
        
        prevT = prof.wetbulb[prevLevel]
        thisT = prof.wetbulb[level]
        
        if (thisT * prevT < 0):
            # print "Crossed wb zero:", prevT, thisT
            
            zeroHt = -9999
            
            if (thisT < prevT):
                zeroHt = np.interp(0,np.flipud(prof.wetbulb[prevLevel:nextLevel]),np.flipud(prof.hght[prevLevel:nextLevel]))
            else:
                zeroHt = np.interp(0,prof.wetbulb[prevLevel:nextLevel],prof.hght[prevLevel:nextLevel])
            
            zeroPres = interp.pres(prof,zeroHt)
            
            prof.wetbulb = np.insert(prof.wetbulb,level,[0])
            
            prof.tmpc = np.insert(prof.tmpc,level,
                          [interp.temp(prof,zeroPres)])
             
            prof = insertLevels(prof,zeroHt,zeroPres,level)                       
            
            tLength += 1
            
            # Double increment j since you just added a new element
            level = level + 2
        else:
            level = level + 1
            
    return prof
    
def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)
        
pathname = os.path.normpath("C:/Users/eric.lenning/Python/Soundings/")
filename = "spc*"
spcfiles = os.path.normpath(pathname + "/2 - spcFormat\\" + filename)

for spcFilename in glob.glob(spcfiles): # pathname + filename):
    print "Reading from", spcFilename
    spc_path, spc_name = os.path.split(os.path.abspath(spcFilename))
    outFilename1 = os.path.normpath(pathname + "/4 - pklFormat\\pkl_" + spc_name)
    outFilename2 = os.path.normpath(pathname + "/3 - expanded spcFormat\\exp_" + spc_name)
    outFilename3 = os.path.normpath(pathname + "/5 - dgzStats\\dgz_" + spc_name)
    print "Writing to", outFilename1, outFilename2, outFilename3
    
    expFile = open(outFilename2, "w")
    dgzFile = open(outFilename3, "w")
    
    with open(spcFilename, "r") as myfile:

        spc_data=myfile.read()
        
        data = np.array([l.strip() for l in spc_data.split('\n')])    
        
        sounding_text = ""
        
        profiles = []
        
        levels = 0               
        
#        prof = spcBasicProfile()
        
        for d in data:
        
            levels += 1
            
            sounding_text += d + "\n"
            
            m = re.match("%END%",d)
            
            if m:
                if (levels > 30):
                    #print sounding_text
                    
                    pres, hght, tmpc, dwpc, wdir, wspd, loc, datestring = parseSPC(sounding_text)                                       
                    
                    # Somewhat arbitrary choice for minimum number of levels in a sounding.
                    # But also require the top-level temperature to be below -20.
                    
                    if (tmpc[-1] < -20) and (tmpc[0] > -12): #  and (m15_rh > 50):
                        
                        # Do a rough estimate of RH in the DGZ to weed out dry-ish soundings.
                        m15 = next(x[0] for x in enumerate(tmpc) if x[1] <= -15)
                        
                        m15_rh = thermo.relh(pres[m15], tmpc[m15], dwpc[m15])     

                        if m15_rh > 50:                                   
                            
                            date = datetime.datetime.strptime(datestring, "%m%d%y/%H%M")
                            
                            # this returns a SHARPPy profile, including theta-e and wetbulb                        
                            prof = profile.create_profile(profile='default', pres=pres, hght=hght, tmpc=tmpc, \
                                            dwpc=dwpc, wspd=wspd, wdir=wdir, missing=-9999, date=date, \
                                            location=loc, strictQC=True)
                        
                            prof = interpZeroTemps(prof)
                            
                            prof = interpZeroWetbulbs(prof)
                            
                            #mean_rh_tmp = params.mean_relh(prof)
                            mean_rh = params.mean_relh(prof, ptop=400)
                            
                            # print mean_rh
                            
                            #if prof.tmpc[0] > -12:
                                
                            #dgz_pbot, dgz_ptop = params.dgz(prof)
                            dgz_pbot, dgz_ptop = dgz(prof)
                            
                            # If DGZ properly decreases in pressure, meaning both 
                            # the top and bottom were found, proceed.
                            
                            if dgz_pbot > dgz_ptop:
                                
                                dgz_meanrh = params.mean_relh(prof, pbot=dgz_pbot, ptop=dgz_ptop)
                                
                                dgzFile.write("{} {} {} {} {}\n".format(
                                        loc, datestring, 
                                        dgz_pbot, dgz_ptop, dgz_meanrh))
                                
                                # print dgz_meanrh
                                
                                if (mean_rh > 80) and (dgz_meanrh > 75):
                                    
                                    #if m15_rh <= 50:
                                    #    print "*** DRY m15_rh PASSED: ", m15_rh
                                    
                                    profiles.append(prof)
                                
                                    sProf = spcExpProfile(prof)
                                
                                    sProf.toExpandedFile(expFile)
                                    
                                #else:
                                    #print "   Failed: Sounding RH {} or DGZ RH {} too dry. (m15_rh = {})".format(
                                    #        mean_rh, dgz_meanrh, m15_rh)
                                    
                            else:
                                print "   Failed: Invalid DGZ."
                            
                    #else:
                        #print "   Failed: SfcT {} TopT {} too warm or minus 15 RH {} too dry.".format(
                        #        tmpc[0], tmpc[-1], m15_rh)
                        
                #else:
                    #print "   Failed: Only {} levels.".format(levels)
                    
                sounding_text = ""

                levels = 0            
    
    expFile.close()
    
    dgzFile.close()
    
    save_object(profiles,outFilename1)
    
