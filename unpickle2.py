# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 08:55:06 2017

@author: eric.lenning
"""
import glob
import csv
import os
from collections import defaultdict
import cPickle as pickle
#import sharppy.sharptab.profile
import sounding_functions as sf
# import numpy as np
import plotting
# import matplotlib.pyplot as plt
from EnergyProfile import EnergyProfile

def getStations(stns, idCol, numCol, stnfile, headerCount, splitNeeded):
    
    with open(stnfile) as csvfile:
        for i in range(headerCount):
            csvfile.next()  
        reader = csv.reader(csvfile)
        for nextrow in reader:
            row = nextrow
            if splitNeeded:
                row = nextrow[0].split()
            if len(row) > 1:
                #print row[0], row
                #if row[1].startswith('71') or row[numCol].startswith('72'):
                num = row[numCol].strip()
                id = row[idCol].strip()
                if ("-" not in id) and (num.startswith('7')):
                    #print id, row
                    if len(id) > 3:
                        if id.startswith('C') or id.startswith('K'):
                            if id[1:] in stns:
                                if (stns[id[1:]] != num):
                                    print id, row
                                    print "Inconsistent Duplicate:", id[1:], stns[id[1:]], num
                            stns[id[1:]] = [num, id]
                    else:
                        if id in stns:
                            if (stns[id] != num):
                                print id, row
                                print "Inconsistent Duplicate:", id, stns[id], num
                        stns[id] = [num, id]
                    
    return stns
"""
pathname = os.path.normpath("C:/Users/eric.lenning/Python/Soundings/")
filename = "/raob*"

for fslFilename in glob.glob(pathname + filename):
    print "Reading from", fslFilename
    fsl_path, fsl_name = os.path.split(os.path.abspath(fslFilename))
    outFilename = os.path.normpath(fsl_path + "\\spc" + fsl_name)
    print "Writing to", outFilename
    parseFSL(fslFilename,outFilename)
"""

pathname = os.path.normpath("C:/Users/eric.lenning/Python/Soundings/4 - pklFormat")

# First get list of RAOB IDs and names since that's not part of the SPC format.
#
# The WMO identifier, often called the "index number" relies on a 5-digit
# numeric code to identify a land weather station. The first two digits are
# referred to as the "block number" and refer to the geographic area 
# (00-29 Europe, 30-59 Asia, 60-68 Africa, 69 special use, 70-79 North America,
# 80-89 South America, 90-99 Oceania).

stns = defaultdict()

stnfile1 = os.path.normpath(pathname + "/raob-short.txt")

stns = getStations(stns, 0, 1, stnfile1, 4, True)

stnfile2 = os.path.normpath(pathname + "/raob-stn.txt")

stns = getStations(stns, 1, 0, stnfile2, 9, False)

# headers = next(reader)[8:]

#stns = dict((rows[1],rows[0]) for rows in reader) 
#if (rows[1].startswith('71') or rows[1].startswith['72']))

# Now get to the pickled files

filename = "/pkl*"

#freezingSfcWb = 0
#freezingSfcTmpc = 0

#maxtmpcWarmLayerAloftCount = 0
#maxWbWarmLayerAloftCount = 0

maxWarmLayerAloftCount = defaultdict(int)

profileCount = 0

# profileLevels = 0
minLevels = 9999
maxLevels = 0

tmpcLayerCounts = defaultdict(int)
wbLayerCounts = defaultdict(int)

"""
freezingLevels = []
wbFreezingLevels = []
sfcTmpcEnergies = []
sfcWbEnergies = []
sfcTemps = []
sfcWbTemps = []
"""
# minEnergies = []
"""
aloftMaxTs = []
aloftMaxWbs = []
minRefreezeTs = []
minRefreezeTws = []
posTmpcEnergies = []
posWbEnergies = []
refreezeTmpcEnergies = []
refreezeWbEnergies = []
warmTAloftDepths = []
warmTwAloftDepths = []
refreezeTDepths = []
refreezeTwDepths = []
tmpcPtypes = []
wbPtypes = []
"""

# Read each pickled file, with pkl prefix.  Each file will contain a list of
# RAOB profiles.

tmpcStats = defaultdict(list)
wbStats = defaultdict(list)
# stats = defaultdict(list)
                    
for pklFilename in glob.glob(pathname + filename):
    print "Reading from", pklFilename
    pathname, pkl_name = os.path.split(os.path.abspath(pklFilename))
    # outFilename = os.path.normpath(pathname + "\\stats" + pkl_name)
    # print "Writing to", outFilename
        
    with open(pklFilename, 'rb') as input:

        # outfh = open(outFilename, 'w')
        
        new_profiles = pickle.load(input)
        
        # Process each profile in the file.
        
        for p in new_profiles:
            
            # Only keep those with a surface temperature less than 4 C since
            # those are of greatest interest for mixed precipitation.  Any
            # warmer and you're very likely just dealing with rain.
            
            if (p.tmpc[0] < 4):
                
                # Keep track of the number of profiles processed.
                
                profileCount += 1
                    
                # levels = len(p.pres)
    
                # Track basic information about all profiles
                
                # profileLevels += len(p.pres)
    
                # Track the least and most number of levels found in a profile.
                # Note that no profile was stored in the spc2profiles script
                # if there were less than 30 levels.
                
                # if levels < minLevels:
                #    minLevels = levels
    
                # if levels > maxLevels:
                #    maxLevels = levels
                
                tmpcProfile = EnergyProfile(pres=p.pres, hght=p.hght, temp=p.tmpc,
                                            sfct=p.tmpc[0], sfctw=p.wetbulb[0], tstr="T")
                
                wbProfile = EnergyProfile(pres=p.pres, hght=p.hght, temp=p.wetbulb,
                                          sfct=p.tmpc[0], sfctw=p.wetbulb[0], tstr="Tw")
    
                # Build lists for properties that all profiles have.
                
                for prof in 'T', 'Tw':                    
                    
                    tempProfile = tmpcProfile                    
                    stats = tmpcStats
                    
                    if prof == 'Tw':
                        tempProfile = wbProfile
                        stats = wbStats
                    
                    stats['stns'].append([stns[p.location][0],stns[p.location][1]])
                    
                    stats['dates'].append(p.date)
                    
                    stats['freezingLevels'].append(tempProfile.freezingLevel)                    
        
                    stats['sfcTemps'].append(tempProfile.sfcTemp)
                    
                    # print p.location, p.date
                        
                    stats['ptypes'].append(tempProfile.ptype)
                    
                    stats['ptypest'].append(tempProfile.ptypet)
                                        
                    stats['sfcEnergies'].append(tempProfile.sfcEnergy)                    
                    
                    # Now build lists for properties that only matter if there is
                    # warm air at the sfc or aloft in a profile (1 or more energy layers).
                    
                    """    
                    if tempProfile.energyLayerCount > 1:
                                                
                        warmAloftDepth = tempProfile.warmAloftDepth
    
                        refreezeLayerDepthMeters = tempProfile.refreezeLayerDepthMeters
                    """
                    
                    # Now build lists for properties that only matter if there is
                    # warm air aloft in a profile (1 or more energy layers).
                    
                    #if tempProfile.maxTempAloft > 0:
                        
                    stats['aloftMaxTemps'].append(tempProfile.maxTempAloft)
                    
                    stats['above900MaxTemps'].append(tempProfile.maxTempAbove900)
                    
                    stats['above2kftMaxTemps'].append(tempProfile.maxTempAbove2kft)
                    
                    stats['posEnergyAbove900'].append(tempProfile.posEnergyAbove900)
                    
                    stats['posEnergyAbove2kft'].append(tempProfile.posEnergyAbove2kft)                    
                    
                    stats['posEnergyBelow2kft'].append(tempProfile.posEnergyBelow2kft)
                    
                    stats['negEnergyBelow850'].append(tempProfile.negEnergyBelow850)
                    
                    stats['posEnergies'].append(tempProfile.posEnergyAloft)
                    
                    stats['posTotals'].append(tempProfile.posEnergyTotal)
                    
                    stats['warmAloftDepths'].append(tempProfile.warmAloftDepth)
                
                    stats['minRefreezeTemps'].append(tempProfile.minRefreezeTemp)
                    
                    stats['refreezeEnergiesLowest'].append(tempProfile.refreezeEnergyLowest)
                    
                    stats['refreezeEnergiesAll'].append(tempProfile.refreezeEnergyAll)

                    stats['refreezeEnergiesSfc'].append(tempProfile.refreezeEnergySfc)                   
                    
                    stats['refreezeDepths'].append(tempProfile.refreezeLayerDepthMeters)
                                    
                    if tempProfile.sfcTemp > 0:
                        stats['minEnergies'].append(sf.calcMinEnergy(tempProfile.sfcTemp))
                    else:
                        stats['minEnergies'].append(0)
                        
                    if prof == 'T':
                        tmpcStats = stats
                    else:
                        wbStats = stats                        
                   
                    warmLayerAloftCount = tempProfile.energyLayerCount / 2
                    
                    if warmLayerAloftCount > maxWarmLayerAloftCount[prof]:
                        maxWarmLayerAloftCount[prof] = warmLayerAloftCount
                
        
        """        
        for i, v in enumerate(stats['sfcTemps']):
            outfh.write("{} {} {} {} {} {} {}\n".format(
                        stats['sfcTemps'][i], stats['sfcTmpcEnergies'][i], 
                        stats['minEnergies'][i], stats['freezingLevels'][i],
                        stats['sfcWbTemps'][i], stats['sfcWbEnergies'][i], 
                        stats['wbFreezingLevels'][i]))

        outfh.close()
        """
        
print profileCount, "profiles"

for t in 'T', 'Tw':

    tempLabel = "Temperature"
    stats = tmpcStats    
    
    if t == "Tw":
        tempLabel = "Wet Bulb Temperature"
        stats = wbStats
    
#    plotting.ptypes_on_sfct_vs_fzlevel(stats['ptypes'], stats['sfcTemps'], stats['freezingLevels'], t)
    
    #tFilt = sTemps[sTemps > 32]
    #eFilt = sEnergies[sTemps > 32]
    
    #tFilt = tFilt[tFilt < 40]
    #eFilt = eFilt[tFilt < 40]
    
    #plotting.sfct_vs_energy(stats['sfcTemps'], stats['sfcEnergies'], stats['minEnergies'], tempLabel, stats['stns'], stats['dates'])    

    plotting.hist(stats['posEnergies'], "Positive {}-Energy Aloft".format(t), t)
    
    plotting.hist(stats['aloftMaxTemps'], "Max{} Aloft".format(t), t)
    
    plotting.hist(stats['above900MaxTemps'], "Max{} Above 900mb".format(t), t)

    plotting.hist(stats['above2kftMaxTemps'], "Max{} Above 2kft".format(t), t)    
    
    plotting.hist(stats['refreezeEnergiesLowest'], "Refreezing {}-Energy".format(t), t)
    
    plotting.hist(stats['minRefreezeTemps'], "Min {}-Energy".format(t), t)
    
    """
    plotting.t_vs_posEnergy(stats['sfcTemps'], stats['sfcEnergies'], stats['minEnergies'], tempLabel, stats['stns'], stats['dates'], xMax=4.44, yMax=80, snowRainTemp=1.11, rainOnlyTemp=2.78, tempLayer="Sfc", energyLayer="Sfc", ptypes=stats['ptypes'])
    """
    
    plotting.maxt_vs_energy(stats['aloftMaxTemps'], stats['posEnergies'], stats['sfcTemps'], t, stats['stns'], stats['dates'])
    

    #plotting.t_vs_posEnergy(stats['aloftMaxTemps'], stats['posEnergies'], [0]*len(stats['aloftMaxTemps']), t, stats['stns'], stats['dates'], xMax=4, yMax=80, snowRainTemp=1, rainOnlyTemp=3, tempLayer="Aloft", energyLayer="Aloft", ptypes=stats['ptypes'])
    
    plotting.t_vs_posEnergy(stats['aloftMaxTemps'], stats['posEnergies'], [0]*len(stats['aloftMaxTemps']), t, stats['stns'], stats['dates'], xMax=4, yMax=80, snowRainTemp=1, rainOnlyTemp=3, tempLayer="Aloft", energyLayer="Aloft", ptypes=stats['ptypes'], ptypest=stats['ptypest'])
      
    plotting.t_vs_lowEnergy(stats['sfcTemps'], stats['posEnergyBelow2kft'], stats['posEnergyAbove2kft'], [0]*len(stats['sfcTemp']), t, stats['stns'], stats['dates'], xMax=4, yMax=80, snowRainTemp=1, rainOnlyTemp=3, tempLayer="Sfc", energyLayer="Below2kft", ptypes=stats['ptypes'], ptypest=stats['ptypest'])

    """
    plotting.t_vs_posEnergy(stats['sfcTemps'], stats['posTotals'], [0]*len(stats['sfcTemps']), t, stats['stns'], stats['dates'], xMax=4, yMax=80, snowRainTemp=1, rainOnlyTemp=3, tempLayer="Sfc", energyLayer="Total", ptypes=stats['ptypes'])

    plotting.aloft_vs_above(stats['aloftMaxTemps'], stats['above900MaxTemps'], "Max {}".format(t), "900mb", stats['stns'], stats['dates'])
    """
    plotting.aloft_vs_above(stats['aloftMaxTemps'], stats['above2kftMaxTemps'], "Max {}".format(t), "2000ft", stats['stns'], stats['dates'], stats['ptypes'], stats['ptypest'])
    
    plotting.aloft_vs_total(stats['posEnergies'], stats['posTotals'], stats['refreezeEnergiesLowest'], stats['sfcTemps'], "Aloft vs Total", stats['stns'], stats['dates'])
    
    plotting.aloft_vs_above(stats['posEnergies'], stats['posEnergyAbove900'], "{} Energy".format(t), "900mb", stats['stns'], stats['dates'], stats['ptypes'], stats['ptypest'])
    
    plotting.aloft_vs_above(stats['posEnergies'], stats['posEnergyAbove2kft'], "{} Energy".format(t), "2000ft", stats['stns'], stats['dates'], stats['ptypes'], stats['ptypest'])    

    plotting.aloft_not_above(stats['posEnergies'], stats['posEnergyAbove2kft'], "{} Energy".format(t), "2000ft", stats['stns'], stats['dates'], stats['ptypes'], stats['ptypest'])    
 
    plotting.total_vs_below(stats['refreezeEnergiesLowest'], stats['negEnergyBelow850'], "{} Energy".format(t), "850mb", stats['stns'], stats['dates'], stats['ptypes'], stats['ptypest'])
    
    plotting.total_vs_below(stats['refreezeEnergiesLowest'], stats['refreezeEnergiesAll'], "{} Energy".format(t), "Highest Melting Layer", stats['stns'], stats['dates'], stats['ptypes'], stats['ptypest'])    
    
    plotting.maxt_vs_depth(stats['aloftMaxTemps'], stats['warmAloftDepths'], tempLabel, stats['stns'], stats['dates'])
    """
    maxTs, depths = zip(*[(i,j) for i,j in zip(stats['aloftMaxTemps'], stats['warmAloftDepths']) if i <= 4])
    
    plotting.maxt_vs_depth(maxTs, depths, tempLabel)
    
    plotting.maxtXdepth_vs_energies(stats['aloftMaxTemps'], stats['warmAloftDepths'], stats['posEnergies'], t)
    
    # plotting.energy_vs_depth(posTmpcEnergies, warmTAloftDepths, tempLabel)
    
    # plotting.energy_vs_depth(posWbEnergies, warmTwAloftDepths, "Wet Bulb")
    
    plotting.depth_vs_energy(stats['warmAloftDepths'], stats['posEnergies'], tempLabel)
    """
    plotting.t_vs_refreeze_depth(stats['minRefreezeTemps'], stats['refreezeDepths'], stats['refreezeEnergiesLowest'], t, stats['aloftMaxTemps'], stats['stns'], stats['dates'],ptypes=stats['ptypes'])
    """
    plotting.maxt_vs_mint(stats['aloftMaxTemps'], stats['minRefreezeTemps'], t)
    
    """
    plotting.maxenergy_vs_minenergy(stats['posEnergies'], stats['refreezeEnergiesLowest'], t)
    """
    plotting.sfct_vs_energyAloft(stats['sfcTemps'], stats['posEnergies'], t)    
    
    plotting.aloft_vs_sfc(stats['posEnergies'], stats['sfcEnergies'], stats['refreezeEnergiesLowest'], t)
    """ 
    
    plotting.rare_sleet(stats['sfcTemps'], stats['refreezeEnergiesLowest'], stats['aloftMaxTemps'], t, stats['stns'], stats['dates'], ptypes=stats['ptypes'])
    
    plotting.cold_not_snow(stats['sfcTemps'], stats['posEnergies'], stats['above2kftMaxTemps'], t, stats['stns'], stats['dates'], ptypes=stats['ptypes'])
    
    plotting.max_sleet(stats['posEnergies'], stats['refreezeEnergiesLowest'], stats['aloftMaxTemps'], stats['sfcTemps'], t, stats['stns'], stats['dates'], ptypes=stats['ptypes'], ptypest=stats['ptypest'])
    
    # print "Overall temp and wb layer counts:", tmpcLayerCounts, wbLayerCounts
    
    # print "Most levels", maxLevels
    # print "Least levels", minLevels
    
    # print "Average number of levels per profile:", (profileLevels / profileCount)
    print "Max {} warm layer count: {}".format(t, maxWarmLayerAloftCount[t])
