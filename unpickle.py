# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 08:55:06 2017

@author: eric.lenning
"""
import glob
import os
from collections import defaultdict
import cPickle as pickle
#import sharppy.sharptab.profile
import sounding_functions as sf
# import numpy as np
import plotting
# import matplotlib.pyplot as plt

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

pathname = os.path.normpath("C:/Users/eric.lenning/Python/Soundings/")
filename = "/pkl*"

#freezingSfcWb = 0
#freezingSfcTmpc = 0

maxtmpcWarmLayerAloftCount = 0
maxWbWarmLayerAloftCount = 0

profileCount = 0
# profileLevels = 0
minLevels = 9999
maxLevels = 0

tmpcLayerCounts = defaultdict(int)
wbLayerCounts = defaultdict(int)

freezingLevels = []
wbFreezingLevels = []
sfcTempEnergies = []
sfcWbEnergies = []
sfcTemps = []
sfcWbTemps = []
minEnergies = []
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
ptypes = []

# Read each pickled file, with pkl prefix.  Each file will contain a list of
# RAOB profiles.

for pklFilename in glob.glob(pathname + filename):
    print "Reading from", pklFilename
    pathname, pkl_name = os.path.split(os.path.abspath(pklFilename))
    outFilename = os.path.normpath(pathname + "\\stats" + pkl_name)
    print "Writing to", outFilename
        
    with open(pklFilename, 'rb') as input:

        outfh = open(outFilename, 'w')
        
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
                
                freezingLevel = sf.get_freezing_level(p.hght, p.tmpc)
                wbFreezingLevel = sf.get_freezing_level(p.hght, p.wetbulb)
    
                freezingLevels.append(freezingLevel)
                wbFreezingLevels.append(wbFreezingLevel)
    
                sfcTemps.append(p.tmpc[0])
                sfcWbTemps.append(p.wetbulb[0])
    
                # Calculate the magnitude of each energy layer in a sounding.
                # - A sounding with zero energy layers starts and stays below
                #   freezing.
                # - A sounding with an odd number of energy layers starts above
                #   freezing and crosses the freezing level an odd number of
                #   times. The even indices (0,2,4, etc.) are warm layers.
                # - A sounding with an even number of energy layers starts below
                #   freezing and crosses the freezing level an even number of
                #   times. The odd indices (1,3,5, etc.) are warm layers.
                
                tmpcEnergy = sf.calculate_energy(p.hght, p.tmpc, p.pres)                
                
                # print p.location, p.date
                    
                ptypes.append(sf.calc_precip(tmpcEnergy, p.pres[-1]))
    
                wbEnergy = sf.calculate_energy(p.hght, p.wetbulb, p.pres)
                
                # Build lists to record the energy in the surface layer
                # for each profile
                
                tmpcLayers = len(tmpcEnergy)
                maxTAloft = -9999                
                tmpcPosEnergyAloft = 0
                
                if tmpcLayers > 0:
                    sfcTempEnergies.append(tmpcEnergy[0])
                    
                    if tmpcLayers > 1:
                        maxTAloft = sf.get_maxt_aloft(p.tmpc)
                        minRefreezeT = sf.get_min_refreeze_t(p.tmpc)
                        
                        tmpcPosEnergyAloft = sum(tmpcEnergy[1::2])
                        tmpcRefreezeEnergy = tmpcEnergy[0]
                        
                        if tmpcLayers % 2 == 1:
                            tmpcPosEnergyAloft = sum(tmpcEnergy[2::2])
                            tmpcRefreezeEnergy = tmpcEnergy[1]
                        
                        warmTAloftDepth = sf.getWarmLayerDepth(p.tmpc, p.hght)

                        refreezeTLayerDepth = sf.getRefreezeLayerDepth(p.tmpc, p.hght)
                    
                else:
                    sfcTempEnergies.append(0)
    
                if maxTAloft > 0:
                    
                    aloftMaxTs.append(maxTAloft)
                    
                    posTmpcEnergies.append(tmpcPosEnergyAloft)
                    
                    warmTAloftDepths.append(warmTAloftDepth)
                    
                    if tmpcPosEnergyAloft < 0:
                        print "Neg T Energy:", tmpcLayers, maxTAloft, tmpcPosEnergyAloft, p.tmpc
                
                    minRefreezeTs.append(minRefreezeT)
                    
                    refreezeTmpcEnergies.append(tmpcRefreezeEnergy)
                    
                    refreezeTDepths.append(refreezeTLayerDepth)
                    
                wbLayers = len(wbEnergy)
                maxWbAloft = -9999
                wbPosEnergyAloft = 0
                
                if wbLayers > 0:
                    sfcWbEnergies.append(wbEnergy[0])
                    
                    if wbLayers > 1:
                        maxWbAloft = sf.get_maxt_aloft(p.wetbulb)
                        minRefreezeWb = sf.get_min_refreeze_t(p.wetbulb)
                            
                        wbPosEnergyAloft = sum(wbEnergy[1::2])
                        wbRefreezeEnergy = wbEnergy[0]
                        
                        if wbLayers % 2 == 1:                        
                            wbPosEnergyAloft = sum(wbEnergy[2::2])
                            wbRefreezeEnergy = wbEnergy[1]
                            
                        warmTwAloftDepth = sf.getWarmLayerDepth(p.wetbulb, p.hght)
                        
                        refreezeTwLayerDepth = sf.getRefreezeLayerDepth(p.wetbulb, p.hght)
                else:
                    sfcWbEnergies.append(0)
    
                if maxWbAloft > 0:
                    
                    aloftMaxWbs.append(maxWbAloft)
                    
                    posWbEnergies.append(wbPosEnergyAloft)           
                    
                    warmTwAloftDepths.append(warmTwAloftDepth)
                    
                    if wbPosEnergyAloft < 0:
                        print "Neg Wb Energy:", wbLayers, maxWbAloft, wbPosEnergyAloft, p.tmpc
    
                    minRefreezeTws.append(minRefreezeWb)
                    
                    refreezeWbEnergies.append(wbRefreezeEnergy)
                    
                    refreezeTwDepths.append(refreezeTwLayerDepth)
                    
                #print "T", p.tmpc, maxTAloft, tmpcPosEnergyAloft
    
                #print "Tw", p.wetbulb, maxWbAloft, wbPosEnergyAloft
                
                tmpcLayerCounts[tmpcLayers] += 1
    
                wbLayerCounts[wbLayers] += 1
                                        
                tmpcWarmLayerAloftCount = tmpcLayers / 2
                
                if tmpcWarmLayerAloftCount > maxtmpcWarmLayerAloftCount:
                    maxtmpcWarmLayerAloftCount = tmpcWarmLayerAloftCount
    
                wbWarmLayerAloftCount = wbLayers / 2
                
                if wbWarmLayerAloftCount > maxWbWarmLayerAloftCount:
                    maxWbWarmLayerAloftCount = wbWarmLayerAloftCount
                
                if p.tmpc[0] > 0:
                    minEnergies.append(sf.calcMinEnergy(p.tmpc[0]))
                else:
                    minEnergies.append(0)
                
        for i, v in enumerate(sfcTemps):
            outfh.write("{} {} {} {} {} {} {}\n".format(
                        sfcTemps[i], sfcTempEnergies[i], minEnergies[i], freezingLevels[i],
                        sfcWbTemps[i], sfcWbEnergies[i], wbFreezingLevels[i]))

        outfh.close()
        
plotting.ptypes_on_sfct_vs_fzlevel(ptypes, sfcTemps, freezingLevels, "T")

#tFilt = sTemps[sTemps > 32]
#eFilt = sEnergies[sTemps > 32]

#tFilt = tFilt[tFilt < 40]
#eFilt = eFilt[tFilt < 40]

plotting.sfct_vs_energy(sfcTemps, sfcTempEnergies, minEnergies,"Temperature")

plotting.sfct_vs_energy(sfcWbTemps, sfcWbEnergies, minEnergies,"Wet Bulb Temperature")

plotting.maxt_vs_energy(aloftMaxTs, posTmpcEnergies, "MaxT")

plotting.maxt_vs_energy(aloftMaxWbs, posWbEnergies, "MaxWb")

plotting.maxt_vs_depth(aloftMaxTs, warmTAloftDepths, "Temperature")

plotting.maxt_vs_depth(aloftMaxWbs, warmTwAloftDepths, "Wet Bulb Temperature")

plotting.maxtXdepth_vs_energies(aloftMaxTs, warmTAloftDepths, posTmpcEnergies, "T")

# plotting.energy_vs_depth(posTmpcEnergies, warmTAloftDepths, "Temperature")

# plotting.energy_vs_depth(posWbEnergies, warmTwAloftDepths, "Wet Bulb")

plotting.depth_vs_energy(warmTAloftDepths, posTmpcEnergies, "Temperature")

plotting.depth_vs_energy(warmTwAloftDepths, posWbEnergies, "Wet Bulb")

plotting.t_vs_refreeze_depth(minRefreezeTs, refreezeTDepths, refreezeTmpcEnergies, "T")

plotting.t_vs_refreeze_depth(minRefreezeTws, refreezeTwDepths, refreezeWbEnergies, "Tw")

plotting.maxt_vs_mint(aloftMaxTs, minRefreezeTs, "T")

plotting.maxt_vs_mint(aloftMaxWbs, minRefreezeTws, "Tw")

plotting.maxenergy_vs_minenergy(posTmpcEnergies, refreezeTmpcEnergies, "T")

plotting.maxenergy_vs_minenergy(posWbEnergies, refreezeWbEnergies, "Tw")

plotting.sfct_vs_energyAloft(sfcTemps, posTmpcEnergies, "T")

print profileCount, "profiles"

print "Overall temp and wb layer counts:", tmpcLayerCounts, wbLayerCounts

# print "Most levels", maxLevels
# print "Least levels", minLevels

# print "Average number of levels per profile:", (profileLevels / profileCount)
print "Max tmpc warm layer count", maxtmpcWarmLayerAloftCount
print "Max wb warm layer count", maxWbWarmLayerAloftCount
