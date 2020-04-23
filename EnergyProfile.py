# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 10:03:45 2017

@author: eric.lenning
"""
import numpy as np
import sounding_functions as sf

"""
A class containing an EnergyProfile meaning a profile represented in terms
of energy layers and related properties.  These properties can be based on
either Temperature or Wet Bulb, so the descriptions below are generic and do 
not assume one or the other.  In the code, 'temp' stands for either a T or Tw
profile.

These properties include:
    
    energyLayers -- list of each enery layer in the profile - could be empty
        
    energyLayerCount -- number of energy layers - could be Zero if whole
                        profile is subfreezing
    
    sfcTemp -- temperature (C) at surface
    
    freezingLevel -- height AGL (m) of the freezing level - could be Zero if whole
                     profile is subfreezing    
    
    ptype -- tuple representing chance of RA, ZR, IP, SN 
             with 0 meaning unlikely, 1 meaning chance, and 2 meaning likely.
             So (0, 1, 0, 2) means chance of ZR with SN likely. 
             Shouldn't be (0, 0, 0, 0).
    
    maxTempAloft -- maximum temperature (C) in a warm layer aloft - could be
                    Zero if entire profile is subZero
    
    posEnergyAloft -- total amount of positive energy in layers other than the
                      surface layer - could be Zero
    
    posEnergyTotal -- total amount of positive energy in all layers - could be
                      Zero
                      
    refreezeEnergy -- total amount of negative/refreezing energy in the profile -
                      could ironically be Zero if entire sounding is subfreezing
    
    sfcEnergy -- amount of energy (pos or neg) in the lowest layer - could be 
                 Zero if entire sounding is subfreezing
    
    minRefreezeTemp -- coldest temperature in a refreezing layer - could be
                       Zero if entire sounding is subfreezing
    
    warmAloftDepth -- depth (m) of the deepest elevated warm layer - could be
                      Zero if no elevated warm layers
    
    refreezeLayerDepthMeters -- depth (m) of the deepest refreezing layer - could be
                          Zero if no elevated warm layers
    

"""

class EnergyProfile(object):
    def __init__(self, **kwargs):
        
        ## set the missing variable
        
        ## get the data and turn some of them into arrays
        self.pres = np.array(kwargs.get('pres'), dtype=float)
        self.hght = np.array(kwargs.get('hght'), dtype=float)
        self.temp = np.array(kwargs.get('temp'), dtype=float)
        
        self.sfct = kwargs.get('sfct')
        self.sfctw = kwargs.get('sfctw')
        self.tstr =  kwargs.get('tstr')
#        self.dwpc = np.array(kwargs.get('dwpc'), dtype=float)
        
        self.freezingLevel = sf.get_freezing_level(self.hght, self.temp)

        # External
        # freezingLevels.append(freezingLevel)

        # External
        # sfcTemps.append(p.tmpc[0])
        

        # Calculate the magnitude of each energy layer in a sounding.
        # - A sounding with zero energy layers starts and stays below
        #   freezing.
        # - A sounding with an odd number of energy layers starts above
        #   freezing and crosses the freezing level an odd number of
        #   times. The even indices (0,2,4, etc.) are warm layers.
        # - A sounding with an even number of energy layers starts below
        #   freezing and crosses the freezing level an even number of
        #   times. The odd indices (1,3,5, etc.) are warm layers.
        
        self.energyLayers = sf.calculate_energy(self.hght, self.temp, self.pres, False)                
        
        # print self.location, self.date
          
        self.ptype = sf.calc_precip_energy(self.energyLayers, self.pres[-1], self.sfct, self.sfctw)
        
        # External
        # ptypes.append(sf.calc_precip(energyLayers, self.pres[-1]))
        
        self.sfcTemp = self.temp[0]
        
        self.maxTempAloft = 0                
        
        self.maxTempAbove900 = 0
        
        self.maxTempAbove2kft = 0
        
        self.posEnergyAloft = 0
        
        self.posEnergyAbove900 = 0
        
        self.posEnergyAbove2kft = 0
        
        self.negEnergyBelow850 = 0
        
        self.posEnergyTotal = 0
        
        self.sfcEnergy = 0
        
        self.minRefreezeTemp = 0
        
        self.refreezeEnergyLowest = 0
        
        self.refreezeEnergyAll = 0
        
        self.refreezeEnergySfc = 0
        
        self.warmAloftDepth = 0
        
        self.refreezeLayerDepthMeters = 0
        
        self.energyLayerCount = len(self.energyLayers)

        # These next 'Above' calculations apply to all soundings...sorta.
        
        self.maxTempAbove900 = sf.get_maxt_above900(self.temp, self.pres)
        
        self.maxTempAbove2kft = sf.get_maxt_above2kft(self.temp, self.hght)
        
        self.posEnergyAbove900 = sf.get_posenergy_above900(self.hght, self.temp, self.pres)
        
        self.posEnergyAbove2kft = sf.get_posenergy_above2kft(self.hght, self.temp, self.pres)

        self.negEnergyBelow850 = sf.get_negenergy_below850(self.hght, self.temp, self.pres)                
        
        self.posEnergyBelow2kft = 0
        
        # Do only if we have at least one energy layer
        
        if self.energyLayerCount > 0:
            
            self.sfcEnergy = self.energyLayers[0]
            
            # External
            # sfcTempEnergies.append(energyLayers[0])
            
            # Do only if we have at least two energy layers, meaning we have
            # at least one warm layer aloft and at least one refreezing layer
            # below that.
            
            if self.energyLayerCount > 1:
                
                self.maxTempAloft = sf.get_maxt_aloft(self.temp)
                
                #self.negEnergyBelow850 = sf.get_negenergy_below850(self.hght, self.temp, self.pres)
                
                self.minRefreezeTemp = sf.get_min_refreeze_t(self.temp)
                
                # Don't count a surface-based warm layer in the posEnergyAloft
                # calculation, but count all subfreezing layers in refreezeEnergy.
                
                if self.energyLayerCount % 2 == 1:
                    self.posEnergyAloft = sum(self.energyLayers[2::2])
                    self.posEnergyTotal = sum(self.energyLayers[0::2])
                    self.refreezeEnergyLowest = self.energyLayers[1]
                    self.refreezeEnergyAll = sum(self.energyLayers[1::2])
                    self.refreezeEnergySfc = 0
                    
                else:
                    self.posEnergyAloft = sum(self.energyLayers[1::2])
                    self.posEnergyTotal = self.posEnergyAloft
                    self.refreezeEnergyLowest = self.energyLayers[0]
                    self.refreezeEnergyAll = sum(self.energyLayers[0::2])
                    self.refreezeEnergySfc = self.refreezeEnergyLowest
                
                self.posEnergyBelow2kft = self.posEnergyTotal - self.posEnergyAbove2kft
                
                self.warmAloftDepth = sf.getWarmLayerDepth(self.temp, self.hght)

                self.refreezeLayerDepthMeters = sf.getRefreezeLayerDepthMeters(self.temp, self.hght)

        self.ptypet = sf.calc_precip_traditional(
                self.maxTempAbove2kft, self.minRefreezeTemp, 
                self.refreezeLayerDepthMeters, self.sfct, self.sfctw, 
                self.temp, self.hght, self.tstr
                )
        
        # Do this quick check to make sure all is consistent.
        
        if self.maxTempAloft > 0:
            
            # External
            # aloftMaxTs.append(self.maxTempAloft)
            
            # External
            # posTmpcEnergies.append(self.posEnergyAloft)
            
            # External
            # warmAloftDepths.append(self.warmAloftDepth)
            
            if self.posEnergyAloft < 0:
                print "Neg T Energy:", self.energyLayerCount, self.maxTempAloft, self.posEnergyAloft, self.tmpc
        
            # External
            # minRefreezeTemps.append(self.minRefreezeTemp)            
            # refreezeTmpcEnergies.append(self.refreezeEnergy)            
            # refreezeTDepths.append(self.refreezeLayerDepthMeters)           
            
        #print "T", self.tmpc, self.maxTempAloft, self.posEnergyAloft

        #print "Tw", self.wetbulb, maxWbAloft, wbPosEnergyAloft
        
        # External
        # tmpcLayerCounts[self.energyLayerCount] += 1
        # wbLayerCounts[wbLayers] += 1                                
        # tmpcWarmLayerAloftCount = self.energyLayerCount / 2        
        # if tmpcWarmLayerAloftCount > maxtmpcWarmLayerAloftCount:
        #    maxtmpcWarmLayerAloftCount = tmpcWarmLayerAloftCount
        # wbWarmLayerAloftCount = wbLayers / 2        
        #if wbWarmLayerAloftCount > maxWbWarmLayerAloftCount:
        #    maxWbWarmLayerAloftCount = wbWarmLayerAloftCount