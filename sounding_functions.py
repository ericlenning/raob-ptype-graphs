# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 12:02:26 2017

@author: eric.lenning
"""
# import sharppy.sharptab.profile
import numpy.ma as ma
import numpy as np
import math

clamp = lambda n, minn, maxn: max(min(maxn, n), minn)

def count_soundings(fp):
    num_soundings = 0
    for line in fp:
        if line.lstrip().startswith('254'):
            num_soundings+=1
    return num_soundings

def get_sfc_tmpc(prof):
    tmpc = prof.tmpc
    
    return tmpc[0]

def get_sfc_wetbulb(prof):
    wetbulb = prof.wetbulb
    
    return wetbulb[0]

def get_freezing_level(hghts, temps):
    i = 0
    while i < len(temps):
        if (temps[i] <= 0):
            return hghts[i]-hghts[0]
        i += 1
"""
def get_positive_energy(hghts, temps):
    sfcTemp = temps[0]
    inWarmLayer = False
        
    i = 1
    
    while i < len(temps):
"""
def get_maxt_aloft(temps):

    firstNegative = next(x[0] for x in enumerate(temps) if x[1] < 0)
    
    result = max(temps[firstNegative:])
    
    # Return -9999 if entire sounding is below freezing.
    
    if result > 0:
        return result
    else:
        return 0

def get_maxt_above900(temps, pres):

    firstOver900 = next(x[0] for x in enumerate(pres) if x[1] <= 900)
    
    result = max(temps[firstOver900:])
    
    # Return -9999 if entire sounding is below freezing.
    
    if result > 0:
        return result
    else:
        return 0

def get_maxt_above2kft(temps, hght):

    hght = hght - hght[0]
    
    # 2000 ft is 609.6 meters
    
    firstOver2kft = next(x[0] for x in enumerate(hght) if x[1] >= 609.6)
    
    result = max(temps[firstOver2kft:])
    
    # Return -9999 if entire sounding is below freezing.
    
    if result > 0:
        return result
    else:
        return 0
    
def get_posenergy_above900(hght, temps, pres):

    firstOver900 = next(x[0] for x in enumerate(pres) if x[1] <= 900)
    
    energiesOver900 = calculate_energy(hght[firstOver900:], temps[firstOver900:], pres[firstOver900:], False)
    
    result = sum(x for x in energiesOver900 if x > 0)
    
    # Return -9999 if entire sounding is below freezing.
    
    if result > 0:
        return result
    else:
        return 0    

def get_posenergy_above2kft(hght, temps, pres):

    firstOver2kft = next(x[0] for x in enumerate(hght) if x[1] >= 609.6)
    
    energiesOver2kft = calculate_energy(hght[firstOver2kft:], temps[firstOver2kft:], pres[firstOver2kft:], False)
    
    result = sum(x for x in energiesOver2kft if x > 0)
    
    # Return -9999 if entire sounding is below freezing.
    
    if result > 0:
        return result
    else:
        return 0    
    
def get_negenergy_below850(hght, temps, pres):
    firstOver850 = next(x[0] for x in enumerate(pres) if x[1] <= 850)
    
    negEnergyBelow850 = calculate_energy(hght[:firstOver850], temps[:firstOver850], pres[:firstOver850:], True)
    
    result = sum(x for x in negEnergyBelow850 if x < 0)
    
    # Return -9999 if entire sounding is below freezing.
    
    if result < 0:
        return result
    else:
        return 0    

# This doesn't assume/require a surface based layer.  If not a surface based
# layer then this gets ignored by the traditional top-down method, I think.
        
def get_min_refreeze_t(temps):
    
    # print "gmrt: temps", temps
    
    # Get index of first negative temperature in sounding.  Could be at surface.
    
    firstNegative = next(x[0] for x in enumerate(temps) if x[1] < 0)
    
    # Starting from the first negative, get index of first zero above this.
    # With no warm layer aloft, there won't be a zero above this level.
    
    refreezeTop = firstNegative + next(x[0] for x in enumerate(temps[firstNegative:]) if x[1] == 0)
    
    #print "gmrt: firstNeg refreezeTop", firstNegative, refreezeTop
    
    result = min(temps[firstNegative:refreezeTop])
    
    #print "gmrt: result", result
    
    return result

"""
RefreezeLayerDepth (m) -- This could be calculated two ways.  
    -- In the traditional top-down technique, this is required to be a surface
       based layer, with a SfcT < 0C.  Otherwise refreezing is not considered.
    -- It also could be calculated using the first subzero layer in the sounding.
    
    Since this is only used by the traditional top-down technique, we will
    do the calculation according to their method, setting it to zero if not at
    the surface.
"""
def getRefreezeLayerDepthMeters(temps, hghts):
    
    
    if temps[0] >= 0:
        return 0
    
    else:
        """
        # Get index of first negative temperature in sounding.  Could be at surface.
        
        firstNegative = next(x[0] for x in enumerate(temps) if x[1] < 0)
            
        # Only look at portion of sounding above first negative value.
        # What if we never cross the zero line (i.e. rest of sounding is freezing)?
        
        zeroAbove = next(x[0] for x in enumerate(temps[firstNegative:]) if x[1] == 0)
        
        refreezeTop = firstNegative + zeroAbove
        
        # We want level where sounding crosses zero, or the surface.
        
        if firstNegative > 0:
            firstNegative -= 1    
            
        # print "grld temps zeroAbove firstNeg refreezeTop:", \
           # temps, zeroAbove, firstNegative, refreezeTop
        
        thickness = hghts[refreezeTop] - hghts[firstNegative]
        
        # print "grld hghts thickness:", hghts, thickness
        """
        refreezeTop = next(x[0] for x in enumerate(temps) if x[1] == 0)
        
        thickness = hghts[refreezeTop] - hghts[0]
        
        # print "grld hghts thickness:", hghts, thickness        
        
        return thickness
    
def getWarmLayerDepth(temps, hghts):
    
    # Only look at portion of sounding above first negative value
    
    if (temps[0] < 0):
        firstNegative = 0
    else:
        firstNegative = next(x[0] for x in enumerate(temps) if x[1] < 0)
    
    # Find the warmest temperature above the first negative value
    
    tmax = max(temps[firstNegative:])
    
    if tmax > 0:
        
        # Shorten the heights list to match the shortened temps list
        
        tmpHghts = hghts[firstNegative:]
        
        # Figure out the index of the first zero value above and below the maxT
        
        bottomIndex = 0
        topIndex = 0
        tMaxFound = False
        
        for i,t in enumerate(temps[firstNegative:]):
            if t == 0:
                if not tMaxFound:
                    bottomIndex = i
                else:
                    topIndex = i
                    depth = tmpHghts[topIndex] - tmpHghts[bottomIndex]
                    return depth
            elif t == tmax:
                tMaxFound = True

    else:
        return -9998    
    
# There is a better/faster/easier way to do this.  Don't use this method.

def countWarmLayersAloft(temperatures):
    
    temperatures[temperatures == -9999] = ma.masked
    
    temps = temperatures.compressed()
    
    warmLayersAloft = 0    
    
    notInWarmLayer = True
    
    if (temps[0] >= 0) and (temps[1] >= 0):
        notInWarmLayer = False
        
    # reaching 0.0 doesn't start a new layer, only passing it does
    
    for t in temps[2:]:
        #print t
        if notInWarmLayer:
            if t > 0: # new warm layer aloft
                notInWarmLayer = False
                warmLayersAloft += 1
        else:
            if t < 0:
                notInWarmLayer = True        
    
    return warmLayersAloft

# A very basic calculation of the minimum possible positive energy in a
# surface layer based on a dry adiabatic lapse rate from the starting surface
# temperature (C).
    
def calcMinEnergy(sfcTemp):
    dalr = 9.8 # dry adiabatic lapse rate
    
    fzLevelMeters = 1000 * sfcTemp / dalr
    
    avgLayerTemp = 0.5*sfcTemp
    
    minEnergy = dalr * fzLevelMeters * avgLayerTemp / 273.15
    
    # minEnergy = dalr * 1000 * (sfcTemp / dalr) * 0.5 * sfcTemp / 273.15
    
    return minEnergy

"""        
Energy calculation for a layer (T in K and hghts in m) is:
(9.8*(thickness)*(avg temp))-273.15) / 273.15

Calculate and return the energies in each layer of a sounding.  The layer
above the last zero will not be included unless truncated is True meaning we
are only examining up to a certain pressure or height level.

"""
 
def calculate_energy(h, T, p, truncated):

    lastT = len(T) - 1
    
    #print "len(p) and len(T/Tw) are:",len(p),len(T)
    #print "p and T/Tw are:",p,T
    
    energies = []

    # Keep track of energy in a single layer
    
    layersum = 0
 
    layercount = 0

    # Keep track of coldest or warmest temperature in each layer.
    
    # layerExtreme = []

    # Our surface temperature is by default the max or min T so far in its layer.
    
    # layerExtreme.append([p[0], T[0]])

    # Keep track if we're presently examining a layer.
    
    examiningLayer = False
    
    # Each time we CROSS zero, start a new layer.
    # Our data will have interpolated zero levels in T and Tw.
    # So upon encountering a zero, start a new layer.
    # Above our last zero we won't store any data.
    # Sometimes a zero occurs in multiple adjacent layers.

    for i, v in enumerate(h):
        
        if i > 0:  # can't compute thickness of layer below the surface
            
            examiningLayer = True
            
            thickness = h[i] - h[i-1]
            
            avgtemp = (T[i] + T[i-1])/2.0

            layerenergy = 9.8*thickness*avgtemp/273.15

            layersum += layerenergy

            #print "Checking T[i] of {} with i of {} against {}".format(T[i],i,layerExtreme[layercount][1])

            #if layerenergy < 0:  # In a cold layer
            #    if T[i] < layerExtreme[layercount][1]:
            #        #print "In COLD layer {} with topT {} colder than {}".format(layercount,T[i],layerExtreme[layercount][1])
            #        layerExtreme[layercount] = [p[i],T[i]]
            #else:         # In a warm layer
            #    if T[i] > layerExtreme[layercount][1]:
            #        #print "In WARM layer {} with topT {} warmer than {}".format(layercount,T[i],layerExtreme[layercount][1])
            #        layerExtreme[layercount] = [p[i],T[i]]
                    
            # If we've reached a zero level and the next temperature is of
            # opposite sign of the energy in this layer, we are done with this
            # layer. And we could have zeros in adjacent layers.

            if T[i] == 0:
                
                if i < lastT and T[i+1] and (layersum * T[i+1] < 0):
                    
                    energies.append(layersum)
    
                    #print "Finished layer {} with energy of {}".format(layercount,layersum)
    
                    layersum = 0
     
                    layercount = layercount + 1
    
                    # We are done examining this layer and haven't started the
                    # next one yet.
                    
                    examiningLayer = False
                    
                    #layerExtreme.append([p[i],T[i]])
                
    # In case we never got to a zero layer

    if (examiningLayer & truncated):
        energies.append(layersum)

    #energies.append( (9.8*thickness*avgtemp-273.15)/273.15 )

            
        #plt.plot(p,earray)
        #plt.fill_betweenx(p, 0, earray, (earray >= 0), color='red', alpha=.25)
        #plt.fill_betweenx(p, 0, earray, (earray <= 0), color='blue',  alpha=.25)
    
    #esum = sum([i for i in energy if i > 0])

    #print energies,esum
    """
    finalLayerCount = len(energies)
    
    if finalLayerCount > 0:
        print finalLayerCount,"energy layers in",energies,"and maxes/mins are:",layerExtreme[:finalLayerCount]
    else:
        print "Zero layers -- all below freezing."
    """
    
    # So we have the list called layerExtreme and the list called energies.
    # Each entry contains info for one positive or negative energy layer

    #for i, m in enumerate(energies):
    #    #print i
    #    maxp = layerExtreme[i][0]
    #    maxt = layerExtreme[i][1]

    return energies

"""
We only analyze soundings with a top temperature (at highest height) that is below zero.

If the bottom layer has negative energy, zero should be crossed an even number of times.
    Ptype options are PL, SN, ZR
If the bottom layer has positive energy, zero should be crossed an odd number of times.
    Ptype options are PL, SN, RA
 
If the highest pressure level (ptop) is below 700 we don't try to forecast sleet.

We return a 4-element tuple with each place having one of:
 0 (low/no chance) -- less than 25%
 1 (chance)        -- 25% to 74%
 2 (likely)        -- 75% and higher
 
The four elements represent (RA, ZR, IP, SN)
"""
 
def calc_precip_energy(energylayers, ptop, sfct, sfctw): # ,fig,offset):

    probRA = 0
    probSN = 100  # Assume snow.
    probZR = 0
    probPL = 0

    layers = np.array(energylayers)

#    print "Calculating energy with these layers:", layers

    soundingType = "All below freezing"
    
    result = (0, 0, 0, 2)

    # no layers; all below freezing
    
    if len(layers) == 0:
        return result
    
    if (all(i <=0 for i in layers)): # Ptype options are SN or FZRA/FZDZ but say SN for now.

        return result

        """
        if ptop <= 600:
            probPL = 0
            probZR = 0
        """
        
    else:  # some layers above freezing, though perhaps just sfc

        sfcEnergy = energylayers[0]

        # Maybe liquid or sleet is being introduced into sfc layer.
        # For simplicity compare combination of all pos layers and all neg layers above sfc.
        
        negativeEnergy = 0

        negativeEnergy = -1 * sum([e for e in energylayers if e < 0])      # All cold air in sounding

        positiveEnergyAloft = 0

        positiveEnergyAloft = sum([e for e in energylayers[1:] if e > 0])  # Any warm air above surface

        if ptop >= 700:
            print "Sleet not forecast: ptop of {} is too low.".format(ptop)
            probPL = 0

        # The following sleet/liquid relationships were derived by Pierre Bourgouin and revised by Kevin Birk.
        # Original line to be exceeded for ProbPL to equal 100 was:
        #    y = 65 + 0.66x
        # Revised line from Birk's work is flatter and has a larger y-intercept:
        #    y = 105 + 0.128x
        # The flatter slope means that overall the amount of cold air required to refreeze partially 
        # or fully melted hydrometeors does not depend quite as much on the degree of warming aloft.
        # This change makes it easier to get freezing rain when there is a relatively small amount of warm air
        # aloft (<75 J/kg), but easier to get sleet when there is relatively more warm air aloft.

        # Former equations to see if there is enough cold air to produce sleet below warm air aloft.
        #if (negativeEnergy >= 65+(0.66*positiveEnergyAloft)):
            #probPL = 100
        #elif (negativeEnergy >= 56+(0.66*positiveEnergyAloft)):
            #probPL = 75
        #elif (negativeEnergy >= 46+(0.66*positiveEnergyAloft)):
            #probPL = 50
        #elif (negativeEnergy >= 30+(0.66*positiveEnergyAloft)):
            #probPL = 24

        # New equations to see if there is enough cold air to produce sleet below warm air aloft.
        """
        if (negativeEnergy >= 105+(0.128*positiveEnergyAloft)):
            probPL = 100
        elif (negativeEnergy >= 85+(0.66*positiveEnergyAloft)):
            probPL = 75
        elif (negativeEnergy >= 70+(0.66*positiveEnergyAloft)):
            probPL = 50
        elif (negativeEnergy >= 58+(0.66*positiveEnergyAloft)):
            probPL = 24
        """
        
        # Just use three categories, and recall probPL defaults to 0
        
        if (negativeEnergy >= 85+(0.128*positiveEnergyAloft)):
            probPL = 2
        elif (negativeEnergy >= 58+(0.128*positiveEnergyAloft)):
            probPL = 1

        # probRA = 0

        # Former equations to see if there is enough cold air under warm air to reduce possibility of liquid to surface.
        #if (negativeEnergy > 85+(0.66*positiveEnergyAloft)):
            #probRA = 20   # Always allow for a small chance of liquid.
        #elif (negativeEnergy <= 66+(0.66*positiveEnergyAloft)):
            #probRA = 100  
        #elif (negativeEnergy <= 75+(0.66*positiveEnergyAloft)):
            #probRA = 75   
        #elif (negativeEnergy <= 85+(0.66*positiveEnergyAloft)):
            #probRA = 50   

        # New equations to see if there is enough cold air under warm air to reduce possibility of liquid to surface.
        """
        if (negativeEnergy > 105+(0.128*positiveEnergyAloft)):
            probRA = 20   # Always allow for a small chance of liquid.
        elif (negativeEnergy <= 58+(0.128*positiveEnergyAloft)):
            probRA = 100  
        elif (negativeEnergy <= 85+(0.128*positiveEnergyAloft)):
            probRA = 75   
        else:  # if (negativeEnergy <= 105+(0.128*positiveEnergyAloft)):
            probRA = 50   
        """
        
        # Just use three categories, and recall probRA defaults to 0
        
        if (negativeEnergy > 105+(0.128*positiveEnergyAloft)):
            probRA = 1   # Always allow for a small chance of liquid.
        elif (negativeEnergy <= 58+(0.128*positiveEnergyAloft)):
            probRA = 2  
        elif (negativeEnergy <= 85+(0.128*positiveEnergyAloft)):
            probRA = 2   
        else:  # if (negativeEnergy <= 105+(0.128*positiveEnergyAloft)):
            probRA = 1               

        # Sanity check for not enough warm air aloft to melt snow at all.
        # Note: Birk's GFE tools use 0.5 J/kg for the sanity check.
        # Assuming no liquid with < 5 J/kg makes it consistent with calcs in sfc layer.
        # Also, with < 5 J/kg aloft there could be sleet but the impact likely is minimal.

        if ((probPL > 0) & (positiveEnergyAloft > 0) & (positiveEnergyAloft < 5)):    #   OR  0.5)):
            #probPL = 0
            probPL = 1 # int(round(probPL * positiveEnergyAloft/5.0,-1))

        # Sanity check for not enough warm air aloft to melt snow much.

        if positiveEnergyAloft < 5:
            probRA = 0
        elif positiveEnergyAloft < 10:
            probRA = 1 # min(int(round(16.5*positiveEnergyAloft - 74.25,-1)),probRA)
        else:
            probRA = probRA

        # Now just look at melting of snow by warm air aloft.

        if positiveEnergyAloft > 20:
            probSN = 0
        elif positiveEnergyAloft > 13:
            probSN = 1 # min(int(round(-14.286*positiveEnergyAloft + 292.86,-1)),100)
        else:
            probSN = 2
            
        if sfcEnergy < 0: # Ptype options are PL, SN, ZR
            soundingType = "Some above zero, but not sfc"
            probZR = probRA
            probRA = 0
        else:             # Ptype options are PL, SN, RA
            soundingType = "Some above zero, including sfc"
            probZR = 0

            # The basic probRA/probSN assume snow is being introduced into sfc layer.

            if sfcEnergy < 5:
                probRA = probRA # 0
            elif sfcEnergy < 10:
                probRA = max(1, probRA) # max(int(round(16.5*sfcEnergy - 74.25,-1)),probRA)
            else:
                probRA = 2

            # Look at melting of snow and ice in warm sfc layer.

            if sfcEnergy > 20:
                probSN = 0
                probPL = 0
            elif sfcEnergy > 13:
                probSN = min(1, probSN) # min(int(round(-14.286*sfcEnergy + 292.86,-1)),probSN)
                probPL = min(probSN, probPL)
            else:
                probSN = min(2, probSN)
                
        return (probRA, probZR, probPL, probSN)
            
    #ptypeString = "RA: {} SN: {} ZR: {} PL: {}".format(probRA,probSN,probZR,probPL)

    #ax.annotate(ptypeString,xy=(-60,125+offset_mb),xytext=(-60,125+offset_mb),size='x-small')
    # fig.text(.55,.8-offset,ptypeString,size='x-small',backgroundcolor='w')

    # fig.text(.55,.8-offset-0.02,soundingType,size='x-small',backgroundcolor='w')

    #print energylayers,"means",ptypeString

"""
Determine precip type using the traditional top-down method.  Assume for simplicity
we have ice nucleation, since soundings were filtered for saturation in the DGZ.

RefreezeDepth is in meters.

Also, MaxTAloft may be MaxTwAloft, and MinRefreezeT may be MinRefreezeTw.

We return a 4-element tuple with each place having one of:
 0 (low/no chance) -- less than 25%
 1 (chance)        -- 25% to 74%
 2 (likely)        -- 75% and higher
"""

def calc_precip_traditional(MaxTAloft, MinRefreezeT, RefreezeDepth, 
                            SfcT, SfcTw, temps, hghts, tStr):
    
    # Arbitrary cutoff temperature for no mention of snow?  37F = 2.78C, 38F = 3.33C
    # Arbitrary cutoff temperature for no mention of liquid?  35F = 5J = 1.67C
    sThreshold = 1.67
    rThreshold = 3.33
    
    # Starting assumptions for when MaxTAloft < 1    
    probRA = 0
    probZR = 0
    probPL = 0
    probSN = 100  # Assume snow.
    
    """
    Rules for ProbRefreezeSleet, which you'll need later:
        Cold layer must be surface based and >= 2500 ft (762 m) deep.
        Only considered when MaxTAloft >= 3C.
    """
    ProbRefreezeSleet = 0

    if RefreezeDepth >= 762:
        ProbRefreezeSleet = 2.3036 * MinRefreezeT * MinRefreezeT + \
                            5.979 * MinRefreezeT + 0.8295
        ProbRefreezeSleet = clamp(ProbRefreezeSleet, 0, 100)

    """
    ProbLiquid -- before we know ProbRA or ProbZR
    Now figure out the chance of getting liquid to the surface.
    This increases with MaxTAloft >= 2.5 C or with SfcT > 0.
    Liquid to the surface can end up as RA, ZR, or a mix.
    """
    ProbLiquid = 0

    # ProbLiquid chances start at 2.5C aloft and are 100% at 4C aloft.
    # Otherwise they remain zero.
    
    if MaxTAloft >= 4:
        ProbLiquid = 100
    elif MaxTAloft > 2.4:
        ProbLiquid = 6.4286 * MaxTAloft * MaxTAloft - 1.7143 * MaxTAloft + 1.8571

    # Lower the chance of liquid if you have a refreezing layer.
    # Not sure why we subtract instead of multiplying.
    
    ProbLiquid = ProbLiquid - ProbRefreezeSleet
    ProbLiquid = clamp(ProbLiquid, 0, 100)

    # Increase chance of liquid to 100% if SfcT > 0 and MaxTAloft > 0 (or 1?)
    # This is supposed to address isothermal layers.
    
    if ((SfcT > 0) & (MaxTAloft > 0)):
        ProbLiquid = 100

    # Where T > MaxTwAloft and MaxTwAloft < 0.5, use the SfcT snow threshold.
    # This is basically a sanity check?
    
    if (SfcT > MaxTAloft) & (MaxTAloft <= 0):
        if SfcT > sThreshold:
            ProbLiquid = 100
        else:
            ProbLiquid = 0

    # Separate out the true rain/drizzle from the freezing rain/drizzle.
    #PotTrueLiquid = 0
    #PotFreezingLiquid = 0

    if SfcT > 0:
        probRA = ProbLiquid
    else:
        probZR = ProbLiquid

    if MaxTAloft > 0:
        #Need MaxTAloft > 0, T between 32 and 35, and Tw <= 0 to consider assigning freezing
        wetBulbFreezing = ((MaxTAloft > 0) & (SfcT > 0) & (SfcT < 2.22) & (SfcTw <= 0))
        zrTrue = ((ProbLiquid > 0) & wetBulbFreezing)

        """
        if MaxTwAloft is not None:   
            #Need MaxTwAloft > 0, T between 32 and 35, and Tw <= 0 to consider assigning freezing
            wetBulbFreezingMask = ((MaxTwAloft > 0) & (T > 32) & (T < 36) & (Tw <= 0))
            Mask = ((PotLiquid > 0) & wetBulbFreezingMask)         

            # In areas of special SfcT range, set ProbFreezing using equation 
            # and ProbLiquid as its inverse
            
            PotFreezingLiquid[Mask] = (100 * (5572412337462.57 * exp(-0.930254926892308 * T)))[Mask]
            PotTrueLiquid[Mask] = (100 - PotFreezingLiquid)[Mask]
    
            # In areas of special SfcT range, and areas where PotFreezing > PotLiquid,
            # set the PotFreezing to 100.  In those same areas...?
            
            greaterFreezingMask = (Mask & (PotFreezingLiquid > PotTrueLiquid))
            PotFreezingLiquid[greaterFreezingMask] = 100
            PotTrueLiquid[greaterFreezingMask] = ((100/PotFreezingLiquid) * PotTrueLiquid)[greaterFreezingMask]
                             
            greaterLiquidMask = (Mask & (PotTrueLiquid > PotFreezingLiquid))
            PotTrueLiquid[greaterLiquidMask] = 100
            PotFreezingLiquid[greaterLiquidMask] = ((100/PotTrueLiquid) * PotFreezingLiquid)[greaterLiquidMask]

        """

        # Here T seems to be in F, so do math.
        T = 32 + 1.8 * SfcT
        
        # There's a thin range where ZR is possible above 0 C at the sfc.       
        if zrTrue:
            probZR = (100 * (5572412337462.57 * math.exp(-0.930254926892308 * T)))
            probRA = (100 - probZR)

            # If the calculated chance of freezing liquid is larger than for 
            # pure liquid, then...
            
            if probZR > probRA:
                probZR = 100
                probRA = ((100/probZR) * probRA)                            
            
            # ...otherwise...
            
            if probRA > probZR:
                probRA = 100
                probZR = ((100/probRA) * probZR)

        """
        Now for sleet...
        """
        
        MaskZeroTwo = ((MaxTAloft > 0) & (MaxTAloft < 2))
        MaskTwoThree = ((MaxTAloft >= 2) & (MaxTAloft <= 3))
        MaskThreeFour = ((MaxTAloft > 3) & (MaxTAloft < 4))
        
        if MaskZeroTwo:
            probPL = (30.857 * MaxTAloft * MaxTAloft - 12.114 * MaxTAloft + 0.2286)
        elif MaskTwoThree:
            probPL = 100
        elif MaskThreeFour:
            probPL = (40 * MaxTAloft * MaxTAloft - 380 * MaxTAloft + 880)

        # ProbRefreezeSleet is only considered when MaxTAloft > 3
        # Ensure probPL is at least equal to the ProbRefreezeSleet
        
        if MaxTAloft > 3:
            if probPL < ProbRefreezeSleet:
                probPL = ProbRefreezeSleet
            
        if ((SfcT > 0) & (MaxTAloft >= 1.0)):
            probPL = 0
                
        """
        Now for snow...which by default is 100.
        """

        if (MaxTAloft >= 3):
            
            probSN = 0
            
        elif (MaxTAloft > 0):
            #PotSnow uses the formula PotSnow = 8.3333*(MaxTwAloft)^3 - 40 *(MaxTwAloft)^2 + 11.667*MaxTwAloft + 100
            #PotSnow uses the formula PotSnow = 10.889*(MaxTwAloft)^3 - 49.041 *(MaxTwAloft)^2 + 16.138*MaxTwAloft + 99.599

            probSN = (10.889*MaxTAloft*MaxTAloft*MaxTAloft - 49.041*MaxTAloft*MaxTAloft + 16.138*MaxTAloft + 99.599)

            #Sep 7, 2017 - commented out the following code for Mask1 which results in a smoother transition from 100 > 0.
            #Mask1 = ((T > 32) & (MaxTwAloft > 0) & (MaxTwAloft < 1.0) & (T < rthreshold))
            #PotSnow[Mask1] = 100
            
            if (SfcT > 0) & (MaxTAloft >= 1.0):
                probSN = 0
                
            # Where the hourly T > MaxTwAloft and MaxTwAloft < 1, use the forecaster entry surface T for all liquid.
            # First have to convert hourly T to C.
            #TinC = (T-32) * 5/9
            
            if ((SfcT > MaxTAloft) & (MaxTAloft < 1) & (SfcT < rThreshold)):
                #Mask1 = T < rthreshold
                probSN = 100
                
            #PotSnow = where(logical_and(greater(TinC, MaxTwAloft), less(MaxTwAloft, 1)), where(less(T,rthreshold),100.,PotSnow), PotSnow)
            if ((SfcT > MaxTAloft) & (MaxTAloft < 1) & (SfcT >= rThreshold)):
                #Mask2 = T >= rthreshold
                probSN = 0
            #PotSnow = where(logical_and(greater(TinC, MaxTwAloft), less(MaxTwAloft, 1)), where(greater_equal(T,rthreshold),0.,PotSnow), PotSnow)
           
        else:
            # Adjust PotSnow to 0 where surface temperatures are greater than the forecaster entry for all liquid
            if SfcT > rThreshold:
                probSN = 0

    probs = []
    
    for prob in probRA, probZR, probPL, probSN:            
        if prob >= 75:
            probs.append(2)
        elif prob >= 25:
            probs.append(1)
        else:
            probs.append(0)
            
    return tuple(probs)