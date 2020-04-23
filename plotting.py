# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 10:54:53 2017

@author: eric.lenning
"""
#import heapq
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LinearSegmentedColormap
from PtypeMarker import PtypeMarker
import os 
import itertools

widthInches = 10
heightInches = 8
dpi = 80
sounding_width = 800
sounding_height = 640    

def get_url(yyyy, mm, dd, hh, stnid):
    snd_base = "http://weather.uwyo.edu/cgi-bin/sounding?region=naconf&" + \
               "TYPE=GIF%3ASKEWT&"
    snd_year = "YEAR={yyyy}&".format(yyyy=yyyy)
    snd_mm = "MONTH={mm:02d}&".format(mm=int(mm))
    #snd_ddhh = "FROM={dd:02d}{hh:02d}&TO={dd:02d}{hh:02d}&".format(dd=dd,hh=hh)
    snd_ddhh = "FROM={dd}{hh:02d}&TO={dd:02d}{hh:02d}&".format(dd=int(dd),hh=int(hh))
    snd_stnid = "STNM={stnid}".format(stnid=stnid)
               
    snd_url = snd_base + snd_year + snd_mm + snd_ddhh + snd_stnid
    
    return snd_url

def get_texturl(yyyy, mm, dd, hh, stnid):
    snd_base = "http://weather.uwyo.edu/cgi-bin/sounding?region=naconf&" + \
               "TYPE=TEXT%3ALIST&"
    snd_year = "YEAR={yyyy}&".format(yyyy=yyyy)
    snd_mm = "MONTH={mm:02d}&".format(mm=int(mm))
    #snd_ddhh = "FROM={dd:02d}{hh:02d}&TO={dd:02d}{hh:02d}&".format(dd=dd,hh=hh)
    snd_ddhh = "FROM={dd}{hh:02d}&TO={dd:02d}{hh:02d}&".format(dd=int(dd),hh=int(hh))
    snd_stnid = "STNM={stnid}".format(stnid=stnid)
               
    snd_url = snd_base + snd_year + snd_mm + snd_ddhh + snd_stnid
    
    return snd_url

# https://vortex.plymouth.edu/cgi-bin/sfc/gen-statlog-a.cgi?ident=Kdvn&pl=rawspec&yy=15&mm=12&dd=28&pg=print      
# https://vortex.plymouth.edu/cgi-bin/sfc/gen-statlog-a.cgi?ident=Kdvn&pl=rawspec&yy=15&mm=12&dd=28&pg=web    
def get_metarLink(yyyy, mm, dd, hh, stnid):
    mtr_base = "https://vortex.plymouth.edu/cgi-bin/sfc/gen-statlog-a.cgi?"
    mtr_stnid = "ident={stnid}&".format(stnid=stnid)
    mtr_pl = "pl=rawspec&"
    mtr_year = "yy={yyyy}&".format(yyyy=yyyy[2:])
    mtr_mm = "mm={mm:02d}&".format(mm=int(mm))
    #mtr_ddhh = "FROM={dd:02d}{hh:02d}&TO={dd:02d}{hh:02d}&".format(dd=dd,hh=hh)
    mtr_dd = "dd={dd:02d}&".format(dd=int(dd),hh=int(hh))
    mtr_pg = "pg=print"
               
    mtr_url = mtr_base +  mtr_stnid + mtr_pl + mtr_year + mtr_mm + mtr_dd + mtr_pg
    
    return mtr_url
    
def createWebpage(xcoords, ycoords, img_width, img_height, figname, title,
                  stns, dates, ptypes=None, ptypest=None):
    
    urls = []
    mtr_urls = []
    text_urls = []
    
    if ptypes is None:
        ptypes = [''] * len(xcoords)

    if ptypest is None:
        ptypest = [''] * len(xcoords)
        
    # The date fields will look like:  mmddyy/hh00 or 120114/0000
    for s, d in zip(stns, dates):
        mm = d.strftime("%m")
        dd = d.strftime("%d")
        yyyy = d.strftime("%Y")
        hh = d.strftime("%H")
        
        urls.append(get_url(yyyy, mm, dd, hh, s[0])) 
        
        mtr_urls.append(get_metarLink(yyyy, mm, dd, hh, s[1]))
        
        text_urls.append(get_texturl(yyyy, mm, dd, hh, s[0]))
    """
    div {
            display: inline-block;
            width: 1000px;
            overflow: hidden;
            border-style: solid;
            border-color: black;
        }        
    """
    f = open(figname + ".html", "w")
    f.write('''<HTML>
    <HEAD>
     <TITLE>%s</TITLE>
    </HEAD>
    <BODY>
    <STYLE>
    
    iframe#sounding {
        width: 900px;
        height: 700px;
    }
    iframe#metar {
        width: 600px;
        height: 200px;
    }
    iframe#soundingtext {
        width: 900px;
        height: 200px;
    }
    </STYLE>
    <SCRIPT>
    function mouseover(name) {
      var cid = document.getElementById("cid");
      cid.innerHTML = name;
    }
    function show(filename, mtrlink, textlink) {
      var sounding = document.getElementById("sounding");
      sounding.src = filename;
      //sounding.innerHTML='<object type="text/html" data="'+filename+'"></object>';
      //alert(filename);
      //sounding.width = "800px";

      var metar = document.getElementById("metar");
      metar.src = mtrlink;
      
      var text = document.getElementById("soundingtext");
      soundingtext.src = textlink;
    }
    
    </SCRIPT>    
    <table>
    <tr>
    <td>
    Mouse is over: <SPAN id="cid"></SPAN><BR>
    Pick a point to see the depiction<BR>
    <IMG SRC="%s" ismap usemap="#points"
      WIDTH="%d" HEIGHT="%d" style="vertical-align: text-top; float: left;">
    </td>
    <td rowspan="1" valign="top">
    <iframe ID="sounding">Sounding Plot</iframe>
    </td>
    </tr>
    <tr>
    <td align="center">
    <iframe ID="metar">METAR Plot</iframe>
    </td>
    <td>
    <iframe ID="soundingtext">Sounding Text</iframe>
    </td>
    </tr>
    </table>
    
    <MAP name="points">
    ''' % (title, os.path.basename(figname) + ".png", img_width, img_height))
    
    # Figure height in pixels.
    # HTML image coordinates have y=0 on the top.  Matplotlib
    # has y=0 on the bottom.  We'll need to flip the numbers
    #for cid, x, y, url in zip(stn_ids, xcoords, ycoords, urls):
    for x, y, url, stn, date, m_url, t_url, ptype, ptypet in zip(
            xcoords, ycoords, urls, stns, dates, mtr_urls, text_urls, ptypes, ptypest):
        f.write('''<AREA shape="circle" coords="%d,%d,5"
    onmouseover="javascript:mouseover('%s');"
    href="javascript:show('%s','%s','%s');">\n''' %
                (x, img_height-y, "{}/{} on {}</br>[R Z P S] Layer-Energy: {} Traditional: {}".format(
                        stn[0], stn[1], date, ptype, ptypet), url, m_url, t_url))
    
    f.write("</MAP>\n</BODY></HTML>\n")
    
    f.close()    
# import scipy.optimize as optimize

def expFunc(x, a, b, c):
    
    return a*np.exp(-b*x)-c

def plotMinEnergyLine(ax):
    
    xCrit = np.linspace(0, 4, num=5)

    ctok = np.full_like(xCrit, 273.15)

    yCrit = 1000*9.8*(xCrit/9.8)*(((xCrit+273.15+ctok)/2)-273.15)/273.15

    print xCrit + 273.15, yCrit

    ax.plot(xCrit, yCrit)
    
# Determine the chance that a melted or partially melted hydrometeor refreezes.
# Return 0 for unlikely (<20%), 1 for chance (>= 20 to < 75%), and 2 for likely.

def getProbRefreeze(positiveEnergiesAloft, negativeEnergies):

    negativeEnergies = [ -x for x in negativeEnergies ]
    
    probPL = np.zeros_like(positiveEnergiesAloft, dtype=int)

    for i in range(0, len(positiveEnergiesAloft)):
        
        if (negativeEnergies[i] >= 75+(0.128*positiveEnergiesAloft[i])):
            probPL[i] = 2
        elif (negativeEnergies[i] >= 20+(0.128*positiveEnergiesAloft[i])):
            probPL[i] = 1
        
    return probPL
    
def addBestFitLine(x, y, plt):
    
    ax = plt.gca()
    
    ax.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 2))(np.unique(x)),c='k')
    
    ax.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)),c='k')
    
    # R-Square
 
    correlation = np.corrcoef(x, y)[0,1]
    rsq = correlation**2
    
    """
    plt.text(0.8, 0.2, r'$R^2$'+'Value = %0.4f' %rsq, transform=ax.transAxes,
         horizontalalignment='center',
         verticalalignment='center')
    """

def annotateRainSnow(snd_count, maxSnowT, minRainT, plt):

    ax = plt.gca()

    xMin, xMax = ax.get_xlim()
    
    snPosX = 0.5 * (xMin + maxSnowT)
    mixPosX = 0.5 * (maxSnowT + minRainT)
    raPosX = 0.5 * (minRainT + xMax)
    
    snLength = 0.25 * (maxSnowT - xMin)
    mixLength = 0.25 * (minRainT - maxSnowT)
    raLength = 0.25 * (xMax - minRainT)
    
    snArrow = .7 * snLength
    mixArrow = .7 * mixLength
    raArrow = .7 * raLength
    
    snColor = 'steelblue'
    mixColor = '#386092'
    raColor = 'limegreen'
    
    ax.text(snPosX, 61, "SNOW", fontsize=10, 
            color=snColor, fontweight='bold', ha='center',
            bbox=dict(alpha=0.9, facecolor='w', edgecolor='none', pad=0.0))
    ax.text(mixPosX, 61, "MIX", fontsize=10, 
            color=mixColor, fontweight='bold', ha='center',
            bbox=dict(alpha=0.9, facecolor='w', edgecolor='none', pad=0.0))
    ax.text(raPosX, 61, "RAIN", fontsize=10, 
             color=raColor, fontweight='bold', ha='center',
             bbox=dict(alpha=0.9, facecolor='w', edgecolor='none', pad=0.0))

    ax.arrow(xMin + snLength, 62, -snArrow, 0, head_width=3, head_length=0.05, fc=snColor, ec=snColor)
    ax.arrow(maxSnowT - snLength, 62, snArrow, 0, head_width=3, head_length=0.05, fc=snColor, ec=snColor)

    ax.arrow(maxSnowT + mixLength, 62, -mixArrow, 0, head_width=3, head_length=0.05, fc=mixColor, ec=mixColor)
    ax.arrow(minRainT - mixLength, 62, mixArrow, 0, head_width=3, head_length=0.05, fc=mixColor, ec=mixColor)
    
    ax.arrow(minRainT + raLength, 62, -raArrow, 0, head_width=3, head_length=0.05, fc=raColor, ec=raColor)
    ax.arrow(xMax - raLength, 62, raArrow, 0, head_width=3, head_length=0.05, fc=raColor, ec=raColor)
    
    ax.axvline(maxSnowT, c='b')
    ax.axvline(minRainT, c='g')

    ax.text(0.08, 0.9, "{} soundings".format(snd_count), fontsize=10, 
             transform=ax.transAxes, backgroundcolor='w', 
             bbox=dict(alpha=0.9, facecolor='w'))
    
def addMeltingRects(xMax, plt):
    
    #fig, ax = plt.subplots(figsize=(widthInches, heightInches), dpi=dpi) # plt.gca()
    
    ax = plt.gca()
    
    poly1 = plt.Rectangle((0, 0), xMax, 5, color='b', alpha=.5)
    poly2 = plt.Rectangle((0, 5), xMax, 5, color='dodgerblue', alpha=.5)
    poly3 = plt.Rectangle((0, 10), xMax, 3.5, color='c', alpha=.5)
    poly4 = plt.Rectangle((0, 13.5), xMax, 6.5, color='lightgreen', alpha=.5)

    ax.add_patch(poly1)
    ax.add_patch(poly2)
    ax.add_patch(poly3)
    ax.add_patch(poly4)

    ax.text(0.87, 0.01, "SNOW", fontsize=10, 
             transform=ax.transAxes, color='#C6D9F1', fontweight='bold', ha='center')
    ax.text(0.87, 0.076, "SN/CHC RA", fontsize=10, 
             transform=ax.transAxes, color='#4F81BD', fontweight='bold', ha='center')
    ax.text(0.87, 0.13, "RAIN/SNOW", fontsize=10, 
             transform=ax.transAxes, color='#386092', fontweight='bold', ha='center')    
    ax.text(0.87, 0.19, "RA/CHC SN", fontsize=10, 
             transform=ax.transAxes, color='darkgreen', fontweight='bold', ha='center')    
    ax.text(0.87, 0.26, "RAIN", fontsize=10, 
             transform=ax.transAxes, color='limegreen', fontweight='bold', ha='center')
    
def ptypes_on_sfct_vs_fzlevel(ptypes, temps, fzlevels, tempLabel):
    
#    for i, t in enumerate(temps):
#        print temps[i], fzlevels[i], ptypes[i]
    
    return
        
def aloft_vs_total(posAloft, posTotal, refreezeLowest, sfcTemps, tempLabel, stns, dates):

    ax = plt.gca()
    
    ax.set_xlim(0, 30)    
    ax.set_ylim(0, 30)
    
    pAloft = []
    pTotal = []
    
    # Filter out boring duplicates
    
    for x in xrange(0, len(posAloft)):        
        if posAloft[x] != posTotal[x]:
            pAloft.append(posAloft[x])
            pTotal.append(posTotal[x])
            print "{:.1f} Aloft and {:.1f} Total with Refreeze {:.1f} and Sfc {:.1f} at {} on {}".format(
                    posAloft[x], posTotal[x], refreezeLowest[x], sfcTemps[x],
                    stns[x], dates[x])
    
    addBestFitLine(pAloft, pTotal, plt)
    
    ax.scatter(pAloft, pTotal)
    
    ax.set_ylabel("Melting {}-Energy [J] Total".format(tempLabel))
    ax.set_xlabel("Melting {}-Energy [J] Aloft".format(tempLabel))  
    
    plt.show()

def aloft_vs_sfc(posAloft, posSfc, refreezeEnergies, tempLabel):

    pAloft = np.array(posAloft)
    pSfc = np.array(posSfc)
    rEnergies = np.array(refreezeEnergies)
    
    probRefreeze = getProbRefreeze(pAloft, rEnergies)
    
    xMax = 50
    yMax = 50
    
    mask = np.logical_and(np.logical_and(pAloft <= yMax, pAloft > 0),
                          np.logical_and(pSfc <= xMax, pSfc > 0))
    
    pAloft = pAloft[mask]
    pSfc = pSfc[mask]
    probRefreeze = probRefreeze[mask]
    
    cmap = LinearSegmentedColormap.from_list('name', ['green', 'cyan', 'blue'])

    ax = plt.gca()
    
    ax.scatter(pSfc, pAloft, c=probRefreeze, s=10 + 50*probRefreeze, cmap=cmap)
    
    ax.set_xlabel("Melting {}-Energy [J] Surface".format(tempLabel))
    ax.set_ylabel("Melting {}-Energy [J] Aloft".format(tempLabel))  

    ax.text(0.1, 0.9, "{} soundings".format(len(pAloft)), fontsize=10, 
             transform=ax.transAxes, backgroundcolor='w', 
             bbox=dict(alpha=0.9, facecolor='w'))
    
    ax.plot([0,5],[5,0],[0,10],[10,0],[0,13.5],[13.5,0],[0,20],[20,0])
    
    plt.show()
    
def maxenergy_vs_minenergy(maxEnergy, minEnergy, tempLabel):
    
    ax = plt.gca()
    
    ax.scatter(maxEnergy, minEnergy, s=10)
    
    ax.set_ylim(-210, 0)
    ax.set_xlim(0, 500)
    
    # ax = plt.gca()

    ax.invert_yaxis()
    
    # ax.set_xlim(0,4)
    # ax.set_ylim(0,80)
    
    ax.set_ylabel("Refreeze {}-Energy [J]".format(tempLabel))
    ax.set_xlabel("Melting {}-Energy [J] Aloft".format(tempLabel))
    
    #addMeltingRects(40, plt)
    
    #plt.axvline(1, c='b')
    #plt.axvline(3, c='g')

    ax.text(0.1, 0.9, "{} soundings".format(len(maxEnergy)), fontsize=10, 
             transform=ax.transAxes, backgroundcolor='w', 
             bbox=dict(alpha=0.9, facecolor='w'))
    
    xRange = np.arange(0, 500, 10)
    
    yRange1 = -105 - 0.128 * xRange
    
    yRange2 = -85 - 0.128 * xRange
    
    yRange3 = -58 - 0.128 * xRange
    
    ax.plot(xRange, yRange1)
    ax.plot(xRange, yRange2)
    ax.plot(xRange, yRange3)
    
    plt.show()
    
def maxt_vs_mint(maxTs, minTs, tempLabel):
    
    ax = plt.gca()
    
    ax.scatter(maxTs, minTs, s=10)   

    ax.invert_yaxis()
    
    # ax.set_xlim(0,4)
    # ax.set_ylim(0,80)
    
    ax.set_ylabel("Min Refreeze {} [C]".format(tempLabel))
    ax.set_xlabel("Max{} Aloft [C]".format(tempLabel))
    
    #addMeltingRects(40, plt)
    
    #plt.axvline(1, c='b')
    #plt.axvline(3, c='g')

    ax.text(0.1, 0.9, "{} soundings".format(len(maxTs)), fontsize=10, 
             transform=ax.transAxes, backgroundcolor='w', 
             bbox=dict(alpha=0.9, facecolor='w'))
    
    plt.show()
    
def maxt_vs_energy(maxTs, posEnergies, sfcTmps, tempLabel, stns, dates):

    # From lists, make np arrays for easier manipulation    
    sTemps = np.array(maxTs)
    sEnergies = np.array(posEnergies)
    sfTemps = np.array(sfcTmps)
    sStns = np.array(stns)
    sDates = np.array(dates)   
    
    tmp_indices = (sTemps > 0) & (sfTemps <= 0)
    
    sEnergies = sEnergies[tmp_indices]
    sStns = sStns[tmp_indices]
    sfTemps = sfTemps[tmp_indices]
    sDates = sDates[tmp_indices]
    sTemps = sTemps[tmp_indices]   
    
    # Create plot and retrieve fig and ax objects    
    fig, ax = plt.subplots(figsize=(widthInches, heightInches), dpi=dpi)
    
    # Filter points by some criteria, if desired.
    colors = np.where(sEnergies < 20, 'b', 'g')
    markers = np.where(sEnergies < 20, r'$\ast$', '.')

    unique_markers = set(markers)

    # Plot points
    for um in unique_markers:
        mask = markers == um
        plt.scatter(sTemps[mask], sEnergies[mask], c=colors[mask], marker=um)
    
    # Create titles and labels
    axTitle = "{} Aloft vs Positive {}-Energy Aloft".format(tempLabel, tempLabel)    
    plt.title(axTitle)
    ax.set_ylabel("Positive Energy Aloft (J)")
    ax.set_xlabel("{} Aloft [C]".format(tempLabel))
    
    xMax = 4
    yMax = 80
    ax.set_xlim(0,xMax)
    ax.set_ylim(0,yMax)
    
    # Add some annotations to the figure
    
    addMeltingRects(xMax, plt)
    
    snowRainTemp = 1
    rainOnlyTemp = 3
    
    snd_count = len(sEnergies[sEnergies < yMax])
    
    annotateRainSnow(snd_count, snowRainTemp, rainOnlyTemp, plt)

    # Save the figure to a file
    figname= "plots\\max{}vsEnergy".format(tempLabel)    
    fig.savefig(figname + '.png', dpi=dpi)
    
    img_width = fig.get_figwidth() * dpi
    img_height = fig.get_figheight() * dpi
    
    # Convert the data set points into screen space coordinates
    trans = ax.transData 
    xcoords, ycoords = zip(*trans.transform(zip(sTemps, sEnergies)))
    
    # Create the webpage for displaying this information
    createWebpage(xcoords, ycoords, img_width, img_height, figname, axTitle, sStns, sDates)
    
    plt.show()
    
def sfct_vs_energyAloft(sfcTemps, energiesAloft, tempLabel):
    
    temps = np.array(sfcTemps)

    energies = np.array(energiesAloft)
    
    temps = temps[energies > 0]
    
    energies = energies[energies > 0]
    
    ax = plt.gca()
    
    ax.scatter(temps, energies)    
    
    ax.set_ylabel("Positive Energy Aloft (J)")
    
    ax.set_xlabel("Sfc{} [C]".format(tempLabel))
    
    plt.show()
    
def sfct_vs_energy(sfcTemps, sfcTempEnergies, minEnergies, tempLabel, stns, dates):
    
    # Create nparrays from lists for easier handling
    
    sTemps = np.array(sfcTemps)
    sEnergies = np.array(sfcTempEnergies)
    mEnergies = np.array(minEnergies)
    sStns = np.array(stns)
    sDates = np.array(dates)
    
    # Filter to only plot sfc energies > minEnergy for a given temperature
    indices = sEnergies > mEnergies
    sTemps = sTemps[indices]
    sEnergies = sEnergies[indices]
    sStns = sStns[indices]
    sDates = sDates[indices]    
    
    # Create plot and retrieve fig and ax objects    
    fig, ax = plt.subplots(figsize=(widthInches, heightInches), dpi=dpi)
    
    colors = np.where(sEnergies < 20, 'b', 'g')
    markers = np.where(sEnergies < 20, r'$\ast$', '.')

    #plt.scatter(sfcTemps, sfcTempEnergies)
    unique_markers = set(markers)

    for um in unique_markers:
        mask = markers == um
        plt.scatter(sTemps[mask], sEnergies[mask], c=colors[mask], marker=um)
        
    xMax = 4.44
    ax.set_xlim(0, xMax)
    ax.set_ylim(0, 80)

    # ax = plt.gca()
    
    addMeltingRects(xMax, plt)

    snd_count = len(sTemps)
    
    snowRainTemp = 1.11
    rainOnlyTemp = 3.33
    
    annotateRainSnow(snd_count, snowRainTemp, rainOnlyTemp, plt)
    
    xCrit = np.linspace(0, 4, num=5)

    ctok = np.full_like(xCrit, 273.15)

    yCrit = 1000*9.8*(xCrit/9.8)*(((xCrit+273.15+ctok)/2)-273.15)/273.15

    print xCrit + 273.15, yCrit

    ax.plot(xCrit, yCrit)

    # Create titles and labels
    axTitle = "Sfc{} vs Positive {}-Energy Sfc".format(tempLabel, tempLabel)    
    ax.set_title(axTitle, y = 1.08)
    ax2 = ax.twiny()
    ax.set_ylabel("Sfc-Based Energy (J)")
    ax.set_xlabel("Sfc {} [C]".format(tempLabel))
    ax2.set_xlabel("Sfc {} [F]".format(tempLabel))
    #plt.title(title)
    
    # set twin scale (convert degree celsius to fahrenheit)
    T_f = lambda T_c: T_c*1.8 + 32.
    # get bottom axis limits
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    # apply function and set transformed values to top axis limits
    ax2.set_xlim((T_f(xmin),T_f(xmax)))
    ax2.set_ylim(ymin,ymax)
    # set an invisible artist to twin axes 
    # to prevent falling back to initial values on rescale events
    ax2.plot([],[])
    
    # addBestFitLine(sTemps, sEnergies, plt)    
    
    #plt.show()

    # Save the figure to a file
    figname= "plots\\sfc{}vsSfcEnergy".format(tempLabel)    
    fig.savefig(figname + '.png', dpi=dpi)
    
    img_width = fig.get_figwidth() * dpi
    img_height = fig.get_figheight() * dpi
    
    # Convert the data set points into screen space coordinates
    trans = ax.transData 
    xcoords, ycoords = zip(*trans.transform(zip(sTemps, sEnergies)))
    
    # Create the webpage for displaying this information
    createWebpage(xcoords, ycoords, img_width, img_height, figname, axTitle, sStns, sDates)
    
    plt.show()
    
def t_vs_refreeze_depth(minRefreezeTs, refreezeTDepths, energies, tempLabel, maxTemps, stns, dates, ptypes):
    
    # From lists make np arrays for easier handling.    
    minTs = np.array(minRefreezeTs)    
    depths = np.array(refreezeTDepths)    
    invEnergies = np.array(energies) * -1    
    #minTs = np.array(minRefreezeTs)
    #depths = np.array(refreezeTDepths)
    maxTs = np.array(maxTemps)
    sStns = np.array(stns)
    sDates = np.array(dates)
    pTypes = np.array(ptypes)
    
    indices = maxTs > 3
    minTs = minTs[indices]
    depths = depths[indices]
    invEnergies = invEnergies[indices]
    sStns = sStns[indices]
    sDates = sDates[indices]
    pTypes = pTypes[indices]

    # Create plot and retrieve fig and ax objects    
    fig, ax = plt.subplots(figsize=(widthInches, heightInches), dpi=dpi)
    
    xMax = 0
    ax.set_xlim(-20, xMax)
    ax.set_ylim(0, 3000)
    
    #ax = plt.gca()   
    
    # ax.set_xlim(0, 4.44)
    # ax.set_ylim(0, 160)
    
    ax.invert_xaxis()

    # Plot all points using a gradient color scale and dot size indicating the
    # amount of negative energy
    
    plt.scatter(minTs, depths, s=invEnergies, c=invEnergies)

    """
    plt.scatter(mrts[invEnergies > 120], rtds[invEnergies > 120], 
                s=invEnergies[invEnergies > 120], facecolors='none', 
                edgecolors='k', linewidth=3)
    """
    
    # Find points where pType is likely sleet. Plot these with a black outline.
    
    indices = [x[2] > 1 for x in pTypes]
    
    plt.scatter(minTs[indices], depths[indices], 
                s=invEnergies[indices], facecolors='none', 
                edgecolors='k', linewidth=3)
    
    # Find points where pType is not sleet. Plot these as white but a little
    # smaller than existing points so you still see the outline of the original
    # points at these locations.
    
    #plt.scatter(minTs, depths, s=invEnergies[invEnergies < 60]*.8, color='w')
    indices = [x[2] == 0 for x in pTypes]
    
    plt.scatter(minTs[indices], depths[indices], s=invEnergies[indices]*.8, color='w')
    
    ax2 = ax.twinx()

    # Create titles and labels
    axTitle = "MinRefreeze{} vs Refreeze {}-Energy (Max{}Aloft > 3)".format(tempLabel, tempLabel, tempLabel)
    plt.title(axTitle)
    ax.set_ylabel("Refreeze Layer Depth (m)")
    ax.set_xlabel("Min{} in Refreeze Layer [C]".format(tempLabel))
    ax2.set_ylabel("Refreeze Layer Depth (ft)")
    
    H_ft = lambda H_m: H_m*3.28

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    
    # apply function and set transformed values to top axis limits
    ax2.set_ylim((H_ft(ymin),H_ft(ymax)))
    #ax2.set_xlim(ymin,ymax)
    
    ax2.set_yticks(np.arange(ax2.get_ylim()[0],ax2.get_ylim()[1],1000))
    
    ax.set_xticks(np.arange(-20,2,2))
    
    # set an invisible artist to twin axes 
    # to prevent falling back to initial values on rescale events
    ax2.plot([],[])
    
    poly1 = plt.Rectangle((-4, 2500), -2, ax2.get_ylim()[1], color='c', alpha=.5, fill=None)
    poly2 = plt.Rectangle((-6, 2500), -2, ax2.get_ylim()[1], color='dodgerblue', alpha=.5, fill=None)
    poly3 = plt.Rectangle((-8, 2500), xmax, ax2.get_ylim()[1], color='b', alpha=.5, fill=None)
    
    ax2.add_patch(poly1)
    ax2.add_patch(poly2)
    ax2.add_patch(poly3)
    
    snd_count = len(minTs)
    
    plt.text(0.1, 0.9, "{} soundings".format(snd_count), fontsize=10, 
             transform=ax.transAxes, backgroundcolor='w', 
             bbox=dict(alpha=0.9, facecolor='w'))
    
    #plt.show()
    
    # Save the figure to a file
    figname= "plots\\Min{}vsRefreezeEnergy".format(tempLabel)    
    fig.savefig(figname + '.png', dpi=dpi)
    
    img_width = fig.get_figwidth() * dpi
    img_height = fig.get_figheight() * dpi
    
    # Convert the data set points into screen space coordinates
    trans = ax.transData 
    xcoords, ycoords = zip(*trans.transform(zip(minTs, depths)))
    
    # Create the webpage for displaying this information
    createWebpage(xcoords, ycoords, img_width, img_height, figname, axTitle, sStns, sDates, pTypes)
    
    plt.show()
    
def maxt_vs_depth(temps, depths, tempLabel, stns, dates):
    
    # Create nparrays from lists for easier handling
    
    sTemps = np.array(temps)

    sDepths = np.array(depths)
    
    plt.scatter(sTemps, sDepths, s=10)
    
    ax = plt.gca()
    
    #ax.set_xlim(0,4)
    #ax.set_ylim(0,2000)
    
    ax.set_xlabel("Max {} Aloft [C]".format(tempLabel))
    ax.set_ylabel("Depth of Elevated Warm Layer (m)")
    
    plt.text(0.1, 0.9, "{} soundings".format(len(temps)), fontsize=10, 
             transform=ax.transAxes, backgroundcolor='w', 
             bbox=dict(alpha=0.9, facecolor='w'))
        
    addBestFitLine(sTemps, sDepths, plt)
    
    plt.show()

def maxtXdepth_vs_energies(temps, depths, energies, tempLabel):
    
    t = np.array(temps)
    d = np.array(depths)
    
    txd = t * d
    
    plt.scatter(txd, energies)
    
    ax = plt.gca()
    
    ax.set_xlabel("{} Aloft [C] X Depth of Melting {}-Layer [m]".format(tempLabel, tempLabel))
    ax.set_ylabel("{}-based Melting Energy [J]".format(tempLabel))
    
    plt.text(0.1, 0.9, "{} soundings".format(len(temps)), fontsize=10, 
             transform=ax.transAxes, backgroundcolor='w', 
             bbox=dict(alpha=0.9, facecolor='w'))
    
    addBestFitLine(txd, energies, plt)
    
    plt.show()
    
def energy_vs_depth(energies, depths, tempLabel):

    sEnergies = np.array(energies)
    
    sDepths = np.array(depths)
    
    plt.scatter(sEnergies, sDepths, s=10)
    
    ax = plt.gca()
    
    ax.set_xlabel("{}-Based Energy Aloft [J]".format(tempLabel))
    ax.set_ylabel("Depth of Elevated Warm Layer (m)")

    plt.text(0.1, 0.9, "{} soundings".format(len(energies)), fontsize=10, 
             transform=ax.transAxes, backgroundcolor='w', 
             bbox=dict(alpha=0.9, facecolor='w'))
    
    # addBestFitLine(sEnergies, sDepths, plt)
    
    plt.show()    
    
def depth_vs_energy(depths, energies, tempLabel):

    sEnergies = np.array(energies)
    
    sDepths = np.array(depths)
    
    #sTemps = np.array(maxTs)

    #sEnergies = np.array(posEnergies)

    #mEnergies = np.array(minEnergies)
    
    # Filter to only plot sfc energies > minEnergy for a given temperature
    
    #sTemps = sTemps[sEnergies > mEnergies]    
    
    #sEnergies = sEnergies[sEnergies > mEnergies]
    
    colors = np.where(sEnergies < 20, 'b', 'g')
    markers = np.where(sEnergies < 20, r'$\ast$', '.')

    #plt.scatter(sfcTemps, sfcTempEnergies)
    unique_markers = set(markers)

    for um in unique_markers:
        mask = markers == um
        plt.scatter(sDepths[mask], sEnergies[mask], c=colors[mask], marker=um)
        
    # plt.scatter(sDepths, sEnergies, s=10)
    
    xmax = 1000

    ax = plt.gca()
    
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, 80)
    
    ax = plt.gca()
    
    ax.set_ylabel("{}-Based Energy Aloft [J]".format(tempLabel))
    ax.set_xlabel("Depth of Elevated {} Warm Layer (m)".format(tempLabel))
    
    addMeltingRects(xmax, plt)

    plt.text(0.1, 0.9, "{} soundings".format(len(depths)), fontsize=10,
             transform=ax.transAxes, backgroundcolor='w', 
             bbox=dict(alpha=0.9, facecolor='w'))
    
    # addBestFitLine(sDepths, sEnergies, plt)
    
    plt.show()        

def t_vs_lowEnergy(temps, energies, eAloft, minValues, tempLabel, stns, dates, xMax, yMax, snowRainTemp, rainOnlyTemp, tempLayer, energyLayer, ptypes, ptypest):
    sTemps = np.array(temps)
    sEnergies = np.array(energies)
    mValues = np.array(minValues)
    sStns = np.array(stns)
    sDates = np.array(dates)
    sTypes = np.array(ptypes)
    sTypest = np.array(ptypest)
    eAlofts = np.array(eAloft)
    
    indices = (sEnergies > 5) & (eAlofts <= 5)
    sTemps = sTemps[indices]
    sEnergies = sEnergies[indices]
    sStns = sStns[indices]
    sDates = sDates[indices]    
    sTypes = sTypes[indices]
    sTypest = sTypest[indices]   
    eAlofts = eAlofts[indices]
    
    # Create plot and retrieve fig and ax objects    
    fig, ax = plt.subplots(figsize=(widthInches, heightInches), dpi=dpi)
    
    colors = np.where(sEnergies < 20, 'b', 'g')
    markers = np.where(sEnergies < 20, r'$\ast$', '.')

    unique_markers = set(markers)

    for um in unique_markers:
        mask = markers == um
        plt.scatter(sTemps[mask], sEnergies[mask], c=colors[mask], marker=um)
        
    #xMax = 4.44
    #yMax = 80
    ax.set_xlim(-10, xMax)
    ax.set_ylim(0, yMax)

    # ax = plt.gca()
    
    addMeltingRects(xMax, plt)

    snd_count = len(sTemps)
    
    #snowRainTemp = 1.11
    #rainOnlyTemp = 3.33
    
    annotateRainSnow(snd_count, snowRainTemp, rainOnlyTemp, plt)
    
    #tempLayer = "Aloft"
    titlePosition = 1
    
    if tempLayer == "Sfc":
        plotMinEnergyLine(ax)
        titlePosition = 1.08        
        ax2 = ax.twiny()
        ax2.set_xlabel("{} {} [F]".format(tempLabel, tempLayer))
    
    # set twin scale (convert degree celsius to fahrenheit)
        T_f = lambda T_c: T_c*1.8 + 32.
        # get bottom axis limits
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        # apply function and set transformed values to top axis limits
        ax2.set_xlim((T_f(xmin),T_f(xmax)))
        ax2.set_ylim(ymin,ymax)
        # set an invisible artist to twin axes 
        # to prevent falling back to initial values on rescale events
        ax2.plot([],[])
                
    # Create other titles and labels
    axTitle = "{} {} vs Positive {}-Energy {}".format(tempLayer, tempLabel, tempLabel, energyLayer)
    plt.title(axTitle, y = titlePosition)
    ax.set_ylabel("{}-Based Positive Energy (J)".format(energyLayer))
    ax.set_xlabel("{} {} [C]".format(tempLayer, tempLabel))
    #ax2.set_xlabel("{} {} [F]".format(tempLabel))
    #plt.title(title)
    
    # addBestFitLine(sTemps, sEnergies, plt)    
    
    #plt.show()

    # Save the figure to a file
    figname= "plots\{}{}vs{}LowEnergy{}".format(tempLayer, tempLabel, tempLabel, energyLayer)
    fig.savefig(figname + '.png', dpi=dpi)
    
    img_width = fig.get_figwidth() * dpi
    img_height = fig.get_figheight() * dpi
    
    # Convert the data set points into screen space coordinates
    trans = ax.transData 
    xcoords, ycoords = zip(*trans.transform(zip(sTemps, sEnergies)))
    
    # Create the webpage for displaying this information
    createWebpage(xcoords, ycoords, img_width, img_height, figname, axTitle, sStns, sDates, sTypes, sTypest)

    plt.show()    
def t_vs_posEnergy(temps, energies, minValues, tempLabel, stns, dates, xMax, yMax, snowRainTemp, rainOnlyTemp, tempLayer, energyLayer, ptypes, ptypest):
#def t_vs_posEnergy(**chartData):
    # Create nparrays from lists for easier handling    
    sTemps = np.array(temps)
    sEnergies = np.array(energies)
    mValues = np.array(minValues)
    sStns = np.array(stns)
    sDates = np.array(dates)
    sTypes = np.array(ptypes)
    sTypest = np.array(ptypest)
    
    # Filter to only plot sfc energies > minEnergy for a given temperature
    indices = sEnergies > mValues
    sTemps = sTemps[indices]
    sEnergies = sEnergies[indices]
    sStns = sStns[indices]
    sDates = sDates[indices]    
    sTypes = sTypes[indices]
    sTypest = sTypest[indices]
    
    # Create plot and retrieve fig and ax objects    
    fig, ax = plt.subplots(figsize=(widthInches, heightInches), dpi=dpi)
    
    colors = np.where(sEnergies < 20, 'b', 'g')
    markers = np.where(sEnergies < 20, r'$\ast$', '.')

    unique_markers = set(markers)

    for um in unique_markers:
        mask = markers == um
        plt.scatter(sTemps[mask], sEnergies[mask], c=colors[mask], marker=um)
        
    #xMax = 4.44
    #yMax = 80
    ax.set_xlim(-10, xMax)
    ax.set_ylim(0, yMax)

    # ax = plt.gca()
    
    addMeltingRects(xMax, plt)

    snd_count = len(sTemps)
    
    #snowRainTemp = 1.11
    #rainOnlyTemp = 3.33
    
    annotateRainSnow(snd_count, snowRainTemp, rainOnlyTemp, plt)
    
    #tempLayer = "Aloft"
    titlePosition = 1
    
    if tempLayer == "Sfc":
        plotMinEnergyLine(ax)
        titlePosition = 1.08        
        ax2 = ax.twiny()
        ax2.set_xlabel("{} {} [F]".format(tempLabel, tempLayer))
    
    # set twin scale (convert degree celsius to fahrenheit)
        T_f = lambda T_c: T_c*1.8 + 32.
        # get bottom axis limits
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        # apply function and set transformed values to top axis limits
        ax2.set_xlim((T_f(xmin),T_f(xmax)))
        ax2.set_ylim(ymin,ymax)
        # set an invisible artist to twin axes 
        # to prevent falling back to initial values on rescale events
        ax2.plot([],[])
                
    # Create other titles and labels
    axTitle = "{} {} vs Positive {}-Energy {}".format(tempLayer, tempLabel, tempLabel, energyLayer)
    plt.title(axTitle, y = titlePosition)
    ax.set_ylabel("{}-Based Positive Energy (J)".format(energyLayer))
    ax.set_xlabel("{} {} [C]".format(tempLayer, tempLabel))
    #ax2.set_xlabel("{} {} [F]".format(tempLabel))
    #plt.title(title)
    
    # addBestFitLine(sTemps, sEnergies, plt)    
    
    #plt.show()

    # Save the figure to a file
    figname= "plots\{}{}vs{}Energy{}".format(tempLayer, tempLabel, tempLabel, energyLayer)
    fig.savefig(figname + '.png', dpi=dpi)
    
    img_width = fig.get_figwidth() * dpi
    img_height = fig.get_figheight() * dpi
    
    # Convert the data set points into screen space coordinates
    trans = ax.transData 
    xcoords, ycoords = zip(*trans.transform(zip(sTemps, sEnergies)))
    
    # Create the webpage for displaying this information
    createWebpage(xcoords, ycoords, img_width, img_height, figname, axTitle, sStns, sDates, sTypes, sTypest)

    plt.show()

def total_vs_below(totalAloft, totalBelow, tempLabel, levelLabel, sStns, sDates, pTypes, pTypest):
       
    # Create plot and retrieve fig and ax objects    
    fig, ax = plt.subplots(figsize=(widthInches, heightInches), dpi=dpi)
    
    totalAloft, totalBelow, sStns, sDates = zip(*[(i,j,k,l) for i,j,k,l 
                                                  in zip(totalAloft, totalBelow, sStns, sDates) if ((i < 0) and (i != j))])
    
    plt.scatter(totalAloft, totalBelow, s=5)

    axTitle = "Refreezing {} Below Lowest Melting Layer vs {} Below {}".format(tempLabel, tempLabel, levelLabel)
    
    ax.invert_xaxis()
    ax.invert_yaxis()
    
    ax.set_xlim(0, -500)
    ax.set_ylim(0, -500)
    
    plt.title(axTitle)
    
    # ax = plt.gca()
    
    ax.set_xlabel("Refreezing {} below Lowest Melting Layer".format(tempLabel))
    ax.set_ylabel("Refreezing {} below {}".format(tempLabel, levelLabel))
    
    snd_count = len(totalAloft)
    
    plt.text(0.08, 0.9, "{} soundings".format(snd_count), fontsize=10, 
             transform=ax.transAxes, backgroundcolor='w', 
             bbox=dict(alpha=0.9, facecolor='w'))    
    
    # Save the figure to a file
    figname= "plots\\{}TotalvsBelow{}".format(tempLabel, levelLabel)
    
    fig.savefig(figname + '.png', dpi=dpi)
    
    img_width = fig.get_figwidth() * dpi
    img_height = fig.get_figheight() * dpi
    
    # Convert the data set points into screen space coordinates
    trans = ax.transData 
    xcoords, ycoords = zip(*trans.transform(zip(totalAloft, totalBelow)))
    
    # Create the webpage for displaying this information
    createWebpage(xcoords, ycoords, img_width, img_height, figname, axTitle, sStns, sDates, pTypes, pTypest)
    
    plt.show()
    
def aloft_vs_above(alofts, aboves, tempLabel, levelLabel, sStns, sDates, pTypes, pTypest):
    
    # Create plot and retrieve fig and ax objects    
    fig, ax = plt.subplots(figsize=(widthInches, heightInches), dpi=dpi)
    
    alofts, aboves, sStns, sDates = zip(*[(i,j,k,l) for i,j,k,l 
                                                  in zip(alofts, aboves, sStns, sDates)]) # if (i > 0)]) # & (j > 0))])
    
    ax.scatter(alofts, aboves, s=5)

    axTitle = "{} Aloft vs {} Above {}".format(tempLabel, tempLabel, levelLabel)
    
    ax.set_title(axTitle)
    
    # ax = plt.gca()
    
    ax.set_xlabel("{} above First Freezing Layer".format(tempLabel))
    ax.set_ylabel("{} above {}".format(tempLabel, levelLabel))
    
    snd_count = len(alofts)
    
    ax.text(0.08, 0.9, "{} soundings".format(snd_count), fontsize=10, 
             transform=ax.transAxes, backgroundcolor='w', 
             bbox=dict(alpha=0.9, facecolor='w'))
    
    unmatchedCount = 0
    
    for i,j in itertools.izip(alofts, aboves):
        if (i != j):
            unmatchedCount += 1
    
    ax.text(0.08, 0.80, "{} with differences".format(unmatchedCount), fontsize=10, 
             transform=ax.transAxes, backgroundcolor='w', 
             bbox=dict(alpha=0.9, facecolor='w'))
    
    # Save the figure to a file
    figname= "plots\\{}AloftvsAbove{}".format(tempLabel, levelLabel)
    
    fig.savefig(figname + '.png', dpi=dpi)
    
    img_width = fig.get_figwidth() * dpi
    img_height = fig.get_figheight() * dpi
    
    # Convert the data set points into screen space coordinates
    trans = ax.transData 
    xcoords, ycoords = zip(*trans.transform(zip(alofts, aboves)))
    
    # Create the webpage for displaying this information
    createWebpage(xcoords, ycoords, img_width, img_height, figname, axTitle, sStns, sDates, pTypes, pTypest)
    
    plt.show()
    
def aloft_not_above(alofts, aboves, tempLabel, levelLabel, sStns, sDates, pTypes, pTypest):
    
    # Create plot and retrieve fig and ax objects    
    fig, ax = plt.subplots(figsize=(widthInches, heightInches), dpi=dpi)
    
    alofts, aboves, sStns, sDates = zip(*[(i,j,k,l) for i,j,k,l 
                                                  in zip(alofts, aboves, sStns, sDates) if (i > 0) & (j <= 0)])
    
    ax.scatter(alofts, aboves, s=5)

    axTitle = "{} Aloft not {} Above {}".format(tempLabel, tempLabel, levelLabel)
    
    ax.set_title(axTitle)
    
    # ax = plt.gca()
    
    ax.set_xlabel("{} above First Freezing Layer".format(tempLabel))
    ax.set_ylabel("{} above {}".format(tempLabel, levelLabel))
    
    snd_count = len(alofts)
    
    ax.text(0.08, 0.9, "{} soundings".format(snd_count), fontsize=10, 
             transform=ax.transAxes, backgroundcolor='w', 
             bbox=dict(alpha=0.9, facecolor='w'))
    
    unmatchedCount = 0
    
    for i,j in itertools.izip(alofts, aboves):
        if (i != j):
            unmatchedCount += 1
    
    ax.text(0.08, 0.80, "{} with differences".format(unmatchedCount), fontsize=10, 
             transform=ax.transAxes, backgroundcolor='w', 
             bbox=dict(alpha=0.9, facecolor='w'))
    
    # Save the figure to a file
    figname= "plots\\{}AloftNotAbove{}".format(tempLabel, levelLabel)
    
    fig.savefig(figname + '.png', dpi=dpi)
    
    img_width = fig.get_figwidth() * dpi
    img_height = fig.get_figheight() * dpi
    
    # Convert the data set points into screen space coordinates
    trans = ax.transData 
    xcoords, ycoords = zip(*trans.transform(zip(alofts, aboves)))
    
    # Create the webpage for displaying this information
    createWebpage(xcoords, ycoords, img_width, img_height, figname, axTitle, sStns, sDates, pTypes, pTypest)
    
    plt.show()    
    
def hist(element, xLabel, tempLabel):
    
    element = [abs(x) for x in element if abs(x) > 0] 
    
    maxValue = int(max(element))
    
    tickInterval = 1
    
    if maxValue > 200:
        tickInterval = 20
        maxValue = int(math.ceil(maxValue/float(tickInterval))*tickInterval)
    elif maxValue > 100:
        maxValue = int(round(max(element),-1))
        tickInterval = 10
    elif maxValue > 25:
        maxValue = int(round(max(element),-1))
        tickInterval = 5

    bins = np.arange(0, maxValue, tickInterval)
    
    ax = plt.gca()      
    
    ax.hist(element, bins=bins)
    
    ax.set_xlabel(xLabel)
    
    ax.set_ylabel("Count")
    

    
    ax.text(0.65, 0.9, "{} soundings".format(len(element)), fontsize=10, 
             transform=ax.transAxes, backgroundcolor='w', 
             bbox=dict(alpha=0.9, facecolor='w'))

    ax.xaxis.set_minor_locator(ticker.MultipleLocator(tickInterval))
    
    ax.set_xlim(0, maxValue)
    
    plt.show()

def rare_sleet(sTemps, negEnergies, maxTempsAloft, t, stns, dates, ptypes=None):
    
    sTemps = np.array(sTemps)
    negEnergies = np.array([-x for x in negEnergies])
    maxTempsAloft = np.array(maxTempsAloft)
    stns = np.array(stns)
    dates = np.array(dates)
    ptypes = np.array(ptypes)
    
    # Filter to only plot sfc temps > 0
    indices = [x > 0 for x in sTemps]    
    sTemps = sTemps[indices]
    negEnergies = negEnergies[indices]
    maxTempsAloft = maxTempsAloft[indices]
    stns = stns[indices]
    dates = dates[indices]    
    ptypes = ptypes[indices]
    
    indices = [x >= 1 for x in maxTempsAloft]
    sTemps = sTemps[indices]
    negEnergies = negEnergies[indices]
    maxTempsAloft = maxTempsAloft[indices]
    stns = stns[indices]
    dates = dates[indices]    
    ptypes = ptypes[indices]
    
    indices = [x >= 60 for x in negEnergies]
    sTemps = sTemps[indices]
    negEnergies = negEnergies[indices]
    maxTempsAloft = maxTempsAloft[indices]
    stns = stns[indices]
    dates = dates[indices]    
    ptypes = ptypes[indices]
    
    # Create plot and retrieve fig and ax objects    
    fig, ax = plt.subplots(figsize=(widthInches, heightInches), dpi=dpi)
    
    ax.scatter(sTemps, maxTempsAloft, s=negEnergies, c=negEnergies)
    
    axTitle = "Sfc{} vs Max{} Aloft".format(t, t)
    
    ax.set_title(axTitle)
    
    # plt.gca() 
    
    ax.set_ylabel("Max{} Aloft (C)".format(t))
    ax.set_xlabel("Sfc{} (C)".format(t))
    
    snd_count = len(maxTempsAloft)
    
    ax.text(0.08, 0.9, "{} soundings".format(snd_count), fontsize=10, 
             transform=ax.transAxes, backgroundcolor='w', 
             bbox=dict(alpha=0.9, facecolor='w'))

    # Save the figure to a file
    figname= "plots\\Sfc{}vsMax{}Aloft".format(t, t)
    
    fig.savefig(figname + '.png', dpi=dpi)
    
    img_width = fig.get_figwidth() * dpi
    img_height = fig.get_figheight() * dpi
    
    # Convert the data set points into screen space coordinates
    trans = ax.transData 
    xcoords, ycoords = zip(*trans.transform(zip(sTemps, maxTempsAloft)))
    
    # Create the webpage for displaying this information
    createWebpage(xcoords, ycoords, img_width, img_height, figname, axTitle, stns, dates, ptypes)
    
    plt.show()
    
def cold_not_snow(sTemps, posEnergies, maxTempsAloft, t, stns, dates, ptypes=None):
    
    sTemps = np.array(sTemps)
    posEnergies = np.array(posEnergies)
    maxTempsAloft = np.array(maxTempsAloft)
    stns = np.array(stns)
    dates = np.array(dates)
    ptypes = np.array(ptypes)
    
    # Filter to only plot sfc temps > 0
    indices = [(x < 1) for x in sTemps]    
    sTemps = sTemps[indices]
    posEnergies = posEnergies[indices]
    maxTempsAloft = maxTempsAloft[indices]
    stns = stns[indices]
    dates = dates[indices]    
    ptypes = ptypes[indices]
    
    indices = [(x > -999) and (x < 1) for x in maxTempsAloft]
    sTemps = sTemps[indices]
    posEnergies = posEnergies[indices]
    maxTempsAloft = maxTempsAloft[indices]
    stns = stns[indices]
    dates = dates[indices]    
    ptypes = ptypes[indices]
    
    indices = [x > 5 for x in posEnergies]
    sTemps = sTemps[indices]
    posEnergies = posEnergies[indices]
    maxTempsAloft = maxTempsAloft[indices]
    stns = stns[indices]
    dates = dates[indices]    
    ptypes = ptypes[indices]
    
    # Create plot and retrieve fig and ax objects    
    fig, ax = plt.subplots(figsize=(widthInches, heightInches), dpi=dpi)
    
    ax.scatter(sTemps, maxTempsAloft, s=posEnergies*10, c=posEnergies)
    
    axTitle = "Sfc{} vs Max{} Above 2kft".format(t, t)
    
    ax.set_title(axTitle)
    
    # plt.gca() 
    
    ax.set_ylabel("Max{} Aloft (C)".format(t))
    ax.set_xlabel("Sfc{} (C)".format(t))
    
    snd_count = len(maxTempsAloft)
    
    ax.text(0.08, 0.9, "{} soundings".format(snd_count), fontsize=10, 
             transform=ax.transAxes, backgroundcolor='w', 
             bbox=dict(alpha=0.9, facecolor='w'))

    # Save the figure to a file
    figname= "plots\\Sfc{}vsMax{}AloftCold".format(t, t)
    
    fig.savefig(figname + '.png', dpi=dpi)
    
    img_width = fig.get_figwidth() * dpi
    img_height = fig.get_figheight() * dpi
    
    # Convert the data set points into screen space coordinates
    trans = ax.transData 
    xcoords, ycoords = zip(*trans.transform(zip(sTemps, maxTempsAloft)))
    
    # Create the webpage for displaying this information
    createWebpage(xcoords, ycoords, img_width, img_height, figname, axTitle, stns, dates, ptypes)
    
    plt.show()    
    
def max_sleet(posEnergies, negEnergies, maxTempsAloft, sfcTemps, t, stns, dates, ptypes=None, ptypest=None):
    
    posEnergies = np.array(posEnergies)
    negEnergies = np.array([-x for x in negEnergies])
    maxTempsAloft = np.array(maxTempsAloft)
    sfcTemps = np.array(sfcTemps)
    stns = np.array(stns)
    dates = np.array(dates)
    ptypes = np.array(ptypes)
    ptypest = np.array(ptypest)
    
    # Filter to only plot sfc temps > 0
    
    indices = [x > 0 for x in posEnergies]    
    posEnergies = posEnergies[indices]
    negEnergies = negEnergies[indices]
    maxTempsAloft = maxTempsAloft[indices]
    sfcTemps = sfcTemps[indices]
    stns = stns[indices]
    dates = dates[indices]    
    ptypes = ptypes[indices]
    ptypest = ptypest[indices]
    
    
    indices = [((x >= 1) and (x <= 3)) for x in maxTempsAloft]
    posEnergies = posEnergies[indices]
    negEnergies = negEnergies[indices]
    maxTempsAloft = maxTempsAloft[indices]
    sfcTemps = sfcTemps[indices]
    stns = stns[indices]
    dates = dates[indices]    
    ptypes = ptypes[indices]
    ptypest = ptypest[indices]
    
    indices = [x <= 0 for x in sfcTemps]
    posEnergies = posEnergies[indices]
    negEnergies = negEnergies[indices]
    maxTempsAloft = maxTempsAloft[indices]
    sfcTemps = sfcTemps[indices]
    stns = stns[indices]
    dates = dates[indices]    
    ptypes = ptypes[indices]
    ptypest = ptypest[indices]
    
    """
    indices = [x >= 60 for x in negEnergies]
    posEnergies = posEnergies[indices]
    negEnergies = negEnergies[indices]
    maxTempsAloft = maxTempsAloft[indices]
    stns = stns[indices]
    dates = dates[indices]    
    ptypes = ptypes[indices]
    """

    # Create plot and retrieve fig and ax objects    
    fig, ax = plt.subplots(figsize=(widthInches, heightInches), dpi=dpi)
    
    pm = PtypeMarker()
    
    ra = [x[0] > 0 for x in ptypes]
    zr = [x[1] > 0 for x in ptypes]
    pl = [x[2] > 0 for x in ptypes]
    sn = [x[3] > 0 for x in ptypes]
    
#    plt.scatter(posEnergies[ra], negEnergies[ra], s=80, marker=pm.ra, 
#            edgecolors="green", facecolors='none', linewidth=[x[0] for x in ptypes[ra]])
    ax.scatter(posEnergies[zr], negEnergies[zr], s=80, marker=pm.zr, 
            edgecolors="red", facecolors='none', linewidth=[x[1] for x in ptypes[zr]])
    ax.scatter(posEnergies[pl], negEnergies[pl], s=80, marker=pm.pl, 
            edgecolors="orange", facecolors='none', linewidth=[x[2] for x in ptypes[pl]])
    
    if (sn):
        ax.scatter(posEnergies[sn], negEnergies[sn], s=80, marker=pm.sn, 
            edgecolors="blue", facecolors='none', linewidth=[x[3] for x in ptypes[sn]])        
    
    #colors = np.where(sEnergies < 20, 'b', 'g')
    #markers = np.where(sEnergies < 20, r'$\ast$', '.')

    #unique_markers = set(markers)
    
    #for ptype in ptypes:
    #    mask = markers == um
    #    plt.scatter(sTemps[mask], sEnergies[mask], c=colors[mask], marker=um)
    
#    plt.scatter(posEnergies, negEnergies, s=negEnergies, c=negEnergies)

    xRange = np.arange(0, 150, 10)
    
    yRange1 = 105 + 0.128 * xRange
    
    yRange2 = 85 + 0.128 * xRange
    
    yRange3 = 58 + 0.128 * xRange
    
    ax.plot(xRange, yRange1)
    ax.plot(xRange, yRange2)
    ax.plot(xRange, yRange3)
    
    axTitle = "Positive Energy vs Refreeze Energy ({}-Based) with Max{} Aloft 1-3 C".format(t,t)
    
    ax.set_title(axTitle)
    
    # plt.gca() 
    
    ax.set_xlabel("Positive {}-Energy Aloft (J/kg)".format(t))    
    
    ax.set_ylabel("Refreezing {}-Energy (J/kg)".format(t))    
    
    snd_count = len(maxTempsAloft)
    
    ax.text(0.08, 0.9, "{} soundings".format(snd_count), fontsize=10, 
             transform=ax.transAxes, backgroundcolor='w', 
             bbox=dict(alpha=0.9, facecolor='w'))

    # Save the figure to a file
    figname= "plots\\{}Aloft1to3".format(t)
    
    fig.savefig(figname + '.png', dpi=dpi)
    
    img_width = fig.get_figwidth() * dpi
    img_height = fig.get_figheight() * dpi
    
    # Convert the data set points into screen space coordinates
    trans = ax.transData 
    xcoords, ycoords = zip(*trans.transform(zip(posEnergies, negEnergies)))
    
    # Create the webpage for displaying this information
    createWebpage(xcoords, ycoords, img_width, img_height, figname, axTitle, stns, dates, ptypes, ptypest)
    
    plt.show()    