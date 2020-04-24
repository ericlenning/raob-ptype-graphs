# raob-ptype-graphs
Examine properties of RAOB soundings relevant to evaluating precipitation-type.

The purpose of this project is to test and qualitatively evaluate various methods of forecasting precipitation type using archived RAOB data.  Output can be compared to representative surface observations.  It allows for direct comparisons of different methods (MaxTAloft, original and revised Bourgouin layer-enegy approaches, temperature-based vs wetbulb based, etc.).  

It also allows for finding and examining specific situations.  Just a few examples include:
* Surface temperatures in certain ranges
* MaxTAloft or MaxTwAloft in certain ranges
* Sfc-based energy in certain ranges
* Soundings with elevated warm layers below certain heights
* Soundings above certain elevations

This project also has the potential to provide statistical or climatological information about RAOBs.

Main Programs
-------------
The programs need to run in the order listed, but once the first two prepare the data, the third can be execuated over and over with various additions and tweaks as needed.  So the first two wrangle all the data one time and the third does stuff with the data.

* fsl2spc.py         
  * Translates FSL-formatted soundings into SPC format.
  * Utilizes SPCSounding.py.
  * Reads data from "Soundings/raob/1 - fslFormat".
  * Writes data to "Soundings/raob/2 - spcFormat".
                
* spc2profiles.py
  * Reads SPC-formatted soundings.
  * Utilizes sharpy and spcProfile.py.
  * Does some QC and filtering, only keeping soundings:
    * Saturated enough to support precip;
    * Saturated enough in DGZ to support ice nucleation;
    * Colder than -20C at top and warmer than -12C at sfc;
    * With more than 30 vertical levels.
  * We also interpolate and add levels crossing 0C T and Tw.
  * Reads data from "Soundings/raob/2 - spcFormat".
  * Writes output to:
    * "Soundings/3 - expanded spcFormat/"
    * "Soundings/4 - pklFormat/"
    * "Soundings/5 - dgzStats/"

* unpickle2.py
  * Reads prepared soundings, calculates various ptype related parameters, and creates graphs to display these.
  * Utilizes:
    * EnergyProfile.py
    * sounding_functions.py
    * plotting.py
  * Reads data from "Soundings/4 - pklFormat/".
  * Creates plots and html pages under "plots/".

Classes and Helper Programs
---------------------------
* SPCSounding.py -- Class representing a basic SPC six-column sounding.

* spcProfile.py -- Contains the class spcExpProfile which represents an expanded version of an SPC profile with more columns.

* EnergyProfile.py -- Represents a vertical profile in terms of energy layers and related properties. 

* sounding_functions.py -- Calculates various properties of soundings.

* plotting.py -- Utilizes matplotlib and PtypeMarker.py to create various plots and the webpages to display them along with the skew-t diagram and representative METAR (if available).
                    
External Packages
-----------------
numpy
sharpy
matplotlib
