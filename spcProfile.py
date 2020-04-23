# -*- coding: utf-8 -*-
"""
Created on Thu Oct 05 10:03:14 2017

@author: eric.lenning
"""
#import sharppy.sharptab.profile as profile
import sharppy.sharptab.utils as utils
import getpass
from datetime import datetime

class spcExpProfile(object):

    def __init__(self, prof):
        self._prof = prof # using composition
        
    def toExpandedFile(self, snd_file):
        #snd_file = open(file_name, 'w')
        def qc(val):
            return -9999. if not utils.QC(val) else val

        snd_loc = (" " * (4 - len(self._prof.location))) + self._prof.location

        now = datetime.utcnow()
        user = getpass.getuser()
        snd_file.write("%TITLE%\n")
        snd_file.write("%s   %s\n Saved by user: %s on %s UTC\n" % (snd_loc, self._prof.date.strftime("%y%m%d/%H%M"), user, now.strftime('%Y%m%d/%H%M')))
        snd_file.write("   LEVEL       HGHT       TEMP       DWPT       WDIR       WSPD       WBLB       THTE       LR\n")
        snd_file.write("--------------------------------------------------------------------------------------------------\n")
        snd_file.write("%RAW%\n")
        for idx in xrange(self._prof.pres.shape[0]):
            str = ""
            for col in ['pres', 'hght', 'tmpc', 'dwpc', 'wdir', 'wspd', 'wetbulb','thetae']:
                str += "%8.2f,  " % qc(self._prof.__dict__[col][idx])
                
            if idx > 0:
                tdif = self._prof.__dict__['tmpc'][idx] - self._prof.__dict__['tmpc'][idx-1]
                hdif = self._prof.__dict__['hght'][idx] - self._prof.__dict__['hght'][idx-1]
                lr = 1000 * tdif/hdif
                
                if lr < -9.8:
                    str += "%8.2f *,  " % lr
                else:
                    str += "%8.2f,  " % lr
                
            snd_file.write(str[:-3] + "\n")
        snd_file.write("%END%\n")
        #snd_file.close()        