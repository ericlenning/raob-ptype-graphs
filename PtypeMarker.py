# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 11:51:23 2017

@author: eric.lenning
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath

class PtypeMarker(object):
    
    def __init__(self, **kwargs):
        self.ra = self.get_ra()
        self.sn = self.get_sn()
        self.pl = self.get_pl()
        self.zr = self.get_zr()
        
        self.razr = self.ra.make_compound_path(self.ra, self.zr)
        self.rasn = self.ra.make_compound_path(self.ra, self.sn)
        self.rapl = self.ra.make_compound_path(self.ra, self.pl)
        self.zrpl = self.ra.make_compound_path(self.zr, self.pl)
        self.zrsn = self.ra.make_compound_path(self.zr, self.sn)
        self.plsn = self.ra.make_compound_path(self.pl, self.sn)
        
        self.raplsn = self.ra.make_compound_path(self.ra, self.pl, self.sn)
        self.razrsn = self.ra.make_compound_path(self.ra, self.zr, self.sn)
        self.zrplsn = self.ra.make_compound_path(self.pl, self.zr, self.sn)
        self.razrpl = self.ra.make_compound_path(self.ra, self.pl, self.zr)
        
        self.razrplsn = self.ra.make_compound_path(self.ra, self.pl, self.zr, self.sn)
        
    def get_ra(self):
        u = np.array([  [2.0,2.0],
                        [-2.0,2.0],
                        [-2.0,-2.0]  ])
        
        codes = [1] + [2]*(len(u)-2) + [2] 
    
        return mpath.Path(u, codes, closed=False)
    
    def get_sn(self):
        u = np.array([  [-2.0,-2.0],
                        [2.0,-2.0],
                        [2.0,2.0]  ])
        
        codes = [1] + [2]*(len(u)-2) + [2] 
    
        return mpath.Path(u, codes, closed=False)
    
    def get_pl(self):
        u = np.array([  [2,-0.5],
                        [-2,-0.5],
                        [0.5,2],
                        [0.5,-2]  ])
        
        codes = [1, 2, 1, 2] 
    
        return mpath.Path(u, codes, closed=False)
    
    def get_zr(self):
        u = np.array([  [-0.5,2.0],
                        [-0.5,-2.0],
                        [-2.0,0.5],
                        [2.0,0.5]  ])
        
        codes = [1, 2, 1, 2] 
    
        return mpath.Path(u, codes, closed=False)

    def demo(self):
        
        pm = PtypeMarker()
        
        print pm.ra, pm.sn, pm.pl, pm.zr, pm.razr
        
        size = 70
        
        plt.scatter([2,2,3],[1.4,2.3,2.8], s=size, marker=pm.ra, 
                    edgecolors="crimson", facecolors='none', linewidth=1)
        plt.scatter([2,2,3],[1.4,2.3,2.8], s=size, marker=pm.sn, 
                    edgecolors="blue", facecolors='none', linewidth=1)
        plt.scatter([1,1,2],[1.4,2.3,2.8], s=size, marker=pm.zr, 
                    edgecolors="crimson", facecolors='none', linewidth=2)
        
        #plt.scatter([1,1,2],[1.4,2.3,2.8], s=size, marker=pm.ra, 
        #            edgecolors="crimson", facecolors='none') # , linewidth=2)
        plt.scatter([0,1,2],[1,3,1], s=size, marker=pm.sn, 
                    edgecolors="k", facecolors='none')
        
        plt.scatter([0,2,5],[2,3,2], s=size, marker=pm.sn, 
                    edgecolors="k", facecolors='none')
        plt.scatter([0,2,5],[2,3,2], s=size, marker=pm.pl, 
                    edgecolors="k", facecolors='none', linewidth=0)
        plt.scatter([0,2,5],[2,3,2], s=size, marker=pm.zr, 
                    edgecolors="k", facecolors='none', linewidth=2)
        
        plt.scatter([1,3,5],[5,3,1], s=size, marker=pm.zr, 
                    edgecolors="k", facecolors='none')
        plt.scatter([2,4,6],[6,4,2], s=size, marker=pm.razr, 
                    edgecolors="k", facecolors='none')
        plt.scatter([2.7,4.7,6.7],[6.7,4.7,2.7], s=size, marker=pm.rapl, 
                    edgecolors="olive", facecolors='none')
        plt.scatter([2,3,4],[2,2,2], s=size, marker=pm.rasn, 
                    edgecolors="b", facecolors='none')
        plt.scatter([1,2,3],[3,2,1], s=size, marker=pm.zrpl, 
                    edgecolors="g", facecolors='none')
        #plt.scatter([0,1.8,3],[0,2,4], s=size, marker=pm."o", 
        #            edgecolors="k", facecolors='none')
        plt.scatter([0.5,2.8,3.4],[0,2,4], s=size, marker=pm.zrsn, 
                    edgecolors="m", facecolors='none')
        plt.scatter([0.5,2.8,3.4],[0.6,2.6,4.9], s=size, marker=pm.plsn, 
                    edgecolors="y", facecolors='none')
        plt.scatter([2.5,4.8,5.4],[0.6,2.6,4.9], s=size, marker=pm.raplsn, 
                    edgecolors="brown", facecolors='none')
        plt.scatter([2.5,4.8,5.4],[1.6,3.6,5.9], s=size, marker=pm.razrsn, 
                    edgecolors="grey", facecolors='none')
        plt.scatter([2.1,4.1,5.1],[1.1,3.1,5.1], s=size, marker=pm.zrplsn, 
                    edgecolors="darkcyan", facecolors='none')
        plt.scatter([5.1,3.1,2.1],[1.1,3.1,5.1], s=size, marker=pm.razrpl, 
                    edgecolors="orange", facecolors='none')
        plt.scatter([5.1,3.1,2.1],[2.1,2.1,2.1], s=size, marker=pm.razrplsn, 
                    edgecolors="violet", facecolors='none')
        
        plt.show()
        
        for i in 't','v','a':
            print i
            
if __name__ == "__main__":
        
    PtypeMarker().demo()