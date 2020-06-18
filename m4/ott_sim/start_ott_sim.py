# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 10:06:35 2020

@author: Runa
list of startup instruction to startup the simulated OTT and get interferometer frames
usage:
    import everything as follows
    setup the OTT:
        ott.slide, ott.rslide, ott.angle
            slide=0.8 for segment, rslide =0.3 for alignment, angle 0,60,120... for segments
    setup the mirrors in ott:
        ott.param, ott.m4, ott.refflat
        
    acquire an image 
"""


from importlib import reload
import os
from matplotlib import pyplot as plt
a='D:\Astro\ARCETRI\Python\M4-master'
os.chdir(a)
from m4.configuration.create_ott import *
from m4.configuration import start
ott=start.create_ott()
from m4.ott_sim.ott_images import *
from m4.utils.Interface_4d import comm4d
#test
plt.imshow(ott.mask)



ott.slide(0.75)
ott.angle(30.)
ott.rslide(0.3)
ott.parab([0,0,1e-5,0,0,0])
ott.m4([0,0,0,1e-6,0,0])
ott.refflat([0,0,0,0,1e-6,0])


p,m = comm4d.acq4d(1,show=1)
#which is a call to the simulated interf. acquisition
p,m=ott_smap(show=1)



