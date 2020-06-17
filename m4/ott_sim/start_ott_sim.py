# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 10:06:35 2020

@author: Runa
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
#test
plt.imshow(ott.mask)


ott.slide(0.75)
ott.angle(30.)
ott.rslide(0.3)
ott.parab([0,0,1e-5,0,0,0])
ott.m4([0,0,0,1e-6,0,0])
ott.refflat([0,0,0,0,1e-6,0])

p,m=ott_smap()



