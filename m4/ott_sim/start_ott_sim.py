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

#test
plt.imshow(ott.mask)
