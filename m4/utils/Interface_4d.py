# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 06:20:03 2020

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
from m4.configuration import config as conf
from m4.configuration.ott_parameters import *


class comm4d():
    
    phcam_ip = ''
    phcam_datapath = ''
    phcam_conf = ''
    
    def acq4d(nframes = 1, show=0):
        
        if conf.simulated ==1:
            opd,mask = ott_smap(show=show)
            
        else:
            print('some function to acquire the interferometer....')
            opd = np.array(OTT_Parameters.N_PIXEL)
            mask = opd
        
        return(opd, mask)