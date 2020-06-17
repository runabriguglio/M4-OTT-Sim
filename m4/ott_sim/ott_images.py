# -*- coding: utf-8 -*-
"""
Created on Thu May 21 18:06:09 2020

@author: Runa
"""
"""

"""
import numpy as np
from m4.configuration import start
ott = start.create_ott()
from m4.configuration.ott_parameters import *
import geo

def ott_smap():
    npix = Interferometer.N_PIXEL
    pmap    = ott_parab_ima()
    m4map   =  ott_m4_ima1(bad=maskm, idm=idm, idtot=idtot, optpupil=optpupil)
    smap    = (pmap + m4map) * maskm
    if offset != None:
        smap = smap+offset 
    
    fullmask  = maskm    
    roffset = (ott_rslide()-ott_slide())*OttParameters.pscale
    draw_ref = (abs(roffset)-ott.rflat_radius*OttParameters.pscale <  npix[1]/2)
    if draw_ref ==1:
        rmap                = ott_rflat_ima(bad=maskr, idr=idr, idc=idc,deshape=rflatshape)
        smap[ott.idx[idr]]  = rmap[ott.idx[idr]]+pmap[ott.idx[idr]]
        smap[idc]           = 0
        fullmask            = maskm+maskr
        fullmask[np.where(fullmask)]  = 1
        fullmask[idc]                   = 0
   
   
    smap = smap + ott.offset  
   
    if quant ==1:
        print('add quantization')
       #ToDo smap  = add_quantization(smap, bad=fullmask, quant=ott.interf.quantization):
           
    mask  = fullmask
    cir   = qpupil(fullmask, xx=xx, yy=yy,nocircle=1)
    """ 
    ToDo
    if keyword_set(show) then begin
      smap1 = smap
      
      out = where(fullmask eq 0)
      smap1[out] = max(smap1[where(fullmask)])
      wset, 0
      image_show, /As,  smap1;, titl='WF (?) map'
      wset, 1
      gg = ott_geometry(/sh)
      wset, 3
      display_data,  mirror.gpos[*], mirror.coord, spot_mag=2, /no_n,/Sh, back=max(mirror.gpos);image_Show, /As, /sh, mirror.m4ima
      if keyword_set(fringes) then begin
        wset, 2
        loadct, 0,/silent
        interf = pwrap(smap, lambda=ott.interf.lambda, optfact=1, bad=mask, detector_mask=detmask)*fullmask
        
        out = where(fullmask eq 0)
        
        interf[0:ott.interf.npix[0]/5,4.5*ott.interf.npix[1]/5:ott.interf.npix[1]-1]=-1
        interf[4*ott.interf.npix[0]/5:ott.interf.npix[0]-1,4.5*ott.interf.npix[1]/5:ott.interf.npix[1]-1]=1
       
        image_show, /As, interf, min_v=-1, max_v=1
        loadct, 3,/silent
      endif
  endif
"""   
    return(smap)
    
    
def ott_parab_ima():
    npix = Interferometer.N_PIXEL
    smap    = ott.smap.copy
    ww      = numpy.dot(ott.zmat(), ott.zmx_parpos2z()) 
    for i in range(0,5):
        smap[ott.idx] = smap[ott.idx] + ww[i,:]* ott.parab[i]
        
    mask = geo.draw_mask(ott.mask(),npix[0]/2,npix[1]/2,ott.fold_radius*OttParameters.pscale)  
    smap = smap * mask
    return(smap)
    
def ott_m4_ima1(): #bad=mask, idm=idm, idtot=idtot, optpupil=optpupil

    theta   = ott.angle()*np.pi/180.
    rmat    = [[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]]
    ss      = ott.m4pupil.shape
  
    mirror.m4ima = mirror.m4ima + ott.m4offset
    usepupil = ott.m4pupil
    #if keyword_set(optpupil) then usepupil = usepupil*ott.m4optpupil   ;mod 20161010
  
    m4smask   = self._ott_map2ima(usepupil)
    m4img     = self._ott_map2ima(mirror.m4ima)
    mask      = m4smask*ott.mask
    idtot     = where(mask[ott.idx])
    m4smap  = ott.smap
    ww      = np.dot(ott.zmat(), ott.zmx_m4pos2z)
    for i in range(0,5): 
        m4smap[ott.idx[idtot]] = m4smap[ott.idx[idtot]] + ww[i,idtot]* ott.m4[i]
    
    m4img[ott.idx[idtot]] += m4smap[ott.idx[idtot]]
    return( m4img)
  

def ott_map2ima(w):   #debugged
    npix =  Interferometer.N_PIXEL
    theta   = ott.angle()*np.pi/180.
    rmat    = [[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]]
    ss      = w.shape
    theangle = -90.-ott.angle()
    simg    = geo.rotate(w, theangle)
    parxy   = ott.slide() * OttParameters.pscale
    x0      = np.fix(ss[0]/2-npix[0]/2)
    x1      = np.fix(ss[0]/2+npix[0]/2-1)
    y0      = np.fix(ss[1]/2+parxy-npix[1]/2)
    y1      = np.fix(ss[1]/2+parxy+npix[1]/2-1)
    x0=x0.astype(int)
    x1=x1.astype(int)
    y0=y0.astype(int)
    y1=y1.astype(int)
    simg    = simg[x0:x1+1, y0:y1+1]
    return(simg)