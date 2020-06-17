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
import matplotlib as plt

def ott_smap(offset=None, quant=None):
    npix = Interferometer.N_PIXEL
    pmap    = ott_parab_ima()
    m4map,m4smask   =  ott_m4_ima()
    mmask = m4smask*ott.mask
    smap    = (pmap + m4map) * mmask
    if offset != None:
        smap = smap+offset 
    
    #fullmask  = maskm    
    roffset = (ott.rslide()-ott.slide())*OttParameters.pscale
    draw_ref = (abs(roffset)-OttParameters.rflat_radius*OttParameters.pscale <  npix[1]/2)
    if draw_ref ==1:
        rmap,rmask                = ott_rflat_ima()
        #smap[ott.idx[idr]]  = rmap[ott.idx[idr]]+pmap[ott.idx[idr]]
        #smap[idc]           = 0
        #fullmask            = maskm+maskr
        #fullmask[np.where(fullmask)]  = 1
        #fullmask[idc]                   = 0
   
   
    rmap = (rmap + pmap)*rmask
    romask = geo.draw_mask(ott.mask*0, npix[0]/2,npix[1]/2+roffset,OttParameters.rflat_radius*OttParameters.pscale-OttParameters.rflat_cell*OttParameters.pscale)
    romask = romask+rmask
    irmask = -rmask+1
    
    smap1 = smap*irmask+rmap
    smap1 = smap1 + ott.offset  #this is to be verified
    smask = rmap +mmask
    smask[smask != 0] =1
    smask[romask == 1] =0
    smap1 = smap1*smask
    if quant ==1:
        print('add quantization')
       #ToDo smap  = add_quantization(smap, bad=fullmask, quant=ott.interf.quantization):
           
    #mask  = fullmask
    #cir   = qpupil(fullmask, xx=xx, yy=yy,nocircle=1)
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
    plt.imshow(smap1, cmap='hot');plt.colorbar()    
    return(smap1, smask)
    
def pwrap(img):
    wav = Interferometer.WAVEL
    img1 = img.copy()
    optfact = 1
    img1[img1 != 0] = np.sin(2*!pi*img1[img1 != 0]*optfact/(wav))
    return(img1)    
    
def ott_parab_ima():
    npix = Interferometer.N_PIXEL
    smap    = ott.smap.copy()
    ww      = np.dot(ott.zmat(), ott.zmx_parpos2z()) 
    
    for i in range(0,5):
        smap[ott.parmask == True] = smap[ott.parmask == True] + ww[:,i]* (ott.parab())[i]
    
    mask = geo.draw_mask(ott.mask,npix[0]/2,npix[1]/2,ott.fold_radius*OttParameters.pscale)  
    smap = smap * mask
    return(smap)
    
    
    
def ott_rflat_ima( deshape=0):
    npix = Interferometer.N_PIXEL
    pscale = OttParameters.pscale
    rmap    = ott.smap
    roffset = (ott.rslide()-ott.slide())*pscale
    rflat_radius = OttParameters.rflat_radius
    if abs(roffset)-rflat_radius*pscale <  npix[1]/2 :
        rmask = np.ones([npix[0],npix[1]])
        rmask = geo.draw_mask(ott.mask*0, npix[0]/2,npix[1]/2+roffset,rflat_radius*pscale-1)
        if deshape==1:
            #rima = read_refflat_deformation( 1.)
            rmap[npix[0]/2-rflat_radius*pscale-1,npix[1]/2-rflat_radius*pscale-1]= rima 
        
        #idr   = where(rmask[ott.idx] )
        ww      = np.dot(ott.zmat(), ott.zmx_refflatpos2z()) 
        for i in range(0,5):
            rmap[ott.parmask == True] = rmap[ott.parmask == True] + ww[:,i]*(ott.refflat())[i]
            #rmap[ott.idx[idr]] += ww[i,idr]* ott.refflat[i]
            
        rmap = rmap*rmask
        #cm      = edge_erode(erode=2, rmask)
        #idc     = where((cm+rmask) eq 1)
        
    return( rmap, rmask  )
    
    
    
def ott_m4_ima(): #bad=mask, idm=idm, idtot=idtot, optpupil=optpupil

    theta   = ott.angle()*np.pi/180.
    rmat    = [[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]]
    ss      = ott.m4pupil.shape
  
    m4ima = ott.m4ima + ott.m4offset
    usepupil = ott.m4pupil.copy()
    #if keyword_set(optpupil) then usepupil = usepupil*ott.m4optpupil   ;mod 20161010
  
    m4smask   = ott_map2ima(usepupil)
    m4img     = ott_map2ima(m4ima)
    mask      = m4smask*ott.mask
    #idtot     = where(mask[ott.idx])
    m4smap  = ott.smap.copy()
    ww      = np.dot(ott.zmat(), ott.zmx_m4pos2z())
    for i in range(0,5): 
        m4smap[ott.parmask == True] = m4smap[ott.parmask == True] + ww[:,i]* (ott.m4())[i]
        #m4smap[ott.idx[idtot]] = m4smap[ott.idx[idtot]] + ww[i,idtot]* ott.m4[i]
    
    #m4img[ott.idx[idtot]] += m4smap[ott.idx[idtot]]
    m4img = (m4img + m4smap)*m4smask
    return( m4img, m4smask)
  

def ott_map2ima(w):   #debugged
    npix =  Interferometer.N_PIXEL
    theta   = ott.angle()*np.pi/180.
    rmat    = [[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]]
    ss      = w.shape
    theangle = -30.-ott.angle()
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