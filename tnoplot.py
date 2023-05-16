from astroquery.mpc import MPC
from astropy.table import Table
import numpy as np
from convmf import *
import matplotlib.pyplot as pl
from astropy.time import Time

# All transneptunian objects
result = MPC.query_objects('asteroid', semimajor_axis_min=30.1, tisserand_jupiter_min=3.05) 
#result = MPC.query_objects('asteroid', perihelion_distance_min=30.1) 

tnotab=Table(result)

ms=np.radians(tnotab['mean_anomaly'].data.astype(np.float))
es=tnotab['eccentricity'].data.astype(np.float)
omegas=np.radians(tnotab['ascending_node'].data.astype(np.float))
ws=np.radians(tnotab['argument_of_perihelion'].data.astype(np.float))
incs=np.radians(tnotab['inclination'].data.astype(np.float))
sas=tnotab['semimajor_axis'].data.astype(np.float)
epochjds=tnotab['epoch_jd'].data.astype(np.float)

# Shift to current epoch
dts=(Time.now().jd-epochjds)/365.25 # convert to years
mu=4*np.pi**2 						# mu = G(M+m), G=4pi**2 au**3/MSun/yr**2, M=1MSun, m << M
mmotion=np.sqrt(mu/sas**3) 			# mean motion
mcurrent=ms+mmotion*dts

fs=np.squeeze(convmf(mcurrent,es))
cfs   = np.cos(fs)
rs    = sas*(1-es**2)/(1+es*cfs)

# 3D
xs = rs*(np.cos(omegas)*np.cos(ws+fs)-np.sin(omegas)*np.sin(ws+fs)*np.cos(incs))
ys = rs*(np.sin(omegas)*np.cos(ws+fs)+np.cos(omegas)*np.sin(ws+fs)*np.cos(incs))
zs = rs*(np.sin(ws+fs)*np.sin(incs))

# Pluto's orbit
pid=np.where(tnotab['name']=='Pluto')[0]
morb=np.linspace(0,2*np.pi,num=100)
forb=np.squeeze(convmf(morb,es[pid]))
xorb = rs[pid]*(np.cos(omegas[pid])*np.cos(ws[pid]+forb)-np.sin(omegas[pid])*np.sin(ws[pid]+forb)*np.cos(incs[pid]))
yorb = rs[pid]*(np.sin(omegas[pid])*np.cos(ws[pid]+forb)+np.cos(omegas[pid])*np.sin(ws[pid]+forb)*np.cos(incs[pid]))

#fig,ax=pl.figure(figsize=[13.3,13.3],dpi=100,facecolor='k',edgecolor='k')
#pl.plot(xs,ys,'.',color='w')
#pl.ylim(-100,100)
#pl.xlim(-100,100)

fig1, ax1 = pl.subplots(figsize=[10.5,10.5],dpi=100,facecolor='k')
ax1.set_aspect('equal')
ax1.set_facecolor('k')
ax1.plot(xs,ys,'.',color='w',markersize=6)
ax1.plot(0,0,'y*',markersize=10)
ax1.plot(xorb,yorb,'xkcd:sky blue',lw=5)
#ax1.arrow(arrowco[j,0],arrowco[j,1],0,arrowco[j,3], fc='k', ec='k',head_width=3)
ax1.annotate('Pluto',xy=(xs[pid],ys[pid]),xytext=(10,-55),arrowprops=dict(facecolor='#75bbfd', shrink=0.05),c='#75bbfd')
pl.ylim(-100,100)
pl.xlim(-100,100)
pl.axis('off')
pl.savefig("tnosnow.svg",facecolor='k', bbox_inches='tight',dpi=300)

#fig1, ax1 = pl.subplots(figsize=[10.5,10.5],dpi=100,facecolor='k')
#ax1.set_aspect('equal')
#ax1.set_facecolor('k')
#ax1.plot(xs,ys,'.',color='w',markersize=2)
#ax1.plot(0,0,'y*',markersize=10)
##ax1.plot(xorb,yorb,'xkcd:sky blue',lw=5)
###ax1.arrow(arrowco[j,0],arrowco[j,1],0,arrowco[j,3], fc='k', ec='k',head_width=3)
##ax1.annotate('Pluto',xy=(xs[pid],ys[pid]),xytext=(10,-55),arrowprops=dict(facecolor='#75bbfd', shrink=0.05),c='#75bbfd')
#pl.ylim(-100,100)
#pl.xlim(-100,100)
#pl.axis('off')
#pl.savefig("tnos-nop2.png",facecolor='k', bbox_inches='tight',dpi=300)
