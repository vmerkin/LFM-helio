###########################
masdir = '/Users/merkivg1/work/LFM-helio_2.0/PSI/Polytropic_run/tdmpoly01f'
#masdir = '/Users/merkivg1/work/LFM-helio_2.0/PSI/parker01'
#lfmfile = './mas2lfm_106x96x128_mhd_0003000.hdf'
#lfmfile = '/Users/merkivg1/work/LFM-helio_2.0/PSI/results/CME/mas2lfm_tdmpoly01R_012_mhd_0003000.hdf'
#lfmfile = './mas2lfm_parker_106x96x128_tx2_mhd_0003000.hdf'
#lfmfile = './mas2lfm_parker_106x96x128_tx2_gravity_mhd_0003000.hdf'
#lfmfile = './new_relaxed/mas2lfm_106x96x128_tx2_gravity_mhd_0003000.hdf'
lfmfile = './mas2lfm_106x192x256_mhd_0005560.hdf'
#lfmfile = './mas2lfm_106x96x128_tx2_gravity_mhd_0004320.hdf'
lfmfile = '/Users/merkivg1/work/LFM-helio_2.0/PSI/mas2lfm_helio_g105.hdf'
kslice=60
variable = 'br'
time_label = 251
#time_label = 7
gam = 1.05
###########################




import os,sys,glob,lfm
from pyhdf.SD import SD, SDC
import mas

from numpy import linspace,pi,meshgrid,sin,cos,zeros,ones,dstack,diff,sqrt,array,savetxt,flatnonzero,insert
from matplotlib import pyplot as plt
from matplotlib import cm as cm
from scipy import interpolate




#############################
# init variable dictionaires
vars={}
vars['lfm']={}

#############################
### LFM READ STUFF ####
#############################
from pyhdf.SD import SD, SDC
from pylab import * #sqrt,pcolor,figure
hdffile = SD(lfmfile,SDC.READ)
x=hdffile.select('X_grid').get()/6.96e10
y=hdffile.select('Y_grid').get()/6.96e10
z=hdffile.select('Z_grid').get()/6.96e10
bx=hdffile.select('bx_').get()[:-1,:-1,:-1]
by=hdffile.select('by_').get()[:-1,:-1,:-1]
bz=hdffile.select('bz_').get()[:-1,:-1,:-1]
vx=hdffile.select('vx_').get()[:-1,:-1,:-1]
vy=hdffile.select('vy_').get()[:-1,:-1,:-1]
vz=hdffile.select('vz_').get()[:-1,:-1,:-1]
rho=hdffile.select('rho_').get()[:-1,:-1,:-1]
c=hdffile.select('c_').get()[:-1,:-1,:-1]
hdffile.end()

# =========== Cell centers ==============
x=0.125*( x[:-1,:-1,:-1]+x[1:,:-1,:-1]+x[:-1,1:,:-1]+x[:-1,:-1,1:]+
          x[1:,1:,:-1]+x[1:,:-1,1:]+x[:-1,1:,1:]+x[1:,1:,1:] )
y=0.125*( y[:-1,:-1,:-1]+y[1:,:-1,:-1]+y[:-1,1:,:-1]+y[:-1,:-1,1:]+
          y[1:,1:,:-1]+y[1:,:-1,1:]+y[:-1,1:,1:]+y[1:,1:,1:] )
z=0.125*( z[:-1,:-1,:-1]+z[1:,:-1,:-1]+z[:-1,1:,:-1]+z[:-1,:-1,1:]+
          z[1:,1:,:-1]+z[1:,:-1,1:]+z[:-1,1:,1:]+z[1:,1:,1:] )
# =======================================

R=sqrt(x**2+y**2+z**2)

theta=arccos(z[kslice,:,:]/R[kslice,:,:])
phi=arctan2(y[kslice,:,:],x[kslice,:,:])
phi[phi<0]+=2*pi

br     = bx[kslice,:,:]*cos(phi)*sin(theta) + by[kslice,:,:]*sin(phi)*sin(theta) + bz[kslice,:,:]*cos(theta)
btheta = bx[kslice,:,:]*cos(phi)*cos(theta) + by[kslice,:,:]*sin(phi)*cos(theta) - bz[kslice,:,:]*sin(theta)
bphi   =-bx[kslice,:,:]*sin(phi)            + by[kslice,:,:]*cos(phi)

vr     = vx[kslice,:,:]*cos(phi)*sin(theta) + vy[kslice,:,:]*sin(phi)*sin(theta) + vz[kslice,:,:]*cos(theta)
vtheta = vx[kslice,:,:]*cos(phi)*cos(theta) + vy[kslice,:,:]*sin(phi)*cos(theta) - vz[kslice,:,:]*sin(theta)
vphi   =-vx[kslice,:,:]*sin(phi)            + vy[kslice,:,:]*cos(phi)


(vars['lfm']['bp'],
 vars['lfm']['bt'],
 vars['lfm']['br'],
 vars['lfm']['vp'],
 vars['lfm']['vt'],
 vars['lfm']['vr'],
 vars['lfm']['rho'],
 vars['lfm']['c']) = (bphi,btheta,br,vphi,vtheta,vr,rho[kslice,:,:],c[kslice,:,:])

vars['lfm']['t'] =  vars['lfm']['c']**2/gam*1.67e-8/1.38
#############################


#############################
vars['mas'] = mas.read_all_vars(masdir,time_label)
mas.set_plot_limits(vars['mas'],'ConfigScripts/plot.config')
#############################

nk=x.shape[0]
mas_kslice=where(vars['mas'][variable]['phi']>=float(kslice)/nk*2.*pi)[0][0]

# interpolate MAS onto LFM grid
f=interpolate.RectBivariateSpline(vars['mas'][variable]['theta'],
                                  vars['mas'][variable]['r']/6.96e10,
                                  vars['mas'][variable]['data'][mas_kslice,:,:],kx=1,ky=1)
mas_interp=f(theta[:,0],R[mas_kslice,0,:])

#############################
###### PLOTTING ######
#############################
fig = plt.figure(figsize = (18,16))

vmin,vmax = vars['mas'][variable]['lims']
fmt = vars['mas'][variable]['fmt']

ax1 = plt.subplot(311)

p1=ax1.pcolormesh(R[kslice,:,:],pi-theta,mas_interp,vmin=vmin,vmax=vmax)
ax1.xaxis.set_visible(False)
ax1.set_title( ('MAS: Min/Max = '+fmt+'/'+fmt) % 
               (mas_interp.min(),mas_interp.max()))

ax2 = plt.subplot(312)
p2 = ax2.pcolormesh(R[kslice,:,:],pi-theta,vars['lfm'][variable],vmin=vmin,vmax=vmax)
ax2.xaxis.set_visible(False)
ax2.set_title( ('LFM: Min/Max = '+fmt+'/'+fmt) % 
               (vars['lfm'][variable].min(),vars['lfm'][variable].max()))

ax3 = plt.subplot(313)
p3 = ax3.pcolormesh(R[kslice,:,:],pi-theta,vars['lfm'][variable]/mas_interp-1,vmin=-0.05,vmax=0.05,cmap=cm.RdBu_r)
ax3.set_xlabel('R')

ax3.set_title( ('Ratio-1: Min/Max = '+fmt+'/'+fmt) % 
               ((vars['lfm'][variable]/mas_interp-1).min(),(vars['lfm'][variable]/mas_interp-1).max()))



for (ax,p) in zip([ax1,ax2,ax3],[p1,p2,p3]):
    ax.set_xlim(20,50)
    ax.set_ylim(0,pi)
    ax.set_ylabel(r'$\theta$')
    plt.colorbar(p,ax=ax).set_label(variable)




plt.show()

    
