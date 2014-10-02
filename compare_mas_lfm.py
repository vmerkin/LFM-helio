###########################
masdir = '/Users/merkivg1/work/LFM-helio_2.0/PSI/Polytropic_run/tdmpoly01R'
#lfmfile = './mas2lfm_106x96x128_mhd_0003000.hdf'
lfmfile = '/Users/merkivg1/work/LFM-helio_2.0/PSI/results/CME/mas2lfm_tdmpoly01R_012_mhd_0003000.hdf'
Rslice = 30.
variable = 'br'
time_label = 12
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

#lfm_islice = flatnonzero(R[0,0,:]>=Rslice)[0]

theta=arccos(z[0,:,:]/R[0,:,:])
phi=arctan2(y[0,:,:],x[0,:,:])
phi[phi<0]+=2*pi

br     = bx[0,:,:]*cos(phi)*sin(theta) + by[0,:,:]*sin(phi)*sin(theta) + bz[0,:,:]*cos(theta)
btheta = bx[0,:,:]*cos(phi)*cos(theta) + by[0,:,:]*sin(phi)*cos(theta) - bz[0,:,:]*sin(theta)
bphi   =-bx[0,:,:]*sin(phi)            + by[0,:,:]*cos(phi)

vr     = vx[0,:,:]*cos(phi)*sin(theta) + vy[0,:,:]*sin(phi)*sin(theta) + vz[0,:,:]*cos(theta)
vphi   =-vx[0,:,:]*sin(phi)            + vy[0,:,:]*cos(phi)


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
vars['mas']['br']['lims'] = -1.e-3,1.e-3
vars['mas']['vr']['lims'] = 3.e7,4.e7
vars['mas']['vt']['lims'] = -2.e6,2.e6
vars['mas']['vp']['lims'] = -2.e6,2.e6
vars['mas']['rho']['lims'] = 5.e-21,1.e-20
vars['mas']['t']['lims'] = 9.7e5,1.05e6

vars['mas']['bt']['fmt'] = '%.8f'
vars['mas']['bp']['fmt'] = '%.8f'
vars['mas']['br']['fmt'] = '%.4f'
vars['mas']['vr']['fmt'] = '%.4e'
vars['mas']['vt']['fmt'] = '%.4e'
vars['mas']['vp']['fmt'] = '%.4e'
vars['mas']['rho']['fmt'] = '%.4e'
vars['mas']['t']['fmt'] = '%.4e'

#mas_islice = flatnonzero(vars['mas'][variable]['r']>=Rslice)[0]
#############################

# interpolate MAS onto LFM grid
f=interpolate.RectBivariateSpline(vars['mas'][variable]['theta'],
                                  vars['mas'][variable]['r']/6.96e10,
                                  vars['mas'][variable]['data'][0,:,:],kx=1,ky=1)
mas_interp=f(theta[:,0],R[0,0,:])

#############################
###### PLOTTING ######
#############################
fig = plt.figure(figsize = (18,16))

vmin,vmax = vars['mas'][variable]['lims']
fmt = vars['mas'][variable]['fmt']

ax1 = plt.subplot(311)

p1=ax1.pcolormesh(R[0,:,:],pi-theta,mas_interp,vmin=vmin,vmax=vmax)
ax1.xaxis.set_visible(False)
ax1.set_title( ('MAS: Min/Max = '+fmt+'/'+fmt) % 
               (mas_interp.min(),mas_interp.max()))

ax2 = plt.subplot(312)
p2 = ax2.pcolormesh(R[0,:,:],pi-theta,vars['lfm'][variable],vmin=vmin,vmax=vmax)
ax2.xaxis.set_visible(False)
ax2.set_title( ('LFM: Min/Max = '+fmt+'/'+fmt) % 
               (vars['lfm'][variable].min(),vars['lfm'][variable].max()))

ax3 = plt.subplot(313)
p3 = ax3.pcolormesh(R[0,:,:],pi-theta,vars['lfm'][variable]/mas_interp-1,vmin=-0.05,vmax=0.05,cmap=cm.RdBu_r)
ax3.set_xlabel('R')

ax3.set_title( ('Ratio-1: Min/Max = '+fmt+'/'+fmt) % 
               ((vars['lfm'][variable]/mas_interp-1).min(),(vars['lfm'][variable]/mas_interp-1).max()))



for (ax,p) in zip([ax1,ax2,ax3],[p1,p2,p3]):
    ax.set_xlim(20,50)
    ax.set_ylim(0,pi)
    ax.set_ylabel(r'$\theta$')
    plt.colorbar(p,ax=ax).set_label(variable)




plt.show()

    
