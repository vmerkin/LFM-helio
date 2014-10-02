###########################
masdir = '/Users/merkivg1/work/LFM-helio_2.0/PSI/Polytropic_run/tdmpoly01f'
lfmfile = '/Users/merkivg1/work/LFM-helio_2.0/PSI/results/CME/mas2lfm_tdmpoly01R_012_mhd_0006290.hdf'
#lfmfile = '/Users/merkivg1/GitHub/LFM-helio/new_relaxed/mas2lfm_106x96x128_tx2_gravity_mhd_0004320.hdf'
#lfmfile = '/Users/merkivg1/GitHub/LFM-helio/mas2lfm_106x192x256_mhd_0005560.hdf'

Rslice = 20.
equatorial_slice = True
variable = 'bt'
time_label = 251

gam = 1.05
###########################







import os,sys,glob,lfm
from pyhdf.SD import SD, SDC
import mas

from numpy import linspace,pi,meshgrid,sin,cos,zeros,ones,dstack,diff,sqrt,array,savetxt,flatnonzero,insert
from matplotlib import pyplot as plt



#############################
# init variable dictionaires
vars={}
vars['mas']={}
vars['lfm']={}
#############################





#############################
### MAR READ STUFF ###
#############################
mas_var_names = ['t','rho','vt','vp','vr','bt','bp','br']
mas_var_units = ['temperature','mass density',
                 'velocity','velocity','velocity',
                 'magnetic field','magnetic field','magnetic field']

for var_name,var_unit in zip(mas_var_names,mas_var_units):
    vars['mas'][var_name]={}
    (vars['mas'][var_name]['phi'],
     vars['mas'][var_name]['theta'],
     vars['mas'][var_name]['r'],
     vars['mas'][var_name]['data']) = mas.read_var(os.path.join(masdir,var_name+'%03d.hdf'%time_label),var_unit)
r_mas = vars['mas'][variable]['r']

vars['mas']['bt']['lims'] = -5.e-4,5.e-4
vars['mas']['bp']['lims'] = -5.e-4,5.e-4
vars['mas']['br']['lims'] = -1.e-3,1.e-3
vars['mas']['vr']['lims'] = 3.e7,4.e7
vars['mas']['vt']['lims'] = -2.e6,2.e6
vars['mas']['vp']['lims'] = -2.e6,2.e6
vars['mas']['rho']['lims'] = 2.e-21,2.e-20  #5.e-21,1.e-20
vars['mas']['t']['lims'] = 9.7e5,1.05e6

vars['mas']['bt']['fmt'] = '%.4f'
vars['mas']['bp']['fmt'] = '%.4f'
vars['mas']['br']['fmt'] = '%.4f'
vars['mas']['vr']['fmt'] = '%.4e'
vars['mas']['vt']['fmt'] = '%.4e'
vars['mas']['vp']['fmt'] = '%.4e'
vars['mas']['rho']['fmt'] = '%.4e'
vars['mas']['t']['fmt'] = '%.4e'
#############################






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
theta=arccos(z/R)
phi=arctan2(y,x)
phi[phi<0]+=2*pi


mas_islice = flatnonzero(vars['mas'][variable]['r']>=Rslice)[0]
lfm_islice = flatnonzero(R[0,0,:]>=Rslice)[0]
if not equatorial_slice:
    lfm_selection = (slice(None),slice(None),lfm_islice)
    mas_selection = (slice(None),slice(None),mas_islice)
else:
    mas_jslice = flatnonzero(vars['mas'][variable]['theta']>=pi/2.)[0]
    lfm_jslice = flatnonzero(theta[0,:,0]>=pi/2.)[0]
    lfm_selection = (slice(None),lfm_jslice,slice(lfm_islice,None))
    mas_selection = (slice(None),mas_jslice,slice(mas_islice,None))
    


br     = bx[lfm_selection]*cos(phi[lfm_selection])*sin(theta[lfm_selection]) + by[lfm_selection]*sin(phi[lfm_selection])*sin(theta[lfm_selection]) + bz[lfm_selection]*cos(theta[lfm_selection])
btheta = bx[lfm_selection]*cos(phi[lfm_selection])*cos(theta[lfm_selection]) + by[lfm_selection]*sin(phi[lfm_selection])*cos(theta[lfm_selection]) - bz[lfm_selection]*sin(theta[lfm_selection])
bphi   =-bx[lfm_selection]*sin(phi[lfm_selection])            + by[lfm_selection]*cos(phi[lfm_selection])

vr     = vx[lfm_selection]*cos(phi[lfm_selection])*sin(theta[lfm_selection]) + vy[lfm_selection]*sin(phi[lfm_selection])*sin(theta[lfm_selection]) + vz[lfm_selection]*cos(theta[lfm_selection])
vtheta = vx[lfm_selection]*cos(phi[lfm_selection])*cos(theta[lfm_selection]) + vy[lfm_selection]*sin(phi[lfm_selection])*cos(theta[lfm_selection]) - vz[lfm_selection]*sin(theta[lfm_selection])
vphi   =-vx[lfm_selection]*sin(phi[lfm_selection])            + vy[lfm_selection]*cos(phi[lfm_selection])


(vars['lfm']['bp'],
 vars['lfm']['bt'],
 vars['lfm']['br'],
 vars['lfm']['vp'],
 vars['lfm']['vt'],
 vars['lfm']['vr'],
 vars['lfm']['rho'],
 vars['lfm']['c']) = (bphi,btheta,br,vphi,vtheta,vr,rho[lfm_selection],c[lfm_selection])

vars['lfm']['t'] =  vars['lfm']['c']**2/gam*1.67e-8/1.38
#############################





#############################
###### PLOTTING ######
#############################
p_mas,t_mas = vars['mas'][variable]['phi'],vars['mas'][variable]['theta']

mas_x = p_mas
lfm_x = phi[lfm_selection]

if not equatorial_slice:
    mas_y = pi-t_mas
    lfm_y = pi-theta[lfm_selection]
else:
    mas_y = r_mas[mas_islice:]
    lfm_y = R[lfm_selection]


fig = plt.figure(figsize = (18,16))

vmin,vmax = vars['mas'][variable]['lims']
fmt = vars['mas'][variable]['fmt']

if not equatorial_slice:
    ax1 = plt.subplot(211,aspect='equal')
    ax2 = plt.subplot(212,aspect='equal')
else:
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)

p1=ax1.pcolormesh(mas_x,mas_y,vars['mas'][variable]['data'][mas_selection].T,vmin=vmin,vmax=vmax)
ax1.xaxis.set_visible(False)
ax1.set_title( ('MAS: Min/Max = '+fmt+'/'+fmt) % 
               (vars['mas'][variable]['data'][mas_selection].min(),vars['mas'][variable]['data'][mas_selection].max()))

p2 = ax2.pcolormesh(lfm_x,lfm_y,vars['lfm'][variable],vmin=vmin,vmax=vmax)
ax2.set_xlabel(r'$\phi$')
ax2.set_title( ('LFM: Min/Max = '+fmt+'/'+fmt) % 
               (vars['lfm'][variable].min(),vars['lfm'][variable].max()))



for ax in [ax1,ax2]:
    ax.set_xlim(0,2*pi)
    if not equatorial_slice:
        ax.set_ylim(0,pi)
        ax.set_ylabel(r'$\theta$')
    else:
        ax.set_ylim(Rslice,r_mas.max())
        ax.set_ylabel('R')
    plt.colorbar(p1,ax=ax).set_label(variable)

plt.show()

    
