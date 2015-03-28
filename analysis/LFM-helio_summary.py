#!/usr/bin/python
import sys,os
from numpy import argmin,abs,pi,log10
import matplotlib
import matplotlib.pyplot as plt
home = os.path.expanduser("~");
if os.path.join(home,'GitHub/Python/LFM-helio/lib/') not in sys.path: sys.path.append(os.path.join(home,'GitHub/Python/LFM-helio/lib/'))
import lfmhlib

#----------- PARSE ARGUMENTS ---------#
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('hdfFileName',help='LFM-helio file to use')
parser.add_argument('shell',help='Radial shell in code units of distance.',type=float)
 #parser.add_argument('limits',help='',type=float)
parser.add_argument('Tsolar',help='T solar rotation in code units.',type=float,default=25.38*24*3600,nargs='?')
args = parser.parse_args()
#----------- PARSE ARGUMENTS ---------#


# get LFM stuff out
simtime = lfmhlib.get_time(args.hdfFileName)
R,theta,phi = lfmhlib.r_theta_phi_uniform(args.hdfFileName)
R/=6.96e10   # FIX ME HARD CODED UNITS
islice=argmin(abs(R-args.shell))   # index of the desired radial distance


br = lfmhlib.read_var(args.hdfFileName,'br')
vr = lfmhlib.read_var(args.hdfFileName,'vr')
rho = lfmhlib.read_var(args.hdfFileName,'rho')
c = lfmhlib.read_var(args.hdfFileName,'c')

# shift to start of CR
br = lfmhlib.time_shift(simtime,br[:,:,islice],args.Tsolar)
vr = lfmhlib.time_shift(simtime,vr[:,:,islice],args.Tsolar)
rho = lfmhlib.time_shift(simtime,rho[:,:,islice],args.Tsolar)
c = lfmhlib.time_shift(simtime,c[:,:,islice],args.Tsolar)


# Plot
matplotlib.rc('mathtext',fontset='stixsans',default='regular')
fig=plt.figure(figsize=(12,12))
fig.suptitle('Radial distance= %.2f' % args.shell)

p1c = 0.5*(phi[1:]+phi[:-1])
t1c = 0.5*(theta[1:]+theta[:-1])

# Br
ax1=plt.subplot(4,1,1,aspect='equal')
p1=plt.pcolormesh(phi*180/pi,(pi/2.-theta)*180./pi,br.T*1.e5,vmin=-5.,vmax=5.,cmap='RdBu_r',rasterized=True)
plt.colorbar(p1,ax=ax1).set_label(r'B$_r$, nT')
ax1.contour(p1c*180./pi,(pi/2-t1c)*180/pi,br.T,[0.],colors='black')

# Vr
ax2=plt.subplot(4,1,2,aspect='equal')
p2=plt.pcolormesh(phi*180/pi,(pi/2.-theta)*180./pi,vr.T*1.e-5,vmin=200.,vmax=800.,rasterized=True)
plt.colorbar(p2,ax=ax2).set_label(r'V$_r$, km/s')
plt.axhline(0)
ax2.contour(p1c*180./pi,(pi/2-t1c)*180/pi,br.T,[0.],colors='white')

# rho
ax3=plt.subplot(4,1,3,aspect='equal')
p3=plt.pcolormesh(phi*180/pi,(pi/2-theta)*180./pi,log10(rho.T/1.67e-24),vmin=0.4,vmax=1.5,rasterized=True)#vmin=log10(rho.min()/1.67e-24),vmax=log10(rho.max()/1.67e-24),rasterized=True)
plt.colorbar(p3,ax=ax3).set_label(r'log$_{10}(\rho\,\mathrm{cm}^{-3}$)')
ax3.contour(p1c*180./pi,(pi/2-t1c)*180/pi,br.T,[0.],colors='white')

# temperature
T=3./5.*((c+1.)*1.e-2)**2*1.67e-27/1.38e-23 # T in K  # FIX ME HARD CODE UNITS AND GAMMA
ax4=plt.subplot(4,1,4,aspect='equal')
p4=plt.pcolormesh(phi*180/pi,(pi/2-theta)*180./pi,log10(T.T),vmin=4.5,vmax=5.4,rasterized=True)#log10(T.min()),vmax=log10(T.max()),rasterized=True)
plt.colorbar(p4,ax=ax4).set_label(r'log$_{10}$(T, K)')
ax4.contour(p1c*180./pi,(pi/2-t1c)*180/pi,br.T,[0.],colors='white')


for ax in [ax1,ax2,ax3,ax4]:
    if ax!=ax4:
        ax.set_xticklabels([])
    else:
        ax.set_xlabel(r'$\phi$')
        ax.set_ylabel(r'$\theta$')
        ax.set_xlim(0,360)
        ax.set_ylim(-70,70)

#plt.show()
plot.savefig(args.hdfFileName[:-4]+'.png')
 
 
