#!/usr/bin/python
import sys,os
from numpy import argmin,abs,pi
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
parser.add_argument('variable',help='What to plot')
parser.add_argument('Tsolar',help='T solar rotation in code units.',type=float,default=25.38*24*3600,nargs='?')
args = parser.parse_args()
#----------- PARSE ARGUMENTS ---------#

# get LFM stuff out
simtime = lfmhlib.get_time(args.hdfFileName)
R,theta,phi = lfmhlib.r_theta_phi_uniform(args.hdfFileName)
R/=6.96e10
islice=argmin(abs(R-args.shell))   # index of the desired radial distance
var = lfmhlib.read_var(args.hdfFileName,args.variable)
var = lfmhlib.time_shift(simtime,var[:,:,islice],args.Tsolar)


# Plot
matplotlib.rc('mathtext',fontset='stixsans',default='regular')
fig=plt.figure(figsize=(19.5,5.5))
fig.suptitle('Radial distance= %.2f' % args.shell)

ax1=plt.subplot(1,1,1,aspect='equal')
#plot(UT,nl1)
#plot(time,rho[:,nj/2,1]/1.67e-24/1.16/100.)
p=plt.pcolormesh(phi*180/pi,90.-theta*180./pi,var.T*1.e5,vmin=-1.,vmax=1.,cmap='RdBu_r',rasterized=True)
plt.colorbar(p,ax=ax1).set_label(args.variable)
#setp( ax1.get_xticklabels(), visible=False)
ax1.set_ylabel(r'$\theta$')
ax1.set_xlabel(r'$\phi$')
ax1.set_xlim(0,360)
ax1.set_ylim(-70,70)

plt.show()
