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
parser.add_argument('crStartTime',help='CR start time (YYYY MM DD HH MM)')
parser.add_argument('shell',help='Radial shell in code units of distance.',type=float)
parser.add_argument('variable',help='What to plot')
args = parser.parse_args()
#----------- PARSE ARGUMENTS ---------#


# CR start time from umtof.umd.edu/pm/crn
import datetime
start_time=datetime.datetime.strptime(args.crStartTime,'%Y %m %d %H %M')


# get LFM stuff out
simtime = lfmhlib.get_time(args.hdfFileName)
R,theta,phi = lfmhlib.r_theta_phi_uniform(args.hdfFileName)
islice=argmin(abs(R-args.shell))   # index of the desired radial distance
var = lfmhlib.read_var(args.hdfFileName,args.variable)
var = 0.5*(var[:,(theta.shape[0]-1)/2-1,islice] +
           var[:,(theta.shape[0]-1)/2,islice])   # extract equatorial value  # extract equatorial value
lfmhlib.time_shift(simtime,var,25.38*24*3600)


# Plot
# ========== time dependence ==========
ut = lfmhlib.get_UT(simtime,phi)
ut=[start_time+datetime.timedelta(seconds=float(ttt)) for ttt in ut]

plt.plot(ut,var[::-1])
"""
matplotlib.rc('mathtext',fontset='stixsans',default='regular')
fig=plt.figure(figsize=(19.5,5.5))
ax1=plt.subplot(1,1,1,aspect='equal')
#plot(UT,nl1)
#plot(time,rho[:,nj/2,1]/1.67e-24/1.16/100.)
p=plt.pcolormesh(phi*180/pi,90.-theta*180./pi,var[:,:,islice].T)
plt.colorbar(p,ax=ax1).set_label(args.variable)
#setp( ax1.get_xticklabels(), visible=False)
ax1.set_ylabel(r'$\theta$')
ax1.set_xlabel(r'$\phi$')
ax1.set_xlim(0,360)
ax1.set_ylim(-70,70)

ax2=subplot(3,1,2)
plot(UT,vl1)
setp( ax2.get_xticklabels(), visible=False)
ylabel('V, km/s')
#xlim(0,28)

ax2=subplot(3,1,3)
plot(UT,bl1)
setp( ax1.get_xticklabels(), visible=False)
ylabel('|B|, nT')
xlabel('Time, days')
#xlim(0,28)


show()

"""
