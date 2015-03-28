#!/usr/bin/python
import sys,os
from numpy import argmin,abs,pi,log10
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
home = os.path.expanduser("~");
if os.path.join(home,'GitHub/Python/LFM-helio/lib/') not in sys.path: sys.path.append(os.path.join(home,'GitHub/Python/LFM-helio/lib/'))
import lfmhlib

def plots(fName,shell,tsolar,gamma):
    # get LFM stuff out
    simtime = lfmhlib.get_time(fName)
    R,theta,phi = lfmhlib.r_theta_phi_uniform(fName)
    R/=6.96e10   # FIX ME HARD CODED UNITS
    islice=argmin(abs(R-shell))   # index of the desired radial distance


    br = lfmhlib.read_var(fName,'br')
    vr = lfmhlib.read_var(fName,'vr')
    rho = lfmhlib.read_var(fName,'rho')
    c = lfmhlib.read_var(fName,'c')

    # shift to start of CR
    br = lfmhlib.time_shift(simtime,br[:,:,islice],tsolar)
    vr = lfmhlib.time_shift(simtime,vr[:,:,islice],tsolar)
    rho = lfmhlib.time_shift(simtime,rho[:,:,islice],tsolar)
    c = lfmhlib.time_shift(simtime,c[:,:,islice],tsolar)


    # Plot
    gs = gridspec.GridSpec(4,2,width_ratios=[14,1])

    p1c = 0.5*(phi[1:]+phi[:-1])
    t1c = 0.5*(theta[1:]+theta[:-1])

    # Br
    ax1=plt.subplot(gs[0,0],aspect='equal')
    ax1_cb=plt.subplot(gs[0,1])
    p1=ax1.pcolormesh(phi*180/pi,(pi/2.-theta)*180./pi,br.T*1.e5,vmin=-120.,vmax=120.,cmap='RdBu_r',rasterized=True)
    plt.colorbar(p1,cax=ax1_cb).set_label(r'B$_r$, nT')
    ax1.contour(p1c*180./pi,(pi/2-t1c)*180/pi,br.T,[0.],colors='black')

    # Vr
    ax2=plt.subplot(gs[1,0],aspect='equal',sharex=ax1)
    ax2_cb=plt.subplot(gs[1,1])
    p2=ax2.pcolormesh(phi*180/pi,(pi/2.-theta)*180./pi,vr.T*1.e-5,vmin=200.,vmax=800.,rasterized=True)
    plt.colorbar(p2,cax=ax2_cb).set_label(r'V$_r$, km/s')
    plt.axhline(0)
    ax2.contour(p1c*180./pi,(pi/2-t1c)*180/pi,br.T,[0.],colors='white')

    # rho
    ax3=plt.subplot(gs[2,0],aspect='equal',sharex=ax1)
    ax3_cb=plt.subplot(gs[2,1])
#    p3=ax3.pcolormesh(phi*180/pi,(pi/2-theta)*180./pi,log10(rho.T/1.67e-24),vmin=0.4,vmax=1.5,rasterized=True)
    p3=ax3.pcolormesh(phi*180/pi,(pi/2-theta)*180./pi,rho.T/1.67e-24,vmin=200,vmax=1600,rasterized=True)
#    plt.colorbar(p3,cax=ax3_cb).set_label(r'log$_{10}(n\,\mathrm{cm}^{-3}$)')
    plt.colorbar(p3,cax=ax3_cb).set_label(r'n, cm$^{-3}$')
    ax3.contour(p1c*180./pi,(pi/2-t1c)*180/pi,br.T,[0.],colors='white')

    # temperature
    T=1./gamma*c**2*1.67e-24/1.38e-16 # T in K  # FIX ME HARD CODE UNITS AND GAMMA
    ax4=plt.subplot(gs[3,0],aspect='equal',sharex=ax1)
    ax4_cb=plt.subplot(gs[3,1])
#    p4=ax4.pcolormesh(phi*180/pi,(pi/2-theta)*180./pi,log10(T.T),vmin=4.5,vmax=5.4,rasterized=True)
    p4=ax4.pcolormesh(phi*180/pi,(pi/2-theta)*180./pi,T.T*1.e-6,vmin=0.1,vmax=2,rasterized=True)
#    plt.colorbar(p4,cax=ax4_cb).set_label(r'log$_{10}$(T, K)')
    plt.colorbar(p4,cax=ax4_cb).set_label(r'T, MK')
    ax4.contour(p1c*180./pi,(pi/2-t1c)*180/pi,br.T,[0.],colors='white')


    for ax in [ax1,ax2,ax3,ax4]:
        ax.set_ylabel(r'$\theta$')
        ax.set_xlim(0,360)
        ax.set_ylim(-90,90)

        if ax!=ax4:
            plt.setp(ax.get_xticklabels(), visible=False) 
        else:
            ax.set_xlabel(r'$\phi$')

    return(ax1,ax2,ax3,ax4)

#plt.show()
#    plt.savefig(fName[:-4]+'.png',bbox_inches='tight')
#    plt.clf()

def plotbrvr(fName,shell,tsolar,outer_grid=None,lims={'br':[-5.,5.],
                                      'vr':[200.,800]}):
    # get LFM stuff out
    simtime = lfmhlib.get_time(fName)
    R,theta,phi = lfmhlib.r_theta_phi_uniform(fName)
    R/=6.96e10   # FIX ME HARD CODED UNITS
#    islice=argmin(abs(R-shell))   # index of the desired radial distance
    islice = shell


    br = lfmhlib.read_var(fName,'br')
    vr = lfmhlib.read_var(fName,'vr')

    # shift to start of CR
    br = lfmhlib.time_shift(simtime,br[:,:,islice],tsolar)
    vr = lfmhlib.time_shift(simtime,vr[:,:,islice],tsolar)

    # Plot
    if outer_grid==None:
        gs = gridspec.GridSpec(2,2,width_ratios=[14,1],hspace=0.1)
    else:
        gs = gridspec.GridSpecFromSubplotSpec(2,2,subplot_spec=outer_grid,width_ratios=[14,1],hspace=0.1)

    p1c = 0.5*(phi[1:]+phi[:-1])
    t1c = 0.5*(theta[1:]+theta[:-1])

    # Br
    ax1=plt.subplot(gs[0,0],aspect='equal')
    ax1_cb=plt.subplot(gs[0,1])
    p1=ax1.pcolormesh(phi*180/pi,(pi/2.-theta)*180./pi,br.T*1.e5,vmin=lims['br'][0],vmax=lims['br'][1],cmap='RdBu_r',rasterized=True)
    plt.colorbar(p1,cax=ax1_cb).set_label(r'B$_r$, nT')
    ax1.contour(p1c*180./pi,(pi/2-t1c)*180/pi,br.T,[0.],colors='black')
    ax1.axhline(0.,linestyle='--',color='black')

    # Vr
    ax2=plt.subplot(gs[1,0],aspect='equal',sharex=ax1)
    ax2_cb=plt.subplot(gs[1,1])
    p2=ax2.pcolormesh(phi*180/pi,(pi/2.-theta)*180./pi,vr.T*1.e-5,vmin=lims['vr'][0],vmax=lims['vr'][1],rasterized=True)
    plt.colorbar(p2,cax=ax2_cb).set_label(r'V$_r$, km/s')
    plt.axhline(0)
    ax2.contour(p1c*180./pi,(pi/2-t1c)*180/pi,br.T,[0.],colors='white')
    ax2.axhline(0.,linestyle='--',color='white')

    for ax in [ax1,ax2]:
        ax.set_ylabel(r'$\theta$')
        ax.set_xlim(0,360)
        ax.set_ylim(-70,70)

        if ax!=ax2:
            plt.setp(ax.get_xticklabels(), visible=False) 
        else:
            ax.set_xlabel(r'$\phi$')

    return(ax1,ax2)

 
 
if __name__ == "__main__":
    #----------- PARSE ARGUMENTS ---------#
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('hdfFileName',help='LFM-helio file to use')
    parser.add_argument('shell',help='Radial shell in code units of distance.',type=float)
    #parser.add_argument('limits',help='',type=float)
    parser.add_argument('Tsolar',help='T solar rotation in code units.',type=float,default=25.38*24*3600,nargs='?')
    parser.add_argument('gamma',help='Gamma.',type=float,default=1.5,nargs='?')
    args = parser.parse_args()
    #----------- PARSE ARGUMENTS ---------#

    matplotlib.rc('mathtext',fontset='stixsans',default='regular')
    fig=plt.figure(figsize=(8,8))
    fig.suptitle('Radial distance= %.2f' % args.shell,y=0.95)
    plots(args.hdfFileName,args.shell,args.Tsolar,args.gamma)
    plt.show()
