from pylab import *
import glob,os,sys,datetime
import summary
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.dates import DateFormatter

home = os.path.expanduser("~");
if os.path.join(home,'GitHub/Python/LFM-helio/lib/') not in sys.path: sys.path.append(os.path.join(home,'GitHub/Python/LFM-helio/lib/'))
import lfmhlib,ace,cspice

####################################################################################################
########### CR 2094 ################################################################################
####################################################################################################
"""
lfmFileDir = '/glade/scratch/vgm/LFM-helio/ADAPT/2094rf'
lfmFiles = sorted(glob.glob(os.path.join(lfmFileDir,'*v1*.hdf')))
start_step = 25000
#start_file = 'cr2094_192x96x192_v1_mhd_0025000.hdf'
start_dati = datetime.datetime(2010,1,1,23,0,0)
nk=192
Omega = 2*pi/27.27/24./3600.

# get ACE position
ACE_dati,ACE_x,ACE_y,ACE_z,ACE_r,ACE_theta,ACE_phi=cspice.cspice('/glade/u/home/vgm/LFM-helio_2.0/SCDATA/Earth_cspice_HCI_position_2010_001_060')

# get ACE data
ACE_mf_time,ACE_bx,ACE_by,ACE_bz=ace.cdaweb_mfi('/glade/u/home/vgm/LFM-helio_2.0/SCDATA/ACE/2094/AC_H2_MFI_14914.txt')
ACE_plasma_time,ACE_v,ACE_n=ace.cdaweb_swe('/glade/u/home/vgm/LFM-helio_2.0/SCDATA/ACE/2094/AC_H2_SWE_14914.txt')
"""
####################################################################################################

####################################################################################################
########### CR 2066  ###############################################################################
####################################################################################################
lfmFileDir = '/glade/scratch/vgm/LFM-helio/ADAPT/2066rf'
lfmFiles = sorted(glob.glob(os.path.join(lfmFileDir,'*.hdf')))
start_step = 25000
#start_file = 'cr2094_192x96x192_v1_mhd_0025000.hdf'
start_dati = datetime.datetime(2008,1,1,23,0,0)
nk=192
Omega = 2*pi/27.27/24./3600.


# get ACE position
ACE_dati,ACE_x,ACE_y,ACE_z,ACE_r,ACE_theta,ACE_phi=cspice.cspice('/glade/u/home/vgm/LFM-helio_2.0/SCDATA/ACE/Earth_cspice_HCI_position_2008_001_060')

# get MSGR position
MSGR_dati,MSGR_x,MSGR_y,MSGR_z,MSGR_r,MSGR_theta,MSGR_phi=cspice.cspice('/glade/u/home/vgm/LFM-helio_2.0/SCDATA/MSGR/MSGR_cspice_HCI_position_2008_001_060')

# get ACE data
ACE_mf_time,ACE_bx,ACE_by,ACE_bz=ace.cdaweb_mfi('/glade/u/home/vgm/LFM-helio_2.0/SCDATA/ACE/2066/AC_H2_MFI_20155.txt')
ACE_plasma_time,ACE_v,ACE_n=ace.cdaweb_swe('/glade/u/home/vgm/LFM-helio_2.0/SCDATA/ACE/2066/AC_H2_SWE_20155.txt')
####################################################################################################

####################################################################################################
########### YEAR 2012###############################################################################
####################################################################################################
# lfmFileDir = '/glade/scratch/vgm/LFM-helio/ADAPT/2012rf'
# lfmFiles = sorted(glob.glob(os.path.join(lfmFileDir,'*.hdf')))
# start_step = 25000
# #start_file = 'cr2094_192x96x192_v1_mhd_0025000.hdf'
# start_dati = datetime.datetime(2012,2,15,0,0,0)
# nk=192
# Omega = 2*pi/27.27/24./3600.


# # get ACE position
# ACE_dati,ACE_x,ACE_y,ACE_z,ACE_r,ACE_theta,ACE_phi=cspice.cspice('/glade/u/home/vgm/LFM-helio_2.0/SCDATA/ACE/Earth_cspice_HCI_position_2012_001_120')

# # get ACE data
# ACE_mf_time,ACE_bx,ACE_by,ACE_bz=ace.cdaweb_mfi('/glade/u/home/vgm/LFM-helio_2.0/SCDATA/ACE/AC_H2_MFI_7837.txt')
# ACE_plasma_time,ACE_v,ACE_n=ace.cdaweb_swe('/glade/u/home/vgm/LFM-helio_2.0/SCDATA/ACE/AC_H2_SWE_7837.txt')
####################################################################################################



#start_ind = lfmFiles.index(os.path.join(lfmFileDir,start_file))
start_ind = [start_step==int(f[-11:-4]) for f in lfmFiles].index(True)
dates = []
vr = []
br = []
vr0= []
br0= []

matplotlib.rc('mathtext',fontset='stixsans',default='regular')
close('all')
fig=plt.figure(figsize=(12,15))
gs = gridspec.GridSpec(2,1,height_ratios=[6,4],hspace=0.2)
gs1 = gridspec.GridSpecFromSubplotSpec(2,1,subplot_spec=gs[1],hspace=0.1)
for i,f in enumerate(lfmFiles[start_ind:]):
    if i==0: 
        simtime0 = lfmhlib.get_time(f)
        R,theta,phi = lfmhlib.r_theta_phi_uniform(f)
        R/=6.96e10   # FIX ME HARD CODED UNITS

    simtime = lfmhlib.get_time(f)

    d = start_dati+datetime.timedelta(seconds=simtime-simtime0)
    dates.append(d)
    print(d)

    # look for the corresponding ACE record
    aceind = argmin([abs((ddd-d).total_seconds()) for ddd in ACE_dati])

    pace = (pi-Omega*(simtime-simtime0))%(2*pi)
    tace = 90.-ACE_theta[aceind]

    # look for the corresponding MSGR record
    msgrind = argmin([abs((ddd-d).total_seconds()) for ddd in MSGR_dati])
    pace1 = (pi-Omega*(simtime-simtime0))
    pmsgr = (pace1 - (ACE_phi[aceind]-MSGR_phi[msgrind])*pi/180.)%(2*pi)
    tmsgr = 90.-MSGR_theta[msgrind]


    iace=argmin(abs(R-ACE_r[aceind]*215.))   # index of the desired radial distance
    jace=argmin(abs(pi/2.-theta-tace*pi/180.))   # index of the desired theta-angle
    kace = floor(pace/(2*pi)*nk)

    v = lfmhlib.read_var(f,'vr')
    b = lfmhlib.read_var(f,'br')
    vr.append(v[kace,jace,iace])
    br.append(b[kace,jace,iace])

    v= lfmhlib.read_var(lfmFiles[start_ind],'vr')
    b= lfmhlib.read_var(lfmFiles[start_ind],'br')
    vr0.append(v[kace,jace,iace])
    br0.append(b[kace,jace,iace])

    fig.suptitle(d.strftime('%Y:%m:%d %H:%M'),y=0.95)
    
#    (ax1,ax2,ax3,ax4) = summary.plots(f,iace,1.e25)
    (ax1,ax2) = summary.plotbrvr(f,iace,1.e25,outer_grid=gs[0])
    ax1.scatter(pace*180/pi,tace)
    ax2.scatter(pace*180/pi,tace)
    ax1.scatter(pmsgr*180/pi,tmsgr,color='purple')
    ax2.scatter(pmsgr*180/pi,tmsgr,color='purple')


    # ACE data
    ax3 = plt.subplot(gs1[0])
    ax4 = plt.subplot(gs1[1],sharex=ax3)

    ax3.plot(ACE_plasma_time,ACE_v,'r',label='ACE')
    ax3.legend()
    ax3.plot(dates,array(vr0)*1.e-5,color='gray')
    ax3.plot(dates,array(vr)*1.e-5,'k-')
    ax3.axvline(d,color='black')
    ax3.set_ylim(200,800)
    ax3.set_ylabel(r'V$_r$, km/s')

    ax4.plot(ACE_mf_time,-ACE_bx,'r')
    ax4.plot(dates,array(br0)*1.e5,'gray')
    ax4.plot(dates,array(br)*1.e5,'k-')
    ax4.axvline(d,color='black')
    ax4.xaxis.set_major_formatter(DateFormatter('%m/%d'))
    ax4.set_ylim(-10,10)
    ax4.set_ylabel(r'B$_r$, nT')
    ax4.axhline(0.,color='black',linestyle='--')
    plt.setp(ax3.get_xticklabels(),visible=False)

    
    savefig(f[:-4]+'.png',bbox_inches='tight')
    clf()

    
