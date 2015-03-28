from pylab import *
import glob,os,sys,datetime
import summary
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
from matplotlib.dates import DateFormatter

home = os.path.expanduser("~");
if os.path.join(home,'GitHub/Python/LFM-helio/lib/') not in sys.path: sys.path.append(os.path.join(home,'GitHub/Python/LFM-helio/lib/'))
import lfmhlib,ace,cspice,msgr



####################################################################################################
########### CR 2094 ################################################################################
"""
lfmFileDir = '/glade/scratch/vgm/LFM-helio/ADAPT/2094rf'
lfmFiles = sorted(glob.glob(os.path.join(lfmFileDir,'cr2094_192x96x192_v1_mhd*.hdf')))
start_step = 25000
start_dati = datetime.datetime(2010,1,1,23,0,0)
nk=192
Omega = 2*pi/27.27/24./3600.

# get MSGR position
MSGR_dati,MSGR_x,MSGR_y,MSGR_z,MSGR_r,MSGR_theta,MSGR_phi=cspice.cspice('/glade/u/home/vgm/LFM-helio_2.0/SCDATA/MSGR/MSGR_cspice_HCI_position_2010_001_060')

# get ACE position
ACE_dati,ACE_x,ACE_y,ACE_z,ACE_r,ACE_theta,ACE_phi=cspice.cspice('/glade/u/home/vgm/LFM-helio_2.0/SCDATA/ACE/Earth_cspice_HCI_position_2010_001_060')

# get MSGR data
#MSGR_dati,MSGR_br,MSGR_bt,MSGR_bn=msgr.msgr('/glade/u/home/vgm/LFM-helio_2.0/SCDATA/MSGR/MAGRTNSCIAVG10_001_035_60_V06.TAB')
MSGR_dati,MSGR_br,MSGR_bt,MSGR_bn= msgr.msgr_v1('/glade/u/home/vgm/LFM-helio_2.0/SCDATA/MSGR/mag_10_001-059_1hour_promitjos.dat',start_dati.year)
"""
####################################################################################################

####################################################################################################
########### CR 2066 ################################################################################
lfmFileDir = '/glade/scratch/vgm/LFM-helio/ADAPT/2066rf'
lfmFiles = sorted(glob.glob(os.path.join(lfmFileDir,'cr2066_192x96x192_mhd*.hdf')))
start_step = 25000
start_dati = datetime.datetime(2008,1,1,23,0,0)
nk=192
Omega = 2*pi/27.27/24./3600.

# get MSGR position
MSGR_dati,MSGR_x,MSGR_y,MSGR_z,MSGR_r,MSGR_theta,MSGR_phi=cspice.cspice('/glade/u/home/vgm/LFM-helio_2.0/SCDATA/MSGR/MSGR_cspice_HCI_position_2008_001_060')

# get ACE position
ACE_dati,ACE_x,ACE_y,ACE_z,ACE_r,ACE_theta,ACE_phi=cspice.cspice('/glade/u/home/vgm/LFM-helio_2.0/SCDATA/ACE/Earth_cspice_HCI_position_2008_001_060')


# get MSGR data
#MSGR_dati,MSGR_br,MSGR_bt,MSGR_bn=msgr.msgr('/glade/u/home/vgm/LFM-helio_2.0/SCDATA/MSGR/MAGRTNSCIAVG10_001_035_60_V06.TAB')
MSGR_dati,MSGR_br,MSGR_bt,MSGR_bn= msgr.msgr_v1('/glade/u/home/vgm/LFM-helio_2.0/SCDATA/MSGR/mag_08_001-055_1hour_promitjos.dat',start_dati.year)
####################################################################################################


start_ind = [start_step==int(f[-11:-4]) for f in lfmFiles].index(True)
dates = []
vr = []
br = []


matplotlib.rc('mathtext',fontset='stixsans',default='regular')
fig=plt.figure(figsize=(10,10))
gs = gridspec.GridSpec(2,3,width_ratios=[20,1,6],height_ratios=[2,1])
for i,f in enumerate(lfmFiles[start_ind:]):
    simtime = lfmhlib.get_time(f)
    if i==0: 
        simtime0 = simtime
        R,theta,phi = lfmhlib.r_theta_phi_uniform(f)
        R/=6.96e10   # FIX ME HARD CODED UNITS
    
    d = start_dati+datetime.timedelta(seconds=simtime-simtime0)
    dates.append(d)
    print(d)

    # look for the corresponding ACE record
    aceind = argmin([abs((ddd-d).total_seconds()) for ddd in ACE_dati])

    # look for the corresponding MSGR record
    msgrind = argmin([abs((ddd-d).total_seconds()) for ddd in MSGR_dati])

    # Need ACE longitude to figure out that of MSGR
    pace = (pi-Omega*(simtime-simtime0)) #%(2*pi)

    pmsgr = (pace - (ACE_phi[aceind]-MSGR_phi[msgrind])*pi/180.)%(2*pi)
    tmsgr = 90.-MSGR_theta[msgrind]


    imsgr=argmin(abs(R-MSGR_r[msgrind]*215.))   # index of the desired radial distance
    jmsgr=argmin(abs(pi/2.-theta-tmsgr*pi/180.))   # index of the desired theta-angle
    kmsgr = floor(pmsgr/(2*pi)*nk)

    b = lfmhlib.read_var(f,'br')
    br.append(b[kmsgr,jmsgr,imsgr])

    fig.suptitle(dates[i].strftime('%Y:%m:%d %H:%M'),y=0.95)
    
    rc=0.5*(R[:-1]+R[1:])
    phic=0.5*(phi[:-1]+phi[1:])
    phic=append(phic,phic[0])
    rr,pp = meshgrid(rc*sin((90.-tmsgr)*pi/180.),phic)
    ax1=subplot(gs[0,0],aspect='equal')
    bplot = vstack([b[:,jmsgr,:],b[0,jmsgr,:]])*(rc/rc[0])**2*1.e5
    p=ax1.pcolormesh(rr*cos(pp),rr*sin(pp),bplot,shading='gouraud',vmin=-160,vmax=160,cmap=cm.RdBu_r)
    ax1.set_xlim(-230,230)
    ax1.set_ylim(-230,230)
    ax1.set_xlabel('X, Rs')
    ax1.set_ylabel('Y, Rs')
    cb_ax = subplot(gs[0,1])
    colorbar(p,cax=cb_ax).set_label(r'B$_r\times (r/r_0)^2$, nT')

    ax1.scatter(MSGR_r[msgrind]*215.*sin((90.-tmsgr)*pi/180.)*cos(pmsgr),
                MSGR_r[msgrind]*215.*sin((90.-tmsgr)*pi/180.)*sin(pmsgr))
    ax1.scatter(ACE_r[aceind]*215.*sin((90.-tmsgr)*pi/180.)*cos(pace),
                ACE_r[aceind]*215.*sin((90.-tmsgr)*pi/180.)*sin(pace))

    # MSGR data
    ax2 = plt.subplot(gs[1,:])
    ax2.plot(MSGR_dati,MSGR_br,'r',label='MSGR')
    ax2.set_xlim((MSGR_dati[0],MSGR_dati[-1]))  # this fixes the size of axis so they don't change in the movie
    ax2.legend()
    ax2.plot(dates,array(br)*1.e5,'k-')
    ax2.axvline(d,color='black')
    ax2.xaxis.set_major_formatter(DateFormatter('%m/%d'))
    ax2.set_ylim(-50,50)
    ax2.set_ylabel(r'B$_r$, nT')
    ax2.axhline(0.,color='black',linestyle='--')
    
    savefig(f[:-4]+'.png',bbox_inches='tight')
    clf()

    
