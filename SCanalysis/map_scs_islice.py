import matplotlib
from matplotlib import pyplot as plt
from matplotlib.dates import DateFormatter
import glob,os,sys,datetime
from numpy import pi,argmin
import time

home = os.path.expanduser("~");
if os.path.join(home,'GitHub/Python/LFM-helio/lib/') not in sys.path: sys.path.append(os.path.join(home,'GitHub/Python/LFM-helio/lib/'))
if os.path.join(home,'GitHub/Python/LFM-helio/analysis/') not in sys.path: sys.path.append(os.path.join(home,'GitHub/Python/LFM-helio/analysis/'))
import lfmhlib,ace,cspice,stereo
import summary

####################################################################################################
############## SET UP PARAMETERS      ##############################################################
####################################################################################################
save_pdf = True




####################################################################################################
##############  THIS IS FOR 2008      ##############################################################
####################################################################################################
#lfmFileDir = '/glade/p/ujhb0003/LFM-helio/ADAPT/2066rf_fullpb_192x192x768_converged'
#lfmFileDir = '/glade/scratch/vgm/LFM-helio/ADAPT/2066rf_fullpb_192x192x768_converged-smoothed'
#lfmFileDir = '/glade/scratch/vgm/LFM-helio/ADAPT/2066rf_fullpb_192x192x768_converged-vrsmoothed'
#lfmFileDir = '/glade/scratch/vgm/LFM-helio/ADAPT/2066rf_fullpb_192x192x768_converged-vrsmoothed16'
#lfmFileDir = '/glade/scratch/vgm/LFM-helio/ADAPT/2066rf_fullpb_192x192x768_converged-rhovrsmoothed'
#lfmFileDir = '/glade/scratch/vgm/LFM-helio/ADAPT/2066rf_fullpb_192x192x768_converged_fixedboundaries'
#lfmFileDir = '/glade/scratch/vgm/LFM-helio/ADAPT/2066rf_fullpb_192x192x768_converged-vrsmoothed_fixedboundaries'
lfmFileDir = '/glade/scratch/vgm/LFM-helio/ADAPT/2066rf_fullpb_192x192x768_converged_fixedboundaries_stairfix'
#lfmFileDir = '/glade/scratch/vgm/LFM-helio/ADAPT/2066rf_fullpb_192x192x768_converged-vrsmoothed_fixedboundaries_stairfix'
#lfmFileDir = '/glade/scratch/vgm/LFM-helio/ADAPT/2066rf_fullpb_192x192x768_converged-rhovrsmoothed_fixedboundaries_stairfix'
start_step = 64000
start_dati = datetime.datetime(2008,1,1,23,0,0)

ace_loc_file = '/glade/u/home/vgm/LFM-helio_2.0/SCDATA/ACE/Earth_cspice_HCI_position_2008_001_060'
sta_loc_file = '/glade/u/home/vgm/LFM-helio_2.0/SCDATA/STEREO/STA_scpice_HCI_position_2008_001_060'
stb_loc_file = '/glade/u/home/vgm/LFM-helio_2.0/SCDATA/STEREO/STB_scpice_HCI_position_2008_001_060'
msgr_loc_file = '/glade/u/home/vgm/LFM-helio_2.0/SCDATA/MSGR/MSGR_cspice_HCI_position_2008_001_060'
####################################################################################################

####################################################################################################
##############  THIS IS FOR 2010      ##############################################################
####################################################################################################
# lfmFileDir = '/glade/scratch/vgm/LFM-helio/ADAPT/2094rf_fullpb_v1_192x192x768_converged'
# start_step = 61600
# start_dati = datetime.datetime(2010,1,1,23,0,0)

# ace_loc_file = '/glade/u/home/vgm/LFM-helio_2.0/SCDATA/ACE/Earth_cspice_HCI_position_2010_001_060'
# sta_loc_file = '/glade/u/home/vgm/LFM-helio_2.0/SCDATA/STEREO/STA_scpice_HCI_position_2010_001_060'
# stb_loc_file = '/glade/u/home/vgm/LFM-helio_2.0/SCDATA/STEREO/STB_scpice_HCI_position_2010_001_060'
# msgr_loc_file = '/glade/u/home/vgm/LFM-helio_2.0/SCDATA/MSGR/MSGR_cspice_HCI_position_2010_001_060'
####################################################################################################


Omega = 2*pi/27.27/24./3600.   #  note, the sun rotation rate used in
                               #  the runs was 25.38 days. But here we
                               #  have to use 27.27 because this is
                               #  the rate at which Earth moves
                               #  relative to the frame rotating with
                               #  the sun
####################################################################################################
home = os.path.expanduser("~");
if os.path.join(home,'GitHub/Python/LFM-helio/lib/') not in sys.path: sys.path.append(os.path.join(home,'GitHub/Python/LFM-helio/lib/'))
if os.path.join(home,'GitHub/Python/LFM-helio/analysis/') not in sys.path: sys.path.append(os.path.join(home,'GitHub/Python/LFM-helio/analysis/'))
import lfmhlib,ace,cspice

# get ACE position
ACE_dati,ACE_x,ACE_y,ACE_z,ACE_r,ACE_theta,ACE_phi=cspice.cspice(ace_loc_file)

# GET THIS FOR OTHER SCS FROM MAP_SCS.PY
# get STEREO A position
STA_dati,STA_x,STA_y,STA_z,STA_r,STA_theta,STA_phi=cspice.cspice(sta_loc_file)

STB_dati,STB_x,STB_y,STB_z,STB_r,STB_theta,STB_phi=cspice.cspice(stb_loc_file)

# get MSGR position
MSGR_dati,MSGR_x,MSGR_y,MSGR_z,MSGR_r,MSGR_theta,MSGR_phi=cspice.cspice(msgr_loc_file)
####################################################################################################
lfmFiles = sorted(glob.glob(os.path.join(lfmFileDir,'*.hdf')))
start_ind = [start_step==int(f[-11:-4]) for f in lfmFiles].index(True)

plt.style.use('publication')
plt.close('all')

## Note, we use analysis/summary below to make the plots. That
## function defines y-limits, but they may get overriden by the figure
## size below as mpl is trying to plot them with equal aspect ratio
## within the figure size constraints that we've given
fig = plt.figure(figsize=(30/6.,45./6.)) 

for i,f in enumerate(lfmFiles[start_ind:]):
    if i==0: 
        simtime0 = lfmhlib.get_time(lfmFiles[0])
        R,theta,phi,Rc,thetac,phic = lfmhlib.r_theta_phi_uniform(f)
        nk = phic.size

    simtime = lfmhlib.get_time(f)

    d = start_dati+datetime.timedelta(seconds=simtime-simtime0)
    print(d)
#    if ( abs((d-datetime.datetime(2008,1,12,3,51)).total_seconds())>300.): continue    # for staircase fix
    if ( abs((d-datetime.datetime(2008,1,25,5,7)).total_seconds())>300.): continue    # for staircase fix
#    if ( abs((d-datetime.datetime(2008,1,12,3,21)).total_seconds())>300.): continue
#    if ( abs((d-datetime.datetime(2008,1,25,5,16)).total_seconds())>300.): continue

    # look for the corresponding ACE record
    aceind = argmin([abs((ddd-d).total_seconds()) for ddd in ACE_dati])

    pace = (pi-Omega*(simtime-simtime0))%(2*pi)
    pace1 = (pi-Omega*(simtime-simtime0))
    tace = 90.-ACE_theta[aceind]

    # look for the corresponding STA Record
    staind = argmin([abs((ddd-d).total_seconds()) for ddd in STA_dati])
    psta = (pace1 - (ACE_phi[aceind]-STA_phi[staind])*pi/180.)%(2*pi)
    tsta = 90.-STA_theta[staind]

    # look for the corresponding STB Record
    stbind = argmin([abs((ddd-d).total_seconds()) for ddd in STB_dati])
    pstb = (pace1 - (ACE_phi[aceind]-STB_phi[stbind])*pi/180.)%(2*pi)
    tstb = 90.-STB_theta[stbind]

    # look for the corresponding MSGR record
    msgrind = argmin([abs((ddd-d).total_seconds()) for ddd in MSGR_dati])
    pmsgr = (pace1 - (ACE_phi[aceind]-MSGR_phi[msgrind])*pi/180.)%(2*pi)
    tmsgr = 90.-MSGR_theta[msgrind]

    islice=argmin(abs(Rc-ACE_r[aceind]*215.))   # index of the desired radial distance
    axs=summary.plots(f,islice,1.e25,1.5)

    for ax in axs:
        ax.scatter(pace*180/pi,tace,color='black',edgecolor='black',s=20,zorder=10)
        ax.scatter(pmsgr*180/pi,tmsgr,color='purple',edgecolor='black',s=20,zorder=10)
        ax.scatter(psta*180/pi,tsta,color='red',edgecolor='black',s=20,zorder=10)
        ax.scatter(pstb*180/pi,tstb,color='dodgerblue',edgecolor='black',s=20,zorder=10)

    plt.tight_layout(h_pad=0.25)
    fig.suptitle(d.strftime('%Y:%m:%d %H:%M'),y=1)
#    axs[0].set_title(d.strftime('%Y:%m:%d %H:%M'))
    if save_pdf:
        plt.savefig(lfmFileDir+'/figs/islice/summary_%05d.pdf'%i,bbox_inches='tight',dpi=300)
    else:
        plt.savefig(lfmFileDir+'/figs/islice/summary_%05d.png'%i,bbox_inches='tight',dpi=300)
    plt.clf()



    
