import matplotlib
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.dates import DateFormatter
import glob,os,sys,datetime
from numpy import pi,argmin,sin,cos,vstack,linspace,ones_like,arcsinh
import time

home = os.path.expanduser("~");
if os.path.join(home,'GitHub/Python/LFM-helio/lib/') not in sys.path: sys.path.append(os.path.join(home,'GitHub/Python/LFM-helio/lib/'))
if os.path.join(home,'GitHub/Python/LFM-helio/analysis/') not in sys.path: sys.path.append(os.path.join(home,'GitHub/Python/LFM-helio/analysis/'))
import lfmhlib,ace,cspice,stereo
import summary

####################################################################################################
############## SET UP PARAMETERS      ##############################################################
####################################################################################################
save_pdf = False




####################################################################################################
##############  THIS IS FOR 2008      ##############################################################
####################################################################################################
#lfmFileDir = '/glade/p/ujhb0003/LFM-helio/ADAPT/2066rf_fullpb_192x192x768_converged'
#lfmFileDir = '/glade/scratch/vgm/LFM-helio/ADAPT/2066rf_fullpb_192x192x768_converged-smoothed'
#lfmFileDir = '/glade/scratch/vgm/LFM-helio/ADAPT/2066rf_fullpb_192x192x768_converged_fixedboundaries'
lfmFileDir = '/glade/scratch/vgm/LFM-helio/ADAPT/2066rf_fullpb_192x192x768_converged_fixedboundaries_stairfix'
start_step = 40000 #64000
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
lfmFiles = sorted(glob.glob(os.path.join(lfmFileDir,'cr2066_192x192x768_mhd_*.hdf')))
start_ind = [start_step==int(f[-11:-4]) for f in lfmFiles].index(True)

plt.style.use('publication')
plt.close('all')

fig = plt.figure(figsize=(20/6.,16./6.)) 

for i,f in enumerate(lfmFiles[start_ind:]):
    if i==0: 
        simtime0 = lfmhlib.get_time(lfmFiles[0])
        R,theta,phi,Rc,thetac,phic = lfmhlib.r_theta_phi_uniform(f)
        nk = phic.size

#    if i<>1052: continue
    if i<>0: continue

    simtime = lfmhlib.get_time(f)

    d = start_dati+datetime.timedelta(seconds=simtime-simtime0)
    print(d)

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

    jsc=argmin(abs(pi/2.-thetac-tmsgr*pi/180.))
    br=lfmhlib.read_var_jslice(f,'br',jsc,thetac,phic)
    xc=lfmhlib.read_var_jslice(f,'X_grid',jsc,thetac,phic)/6.96e10
    yc=lfmhlib.read_var_jslice(f,'Y_grid',jsc,thetac,phic)/6.96e10

    # fix up periodic axis
    xc=vstack([xc,xc[0,:]])
    yc=vstack([yc,yc[0,:]])
    br=vstack([br,br[0,:]])

    # import astropy
    # from astropy.convolution import convolve,Gaussian2DKernel

    # gauss=Gaussian2DKernel(width=1)
    # br   =astropy.convolution.convolve(br,gauss,boundary='extend')


    ax = fig.add_subplot(111,aspect='equal')
    p=ax.pcolormesh(xc,yc,br*(Rc/Rc[0])**2*1.e5,vmin=-160,vmax=160,shading='gouraud',cmap=cm.RdBu_r)
#    p=ax.pcolormesh(xc,yc,br,shading='gouraud',cmap=cm.RdBu_r)
#    p=ax.contourf(xc,yc,br*(Rc/Rc[0])**2*1.e5,linspace(-160,160,30),extend='both',shading='gouraud',cmap=cm.RdBu_r)
    cb=plt.colorbar(p,ax=ax)
    cb.set_label(r'B$_r\times (r/r_0)^2$, nT')
    cb.solids.set_rasterized(True)
    ax.scatter(MSGR_r[msgrind]*215.*sin((90.-tmsgr)*pi/180.)*cos(pmsgr),
               MSGR_r[msgrind]*215.*sin((90.-tmsgr)*pi/180.)*sin(pmsgr),
               color='purple',edgecolor='black',s=10,zorder=10,linewidth=0.5)
    ax.scatter(ACE_r[aceind]*215.*sin((90.-tmsgr)*pi/180.)*cos(pace),
                ACE_r[aceind]*215.*sin((90.-tmsgr)*pi/180.)*sin(pace),
                color='black',edgecolor='black',s=10,zorder=10,linewidth=0.5)
    ax.scatter(ACE_r[staind]*215.*sin((90.-tmsgr)*pi/180.)*cos(psta),
                ACE_r[staind]*215.*sin((90.-tmsgr)*pi/180.)*sin(psta),
                color='red',edgecolor='black',s=10,zorder=10,linewidth=0.5)
    ax.scatter(ACE_r[stbind]*215.*sin((90.-tmsgr)*pi/180.)*cos(pstb),
                ACE_r[stbind]*215.*sin((90.-tmsgr)*pi/180.)*sin(pstb),
                color='dodgerblue',edgecolor='black',s=10,zorder=10,linewidth=0.5)

    ax.set_xlim(-230,230)
    ax.set_ylim(-230,230)
    ax.set_xlabel('X, Rs')
    ax.set_ylabel('Y, Rs')

#    plt.tight_layout(h_pad=0.25)
    fig.suptitle(d.strftime('%Y:%m:%d %H:%M'),y=1)
#    axs[0].set_title(d.strftime('%Y:%m:%d %H:%M'))
    if save_pdf:
        plt.savefig(lfmFileDir+'/figs/jslice/summary__%05d.pdf'%i,bbox_inches='tight',dpi=300)
    else:
        plt.savefig(lfmFileDir+'/figs/jslice/summary_%05d.png'%i,bbox_inches='tight',dpi=300)
    plt.clf()



    
