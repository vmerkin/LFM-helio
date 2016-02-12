from pylab import *
import glob,os,sys,datetime
if len(sys.argv)!=2: sys.exit('Please, provide the spacecraft name (ACE, STA, STB, MSGR) and this should be the only parameter')
if sys.argv[1] not in ['ACE','STA','STB','MSGR']: sys.exit('Spacecraft name should be one of these: (ACE, STA, STB, MSGR)')

####################################################################################################
############## SET UP PARAMETERS      ##############################################################
####################################################################################################
SC = sys.argv[1] #'ACE'  'STA', 'STB', 'MSGR'     # which spacecraft to map

####################################################################################################
##############  THIS IS FOR 2008      ##############################################################
####################################################################################################
#lfmFileDir = '/glade/p/ujhb0003/LFM-helio/ADAPT/2066rf_fullpb_192x192x768_converged'
#'/glade/p/ujhb0003/LFM-helio/ADAPT/2066rf'
lfmFileDir = '/glade/scratch/vgm/LFM-helio/ADAPT/2066rf_fullpb_192x192x768_converged_fixedboundaries_stairfix'
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
dates = []
vr = []
br = []
n  = []
t  = []
vr0= []
br0= []
n0 = []
t0 = []
i_ind = []
j_ind = []
k_ind = []
for i,f in enumerate(lfmFiles[start_ind:]):
    if i==0: 
        simtime0 = lfmhlib.get_time(f)
        R,theta,phi,Rc,thetac,phic = lfmhlib.r_theta_phi_uniform(f)
        nk = phic.size
        
    simtime = lfmhlib.get_time(f)

    d = start_dati+datetime.timedelta(seconds=simtime-simtime0)
    dates.append(d)
    print(d)

    # look for the corresponding ACE record
    aceind = argmin([abs((ddd-d).total_seconds()) for ddd in ACE_dati])

    pace = (pi-Omega*(simtime-simtime0))%(2*pi)
    pace1 = (pi-Omega*(simtime-simtime0))
    tace = 90.-ACE_theta[aceind]

    if SC=='ACE':
        isc=argmin(abs(Rc-ACE_r[aceind]*215.))   # index of the desired radial distance
        jsc=argmin(abs(pi/2.-thetac-tace*pi/180.))   # index of the desired theta-angle
        ksc = int(floor(pace/(2*pi)*nk))
    elif SC=='STA':
        # look for the corresponding STA Record
        staind = argmin([abs((ddd-d).total_seconds()) for ddd in STA_dati])
        psta = (pace1 - (ACE_phi[aceind]-STA_phi[staind])*pi/180.)%(2*pi)
        tsta = 90.-STA_theta[staind]

        isc=argmin(abs(Rc-STA_r[staind]*215.))   # index of the desired radial distance
        jsc=argmin(abs(pi/2.-thetac-tsta*pi/180.))   # index of the desired theta-angle
        ksc = int(floor(psta/(2*pi)*nk))
    elif SC=='STB':
        # look for the corresponding STB Record
        stbind = argmin([abs((ddd-d).total_seconds()) for ddd in STB_dati])
        pstb = (pace1 - (ACE_phi[aceind]-STB_phi[stbind])*pi/180.)%(2*pi)
        tstb = 90.-STB_theta[stbind]

        isc=argmin(abs(Rc-STB_r[stbind]*215.))   # index of the desired radial distance
        jsc=argmin(abs(pi/2.-thetac-tstb*pi/180.))   # index of the desired theta-angle
        ksc = int(floor(pstb/(2*pi)*nk))
    elif SC=='MSGR':
        # look for the corresponding MSGR record
        msgrind = argmin([abs((ddd-d).total_seconds()) for ddd in MSGR_dati])
        pmsgr = (pace1 - (ACE_phi[aceind]-MSGR_phi[msgrind])*pi/180.)%(2*pi)
        tmsgr = 90.-MSGR_theta[msgrind]

        isc=argmin(abs(Rc-MSGR_r[msgrind]*215.))   # index of the desired radial distance
        jsc=argmin(abs(pi/2.-thetac-tmsgr*pi/180.))   # index of the desired theta-angle
        ksc = int(floor(pmsgr/(2*pi)*nk))

    v = lfmhlib.read_var_point(f,'vr',isc,jsc,ksc,thetac,phic)
#    vx = lfmhlib.read_var_point(f,'vx',isc,jsc,ksc,thetac,phic)
#    vy = lfmhlib.read_var_point(f,'vy',isc,jsc,ksc,thetac,phic)
#    vz = lfmhlib.read_var_point(f,'vz',isc,jsc,ksc,thetac,phic)
#    v = sqrt(vx**2+vy**2+vz**2)
    b = lfmhlib.read_var_point(f,'br',isc,jsc,ksc,thetac,phic)
    nn = lfmhlib.read_var_point(f,'rho',isc,jsc,ksc,thetac,phic)/1.67e-24
    tt = lfmhlib.read_var_point(f,'c',isc,jsc,ksc,thetac,phic)**2*1.67e-8/1.38/1.5
    vr.append(v)
    br.append(b)
    n.append(nn)
    t.append(tt)


    v= lfmhlib.read_var_point(lfmFiles[start_ind],'vr',isc,jsc,ksc,thetac,phic)
#    vx = lfmhlib.read_var_point(lfmFiles[start_ind],'vx',isc,jsc,ksc,thetac,phic)
#    vy = lfmhlib.read_var_point(lfmFiles[start_ind],'vy',isc,jsc,ksc,thetac,phic)
#    vz = lfmhlib.read_var_point(lfmFiles[start_ind],'vz',isc,jsc,ksc,thetac,phic)
#    v = sqrt(vx**2+vy**2+vz**2)
    b= lfmhlib.read_var_point(lfmFiles[start_ind],'br',isc,jsc,ksc,thetac,phic)
    nn = lfmhlib.read_var_point(lfmFiles[start_ind],'rho',isc,jsc,ksc,thetac,phic)/1.67e-24
    tt = lfmhlib.read_var_point(lfmFiles[start_ind],'c',isc,jsc,ksc,thetac,phic)**2*1.67e-8/1.38/1.5
    vr0.append(v)
    br0.append(b)
    n0.append(nn)
    t0.append(tt)

    i_ind.append(isc)
    j_ind.append(jsc)
    k_ind.append(ksc)

out = {'Time':{'data':dates,'meta':''},
       'vr0':{'data':array(vr0),'meta':'Vr[cm/s] time series extracted from the relaxed start file'},
       'vr':{'data':array(vr),'meta':'Vr[cm/s] time series extracted from the time-dependent files'},
       'br0':{'data':array(br0),'meta':'Br[Gs] time series extracted from the relaxed start file'},
       'br':{'data':array(br),'meta':'Br[Gs] time series extracted from the time-dependent files'},
       'n0':{'data':array(n0),'meta':'n[/cc] time series extracted from the relaxed start file'},
       'n':{'data':array(n),'meta':'n[/cc] time series extracted from the time-dependent files'},
       'T0':{'data':array(t0),'meta':'T[K] time series extracted from the relaxed start file'},
       'T':{'data':array(t),'meta':'T[K] time series extracted from the time-dependent files'},
       'i index list':{'data':i_ind,'meta':'list of i indices marking the sc position in R'},
       'j index list':{'data':j_ind,'meta':'list of j indices marking the sc position in theta'},
       'k index list':{'data':k_ind,'meta':'list of k indices marking the sc position in phi'}}
import pickle
pickle.dump(out,open(SC+'_cr2066rf_192x192x768.pkl','wb'))
#pickle.dump(out,open(SC+'_cr2094rf_192x192x768.pkl','wb'))
