from numpy import linspace,pi,meshgrid,sin,cos,zeros,ones,dstack,diff,sqrt,array,savetxt,flatnonzero,insert,asarray,zeros_like
import os,sys,glob
from scipy import interpolate
import time
import wsa2lfm
sys.path.append('../lib')
import wsa
import pyLTR

#----------- PARSE ARGUMENTS ---------#
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('ConfigFileName',help='The name of the configuration file to use',default='startup.config')
args = parser.parse_args()
#----------- PARSE ARGUMENTS ---------#

# Read params from config file
prm = wsa2lfm.params.params(args.ConfigFileName)
(ni,nj,nk) = (prm.ni,prm.nj,prm.nk)
gamma = prm.gamma


############### WSA STUFF #####################
phi_wsa_v,theta_wsa_v,phi_wsa_c,theta_wsa_c,bi_wsa,v_wsa,n_wsa,T_wsa = wsa.read(prm.wsaFile,prm.densTempInfile)
############### WSA STUFF #####################

############### LFM STUFF #####################

# LFM GRID
sg = pyLTR.Grids.SphereGrid((ni,nj,nk))
(P,T,R) = sg.ptrCorner(rmin=prm.rmin*prm.scale,rmax=prm.rmax*prm.scale,thetamin=prm.thetamin)
(Pc,Tc,Rc) = sg.ptrCenter()
(x,y,z) = sg.xyzCorner()
(xc,yc,zc) = sg.xyzCenter()

# this is fast and better than griddata in that it nicely extrapolates boundaries:
fbi      = interpolate.RectBivariateSpline(phi_wsa_c,theta_wsa_c,bi_wsa.T,kx=1,ky=1)  
br = fbi(Pc[:,0,0],Tc[0,:,0])


############### SMOOTHING #####################
if not prm.gaussSmoothWidth==0:
    import astropy
    from astropy.convolution import convolve,Gaussian2DKernel

    gauss=Gaussian2DKernel(width=prm.gaussSmoothWidth)
    br   =astropy.convolution.convolve(br,gauss,boundary='extend')

fv      = interpolate.RectBivariateSpline(phi_wsa_c,theta_wsa_c,v_wsa.T,kx=1,ky=1)  
vr = fv(Pc[:,0,0],Tc[0,:,0])

f      = interpolate.RectBivariateSpline(phi_wsa_c,theta_wsa_c,n_wsa.T,kx=1,ky=1)  
rho = f(Pc[:,0,0],Tc[0,:,0])

f      = interpolate.RectBivariateSpline(phi_wsa_c,theta_wsa_c,T_wsa.T,kx=1,ky=1)  
temp = f(Pc[:,0,0],Tc[0,:,0])


cs = sqrt(gamma*1.38e-23*temp/1.67e-27)*1.e2   # FIX ME. in cm/s

############### LFM STUFF #####################
if prm.dumpInit:
    lfmh = pyLTR.Tools.lfmstartup.lfmstartup(os.path.join(prm.dirInitLFMfile,prm.initLFMfile),(ni,nj,nk))
    lfmh.open()
    # Set variables
    lfmh.writeVar('X_grid',x)
    lfmh.writeVar('Y_grid',y)
    lfmh.writeVar('Z_grid',z)
        
    # lfmh.setVar('rho_',400.*1.67e-24*ones((nk+1,nj+1,ni+1)))
    # lfmh.setVar('c_',  5.e6*ones((nk+1,nj+1,ni+1)))
    # lfmh.setVar('vx_',4.e7*sin(Tc)*cos(Pc))
    # lfmh.setVar('vy_',4.e7*sin(Tc)*sin(Pc))
    # lfmh.setVar('vz_',4.e7*cos(Tc))

    lfmh.writeVar('rho_',rho[:,:,None]*(R[0,0,0]/Rc)**2)  # note, we map to the bottom boundary, to the first cell center. Could use Rc in the numerator.
    lfmh.writeVar('c_',  cs[:,:,None]*ones((nk,nj,ni)))
    lfmh.writeVar('vx_',vr[:,:,None]*sin(Tc)*cos(Pc))
    lfmh.writeVar('vy_',vr[:,:,None]*sin(Tc)*sin(Pc))
    lfmh.writeVar('vz_',vr[:,:,None]*cos(Tc))


    tmp    = diff(T[:,:,0],axis=1)
    dtheta = 0.5*(tmp[:-1,:]+tmp[1:,:])
    tmp    = diff(P[:,:,0],axis=0)
    dphi   = 0.5*(tmp[:,:-1]+tmp[:,1:])

    # note the fancy dstack to fill in the last i-face with correct bi values
    #
    # note we use R[0,0,0] instead of Rc
    # we thus assume that WSA data are at rmin as above for the density
    bi = dstack(R.shape[2]*[br*R[0,0,0]**2*sin(Tc[:,:,0])*dtheta*dphi])
    lfmh.writeVar('bi_',bi)
    lfmh.close()


if prm.dumpBC:
    # ghost cell coordinates
    (xg,yg,zg,xcg,ycg,zcg) = [2*arr[:,:,[0]]-arr[:,:,1:prm.NO2+1] for arr in (x,y,z,xc,yc,zc)]
    # Here we assume grid radii don't vary with theta and phi, at least in the ghost region
    rcg = sqrt(xcg[0,0,:]**2+ycg[0,0,:]**2+zcg[0,0,:]**2)
    rg = sqrt(xg[0,0,:]**2+yg[0,0,:]**2+zg[0,0,:]**2)

    # note, redefining interpolation functions we could also
    # interpolate from bi_wsa as above, but then we would have to
    # smooth bk, if necessary. The way we're doing it here, bk will be
    # smoothed or not, dependent on whether br has been smoothed.
    # note also, this has to extrapolate
    fbi = interpolate.RectBivariateSpline(Pc[:,0,0],Tc[0,:,0],br,kx=1,ky=1)
    fv  = interpolate.RectBivariateSpline(Pc[:,0,0],Tc[0,:,0],vr,kx=1,ky=1)
    br_kface = fbi(P[:,0,0],Tc[0,:,0])
    vr_kface  = fv (P[:,0,0],Tc[0,:,0])

    bp_kface = -br_kface*2*pi/(prm.Tsolar*24.*3600.)*R[0,0,0]*sin(Tc[0,:,0])/vr_kface  # beautiful numpy broadcasting
    bp       = -br      *2*pi/(prm.Tsolar*24.*3600.)*R[0,0,0]*sin(Tc[0,:,0])/vr          

    # Scale inside ghost region
    (vr,rho,cs,br,bp,bp_kface) = [dstack(prm.NO2*[var]) for var in (vr,rho,cs,br,bp,bp_kface)]
    rho*=(R[0,0,0]/rcg)**2
    br*=(R[0,0,0]/rcg)**2
    bp*=(R[0,0,0]/rcg)
    bp_kface*=(R[0,0,0]/rcg)

    bt = zeros_like(bp)
    vp = zeros_like(bp)
    vt = zeros_like(bp)

    for i in range(prm.NO2):
        # note, by indexing below, we're making sure we are dumping arrays with shape (nk,nj)
        # instead of thinking whether they are already of this size
        out = array([ bt[:nk,:nj,i].T.ravel(),              # Note, bt (1st column) is a place holder used only for tiem-dependent runs (should be bt_jface)
                      bp_kface[:nk,:nj,i].T.ravel(),        # the 1st two columns are used for e-field calculation; the next 3 for bx,by,bz calc in ghost region
                      br[:nk,:nj,i].T.ravel(),
                      bt[:nk,:nj,i].T.ravel(),
                      bp[:nk,:nj,i].T.ravel(),
                      vr[:nk,:nj,i].T.ravel(),
                      vt[:nk,:nj,i].T.ravel(),
                      vp[:nk,:nj,i].T.ravel(),
                      rho[:nk,:nj,i].T.ravel(),
                      cs[:nk,:nj,i].T.ravel()])
                  
        savetxt('innerbc_000_%d.dat'%i,out.T,
                fmt=['%13.8f','%13.8f','%13.8f','%13.8f','%13.8f',
                     '%15.5e','%15.5e','%15.5e',
                     '%14.5e','%14.5e'],
                delimiter='')

"""
if prm.plots:

    import matplotlib.pyplot as plt
    plt.figure(); plt.pcolormesh(phi_wsa_v,theta_wsa_v,bi_wsa); plt.colorbar()
    plt.figure(); plt.pcolormesh(phi_lfm_v,theta_lfm_v,bi_lfm); plt.colorbar()
    show()

    if smooth: 
        figure(); pcolormesh(phi_lfm_v,theta_lfm_v,bi_lfm_smooth); colorbar()
        contour(phi_lfm_c,theta_lfm_c,bi_lfm_smooth,[0.])
        contour(phi_lfm_v[:-1],theta_lfm_c,bk_lfm_smooth,[0.],colors='purple',linewidth=5)
"""    




