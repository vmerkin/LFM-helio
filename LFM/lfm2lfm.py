from numpy import linspace,pi,meshgrid,sin,cos,zeros,ones,dstack,diff,sqrt,array,savetxt,flatnonzero,insert,asarray,zeros_like,arccos,arctan2
import os,sys,glob
from scipy import interpolate
import time
import lfm2lfm
#sys.path.append('../lib')
import pyLTR

#----------- PARSE ARGUMENTS ---------#
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('ConfigFileName',help='The name of the configuration file to use',default='startup.config')
args = parser.parse_args()
#----------- PARSE ARGUMENTS ---------#

# Read params from config file
prm = lfm2lfm.params.params(args.ConfigFileName)
(ni,nj,nk) = (prm.ni,prm.nj,prm.nk)
gamma = prm.gamma

b0 = prm.bscale
n0 = prm.nscale
l0 = prm.xscale
v0 = b0*1.e12/sqrt(1.67*n0)
t0 = l0/v0
mp = 1.67e-24
islice = -10


############### inLFM STUFF #####################
from pyhdf.SD import SD,SDC
inLFM = SD(prm.inFile)
vx = inLFM.select('vx_').get()[:-1,:-1,islice]/v0
vy = inLFM.select('vy_').get()[:-1,:-1,islice]/v0
vz = inLFM.select('vz_').get()[:-1,:-1,islice]/v0
bx = inLFM.select('bx_').get()[:-1,:-1,islice]/b0
by = inLFM.select('by_').get()[:-1,:-1,islice]/b0
bz = inLFM.select('bz_').get()[:-1,:-1,islice]/b0
rho = inLFM.select('rho_').get()[:-1,:-1,islice]/n0/mp
cs  = inLFM.select('c_').get()[:-1,:-1,islice]/v0
bi  = inLFM.select('bi_').get()[:-1,:-1,islice]/b0/l0**2

x = inLFM.select('X_grid').get()
y = inLFM.select('Y_grid').get()
z = inLFM.select('Z_grid').get()
# =========== Cell centers ==============
xc=0.125*( x[:-1,:-1,:-1]+x[1:,:-1,:-1]+x[:-1,1:,:-1]+x[:-1,:-1,1:]+
           x[1:,1:,:-1]+x[1:,:-1,1:]+x[:-1,1:,1:]+x[1:,1:,1:] )
yc=0.125*( y[:-1,:-1,:-1]+y[1:,:-1,:-1]+y[:-1,1:,:-1]+y[:-1,:-1,1:]+
           y[1:,1:,:-1]+y[1:,:-1,1:]+y[:-1,1:,1:]+y[1:,1:,1:] )
zc=0.125*( z[:-1,:-1,:-1]+z[1:,:-1,:-1]+z[:-1,1:,:-1]+z[:-1,:-1,1:]+
           z[1:,1:,:-1]+z[1:,:-1,1:]+z[:-1,1:,1:]+z[1:,1:,1:] )
# =======================================
Rc=sqrt(xc[:,:,islice+1]**2+yc[:,:,islice+1]**2+zc[:,:,islice+1]**2)
thetac=arccos(zc[:,:,islice+1]/Rc)
phic=arctan2(yc[:,:,islice+1],xc[:,:,islice+1])
phic[phic<0]+=2*pi

br     = bx*cos(phic)*sin(thetac) + by*sin(phic)*sin(thetac) + bz*cos(thetac)
vr     = vx*cos(phic)*sin(thetac) + vy*sin(phic)*sin(thetac) + vz*cos(thetac)

#br = sqrt(bx**2+by**2+bz**2)
#vr = sqrt(vx**2+vy**2+vz**2)
inLFM.end()

# ############### LFM STUFF #####################

# LFM GRID
sg = pyLTR.Grids.SphereGrid((ni,nj,nk))
(P,T,R) = sg.ptrCorner(rmin=prm.rmin,rmax=prm.rmax,thetamin=prm.thetamin)
(Pc,Tc,Rc) = sg.ptrCenter()
(x,y,z) = sg.xyzCorner()
(xc,yc,zc) = sg.xyzCenter()

# ############### LFM STUFF #####################
if prm.dumpInit:
    lfmh = pyLTR.Tools.lfmstartup.lfmstartup(os.path.join(prm.dirInitLFMfile,prm.initLFMfile),(ni,nj,nk))
    lfmh.open()
    # Set variables
    lfmh.writeVar('X_grid',x)
    lfmh.writeVar('Y_grid',y)
    lfmh.writeVar('Z_grid',z)
        
    lfmh.writeVar('rho_',rho[:,:,None]*(R[0,0,0]/Rc)**2)
    lfmh.writeVar('c_',  cs[:,:,None]*ones((nk,nj,ni)))
    lfmh.writeVar('vx_',vr[:,:,None]*sin(Tc)*cos(Pc))
    lfmh.writeVar('vy_',vr[:,:,None]*sin(Tc)*sin(Pc))
    lfmh.writeVar('vz_',vr[:,:,None]*cos(Tc))

#    lfmh.writeVar('vx_',vx[:,:,None]*ones((nk,nj,ni)))
#    lfmh.writeVar('vy_',vy[:,:,None]*ones((nk,nj,ni)))
#    lfmh.writeVar('vz_',vz[:,:,None]*ones((nk,nj,ni)))


    # note the fancy dstack to fill in the last i-face with correct bi values
    #
    # note we use R[0,0,0] instead of Rc
    # we thus assume that WSA data are at rmin as above for the density
    bi = dstack(R.shape[2]*[bi])
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

# DEAL WITH UNITS HERE!!!

    bp_kface = -br_kface*2*pi/(prm.Tsolar*24.*3600./t0)*R[0,0,0]*sin(Tc[0,:,0])/vr_kface  # beautiful numpy broadcasting
    bp       = -br      *2*pi/(prm.Tsolar*24.*3600./t0)*R[0,0,0]*sin(Tc[0,:,0])/vr          
#    et   = -br_kface*2*pi/(prm.Tsolar*24.*3600.)*R[0,0,0]*sin(Tc[0,:,0])   # etheta electric field. Note, only defined in 2D

     # Scale inside ghost region
    (vr,rho,cs,br,bp,bp_kface) = [dstack(prm.NO2*[var]) for var in (vr,rho,cs,br,bp,bp_kface)]
    rho*=(R[0,0,0]/rcg)**2
    br*=(R[0,0,0]/rcg)**2
    bp*=(R[0,0,0]/rcg)
    bp_kface*=(R[0,0,0]/rcg)
    
    bt = zeros_like(bp)
    vp = zeros_like(bp)
    vt = zeros_like(bp)
    et = zeros_like(bp)
    ep = zeros_like(bp)
     

    for i in range(prm.NO2):
        # note, by indexing below, we're making sure we are dumping arrays with shape (nk,nj)
        # instead of thinking whether they are already of this size
        out = array([ bt[:nk,:nj,i].T.ravel(),              # Note, bt (1st column) is a place holder used only for tiem-dependent runs (should be bt_jface)
                      bp_kface[:nk,:nj,i].T.ravel(),        # the 1st two columns are used for e-field calculation; the next 3 for bx,by,bz calc in ghost region
                      et[:nk,:nj,i].T.ravel(),
                      ep[:nk,:nj,i].T.ravel(),
                      br[:nk,:nj,i].T.ravel(),
                      bt[:nk,:nj,i].T.ravel(),
                      bp[:nk,:nj,i].T.ravel(),
                      vr[:nk,:nj,i].T.ravel(),
                      vt[:nk,:nj,i].T.ravel(),
                      vp[:nk,:nj,i].T.ravel(),
                      rho[:nk,:nj,i].T.ravel(),
                      cs[:nk,:nj,i].T.ravel()])
                  
        savetxt('innerbc_001_%d.dat'%i,out.T,
                # fmt=['%13.8f','%13.8f','%15.8f','%13.8f','%13.8f','%13.8f',
                #      '%15.5e','%15.5e','%15.5e',
                #      '%14.5e','%14.5e'],
                fmt=['%15.5e']*12,
                delimiter='')

# """
# if prm.plots:

#     import matplotlib.pyplot as plt
#     plt.figure(); plt.pcolormesh(phi_wsa_v,theta_wsa_v,bi_wsa); plt.colorbar()
#     plt.figure(); plt.pcolormesh(phi_lfm_v,theta_lfm_v,bi_lfm); plt.colorbar()
#     show()

#     if smooth: 
#         figure(); pcolormesh(phi_lfm_v,theta_lfm_v,bi_lfm_smooth); colorbar()
#         contour(phi_lfm_c,theta_lfm_c,bi_lfm_smooth,[0.])
#         contour(phi_lfm_v[:-1],theta_lfm_c,bk_lfm_smooth,[0.],colors='purple',linewidth=5)
# """    




