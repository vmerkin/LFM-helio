from numpy import linspace,pi,meshgrid,sin,cos,zeros,ones,dstack,diff,sqrt,array,savetxt,flatnonzero,insert,asarray,zeros_like,where,roll,ediff1d,empty,tanh,insert,append
import os,sys,glob
from scipy import interpolate
from scipy.optimize import newton_krylov,anderson

import time
import adapt2lfm
if '../lib' not in sys.path: sys.path.append('../lib')
import wsa
import poisson
import pyLTR

#----------- PARSE ARGUMENTS ---------#
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('ConfigFileName',help='The name of the configuration file to use',default='startup.config')
args = parser.parse_args()
#----------- PARSE ARGUMENTS ---------#

# Read params from config file
prm = adapt2lfm.params.params(args.ConfigFileName)
(ni,nj,nk) = (prm.ni,prm.nj,prm.nk)




wsaFiles = sorted(glob.glob(os.path.join(prm.adaptdir,prm.adaptWildCard)))

et_save = zeros( (nj,nk+1) )
ep_save = zeros( (nj+1,nk) )

for (fcount,wsaFile) in enumerate(wsaFiles):
    print(fcount)
    ############### WSA STUFF #####################
    if wsaFile==wsaFiles[0]: 
        verbose=True
    else:
        verbose=False

    phi_wsa_v,theta_wsa_v,phi_wsa_c,theta_wsa_c,bi_wsa,v_wsa,n_wsa,T_wsa = wsa.read(wsaFile,prm.densTempInfile,prm.normalized,verbose=verbose)
    ############### WSA STUFF #####################


    if wsaFile == wsaFiles[0]:
        ############### LFM STUFF #####################
        # LFM GRID
        sg = pyLTR.Grids.SphereGrid((ni,nj,nk))
        (P,T,R) = sg.ptrCorner(rmin=prm.rmin*prm.scale,rmax=prm.rmax*prm.scale,thetamin=prm.thetamin)
        (Pc,Tc,Rc) = sg.ptrCenter()
        (x,y,z) = sg.xyzCorner()
        (xc,yc,zc) = sg.xyzCenter()
        phi = Pc[:,0,0]
        theta = Tc[0,:,0]
        pois = poisson.poisson(theta,phi)
    ####################################################################################################
    #
    # We do this below to align WSA solutions with the very first file
    # -- i.e., rotate into the rotating frame
    #
    # we need the sidereal
    # period here in the denominator, but we want to make it general
    # enough that Tsolar is in arbitrary time units. Thus we just
    # multiply by the factor in days, leaving Tsolar alone.

    omega=2*pi/prm.Tsolar*(25.38/27.27)
    phi_prime=(phi_wsa_c-omega*prm.adaptCadence*fcount)%(2*pi)
    if where(ediff1d(phi_prime)<0)[0].size!=0:
        ind0=where(ediff1d(phi_prime)<0)[0][0]+1
    else:
        ind0=0
    phi_prime=roll(phi_prime,-ind0)
    bi_wsa_rolled=roll(bi_wsa,-ind0,axis=1)
    v_wsa_rolled=roll(v_wsa,-ind0,axis=1)
    n_wsa_rolled=roll(n_wsa,-ind0,axis=1)
    T_wsa_rolled=roll(T_wsa,-ind0,axis=1)
    # ####################################################################################################

    # this is fast and better than griddata in that it nicely extrapolates boundaries:
    fbi      = interpolate.RectBivariateSpline(phi_wsa_c,theta_wsa_c,bi_wsa_rolled.T,kx=1,ky=1)  
    br = fbi(phi,theta)

    ############### SMOOTHING #####################
    if not prm.gaussSmoothWidth==0:
        import astropy
        from astropy.convolution import convolve,Gaussian2DKernel

        gauss=Gaussian2DKernel(width=prm.gaussSmoothWidth)
        br   =astropy.convolution.convolve(br,gauss,boundary='extend')


    fv      = interpolate.RectBivariateSpline(phi_wsa_c,theta_wsa_c,v_wsa_rolled.T,kx=1,ky=1)  
    vr = fv(Pc[:,0,0],Tc[0,:,0])

    f      = interpolate.RectBivariateSpline(phi_wsa_c,theta_wsa_c,n_wsa_rolled.T,kx=1,ky=1)  
    rho = f(Pc[:,0,0],Tc[0,:,0])

    f      = interpolate.RectBivariateSpline(phi_wsa_c,theta_wsa_c,T_wsa_rolled.T,kx=1,ky=1)  
    temp = f(Pc[:,0,0],Tc[0,:,0])


    cs = sqrt(prm.gamma*1.38e-23*temp/1.67e-27)*1.e2   # FIX ME. in cm/s


    if wsaFile == wsaFiles[0]:
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


    # Poisson solve after interpolation onto LFM grid
    # This way we get better accuracy when LFM grid is finer than WSA
    #
    # Note also, the poisson function is written assuming dimensions
    # are (nj,nk) while to conform to the LFM grid allignement
    # (nk,nj,ni), we have transposed the WSA arrays above so br by now
    # is (nk,nj)
    if fcount>0:
        pois.setRHS( (br-br_save).T )
        guess=zeros_like(br.T)
        Psi = newton_krylov(pois.residual,guess, method='lgmres')#,f_rtol=1.e-6) #iter=100
        print('Residual: %g' % abs(pois.residual(Psi)).max())


        """
        Now we have our electric field along edges -- just like LFM
        wants the suffix _a denotest that this is adapt field, as
        opposed to corotation field et that we calculate below
        """
        et_a = zeros( (Psi.shape[0],Psi.shape[1]+1) )
        et_a[:,1:-1] = diff(Psi,axis=1)/diff(phi)
        et_a[:,0] = (Psi[:,0] - Psi[:,-1])/(phi[0]-phi[-1]+2*pi); et_a[:,-1]=et_a[:,0]
        et_a/= sin(theta[:,None]) 
        """
        note, assuming theta constant along phi and same theta on the
        boundary and in the center of the LFM cell
        """
        ep_a=zeros((Psi.shape[0]+1,Psi.shape[1]))
        ep_a[1:-1,:] = -diff(Psi,axis=0)/diff(theta)[:,None]
        ep_a[0,:]=ep_a[1,:]   # used to set these to zero, but more appropriate to repeat from next theta, since Ephi does not depend on theta at the pole.
        ep_a[-1,:]=ep_a[-2,:] # Ek(j=0,njp1) is set to zero anyway inside LFM (because edge length is zero)

        """
        OK, now we do the trick with e-field definition as a function
        of time. We assume that e-field is changing linearly with time
        (which means Br is changing quadratically.  If we denote the
        E-field calculated above as E(T+1/2) (i.e., defined in between
        times T when Br is defined by ADAPT), the simple integration gives

        Q(T+1)+Q(T) = 2*dB/dT, where Q is curl(E) and dB/dT is the RHS
        of the equation above. Assuming E=0 at T=0, we then get

        Q(0) = 0
        Q(1) = 2*dB/dT(1/2)
        Q(2) = 2*dB/dT(3/2)-dB/dT(1)
        Q(3) = 2*dB/dT(5/2)-dB/dT(2)
        ...

        """
        # Convert to CGS. FIX ME!!! UNITS HARD CODED
        et_a*= prm.rmin*prm.scale/prm.adaptCadence/24./3600.
        ep_a*= prm.rmin*prm.scale/prm.adaptCadence/24./3600.

#        et_save = 2*et_a - et_save
#        ep_save = 2*ep_a - ep_save

        # ramp = lambda t: 0.5*(tanh(20/pi*(t-pi/10))-tanh(20/pi*(t-pi+pi/10)))
        # et_save = et_a*ramp(theta[:,None])
        # q=0.5*(theta[1:]+theta[:-1])
        # q=insert(q,0,0.)
        # q=append(q,pi)
        # ep_save = ep_a*ramp(q[:,None])

        et_save = et_a
        ep_save = ep_a

    br_save = br

    """
     The potential derived above is defined at LFM grid cell
     centers. For LFM-electric field, however, we need to interpolate
     it to cell vortices:
    """
#    Psiv=zeros((Psi.shape[0]+1,Psi.shape[1]+1))
#    Psiv[1:-1,1:-1] = 0.25*(Psi[:-1,:-1]+Psi[1:,:-1]+Psi[:-1,1:]+Psi[1:,1:])
    """
    Fix up periodic boundary
    """
#    Psiv[1:-1,0] = 0.25*(Psi[:-1,0]+Psi[1:,0]+Psi[:-1,-1]+Psi[1:,-1]); Psiv[1:-1,-1]=Psiv[1:-1,0]
    
    """
    Note, the above code fixes the potential at 0 on the pole. Our
    boundary condition on the Poisson solve assumed zero potential
    at the pole, so this is consistent.
    """


    """
    The above are our electric fields defined at times of ADAPT
    cadence. We now need the mag field corrections at the same
    time. Note the interpolation of the e-field to cell-centers.
    """

    # Note the ugliness with the transposes. This is all a result of
    # our trying to conform to the LFM array dimensions (nk,nj,ni) and
    # keep everything consistent with the wsa2lfm script
    vrt = vr.T
    bp_a = zeros_like(vrt)
    bt_a = zeros_like(vrt)
    bp_kface_a = zeros( (vrt.shape[0],vrt.shape[1]+1) )
    bt_jface_a = zeros( (vrt.shape[0]+1,vrt.shape[1]) )
    vrt_kface  = zeros( (vrt.shape[0],vrt.shape[1]+1) )
    vrt_jface  = zeros( (vrt.shape[0]+1,vrt.shape[1]) )
    if fcount >0:
        bp_a = 0.5*(et_save[:,:-1]+et_save[:,1:])/vrt
        bt_a = -0.5*(ep_save[:-1,:]+ep_save[1:,:])/vrt

        # the above are at cell centers, also need at the
        # corresponding faces
        #
        # First interpolate velocity to the edges
        vrt_kface[:,1:-1] = 0.5*(vrt[:,:-1]+vrt[:,1:]); vrt_kface[:,0] = 0.5*(vrt[:,-1]+vrt[:,0]); vrt_kface[:,-1] = vrt_kface[:,0]
        vrt_jface[1:-1,:] = 0.5*(vrt[1:,:]+vrt[:-1,:]) ; vrt_jface[0,:]=vrt[1,:].mean(); vrt_jface[-1,:]=vrt[-2,:].mean(); 
        
        bp_kface_a = et_save/vrt_kface
        bt_jface_a = -ep_save/vrt_jface
    # Note, these are defined at cell centers on the boundary (at rmin)
    bp_a = bp_a.T  
    bt_a = bt_a.T
    bt_jface_a = bt_jface_a.T
    bp_kface_a = bp_kface_a.T


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

        bp_kface = bp_kface_a
        #-br_kface*2*pi/(prm.Tsolar*24.*3600.)*R[0,0,0]*sin(Tc[0,:,0])/vr_kface + bp_kface_a # beautiful numpy broadcasting
        bp       = bp_a
        #-br      *2*pi/(prm.Tsolar*24.*3600.)*R[0,0,0]*sin(Tc[0,:,0])/vr       + bp_a
        et   = et_save.T
        #-br_kface*2*pi/(prm.Tsolar*24.*3600.)*R[0,0,0]*sin(Tc[0,:,0])              + et_save.T  # etheta electric field. Note, only defined in 2D


        # Scale inside ghost region
        (vr,rho,cs,br,bp,bp_kface) = [dstack(prm.NO2*[var]) for var in (vr,rho,cs,br,bp,bp_kface)]
        rho*=(R[0,0,0]/rcg)**2
        br*=(R[0,0,0]/rcg)**2
#        bp*=(R[0,0,0]/rcg)
#        bp_kface*=(R[0,0,0]/rcg)

#        bt = zeros_like(bp)
        vp = zeros_like(bp)
        vt = zeros_like(bp)

        for i in range(prm.NO2):
            # note, by indexing below, we're making sure we are dumping arrays with shape (nk,nj)
            # instead of thinking whether they are already of this size
            out = array([ bt_jface_a[:nk,:nj].T.ravel(),             
                          bp_kface[:nk,:nj,i].T.ravel(),        
                          et[:nk,:nj].T.ravel(),
                          ep_save[:nj,:nk].ravel(),
                          br[:nk,:nj,i].T.ravel(),
                          bt_a[:nk,:nj].T.ravel(),
                          bp[:nk,:nj,i].T.ravel(),
                          vr[:nk,:nj,i].T.ravel(),
                          vt[:nk,:nj,i].T.ravel(),
                          vp[:nk,:nj,i].T.ravel(),
                          rho[:nk,:nj,i].T.ravel(),
                          cs[:nk,:nj,i].T.ravel()])
                  
            savetxt(os.path.join(prm.dirInitLFMfile,'innerbc_%03d_%d.dat'% (fcount+1,i)),out.T,
                    # fmt=['%13.8f','%13.8f','%15.8f','%13.8f','%13.8f','%13.8f',
                    #      '%15.5e','%15.5e','%15.5e',
                    #      '%14.5e','%14.5e'],
                    fmt=['%15.5e']*12,
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




