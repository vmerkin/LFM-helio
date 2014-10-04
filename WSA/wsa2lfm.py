def sph2cart(a,theta,phi):
    ar,at,ap = a

    ax = ar*cos(phi)*sin(theta)+at*cos(phi)*cos(theta)-ap*sin(phi)
    ay = ar*sin(phi)*sin(theta)+at*sin(phi)*cos(theta)+ap*cos(phi)
    az = ar*cos(theta)-at*sin(theta)

    return(ax,ay,az)




import ConfigParser
config = ConfigParser.ConfigParser()
config.read('wsa.config')

ni = config.getint('Dimensions','NI')
nj = config.getint('Dimensions','NJ')
nk = config.getint('Dimensions','NK')
fileName = config.get('OutputFileName','Prefix')+'_%dx%dx%d_mhd_0000000.hdf'%(ni,nj,nk)

scale = config.getfloat('Normalization','Xscale')
rmin = config.getfloat('GridSpecs','RMIN')
rmax = config.getfloat('GridSpecs','RMAX')
thetamin = config.getfloat('GridSpecs','THETAMIN')


gamma = config.getfloat('Constants','gamma')
NO2   = config.getint('Constants','NO2')
Tsolar = config.getfloat('Constants','Tsolar')

wsa_file = config.get('WSA','wsafile')
gauss_smooth_width = config.getint('WSA','gauss_smooth_width')

dumpInit = True
dumpBC   = False
vis      = False





from numpy import linspace,pi,meshgrid,sin,cos,zeros,ones,dstack,diff,sqrt,array,savetxt,flatnonzero,insert,asarray,zeros_like
import os,sys,glob
import lfmpy5 as lfmpy
from scipy import interpolate
import time

############### WSA STUFF #####################
import pyfits
hdulist = pyfits.open(wsa_file)
print hdulist[0].header.ascardlist()

n_phi_wsa_v = hdulist[0].header['NAXIS1']+1  # number of cell vertices
n_phi_wsa_c = hdulist[0].header['NAXIS1']    # number of cell centers
phi_wsa_v     = linspace(0,360,n_phi_wsa_v)/180.*pi
phi_wsa_c     = 0.5*(phi_wsa_v[:-1]+phi_wsa_v[1:])

n_theta_wsa_v = hdulist[0].header['NAXIS2']+1  # number of cell vertices
n_theta_wsa_c = hdulist[0].header['NAXIS2']    # number of cell centers
theta_wsa_v     = linspace(0,180,n_theta_wsa_v)/180.*pi
theta_wsa_c     = 0.5*(theta_wsa_v[:-1]+theta_wsa_v[1:])

bi_wsa = hdulist[0].data[0,::-1,:]  # note the theta reversal to convert from wsa theta to lfm theta definition
v_wsa  = hdulist[0].data[1,::-1,:]  # note the theta reversal to convert from wsa theta to lfm theta definition
n_wsa  = hdulist[0].data[2,::-1,:]
T_wsa  = hdulist[0].data[3,::-1,:]

hdulist.close()
############### WSA STUFF #####################

############### LFM STUFF #####################
# LFM GRID
grid = lfmpy.grids((ni,nj,nk))
(P,T,R,Pc,Tc,Rc,x,y,z,xc,yc,zc) = grid.getGrid(rmin=rmin*scale,rmax=rmax*scale,thetamin=thetamin)

# this is fast and better than griddata in that it nicely extrapolates boundaries:
fbi      = interpolate.RectBivariateSpline(phi_wsa_c,theta_wsa_c,bi_wsa.T,kx=1,ky=1)  
br = fbi(Pc[:,0,0],Tc[0,:,0])

############### SMOOTHING #####################
if not gauss_smooth_width==0:
    import astropy
    from astropy.convolution import convolve,Gaussian2DKernel

    gauss=Gaussian2DKernel(width=gauss_smooth_width)
    br   =astropy.convolution.convolve(br,gauss,boundary='extend')

fv      = interpolate.RectBivariateSpline(phi_wsa_c,theta_wsa_c,v_wsa.T,kx=1,ky=1)  
vr = fv(Pc[:,0,0],Tc[0,:,0])

f      = interpolate.RectBivariateSpline(phi_wsa_c,theta_wsa_c,n_wsa.T,kx=1,ky=1)  
rho = f(Pc[:,0,0],Tc[0,:,0])*1.67e-24

f      = interpolate.RectBivariateSpline(phi_wsa_c,theta_wsa_c,T_wsa.T,kx=1,ky=1)  
temp = f(Pc[:,0,0],Tc[0,:,0])


cs = sqrt(gamma*1.38e-23*temp/1.67e-27)*1.e2   # in cm/s

############### LFM STUFF #####################
if dumpInit:
    with lfmpy.startup(fileName,(ni,nj,nk)) as lfmh:
        # Set variables
        lfmh.setVar('X_grid',x)
        lfmh.setVar('Y_grid',y)
        lfmh.setVar('Z_grid',z)
        
        # lfmh.setVar('rho_',400.*1.67e-24*ones((nk+1,nj+1,ni+1)))
        # lfmh.setVar('c_',  5.e6*ones((nk+1,nj+1,ni+1)))
        # lfmh.setVar('vx_',4.e7*sin(Tc)*cos(Pc))
        # lfmh.setVar('vy_',4.e7*sin(Tc)*sin(Pc))
        # lfmh.setVar('vz_',4.e7*cos(Tc))

        lfmh.setVar('rho_',rho[:,:,None]*(Rc[0,0,0]/Rc)**2)
        lfmh.setVar('c_',  cs[:,:,None]*ones((nk,nj,ni)))
        lfmh.setVar('vx_',vr[:,:,None]*sin(Tc)*cos(Pc))
        lfmh.setVar('vy_',vr[:,:,None]*sin(Tc)*sin(Pc))
        lfmh.setVar('vz_',vr[:,:,None]*cos(Tc))


        tmp    = diff(T[:,:,0],axis=1)
        dtheta = 0.5*(tmp[:-1,:]+tmp[1:,:])
        tmp    = diff(P[:,:,0],axis=0)
        dphi   = 0.5*(tmp[:,:-1]+tmp[:,1:])

        # note the fancy dstack to fill in the last i-face with correct bi values
        #
        # note we use Rc[0,0,0] instead of rmin*scale (half-cell shift)
        # we thus assume that WSA data are not at rmin but at the center of the first cell 
        # and scale bi appropriately
        bi = dstack(R.shape[2]*[br*Rc[0,0,0]**2*sin(Tc[:,:,0])*dtheta*dphi])
        lfmh.setVar('bi_',bi)

if dumpBC:
    # ghost cell coordinates
    (xg,yg,zg,xcg,ycg,zcg) = [2*arr[:,:,[0]]-arr[:,:,1:NO2+1] for arr in (x,y,z,xc,yc,zc)]
    # Here we assume grid radii don't vary with theta and phi, at least in the ghost region
    rcg = sqrt(xcg[0,0,:]**2+ycg[0,0,:]**2+zcg[0,0,:]**2)/scale
    rg = sqrt(xg[0,0,:]**2+yg[0,0,:]**2+zg[0,0,:]**2)/scale

    # note, redefining interpolation functions we could also
    # interpolate from bi_wsa as above, but then we would have to
    # smooth bk, if necessary. The way we're doing it here, bk will be
    # smoothed or not, dependent on whether br has been smoothed.
    # note also, this has to extrapolate
    fbi = interpolate.RectBivariateSpline(Pc[:,0,0],Tc[0,:,0],br,kx=1,ky=1)
    fv  = interpolate.RectBivariateSpline(Pc[:,0,0],Tc[0,:,0],vr,kx=1,ky=1)
    br_kface = fbi(P[:,0,0],Tc[0,:,0])
    vr_kface  = fv (P[:,0,0],Tc[0,:,0])

    # note we use Rc[0,0,0] instead of rmin*scale (half-cell shift)
    # we thus assume that WSA data are not at rmin but at the center of the first cell 
    # and scale bi appropriately
    bp_kface = -br_kface*2*pi/(Tsolar*24.*3600.)*Rc[0,0,0]*sin(Tc[0,:,0])/vr_kface  # beautiful numpy broadcasting
    bp       = -br      *2*pi/(Tsolar*24.*3600.)*Rc[0,0,0]*sin(Tc[0,:,0])/vr          




    (vr,rho,cs,br,bp,bp_kface) = [dstack(NO2*[var]) for var in (vr,rho,cs,br,bp,bp_kface)]
    rho*=(Rc[0,0,0]/scale/rcg)**2
    br*=(Rc[0,0,0]/scale/rcg)**2
    bp*=(Rc[0,0,0]/scale/rcg)
    bp_kface*=(Rc[0,0,0]/scale/rcg)

    bt = zeros_like(bp)
    vp = zeros_like(bp)
    vt = zeros_like(bp)

    bx,by,bz = sph2cart( (br,bt,bp),Tc[:,:,[0]],Pc[:,:,[0]] )
    vx,vy,vz = sph2cart( (vr,vt,vp),Tc[:,:,[0]],Pc[:,:,[0]] )


    for i in range(NO2):
        # note, by indexing below, we're making sure we are dumping arrays with shape (nk,nj)
        # instead of thinking whether they are already of this size
        out = array([ bt[:nk,:nj,i].T.ravel(),
                      bp_kface[:nk,:nj,i].T.ravel(),
                      # bx[:nk,:nj,i].T.ravel(),
                      # by[:nk,:nj,i].T.ravel(),
                      # bz[:nk,:nj,i].T.ravel(),
                      # vx[:nk,:nj,i].T.ravel(),
                      # vy[:nk,:nj,i].T.ravel(),
                      # vz[:nk,:nj,i].T.ravel(),
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


if vis:
    figure(); pcolormesh(phi_wsa_v,theta_wsa_v,bi_wsa); colorbar()
    figure(); pcolormesh(phi_lfm_v,theta_lfm_v,bi_lfm); colorbar()

    if smooth: 
        figure(); pcolormesh(phi_lfm_v,theta_lfm_v,bi_lfm_smooth); colorbar()
        contour(phi_lfm_c,theta_lfm_c,bi_lfm_smooth,[0.])
        contour(phi_lfm_v[:-1],theta_lfm_c,bk_lfm_smooth,[0.],colors='purple',linewidth=5)
    

    show()


