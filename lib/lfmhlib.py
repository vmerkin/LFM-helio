from pyhdf.SD import SD, SDC
from numpy import sqrt,arccos,arctan2,cos,sin,pi,roll

def read(filename):
    """
    Returns R,theta,phi,Rc,thetac,phic,br,btheta,bphi,vr,vtheta,vphi,rho,cs
    """

    hdffile = SD(filename,SDC.READ)

    x=hdffile.select('X_grid').get()
    y=hdffile.select('Y_grid').get()
    z=hdffile.select('Z_grid').get()
    bx=hdffile.select('bx_').get()[:-1,:-1,:-1]
    by=hdffile.select('by_').get()[:-1,:-1,:-1]
    bz=hdffile.select('bz_').get()[:-1,:-1,:-1]
    vx=hdffile.select('vx_').get()[:-1,:-1,:-1]
    vy=hdffile.select('vy_').get()[:-1,:-1,:-1]
    vz=hdffile.select('vz_').get()[:-1,:-1,:-1]
    rho=hdffile.select('rho_').get()[:-1,:-1,:-1]
    cs=hdffile.select('c_').get()[:-1,:-1,:-1]

    t=hdffile.time
    hdffile.end()

    # =========== Cell centers ==============
    xc=0.125*( x[:-1,:-1,:-1]+x[1:,:-1,:-1]+x[:-1,1:,:-1]+x[:-1,:-1,1:]+
              x[1:,1:,:-1]+x[1:,:-1,1:]+x[:-1,1:,1:]+x[1:,1:,1:] )
    yc=0.125*( y[:-1,:-1,:-1]+y[1:,:-1,:-1]+y[:-1,1:,:-1]+y[:-1,:-1,1:]+
              y[1:,1:,:-1]+y[1:,:-1,1:]+y[:-1,1:,1:]+y[1:,1:,1:] )
    zc=0.125*( z[:-1,:-1,:-1]+z[1:,:-1,:-1]+z[:-1,1:,:-1]+z[:-1,:-1,1:]+
              z[1:,1:,:-1]+z[1:,:-1,1:]+z[:-1,1:,1:]+z[1:,1:,1:] )
    # =======================================
    
    R=sqrt(x**2+y**2+z**2)
    theta=arccos(z/R)
    phi=arctan2(y,x)
    phi[phi<0]+=2*pi

    Rc=sqrt(xc**2+yc**2+zc**2)
    thetac=arccos(zc/Rc)
    phic=arctan2(yc,xc)
    phic[phic<0]+=2*pi


    br     = bx*cos(phic)*sin(thetac) + by*sin(phic)*sin(thetac) + bz*cos(thetac)
    btheta = bx*cos(phic)*cos(thetac) + by*sin(phic)*cos(thetac) - bz*sin(thetac)
    bphi   =-bx*sin(phic)            + by*cos(phic)

    vr     = vx*cos(phic)*sin(thetac) + vy*sin(phic)*sin(thetac) + vz*cos(thetac)
    vtheta = vx*cos(phic)*cos(thetac) + vy*sin(phic)*cos(thetac) - vz*sin(thetac)
    vphi   =-vx*sin(phic)            + vy*cos(phic)

    return(t,R,theta,phi,Rc,thetac,phic,br,btheta,bphi,vr,vtheta,vphi,rho,cs)

def r_theta_phi_uniform(filename):
    """
    Return R, theta, phi 1-d arrays, assuming uniform grid
    """

    hdffile = SD(filename,SDC.READ)

    # note, minimal data are read from file

    ni=hdffile.select('X_grid').ni-1
    nj=hdffile.select('X_grid').nj-1
    nk=hdffile.select('X_grid').nk-1


    # first get R
    x=hdffile.select('X_grid').get()[0,0,:] 
    y=hdffile.select('Y_grid').get()[0,0,:] 
    z=hdffile.select('Z_grid').get()[0,0,:] 
    R = sqrt(x**2+y**2+z**2)

    # first get R
    z=hdffile.select('Z_grid').get()[0,:,0] 
    theta = arccos(z/R[0])

    x=hdffile.select('X_grid').get()[:,nj/2,0]
    y=hdffile.select('Y_grid').get()[:,nj/2,0]
    phi=arctan2(y,x)
    phi[phi<0]+=2*pi
    phi[-1]=phi[0]+2*pi
    
    hdffile.end()
    
    return R,theta,phi

def read_var(filename,var_name):
    hdffile = SD(filename,SDC.READ)
    if var_name not in ['br','btheta','bphi','vr','vtheta','vphi']:
        var=hdffile.select(var_name+'_').get()[:-1,:-1,:-1]
    else:
        R,theta,phi=r_theta_phi_uniform(filename)
        thetac = 0.5*(theta[1:]+theta[:-1])
        phic = 0.5*(phi[1:]+phi[:-1])
        
        phic = phic[:,None,None]
        thetac = thetac[None,:,None]

        if var_name in ['br','btheta','bphi']:
            bx=hdffile.select('bx_').get()[:-1,:-1,:-1]
            by=hdffile.select('by_').get()[:-1,:-1,:-1]
            bz=hdffile.select('bz_').get()[:-1,:-1,:-1]

            if var_name=='br':
                var     = bx*cos(phic)*sin(thetac) + by*sin(phic)*sin(thetac) + bz*cos(thetac)
            elif var_name=='btheta':
                var = bx*cos(phic)*cos(thetac) + by*sin(phic)*cos(thetac) - bz*sin(thetac)
            else:
                var   =-bx*sin(phic)            + by*cos(phic)
        else:
            vx=hdffile.select('vx_').get()[:-1,:-1,:-1]
            vy=hdffile.select('vy_').get()[:-1,:-1,:-1]
            vz=hdffile.select('vz_').get()[:-1,:-1,:-1]

            if var_name=='vr':
                var    = vx*cos(phic)*sin(thetac) + vy*sin(phic)*sin(thetac) + vz*cos(thetac)
            elif var_name=='vtheta':
                var = vx*cos(phic)*cos(thetac) + vy*sin(phic)*cos(thetac) - vz*sin(thetac)
            else:
                var   =-vx*sin(phic)            + vy*cos(phic)
    hdffile.end()
    return(var)

def get_time(filename):
    hdffile = SD(filename, SDC.READ)
    t=hdffile.time
    hdffile.end()
    return(t)
    
def time_shift(t,var,Tsolar=25.38*24*3600.):
    """
    Shift variables to the beginning of the CR.

    Parameters:
    t: current time
    var: variable -- assumed to be a 2D i-slice (nk,nj) of the LFM data 
    Tsolar: solar rotation period in units of time used in the code. Default to 25.38 days and seconds as time units.
    """
    time_shift= (t % Tsolar)/Tsolar
    nshift = int(time_shift*(var.shape[0]+1)) # note, using nk+1 here, which gives us full rotation
    return roll(var,-nshift,axis=0)

