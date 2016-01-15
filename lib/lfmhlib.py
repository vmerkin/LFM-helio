from pyhdf.SD import SD, SDC
from numpy import sqrt,arccos,arctan2,cos,sin,pi,roll
import time

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


    # first get R; in principle, could just take z along the axis
    # but are allowing for the possibility that the axis is cut out
    x=hdffile.select('X_grid').get(start=(0,0,0),count=(1,1,ni+1)).squeeze()
    y=hdffile.select('Y_grid').get(start=(0,0,0),count=(1,1,ni+1)).squeeze()
    z=hdffile.select('Z_grid').get(start=(0,0,0),count=(1,1,ni+1)).squeeze()
    R = sqrt(x**2+y**2+z**2)

    z=hdffile.select('Z_grid').get(start=(0,0,0),count=(1,nj+1,1)).squeeze()
    theta = arccos(z/R[0])

    x=hdffile.select('X_grid').get(start=(0,nj/2,0),count=(nk+1,1,1)).squeeze()
    y=hdffile.select('Y_grid').get(start=(0,nj/2,0),count=(nk+1,1,1)).squeeze()
    phi=arctan2(y,x)
    phi[phi<0]+=2*pi
    phi[-1]=phi[0]+2*pi
    
    hdffile.end()

    R/=6.96e10   # FIX ME HARD CODED UNITS
    Rc = 0.5*(R[1:]+R[:-1])
    thetac = 0.5*(theta[1:]+theta[:-1])
    phic = 0.5*(phi[1:]+phi[:-1])
    
    return R,theta,phi,Rc,thetac,phic

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

def read_var_point(filename,var_name,i,j,k,thetac,phic):
    thetac = thetac[j]
    phic   = phic[k]


    hdffile = SD(filename,SDC.READ)
    if var_name not in ['br','btheta','bphi','vr','vtheta','vphi']:
        var=hdffile.select(var_name+'_').get(start=(k,j,i),count=(1,1,1)).squeeze()
    else:
#        R,theta,phi=r_theta_phi_uniform(filename)

        if var_name in ['br','btheta','bphi']:
            bx=hdffile.select('bx_').get(start=(k,j,i),count=(1,1,1)).squeeze()
            by=hdffile.select('by_').get(start=(k,j,i),count=(1,1,1)).squeeze()
            bz=hdffile.select('bz_').get(start=(k,j,i),count=(1,1,1)).squeeze()

            if var_name=='br':
                var     = bx*cos(phic)*sin(thetac) + by*sin(phic)*sin(thetac) + bz*cos(thetac)
            elif var_name=='btheta':
                var = bx*cos(phic)*cos(thetac) + by*sin(phic)*cos(thetac) - bz*sin(thetac)
            else:
                var   =-bx*sin(phic)            + by*cos(phic)
        else:
            vx=hdffile.select('vx_').get(start=(k,j,i),count=(1,1,1)).squeeze()
            vy=hdffile.select('vy_').get(start=(k,j,i),count=(1,1,1)).squeeze()
            vz=hdffile.select('vz_').get(start=(k,j,i),count=(1,1,1)).squeeze()

            if var_name=='vr':
                var    = vx*cos(phic)*sin(thetac) + vy*sin(phic)*sin(thetac) + vz*cos(thetac)
            elif var_name=='vtheta':
                var = vx*cos(phic)*cos(thetac) + vy*sin(phic)*cos(thetac) - vz*sin(thetac)
            else:
                var   =-vx*sin(phic)            + vy*cos(phic)
    hdffile.end()
    return(var)

def read_var_islice(filename,var_name,i,thetac,phic):
    nk = phic.size
    nj = thetac.size
    phic = phic[:,None]
    thetac = thetac[None,:]

    hdffile = SD(filename,SDC.READ)
    if var_name not in ['br','btheta','bphi','vr','vtheta','vphi']:
        var=hdffile.select(var_name+'_').get(start=(0,0,i),count=(nk,nj,1)).squeeze()
    else:
        if var_name in ['br','btheta','bphi']:
            bx=hdffile.select('bx_').get(start=(0,0,i),count=(nk,nj,1)).squeeze()
            by=hdffile.select('by_').get(start=(0,0,i),count=(nk,nj,1)).squeeze()
            bz=hdffile.select('bz_').get(start=(0,0,i),count=(nk,nj,1)).squeeze()

            if var_name=='br':
                var     = bx*cos(phic)*sin(thetac) + by*sin(phic)*sin(thetac) + bz*cos(thetac)
            elif var_name=='btheta':
                var = bx*cos(phic)*cos(thetac) + by*sin(phic)*cos(thetac) - bz*sin(thetac)
            else:
                var   =-bx*sin(phic)            + by*cos(phic)
        else:
            vx=hdffile.select('vx_').get(start=(0,0,i),count=(nk,nj,1)).squeeze()
            vy=hdffile.select('vy_').get(start=(0,0,i),count=(nk,nj,1)).squeeze()
            vz=hdffile.select('vz_').get(start=(0,0,i),count=(nk,nj,1)).squeeze()

            if var_name=='vr':
                var    = vx*cos(phic)*sin(thetac) + vy*sin(phic)*sin(thetac) + vz*cos(thetac)
            elif var_name=='vtheta':
                var = vx*cos(phic)*cos(thetac) + vy*sin(phic)*cos(thetac) - vz*sin(thetac)
            else:
                var   =-vx*sin(phic)            + vy*cos(phic)
    hdffile.end()
    return(var)

def read_var_jslice(filename,var_name,j,thetac,phic):
    hdffile = SD(filename,SDC.READ)
    ni = hdffile.select('X_grid').ni-1
    nk = hdffile.select('X_grid').nk-1
    phic = phic[:,None]
    thetac = thetac[j]

    if var_name not in ['br','btheta','bphi','vr','vtheta','vphi']:
        if var_name in ['X_grid','Y_grid','Z_grid']:
            var=hdffile.select(var_name).get(start=(0,j,0),count=(nk,1,ni)).squeeze()
        else:
            var=hdffile.select(var_name+'_').get(start=(0,j,0),count=(nk,1,ni)).squeeze()
    else:
        if var_name in ['br','btheta','bphi']:
            bx=hdffile.select('bx_').get(start=(0,j,0),count=(nk,1,ni)).squeeze()
            by=hdffile.select('by_').get(start=(0,j,0),count=(nk,1,ni)).squeeze()
            bz=hdffile.select('bz_').get(start=(0,j,0),count=(nk,1,ni)).squeeze()

            if var_name=='br':
                var     = bx*cos(phic)*sin(thetac) + by*sin(phic)*sin(thetac) + bz*cos(thetac)
            elif var_name=='btheta':
                var = bx*cos(phic)*cos(thetac) + by*sin(phic)*cos(thetac) - bz*sin(thetac)
            else:
                var   =-bx*sin(phic)            + by*cos(phic)
        else:
            vx=hdffile.select('vx_').get(start=(0,j,0),count=(nk,1,ni)).squeeze()
            vy=hdffile.select('vy_').get(start=(0,j,0),count=(nk,1,ni)).squeeze()
            vz=hdffile.select('vz_').get(start=(0,j,0),count=(nk,1,ni)).squeeze()

            if var_name=='vr':
                var    = vx*cos(phic)*sin(thetac) + vy*sin(phic)*sin(thetac) + vz*cos(thetac)
            elif var_name=='vtheta':
                var = vx*cos(phic)*cos(thetac) + vy*sin(phic)*cos(thetac) - vz*sin(thetac)
            else:
                var   =-vx*sin(phic)            + vy*cos(phic)
    hdffile.end()
    return(var)

def read_var_ikslice(filename,var_name,i,k,tc,pc):
    hdffile = SD(filename,SDC.READ)
    ni = hdffile.select('X_grid').ni-1
    nj = hdffile.select('X_grid').nj-1

    if var_name not in ['br','btheta','bphi','vr','vtheta','vphi','bt','bp','vt','vp']:
        var=hdffile.select(var_name+'_').get(start=(k,0,i),count=(1,nj,1)).squeeze()
    else:
        if var_name in ['br','btheta','bphi','bt','bp']:
            bx=hdffile.select('bx_').get(start=(k,0,i),count=(1,nj,1)).squeeze()
            by=hdffile.select('by_').get(start=(k,0,i),count=(1,nj,1)).squeeze()
            bz=hdffile.select('bz_').get(start=(k,0,i),count=(1,nj,1)).squeeze()

            if var_name=='br':
                var     = bx*cos(pc[k])*sin(tc) + by*sin(pc[k])*sin(tc) + bz*cos(tc)
            elif (var_name=='btheta' or var_name=='bt'):
                var = bx*cos(pc[k])*cos(tc) + by*sin(pc[k])*cos(tc) - bz*sin(tc)
            else:
                var   =-bx*sin(pc[k])            + by*cos(pc[k])
        else:
            vx=hdffile.select('vx_').get(start=(k,0,i),count=(1,nj,1)).squeeze()
            vy=hdffile.select('vy_').get(start=(k,0,i),count=(1,nj,1)).squeeze()
            vz=hdffile.select('vz_').get(start=(k,0,i),count=(1,nj,1)).squeeze()

            if var_name=='vr':
                var    = vx*cos(pc[k])*sin(tc) + vy*sin(pc[k])*sin(tc) + vz*cos(tc)
            elif (var_name=='vtheta' or var_name=='vt'):
                var = vx*cos(pc[k])*cos(tc) + vy*sin(pc[k])*cos(tc) - vz*sin(tc)
            else:
                var   =-vx*sin(pc[k])            + vy*cos(pc[k])
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

