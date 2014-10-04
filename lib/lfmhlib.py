from pyhdf.SD import SD, SDC
from numpy import sqrt,arccos,arctan2,cos,sin,pi

def read(filename):


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

    return(R,theta,phi,Rc,thetac,phic,br,btheta,bphi,vr,vtheta,vphi,rho,cs)
