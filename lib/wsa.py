import pyfits
from numpy import linspace,pi,meshgrid,sin,cos,zeros,ones,dstack,diff,sqrt,array,savetxt,flatnonzero,insert,asarray,zeros_like,argmin,unravel_index

def read(wsa_file,densTempInfile,normalized=False):
    hdulist = pyfits.open(wsa_file)
    print hdulist[0].header.cards
    
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
    if densTempInfile:
        n_wsa  = hdulist[0].data[2,::-1,:]
        T_wsa  = hdulist[0].data[3,::-1,:]
    else:
        n_wsa = 112.64+9.49e7/v_wsa**2
        T0 = 8e5
        n0 = 300.
        
        B0 = bi_wsa[unravel_index(argmin(abs(n_wsa-n0)),n_wsa.shape)] # this is in nT
        n0 = n_wsa[unravel_index(argmin(abs(n_wsa-n0)),n_wsa.shape)]  # this is in cm^-3

        T_wsa = n0*T0/n_wsa
#        T_wsa = n0*T0/n_wsa + (B0**2-bi_wsa**2)/1.38/8./pi*1.e6/n_wsa
    hdulist.close()

    if normalized:
        return(phi_wsa_v,theta_wsa_v,phi_wsa_c,theta_wsa_c,bi_wsa,v_wsa,n_wsa*1.67e-24,T_wsa)
    else:
        return(phi_wsa_v,theta_wsa_v,phi_wsa_c,theta_wsa_c,bi_wsa*1.e-5,v_wsa*1.e5,n_wsa*1.67e-24,T_wsa)

def info(wsa_file):
    hdulist = pyfits.open(wsa_file)
    print hdulist[0].header.cards
    hdulist.close()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('wsaFile',help='WSA file to use')
    args = parser.parse_args()
    
    info(args.wsaFile)
    
