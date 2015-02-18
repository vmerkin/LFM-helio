from scipy.optimize import newton_krylov,anderson
from numpy import ediff1d,zeros_like,cos,sin,zeros_like,arange

class poisson():
    def __init__(self,theta,phi):
        self.theta  = theta
        self.phi    = phi
        self.nj     = theta.shape[0]
        self.nk     = phi.shape[0]
        self.dtheta = theta[1]-theta[0]   # generalize to the case of arbitrary grid spacing later: ediff1d(theta)
        self.dphi   = phi[1]-phi[0]       # ediff1d(phi)

    def setRHS(self,RHS):
        self.RHS = RHS

    def residual_xsin2(self,P,lhs=False):
        dphi = self.dphi
        dtheta = self.dtheta
        theta = self.theta
        nk = self.nk


        d2x = zeros_like(P)
        d2y = zeros_like(P)
        dy = zeros_like(P)

        d2x[:,1:-1] = (P[:,2:] - 2*P[:,1:-1] + P[:,:-2]) /dphi**2
        d2y[1:-1,:] = (P[2:,:] - 2*P[1:-1,:] + P[:-2,:]) /dtheta**2*sin(theta[1:-1,None])**2

        dy[1:-1,:] = 0.5*(P[2:,:] - P[:-2,:]) /dtheta*sin(theta[1:-1,None])*cos(theta[1:-1,None]) 

        # periodic boundary
        d2x[:,0]  = (P[:,1]-2*P[:,0]+P[:,-1]) /dphi**2
        d2x[:,-1] = (P[:,0]-2*P[:,-1]+P[:,-2]) /dphi**2

        result = self.RHS*sin(theta[:,None])**2

        # Pole boundary condition. Note, the below does not mean that d/dtheta=0 literally.
        # We set sin(theta)*dPsi/dtheta=0 at the pole (see notes for details), 
        # and d2y here represents 1/sin(theta)*d/dtheta (sin(theta)*dPsi/dtheta) in that approxiamtion.
        # NOTE, BELOW ASSUMES CONSTANT SPACING IN THETA!!! (0.5 FACTOR IS BASED ON THAT).
        dy[0,:]=0.
        dy[-1,:]=0.

        d2y[0,:] = 0.5*(P[1,:]-P[0,:])
        d2y[-1,:] = -0.5*(P[-1,:]-P[-2,:])

        ###################################################################################################
        if not lhs:
            return d2x + d2y + dy - result
        else:
            return d2x + d2y + dy

    def residual(self,P,lhs=False):
        dphi = self.dphi
        dtheta = self.dtheta
        theta = self.theta
        nk = self.nk


        d2x = zeros_like(P)
        d2y = zeros_like(P)
        dy = zeros_like(P)

        d2x[:,1:-1] = (P[:,2:] - 2*P[:,1:-1] + P[:,:-2]) /dphi**2/sin(theta[:,None])**2 
        d2y[1:-1,:] = (P[2:,:] - 2*P[1:-1,:] + P[:-2,:]) /dtheta**2

        dy[1:-1,:] = 0.5*(P[2:,:] - P[:-2,:]) /dtheta/sin(theta[1:-1,None])*cos(theta[1:-1,None])

        # periodic boundary
        d2x[:,0]  = (P[:,1]-2*P[:,0]+P[:,-1]) /dphi**2/sin(theta)**2 
        d2x[:,-1] = (P[:,0]-2*P[:,-1]+P[:,-2]) /dphi**2/sin(theta)**2

        result = self.RHS

        # Pole boundary condition. Note, the below does not mean that d/dtheta=0 literally.
        # We set sin(theta)*dPsi/dtheta=0 at the pole (see notes for details), 
        # and d2y here represents 1/sin(theta)*d/dtheta (sin(theta)*dPsi/dtheta) in that approxiamtion.
        # NOTE, BELOW ASSUMES CONSTANT SPACING IN THETA!!! (0.5 FACTOR IS BASED ON THAT).
        dy[0,:]=0.
        dy[-1,:]=0.

        d2y[0,:] = 0.5*(P[1,:]-P[0,:])/sin(theta[0])**2   # Psin=0; 0.5*( 0.5*(P[1,:]+P[0,:]) - Psin)
        d2y[-1,:] = -0.5*(P[-1,:]-P[-2,:])/sin(theta[-1])**2     #-0.5*( 0.5*(P[-2,:]+P[-1,:]) - P[-1,:].mean())

        ###################################################################################################
        if not lhs:
            return d2x + d2y + dy - result
        else:
            return d2x + d2y + dy

if __name__=='__main__':
    from scipy.optimize import newton_krylov,anderson
    from pylab import *

    import wsa
    (phi_wsa_v,theta_wsa_v,phi,theta,br0,v_wsa,n_wsa,T_wsa)=wsa.read('/Users/merkivg1/work/LFM-helio_2.0/WSA/fits/2010_2008_WSA_ADAPT/vel_201001022300R001_ans.fits',False)

    (phi_wsa_v,theta_wsa_v,phi,theta,br1,v_wsa,n_wsa,T_wsa)=wsa.read('/Users/merkivg1/work/LFM-helio_2.0/WSA/fits/2010_2008_WSA_ADAPT/vel_201001032300R001_ans.fits',False)

    omega=2*pi/27.27
    phi_prime=phi-omega*1.
    ind0=where(phi_prime>0)[0][0]
    phi_prime=roll(phi_prime,-ind0)
    phi_prime[phi_prime<0]+=2*pi

    import scipy
    from scipy.interpolate import interpolate
    br1_rolled=roll(br1,-ind0,axis=1)
    fbi      = scipy.interpolate.RectBivariateSpline(theta,phi,br1_rolled,kx=1,ky=1)
    br1_new=fbi(theta,phi_prime)

    import poisson
    pois = poisson.poisson(theta,phi)
    pois.setRHS(br1_new - br0)

    guess=zeros_like(br1_new)
    sol = newton_krylov(pois.residual,guess, method='lgmres', verbose=1,f_rtol=1.e-6) #iter=100
    print('Residual: %g' % abs(pois.residual(sol)).max())

    figure();pcolormesh(sol);colorbar();title('Solution')

    # check curl
    Etheta = -1/sin(theta[:,None])*diff(sol,axis=1)
    Ephi = diff(sol,axis=0)
    figure();pcolormesh(Etheta);colorbar();title('Etheta')
    figure();pcolormesh(Ephi);colorbar();title('Ephi')

