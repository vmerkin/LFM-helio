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

    def residual(self,P):
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

        # poles, set to zero for now
        k    = arange(nk)
        kppi = (nk/2+k)%nk   # k on the opposite side (phi+pi, thus kppi, i.e. "k plus pi")

        ####################################################################################################
        # FLIP AROUND AXIS -- INCORRECT TO USE IN SPHERICAL COORDINATES (NEED
        # TO DO CALCULATIONS IN CARTESIAN FOR THIS TO BE VALID -- LIKE IN LFM
        #
        #    d2y[0,:] = (P[1,:]-2*P[0,:]+P[0,kppi]) /dtheta**2*sin(theta[0])**2
        #    d2y[-1,:] = (P[-1,kppi]-2*P[-1,:]+P[-2,:]) /dtheta**2*sin(theta[-1])**2
        #
        #    dy[0,:] = 0.5*(P[1,:] - P[0,kppi]) /dtheta*sin(theta[0])*cos(theta[0])
        #    dy[-1,:] = 0.5*(P[-1,kppi] - P[-2,:]) /dtheta*sin(theta[-1])*cos(theta[-1])
        ####################################################################################################

        ####################################################################################################
        # FIRST ORDER
        d2y[0,:] = 0.
        d2y[-1,:] = 0.
  
        dy[0,:] = P[0,:]*cos(theta[0])
        dy[-1,:] = P[-1,:]*cos(theta[-1])
        ###################################################################################################

        ####################################################################################################
        # SECOND ORDER
        # d2y[0,:] = P[0,:]+P[0,kppi]
        # d2y[-1,:] = P[-1,:]+P[-1,kppi]
        
        # dy[0,:] = 0.5*(3*P[0,:]+P[0,kppi])*cos(theta[0])
        # dy[-1,:] = 0.5*(3*P[-1,:]+P[-1,kppi])*cos(theta[0])
        ####################################################################################################
        return d2x + d2y + dy - self.RHS*sin(theta[:,None])**2



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

