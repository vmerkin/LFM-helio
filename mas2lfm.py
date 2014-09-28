"""
This script is adapted from ~/work/LFM-helio_2.0/mascme2lfm_helio_coronal_interp.py.
Here I'm going to try to keep things clean and modular.
"""

import os,sys,glob
import mas
from numpy import linspace,pi,meshgrid,sin,cos,zeros,ones,dstack,diff,sqrt,array,savetxt,flatnonzero,insert,asarray,zeros_like
from scipy import interpolate
import time
import pyLTR
import mas2lfm


#----------- PARSE ARGUMENTS ---------#
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('ConfigFileName',help='The name of the configuration file to use',default='startup.config')
args = parser.parse_args()
#----------- PARSE ARGUMENTS ---------#

# Read params from config file
prm = mas2lfm.params.params(args.ConfigFileName)
(ni,nj,nk) = (prm.ni,prm.nj,prm.nk)

if prm.masTimeLabel=='all':
    timeLabels = [int(os.path.basename(s).split('.')[0][-3:]) for s in sorted(glob.glob(os.path.join(prm.masdir,'rho*.hdf')))]
else:
    masTimeLabel = int(prm.masTimeLabel)
    timeLabels = [masTimeLabel]


vars={}
vars['mas']={}
interp={}

for timeLabel in timeLabels:
    print('MAS time label: %d'%timeLabel)
    ######## READ IN MAS DATA #####################
    vars['mas'] = mas.read_all_vars(prm.masdir,timeLabel)

    # We only assume here that the LFM bottom boundary coincides with the shell where MAS br is defined. 
    # We don't have to do this: could define rmin arbitrarily and then linearly interpolate.
    #
    # Do it only for the first MAS file. If we did not have to find rmin here,
    # we could move the LFM grid and ghost defs above the timeLabel loop

    if timeLabel == timeLabels[0]:

        if prm.masFrame=='rotating':
            # transform the azimuthal velocity to inertial frame
            R,T=meshgrid(vars['mas']['vp']['r'],vars['mas']['vp']['theta'])
            # using numpy broadcasting here
            Vphi_solar = 2*pi/(prm.Tsolar*24.*3600./mas.time_unit())*R*sin(T)*481.3711*1.e5 #!!!!!!!! Fix me. Units hard coded.


        # Bottom of the LFM grid. NB we only use prm.rmin here
        # to get the vicinity of the grid. The actual boundary (rmin
        # variable is computed from the MAS grid)
        r_br = vars['mas']['br']['r']
        rmind = flatnonzero(r_br>=prm.rmin*prm.scale)[0]-1
        rmin = r_br[rmind]


        # LFM GRID
        sg = pyLTR.Grids.SphereGrid((ni,nj,nk))
        (P,T,R) = sg.ptrCorner(rmin=rmin,rmax=prm.rmax*prm.scale,thetamin=prm.thetamin)
        (Pc,Tc,Rc) = sg.ptrCenter()
        (x,y,z) = sg.xyzCorner()
        (xc,yc,zc) = sg.xyzCenter()

        # ghost cell coordinates
        (xg,yg,zg,xcg,ycg,zcg) = [2*arr[:,:,[0]]-arr[:,:,1:prm.NO2+1] for arr in (x,y,z,xc,yc,zc)]
        # Here we assume grid radii don't vary with theta and phi, at least in the ghost region
        rcg = sqrt(xcg[0,0,:]**2+ycg[0,0,:]**2+zcg[0,0,:]**2)
        rg = sqrt(xg[0,0,:]**2+yg[0,0,:]**2+zg[0,0,:]**2)


        # interpolate in radius (get coefficients here)
        for var_name in vars['mas']:
            interp[var_name]={};
            interp[var_name]['r']=vars['mas'][var_name]['r']
            mas2lfm.util.radial_interp(interp,var_name,rcg)

        # we also need to interpolate br to i-faces !!! 092514 VGM:
        #                                               no, we don't, since LFM bottom coincides
        #                                               with spherical shell where MAS Br is defined
        # radial_interp(vars['mas'],'br',rg,'r_interp_coefs_iface','r_interp_above_ind_iface')


        if prm.dumpInit:
            # make some aliases for simplicity
            br_mas = vars['mas']['br']['data']
            p_br   = vars['mas']['br']['phi']
            t_br   = vars['mas']['br']['theta']
            br = br_mas[:,:,rmind]
        
            fbi      = interpolate.RectBivariateSpline(p_br,t_br,br,kx=1,ky=1)  
            #        fbi      = interpolate.RectBivariateSpline(p_br,-cos(t_br),br,kx=1,ky=1)  # doing the same on -cos(theta) x phi grid.
            ###############################################

            lfmh = pyLTR.Tools.lfmstartup.lfmstartup(os.path.join(prm.dirInitLFMfile,prm.initLFMfile),(ni,nj,nk))
            lfmh.open()
            bi_lfm = fbi(Pc[:,0,0],Tc[0,:,0])
            #            bi_lfm = fbi(Pc[:,0,0],-cos(Tc[0,:,0]))
                
            # Set variables
            lfmh.writeVar('X_grid',x)
            lfmh.writeVar('Y_grid',y)
            lfmh.writeVar('Z_grid',z)
                
            lfmh.writeVar('rho_',400.*1.67e-24*ones((nk+1,nj+1,ni+1)))
            lfmh.writeVar('c_',  5.e6*ones((nk+1,nj+1,ni+1)))
            lfmh.writeVar('vx_',4.e7*sin(Tc)*cos(Pc))
            lfmh.writeVar('vy_',4.e7*sin(Tc)*sin(Pc))
            lfmh.writeVar('vz_',4.e7*cos(Tc))

            tmp    = diff(T[:,:,0],axis=1)
            dtheta = 0.5*(tmp[:-1,:]+tmp[1:,:])
            tmp    = diff(P[:,:,0],axis=0)
            dphi   = 0.5*(tmp[:,:-1]+tmp[:,1:])

            # note the fancy dstack to fill in the last i-face with correct bi values
            bi = dstack(R.shape[2]*[bi_lfm*rmin**2*sin(Tc[:,:,0])*dtheta*dphi])
            lfmh.writeVar('bi_',bi)
            lfmh.close()

        ### End of stuff done on first iteration only
    # dump cr.bc
    if prm.dumpBC:
        if prm.masFrame=='rotating': 
            vars['mas']['vp']['data']+=Vphi_solar


        vars['lfm']={}

        # some aliases
        for var_name in vars['mas']:
            var = vars['mas'][var_name]['data']
            r  = vars['mas'][var_name]['r']
            p  = vars['mas'][var_name]['phi']
            t  = vars['mas'][var_name]['theta']
            ind_a  = interp[var_name]['r_interp_above_ind']
            q      = interp[var_name]['r_interp_coefs']

            # interpolate lineary in the radial direction
            var_rintrp = (1.-q)*var[:,:,ind_a]+q*var[:,:,ind_a-1]

            # now bilinear interpolation in theta-phi
            tmp=[]
            for i in range(prm.NO2):

                fintrp    = interpolate.RectBivariateSpline(p,t,var_rintrp[:,:,i],kx=1,ky=1) 
#                fintrp    = interpolate.RectBivariateSpline(p,-cos(t),var_rintrp[:,:,i],kx=1,ky=1) 

                tmp.append(fintrp(Pc[:,0,0],Tc[0,:,0]))
#                tmp.append(fintrp(Pc[:,0,0],-cos(Tc[0,:,0])))
            vars['lfm'][var_name] = dstack(tmp)

            # finally interpolate bt and bp to the corresponding faces
            if var_name in ['bt','bp']:
                tmp=[]
                for i in range(prm.NO2):
                    fintrp    = interpolate.RectBivariateSpline(p,t,var_rintrp[:,:,i],kx=1,ky=1) 
#                    fintrp    = interpolate.RectBivariateSpline(p,-cos(t),var_rintrp[:,:,i],kx=1,ky=1) 
                    if var_name == 'bt':
                        tmp.append(fintrp(Pc[:,0,0],T[0,:,0]))
#                        tmp.append(fintrp(Pc[:,0,0],-cos(T[0,:,0])))
                    elif var_name == 'bp':
                        tmp.append(fintrp(P[:,0,0],Tc[0,:,0]))
#                        tmp.append(fintrp(P[:,0,0],-cos(Tc[0,:,0])))
                vars['lfm'][var_name+'_face'] = dstack(tmp)

        vr = vars['lfm']['vr']
        vt = vars['lfm']['vt']
        vp = vars['lfm']['vp']

        br = vars['lfm']['br']
        bt = vars['lfm']['bt']
        bp = vars['lfm']['bp']

        if prm.masFakeRotation:
            bp += -2*pi/(prm.Tsolar*24.*3600.)*rcg[0]*sin(Tc[:,:,[0]])*br[:,:,[0]]/vr[:,:,[0]]

        # note, we're just aliasing the vars here for brevity (should
        # just eliminate the cartesian components altogether because
        # we do the conversion inside boundaries.F now note also we
        # redefine bp and bt below, but that breaks the link between,
        # e.g., bz and bp, so taht bx holds the old bp and bp holds
        # the new one as define below (interpolated to bp_face
        bx,by,bz = br,bt,bp  #sph2cart( (br,bt,bp),Tc[:,:,[0]],Pc[:,:,[0]] )
        vx,vy,vz = vr,vt,vp  #sph2cart( (vr,vt,vp),Tc[:,:,[0]],Pc[:,:,[0]] )

        rho  = vars['lfm']['rho']
        temp = vars['lfm']['t']
        cs = sqrt(prm.gamma*1.38e-23*temp/1.67e-27)*1.e2   # in cm/s

        bp = vars['lfm']['bp_face']
        bt = vars['lfm']['bt_face']

        if prm.masFakeRotation:
            # need to shift Br and Vr to the k-face to add extra Bphi
            # the good news is that we deal with LFM variables now
            # so everything is at least at the same r
            for i in range(prm.NO2):
##### CHANGES
##########                fbr = interpolate.RectBivariateSpline(Pc[:,0,0],Tc[0,:,0],br[:nk,:nj,i])
##########                fvr = interpolate.RectBivariateSpline(Pc[:,0,0],Tc[0,:,0],vr[:nk,:nj,i])
                fbr = interpolate.RectBivariateSpline(Pc[:,0,0],Tc[0,:,0],br[:nk,:nj,0])
                fvr = interpolate.RectBivariateSpline(Pc[:,0,0],Tc[0,:,0],vr[:nk,:nj,0])
#                fbr = interpolate.RectBivariateSpline(Pc[:,0,0],-cos(Tc[0,:,0]),br[:nk,:nj,i])
#                fvr = interpolate.RectBivariateSpline(Pc[:,0,0],-cos(Tc[0,:,0]),vr[:nk,:nj,i])
                br_bp_face = fbr(P[:,0,0],Tc[0,:,0])
                vr_bp_face = fvr(P[:,0,0],Tc[0,:,0])
#                br_bp_face = fbr(P[:,0,0],-cos(Tc[0,:,0]))
#                vr_bp_face = fvr(P[:,0,0],-cos(Tc[0,:,0]))
#                bp[:-1,:,i] += -2*pi/(Tsolar*24.*3600.)*rcg[i]*scale*sin(Tc[:,:,0])*br_bp_face[:-1,:]/vr_bp_face[:-1,:]
##########                bp[:-1,:,i] = -2*pi/(Tsolar*24.*3600.)*rcg[i]*scale*sin(Tc[:,:,0])*br_bp_face[:-1,:]/vr_bp_face[:-1,:]
                bp[:-1,:,i] += -2*pi/(prm.Tsolar*24.*3600.)*rcg[0]*sin(Tc[:,:,0])*br_bp_face[:-1,:]/vr_bp_face[:-1,:]

#        bt=zeros_like(bt)
#        by=zeros_like(by)
#        vy=zeros_like(vy)
#        vz=zeros_like(vz)
        for i in range(prm.NO2):
            # note, by indexing below, we're making sure we are dumping arrays with shape (nk,nj)
            # instead of thinking whether they are already of this size
            out = array([ bt[:nk,:nj,i].T.ravel(),
                          bp[:nk,:nj,i].T.ravel(),
                          bx[:nk,:nj,i].T.ravel(),
                          by[:nk,:nj,i].T.ravel(),
                          bz[:nk,:nj,i].T.ravel(),
                          vx[:nk,:nj,i].T.ravel(),
                          vy[:nk,:nj,i].T.ravel(),
                          vz[:nk,:nj,i].T.ravel(),
                          rho[:nk,:nj,i].T.ravel(),
                          cs[:nk,:nj,i].T.ravel()])

            # Note we always start with 0 in the timelabel; may need to be changed.                  
            savetxt(os.path.join(prm.dirInitLFMfile,'innerbc_%03d_%d.dat'%(timeLabel-timeLabels[0]+1,i)),out.T,
                    fmt=['%13.8f','%13.8f','%13.8f','%13.8f','%13.8f',
                         '%15.5e','%15.5e','%15.5e',
                         '%14.5e','%14.5e'],
                    delimiter='')
