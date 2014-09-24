"""
This script is adapted from ~/work/LFM-helio_2.0/mascme2lfm_helio_coronal_interp.py.
Here I'm going to try to keep things clean and modular.
"""
#----------- PARSE ARGUMENTS ---------#
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('ConfigFileName',help='The name of the configuration file to use',default='startup.config')
args = parser.parse_args()
#----------- PARSE ARGUMENTS ---------#

# Read params from config file
import mas2lfm
import mas2lfm.util
params = mas2lfm.util.configDict(args.ConfigFileName)

# Pretty print config parameters
import pprint
print('=============================\n Configuration parameters \n=============================')
pprint.pprint(params)
print('=============================\n Configuration parameters \n=============================')





import os,sys,glob
import mas
from numpy import linspace,pi,meshgrid,sin,cos,zeros,ones,dstack,diff,sqrt,array,savetxt,flatnonzero,insert,asarray,zeros_like
import lfmpy
from scipy import interpolate
import time


if params['masTimeLabel']=='all':
    timeLabels = [int(os.path.basename(s).split('.')[0][-3:]) for s in sorted(glob.glob(os.path.join(masdir,'rho*.hdf')))]
else:
    masTimeLabel = int(params['masTimeLabel'])
    timeLabels = [masTimeLabel]
    initLFMfile = fileName


import sys
sys.exit(0)

vars={}
vars['mas']={}
mas_var_names = ['t','rho','vt','vp','vr','bt','bp','br']
mas_var_units = ['temperature','mass density',
                 'velocity','velocity','velocity',
                 'magnetic field','magnetic field','magnetic field']
for var_name in mas_var_names: vars['mas'][var_name]={}


for timeLabel in timeLabels:
    print(timeLabel)
    ######## READ IN MAS DATA #####################
    for var_name,var_unit in zip(mas_var_names,mas_var_units):
        (vars['mas'][var_name]['phi'],
         vars['mas'][var_name]['theta'],
         vars['mas'][var_name]['r'],
         vars['mas'][var_name]['data']) = mas.read_var(os.path.join(masdir,var_name+'%03d.hdf'%timeLabel),var_unit)

    # p_t,t_t,r_t,temp_mas      = mas.read_var(os.path.join(masdir,'t%03d.hdf'%timeLabel),'temperature')
    # p_rho,t_rho,r_rho,rho_mas = mas.read_var(os.path.join(masdir,'rho%03d.hdf'%timeLabel),'mass density')
    # p_vr,t_vr,r_vr,vr_mas     = mas.read_var(os.path.join(masdir,'vr%03d.hdf'%timeLabel),'velocity')
    # p_vt,t_vt,r_vt,vt_mas     = mas.read_var(os.path.join(masdir,'vt%03d.hdf'%timeLabel),'velocity')    
    # p_vp,t_vp,r_vp,vp_mas     = mas.read_var(os.path.join(masdir,'vp%03d.hdf'%timeLabel),'velocity')    
    # p_br,t_br,r_br,br_mas     = mas.read_var(os.path.join(masdir,'br%03d.hdf'%timeLabel),'magnetic field')
    # p_bt,t_bt,r_bt,bt_mas     = mas.read_var(os.path.join(masdir,'bt%03d.hdf'%timeLabel),'magnetic field')    
    # p_bp,t_bp,r_bp,bp_mas     = mas.read_var(os.path.join(masdir,'bp%03d.hdf'%timeLabel),'magnetic field')    


    # We only assum here that the LFM bottom boundary coincides with the shell where MAS br is defined. 
    # We don't have to do this: could define rmin arbitrarily and then linearly interpolate.
    #
    # Do it only for the first MAS file. If we did not have to find rmin here,
    # we could move the LFM grid and ghost defs above the timeLabel loop

    if timeLabel == timeLabels[0]:

        # transform the azimuthal velocity to inertial frame
        R,T=meshgrid(vars['mas']['vp']['r'],vars['mas']['vp']['theta'])
        # using numpy broadcasting here
        Vphi_solar = 2*pi/(Tsolar*24.*3600./mas.time_unit())*R*sin(T)*481.3711*1.e5


        # Bottom of the LFM grid
        r_t = vars['mas']['t']['r']
        rmind = flatnonzero(r_t>=rmin)[0]-1
        rmin = r_t[rmind]

        # LFM GRID
        grid = lfmpy.grids((ni,nj,nk))
        (P,T,R,Pc,Tc,Rc,x,y,z,xc,yc,zc) = grid.getGrid(rmin=rmin*scale,rmax=rmax*scale,thetamin=thetamin)

        # ghost cell coordinates
        (xg,yg,zg,xcg,ycg,zcg) = [2*arr[:,:,[0]]-arr[:,:,1:NO2+1] for arr in (x,y,z,xc,yc,zc)]
        # Here we assume grid radii don't vary with theta and phi, at least in the ghost region
        rcg = sqrt(xcg[0,0,:]**2+ycg[0,0,:]**2+zcg[0,0,:]**2)/scale
        rg = sqrt(xg[0,0,:]**2+yg[0,0,:]**2+zg[0,0,:]**2)/scale


        # interpolate in radius (get coefficients here)
        for var_name in mas_var_names:
            radial_interp(vars['mas'],var_name,rcg)

        # we also need to interpolate br to i-faces
        # radial_interp(vars['mas'],'br',rg,'r_interp_coefs_iface','r_interp_above_ind_iface')
    if masFrame=='rotating': 
        vars['mas']['vp']['data']+=Vphi_solar

    if dumpInit:
        br_mas = vars['mas']['br']['data']
        p_br   = vars['mas']['br']['phi']
        t_br   = vars['mas']['br']['theta']
        br = br_mas[:,:,rmind]
#### CHANGES
        fbi      = interpolate.RectBivariateSpline(p_br,t_br,br,kx=1,ky=1)  
#        fbi      = interpolate.RectBivariateSpline(p_br,-cos(t_br),br,kx=1,ky=1)  

        ###############################################
        with lfmpy.startup(os.path.join(masdir,initLFMfile),(ni,nj,nk)) as lfmh:
#### CHANGES
            bi_lfm = fbi(Pc[:,0,0],Tc[0,:,0])
#            bi_lfm = fbi(Pc[:,0,0],-cos(Tc[0,:,0]))
            
            # Set variables
            lfmh.setVar('X_grid',x)
            lfmh.setVar('Y_grid',y)
            lfmh.setVar('Z_grid',z)
            
            lfmh.setVar('rho_',400.*1.67e-24*ones((nk+1,nj+1,ni+1)))
            lfmh.setVar('c_',  5.e6*ones((nk+1,nj+1,ni+1)))
            lfmh.setVar('vx_',4.e7*sin(Tc)*cos(Pc))
            lfmh.setVar('vy_',4.e7*sin(Tc)*sin(Pc))
            lfmh.setVar('vz_',4.e7*cos(Tc))

            tmp    = diff(T[:,:,0],axis=1)
            dtheta = 0.5*(tmp[:-1,:]+tmp[1:,:])
            tmp    = diff(P[:,:,0],axis=0)
            dphi   = 0.5*(tmp[:,:-1]+tmp[:,1:])

            # note the fancy dstack to fill in the last i-face with correct bi values
            bi = dstack(R.shape[2]*[bi_lfm*(rmin*scale)**2*sin(Tc[:,:,0])*dtheta*dphi])
            lfmh.setVar('bi_',bi)


    # dump cr.bc
    if dumpBC:
        vars['lfm']={}
        for var_name in mas_var_names:
            var = vars['mas'][var_name]['data']
            r  = vars['mas'][var_name]['r']
            p  = vars['mas'][var_name]['phi']
            t  = vars['mas'][var_name]['theta']
            ind_a  = vars['mas'][var_name]['r_interp_above_ind']
            q      = vars['mas'][var_name]['r_interp_coefs']

            # interpolate lineary in the radial direction
            var_rintrp = (1.-q)*var[:,:,ind_a]+q*var[:,:,ind_a-1]

            # now bilinear interpolation in theta-phi
            tmp=[]
            for i in range(NO2):
#### CHANGES
                fintrp    = interpolate.RectBivariateSpline(p,t,var_rintrp[:,:,i],kx=1,ky=1) 
#                fintrp    = interpolate.RectBivariateSpline(p,-cos(t),var_rintrp[:,:,i],kx=1,ky=1) 
                tmp.append(fintrp(Pc[:,0,0],Tc[0,:,0]))
#                tmp.append(fintrp(Pc[:,0,0],-cos(Tc[0,:,0])))
            vars['lfm'][var_name] = dstack(tmp)

            # finally interpolate bt and bp to the corresponding faces
            if var_name in ['bt','bp']:
                tmp=[]
                for i in range(NO2):
#### CHANGES
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

        if masFakeRotation:
#            bp += -2*pi/(Tsolar*24.*3600.)*rcg*scale*sin(Tc[:,:,[0]])*br/vr
#            bp = -2*pi/(Tsolar*24.*3600.)*rcg*scale*sin(Tc[:,:,[0]])*br/vr
            bp += -2*pi/(Tsolar*24.*3600.)*rcg[0]*scale*sin(Tc[:,:,[0]])*br[:,:,[0]]/vr[:,:,[0]]

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
        cs = sqrt(gamma*1.38e-23*temp/1.67e-27)*1.e2   # in cm/s

        bp = vars['lfm']['bp_face']
        bt = vars['lfm']['bt_face']

        if masFakeRotation:
            # need to shift Br and Vr to the k-face to add extra Bphi
            # the good news is that we deal with LFM variables now
            # so everything is at least at the same r
            for i in range(NO2):
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
                bp[:-1,:,i] += -2*pi/(Tsolar*24.*3600.)*rcg[0]*scale*sin(Tc[:,:,0])*br_bp_face[:-1,:]/vr_bp_face[:-1,:]

#        bt=zeros_like(bt)
#        by=zeros_like(by)
#        vy=zeros_like(vy)
#        vz=zeros_like(vz)
        for i in range(NO2):
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
                  
            savetxt(os.path.join(masdir,'innerbc_%03d_%d.dat'%(timeLabel,i)),out.T,
                    fmt=['%13.8f','%13.8f','%13.8f','%13.8f','%13.8f',
                         '%15.5e','%15.5e','%15.5e',
                         '%14.5e','%14.5e'],
                    delimiter='')
