import numpy as np

def sph2cart(a,theta,phi):
    ar,at,ap = a

    ax = ar*cos(phi)*sin(theta)+at*cos(phi)*cos(theta)-ap*sin(phi)
    ay = ar*sin(phi)*sin(theta)+at*sin(phi)*cos(theta)+ap*cos(phi)
    az = ar*cos(theta)-at*sin(theta)

    return(ax,ay,az)

def radial_interp(vars_from,var_name,radii_to,coef_label='r_interp_coefs',ind_label='r_interp_above_ind'):
    q=[]
    ind=[]
    for r_to in radii_to:    # interpolating to LFM cell centers
        r_from = vars_from[var_name]['r']
        ind_r_a = np.flatnonzero(r_from>=r_to)[0]     # MAS radius above
        r_a = r_from[ind_r_a]
        r_b = r_from[ind_r_a-1]
        q.append( (r_a-r_to)/(r_a-r_b) )
        ind.append(ind_r_a)
        
    vars_from[var_name][coef_label]=np.asarray(q,dtype='float')
    vars_from[var_name][ind_label]=np.asarray(ind,dtype='int')
