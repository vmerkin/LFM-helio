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
        ind_r_a = flatnonzero(r_from>=r_to)[0]     # MAS radius above
        r_a = r_from[ind_r_a]
        r_b = r_from[ind_r_a-1]
        q.append( (r_a-r_to)/(r_a-r_b) )
        ind.append(ind_r_a)
        
    vars_from[var_name][coef_label]=asarray(q,dtype='float')
    vars_from[var_name][ind_label]=asarray(ind,dtype='int')


def configDict(ConfigFileName):
    import ConfigParser
    config = ConfigParser.ConfigParser()
    config.read(ConfigFileName)

    params = {}

    params['ni'] = config.getint('Dimensions','NI')
    params['nj'] = config.getint('Dimensions','NJ')
    params['nk'] = config.getint('Dimensions','NK')
    params['fileName'] = config.get('OutputFileName','Prefix')+'_%dx%dx%d_mhd_0000000.hdf'%(params['ni'],
                                                                                            params['nj'],
                                                                                            params['nk'])
    
    params['scale'] = config.getfloat('Normalization','Xscale')
    params['rmin'] = config.getfloat('GridSpecs','RMIN')
    params['rmax'] = config.getfloat('GridSpecs','RMAX')
    params['thetamin'] = config.getfloat('GridSpecs','THETAMIN')

    params['gamma'] = config.getfloat('Constants','gamma')
    params['NO2']   = config.getint('Constants','NO2')
    params['Tsolar'] = config.getfloat('Constants','Tsolar')

    params['masdir'] = config.get('MAS','masdir')
    params['masTimeLabel'] = config.get('MAS','masTimeLabel')
    params['masFrame'] = config.get('MAS','masFrame')
    params['masFakeRotation'] = config.getboolean('MAS','masFakeRotation')

    params['dumpInit'] = config.getboolean('DUMPS','init')
    params['dumpBC']   = config.getboolean('DUMPS','BC')

    return(params)
