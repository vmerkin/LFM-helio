from pyhdf.SD import SD,SDC
import os
import ConfigParser

mas_units = {'length':6.96e10,
             'time':1445.87,
             'vr':481.3711*1.e5,
             'vt':481.3711*1.e5,
             'vp':481.3711*1.e5,
             'n':1.e8,
             'rho':1.6726e-16,
             'p':0.3875717,
             't':2.807067e7,
             'br':2.2068908,
             'bt':2.2068908,
             'bp':2.2068908}

mas_var_names = ['t','rho','vt','vp','vr','bt','bp','br']

def read_var(fname,varname,normalized=False):
    f     = SD(fname,SDC.READ)
    phi   = f.select('fakeDim0')[:]
    theta = f.select('fakeDim1')[:]
    r     = f.select('fakeDim2')[:]
    var   = f.select('Data-Set-2')[:]
    f.end()

    if normalized:
        return(phi,theta,r,var)
    else:
        return(phi,theta,r*mas_units['length'],var*mas_units[varname])

def units():
    return(mas_units)


def read_all_vars(dir,timeLabel,normalized=False):
    vars = {}

    for var_name in mas_var_names:
        vars[var_name]={}
        (vars[var_name]['phi'],
         vars[var_name]['theta'],
         vars[var_name]['r'],
         vars[var_name]['data']) = read_var(os.path.join(dir,var_name+'%03d.hdf'%timeLabel),var_name,normalized)

    return(vars)

def set_plot_limits(v,plotConfigFile = 'plot.config'):
    plotConfig = ConfigParser.ConfigParser()
    plotConfig.read(plotConfigFile)

    # assume here v has the necessary keys
    # bascially, this needs to be called after mas.read_all_varaibles
    try:
        v['bt']['lims'] = plotConfig.getfloat('MAS','bt_lim_lo'),plotConfig.getfloat('MAS','bt_lim_hi')
        v['bt']['fmt'] = plotConfig.get('MAS','bt_fmt')

        v['bp']['lims'] = plotConfig.getfloat('MAS','bp_lim_lo'),plotConfig.getfloat('MAS','bp_lim_hi')
        v['bp']['fmt'] = plotConfig.get('MAS','bp_fmt')

        v['br']['lims'] = plotConfig.getfloat('MAS','br_lim_lo'),plotConfig.getfloat('MAS','br_lim_hi')
        v['br']['fmt'] = plotConfig.get('MAS','br_fmt')

        v['vr']['lims'] = plotConfig.getfloat('MAS','vr_lim_lo'),plotConfig.getfloat('MAS','vr_lim_hi')
        v['vr']['fmt'] = plotConfig.get('MAS','vr_fmt')

        v['vt']['lims'] = plotConfig.getfloat('MAS','vt_lim_lo'),plotConfig.getfloat('MAS','vt_lim_hi')
        v['vt']['fmt'] = plotConfig.get('MAS','vt_fmt')

        v['vp']['lims'] = plotConfig.getfloat('MAS','vp_lim_lo'),plotConfig.getfloat('MAS','vp_lim_hi')
        v['vp']['fmt'] = plotConfig.get('MAS','vp_fmt')

        v['rho']['lims'] = plotConfig.getfloat('MAS','rho_lim_lo'),plotConfig.getfloat('MAS','rho_lim_hi')
        v['rho']['fmt'] = plotConfig.get('MAS','rho_fmt')

        v['t']['lims'] = plotConfig.getfloat('MAS','t_lim_lo'),plotConfig.getfloat('MAS','t_lim_hi')
        v['t']['fmt'] = plotConfig.get('MAS','t_fmt')
    except KeyError:
        print('Some keys are missing from the variable dictionary.')
