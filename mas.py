from pyhdf.SD import SD,SDC
import os

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
    mas_var_names = ['t','rho','vt','vp','vr','bt','bp','br']

    vars = {}

    for var_name in mas_var_names:
        vars[var_name]={}
        (vars[var_name]['phi'],
         vars[var_name]['theta'],
         vars[var_name]['r'],
         vars[var_name]['data']) = read_var(os.path.join(dir,var_name+'%03d.hdf'%timeLabel),var_name,normalized)

    return(vars)

    
