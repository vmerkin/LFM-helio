from pyhdf.SD import SD,SDC

units = {'length':6.96e10,
         'time':1445.87,
         'velocity':481.3711*1.e5,
         'number density':1.e8,
         'mass density':1.6726e-16,
         'pressure':0.3875717,
         'temperature':2.807067e7,
         'magnetic field':2.2068908}

def read_var(fname,varname):
    f     = SD(fname,SDC.READ)
    phi   = f.select('fakeDim0')[:]
    theta = f.select('fakeDim1')[:]
    r     = f.select('fakeDim2')[:]
    var   = f.select('Data-Set-2')[:]
    f.end()

    return(phi,theta,r*units['length'],var*units[varname])

def time_unit():
    return(units['time'])
