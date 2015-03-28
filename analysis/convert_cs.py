from pylab import *
import glob

gamma = 1.5
b0 = 0.003

files = sorted(glob.glob('innerbc_*.dat'))

for f in files:
    d=loadtxt(f)

    if f==files[0]:
#        rho = d[:,-2]
#        br  = d[:,4]
        cs  = d[:,-1].copy()
#        d[:,-1] = sqrt(cs**2+gamma*(b0**2-br**2)/8./pi/rho)
    else:
        d[:,-1] = cs
        savetxt('full_pressure_balance_constantcs_radially/'+f,d,fmt=['%15.5e']*12,delimiter='')


