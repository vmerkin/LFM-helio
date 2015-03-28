from pylab import *
import glob

gamma = 1.5
b0 = 0.003

files = sorted(glob.glob('innerbc_*.dat'))

fullp=[]
for f in files:
    d=loadtxt(f)

    rho = d[:,-2]
    br  = d[:,4]
    cs  = d[:,-1]

    fullp.append(1./gamma*rho*cs**2+br**2/8./pi)


