from pylab import *
import datetime
earth = loadtxt('/Users/merkivg1/work/LFM-helio_2.0/SCDATA/helios_13054.lst-earth.txt',skiprows=1)

cr_start = datetime.datetime(2010,1,1)
earth_days = array([(datetime.datetime(int(y),1,1)+datetime.timedelta(days=d)-cr_start).days  for (y,d) in zip(earth[:,0],earth[:,1])],dtype=float) % 27.27  # comment out the % 27.27 part to get monotonic phi

earth_phi = 360.*(1-earth_days/27.27)

phi=linspace(0,360,192)





"""
inside = 0.<=earth_phi<=360.
from scipy import interpolate
f=interpolate.interp1d(earth_phi[inside],earth[inside,3])
earth_theta = f(phi)

outside_right = earth_phi>360.
outside_left  = earth_phi<0.








f=interpolate.UnivariateSpline(earth_phi,earth[:,3],k=1)
earth_theta = f(phi)

"""
