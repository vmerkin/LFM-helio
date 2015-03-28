from pylab import *

def map_helios(filename,cr_start_string):
    """
    Get Earth position based on the beginning time of the carrrington
    map and the time in the correspondning file for the specific
    spacecrat. Output HGI longitude from the file to get the relative
    separation between spacecraft -- this way if we know the EArth
    position, we can get other spacecraft HGI longitude by taking the
    differene. See __main__ below for an example.

    Parameters:
    [In] filename: File in the format of NASA helios http://omniweb.gsfc.nasa.gov/coho/helios/heli.html
    [In] cr_start_string: Data/time start of the CR in the following format "%Y %m %d %H %M %S"
    [Out] phi: longitude in degrees between 0 and 360
    [Out] HGI latitude: from -90 (south) to +90 (north)
    [Out] HGI longitude: use for relative separation between spacecraft
    """
    import datetime
    data = loadtxt(filename,skiprows=1)

    cr_start=datetime.datetime.strptime(cr_start_string,"%Y %m %d %H %M %S")
    days = array([(datetime.datetime(int(y),1,1)+datetime.timedelta(days=d)-cr_start).days  for (y,d) in zip(data[:,0],data[:,1])],dtype=float) % 27.27  # comment out the % 27.27 part to get monotonic phi
    
    phi = 360.*(1-days/27.27)
    
    return(phi,data[:,3],data[:,4])

if __name__ == "__main__":
#    earth = map_helios('/glade/u/home/vgm/LFM-helio_2.0/SCDATA/helios_13054.lst-earth-short.txt',"2010 1 12 7 45 36")
#    mess  = map_helios('/glade/u/home/vgm/LFM-helio_2.0/SCDATA/helios_13066.lst-messenger-short.txt',"2010 1 12 7 45 36")

    earth = map_helios('/glade/u/home/vgm/LFM-helio_2.0/SCDATA/helios_2008_20-50-earth.txt',"2008 1 7 7 45 36")
    mess  = map_helios('/glade/u/home/vgm/LFM-helio_2.0/SCDATA/helios_2008_20-50-messenger.txt',"2008 1 7 7 45 36")


    mess_phi = earth[0] - (earth[2]-mess[2])
    





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
