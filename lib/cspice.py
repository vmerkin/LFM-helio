import numpy as np
import datetime

def cspice(fname):
    dati=[]
    x=[]
    y=[]
    z=[]
    r=[]
    phi=[]
    theta=[]

    for line in open(fname,'r'):
        fields = line.strip().split()
        year  = int(fields[0].strip().split('-')[0])
        doy   = int(fields[0].strip().split('-')[1])
        hr    = int(fields[1].strip().split(':')[0])
        mn   = int(fields[1].strip().split(':')[1])
        dati.append(datetime.datetime(year,1,1,hr,mn)+datetime.timedelta(days=doy-1))

        x.append(float(fields[2]))
        y.append(float(fields[3]))
        z.append(float(fields[4]))

        r.append(float(fields[5]))
        theta.append(float(fields[6]))
        phi.append(float(fields[7]))
        
    return(dati,x,y,z,r,theta,phi)
