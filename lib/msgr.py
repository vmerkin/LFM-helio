import numpy as np
import datetime

def msgr(msgr_file,ave_1h=True):
    r=[]
    lat=[]
    lon=[]
    br=[]
    bt=[]
    bn=[]
    dati=[]
    for line in open(msgr_file,'r'):
        fields = line.strip().split()

        year = int(fields[0])
        doy  = int(fields[1])
        hr   = int(fields[2])
        mn   = int(fields[3])
        r.append(float(fields[7]))
        lat.append(float(fields[8]))
        lon.append(float(fields[9]))
        br.append(float(fields[10]))
        bt.append(float(fields[11]))
        bn.append(float(fields[12]))
        
        dati.append(datetime.datetime(year,1,1,hr,mn)+datetime.timedelta(days=doy-1))

    
    if ave_1h:    
        # assuming 1-min data. Average to 1-h
        n=len(lat)-len(lat)%60
        br_1h=np.mean(np.array(br)[:n].reshape(-1,60),1)
        bt_1h=np.mean(np.array(bt)[:n].reshape(-1,60),1)
        bn_1h=np.mean(np.array(bn)[:n].reshape(-1,60),1)
    
        return (dati[30:n+30:60],br_1h,bt_1h,bn_1h)
    else:
        return (dati,np.array(br),np.array(bt),np.array(bn))

def msgr_v1(msgr_file,year):
    data = np.loadtxt(msgr_file)

    start_dati = datetime.datetime(year,1,1,0)
    dati = [start_dati+datetime.timedelta(days=d-1) for d in data[:,0]]
    
    br=data[:,4]
    bt=data[:,5]
    bn=data[:,6]

    return (dati,br,bt,bn)

