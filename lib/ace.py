import numpy as np
import datetime

def cdaweb_mfi(ACE_mf_file):
    # Mag. field
    ACE_mf_time=[]
    ACE_bx=[]
    ACE_by=[]
    ACE_bz=[]

    i=1
    for line in open(ACE_mf_file,'r'):
        i+=1
        if 'dd-mm-yyyy' in line: break
    start = i

    i = 0
    for line in open(ACE_mf_file,'r'):
        i+=1
        if i<start: continue
        if line[0]=='#': break

        fields = line.strip().split()
        ACE_mf_day  = int(fields[0].strip().split('-')[0])
        ACE_mf_mon  = int(fields[0].strip().split('-')[1])
        ACE_mf_year = int(fields[0].strip().split('-')[2])
        ACE_mf_hr   = int(fields[1].strip().split(':')[0])
        ACE_mf_min  = int(fields[1].strip().split(':')[1])
        ACE_mf_sec  = int(float(fields[1].strip().split(':')[2]))
        ACE_mf_time.append(datetime.datetime(ACE_mf_year, ACE_mf_mon, ACE_mf_day, ACE_mf_hr, ACE_mf_min, ACE_mf_sec))
        ACE_bx.append(float(fields[2]))
        ACE_by.append(float(fields[3]))
        ACE_bz.append(float(fields[4]))

    # invalidate bad data
    ACE_bx=np.array(ACE_bx)
    ACE_by=np.array(ACE_by)
    ACE_bz=np.array(ACE_bz)
    ACE_bx[np.abs(np.array(ACE_bx)) > 1.e20] = np.nan
    ACE_by[np.abs(np.array(ACE_by)) > 1.e20] = np.nan
    ACE_bz[np.abs(np.array(ACE_bz)) > 1.e20] = np.nan
    
    return (ACE_mf_time,ACE_bx,ACE_by,ACE_bz)

def cdaweb_swe(ACE_plasma_file):
    # Plasma
    ACE_plasma_time=[]
    ACE_v=[]
    ACE_n=[]

    i=1
    for line in open(ACE_plasma_file,'r'):
        i+=1
        if 'dd-mm-yyyy' in line: break
    start = i

    i = 0
    for line in open(ACE_plasma_file,'r'):
        i+=1
        if i<start: continue
        if line[0]=='#': break

        fields = line.strip().split()
        ACE_plasma_day  = int(fields[0].strip().split('-')[0])
        ACE_plasma_mon  = int(fields[0].strip().split('-')[1])
        ACE_plasma_year = int(fields[0].strip().split('-')[2])
        ACE_plasma_hr   = int(fields[1].strip().split(':')[0])
        ACE_plasma_min  = int(fields[1].strip().split(':')[1])
        ACE_plasma_sec  = int(float(fields[1].strip().split(':')[2]))
        ACE_plasma_time.append(datetime.datetime(ACE_plasma_year, ACE_plasma_mon, ACE_plasma_day, ACE_plasma_hr, ACE_plasma_min, ACE_plasma_sec))
        ACE_n.append(float(fields[2]))
        ACE_v.append(float(fields[3]))
        

    # invalidate bad data
    ACE_v=np.array(ACE_v)
    ACE_n=np.array(ACE_n)
    ACE_v[np.abs(np.array(ACE_v)) > 1.e20] = np.nan
    ACE_n[np.abs(np.array(ACE_n)) > 1.e20] = np.nan
    
    
    return (ACE_plasma_time,ACE_v,ACE_n)

