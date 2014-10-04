from matplotlib import pyplot as plt
from matplotlib import cm as cm
from numpy import linspace,pi,meshgrid,sin,cos,zeros,ones,dstack,diff,sqrt,array,savetxt,flatnonzero,insert,where
from scipy import interpolate


#############################
###### PLOTTING ######
#############################

def plot_ij(vars,prm):
    nk,nj,ni=vars['lfm']['t'].shape
    R = vars['lfm']['R']
    theta = vars['lfm']['theta']
    Rc = vars['lfm']['Rc']
    thetac = vars['lfm']['thetac']

    kslice = prm.kslice
    mas_kslice=where(vars['mas'][prm.variable]['phi']>=float(kslice)/nk*2.*pi)[0][0]

    fig = plt.figure(figsize = (18,16))

    vmin,vmax = vars['mas'][prm.variable]['lims']
    fmt = vars['mas'][prm.variable]['fmt']

    ax1 = plt.subplot(311)
    ax2 = plt.subplot(312)
    ax3 = plt.subplot(313)

    # interpolate MAS onto LFM grid
    f=interpolate.RectBivariateSpline(vars['mas'][prm.variable]['theta'],
                                      vars['mas'][prm.variable]['r']/6.96e10, #FIX ME. UNITS!
                                      vars['mas'][prm.variable]['data'][mas_kslice,:,:],kx=1,ky=1)
    mas_interp=f(thetac[0,:,0],Rc[mas_kslice,0,:])
    
    p1=ax1.pcolormesh(R[kslice,:,:],pi-theta[kslice,:,:],mas_interp,vmin=vmin,vmax=vmax)
    ax1.xaxis.set_visible(False)
    ax1.set_title( ('MAS: Min/Max = '+fmt+'/'+fmt) % (mas_interp.min(),mas_interp.max()))

    p2 = ax2.pcolormesh(R[kslice,:,:],pi-theta[kslice,:,:],vars['lfm'][prm.variable][kslice,:,:],vmin=vmin,vmax=vmax)
    ax2.xaxis.set_visible(False)
    ax2.set_title( ('LFM: Min/Max = '+fmt+'/'+fmt) % 
                   (vars['lfm'][prm.variable].min(),vars['lfm'][prm.variable].max()))

    ax3 = plt.subplot(313)
    p3 = ax3.pcolormesh(R[kslice,:,:],pi-theta[kslice,:,:],vars['lfm'][prm.variable][kslice,:,:]/mas_interp-1,vmin=-0.05,vmax=0.05,cmap=cm.RdBu_r)
    ax3.set_xlabel('R')

    ax3.set_title( ('Ratio-1: Min/Max = '+fmt+'/'+fmt) % 
                   ((vars['lfm'][prm.variable]/mas_interp-1).min(),(vars['lfm'][prm.variable]/mas_interp-1).max()))

    for (ax,p) in zip([ax1,ax2,ax3],[p1,p2,p3]):
        ax.set_xlim(20,50)
        ax.set_ylim(0,pi)
        ax.set_ylabel(r'$\theta$')
        plt.colorbar(p,ax=ax).set_label(prm.variable)

    plt.show()

    
# This function may need a lot of cleanup. Seems to work but is clunky. Need optimization.
def plot_ik_jk(vars,prm):
    nk,nj,ni=vars['lfm']['t'].shape
    R = vars['lfm']['R']
    theta = vars['lfm']['theta']
    phi = vars['lfm']['phi']
    Rc = vars['lfm']['Rc']
    thetac = vars['lfm']['thetac']
    phic = vars['lfm']['phic']

    mas_islice = flatnonzero(vars['mas'][prm.variable]['r']/6.96e10>=prm.rslice)[0]
    lfm_islice = flatnonzero(Rc[0,0,:]/6.96e10>=prm.rslice)[0]

    if prm.slice == 'equatorial':
        mas_jslice = flatnonzero(vars['mas'][prm.variable]['theta']>=pi/2.)[0]
        lfm_jslice = flatnonzero(vars['lfm']['thetac']>=pi/2.)[0]
        lfm_selection = (slice(None),lfm_jslice,slice(lfm_islice,None))
        mas_selection = (slice(None),mas_jslice,slice(mas_islice,None))
    elif prm.slice == 'i':
        lfm_selection = (slice(None),slice(None),lfm_islice)
        mas_selection = (slice(None),slice(None),mas_islice)


    p_mas,t_mas = vars['mas'][prm.variable]['phi'],vars['mas'][prm.variable]['theta']

    mas_x = p_mas
    lfm_x = phi[lfm_selection]

    if prm.slice=='i':
        mas_y = pi-t_mas
        lfm_y = pi-theta[lfm_selection]
    elif prm.slice=='equatorial':
        mas_y = r_mas[mas_islice:]
        lfm_y = R[lfm_selection]
        

    fig = plt.figure(figsize = (18,16))
    
    vmin,vmax = vars['mas'][prm.variable]['lims']
    fmt = vars['mas'][prm.variable]['fmt']

    if prm.slice=='i':
        ax1 = plt.subplot(211,aspect='equal')
        ax2 = plt.subplot(212,aspect='equal')
    elif prm.slice=='equatorial':
        ax1 = plt.subplot(211)
        ax2 = plt.subplot(212)

    p1=ax1.pcolormesh(mas_x,mas_y,vars['mas'][prm.variable]['data'][mas_selection].T,vmin=vmin,vmax=vmax)
    ax1.xaxis.set_visible(False)
    ax1.set_title( ('MAS: Min/Max = '+fmt+'/'+fmt) % 
                   (vars['mas'][prm.variable]['data'][mas_selection].min(),vars['mas'][prm.variable]['data'][mas_selection].max()))

    p2 = ax2.pcolormesh(lfm_x,lfm_y,vars['lfm'][prm.variable][lfm_selection],vmin=vmin,vmax=vmax)
    ax2.set_xlabel(r'$\phi$')
    ax2.set_title( ('LFM: Min/Max = '+fmt+'/'+fmt) % 
                   (vars['lfm'][prm.variable].min(),vars['lfm'][prm.variable].max()))

    for ax in [ax1,ax2]:
        ax.set_xlim(0,2*pi)
        if prm.slice=='i':
            ax.set_ylim(0,pi)
            ax.set_ylabel(r'$\theta$')
        elif prm.slice=='equatorial':
            ax.set_ylim(prm.rslice,r_mas.max())
            ax.set_ylabel('R')
        plt.colorbar(p1,ax=ax).set_label(prm.variable)

    plt.show()

    
