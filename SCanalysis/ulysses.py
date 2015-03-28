from pylab import *

def ulysses(vr,br,rho,cs,Rc,thetac,phic,filename = '/glade/u/home/vgm/LFM-helio_2.0/SCDATA/Ulysses/hourav_1995.asc'):
    ulysses=loadtxt(filename)
#    ind=logical_and(ulysses[:,1]>=27.5074,ulysses[:,1]<=27.5074+27.27)
    ind=ulysses[:,1]<=140.
#    ind=slice(ulysses[:,0].size)
    ulysses_lon = ulysses[ind,7]
    ulysses_lat = ulysses[ind,6]
    ulysses_rad = ulysses[ind,5]

    lfm_vr=[]
    lfm_br=[]
    lfm_n=[]
    lfm_t=[]

    for (indd,doy) in enumerate(ulysses[ind,1]):
        i=argmin( abs(Rc[0,0,:]/6.955e10/215.-ulysses_rad[indd]))
        j=argmin( abs(90.-180./pi*thetac[0,:,0]-ulysses_lat[indd]) )
        k=argmin( abs(180./pi*phic[:,0,0]-ulysses_lon[indd]) )
        lfm_vr.append(vr[k,j,i])
        lfm_br.append(br[k,j,i])
        lfm_n.append(rho[k,j,i]/1.67e-24)
        lfm_t.append(cs[k,j,i]**2*1.67e-24/1.6667/1.38e-16)


    return(ulysses[ind],lfm_vr,lfm_br,lfm_n,lfm_t)


# figure()
# plot(ulysses[ind,1],ulysses[ind,12])
# plot(ulysses[ind,1],array(lfm_vr)*1.e-5)
# title('Vr, blue-ulysses')

# figure()
# plot(ulysses[ind,1],ulysses[ind,8])
# plot(ulysses[ind,1],ulysses[ind,9],'b--')
# plot(ulysses[ind,1],lfm_n)
# title('Density')

# figure()
# plot(ulysses[ind,1],ulysses[ind,10])
# plot(ulysses[ind,1],ulysses[ind,11],'b--')
# plot(ulysses[ind,1],array(lfm_t))
# axes().set_yscale('log')
# ylim(1e4,1e6)
