
from netCDF4 import Dataset
import matplotlib.pyplot as plt

import glob

fnames=glob.glob("*201705*nc")
fnames=sorted(fnames)
fname='sgparsclkazrcloudsatC1.c1.20170730.000000.nc'

fnames2=['sgparsclkazrcloudsatC1.c1.20170707.000000.nc',
         'sgparsclkazrcloudsatC1.c1.20170708.000000.nc',
         'sgparsclkazrcloudsatC1.c1.20170607.000000.nc',
         'sgparsclkazrcloudsatC1.c1.20170608.000000.nc',
         'sgparsclkazrcloudsatC1.c1.20170614.000000.nc',
         'sgparsclkazrcloudsatC1.c1.20170616.000000.nc',
         'sgparsclkazrcloudsatC1.c1.20170618.000000.nc',
         'sgparsclkazrcloudsatC1.c1.20170627.000000.nc',
         'sgparsclkazrcloudsatC1.c1.20170630.000000.nc',
         'sgparsclkazrcloudsatC1.c1.20170510.000000.nc',
         'sgparsclkazrcloudsatC1.c1.20170511.000000.nc',
         'sgparsclkazrcloudsatC1.c1.20170512.000000.nc',
         'sgparsclkazrcloudsatC1.c1.20170516.000000.nc',
         'sgparsclkazrcloudsatC1.c1.20170517.000000.nc',
         'sgparsclkazrcloudsatC1.c1.20170518.000000.nc',
         'sgparsclkazrcloudsatC1.c1.20170519.000000.nc',
         'sgparsclkazrcloudsatC1.c1.20170520.000000.nc',
         'sgparsclkazrcloudsatC1.c1.20170521.000000.nc',
         'sgparsclkazrcloudsatC1.c1.20170528.000000.nc',
         ]

from kazrRet import *
#'sgparsclkazrcloudsatC1.c1.20170707.000000.nc':[8000,8300],

#'sgparsclkazrcloudsatC1.c1.20170708.000000.nc':[5000,10000],
#'sgparsclkazrcloudsatC1.c1.20170607.000000.nc':[8000,16000],
d={'sgparsclkazrcloudsatC1.c1.20170708.000000.nc':[28000,32000],
   'sgparsclkazrcloudsatC1.c1.20170614.000000.nc':[2000,8000],
   'sgparsclkazrcloudsatC1.c1.20170616.000000.nc':[0,8000],
   'sgparsclkazrcloudsatC1.c1.20170627.000000.nc':[24000,32000],
   'sgparsclkazrcloudsatC1.c1.20170630.000000.nc':[20500,24500],
   'sgparsclkazrcloudsatC1.c1.20170517.000000.nc':[0,10000]}
import numpy as np
ic=-1
for k in d.keys():
    fname=str(k)
    ic+=1
    if ic!=2:
        continue
    fh=Dataset('data/'+fname)
    z=fh['reflectivity_best_estimate'][:]
    t=fh['time'][:]
    h=fh['height'][:]
    t1,t2=d[fname]    
    nx=z.shape[0]

    plt.figure(figsize=(12,6))
    a=np.nonzero((t-t1)*(t-t2)<0)
    c=plt.pcolormesh(t[a[0]],h/1e3,z[a[0],:].T,\
                     cmap='jet',vmin=-20,vmax=15)
    plt.contour(t[a[0]],h/1e3,z[a[0],:].T,\
                levels=[8],colors=['black'])
    plt.ylim(2,15)
    plt.colorbar(c)
    print(fname)
    #plt.show()
    if ic==2:
        break
   
    
from scipy.ndimage import gaussian_filter
dn2d=np.random.randn(a[0].shape[0],400)*2

dn2d=gaussian_filter(dn2d,sigma=3)*2
iwc2d=np.zeros((a[0].shape[0],400))
import lidarSim as lidSim
dn1d=h[100:500].data/1e3*0.1

dmCoeffs=np.polyfit(zKuG,np.log(dmG),1)

zCoeffs=np.polyfit(np.log10(gwc),zKuG,1)
dm_ice=np.zeros((a[0].shape[0],400))
nfreq=3
kexttot_2d=np.zeros((a[0].shape[0],400,nfreq),float)
kscatot_2d=np.zeros((a[0].shape[0],400,nfreq),float)
asymtot_2d=np.zeros((a[0].shape[0],400,nfreq),float)
zKu_2d=np.zeros((a[0].shape[0],400),float)
#zSSRG_2d=np.zeros((a[0].shape[0],400,nfreq),float)
for ik,iprof in enumerate(a[0]):
    z1d=z[iprof,100:500]
    iwc2d
    for k,zka1 in enumerate(z1d):
        if zka1>-12:
            dn1=dn1d[k]+dn2d[ik,k]
            if zka1-10*dn1<-12:
                dn1=(zka1+12)/10.
            dn2d[ik,k]=dn1
            ifind = lidSim.bisection2(zKaS,zka1-10*dn1)
            iwc2d[ik,k]=swc[ifind]*10**dn1
            zKu=zKuS[ifind]
            zKu_2d[ik,k]=zKu+10*dn1
            dm_ice[ik,k]=np.exp(dmCoeffs[0]*(zKu)+dmCoeffs[1])
            ibin2=lidSim.bisection2(dmST[-1,:],dm_ice[ik,k])
            iwc2=iwcST[-1,ibin2]*10**dn1
            #zSSRG[:3,k]=zST[:3,ibin2]+10*dn1
            kexttot_2d[ik,k,:]=kextST[-3:,ibin2]*10**dn1
            kscatot_2d[ik,k,:]=kscaST[-3:,ibin2]*10**dn1
            asymtot_2d[ik,k,:]=gST[-3:,ibin2]

import matplotlib
plt.figure()
plt.pcolormesh(iwc2d.T,cmap='jet',norm=matplotlib.colors.LogNorm())
plt.colorbar()

import xarray as xr
iwc2dx=xr.DataArray(iwc2,dims=['nt','nz'])
dn2dx=xr.DataArray(dn2d,dims=['nt','nz'])
zkax=xr.DataArray(z[a[0],100:500],dims=['nt','nz'])
hx=xr.DataArray(h[100:500],dims=['nz'])
dkazr=xr.Dataset({"iwc":iwc2dx,"dn":dn2dx,"zka":zkax,"height":h})
dkazr.to_netcdf("kazret."+fname[-18:])
