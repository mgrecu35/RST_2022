
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
from procSound import *
presS,hgtS,tempS,qvS=getEnv()
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
from scipy.ndimage import gaussian_filter
import lidarSim as lidSim
dmCoeffs=np.polyfit(zKuG,np.log(dmG),1)
zCoeffs=np.polyfit(np.log10(gwc),zKuG,1)
import matplotlib
import xarray as xr

for k in d.keys():
    fname=str(k)
    ic+=1
    #if ic!=0:
    #    continue
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
  
  
   

    
    dn2d=np.random.randn(a[0].shape[0],400)*2
    dn2d=gaussian_filter(dn2d,sigma=3)*2
    iwc2d=np.zeros((a[0].shape[0],400))

    dn1d=h[100:500].data/1e3*0.1


    dm_ice=np.zeros((a[0].shape[0],400))
    nfreq=4
    kexttot_2d=np.zeros((a[0].shape[0],400,nfreq),float)
    kscatot_2d=np.zeros((a[0].shape[0],400,nfreq),float)
    asymtot_2d=np.zeros((a[0].shape[0],400,nfreq),float)
    zKu_2d=np.zeros((a[0].shape[0],400),float)
    
    for ik,iprof in enumerate(a[0][:]):
        z1d=z[iprof,100:500]
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
                kexttot_2d[ik,k,:]=kextST[-4:,ibin2]*10**dn1
                kscatot_2d[ik,k,:]=kscaST[-4:,ibin2]*10**dn1
                asymtot_2d[ik,k,:]=gST[-4:,ibin2]


    plt.figure()
    plt.pcolormesh(iwc2d.T,cmap='jet',norm=matplotlib.colors.LogNorm())
    plt.colorbar()
    

    iwc2dx=xr.DataArray(iwc2d,dims=['nt','nz'])
    dn2dx=xr.DataArray(dn2d,dims=['nt','nz'])
    zkax=xr.DataArray(z[a[0],100:500],dims=['nt','nz'])
    hx=xr.DataArray(h[100:500],dims=['nz'])
    zKu2d=xr.DataArray(zKu_2d,dims=['nt','nz'])
    dm_icex=xr.DataArray(dm_ice,dims=['nt','nz'])
    kexttot_2dx=xr.DataArray(kexttot_2d,dims=['nt','nz','nfreq'])
    kscatot_2dx=xr.DataArray(kscatot_2d,dims=['nt','nz','nfreq'])
    asymtot_2dx=xr.DataArray(asymtot_2d,dims=['nt','nz','nfreq'])
    
    dkazr=xr.Dataset({"iwc":iwc2dx,"dn":dn2dx,"zka":zkax,"height":h,\
                      "dm_ice":dm_icex,"kext":kexttot_2dx,"kscat":kscatot_2dx,"asymtot":asymtot_2dx})
    dkazr.to_netcdf("kazret."+fname[-18:])


    dz=(h[1]-h[0])*3
    nz=400
    nz3=134
    nx=a[0].shape[0]
    temp=np.zeros((nz3,nx),float)
    pres=np.zeros((nz3,nx),float)
    presf=np.zeros((nz3+1,nx),float)
    q_lsliq=np.zeros((nz3,nx),float)
    q_lsice=np.zeros((nz3,nx),float)
    q_cvliq=np.zeros((nz3,nx),float)
    q_cvice=np.zeros((nz3,nx),float)
    ls_radliq=np.zeros((nz3,nx),float)
    ls_radice=np.zeros((nz3,nx),float)
    cv_radice=np.zeros((nz3,nx),float)
    cv_radliq=np.zeros((nz3,nx),float)
    iwc_lidar=np.zeros((nz3,nx),float)
    
    #presS,hgtS,tempS,qvS
    rho=presS/287.15/(tempS+273.15)*1e2
    z0=h[100]/1e3
    print(dz)
    for i in range(nx):
        temp1=np.interp(h[100:500:3]/1e3,hgtS/1e3,tempS+273.15)
        pres1=np.interp(h[100:500:3]/1e3,hgtS/1e3,presS*1e2)
        presf1=np.interp(z0-dz/2e3+np.arange(nz3+1)*dz/1e3,hgtS/1e3,presS*1e2)
        rho1=np.interp(h[100:500:3]/1e3,hgtS/1e3,rho)
        q_lsice1=iwc2d[i,::3].data/rho1*1e-3
        iwc_lidar[:,i]=iwc2d[i,::3].data
        q_lsice[:,i]=(q_lsice1)
        q_lsliq[:,i]=(np.zeros((nz3),float))
        ls_radice[:,i]=(dm_ice[i,::3]/2*1e-3)
        ls_radliq[:,i]=(np.zeros((nz3),float))
        q_cvice[:,i]=(np.zeros((nz3),float))
        cv_radice[:,i]=(np.zeros((nz3),float))
        q_cvliq[:,i]=(np.zeros((nz3),float))
        cv_radliq[:,i]=(np.zeros((nz3),float))
        temp[:,i]=(temp1)
        pres[:,i]=(pres1)
        presf[:,i]=(presf1)


    temp=np.array(temp)

    q_lsliqX=xr.DataArray(q_lsliq.T,dims=['nt','nz'])
    q_lsiceX=xr.DataArray(q_lsice.T,dims=['nt','nz'])
    q_cvliqX=xr.DataArray(q_cvliq.T,dims=['nt','nz'])
    iwc_lidarX=xr.DataArray(iwc_lidar.T,dims=['nt','nz'])
    q_cviceX=xr.DataArray(q_cvice.T,dims=['nt','nz'])
    ls_radliqX=xr.DataArray(ls_radliq.T,dims=['nt','nz'])
    ls_radiceX=xr.DataArray(ls_radice.T,dims=['nt','nz'])
    cv_radliqX=xr.DataArray(cv_radliq.T,dims=['nt','nz'])
    cv_radiceX=xr.DataArray(cv_radice.T,dims=['nt','nz'])
    tempX=xr.DataArray(temp.T,dims=['nt','nz'])
    presX=xr.DataArray(pres.T,dims=['nt','nz'])
    zKuX=xr.DataArray(zKu_2d[:,::3],dims=['nt','nz'])
    presfX=xr.DataArray(presf.T,dims=['nt','nz1'])
    kexttotX=xr.DataArray(kexttot_2d[:,::3],dims=['nt','nz','nfreq'])
    kscatotX=xr.DataArray(kscatot_2d[:,::3],dims=['nt','nz','nfreq'])
    gtotX=xr.DataArray(asymtot_2d[:,::3],dims=['nt','nz','nfreq'])
    dn2dX=xr.DataArray(dn2d[:,::3],dims=['nt','nz'])
    rhoX=xr.DataArray(rho1,dims=['nz'])
    tempSX=xr.DataArray(tempS+273.15,dims=['nzs'])
    hgtSX=xr.DataArray(hgtS,dims=['nzs'])
    presSX=xr.DataArray(presS,dims=['nzs'])
    qvSX=xr.DataArray(qvS,dims=['nzs'])
    
    dd={"q_lsliq":q_lsliqX,"q_lsice":q_lsiceX,\
        "q_cvliq":q_cvliqX,"q_cvice":q_cviceX,\
        "ls_radliq":ls_radliqX,"ls_radice":ls_radiceX,\
        "cv_radliq":cv_radliqX,"cv_radice":cv_radiceX,\
        "temp":tempX,"pres":presX,"presfX":presfX,"zKu":zKuX,\
        "kext":kexttotX,"kscat":kscatotX,"gtot":gtotX,\
        "dn2d":dn2dX,"rho":rhoX, "z":xr.DataArray(h[100:500:3]),\
        "presS":presSX,"tempS":tempSX,"qvS":qvSX,"presS":presSX,"hgtS":hgtSX,"iwc":iwc_lidarX}
    dX=xr.Dataset(dd)

    comp = dict(zlib=True, complevel=5)
    encoding = {var: comp for var in dX.data_vars}
    dX.to_netcdf("lidarInput_kazr.%s.nc"%fname[-18:-3], encoding=encoding)
    #if ic==0:
    #    break

#npart=4
#nrefl=4

