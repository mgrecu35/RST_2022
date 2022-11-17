#import xarray as xr
import numpy as np
from netCDF4 import Dataset
import lidarSim as lidar
import matplotlib
import matplotlib.pyplot as plt

import numpy as np
import glob

fnames=glob.glob("inputData/*nc")
def procFile(fname):
    fh=Dataset(fname,'a')
    #fh=Dataset("inputData/lidarInput_kazr.20170627.000000.nc","a")
    #fh=Dataset("inputData/lidarInput_kazr.20170517.000000.nc","a")
    #fh=Dataset("inputData/lidarInput_kazr.20170708.000000.nc","a")
    q_lsliq=fh["q_lsliq"][:]
    q_lsice=fh["q_lsice"][:]
    q_cvliq=fh["q_cvliq"][:]
    q_cvice=fh["q_cvice"][:]
    ls_radliq=fh["ls_radliq"][:]
    ls_radice=fh["ls_radice"][:]
    cv_radliq=fh["cv_radliq"][:]
    cv_radice=fh["cv_radice"][:]
    temp=fh["temp"][:]
    pres=fh["pres"][:]
    presf=fh["presfX"][:]
    
    zku=fh["zKu"][:]
    hgtS=fh["hgtS"][:]
    tempS=fh["tempS"][:]
    presS=fh["presS"][:]
    qvS=fh["qvS"][:]
    npart=4
    nrefl=2
    ice_type=1
    undef=0.0

    nx,nz=pres.shape
    kext=fh["kext"][:]
    kscat=fh["kscat"][:]
    gtot=fh["gtot"][:]
    iwc=fh["iwc"][:]
    z=fh["z"][:]
    #pnorm3D=np.zeros((nz,ny,nx),float)
    #betatot3D=np.zeros((nz,ny,nx),float)
    #extinct3D=np.zeros((nz,ny,nx),float)
    #beta_mol3D=np.zeros((nz,ny,nx),float)
    #tau_mol3D=np.zeros((nz,ny,nx),float)
    #alpha_3D=np.zeros((4,nz,ny,nx),float)
    
    pmol,pnorm,pnorm_perp_tot,\
        tautot,betatot_liq,\
        betatot_ice,\
        betatot,refl, \
        zheight,\
        beta_mol, tau_mol,\
        alpha= lidar.lidar_simulator(npart,nrefl,undef,\
                                     pres[:,:],presf[:,:],\
                                     temp[:,:],
                                     q_lsliq[:,:],q_lsice[:,:],\
                                     q_cvliq[:,:],\
                                     q_cvice[:,:],\
                                     ls_radliq[:,:],\
                                     ls_radice[:,:],\
                                    cv_radliq[:,:],cv_radice[:,:],\
                                     ice_type)

    hL=350+np.arange(11)*250
    hgtL=np.concatenate([hL,z-45,np.array([z[-1]+45])])
    hgtm=0.5*(hgtL[1:]+hgtL[:-1])
    pressm=np.interp(hgtm,hgtS,presS)*1e2
    tempm=np.interp(hgtm,hgtS,tempS)
    tempL=np.interp(hgtL,hgtS,tempS)
    qvm=np.interp(hgtm,hgtS,qvS)
    rho1d=pressm/287/tempm
    freqs=['94','183.31', '325.15','660.00']
    ireturn=1
    nz=hgtm.shape[0]
    kexttot_atm=np.zeros((nz,4),float)
    for ik,freq in enumerate(freqs):
        for i in range(nz):
            absair,abswv = lidar.gasabsr98(float(freq),tempm[i],rho1d[i]*qvm[i]*1e-3,pressm[i],ireturn)
            kexttot_atm[i,ik]=(absair+abswv)

    umu=np.cos(53/180*np.pi)
    emis=0.9
    ebar=0.9
    fisot=2.7
    tbL=np.zeros((nx,4))
    tbRadL=np.zeros((nx,4))
    lwp1=np.zeros((nx),float)
    lwp2=np.zeros((nx),float)
    izku=np.zeros((nx),float)
    for i in range(nx):
        for ifreq in range(2,3):
            kext1d=kexttot_atm[:,ifreq].copy()
            kext1d[11:]+=kext[i,:,ifreq]
            salb1d=np.zeros((nz),float)
            salb1d[11:]+=kscat[i,:,ifreq]
            asym1d=np.zeros((nz),float)
            asym1d[11:]+=gtot[i,:,ifreq]
            salb1d/=kext1d
            tbout = lidar.radtran(umu,tempL[0],tempL,hgtL/1e3,kext1d,salb1d,asym1d,fisot,emis,ebar)
            tbL[i,ifreq]=tbout
            
            kext1d=kexttot_atm[:,ifreq].copy()
            a=np.nonzero(zku[i,:]>8)
            salb1d=np.zeros((nz),float)

            asym1d=np.zeros((nz),float)
            if len(a[0])>0:
                kext1d[11:][a]+=kext[i,:,ifreq][a]
                salb1d[11:][a]+=kscat[i,:,ifreq][a]
                asym1d[11:][a]+=gtot[i,:,ifreq][a]
            salb1d/=kext1d
            tbout = lidar.radtran(umu,tempL[0],tempL,hgtL/1e3,kext1d,salb1d,asym1d,fisot,emis,ebar)
            tbRadL[i,ifreq]=tbout
            lwp2[i]=iwc[i,:][a].sum()*(z[1]-z[0])
        a=np.nonzero(zku[i,:]>0)
        izku[i]=len(a[0])
        lwp1[i]=iwc[i,:].sum()*(z[1]-z[0])
        
            
    fh.close()
    a1=np.nonzero(izku>10)
    plt.scatter(lwp1[a1],tbL[a1[0],ifreq],s=1)
    return len(a1[0]),tbL,tbRadL,lwp1,lwp2,zku,pnorm

legL=[]
nt=0
lwpS1,lwpS2=[],[]
for fname in sorted(fnames)[1:2]:
    nt1,tbL,tbRadL,lwp1,lwp2,zku,pnorm=procFile(fname)
    legL.append(fname[-18:-10])
    print((tbRadL[:,2]-tbL[:,2]).max(),fname)
    lwpS1.append(lwp1.sum())
    lwpS2.append(lwp2.sum())
    nt+=nt1
plt.legend(legL)
#plt.pcolormesh(pnorm.T,cmap='jet',norm=matplotlib.colors.LogNorm(vmin=1e-8))
#plt.contour(zku.T,levels=[8,12],colors=['black','black'])
