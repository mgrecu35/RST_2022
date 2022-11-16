#import xarray as xr
import numpy as np
from netCDF4 import Dataset
import lidarSim as lidar

#fh=Dataset("lidarInput_kazr.20170616.000000.nc",'a')
fh=Dataset("lidarInput_kazr.20170627.000000.nc","a")
fh=Dataset("lidarInput_kazr.20170517.000000.nc","a")
fh=Dataset("lidarInput_kazr.20170708.000000.nc","a")
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

import numpy as np
nx,nz=pres.shape
kext=fh["kext"][:]
kscat=fh["kscat"][:]
gtot=fh["gtot"][:]
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
        
for i in range(nx):
    for ifreq in range(4):
        
#fh.createVariable("pnorm3D","f8",("dim_0","dim_1","dim_2"))
#fh.createVariable("betatot3D","f8",("dim_0","dim_1","dim_2"))
#fh.createVariable("tau_mol3D","f8",("dim_0","dim_1","dim_2"))
#fh.createVariable("beta_mol3D","f8",("dim_0","dim_1","dim_2"))
#fh.createDimension("dim_part",4)
#fh.createVariable("alpha_3D","f8",("dim_part","dim_0","dim_1","dim_2"))
#fh.createVariable("extinct3D","f8",("dim_0","dim_1","dim_2"))
#fh["pnorm3D"][:]=pnorm3D
#fh["betatot3D"][:]=betatot3D
#fh["extinct3D"][:]=extinct3D
#fh["tau_mol3D"][:]=tau_mol3D
#fh["beta_mol3D"][:]=beta_mol3D
#fh["alpha_3D"][:]=alpha_3D
#zku=fh["zKu"][:]
#absair,abswv = gasabsr98(f,tk,rhowv,pa,ireturn)
 
fh.close()
import matplotlib
import matplotlib.pyplot as plt
plt.pcolormesh(pnorm.T,cmap='jet',norm=matplotlib.colors.LogNorm(vmin=1e-8))
plt.contour(zku.T,levels=[8,12],colors=['black','black'])
