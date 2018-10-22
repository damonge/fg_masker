import healpy as hp
import numpy as np
import os
import matplotlib.pyplot as plt

#Download Planck 353 maps if not present 
if not os.path.isfile("HFI_SkyMap_353_2048_R2.02_full.fits") :
    import wget
    print("Downloading Planck 353")
    wget.download("https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/HFI_SkyMap_353_2048_R2.02_full.fits","HFI_SkyMap_353_2048_R2.02_full.fits")
    print("\nDone")

#Read 353 maps
print("Reading 353 maps")
mp_i,mp_q,mp_u=hp.read_map("HFI_SkyMap_353_2048_R2.02_full.fits",field=[0,1,2],verbose=False)

#Smooth to signal-dominated scale determined by https://arxiv.org/pdf/1608.02841.pdf (below Eq. 15)
print("Smoothing")
theta_fwhm=np.pi/69.
mp_i,mp_q,mp_u=hp.smoothing([mp_i,mp_q,mp_u],fwhm=theta_fwhm,verbose=False)
mp_i=hp.ud_grade(mp_i,nside_out=512)
mp_q=hp.ud_grade(mp_q,nside_out=512)
mp_u=hp.ud_grade(mp_u,nside_out=512)

#Save to file
print("Saving")
hp.write_map("HFI_353_IQU_smooth.fits",[mp_i,mp_q,mp_u])
