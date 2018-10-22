import numpy as np
from pixell import enmap,reproject
import healpy as hp

def mask_353(en_shape,en_wcs,intensity_cutoff=None,area_cutoff=None,
             fname_353="HFI_353_IQU_smooth.fits",return_healpix=False) :
    """
    Returns binary mask based on total polarized intensity at 353 GHz.
    
    :param en_shape: shape of the output enmap
    :param en_wcs: wcs parameters of the output enmap
    :param intensity_cutoff: if not None, the mask will remove all pixels where the polarized intensity is larger than the maximum polarized intensity across the sky multiplied by this number.
    :param area_cutoff: if not None, this number sets the sky fraction to be removed (e.g. if set to 0.5, the mask will remove the most contaminated 50% of the sky).
    :param fname_353: path to FITS file containing the 353 map.
    :param return_healpix: if set to True, will return both enmap mask and original HEALPix mask.
    """
    #Read 353 and compute polarized intensity
    hp_tqu=np.array(hp.read_map(fname_353,field=[0,1,2]))
    hp_p=np.sqrt(hp_tqu[1]**2+hp_tqu[2]**2);

    #Set intensity cutoff (as a fraction of maximum intensity)
    if intensity_cutoff is not None :
        #This may be set directly as input
        p_cut=intensity_cutoff
    elif area_cutoff is not None :
        #It can be set indirectly by requesting a given total sky fraction
        if area_cutoff>=1. :
            #Just pick everything
            p_cut=1.001
        else :
            #Find cutoff value that selects the right sky fraction
            hist,b=np.histogram(np.log10(hp_p),bins=200)
            b=0.5*(b[1:]+b[:-1])
            hist=np.cumsum(hist)/np.sum(hist+0.)
            p_cut=10.**(b[np.where(hist>area_cutoff)[0][0]])/np.amax(hp_p)
    else :
        #If none are selected, we just return 1 everywhere
        p_cut=1.0001

    #Create healpix mask
    hp_mask=np.ones_like(hp_p); hp_mask[hp_p>=p_cut*np.amax(hp_p)]=0.

    #Transform to enmap
    en_mask=reproject.enmap_from_healpix(hp_mask,en_shape,en_wcs,ncomp=1,rot='gal,equ')
    #Transformation will yield floating point numbers other than 1 and 0. Binarize.
    en_mask[en_mask>0.9]=1; en_mask[en_mask<=0.9]=0;

    if return_healpix :
        return en_mask,hp_mask
    else :
        return en_mask

box=np.deg2rad([[-1.,13.47],[1.,17.47]])
map_en=enmap.read_map("/project/projectdirs/act/data/xlink/deep56_forecast_180223_master_apo_w0.fits",box=box)
mask_353(map_en.shape,map_en.wcs,area_cutoff=0.5)
