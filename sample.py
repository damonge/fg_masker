import numpy as np
from fg_masks import mask_353
from pixell import enmap

box=np.deg2rad([[-1.,13.47],[1.,17.47]])
map_en=enmap.read_map("/project/projectdirs/act/data/xlink/deep56_forecast_180223_master_apo_w0.fits",box=box)
msk=mask_353(map_en.shape,map_en.wcs,area_cutoff=0.5)
