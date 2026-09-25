import tifffile
import aicsimageio as aics
from skimage.morphology import binary_dilation, binary_erosion
from skimage import filters, morphology
from scipy import ndimage as ndi
import sys
import os
import pathlib
import re
import numpy as np

# Global variables
media = ["GluT", "Glu"]
media_titles = ["GluT", "Glu"]
nmedia = 2
nreplicates=3
chan_r = 0
chan_g = 1

# Setup folders 
proj_dir = pathlib.Path(os.getcwd()).parent.parent
raw_data_dir = proj_dir / 'raw_data' / '20260615_D5_GLUCOSE'
proc_data_dir = proj_dir / 'processed_data' / '20260615_D5_GLUCOSE'

if not os.path.exists(proc_data_dir):
   os.makedirs(proc_data_dir)

   
# Import functions
sys.path.insert(0, str(proj_dir / 'code'))
from _image_analysis_functions import *   

if __name__ == "__main__":

    # Image paths for closeups
    paths_closeup = list((raw_data_dir / "CloseUps").glob("*.czi"))
    paths_closeup_fat = [f for f in paths_closeup if re.fullmatch(r".*T\d+\.czi", f.name) and "OaAlone" not in f.name]
    paths_closeup_fa = [f for f in paths_closeup if f not in paths_closeup_fat and "OaAlone" not in f.name]
    if len(paths_closeup_fa) != len(paths_closeup_fat):
        raise ValueError("Number of FA and FAT closeup images do not match.")
    paths_closeup_fa.sort()
    paths_closeup_fat.sort()

    # Calculate and save each mask immediately
    for path_fa, path_fat in zip(paths_closeup_fa, paths_closeup_fat):

        # FA image
        print("Processing FA image: %s" % path_fa.name)
        img_fa = np.squeeze(aics.AICSImage(path_fa).get_image_data())[chan_r]
        focus_mask_fa = find_focus_mask(
            img_fa,
            d_inner=10,
            min_size=1000,
            closing_radius=100
        )

        savename_fa = (
            path_fa.name[:-5] +
            "_focus_mask_Glu_" +
            path_fa.name[-5:]
        ).replace(".czi", ".tiff")

        tifffile.imwrite(proc_data_dir / savename_fa, focus_mask_fa.astype(np.uint8))

        # FAT image
        print("Processing FAT image: %s" % path_fat.name)
        img_fat = np.squeeze(aics.AICSImage(path_fat).get_image_data())[chan_r]
        focus_mask_fat = find_focus_mask(
            img_fat,
            d_inner=10,
            min_size=1000,
            closing_radius=100
        )

        savename_fat = (
            path_fat.name[:-5] +
            "_focus_mask_GluT_" +
            path_fat.name[-5:]
        ).replace(".czi", ".tiff")

        tifffile.imwrite(proc_data_dir / savename_fat, focus_mask_fat.astype(np.uint8))