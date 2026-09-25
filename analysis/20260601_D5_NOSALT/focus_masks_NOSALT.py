import tifffile
import aicsimageio as aics
from skimage.morphology import binary_dilation, binary_erosion
from skimage import filters, morphology
from scipy import ndimage as ndi
import sys
import os
import pathlib
import glob
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import h5py
import re

# Global variables
media = ["FAT_ns", "FA_ns"]
media_titles = ["FATns", "FAns"]
nmedia = 2
nreplicates=3
chan_r = 0
chan_g = 1
thresholds_to_try = np.linspace(0,1,101)
chosen_thr = 80 #67
partners = ['G3','GN4','Ct','AlFa','G4','PsAe','PsFu','S1']

# Setup folders 
proj_dir = pathlib.Path(os.getcwd()).parent.parent
raw_data_dir = proj_dir / 'raw_data' / '20260601_D5_NOSALT'
proc_data_dir = proj_dir / 'processed_data' / '20260601_D5_NOSALT'

if not os.path.exists(proc_data_dir):
   os.makedirs(proc_data_dir)

# Compute parameters
chunk_size = 6  # Adjust based on your system's capabilities

# Define functions for morphological operations
def square_closing(mask, size):
    selem = morphology.square(size)
    return morphology.binary_closing(mask, selem)

def rectangular_closing(mask,rect_height=5, rect_width=700):
    return morphology.binary_closing(mask, morphology.rectangle(rect_height, rect_width))

def rectangular_dilation(mask, rect_height=10, rect_width=800):
    return morphology.binary_dilation(mask, morphology.rectangle(rect_height, rect_width))

def circular_closing(mask, radius=30):
    return morphology.binary_closing(mask, morphology.disk(radius)) 

if __name__ == "__main__":

    # Paths to Ilastik probabilities
    fat_prob_paths = [str(proc_data_dir / f"{partner}_closeup_norm_FATns_{i}_Probabilities.h5") for partner in partners for i in range(1,nreplicates+1)]
    fa_prob_paths = [str(proc_data_dir / f"{partner}_closeup_norm_FAns_{i}_Probabilities.h5") for partner in partners for i in range(1,nreplicates+1)]
    fat_prob_paths.sort()
    fa_prob_paths.sort()
    if len(fat_prob_paths) != len(fa_prob_paths):
        raise ValueError("Number of FA and FAT probability files do not match.")
    else:
        print(f"Found {len(fat_prob_paths)} FAT probability files and {len(fa_prob_paths)} FA probability files.")

    if len(fa_prob_paths) != len(fat_prob_paths):
        print(f"Warning: Found {len(fa_prob_paths)} FA closeup images and {len(fat_prob_paths)} FAT closeup images for the 28 partners in the collection. The numbers do not match.")
    print(f"Found {len(fa_prob_paths)} FA closeup images and {len(fat_prob_paths)} FAT closeup images for the 28 partners in the collection.")

    def sort_key(path,mode="FAT"):
        if mode == "FAT":
            m = re.search(r'/([^/]+)_.*?_FATns_(\d*)_Probabilities\.h5$', path)
        else:
            m = re.search(r'/([^/]+)_.*?_FAns_(\d*)_Probabilities\.h5$', path)
        name = pathlib.Path(path).name
        partner = name.split("_")[0]  # e.g. "AlFa"
        replicate = int(m.group(2)) 
        return (partner, replicate)

    fat_prob_paths.sort(key=lambda x: sort_key(x, mode="FAT"))
    fa_prob_paths.sort(key=lambda x: sort_key(x, mode="FA"))
    print(f"Sorted the FAT and FA probability paths.")
    
    # Load all probabilities
    fat_prob = [h5py.File(p, 'r')['exported_data'][:,:,0] for p in fat_prob_paths]
    fa_prob = [h5py.File(p, 'r')['exported_data'][:,:,0] for p in fa_prob_paths]
    print(f"Loaded probabilities for {len(fat_prob)} FAT replicates and {len(fa_prob)} FA replicates.")

    # Calculate masks
    mask_fat = [(prob > thresholds_to_try[chosen_thr]).astype(int) for prob in fat_prob]
    mask_fa = [(prob > thresholds_to_try[chosen_thr]).astype(int) for prob in fa_prob]
    print(f"Calculated masks for {len(mask_fat)} FAT replicates and {len(mask_fa)} FA replicates.")

    # Define paths to save focus masks
    def prob_path_to_focus_mask_path(path):
        new_path = path.replace("closeup_norm","focus_mask").replace("_Probabilities.h5",".tiff")
        return new_path

    save_paths_fa = [prob_path_to_focus_mask_path(p) for p in fa_prob_paths]
    save_paths_fat = [prob_path_to_focus_mask_path(p) for p in fat_prob_paths]
    if len(save_paths_fa) != len(save_paths_fat):
        print(f"Warning: Defined {len(save_paths_fa)} FA focus mask paths and {len(save_paths_fat)} FAT focus mask paths. The numbers do not match.")
    else:
        print(f"Defined {len(save_paths_fa)} FA focus mask paths and {len(save_paths_fat)} FAT focus mask paths to save the focus masks.")

    print(save_paths_fa)
    print("-----------------------------------------------------------")
    print(save_paths_fat)
    # Morphological operations on FA masks in parallel
    print("-----------------------------------------------------------")
    print("Starting morphological operations on FA masks in parallel...")
    print("-----------------------------------------------------------")
    for i in range(0, len(mask_fa), chunk_size):
        print(f"Processing chunk {i//chunk_size + 1} of {len(mask_fa)//chunk_size + (1 if len(mask_fa) % chunk_size else 0)} for square closing...")
        chunk = mask_fa[i:i + chunk_size]
        save_chunk = save_paths_fa[i:i + chunk_size]
        print("Calculating square closing for the current chunk...")
        with ProcessPoolExecutor() as executor:
            closed_mask_fa_square = list(executor.map(square_closing, chunk, [50]*len(chunk)))
        print("Calculating rectangular dilation for the current chunk...")
        with ProcessPoolExecutor() as executor:
            dilated_mask_rect_fa = list(executor.map(rectangular_dilation, closed_mask_fa_square))
        print("Removing small holes for the current chunk...")
        mask_fa_noholes = [morphology.remove_small_holes(mask, area_threshold=100000) for mask in dilated_mask_rect_fa]
        print("Applying small circular dilation for the current chunk...")
        with ProcessPoolExecutor() as executor:
            focus_mask_fa = list(executor.map(circular_closing, mask_fa_noholes))
        print("Saving focus masks for the current chunk...")
        for mask, save_path in zip(focus_mask_fa, save_chunk):
            tifffile.imwrite(save_path, mask)

    # Morphological operations on FAT masks in parallel
    print("-----------------------------------------------------------")
    print("Starting morphological operations on FAT masks in parallel...")
    print("-----------------------------------------------------------")
    for i in range(0, len(mask_fat), chunk_size):
        print(f"Processing chunk {i//chunk_size + 1} of {len(mask_fat)//chunk_size + (1 if len(mask_fat) % chunk_size else 0)} for square closing...")
        chunk = mask_fat[i:i + chunk_size]
        save_chunk = save_paths_fat[i:i + chunk_size]
        print("Calculating square closing for the current chunk...")
        with ProcessPoolExecutor() as executor:
            closed_mask_fat_square = list(executor.map(square_closing, chunk, [50]*len(chunk)))
        print("Calculating rectangular dilation for the current chunk...")
        with ProcessPoolExecutor() as executor:
            dilated_mask_rect_fat = list(executor.map(rectangular_dilation, closed_mask_fat_square))
        print("Removing small holes for the current chunk...")
        mask_fat_noholes = [morphology.remove_small_holes(mask, area_threshold=100000) for mask in dilated_mask_rect_fat]
        print("Applying small circular dilation for the current chunk...")
        with ProcessPoolExecutor() as executor:
            focus_mask_fat = list(executor.map(circular_closing, mask_fat_noholes))
        print("Saving focus masks for the current chunk...")
        for mask, save_path in zip(focus_mask_fat, save_chunk):
            tifffile.imwrite(save_path, mask)