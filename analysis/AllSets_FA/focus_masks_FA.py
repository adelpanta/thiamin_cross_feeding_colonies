import tifffile
import aicsimageio as aics
from skimage import morphology, filters
from scipy import ndimage as ndi
import sys
import os
import pathlib
import re
import numpy as np
import glob
import h5py
from concurrent.futures import ProcessPoolExecutor

# Global variables
media = ["FAT", "FA"]
media_titles = ["FAT", "FA"]
nmedia = 2
nreplicates=3
chan_r = 0
chan_g = 1
thresholds_to_try = np.linspace(0,1,101)
chosen_thr = 80
all_partners_in_collection = ['At','AeCa','BrVe','CiFr','G3','GN1','GN3','GN4','Ct','G1','G4','G5','G6','G8','AlFa','KlPn','PsFu','PsAe','PsFl','PsPs','PrRe','S1','S2','S5','SC1','SC12','SC14','SC18']

# Setup folders 
proj_dir = pathlib.Path(os.getcwd()).parent.parent
raw_data_dir = proj_dir / 'raw_data' 
proc_data_dir = proj_dir / 'processed_data' 

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

    # Paths for Ilastik probabilities
    all_proc_folders = glob.glob(str(proc_data_dir / '*SET*'))
    all_proc_folders.sort()
    fa_prob_paths = []
    fat_prob_paths = []
    # Loop through files in all_closeup_folders and add to paths_closeup_fa and paths_closeup_fat
    for folder in all_proc_folders:
        # Loop through files in folder and add to paths_closeup_fa and paths_closeup_fat
        for file in glob.glob(str(pathlib.Path(folder) / '*')):
            if any(partner in file for partner in all_partners_in_collection) and "Probabilities" in file:
                if "fa" in file and "fat" not in file:
                    fa_prob_paths.append(file)
                elif "fat" in file:
                    fat_prob_paths.append(file)

    if len(fa_prob_paths) != len(fat_prob_paths):
        print(f"Warning: Found {len(fa_prob_paths)} FA closeup images and {len(fat_prob_paths)} FAT closeup images for the 28 partners in the collection. The numbers do not match.")
    print(f"Found {len(fa_prob_paths)} FA closeup images and {len(fat_prob_paths)} FAT closeup images for the 28 partners in the collection.")

    def sort_key(path,mode="FAT"):
        if mode == "FAT":
            m = re.search(r'/([^/]+)_.*?_fat(\d*)_Probabilities\.h5$', path)
        else:
            m = re.search(r'/([^/]+)_.*?_fa(\d*)_Probabilities\.h5$', path)
        name = pathlib.Path(path).name
        partner = name.split("_")[0]  # e.g. "AlFa"
        replicate = int(m.group(2)) if m.group(2) else 1  # fa -> 1, fa2 -> 2, fa3 -> 3
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
        if "fat" in path:
            mode = "FAT"
        elif "fa" in path:
            mode = "FA"
        else:
            raise ValueError(f"Path {path} does not contain 'fa' or 'fat'.")
        replicate = re.search(r'_(fa|fat)(\d*)_Probabilities\.h5$', path).group(2)
        if not replicate:
            replicate = "1"
            if mode == "FAT":
                new_path = path.replace("closeup_norm","focus_mask").replace("fat","FAT").replace(" Probabilities.h5",f"{replicate}.tiff")
            elif mode == "FA":
                new_path = path.replace("closeup_norm","focus_mask").replace("fa","FA").replace("Probabilities.h5",f"{replicate}.tiff")
        else:
            if mode == "FAT":
                new_path = path.replace("closeup_norm","focus_mask").replace(f"fat{replicate}","_FAT").replace("Probabilities.h5",f"{replicate}.tiff")
            elif mode == "FA":
                new_path = path.replace("closeup_norm","focus_mask").replace(f"fa{replicate}","_FA").replace("Probabilities.h5",f"{replicate}.tiff")
        return new_path

    save_paths_fa = [prob_path_to_focus_mask_path(p) for p in fa_prob_paths]
    save_paths_fat = [prob_path_to_focus_mask_path(p) for p in fat_prob_paths]
    if len(save_paths_fa) != len(save_paths_fat):
        print(f"Warning: Defined {len(save_paths_fa)} FA focus mask paths and {len(save_paths_fat)} FAT focus mask paths. The numbers do not match.")
    else:
        print(f"Defined {len(save_paths_fa)} FA focus mask paths and {len(save_paths_fat)} FAT focus mask paths to save the focus masks.")

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