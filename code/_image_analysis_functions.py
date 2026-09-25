import numpy as np
from skimage import morphology
from skimage import filters
from scipy import ndimage as ndi

def normalize_image_global(img):
    img_norm = (img - np.min(img))/(np.max(img)-np.min(img))
    return img_norm

def remove_shadow_from_mask(mask):
    mask_new = np.zeros(mask.shape)
    for i in range(mask.shape[1]):
        if np.max(mask[:,i]) > 0:
            k = np.where(mask[:,i] > 0)[0][0]
            mask_new[k:,i] = 1
    return mask_new

def find_focus_mask(img, d_inner=10, min_size=1000, closing_radius=100):
    """
    This function takes as input a fluorescent image of a colony edge and returns a mask of the in-focus region.
    The mask is obtained by segmenting the fluorescent image and then applying closing and filling operations to obtain a binary mask of the colony accounting for the presence of the non-fluorescent species.
    Input:
    img: 2D numpy array of the normalized image of the colony edge
    d_inner: distance in pixels from the edge of the mask to consider as in-focus region
    min_size: minimum size of objects to keep in the mask
    closing_radius: radius of the disk used for the closing operation
    """
    # 1. Segment fluorescent region
    img = normalize_image_global(img)
    thr = filters.threshold_otsu(img)
    colony_mask = img > thr

    colony_mask = morphology.remove_small_objects(colony_mask, min_size=min_size)
    colony_mask = morphology.binary_closing(colony_mask, morphology.disk(closing_radius))
    colony_mask = ndi.binary_fill_holes(colony_mask)

    # 2. Distance inward from edge
    dist_from_edge = ndi.distance_transform_edt(colony_mask) # distance_transform_edt gives distance to the nearest background pixel
    
    focus_mask = (colony_mask & (dist_from_edge >= d_inner))

    return focus_mask


def get_arcs_in_focus(focus_mask,oa_mask,ring_width=5):

    """
    Given a binary mask of the in-focus region of a colony, this function identifies the inner boundary of the mask and computes the distance of each pixel to this inner boundary. 
    It then divides the in-focus region into concentric rings of specified width and collects the pixel values and distances for each ring.
    Input:
    - focus_mask: 2D binary numpy array representing the in-focus region of the colony
    - ring_width: Width of the concentric rings (in pixels) to analyze
    Output:
    - ring_distances: List of distances corresponding to the center of each ring
    - pixel_values: List of numpy arrays containing the pixel values for each ring, ordered from left to right
    """

    # 1. Find the inner/bottom boundary of focus_mask
    focus_mask = focus_mask.astype(bool)
    eroded = ndi.binary_erosion(focus_mask)
    boundary = focus_mask & ~eroded

    # keep only bottom-facing boundary pixels
    Y, X = np.indices(focus_mask.shape)

    inner_boundary = np.zeros_like(focus_mask, dtype=bool)

    for x in range(focus_mask.shape[1]):
        ys = np.where(boundary[:, x])[0]
        if len(ys) > 0:
            inner_boundary[ys.max(), x] = True

    # distance to nearest inner-boundary pixel
    dist_from_inner_edge = ndi.distance_transform_edt(~inner_boundary) # by default gives distance to nearest zero pixel, so we invert the inner_boundary mask

    # keep values only inside focus_mask
    dist_from_inner_edge[~focus_mask] = np.nan

    # 2. For each ring, get the pixel values and distances
    pixel_values = []
    ring_distances = []

    min_dist = 50
    max_dist = 0.8*int(np.nanmax(dist_from_inner_edge))

    for d0 in range(min_dist, int(max_dist), ring_width):

        d1 = d0 + ring_width

        ring_mask = (
            focus_mask &
            (dist_from_inner_edge >= d0) &
            (dist_from_inner_edge < d1)
        )

        if ring_mask.sum() < ring_width:
            continue

        vals = oa_mask[ring_mask].astype(int)

        # order along the ring from left to right
        x = X[ring_mask]
        order = np.argsort(x)

        pixel_values.append(vals[order])
        ring_distances.append((d0 + d1) / 2)

    return ring_distances, pixel_values

def intermixing_index(pixel_values):

    """
    The shape of pixel_values is assumed to be (n_rings, n_pixels) where n_rings is the number of rings in the image and n_pixels is the number of pixels in each ring. 
    Each element of pixel_values should be either 0 or 1, indicating the presence or absence of Oa at that pixel.
    The intermixing index is calculated for each ring. 
    """

    # Calculate total abundance of Oa
    pix_flat = np.array([i for j in pixel_values for i in j])
    p = np.count_nonzero(pix_flat)/len(pix_flat)

    if (p <= 0.02) | (p >= 0.98):
        if type(pixel_values) == np.ndarray:
            int_index = [0 for i in range(pixel_values.shape[0])]
        else:
            int_index = [0 for i in range(len(pixel_values))]
    else:
        int_index = []
        for string in pixel_values:
            num_sectors = np.count_nonzero(np.convolve(np.array([-1,1]),string,mode="valid"))
            expected_num_sectors = 2*p*(1-p)*len(string)
            int_index.append(num_sectors/expected_num_sectors) # This way the intermixing index is between 0 and 1

    return int_index

def sector_size_and_freq(pixel_values):

    # Calculate total abundance of Oa
    pix_flat = np.array([i for j in pixel_values for i in j])
    freq = np.count_nonzero(pix_flat)/len(pix_flat)

    # For each row in the image, calculate the mean sector size
    # The mean sector size is obtained by dividing the total abundance of Oa in the row by the number of transitions between Oa and P
    if (freq == 0):
        mean_sector_size_oa = [0 for i in range(pixel_values.shape[0])]
        mean_sector_size_p = [0 for i in range(pixel_values.shape[0])]
    elif (freq == 1):
        mean_sector_size_oa = [1 for i in range(pixel_values.shape[0])]
        mean_sector_size_p = [0 for i in range(pixel_values.shape[0])]
    else:
        mean_sector_size_oa = []
        mean_sector_size_p = []
        for string in pixel_values:
            num_transitions = np.count_nonzero(np.convolve(np.array([-1,1]),string,mode="valid"))
            abund_oa = np.sum(string)
            abund_p = len(string) - abund_oa
            if abund_oa == 0:
                mean_sector_size_oa.append(0)
            elif abund_oa == len(string):
                mean_sector_size_oa.append(1)
            else:
                mean_sector_size_oa.append(2*abund_oa/num_transitions) 
            if abund_p == 0:
                mean_sector_size_p.append(0)
            elif abund_p == len(string):
                mean_sector_size_p.append(1)
            else:
                mean_sector_size_p.append(2*abund_p/num_transitions)

    return mean_sector_size_oa, mean_sector_size_p, freq