# -*- coding: utf-8 -*-
"""
Note:
    - this script performs nuclei matching from confocal to vala
    - and report the nuclei that are identified as positive

"""

import numpy as np
import pandas as pd
import cv2

from skimage import measure, morphology, segmentation, feature

from scipy import ndimage

import tifffile
import matplotlib.pyplot as plt

from statistics import mode
from matplotlib.patches import Circle

import os
import pickle
import json

# -- self define module and functions
import Functions as func

# Pixel size constants
XPIXEL = 0.22  # µm/pixel in xy
ZPIXEL = 0.35  # µm/pixel in z


class spot_class:
    def __init__(self, centroid, ellip_int):
        self.centroid  = centroid
        self.ellip_int = ellip_int

def plot_spots_on_image(image, centroids, nuc_mask, radius=4, color='red', vmin=None, vmax=None, save_path=None):
    """
    Plot the image with circles around each centroid.
    
    Parameters:
    - image: 2D numpy array of the raw nuclei image.
    - centroids: List of tuples, each containing (row, column) of the centroid.
    - radius: Radius of the circles to draw.
    - color: Color of the circles.
    - vmin: Minimum value for image display scaling.
    - vmax: Maximum value for image display scaling.
    - save_path: If provided, save the figure to this path instead of displaying it.
    """
    fig, ax = plt.subplots(dpi=300)
    ax.imshow(image, cmap='gray', vmin=vmin, vmax=vmax)
    ax.contour(nuc_mask, colors='yellow', linewidths=3.5)
    for centroid in centroids:
        circle = Circle((centroid[1], centroid[0]), radius, color=color, fill=False, linewidth=2)
        ax.add_patch(circle)
    plt.xticks([])  # Remove x-axis ticks
    plt.yticks([])  # Remove y-axis ticks
    if save_path:
        plt.savefig(save_path)
        plt.close()
    else:
        plt.show()

def imhmin(I, H):
    """
    Suppress regional minima in the grayscale image I whose depth is less than H.
    
    Parameters:
    - I: 2D numpy array (grayscale image)
    - H: scalar (depth threshold)
    
    Returns:
    - J: 2D numpy array (transformed image)
    """
    I = I.astype(np.float64)
    marker = I + H
    mask = I
    J = morphology.reconstruction(marker, mask, method='erosion')
    return J

def watershed_nuc(bw_nuc):
    
    # Distance transform and watershed for segmentation
    d = -ndimage.distance_transform_edt(bw_nuc)
    d = imhmin(d,2) # Suppress minima less than 2, this function is crucial for avoiding over segmentation
    d_inf = d.copy()
    d_inf[~bw_nuc] = -np.inf
    
    # Watershed segmentation
    return segmentation.watershed(d_inf, mask=bw_nuc, connectivity=2)

def intensity_ellipsoid(raw, position, Rxy=3, Rz=3):
    """
    Return the mean intensity value on the 3D-image raw over the elliptical region
    centered at position with xy-radius Rxy and z-radius Rz.
    
    Parameters:
    - raw (ndarray): 3D image array.
    - position (tuple): (x, y, z) coordinates of the ellipsoid center.
    - Rxy (float): Radius in xy-plane (default=3).
    - Rz (float): Radius in z-direction (default=3).
    
    Returns:
    - float: Mean intensity within the ellipsoid, or np.nan if the region is empty.
    """
    # Translate position to image coordinates (round to nearest integer)
    p = np.round(position).astype(int)
    
    # Get image dimensions
    sz = raw.shape
    
    # Create meshgrid of coordinates
    if Rz == 1:
        # 2D case: only X and Y coordinates
        Y, X = np.ogrid[0:sz[1], 0:sz[2]]  # Note: Python uses 0-based indexing
        Z = None
    else:
        # 3D case: X, Y, Z coordinates
        Z, Y, X = np.ogrid[0:sz[0], 0:sz[1], 0:sz[2]]
    
    # Build mask with ellipse equation
    if Z is None:
        # 2D ellipse
        mask = ((X - p[1]) / Rxy)**2 + ((Y - p[0]) / Rxy)**2 <= 1
    else:
        # 3D ellipsoid
        mask = ((X - p[2]) / Rxy)**2 + ((Y - p[1]) / Rxy)**2 + ((Z - p[0]) / Rz)**2 <= 1
    
    # Extract intensities within the mask
    INT = raw[mask]
    if INT.size == 0:
        return np.nan
    else:
        return np.nanmean(INT)

def spot_detection(img_nuc, mask_nuc, rp3D_nuc, thred_3D, nuc_label, upper_cap=200):
    
    # %%
        
    # img_nuc_blur = ndimage.gaussian_filter(img_nuc, [2,3.5,3.5])
    img_nuc_blur = np.zeros(img_nuc.shape)
    
    img_nuc_subs = img_nuc.astype(float) - img_nuc_blur.astype(float)
    
    # Create 3D nuclear mask and dilate    
    ch2_mask_3d = rp3D_nuc.bw3D
    
    ft = morphology.white_tophat(img_nuc_subs, morphology.ball(2))
    
    # ## -- for the new method, 2025-12-25 
    threshold = np.percentile(ft.ravel(), 99.9)

    
    bw = ft >= threshold
    bw = bw & ch2_mask_3d # Mask out non-nuclear regions
    # Exclude first frame and last frame
    bw[0, :, :] = bw[-1, :, :] = False
    for i in range(bw.shape[0]): # Exclude the edges
        bw[i, 0, :] = bw[i, -1, :] = bw[i, :, 0] = bw[i, :, -1] = False
    
    # Label spots and filter by area (>= 25)
    labeled_bw = measure.label(bw)  # regular segmentation based on threshold
    
    spots = measure.regionprops(labeled_bw, intensity_image=img_nuc)
    spots = [spot for spot in spots if (spot.area >= 3)]
    
    # ##########################
    # This block of code deals with spots that are close to each otehr and can be potentially separated
    #    1. check if there are more than one local maxiums. 
    #    2. if so, use watershed to separate two spots
    # ##########################
    boolean_list = [True] * len(spots)
    spots_new    = []
    for i in range(len(spots)):

        s = spots[i]

        # sigma      = [  1, 0.5, 0.5]
        sigma      = [0.6, 0.6, 0.6]
        local_maxi = feature.peak_local_max(ndimage.gaussian_filter(s.image_intensity, sigma))
        
        if len(local_maxi) > 1:
            
            print('### In spot watershed: Number of local maxium: %d' % len(local_maxi))
            
            boolean_list[i] = False
            
            ## Create mask for selected spots
            mask       = np.zeros(bw.shape)
            coord_list = spots[i].coords
            for k in range(coord_list.shape[0]):
                mask[*coord_list[k]]=1
            bw_select  = mask
            
            ## Create markers for watersheding the spots
            local_maxi2 = feature.peak_local_max(bw_select * ndimage.gaussian_filter(img_nuc, sigma)) # finding the coordinates of peaks
            peaks_mask  = np.zeros_like(bw_select, dtype=bool) # making the peak labels
            for p_coord in local_maxi2:
                peaks_mask[*p_coord] = True
            markers     = ndimage.label(peaks_mask)[0] 
            
            ## watershed spots
            dist         = -ndimage.distance_transform_edt(bw_select)
            # dist         = imhmin(dist,1)
            d_inf        = dist.copy()
            d_inf[~bw]   = -np.inf
            label_select = segmentation.watershed(d_inf, markers, mask=bw_select, connectivity=1)
            
            ## Collect new spots
            s_new = []
            s_new = measure.regionprops(label_select, intensity_image=img_nuc)
            s_new = [s for s in s_new if s.area >= 10]
            spots_new.extend(s_new)
            
            print("###    Got %d spots" % len(s_new))
            
    # Update spots list
    spots2 = []
    for i in range(len(spots)):
        if boolean_list[i]:
            spots2.append(spots[i])
    spots = spots2
    spots.extend(spots_new)
    
    # %% Compute spot intensities using intensity_ellipsoid
    #    keep intensity and centroids for a spot
    spot_new = [ [] ] * len(spots)
    for i, spot in enumerate(spots):
        centroid = spot.weighted_centroid  # (z, y, x)
        ellip_int = intensity_ellipsoid(img_nuc, centroid, Rxy=2, Rz=2)
        # spot.ellip_int = ellip_int
        
        spot_new[i] = spot_class(centroid, ellip_int)
    spots = spot_new
    
    # %% ** Plot each nuclei and detected spots
    # img_nuc_max = np.max(img_nuc[7:29,:,:], axis=0)
    nuc_name = f'nuc_%d' % nuc_label
    
    img_nuc_max = np.max(img_nuc, axis=0)
    centroids   = [ np.round([s.centroid[1],s.centroid[2]]).astype(int) for s in spots ]
    ch2_mask    = morphology.binary_dilation(mask_nuc, morphology.disk(3))
    plot_spots_on_image(img_nuc_max,centroids,ch2_mask, save_path=file_folder+os.sep+r'segmented_nuclei'+os.sep+nuc_name + r'_ch1_ref.png') #,vmin=800,vmax=2000)
    # plot_spots_on_image(img_nuc_max,centroids,np.any(ch2_mask_3d,axis=0)) #,vmin=800,vmax=2000)
    
    # %%
    # Background intensity calculation
    img_nuc_int = img_nuc.astype(float)
    inner_mask  = morphology.binary_dilation(bw, morphology.ball(3)) # was 3
    outer_mask  = morphology.binary_dilation(bw, morphology.ball(4)) # was 6, a bit tighter compared to spot2
    # outer_mask  = morphology.binary_dilation(bw, morphology.ball(9))
    img_nuc_int[inner_mask]       = np.nan  # Exclude spot regions, using nan, because if use False it becomes 0 value
    img_nuc_int[~ (ch2_mask_3d 
                   & outer_mask)] = np.nan  # Keep only nuclear bg region
    spot_bg     = np.nanmean(img_nuc_int)
    
    return spots, spot_bg


def spot_detection2(img, centroid):
 
    ellip_int = intensity_ellipsoid(img, centroid, Rxy=2, Rz=2)
    
    mask       = np.zeros(img.shape)
    mask[int(centroid[0]),int(centroid[1]),int(centroid[2])] = 1
    
    inner_mask  = morphology.binary_dilation(mask, morphology.ball(3)) # was 3
    outer_mask  = morphology.binary_dilation(mask, morphology.ball(6)) # was 6
    
    img_nuc_int = img.astype(float)
    img_nuc_int[inner_mask]  = np.nan
    img_nuc_int[~outer_mask] = np.nan
    spot_bg   = np.nanmean(img_nuc_int)
    
    
    return ellip_int, spot_bg


def nuc_segmentation(img_nuc, percentile):
    thred   = np.percentile(img_nuc, percentile)
        
    bw_nuc = img_nuc >= thred
    bw_nuc = ndimage.binary_fill_holes(bw_nuc)
    bw_nuc = morphology.binary_opening(bw_nuc, morphology.disk(5))
    
    # Watershed segmentation of nuclei
    seg_nuc = watershed_nuc(bw_nuc)
    rps_nuc = measure.regionprops(seg_nuc)
    
    # %
    fig, ax = plt.subplots(dpi=300)
    plt.imshow(seg_nuc)
    plt.xticks([])  # Remove x-axis ticks
    plt.yticks([])  # Remove y-axis ticks
    plt.close()
    
    # % filer by area
    rps_nuc = [rp for rp in rps_nuc if  200 <= rp.area <= 5000] # for 40x objective, using for finding correlation only
    
    # % Nuclear mask after size filtering
    mask = np.zeros_like(bw_nuc)
    seg_nuc_select_label = np.zeros(bw_nuc.shape, dtype=int)
    for rp in rps_nuc:
        cent = np.round(rp.centroid).astype(int)
        mask[*cent] = 1
        for coord in tuple(rp.coords): # recovering label image
            seg_nuc_select_label[*coord] = rp.label
    
    fig, ax = plt.subplots(dpi=300)
    plt.imshow(seg_nuc_select_label)
    plt.xticks([])  # Remove x-axis ticks
    plt.yticks([])  # Remove y-axis ticks
    plt.close()

    return seg_nuc_select_label


def nuc_segmentation3D(img_dapi, percentile=95):
    
    # % Projecting onto a 2D plane
    img_nuc = np.max(img_dapi, axis=0)
    
    thred  = np.percentile(img_nuc, percentile)
        
    bw_nuc = img_nuc >= thred
    bw_nuc = ndimage.binary_fill_holes(bw_nuc)
    bw_nuc = morphology.binary_opening(bw_nuc, morphology.disk(5))
    
    # % Watershed segmentation of nuclei
    seg_nuc = watershed_nuc(bw_nuc)
    rps_nuc = measure.regionprops(seg_nuc)

    # rps_nuc = [rp for rp in rps_nuc if  500 <= rp.area <= 10000]
    rps_nuc = [rp for rp in rps_nuc if  150 <= rp.area <= 5000] # for 40x objective
    
    
    # % Nuclear mask after size filtering
    # Create mask of the original spot
    mask = np.zeros_like(bw_nuc)
    seg_nuc_select_label = np.zeros(bw_nuc.shape, dtype=int)
    for rp in rps_nuc:
        cent = np.round(rp.centroid).astype(int)
        mask[*cent] = 1
        for coord in tuple(rp.coords): # recovering label image
            seg_nuc_select_label[*coord] = rp.label

    fig, ax = plt.subplots(dpi=300)
    plt.imshow(seg_nuc_select_label)
    for rp in rps_nuc:
        cent = np.round(rp.centroid).astype(int)
        plt.text(cent[1], cent[0], rp.label, color="white", fontsize=6)
    plt.xticks([])  # Remove x-axis ticks
    plt.yticks([])  # Remove y-axis ticks
    plt.savefig(file_folder+os.sep+sf+r'_sgGOLDFISH_nuc-labels.png', format='png')
    plt.close()

    # % Turn the mask into a cylinder, carve out the nuclear image
    nuc_label_3D = np.repeat(seg_nuc_select_label[np.newaxis,:,:], img_dapi.shape[0], axis=0)
    rps_nuc_3D   = measure.regionprops(nuc_label_3D, img_dapi)

    # % Creating 3D mask for each nucleus
    for k in range(len(rps_nuc_3D)):
        
        rp = rps_nuc_3D[k]
        
        img_single = ndimage.gaussian_filter(rp.image_intensity, [1.2,3,3])
        thread     = np.percentile(img_single.ravel(), 70)
        bw_single3D  = rp.image_intensity > thread
        for i in range(bw_single3D.shape[0]):
            bw_single3D[i,:,:] = ndimage.binary_fill_holes(bw_single3D[i,:,:])
            bw_single3D[i,:,:] = morphology.binary_opening(bw_single3D[i,:,:], morphology.disk(4))
            bw_single3D[i,:,:] = morphology.binary_erosion(bw_single3D[i,:,:], morphology.disk(4))
        
        rps_nuc_3D[k].bw3D = bw_single3D

    return rps_nuc, rps_nuc_3D, seg_nuc_select_label


def process_file(img, percentile=75, Otsu=0):
    """
    Process a single TIFF file for nuclei segmentation and spot detection.
    
    Parameters:
    - file_path (str): Path to the TIFF file.
    - percentile (float): Percentile for nuclear segmentation threshold (default=75).
    
    Returns:
    - list: List of dictionaries containing nuclei with spots and background intensities.
    """
    
    # Get the nuclei channel
    ch_nuc = img[:,-1, :, :]  # Nuclear channel, usually last stack. Our widefield
    
    pr97 = np.percentile(ch_nuc.ravel(), 97)
    # Get the nuclei segmentaion objects
    # rps_nuc, rps_nuc_3D = NucSeg_SpotDetect.nuc_segmentation(ch_nuc, percentile, Otsu=Otsu)
    rps_nuc, rps_nuc_3D, seg_nuc = nuc_segmentation3D(ch_nuc, percentile)
    print(f"#      Number of nuclei in this FOV: %d" % len(rps_nuc))
        
    ch_range = 1 # In this experiment, we use Ch1 as a reference to check the intensity of Ch2
    
    ch_spot  = img[:,0, :, :] # ch1, the reference channel
    ch_spot2 = img[:,1, :, :] # ch2, the detection channel

    slide_bg = mode(np.ravel(ch_spot))
    print(f"##     Nuclear background: {slide_bg:.2f}")
    
    nuclei = []
 
    # count = 0
    for idx, rp in enumerate(rps_nuc):
        # Crop ch_spot using bounding box
        min_row, min_col, max_row, max_col = rp.bbox
        ch_spot_cropped  = ch_spot[:, min_row:max_row, min_col:max_col]
        ch_spot_cropped2 = ch_spot2[:, min_row:max_row, min_col:max_col]
        mask = rp.image  # 2D binary mask of the nucleus
        
        # Get the spot objects
        print(f"### Nuclei label: {rp.label:3d}")
        upper_cap = 300 # 200
            
        spots , bg_int  = spot_detection(ch_spot_cropped , mask, rps_nuc_3D[idx], pr97, rp.label, upper_cap)
        
        spots_detec  = []
        nuc_bg_detec = []
        for s in range(len(spots)):
            spots2, bg_int2 = spot_detection2(ch_spot_cropped2, spots[s].centroid)
            spots_detec.append(spots2)
            nuc_bg_detec.append(bg_int2)
        
        
        fig, ax = plt.subplots(dpi=300)
        ax.imshow(np.max(ch_spot_cropped2,0), cmap='gray')
        plt.xticks([])  # Remove x-axis ticks
        plt.yticks([])  # Remove y-axis ticks
        nuc_name = f'nuc_%d' % rp.label
        plt.savefig(file_folder+os.sep+r'segmented_nuclei'+os.sep+nuc_name + r'_ch2_detect.png', format='png')
        plt.close()
        
        # Store nucleus data
        nuclei.append({
            'spots_ref'   : spots,                # List of spot objects
            'nuc_bg_ref'  : bg_int,               # Nuclear background intensity
            'spots_detec' : spots_detec,
            'nuc_bg_detec': nuc_bg_detec,
            'slide_bg'  : slide_bg,             # Slide background intensity
            'nuc_img_sz': ch_spot_cropped.shape,# Nuclear image size
            'nuc_label' : rp.label
            })
       
    numb_spot_per_nuc = np.mean([len(nuc['spots_ref']) for nuc in nuclei])
    print(f"### Number of spots per nucleus: {numb_spot_per_nuc:.2f}")

    # nuclei_obj_ch.append(nuclei)
    
    return nuclei, seg_nuc


if __name__ == "__main__":
    
# 2025-11-26 N406K
    # file_folder_root = r'G:\Shared drives\Po-Ta_Chen_Projects\Data\Pota-processed\20251122_N406K_Nikon_Ti2_spinning-disk'    
    file_folder_root = r'G:\Shared drives\Cardiomyocyte project\20251122_N406K_Nikon_Ti2_spinning-disk'
    # sub_folder = [r'D05-1', r'D05-2', r'D05-3', r'D06-1', r'D06-2', r'D06-3', r'D07-1', r'D07-2', r'D07-3',
    #               r'D08-1', r'D08-2', r'D08-3', r'D09-1', r'D09-2', r'D09-3',
    #               r'E03-1', r'E03-2', r'E03-3', r'E04-1', r'E04-2', r'E04-3', r'E05-1', r'E05-2', r'E05-3',
    #               r'E06-1', r'E06-2', r'E06-3', r'E07-1', r'E07-2', r'E07-3', r'E08-1', r'E08-2', r'E08-3',
    #               r'E09-1', r'E09-2', r'E09-3', r'E10-1', r'E10-2', r'E10-3'  ]
    
    sub_folder = [r'E08-2', r'E08-3',
                  r'E09-1', r'E09-2', r'E09-3', r'E10-1', r'E10-2', r'E10-3'  ]
    
    # # sub_folder = [r'D05-1']
    

# # 2025-12-15 20251212-E4950K-2
    # file_folder_root = r'G:\Shared drives\Po-Ta_Chen_Projects\Data\Pota-processed\20251212-E4950K-2'
    # sub_folder = [r'D04-1', r'D04-2', r'D04-3', r'D05-1', r'D05-2', r'D05-3', r'D06-1', r'D06-2', r'D06-3',
    #               r'D07-1', r'D07-2', r'D07-3', r'D08-1', r'D08-2', r'D08-3', r'D09-1', r'D09-2', r'D09-3',
    #               r'D10-1', r'D10-2', r'D10-3',
    #               r'E03-1', r'E03-2', r'E03-3', r'E04-1', r'E04-2', r'E04-3', r'E05-1', r'E05-2', r'E05-3',
    #               r'E06-1', r'E06-2', r'E06-3', r'E07-1', r'E07-2', r'E07-3', r'E08-1', r'E08-2', r'E08-3',
    #               r'E09-1', r'E09-2', r'E09-3', r'E10-1', r'E10-2', r'E10-3',
    #               r'F03-1', r'F03-2', r'F03-3', r'F04-1', r'F04-2', r'F04-3', r'F05-1', r'F05-2', r'F05-3' ]
   
# # 2025-12-16 20251212-E4950K-1
    # file_folder_root = r'G:\Shared drives\Po-Ta_Chen_Projects\Data\Pota-processed\20251212-E4950K-1'
    # sub_folder = [r'D06-1', r'D06-2', r'D06-3', r'D07-1', r'D07-2', r'D07-3', r'D08-1', r'D08-2', r'D08-3',
    #               r'D09-1', r'D09-2', r'D09-3',
    #               r'E04-1', r'E04-2', r'E04-3', r'E05-1', r'E05-2', r'E05-3', r'E06-1', r'E06-2', r'E06-3',
    #               r'E07-1', r'E07-2', r'E07-3', r'E08-1', r'E08-2', r'E08-3', r'E09-1', r'E09-2', r'E09-3',
    #               r'E10-1', r'E10-2', r'E10-3',
    #               r'F03-1', r'F03-2', r'F03-3', r'F04-1', r'F04-2', r'F04-3', r'F05-1', r'F05-2', r'F05-3' ]

# 2025-12-24 20251220-21-A320V-AP
    # # file_folder_root = r'G:\Shared drives\Po-Ta_Chen_Projects\Data\Pota-processed\20251220-21-A320V-AP'
    # file_folder_root = r'G:\Shared drives\Cardiomyocyte project\20251220-21-A320V-AP'
    # sub_folder = [r'D03-1', r'D03-2', r'D03-3', r'D04-1', r'D04-2', r'D04-3', r'D05-1', r'D05-2', r'D05-3',
    #               r'D06-1', r'D06-2', r'D06-3', r'D07-1', r'D07-2', r'D07-3', r'D08-1', r'D08-2', r'D08-3',
    #               r'D09-1', r'D09-2', r'D09-3', r'D10-1', r'D10-2', r'D10-3',
    #               r'E03-1', r'E03-2', r'E03-3', r'E04-1', r'E04-2', r'E04-3', r'E05-1', r'E05-2', r'E05-3',
    #               r'E06-1', r'E06-2', r'E06-3', r'E07-1', r'E07-2', r'E07-3', r'E08-1', r'E08-2', r'E08-3', 
    #               r'E09-1', r'E09-2', r'E09-3', r'E10-1', r'E10-2', r'E10-3',
    #               r'F03-1', r'F03-2', r'F03-3', r'F04-1', r'F04-2', r'F04-3', r'F05-1', r'F05-2', r'F05-3',
    #               r'F06-1', r'F06-2', r'F06-3']
  
# 2025-12-24 20251220-21-A320V-Calcium
    # # file_folder_root = r'G:\Shared drives\Po-Ta_Chen_Projects\Data\Pota-processed\20251220-21-A320V-Calcium'
    # file_folder_root = r'G:\Shared drives\Cardiomyocyte project\20251220-21-A320V-Calcium'
    # # sub_folder = [r'D03-1', r'D03-2', r'D03-3', r'D04-1', r'D04-2', r'D04-3', r'D05-1', r'D05-2', r'D05-3',
    # #               r'D06-1', r'D06-2', r'D06-3', r'D07-1', r'D07-2', r'D07-3', r'D08-1', r'D08-2', r'D08-3',
    # #               r'D09-1', r'D09-2', r'D09-3', r'D10-1', r'D10-2', r'D10-3',
    # #               r'E03-1', r'E03-2', r'E03-3', r'E04-1', r'E04-2', r'E04-3', r'E05-1', r'E05-2', r'E05-3',
    # #               r'E06-1', r'E06-2', r'E06-3', r'E07-1', r'E07-2', r'E07-3', r'E08-1', r'E08-2', r'E08-3', 
    # #               r'E09-1', r'E09-2', r'E09-3', r'E10-1', r'E10-2', r'E10-3',
    # #               r'F04-1', r'F04-2', r'F04-3', r'F05-1', r'F05-2', r'F05-3', r'F06-1', r'F06-2', r'F06-3']    

# 2025-12-28 20251227-N406K
    # file_folder_root = r'G:\Shared drives\Cardiomyocyte project\20251227-N406K'
    # sub_folder = [r'D03-1', r'D03-2', r'D03-3', r'D04-1', r'D04-2', r'D04-3', r'D05-1', r'D05-2', r'D05-3',
    #               r'D07-1', r'D07-2', r'D07-3', r'D08-1', r'D08-2', r'D08-3',
    #               r'D09-1', r'D09-2', r'D09-3', r'D10-1', r'D10-2', r'D10-3',
    #               r'E03-1', r'E03-2', r'E03-3', r'E04-1', r'E04-2', r'E04-3', r'E05-1', r'E05-2', r'E05-3',
    #               r'E07-1', r'E07-2', r'E07-3', r'E08-1', r'E08-2', r'E08-3', 
    #               r'E09-1', r'E09-2', r'E09-3', r'E10-1', r'E10-2', r'E10-3',
    #               r'F03-1', r'F03-2', r'F03-3', r'F04-1', r'F04-2', r'F04-3', r'F05-1', r'F05-2', r'F05-3']
    
    
    for sf in sub_folder:
        
        # print("In folder: %s" % sf)
        
        file_folder = file_folder_root+os.sep+sf
        
        has_file, name = func.file_contains_string(file_folder, 'positive_label', 0)

        
        force_run = 1

        if has_file & ~force_run:
            
            print('Skipping %s' % sf)
            continue
            
        else:
            
            print('Working on %s' % sf)
            os.makedirs(file_folder+os.sep+r'segmented_nuclei', exist_ok=True)
            
            if os.path.exists(file_folder+os.sep+sf+r'_sgGOLDFISH_positive_label.json'):
                os.remove(file_folder+os.sep+sf+r'_positive_label.pkl')
                os.remove(file_folder+os.sep+sf+r'_sgGOLDFISH_mask_selected.pkl')
                os.remove(file_folder+os.sep+sf+r'_sgGOLDFISH_nuc-labels.png')
                os.remove(file_folder+os.sep+sf+r'_sgGOLDFISH_positive_label.json')
                os.remove(file_folder+os.sep+sf+r'_ch2_int_hist.svg')
                
                print(f"The files have been deleted.")
            
            
            try:
            
                # Load TIFF image (assumes shape: (depth, height, width))
                img_Nikon_Ti2_spinning_disk = tifffile.imread(file_folder+os.sep+r'Nikon-Ti2'+os.sep+r'*.tif')
                img_vala_nuc_mask           = tifffile.imread(file_folder+os.sep+r'Vala_mask'+os.sep+r'*.tif')
                img_vala_nuc_mask           = img_vala_nuc_mask[::-1,::-1]
                
            
                # %% Get the nuclei channel
                ch_nuc   = img_Nikon_Ti2_spinning_disk[:,-1, :, :]  # Nuclear channel, usually last stack. Nikon Ti2 spinning disk
                img_nuc1 = np.max(ch_nuc, axis=0)
                # img_nuc2 = img_vala_nuc2
            
                # %%
                # Get the nuclei segmentaion objects
                # nuc_label1 = nuc_segmentation(img_nuc1, 98)
                nuc_label1 = nuc_segmentation(img_nuc1, 90)
                # nuc_label2 = nuc_segmentation(img_nuc2, 95)
                nuc_label2 = img_vala_nuc_mask
                
                
                # %% Finding the peak correlation between 2 images
                # corr = cv2.filter2D(nuc_label1, ddepth=-1, kernel=nuc_label2)
                bw1 = nuc_label1>0
                bw2 = nuc_label2>0
                corr = cv2.filter2D( bw1.astype(int), ddepth=-1, kernel=bw2.astype(int))
                
                fig, ax = plt.subplots(dpi=300)
                plt.imshow(corr)
                plt.xticks([])  # Remove x-axis ticks
                plt.yticks([])  # Remove y-axis ticks
                plt.close()
                
                corr_peak = np.where(corr==np.max(corr))
                shift_x   = corr_peak[0][0]-nuc_label2.shape[0]/2
                shift_y   = corr_peak[1][0]-nuc_label2.shape[1]/2
                shift_x = shift_x.astype(int)
                shift_y = shift_y.astype(int)
                
               
                # %%
                img_Nikon_crop = img_Nikon_Ti2_spinning_disk[:,:,(shift_x+1):(shift_x+nuc_label2.shape[0]+1), (shift_y+1):(shift_y+nuc_label2.shape[1]+1)]
                
                
                # %% Nucleus and spot detection
                # nuclei_perFOV = []
                nuclei, seg_nuc = process_file(img_Nikon_crop, percentile=90)
                
                
                # %% 
                nuc_ref = [ n for n in nuclei if (len(n['spots_ref'])>0) and (len(n['spots_ref'])<3) ]
                
            
                # %%    
                spot_refer_int = []
                spot_detec_int = []
                nuc_label = []
                
                for n in nuc_ref:
                    
                    for i in range(len(n['spots_ref'])):
                        spot_refer_int.append(n['spots_ref'][i].ellip_int-n['nuc_bg_ref'])
                        spot_detec_int.append(n['spots_detec'][i]-n['nuc_bg_detec'][i])
                        nuc_label.append(n['nuc_label'])
            
                nuc_df = pd.DataFrame({'spots_ref': spot_refer_int, 'spots_detec': spot_detec_int, 'nuc_label': nuc_label})
                
                fig, ax = plt.subplots(dpi=300)
                plt.hist(nuc_df['spots_detec'],30)
                plt.savefig(file_folder+os.sep+sf+r'_ch2_int_hist.svg', format='svg')
                plt.close()
                
                # nuc_df_positive = nuc_df[nuc_df['spots_detec']>100]
                nuc_df_positive = nuc_df[nuc_df['spots_detec']>150]
                
                
                # %%
                nuclei_pick_centroids = []
                for i in range(len(nuc_df_positive)):        
                    nuc_select_bw = seg_nuc==nuc_df_positive['nuc_label'].iloc[i]
                    rp = measure.regionprops(nuc_select_bw.astype(int))
                    
                    nuclei_pick_centroids.append(rp[0].centroid)  
                
                label_vala = []
                for i in range(len(nuclei_pick_centroids)):
                    label_vala.append(img_vala_nuc_mask[ np.round(nuclei_pick_centroids[i][0]).astype(int), 
                                                         np.round(nuclei_pick_centroids[i][1]).astype(int) ])
            
                nuc_df_positive = nuc_df_positive.assign(label_vala=pd.Series(label_vala, index=nuc_df_positive.index))
                
                
                select_vale_label_list = nuc_df_positive[nuc_df_positive['label_vala']!=0]['label_vala'].to_list()
                select_nuc_label_list  = nuc_df_positive[nuc_df_positive['label_vala']!=0]['nuc_label'].to_list()
                print('Select list (vala):')
                print(select_vale_label_list)
                print('Select list (sgGF):')
                print(select_nuc_label_list)
                with open(file_folder+os.sep+sf+r'_positive_label.pkl', 'wb') as file:
                    pickle.dump(select_vale_label_list, file)
                with open(file_folder+os.sep+sf+r'_sgGOLDFISH_positive_label.json', 'w') as file:
                    json.dump(select_nuc_label_list, file)
                
                
                label_all_list = np.unique(seg_nuc)
                
                seg_nuc_select = seg_nuc.copy()
                
                for i in [ x for x in label_all_list if x not in select_nuc_label_list]:
                    seg_nuc_select[seg_nuc_select==i]=0
                    
                with open(file_folder+os.sep+sf+r'_sgGOLDFISH_mask_selected.pkl', 'wb') as file:
                    pickle.dump(seg_nuc_select, file)
                
                
                # deleting data instances
                img_Nikon_Ti2_spinning_disk = []
                img_Nikon_crop = []
                img_vala_nuc_mask = []
                ch_nuc = []
                nuc_label1 = []
                nuc_label2 = []
                nuclei = []
                seg_nuc = []
                nuc_ref = []
                    
            except Exception as e: 
                
                print(e) 
                print("Errors in folder: %s" % sf)
    
    
    
    
    
    