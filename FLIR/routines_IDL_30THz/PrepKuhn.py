import os
import sys
import glob
import argparse
import numpy as np
from astropy.io import fits
from skimage.filters import roberts
from skimage.measure import EllipseModel, ransac
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from datetime import datetime
from KLL import mk_flat
import pdb 
import time
from CraamTools.fit import Circle

# --- Define working global variables ---
img_count_for_cube = 100   # Number of images for the data cube
mask_radius        = 330.0 # Radius in pixels to eliminate the solar limb
threshold          = 0.1   # Threshold for valid data in the Kuhn-Lin-Loranz routine

def Test(Show=False):

    iterations = 10
    mask_height, mask_width = 1920, 2560
    cy, cx = mask_height / 2, mask_width / 2
    y, x = np.ogrid[:mask_height, :mask_width]
    dist_from_center = np.sqrt((x - cx)**2 + (y - cy)**2)
    large_mask = (dist_from_center < mask_radius).astype(np.float32)    
    FitsList = glob.glob('./flats/*fits')
    FitsList.sort()
    indices=[0,13,26,28,30]

     # Initialize data containers
    image_cube = np.zeros((480,640,len(indices)), dtype=np.float32)
    displacements = np.zeros((2, len(indices)), dtype=np.float32)
    
    img_idx_in_cube = 0
    img_read_count  = 0
    
    for i in indices:

        FitsFile = FitsList[i]
        with fits.open(FitsFile) as hdul:
            image_data = hdul[0].data.astype(np.float32)
            # Normalize FITS images
        image_data -= np.min(image_data)
        max_val = np.max(image_data)
        if max_val > 0:
            image_data /= max_val

        # --- Find the center of the solar disk ---
        # Identify bright, non-saturated pixels as part of the disk
        max_val = np.max(image_data)
        disk_mask = image_data > (max_val - max_val / 5.0)
        
        if not np.any(disk_mask):
            continue # Skip saturated or dark images

        # Use Roberts filter to find the limb (edge)
        limb_image = roberts(disk_mask)
        limb_points = np.argwhere(limb_image > 0)
        
        if limb_points.shape[0] < 50: # Need enough points to fit a circle
            continue

        # Fit a circle to the limb points
        y_coords, x_coords = limb_points[:, 0], limb_points[:, 1]
        points = np.column_stack([x_coords, y_coords])
         
        try:
            # aistudio.gemini proposed the use of the library skyimage.measure.CircleModel() / EllipseModel()
            # The former cannot fit ANY image, the latter fits with great differences respect to the IDL fitting
            # I decided to use my Circle.py fit with results closer to IDL
            #            model = CircleModel()
            #            model.estimate(points)
            #            xc, yc, r, _ = model.params
            ######################################################################################################
            par,_ = Circle.fit(x_coords,y_coords)
            xc = par[0]
            yc = par[1]
            r = par[2]

        except Exception:
            print(f"  Failed to fit the limb")
            continue # Skip if fitting fails

        # --- Validate the found center and add image to the cube ---
        print(f"Current FITS {FitsFile}, Radius: {r:.0f}  Center: ({xc:.1f}, {yc:.1f})")

        ixc, iyc = int(round(xc)), int(round(yc))
        displacements[:, img_idx_in_cube] = [xc, yc]        
        mask_center_x = int(large_mask.shape[1] / 2)
        mask_center_y = int(large_mask.shape[0] / 2)
        
        start_x = mask_center_x - ixc
        start_y = mask_center_y - iyc
            
        sliced_mask = large_mask[start_y : start_y + 480, start_x : start_x + 640]
            
        # Apply mask and store in the data cube
        masked_image = image_data * sliced_mask
        masked_image[masked_image == 0] = 0.001 # Avoid division by zero
        image_cube[:, :, img_idx_in_cube] = masked_image
        img_idx_in_cube+=1

    # --- Calculate the first-pass flat using the collected data cube ---
    flat_pass = mk_flat(image_cube, threshold, displacements, iterations,Show=Show)
    flat_pass[flat_pass == 0] = 1.0
    
    
    return flat_pass
