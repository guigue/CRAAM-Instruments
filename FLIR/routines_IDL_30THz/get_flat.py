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

import time


# --- Define working global variables ---
img_count_for_cube = 100   # Number of images for the data cube
mask_radius        = 330.0 # Radius in pixels to eliminate the solar limb
threshold          = 0.1   # Threshold for valid data in the Kuhn-Lin-Loranz routine
#iterations         = 10    # Number of iterations for the Kuhn-Lin-Loranz routine

def flatSave(telescope,flat,outdir,fname, date_obs):

    # Create FITS header
    hdu = fits.PrimaryHDU(flat)
    hdr = hdu.header
    hdr['DATE'] = (datetime.now().strftime('%Y-%m-%d'), 'Date of file creation')
    
    # Add telescope-specific keywords
    if telescope == 'AR30T':
        hdr['FILENAME'] = (f'{fname}.fits', 'Output filename')
        hdr['DATE-OBS'] = (date_obs, 'Date of observation')
        hdr['TELESCOP'] = ('AR30T', 'Telescope ID')
        hdr['INSTRUME'] = ('30 THz camera', 'Instrument ID')
        hdr['DETECTOR'] = ('FLIR A645sc', 'Detector ID')
        hdr['EXPTYPE']  = ('FLAT', 'Exposure Type')
        hdr['WAVELNTH'] = (10, 'Wavelength in micrometers')
        hdr['OBSERVAT'] = ('Obs. Astronomico Felix Aguilar', 'Observatory')
        hdr['PLACE']    = ('El Leoncito - San Juan - ARGENTINA', 'Location')
        hdr['LONGITUD'] = ('-69 19.8', 'Observatory Longitude')
        hdr['LATITUDE'] = ('-31 48.1', 'Observatory Latitude')
    elif telescope == 'BR30T':
        hdr['FILENAME'] = (f'{fname}.fits', 'Output filename')
        hdr['DATE-OBS'] = (date_obs, 'Date of observation')
        hdr['TELESCOP'] = ('BR30T', 'Telescope ID')
        hdr['INSTRUME'] = ('30 THz camera', 'Instrument ID')
        hdr['DETECTOR'] = ('FLIR A20', 'Detector ID')
        hdr['EXPTYPE']  = ('FLAT', 'Exposure Type')
        hdr['WAVELNTH'] = (10, 'Wavelength in micrometers')
        hdr['OBSERVAT'] = ('CRAAM', 'Observatory')
        hdr['PLACE']    = ('São Paulo - BRAZIL', 'Location')
        hdr['LONGITUD'] = ('-46 39.1', 'Observatory Longitude')
        hdr['LATITUDE'] = ('-23 32.8', 'Observatory Latitude')
        
    fits_filepath = os.path.join(outdir, f'{fname}.fits')
    hdu.writeto(fits_filepath, overwrite=True)
    print(f"Successfully saved FITS file: {fits_filepath}")

    return

def flatPass(iterations,image_files,semilla,large_mask,Show=False):
    
    # Initialize data containers
    image_cube = np.zeros((480, 640, img_count_for_cube), dtype=np.float32)
    displacements = np.zeros((2, img_count_for_cube), dtype=np.float32)
    
    # Generate a random sequence of images to process
    total_images    = len(image_files)
    np.random.seed(semilla)
    random_indices = np.random.choice(total_images, size=total_images, replace=False)

    img_idx_in_cube = 0
    img_read_count  = 0
    while img_idx_in_cube < img_count_for_cube and img_read_count < total_images:
        
        # --- Read an image from the random sequence ---
        current_file = image_files[random_indices[img_read_count]]
        img_read_count += 1

        with fits.open(current_file) as hdul:
            image_data = hdul[0].data.astype(np.float32)
            # Normalize FITS images
        image_data -= np.min(image_data)
        max_val = np.max(image_data)
        if max_val > 0:
            image_data /= max_val

        # --- Find the center of the solar disk ---
        # Identify bright, non-saturated pixels as part of the disk
        max_val = np.max(image_data)
        disk_mask = image_data > (max_val - max_val / 7.0)
        
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
            model = EllipseModel()
            model.estimate(points)
            xc, yc, rx, ry, _ = model.params
            r = (rx+ry)*0.5

        except Exception:
            print(f"  Failed to fit the limb")
            continue # Skip if fitting fails

        # --- Validate the found center and add image to the cube ---
        if xc > 100 and xc < 359 and yc > 100 and yc < 379 and r > 350:
            
            #print(f"  Center: ({xc:.1f}, {yc:.1f}), Radius: {rx:.0f}, {ry:.0f}, {r:.0f}")
            # Store the displacement (center coordinates)
            displacements[:, img_idx_in_cube] = [xc, yc]

            # Slice the large mask centered on the detected solar disk
            ixc, iyc = int(round(xc)), int(round(yc))
            mask_center_x = int(large_mask.shape[1] / 2)
            mask_center_y = int(large_mask.shape[0] / 2)
            
            start_x = mask_center_x - ixc
            start_y = mask_center_y - iyc
            
            sliced_mask = large_mask[start_y : start_y + 480, start_x : start_x + 640]
            
            # Apply mask and store in the data cube
            masked_image = image_data * sliced_mask
            masked_image[masked_image == 0] = 0.001 # Avoid division by zero
            image_cube[:, :, img_idx_in_cube] = masked_image
            print(f"  Max: {masked_image.max():.1f}, Min: {masked_image.min():.1f}")            
            
            if Show:
                fig  = plt.figure()
                fig.set_size_inches(w=6,h=6)
                Bpos = [0.1,0.1,0.88,0.88]
                ax   = fig.add_subplot(1,1,1, position=Bpos)
                plt.imshow(masked_image,origin='lower',cmap='gray')
                uinput = input(" ")
                print(f"Continue, {uinput}")
                plt.close()

            img_idx_in_cube += 1

    # --- Calculate the first-pass flat using the collected data cube ---
    flat_pass = mk_flat(image_cube, threshold, displacements, iterations)
    flat_pass[flat_pass == 0] = 1.0

    return flat_pass

def apply(indir, telescope, iterations=10,Show=False):
    """
    This process obtains a Flat-Field mask for the 30THz camera using the
    Kuhn, Lin & Loranz (1991) method.

    It reads input images in FITS or FPF format, calculates the flat field
    in a two-pass process, and writes the resulting flat-field image in
    FITS and JPG formats to an 'final' subdirectory.

    Args:
        indir (str): Path to the directory containing the flat-field images.
        telescope (str): Telescope identifier, must be 'AR30T' or 'BR30T'.
    """
    # --- Validate Inputs ---
    if not os.path.isdir(indir):
        raise ValueError(f"*** ERROR: Input directory {indir} not found ***")
        
    if telescope not in ['AR30T', 'BR30T']:
        raise ValueEror(f"*** ERROR: Incorrect telescope specified ('{telescope}'). Must be 'AR30T' or 'BR30T'. ***")

    outdir = os.path.join(indir, 'final/')
    os.makedirs(outdir, exist_ok=True)
    print(f"Output will be saved in: {outdir}")

    # Record the start time
    start_time = time.perf_counter()

    # --- Create a large circular mask to eliminate image edges ---
    # This mask is larger than the images and will be sliced later
    print("Creating circular mask...")
    mask_height, mask_width = 1920, 2560
    cy, cx = mask_height / 2, mask_width / 2
    y, x = np.ogrid[:mask_height, :mask_width]
    dist_from_center = np.sqrt((x - cx)**2 + (y - cy)**2)
    large_mask = (dist_from_center < mask_radius).astype(np.float32)

    # --- Find input images and determine file type ---
    image_files = glob.glob(os.path.join(indir, '*.fits'))+ glob.glob(os.path.join(indir, '*.fts'))
    if not image_files:
        raise ValueError("*** ERROR: No '.fits' images found in {0:s} directory. ***".format(indir))
            
    image_files.sort()
    total_images = len(image_files)
    print(f"Found {total_images} images.")

    # --- Define the output filename based on the date from the first image ---
    first_img_path = image_files[0]
    output_name = ""
    date_obs_str = ""
    with fits.open(first_img_path) as hdul:
        header = hdul[0].header
        date_obs = header.get('DATE-OBS', datetime.now().strftime('%Y-%m-%d'))
        date_parts = date_obs.split('T')[0].split('-')
        output_name = f"{date_parts[0]}{date_parts[1]}{date_parts[2]}"
        date_obs_str = f"{date_parts[0]}-{date_parts[1]}-{date_parts[2]}"

    # =========================================================================
    # --- STEP 1: First Pass Flat-Field Calculation ---
    # =========================================================================
    print("\n--- Starting STEP 1: Initial Flat-Field Calculation ---")
    semilla = 1000
    flat_pass1 = flatPass(iterations,image_files,semilla,large_mask,Show=Show)
    pass1_time = time.perf_counter()
    elapsed_pass1_time = pass1_time - start_time
    print(f"Pass1 Execution time: {elapsed_pass1_time:.4f} seconds")


    # =========================================================================
    # --- STEP 2: Second Pass, using corrected images ---
    # =========================================================================
    print("\n--- Starting STEP 2: Refined Flat-Field Calculation ---")
    semilla = 1001
    flat_pass2 = flatPass(iterations,image_files,semilla,large_mask,Show=Show)    
    pass2_time = time.perf_counter()
    elapsed_pass2_time = pass2_time - pass1_time
    print(f"Pass2 Execution time: {elapsed_pass2_time:.4f} seconds")

    # The final flat is the product of both passes
    final_flat = flat_pass1 * flat_pass2

    # --- Save the final flat-field FITS file ---
    print("\n--- Finalizing and Saving Files ---")
    flatSave(telescope,final_flat,outdir,output_name, date_obs_str)

    # Record the end time
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Total Execution time: {elapsed_time:.4f} seconds")

    return

              
    
