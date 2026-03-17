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

class Flat():

    def __init__(self,indir='./falts/',telescope='AR30T',itoll=1.0E-03):

        self.Circle = Circle()
        self.MetaData = {'indir':indir,
                         'telescope':telescope.upper(),
                         'itoll':itoll,
                         'threshold':0.1,
                         'img_count_for_cube':100,
                         'mask_radius':330}
        self.image_files = []
        
        return

    def apply(self):
    """
    This process obtains a Flat-Field mask for the 30THz camera using the
    Kuhn, Lin & Loranz (1991) method.

    It reads input images in FITS format, creates the flat field
    in a two-pass process, and writes the resulting flat-field image in
    FITS format.

    Args:
        indir (str): Path to the directory containing the flat-field images.
        telescope (str): Telescope identifier, must be 'AR30T' or 'BR30T'.
    """
    # --- Validate Inputs ---
        if not os.path.isdir(self.MetaData['indir']):
            raise ValueError(f"*** ERROR: Input directory {self.MetaData['indir']} not found ***")
        
        if self.MetaData['telescope'] not in ['AR30T', 'BR30T']:
            raise ValueEror(f"*** ERROR: Incorrect telescope specified ({self.MetaData['telescope']}). AR30T | BR30T only allowed values. ***")

    # --- Create a large circular mask to eliminate image edges ---
    # This mask is larger than the images and will be sliced later
        print("Creating circular mask...")
        mask_height, mask_width = 1920, 2560
        cy, cx = mask_height / 2, mask_width / 2
        y, x = np.ogrid[:mask_height, :mask_width]
        dist_from_center = np.sqrt((x - cx)**2 + (y - cy)**2)
        self.large_mask = (dist_from_center < mask_radius).astype(np.float32)

    # --- Find input images and determine file type ---
        image_files = glob.glob(os.path.join(self.MetaData['indir'], '*.fits'))
        if not image_files:
            raise ValueError(f"*** ERROR: No '.fits' images found in {self.MetaData['indir']} directory. ***")
            
        image_files.sort()
        total_images = len(image_files)

    # --- Define the output filename based on the date from the first image ---
        first_img_file = image_files[0]
        output_name = ""
        date_obs_str = ""
        with fits.open(first_img_file) as hdul:
            header = hdul[0].header
            self.date_obs = header.get('DATE-OBS', datetime.now().strftime('%Y-%m-%d'))
            self.date_parts = date_obs.split('T')[0].split('-')
            self.output_name = f"{date_parts[0]}{date_parts[1]}{date_parts[2]}"
            self.date_obs_str = f"{date_parts[0]}-{date_parts[1]}-{date_parts[2]}"

    # =========================================================================
    # --- STEP 1: First Pass Flat-Field Calculation ---
    # =========================================================================
        print("\n--- Starting STEP 1: Initial Flat-Field Calculation ---")
        self.semilla1 = 1000
        flat_pass1 = flatPass(self,1,image_files)

    # =========================================================================
    # --- STEP 2: Second Pass, using corrected images ---
    # =========================================================================
        print("\n--- Starting STEP 2: Refined Flat-Field Calculation ---")
        self.semilla2 = 1001
        flat_pass2 = flatPass(self,2,image_files)    

    # The final flat is the product of both passes
        self.flat = self.flat1 * self.flat2

    # --- Save the final flat-field FITS file ---
        print("\n--- Finalizing and Saving Files ---")
        flatSave(self)

        return

    def flatPass(self,pass,image_files):

     # Initialize data containers
        image_cube = np.zeros((480, 640, img_count_for_cube), dtype=np.float32)
        displacements = np.zeros((2, img_count_for_cube), dtype=np.float32)
    
    # Generate a random sequence of images to process
        total_images    = len(image_files)
        if pass = 1:
            semilla = self.semilla1
        else:
            semilla = self.semilla2
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
                self.image_files.append(current_file)
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
            limb_image =
            roberts(disk_mask)
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
                par,_ = self.Circle.fit(x_coords,y_coords)
                xc = par[0]
                yc = par[1]
                r = par[2]
            
            except Exception:
                print(f"  Failed to fit the limb")
                continue # Skip if fitting fails

        # --- Validate the found center and add image to the cube ---        
        if xc > 100 and xc < 539 and yc > 100 and yc < 379:

            #print(f"Radius: {r:.0f}  Center: ({xc:.1f}, {yc:.1f})")
            
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
            #print(f"  Max: {masked_image.max():.1f}, Min: {masked_image.min():.1f}")

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
    flat_pass = mk_flat(image_cube, threshold, displacements, iterations,Show=Show)
    flat_pass[flat_pass == 0] = 1.0
    
    return flat_pass
    
    def kll(self):
    """
    Make a flat field image using a set of spatially displaced data images.
    Implements the method of Kuhn, Lin, and Loranz, Pub. Astron.
    Soc. Pacific, 103:1097, 1991.

    Parameters
    ----------
    images : numpy.ndarray
        A 3D cube of images (height, width, n_images).
        Temporal variations must be removed prior to using this function.
    thresh : float
        Threshold for selecting valid parts of images.
        If thresh > 0, valid data > thresh.
        If thresh < 0, valid data < |thresh|.
    offsets : numpy.ndarray, optional
        Vector of approximate positional offsets of the images,
        giving x,y shifts for all images. Dimensioned (2, n_images).
        Offset is defined to be the location of a common place in the
        images, in CCD pixels. Offsets should be given to
        fractional pixels, to get shifts between images correct.
        If the IDL Astronomy correlation procedures are used
        to measure the shifts between images, the offsets
        to this routine should be the negatives of the shifts.
        If None, all offsets are assumed to be (0,0).
    niter : int, optional
        Number of iterations to run. Default = 10.

    Returns
    -------
    numpy.ndarray
        The flat field image, normalized. Zero values => no data.

    Notes
    -----
    This Python implementation directly translates the IDL logic.
    IDL's array indexing (e.g., `(dx>0):ncols-1+(0<dx)`) and `WHERE`
    function require careful translation to NumPy slicing and boolean indexing.
    The fixed-pixel shifting approximation for `dx` and `dy` from
    IDL `FIX(dx/ABS(dx)) * FIX(ABS(dx)+0.5)` is replicated.
    """

    # IDL: asize = SIZE(images)
    # IDL array dimensions are (ncols, nrows, nimages)
    # NumPy array dimensions are typically (nrows, ncols, nimages)
    # We need to be careful with transposing or interpreting dimensions.
    # Assuming `images` is (height, width, n_images) in Python.
        if images.ndim != 3:
            raise ValueError("Input 'images' must be a 3D cube (height, width, n_images).")

        nrows, ncols, nimages = images.shape

    # Ensure offsets are provided or create zero offsets
    if offsets is None:
        offsets = np.zeros((2, nimages), dtype=float)
    elif offsets.shape != (2, nimages):
        raise ValueError(f"Offsets must be of shape (2, {nimages}) or None.")

    # Construct the masks array
    # IDL: masks = REPLICATE( BYTE(0), ncols*nrows, nimages)
    # In Python, we can create a boolean mask directly with the image dimensions
    print("Building masks...")
    masks = np.zeros_like(images, dtype=bool) # (nrows, ncols, nimages)

    for i in range(nimages):
        img_slice = images[:, :, i]
        if thresh > 0:
            masks[:, :, i] = (img_slice > thresh)
        else:
            masks[:, :, i] = (img_slice < -thresh)

    # Initialize karray and numarray
    # IDL: karray = FLTARR(ncols, nrows), numarray = INTARR(ncols, nrows)
    # Python: (nrows, ncols)
    karray = np.zeros((nrows, ncols), dtype=float)
    numarray = np.zeros((nrows, ncols), dtype=int)

    print('Building basic image (K array).')
    # Loop over image pairs (i, j)
    for i in range(nimages - 1):
        # IDL: tmp = images(*,*,i), valid = WHERE( masks(*,*,i) GT 0 )
        # IDL: logi = REPLICATE(0.,ncols,nrows), logi(valid) = ALOG(tmp(valid))
        logi = np.zeros((nrows, ncols), dtype=float)
        valid_i = masks[:, :, i]
        logi[valid_i] = np.log(images[:, :, i][valid_i])
        
        for j in range(i + 1, nimages):
            # Calculate integer pixel shifts based on IDL logic
            # IDL: dx = offsets(0,i) - offsets(0,j)
            # IF( dx NE 0. ) THEN dx = FIX(dx/ABS(dx)) * FIX(ABS(dx)+0.5)
            # This is IDL's way of getting a +1, -1, or 0 integer shift.
            # Python equivalent:
            dx_orig = offsets[0, i] - offsets[0, j]
            dy_orig = offsets[1, i] - offsets[1, j]

            # Replicate IDL's FIX(x/abs(x)) * FIX(abs(x)+0.5) behavior
            # This effectively rounds to the nearest integer for dx, dy
            # and takes the sign.
            dx = int(np.sign(dx_orig) * np.round(np.abs(dx_orig)))
            dy = int(np.sign(dy_orig) * np.round(np.abs(dy_orig)))

            logj = np.zeros((nrows, ncols), dtype=float)
            valid_j = masks[:, :, j]
            logj[valid_j] = np.log(images[:, :, j][valid_j])
            
            # --- Process (Image i unshifted, Image j shifted) ---
            # IDL slicing: (dx>0):ncols-1+(0<dx) for X-dim
            # This is a bit tricky to translate directly because IDL's
            # (A>0):B expands to min(A,B):max(A,B) and includes the endpoints.
            # In Python, we need to handle start/end indices carefully for slicing.
            # It defines the overlapping region.

            # Example: IDL (dx>0):ncols-1+(0<dx)
            # if dx > 0, it's dx:ncols-1
            # if dx <= 0, it's 0:ncols-1+dx (which means 0:ncols-1 for dx=0, 0:ncols-1-|dx| for dx<0)
            # This is equivalent to: start = max(0, dx), end = min(ncols, ncols + dx)

            # Define slices for image i (unshifted)
            # Slice for Image i (unshifted in the current iteration)
            slice_i_y = slice(max(0, dy), min(nrows, nrows + dy))
            slice_i_x = slice(max(0, dx), min(ncols, ncols + dx))

            # Define slices for image j (shifted to align with i)
            # The shift for J is -dx, -dy relative to its own coordinate system
            slice_j_y = slice(max(0, -dy), min(nrows, nrows - dy))
            slice_j_x = slice(max(0, -dx), min(ncols, ncols - dx))

            # Extract log values and masks for overlapping regions
            sum_val = np.copy(logi[slice_i_y, slice_i_x])
            
            mask_i_overlap = masks[slice_i_y, slice_i_x, i]
            sum_val -= logj[slice_j_y, slice_j_x]
            mask_j_overlap = np.copy(masks[slice_j_y, slice_j_x, j])

            # Combine masks: only valid if data in *both* images is valid in the overlap
            combined_mask_ij = mask_i_overlap & mask_j_overlap

            # Apply mask to sum_val
            sum_val_masked = sum_val * combined_mask_ij.astype(float)

            # Accumulate K array and N array
            # Karray and Numarray indices are for the "unshifted i" frame of reference
            karray[slice_i_y, slice_i_x] += sum_val_masked
            numarray[slice_i_y, slice_i_x] += combined_mask_ij.astype(int)

            # --- Process (Image j unshifted, Image i shifted) ---
            # Now, for the second pass, we swap roles, so image j is unshifted
            # and image i is shifted by (-dx, -dy) relative to j's frame.

            dx_prime = -dx # The shift for i relative to j's frame
            dy_prime = -dy

            # Slice for Image j (unshifted in this sub-iteration)
            slice_j_prime_y = slice(max(0, dy_prime), min(nrows, nrows + dy_prime))
            slice_j_prime_x = slice(max(0, dx_prime), min(ncols, ncols + dx_prime))

            # Slice for Image i (shifted to align with j)
            slice_i_prime_y = slice(max(0, -dy_prime), min(nrows, nrows - dy_prime))
            slice_i_prime_x = slice(max(0, -dx_prime), min(ncols, ncols - dx_prime))
            
            sum_val_prime = np.copy(logj[slice_j_prime_y, slice_j_prime_x])
            mask_j_overlap_prime = masks[slice_j_prime_y, slice_j_prime_x, j]
            
            sum_val_prime -= logi[slice_i_prime_y, slice_i_prime_x]
            mask_i_overlap_prime = masks[slice_i_prime_y, slice_i_prime_x, i]
            
            combined_mask_ji = mask_j_overlap_prime & mask_i_overlap_prime
            sum_val_prime_masked = sum_val_prime * combined_mask_ji #* combined_mask_ji.astype(float)
            karray[slice_j_prime_y, slice_j_prime_x] += sum_val_prime_masked
            numarray[slice_j_prime_y, slice_j_prime_x] += combined_mask_ji.astype(int)
            
    # Finalize K array
    valid = numarray > 0
    karray[valid] = karray[valid] / numarray[valid]

    if Show:
        # Display initial karray (IDL: TVSCL,karray)
        # This is optional for debugging/visualization
        plt.figure(figsize=(8, 6))
        plt.imshow(karray, origin='lower', cmap='gray')
        plt.title("Initial K array (Before Iterations)")
        plt.colorbar(label="Log Flat Field Difference")
        uinput = input(" ")
        print(f"Continue, {uinput}")
        plt.show()

    # --- Now iterate to a flat ---
    flat = np.copy(karray) # IDL: flat = karray
    print('Starting iterations.')
    print('Iteration, Change (STDEV)')

    for iter_num in range(1, niter + 1):
        correction = np.zeros_like(flat, dtype=float)

        # Sum up the corrections
        for i in range(nimages - 1):
            for j in range(i + 1, nimages):
                dx_orig = offsets[0, i] - offsets[0, j]
                dy_orig = offsets[1, i] - offsets[1, j]

                dx = int(np.sign(dx_orig) * np.round(np.abs(dx_orig)))
                dy = int(np.sign(dy_orig) * np.round(np.abs(dy_orig)))

                # Unshifted I, Shifted J
                slice_i_y = slice(max(0, dy), min(nrows, nrows + dy))
                slice_i_x = slice(max(0, dx), min(ncols, ncols + dx))

                slice_j_y = slice(max(0, -dy), min(nrows, nrows - dy))
                slice_j_x = slice(max(0, -dx), min(ncols, ncols - dx))

                mask_i_overlap = masks[slice_i_y, slice_i_x, i]
                mask_j_overlap = masks[slice_j_y, slice_j_x, j]
                combined_mask_ij = mask_i_overlap & mask_j_overlap

                # Sum for correction is from the *shifted* image's flat field value
                # IDL: sum = flat( (-dx>0):ncols-1-(dx>0), (-dy>0):nrows-1-(dy>0) )
                # This corresponds to the region of `flat` that maps to the shifted J image
                sum_shifted_j_flat = flat[slice_j_y, slice_j_x]
                sum_shifted_j_flat[~combined_mask_ij] = 0.0 # Zero out invalid pixels

                correction[slice_i_y, slice_i_x] += sum_shifted_j_flat

                # Unshifted J, Shifted I
                dx_prime = -dx
                dy_prime = -dy

                slice_j_prime_y = slice(max(0, dy_prime), min(nrows, nrows + dy_prime))
                slice_j_prime_x = slice(max(0, dx_prime), min(ncols, ncols + dx_prime))

                slice_i_prime_y = slice(max(0, -dy_prime), min(nrows, nrows - dy_prime))
                slice_i_prime_x = slice(max(0, -dx_prime), min(ncols, ncols - dx_prime))

                mask_j_overlap_prime = masks[slice_j_prime_y, slice_j_prime_x, j]
                mask_i_overlap_prime = masks[slice_i_prime_y, slice_i_prime_x, i]
                combined_mask_ji = mask_j_overlap_prime & mask_i_overlap_prime

                # Sum for correction is from the *shifted* image's flat field value
                # This corresponds to the region of `flat` that maps to the shifted I image
                sum_shifted_i_flat = flat[slice_i_prime_y, slice_i_prime_x]
                sum_shifted_i_flat[~combined_mask_ji] = 0.0

                correction[slice_j_prime_y, slice_j_prime_x] += sum_shifted_i_flat


        newflat = np.copy(karray) # IDL: newflat = karray
        # IDL: newflat(valid) = newflat(valid) + correction(valid) / numarray(valid)
        newflat[valid] += correction[valid] / numarray[valid]

        change = np.std(newflat - flat) # IDL: STDEV(newflat - flat)
        flat = newflat

        if Show:
        # Display intermediate flat (IDL: TVSCL,flat)
            plt.figure(figsize=(8, 6))
            plt.imshow(flat, origin='lower', cmap='gray')
            plt.title(f"Flat field after iteration {iter_num}")
            plt.colorbar(label="Log Flat Field")
            plt.show() # Uncomment to see each iteration's progress
            uinput = input(" ")
            print(f"Continue, {uinput}")
            plt.show()

        print(f'{iter_num}, {change:.6f}')

    # Final step: Convert back from log space
    # IDL: flat(valid) = EXP(flat(valid))
    final_flat = np.zeros_like(flat)
    final_flat[valid] = np.exp(flat[valid])

    # Zero values => no data, as per IDL
    final_flat[~valid] = 0.0

    return final_flat

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



              

class Circle():

    import numpy as np
    from scipy.optimize import curve_fit
    import pdb

""" ------------------------------------------------------------------

CIRC_FIT

   A series of procedures to fit a circle from a data set.
   It derives from the old and good SST's 'circle_fit.pro'
   The sequence of use is as follows. Assuming you have an image
   M[x,y]

   limb=circ_fit.Get_Limb(M)
   par,cov = circ_fit.fit(limb)

   par is a 3 element list with the results
   Circle Center: Col = par[0], Row = par[1]
   Circle Radius: par[2]

   cov is covariance matrix.

Author:  @Guiguesp
         2018-02-07 after too many mistakes >-(

------------------------------------------------"""

    def __init__(self,x,y):

        self.data = {'x':x,'y':y}
        self.fit()
        return

    def circle_func(self,phi,x0,y0,r):
    """
    CIRCLE

       The fitting function used by curve_fit

    """

        b   = x0 * np.cos(phi) + y0 * np.sin(phi)
        wur = np.sqrt( b*b - x0*x0 - y0*y0 + r*r)
        d   = b + wur

        return d

    def fit(self):

        ierase   = np.array([])
        Niter    = 0
        Nmin     = 50
        MaxNiter = 5
        Continue = True
        x = self.data['x']
        y = self.data['y']

        while Continue:

        # Prepare data for polar fitting
        # First guess of circle center
            if (ierase.shape[0]) > 1 :
                x = np.delete(x,ierase)
                y = np.delete(y,ierase)

            Col_0     = x.mean()
            Row_0     = y.mean()

        # Distance of limb to disc center
            dx = x-Col_0
            dy = y-Row_0
            d  = np.sqrt( dx**2 + dy**2)

        # First guess for the circle radius
            R_0    = d.mean()
            phi = np.arctan2(dy,dx)
            so  = np.argsort(phi)
            phi = phi[so]
            x   = x[so]
            y   = y[so]
            d   = d[so]
            
        # First guess of parameters. Note that, since we subtracted the circle center
        # to the limb, the center first guess is (0,0)
            par0     = (0,0,R_0)
            par, cov = curve_fit(self.circle_func, phi,d,p0=par0)
            dfit     = self.circle_func(phi,par[0],par[1],par[2])
            fsigma   = np.std(d-dfit)
            psigma   = np.sqrt(np.diag(cov))
            ierase   = np.ravel(np.where(np.abs(d-dfit) > fsigma))
            Niter    += 1

            if (len(x) == (x.shape[0])) | (Niter > MaxNiter) | (x.shape[0] < 20)  :
                Continue = False

    # Return the disc center in the matrix coordinates adding what we subtracted
        par[0]   += Col_0
        par[1]   += Row_0
        self.fit = {'x0': par[0], 'y0':par[1],'r':par[2],'cov':cov}
        return
