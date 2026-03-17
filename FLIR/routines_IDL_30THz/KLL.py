import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndi # For potential future smoothing, though not in original
import pdb

def mk_flat(images, thresh, offsets=None, niter=10, Show=False):
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

