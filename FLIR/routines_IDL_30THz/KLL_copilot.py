import numpy as np

def mk_kuhn_flat(images, thresh, niter=10, normalize=True):
    """
    Make a flat-field image using a set of spatially displaced data images
    (Kuhn, Lin & Loranz, PASP 103:1097, 1991), ported from IDL.
    
    Parameters
    ----------
    images : np.ndarray
        Array of shape (nrows, ncols, nimages). Temporal variations should be removed.
    thresh : float
        Threshold for selecting valid parts of images.
        * If thresh > 0: valid pixels satisfy images > thresh
        * If thresh < 0: valid pixels satisfy images < |thresh|
    offsets : np.ndarray
        Approximate positional offsets for all images, shape (2, nimages).
        offsets[0, k] = x (columns), offsets[1, k] = y (rows), in CCD pixels.
        Fractions allowed; they will be rounded to nearest integer (half-up).
        Offsets define the location of a common place in the images.
        If shifts were measured via correlation (IDL Astro lib), feed NEGATIVE of those shifts here.
    niter : int, optional
        Number of iterations (default 10).
    normalize : bool, optional
        If True, normalize flat so that the median of valid (>0) pixels is 1. Default True.
    return_debug : bool, optional
        If True, also return (karray, numarray, changes) for diagnostics.
    
    Returns
    -------
    flat : np.ndarray
        The flat-field image, shape (nrows, ncols). Zero where there is no data.
    (karray, numarray, changes) : tuple (optional)
        If return_debug is True, also returns:
          - karray: the accumulated log-ratio base image before iteration
          - numarray: counts per pixel of contributing pairs
          - changes: list of std(newflat - flat) per iteration (in log space)
    """

    offsets = np.ndarray((2,images.shape[2]))
    images = np.asarray(images, dtype=np.float64)
    assert images.ndim == 3, "images must be (nrows, ncols, nimages)"
    nrows, ncols, nimages = images.shape

    offsets = np.asarray(offsets, dtype=np.float64)
    assert offsets.shape == (2, nimages), "offsets must be shape (2, nimages)"

    # IDL-style half-up rounding (ties go away from zero, not banker's rounding)
    def idl_round(x):
        if x == 0:
            return 0
        return int(np.sign(x) * np.floor(abs(x) + 0.5))

    # Build masks based on threshold rule
    if thresh > 0:
        masks = images > thresh
    else:
        masks = images < abs(thresh)

    # (Optional) avoid log(<=0) by ensuring valid pixels are positive.
    # IDL doesn't add this, but if your data can have <= 0, uncomment:
    # masks &= (images > 0)

    # Precompute logs with zeros where invalid (mirrors IDL's ALOG + fill zeros)
    logs = np.zeros_like(images, dtype=np.float64)
    valid_pos = masks & (images > 0)  # only positive get log
    logs[valid_pos] = np.log(images[valid_pos])
    # pixels valid by threshold but <=0 remain zero in logs, effectively ignored via mask conjunction

    karray = np.zeros((nrows, ncols), dtype=np.float64)
    numarray = np.zeros((nrows, ncols), dtype=np.int32)

    # Helper to compute overlapping slices for two images given integer dx, dy
    # We want indices x,y in the *reference* (i) image such that the corresponding (x-dx,y-dy) is in bounds for image j.
    # dx, dy are pixel shifts: position_i - position_j (rounded)
    def overlap_slices(h, w, dx, dy):
        # x => columns (0..w-1), y => rows (0..h-1)
        # Compute x-range in i:
        x0_i = max(0, dx)
        x1_i = min(w, dx + w)
        y0_i = max(0, dy)
        y1_i = min(h, dy + h)
        if (x1_i <= x0_i) or (y1_i <= y0_i):
            return None  # no overlap
        # Corresponding region in j is shifted by (-dx, -dy)
        x0_j = x0_i - dx
        x1_j = x1_i - dx
        y0_j = y0_i - dy
        y1_j = y1_i - dy
        si = (slice(y0_i, y1_i), slice(x0_i, x1_i))
        sj = (slice(y0_j, y1_j), slice(x0_j, x1_j))
        return si, sj

    # === Build basic image: accumulate karray and numarray from all unordered pairs ===
    for i in range(nimages - 1):
        for j in range(i + 1, nimages):
            dx = idl_round(offsets[0, i] - offsets[0, j])
            dy = idl_round(offsets[1, i] - offsets[1, j])

            # Case 1: unshifted i vs shifted j
            res = overlap_slices(nrows, ncols, dx, dy)
            if res is not None:
                si, sj = res
                mask_ij = masks[si[0], si[1], i] & masks[sj[0], sj[1], j]
                if mask_ij.any():
                    diff = logs[si[0], si[1], i] - logs[sj[0], sj[1], j]
                    # accumulate only where mask true
                    karray[si] += np.where(mask_ij, diff, 0.0)
                    numarray[si] += mask_ij.astype(np.int32)

            # Case 2: unshifted j vs shifted i (reverse)
            dx_rev, dy_rev = -dx, -dy
            res = overlap_slices(nrows, ncols, dx_rev, dy_rev)
            if res is not None:
                si, sj = res
                mask_ji = masks[si[0], si[1], j] & masks[sj[0], sj[1], i]
                if mask_ji.any():
                    diff = logs[si[0], si[1], j] - logs[sj[0], sj[1], i]
                    karray[si] += np.where(mask_ji, diff, 0.0)
                    numarray[si] += mask_ji.astype(np.int32)

    valid_any = numarray > 0
    # Avoid division by zero
    karray_avg = np.zeros_like(karray)
    karray_avg[valid_any] = karray[valid_any] / numarray[valid_any]

    # === Iterate to a flat in log space ===
    flat_log = karray_avg.copy()
    changes = []

    for _ in range(max(1, int(niter))):
        correction = np.zeros_like(flat_log)

        for i in range(nimages - 1):
            for j in range(i + 1, nimages):
                dx = idl_round(offsets[0, i] - offsets[0, j])
                dy = idl_round(offsets[1, i] - offsets[1, j])

                # i (unshifted) receives from j (shifted)
                res = overlap_slices(nrows, ncols, dx, dy)
                if res is not None:
                    si, sj = res
                    mask_ij = masks[si[0], si[1], i] & masks[sj[0], sj[1], j]
                    if mask_ij.any():
                        # In IDL they zero-out where mask <=0
                        add = np.where(mask_ij, flat_log[sj], 0.0)
                        correction[si] += add

                # j (unshifted) receives from i (shifted)
                dx_rev, dy_rev = -dx, -dy
                res = overlap_slices(nrows, ncols, dx_rev, dy_rev)
                if res is not None:
                    si, sj = res
                    mask_ji = masks[si[0], si[1], j] & masks[sj[0], sj[1], i]
                    if mask_ji.any():
                        add = np.where(mask_ji, flat_log[sj], 0.0)
                        correction[si] += add

        new_flat_log = karray_avg.copy()
        new_flat_log[valid_any] += correction[valid_any] / np.maximum(numarray[valid_any], 1)

        change = np.std(new_flat_log - flat_log)  # same metric as IDL's STDEV
        changes.append(change)
        flat_log = new_flat_log

    # Exponentiate to leave log space; zero stays zero for invalid pixels
    flat = np.zeros_like(flat_log)
    flat[valid_any] = np.exp(flat_log[valid_any])

    # Optional normalization to median 1 across valid pixels
    if normalize:
        vals = flat[valid_any]
        if vals.size > 0:
            med = np.median(vals)
            if med > 0:
                flat[valid_any] = flat[valid_any] / med

    if return_debug:
        return flat, karray_avg, numarray, changes
