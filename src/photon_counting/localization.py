import time
import numpy as np
import pandas as pd
import skimage.io as io
from scipy import ndimage
import scipy.optimize as opt
from pathlib import Path
from multiprocessing import Pool

def _two_d_gaussian(xdata_tuple, height, center_x, center_y, width, offset):
    """Flattened 2D Gaussian for curve fitting."""
    (y, x) = xdata_tuple 
    g = height * np.exp(-(((center_x - x) / width)**2 + ((center_y - y) / width)**2) / 2) + offset
    return g.ravel()

def _fit_single_blob(subarray, initial_guess, bounds):
    """Fits a 2D Gaussian to an isolated charge cloud."""
    y, x = np.indices(subarray.shape)
    try:
        popt, pcov = opt.curve_fit(
            _two_d_gaussian, (y, x), subarray.ravel(), 
            p0=initial_guess, bounds=bounds
        )
        perr = np.sqrt(np.diag(pcov))
        return popt, perr
    except RuntimeError:
        return None, None

def process_single_image(image_path, dark_avg_path, threshold, structure_size, save_dir):
    """
    Extracts x-ray interaction coordinates using STORM logic and saves to CSV.
    """
    image_path = Path(image_path)
    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)
    
    # Load and subtract dark average
    im = io.imread(image_path).astype(float) - io.imread(dark_avg_path).astype(float)
    
    # Threshold and apply morphological opening
    binary = (im > threshold).astype(np.uint8)
    structure = np.ones((structure_size, structure_size))
    opened = ndimage.binary_opening(binary, structure=structure)
    
    # Locate charge clouds
    labeled_map, num_features = ndimage.label(opened)
    results = []
    
    # Fit each identified cloud
    for i in range(1, num_features + 1):
        slice_y, slice_x = ndimage.find_objects(labeled_map == i)[0]
        subarray = im[slice_y, slice_x]
        
        # Dynamic bounds based on the sub-region size
        initial_guess = [subarray.max(), subarray.shape[1]/2, subarray.shape[0]/2, 2.0, subarray.min()]
        bounds = ([0, 0, 0, 0, -np.inf], [np.inf, subarray.shape[1], subarray.shape[0], np.inf, np.inf])
        
        popt, perr = _fit_single_blob(subarray, initial_guess, bounds)
        
        if popt is not None:
            height, cx, cy, width, offset = popt
            h_err, cx_err, cy_err, w_err, o_err = perr
            
            # Map local sub-array center back to global image coordinates
            global_cx = cx + slice_x.start
            global_cy = cy + slice_y.start
            area = 2 * np.pi * height * (width ** 2)
            
            results.append({
                'h': height, 'u(h)': h_err,
                'x': global_cx, 'u(x)': cx_err,
                'y': global_cy, 'u(y)': cy_err,
                'width': width, 'u(width)': w_err,
                'offset': offset, 'u(offset)': o_err,
                'area': area
            })
            
    # Save results cleanly using Pandas
    df = pd.DataFrame(results)
    save_path = save_dir / f"{image_path.stem}_storm.csv"
    df.to_csv(save_path, index=False)
    return save_path

def batch_process_storm(image_paths, dark_avg_path, threshold, structure_size, save_dir, workers=None):
    """Multiprocessing wrapper for bulk image analysis."""
    args = [(path, dark_avg_path, threshold, structure_size, save_dir) for path in image_paths]
    
    start_time = time.time()
    # Safely manages the CPU pool for both Windows and Unix
    with Pool(processes=workers) as pool:
        pool.starmap(process_single_image, args)
        
    print(f"Processed {len(image_paths)} images in {time.time() - start_time:.2f} seconds.")