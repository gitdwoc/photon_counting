import numpy as np
import pandas as pd
from pathlib import Path

def construct_super_res_image(csv_paths, original_shape, magnification_factor):
    """
    Constructs a high-resolution array from a list of STORM localization CSVs.
    """
    y_orig, x_orig = original_shape
    super_res_shape = (int(y_orig * magnification_factor), int(x_orig * magnification_factor))
    
    # Create the high-res canvas
    super_res_image = np.zeros(super_res_shape, dtype=np.uint16)
    
    for csv_path in csv_paths:
        try:
            df = pd.read_csv(csv_path)
            if df.empty:
                continue
                
            # Scale coordinates up by the magnification factor
            xs = np.round(df['x'] * magnification_factor).astype(int)
            ys = np.round(df['y'] * magnification_factor).astype(int)
            
            # Filter out coordinates that land outside the bounds
            valid_mask = (xs >= 0) & (xs < super_res_shape[1]) & (ys >= 0) & (ys < super_res_shape[0])
            xs = xs[valid_mask]
            ys = ys[valid_mask]
            
            # Vectorized addition: Adds +1 for every photon that lands on a specific pixel
            np.add.at(super_res_image, (ys, xs), 1)
            
        except Exception as e:
            print(f"Failed to process {csv_path.name}: {e}")
            
    return super_res_image