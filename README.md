# X-Ray STORM: Photon Counting and Super-Resolution Reconstruction

A Python data pipeline designed to surpass the physical spatial resolution limits of visible light sensors used in x-ray imaging. 

By applying principles derived from Stochastic Optical Reconstruction Microscopy (STORM), this package isolates individual x-ray charge clouds, fits their intensities to 2D Gaussians to find sub-pixel coordinates, and reconstructs the accumulated localizations into a high-resolution, zero-noise image.

## 🚀 Features
* **Cross-Platform:** Fully OS-agnostic path handling and multiprocessing.
* **Hardware Optimized:** Uses `multiprocessing.Pool` to distribute the heavy 2D Gaussian curve-fitting workload across all available CPU cores.
* **Vectorized Reconstruction:** Reconstructs massive image canvases (e.g., 4096 x 4096) near-instantly using Numpy vectorization.

## 📁 Repository Structure
* `/src/photon_counting`: Core algorithms for spatial localization and image synthesis.
* `/notebooks`: Contains `01_STORM_Localization_and_Reconstruction.ipynb`, providing a step-by-step tutorial of the processing pipeline.

## ⚙️ Installation & Usage

```bash
git clone [https://github.com/gitdwoc/photon_counting.git](https://github.com/gitdwoc/photon_counting.git)
cd photon_counting
pip install -r requirements.txt