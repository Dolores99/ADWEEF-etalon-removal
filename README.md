# ADWEEF: Automated Dual-Wavelength Raman Spectroscopy Etalon Elimination Framework

This repository contains MATLAB code for **ADWEEF**, a fully automated signal processing framework designed to suppress etalon fringes in dual-wavelength Raman spectroscopy, particularly for biomedical applications (e.g., skin analysis).

> **Status**: This work is currently under review at *Analytical Chemistry* (ACS).

---

## 🔍 Highlights

- Fully automated pipeline, no manual parameter tuning required
- Supports both multi-depth and single-point Raman datasets
- Etalon fringe detection via wavelet transform and statistical masking
- Two reconstruction approaches:
  - Frequency-domain denoising using Continuous Wavelet Transform (CWT)
  - Direct subtraction-based correction
- Optimized for dual-wavelength Raman systems (e.g., 671 nm and 785 nm excitation)

---

## 📁 Project Structure

```
ADWEEF/
│
├── data/                          # Sample datasets
│   ├── example_multiple_671.mat   # Multi-depth spectra (671 nm)
│   ├── example_multiple_785.mat   # Multi-depth spectra (785 nm)
│   ├── example_single_671.mat     # Single-depth spectrum (671 nm)
│   ├── example_single_785.mat     # Single-depth spectrum (785 nm)
│   ├── wl.mat                     # Shared wavelength axis
│   ├── wl_671.mat                 # Wavelengths for 671 nm
│   └── wl_785.mat                 # Wavelengths for 785 nm
│
├── functions/                     # Core processing functions
│   └── *.m                        # Includes wavelet filtering, boundary detection, ICA, visualization, etc.
│
├── main_ADWEEF_multiple.m         # Main script for multi-depth spectra (CWT-based method)
├── main_ADWEEF_single.m           # Main script for single-depth spectra (subtraction-based method)
├── EtalonMaskApp.m                # GUI for single-wavelength excitation spectra (subtraction-based method)
│
├── demo_ADWEEF_multiple.m         # Demonstration script for multi-depth input
├── demo_ADWEEF_single.m           # Demonstration script for single-depth input
│
├── LICENSE.txt                    # Custom academic license (non-commercial use only)
└── README.md                      # This documentation
```

---

## ▶️ Quick Start

```matlab
% Example: Run demo for multi-depth spectra
>> demo_ADWEEF_multiple

% Example: Run demo for single-depth signal
>> demo_ADWEEF_single
```

Make sure the required `.mat` files are in the `/data` folder.

---

## 📌 Dependencies
- MATLAB R2024a or newer
- Signal Processing Toolbox
- Wavelet Toolbox
- 200 colormap
- Bioinformatics Toolbox
- Curve Fitting Toolbox
- Image Processing Toolbox
- Time-reassigned Multisynchrosqueezing Transform

---

## 📄 License

This code is released for **academic and research use only**.  
**Commercial use is strictly prohibited** without prior written permission from the authors.

See [LICENSE.txt](./LICENSE.txt) for full terms.

---

## 📣 Citation

If you use this code or build upon it, please cite our work as:

> Wenyi Xu, Renzhe Bi, Yi Qi, Ruochong Zhang, Poongkulali Rajarahm, Alicia Yap Ann May, Dinish U.S, Qian Cheng, Steven Tien Guan Thng, Malini Olivo. _ADWEEF: Automated Dual-Wavelength Etalon Elimination Framework for Raman Spectroscopy_, under review at _Analytical Chemistry_, 2025.

Once the article is accepted, we will update this section with the DOI and full citation.

---

Thank you for using **ADWEEF**!
