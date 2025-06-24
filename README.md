# ADWEEF: Automated Dual-Wavelength Etalon Elimination Framework

A fully automated algorithm for removing Etalon fringes in dual-wavelength Raman spectroscopy systems, specifically tailored for BI-CCD-based detectors. This method enables robust fringe suppression without the need for manual parameter tuning or prior modeling, and is highly suitable for physiological and biomedical Raman applications such as skin analysis.

---

## 🔍 About

This repository provides the official implementation of **ADWEEF**, the algorithm proposed in our manuscript currently **under review at _Analytical Chemistry_ (ACS)**.

The method is designed for use with dual-wavelength Raman systems and leverages the correlation of fringe patterns across wavelengths to remove Etalon noise in an unsupervised and automated manner.

---

## 🛠 Files

- `main.m` – Core MATLAB implementation of the ADWEEF algorithm  
- `example_data.mat` – Sample synthetic dual-wavelength Raman data  
- `demo_ADWEEF.m` – Example usage script to run the algorithm  
- `LICENSE.txt` – License (academic use only)  
- `README.md` – Project documentation  

---

## 🚀 Getting Started

1. Clone or download this repository  
2. Open `demo_ADWEEF.m` in MATLAB  
3. Run the script to see the algorithm applied on example data  
4. Replace `example_data.mat` with your own data to test the method

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
