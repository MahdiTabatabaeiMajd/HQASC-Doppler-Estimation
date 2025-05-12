# Adaptive Doppler Blood Velocity Estimation

This repository provides MATLAB implementations for adaptive Doppler spectral estimation techniques introduced in the paper:

📄 **[Adaptive Doppler blood flow estimation in ultrasound with enhanced spectral resolution and contrast using limited observation windows](https://doi.org/10.1016/j.ultras.2025.107678)**  
**Authors**: Seyed Mohammad Mahdi Tabatabaei Majd, Leila Eslami, Babak Mohammadzadeh Asl  
**Journal**: *Ultrasonics*, Volume 154, 2025, 107678

---

## 🧠 Overview

Accurate Doppler spectral estimation is essential for evaluating blood flow, particularly during rapid systolic transitions. This toolbox implements **HQASC**, a high-resolution estimator that combines coherence-based post-filtering with eigen-decomposed Capon filtering to enhance contrast and frequency resolution, even with ultra short observation windows (e.g., **N=2**).

Key contributions:
- Implements multiple adaptive estimators: **Welch**, **Capon**, **Pr.Capon**, **MASC**, and **HQASC**
- Designed for **simulations and in vivo Doppler data**
- Supports **short-time Doppler** analysis under challenging flow conditions

---

## 📁 Folder Structure

| Folder       | Description |
|--------------|-------------|
| `datasets/`  | Contains sample or user-supplied Doppler `.mat` files (must contain variable `data`) |
| `scripts/`   | Main scripts for simulation and execution |
| `functions/` | Estimation algorithms and plotting utilities |
| `README.md`  | Overview, usage, and citation |
| `LICENSE`    | Licensing terms (MIT) |

---

## 🔬 Methods Implemented

### 1. Welch Estimator
- Baseline non-adaptive periodogram using rectangular and Hamming windows.

### 2. Capon (Minimum Variance)
- Adaptive filtering to suppress interference and noise, enhancing frequency selectivity.

### 3. Projected Capon (Pr.Capon)
- Eigen-based projection to reduce sidelobe levels and improve estimation.

### 4. MASC (Modified Amplitude Spectrum Capon)  
- Tabatabaei Majd & Mohammadzadeh Asl, *IEEE TUFFC*, 2021  
  [DOI:10.1109/TUFFC.2020.3044774](https://doi.org/10.1109/TUFFC.2020.3044774)

### 5. HQASC (High-Quality Adaptive Spectral Coherence)
- Tabatabaei Majd et al., *Ultrasonics*, 2025  
  [DOI:10.1016/j.ultras.2025.107678](https://doi.org/10.1016/j.ultras.2025.107678)

---

## ▶️ How to Run

### A. Apply Estimators to Real or Simulated Data

```matlab
cd scripts
run_adaptive_doppler_estimation
```

You will be prompted to:
- Select a `.mat` file from the `datasets/` folder (must contain `data`)
- Choose a folder to save results

### B. Monte Carlo Simulation: MSE vs. SNR

```matlab
cd scripts
test_velocity_estimation
```

Edit parameters like `SNR`, `MC`, and `Ns` in the script to match your configuration.

---

## 🛠 Dependencies

- MATLAB R2019b or later
- No external toolboxes required

---

## 📊 Output

- Spectral matrices (`allSpectra_*.mat`)
- Time–velocity spectrogram plots
- Velocity estimation error curves from simulations

---

## 🧾 Citation

If you use these codes, please cite:

```bibtex
@article{tabatabaei2025hqasc,
  title     = {Adaptive Doppler blood flow estimation in ultrasound with enhanced spectral resolution and contrast using limited observation windows},
  author    = {Tabatabaei Majd, Seyed Mohammad Mahdi and Eslami, Leila and Mohammadzadeh Asl, Babak},
  journal   = {Ultrasonics},
  volume    = {154},
  pages     = {107678},
  year      = {2025},
  doi       = {10.1016/j.ultras.2025.107678}
}
```

---

## 📬 Contact

For questions or feedback, please contact:

- **S.M.M. Tabatabaei Majd** — [tabatabaei.majd@gmail.com](mailto:tabatabaei.majd@gmail.com)  
- **Leila Eslami** — [leila.eslami2020@gmail.com](mailto:leila.eslami2020@gmail.com)
