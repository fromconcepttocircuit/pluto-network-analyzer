# ADALM Pluto SDR Network Analyzer (0.1–3 GHz)

**Transform your ADALM Pluto SDR into a budget-friendly 0.1–3 GHz network analyzer using a \$15 RF bridge!**

---

## 📡 See It In Action (Video Demo)

🎬 **Watch the full video demo and tutorial on YouTube**:  

[![Watch the video](https://img.youtube.com/vi/IFquLM-xF30/0.jpg)](https://www.youtube.com/watch?v=IFquLM-xF30)




## 📋 Overview

This interactive Python application turns your ADALM Pluto SDR into a powerful network analyzer capable of measuring S21 (insertion loss/gain) and S11 (return loss) over a 0.1–3 GHz span. Key features include:

* **FFT-based FIR lock-in** for tone detection
* **Interactive S21 & S11 calibration** workflows (LOAD & OPEN)
* **Skip-mask** to reject invalid S11 data
* **Moving-average smoothing** toggle
* **Interactive markers** (left-click to add, right-click to clear)
* **Maximum S11 clamp** at 0 dB

Ideal for RF enthusiasts, hobbyists, and professionals seeking a low-cost network analyzer alternative.

---

## ⚙️ Features

* **Wideband sweep**: 0.1–3 GHz, configurable start/stop frequencies & step count
* **Calibration workflows**: Perform LOAD and OPEN calibration at the press of a button
* **Live plotting**: Dual-panel Matplotlib GUI with S21 (top) & S11 (bottom)
* **Smoothing & raw toggle**: Show raw data points or smoothed curves
* **Markers**: Click to annotate dB & frequency values
* **Excel export**: S11 calibration writes `s11_calibration.xlsx` automatically

---

## 🛠️ Requirements

* Python 3.8+
* [libiio](https://github.com/analogdevicesinc/libiio) & [pyadi-iio](https://github.com/analogdevicesinc/pyadi-iio)
* NumPy
* SciPy
* Matplotlib
* Pandas
* PyQt6

Install dependencies via:

```bash
pip install numpy scipy matplotlib pandas pyqt6 pyadi-iio
```

---

## 📥 Installation

1. **Clone the repo**:

   ```bash
   git clone https://github.com/fromconcepttocircuit/pluto-network-analyzer.git
   cd pluto-network-analyzer
   ```
2. **Install Python dependencies** (see requirements above).
3. **Connect** your ADALM Pluto SDR to your network (default URI: `ip:192.168.2.1`).

---

## 🚀 Usage & Calibration

1. **Run the GUI**:

   ```bash
   python network_analyzer.py
   ```
2. **Set sweep span**: Enter start/stop frequencies (MHz) and number of steps, then click **Apply**.
3. **Calibrate S21**: Click **Cal S21**, wait for the progress dialog to finish. The application will automatically save `cal_s21.npz` in the working directory and load it into the session.
4. **Calibrate S11**: Click **Cal S11**. When prompted, connect a 50 Ω load and click OK, then connect an OPEN termination and click OK. The tool will save `cal_s11_load.npz`, `cal_s11_open.npz`, and export `s11_calibration.xlsx` automatically, then load the new calibration.
5. **Auto-load on startup**: On next launch, if the `.npz` and `.xlsx` files exist, the application will load them automatically. You can also reload them at any time using **Load S21** and **Load S11** buttons.
6. **Plot & analyze**: Toggle smoothing/raw display, and click on the plots to add markers showing dB & frequency.

> **Note:** Calibration files are generated and retrieved at runtime. Recalibrate any time to refresh measurement accuracy.

## 📂 File Structure

```
├── network_analyzer.py    # Main application script
└── README.md             # This document
```

---

## 🎥 YouTube Channel

For more tutorials and deep dives into RF & SDR projects, visit my YouTube channel:

[From Concept to Circuit](http://www.youtube.com/@fromconcepttocircuit)

---


## 📄 License

This project is released under the MIT License. 

---

Happy hacking! 🚀
