#!/usr/bin/env python3
"""
PlutoSDR VNA Pro – interactive edition
======================================
• FFT-based FIR lock-in
• Calibration + skip-mask for S11
• Moving-average smoothing **toggle**
• Add/remove markers with left / right mouse clicks
• Clamp S11 to a maximum of 0 dB
"""

# ───────────── imports ─────────────
import sys, os, time, traceback
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("qtagg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import (
    NavigationToolbar2QT, FigureCanvasQTAgg as FigureCanvas)
import adi
from scipy.signal import kaiserord, firwin, lfilter, fftconvolve, oaconvolve
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QLabel,
    QVBoxLayout, QHBoxLayout, QPushButton, QLineEdit,
    QFrame, QProgressDialog, QMessageBox, QCheckBox)
from PyQt6.QtCore import Qt, QThread, pyqtSignal

# ───────────── configuration constants ─────────────
FILTER_MODE, FFT_METHOD = 'fft', 'fftconvolve'
FILT_RIPPLE_DB, FILT_CUTOFF_HZ, FILT_TRANS_WIDTH_HZ = 70, 500, 100
DEFAULT_STEPS, CAL_POINTS = 500, 500
SMOOTH_WIN21, SMOOTH_WIN11 = 7, 5
S11_DENOM_THRESHOLD, SKIP_MARGIN = 3e-2, 3
SDR_URI = "ip:192.168.2.1"
FS, NUM_S, TONE = 8e6, 2**18, 543e3
RF_F=4e6
DWELL, CLR_READS, EPS = 0.1, 1, 1e-15
MIN_FREQ, MAX_FREQ, SWEEP_INIT_DELAY = 0.3e9, 3e9, 2.0

# ───────────── hardware initialisation ─────────────
sdr = adi.ad9361(uri=SDR_URI)
sdr.sample_rate             = int(FS)
sdr.tx_enabled_channels     = [0]
sdr.rx_enabled_channels     = [0, 1]
sdr.rx_rf_bandwidth         = int(RF_F)
sdr.tx_rf_bandwidth         = int(RF_F)
sdr.rx_buffer_size          = NUM_S
sdr.tx_buffer_size          = NUM_S
sdr.gain_control_mode_chan0 = "manual"
sdr.gain_control_mode_chan1 = "manual"
sdr.tx_cyclic_buffer        = True
sdr.tx_hardwaregain_chan0   = -10
sdr.rx_hardwaregain_chan0   = 50
sdr.rx_hardwaregain_chan1   = 50

_t = np.arange(NUM_S) / FS
sdr.tx((np.exp(2j * np.pi * TONE * _t) * (2**14)).astype(np.complex64))

# ───────────── FIR filter for lock-in ─────────────
nyq     = FS / 2
N, beta = kaiserord(FILT_RIPPLE_DB, FILT_TRANS_WIDTH_HZ/nyq)
b_fir   = firwin(N, FILT_CUTOFF_HZ/nyq, window=('kaiser', beta))

def apply_filter(x):
    if FILTER_MODE == 'direct':
        return lfilter(b_fir, 1.0, x)
    func = fftconvolve if FFT_METHOD == 'fftconvolve' else oaconvolve
    return func(x, b_fir, mode='same')

def lockin(buf: np.ndarray) -> float:
    if np.allclose(buf, 0):
        return 0.0
    t = np.arange(len(buf)) / FS
    y = apply_filter(buf * np.exp(-2j * np.pi * TONE * t))
    y = y[N//2:]  # discard FIR transient
    return np.abs(y).mean()

def to_dB(x):
    return 20 * np.log10(np.maximum(x, EPS))

def smooth_trace(y, k):
    if k <= 1 or y.size < k:
        return y
    mask = np.isfinite(y).astype(float)
    win = np.ones(k)
    num = np.convolve(np.nan_to_num(y), win, 'same')
    den = np.convolve(mask, win, 'same')
    out = np.full_like(y, np.nan)
    good = den > 0
    out[good] = num[good] / den[good]
    return out

# ───────────── worker thread ─────────────
class SweepThread(QThread):
    update  = pyqtSignal(float, float, float)
    started = pyqtSignal()
    error   = pyqtSignal(str)

    def __init__(self, dev):
        super().__init__()
        self.dev       = dev
        self.f0, self.f1, self.n = MIN_FREQ, MAX_FREQ, DEFAULT_STEPS
        self.stop_flag = False
        self.cal21     = None
        self.cal11L    = None
        self.cal11O    = None
        self.skip_mask = None

    def set_span(self, f0, f1, n):
        self.f0, self.f1, self.n = f0, f1, n

    def stop(self):
        self.stop_flag = True

    def load_cal21(self, d):
        self.cal21 = d

    def load_cal11(self, L, O):
        self.cal11L = L
        self.cal11O = O

    def _mask(self, freqs):
        if self.cal11L is None or self.cal11O is None:
            return np.zeros_like(freqs, bool)
        L = np.interp(freqs, self.cal11L['freqs'], self.cal11L['linear'])
        O = np.interp(freqs, self.cal11O['freqs'], self.cal11O['linear'])
        bad = np.abs(O - L) < S11_DENOM_THRESHOLD
        if SKIP_MARGIN:
            bad = np.convolve(bad.astype(int),
                              np.ones(2*SKIP_MARGIN+1, int),
                              'same') > 0
        return bad

    def _safe_rx(self):
        try:
            return self.dev.rx()
        except Exception as e:
            self.error.emit(f"SDR RX error: {e}")
            self.stop_flag = True
            raise

    def run(self):
        time.sleep(SWEEP_INIT_DELAY)
        while not self.stop_flag:
            freqs = np.linspace(self.f0, self.f1, self.n)
            self.skip_mask = self._mask(freqs)
            self.started.emit()
            for i, f in enumerate(freqs):
                if self.stop_flag:
                    break
                NUM_R = 4 if f < 1e9 else 1
                try:
                    self.dev.tx_lo = self.dev.rx_lo = int(f)
                except Exception as e:
                    self.error.emit(f"SDR tune error: {e}")
                    self.stop_flag = True
                    break

                time.sleep(DWELL)
                for _ in range(CLR_READS):
                    self._safe_rx()

                acc0 = np.zeros(NUM_S * NUM_R, np.complex64)
                acc1 = np.zeros_like(acc0)
                for j in range(NUM_R):
                    r = self._safe_rx()
                    acc0[j*NUM_S:(j+1)*NUM_S] = (r[0]/2**12) * 7
                    acc1[j*NUM_S:(j+1)*NUM_S] = (r[1]/2**12) * 7

                s21 = to_dB(lockin(acc0))
                if self.cal21 is not None:
                    s21 -= np.interp(f, self.cal21['freqs'], self.cal21['db'])

                if self.skip_mask[i]:
                    s11 = np.nan
                else:
                    A11 = lockin(acc1)
                    if self.cal11L is not None:
                        L = np.interp(f, self.cal11L['freqs'], self.cal11L['linear'])
                        O = np.interp(f, self.cal11O['freqs'], self.cal11O['linear'])
                        s11 = to_dB(abs((A11 - L) / (O - L)))
                    else:
                        s11 = to_dB(A11)

                # clamp S11 at max 0 dB
                if not np.isnan(s11) and s11 > 0:
                    s11 = 0.0

                self.update.emit(f, s21, s11)

            if not self.stop_flag:
                time.sleep(1)

# ───────────── GUI with toggles + markers ─────────────
class VNA(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("PlutoSDR VNA Pro – interactive")
        self._build_ui()
        self._init_plot()
        self._spawn_worker()

        if os.path.exists('cal_s21.npz'):
            self.load_s21()
        if os.path.exists('cal_s11_load.npz') and os.path.exists('cal_s11_open.npz'):
            self.load_s11()

    # --- UI bar + checkboxes ---
    def _build_ui(self):
        cw = QWidget()
        self.setCentralWidget(cw)
        v = QVBoxLayout(cw)

        top = QFrame()
        top.setFrameStyle(QFrame.Shape.Box | QFrame.Shadow.Raised)
        h = QHBoxLayout(top)

        for lbl, attr, val in [
            ("Start (MHz):","le0", int(MIN_FREQ/1e6)),
            ("Stop (MHz):", "le1", int(MAX_FREQ/1e6)),
            ("Steps:",       "leN", DEFAULT_STEPS)
        ]:
            h.addWidget(QLabel(lbl))
            le = QLineEdit(str(val))
            setattr(self, attr, le)
            h.addWidget(le)

        pb = QPushButton("Apply")
        pb.clicked.connect(self.apply_span)
        h.addWidget(pb)

        for txt, fn in [
            ("Cal S21", self.cal_s21),
            ("Load S21", self.load_s21),
            ("Cal S11", self.cal_s11),
            ("Load S11", self.load_s11)
        ]:
            btn = QPushButton(txt)
            btn.clicked.connect(fn)
            h.addWidget(btn)

        h.addWidget(QLabel("Show:"))
        self.cbSmooth = QCheckBox("Smoothed")
        self.cbSmooth.setChecked(True)
        h.addWidget(self.cbSmooth)
        self.cbRaw = QCheckBox("Raw")
        self.cbRaw.setChecked(True)
        h.addWidget(self.cbRaw)
        self.cbSmooth.stateChanged.connect(self._vis_toggle)
        self.cbRaw.stateChanged.connect(self._vis_toggle)

        v.addWidget(top)

        self.fig, (self.ax21, self.ax11) = plt.subplots(2,1, figsize=(12,9))
        self.canvas = FigureCanvas(self.fig)
        v.addWidget(self.canvas)
        v.addWidget(NavigationToolbar2QT(self.canvas, self))

    # --- traces + markers storage + titles ---
    def _init_plot(self):
        self.x21, self.y21r = [], []
        self.x11, self.y11r = [], []

        self.l21, = self.ax21.plot([], [], 'b-', lw=1.2, label='S21 smoothed')
        self.d21, = self.ax21.plot([], [], 'ro', ms=4,   label='S21 raw')
        self.l11, = self.ax11.plot([], [], 'g-', lw=1.2, label='S11 smoothed')
        self.d11, = self.ax11.plot([], [], 'mo', ms=4,   label='S11 raw')

        self.ax21.set_title("S21 (dB)")
        self.ax11.set_title("S11 (dB)")

        self.ax21.set_xlim(MIN_FREQ/1e9, MAX_FREQ/1e9)
        self.ax11.set_xlim(MIN_FREQ/1e9, MAX_FREQ/1e9)
        self.ax21.set_ylabel("dB")
        self.ax11.set_ylabel("dB")
        self.ax11.set_xlabel("Frequency (GHz)")

        for ax in (self.ax21, self.ax11):
            ax.grid(True)

        self.canvas.mpl_connect('button_press_event', self._on_click)
        self.markers = []
        self._vis_toggle()

    # --- worker startup ---
    def _spawn_worker(self):
        self.wk = SweepThread(sdr)
        self.wk.update.connect(self._update)
        self.wk.started.connect(self._reset)
        self.wk.error.connect(lambda m: QMessageBox.critical(self, "Worker error", m))
        self.wk.start()

    # --- span handling ---
    def apply_span(self):
        try:
            f0 = float(self.le0.text())*1e6
            f1 = float(self.le1.text())*1e6
            n  = int(self.leN.text())
            if not (MIN_FREQ <= f0 < f1 <= MAX_FREQ and n >= 2):
                raise ValueError
        except ValueError:
            print("Span input error")
            return
        self.wk.stop(); self.wk.wait()
        self.wk.set_span(f0, f1, n)
        self.wk.stop_flag = False
        self.wk.start()
        self.ax21.set_xlim(f0/1e9, f1/1e9)
        self.ax11.set_xlim(f0/1e9, f1/1e9)
        self.canvas.draw()

    # --- plotting updates ---
    def _reset(self):
        self.x21.clear(); self.y21r.clear()
        self.x11.clear(); self.y11r.clear()
        for ln in (self.l21, self.d21, self.l11, self.d11):
            ln.set_data([], [])
        self._clear_markers()
        self.canvas.draw()

    def _update(self, f, s21, s11):
        fGHz = f/1e9
        self.x21.append(fGHz); self.y21r.append(s21)
        self.x11.append(fGHz); self.y11r.append(s11)

        self.l21.set_data(self.x21, smooth_trace(np.array(self.y21r), SMOOTH_WIN21))
        self.d21.set_data(self.x21, self.y21r)

        # S11 processing (already clamped in worker)
        x = np.array(self.x11)
        y_raw = np.array(self.y11r)
        valid = np.isfinite(y_raw)
        # fill gaps
        if valid.sum() >= 2:
            y_filled = y_raw.copy()
            y_filled[~valid] = np.interp(x[~valid], x[valid], y_raw[valid])
        else:
            y_filled = y_raw.copy()
        # smoothing
        if self.cbSmooth.isChecked():
            y_plot = smooth_trace(y_filled, SMOOTH_WIN11)
        else:
            y_plot = y_filled

        # clamp smoothed as well
        y_plot = np.minimum(y_plot, 0.0)

        self.l11.set_data(x, y_plot)
        if self.cbRaw.isChecked():
            self.d11.set_data(x[valid], y_raw[valid])
        else:
            self.d11.set_data([], [])

        for ax in (self.ax21, self.ax11):
            ax.relim(); ax.autoscale_view(scalex=False)
        self.canvas.draw_idle()

    # --- visibility toggles ---
    def _vis_toggle(self):
        showS = self.cbSmooth.isChecked()
        showR = self.cbRaw.isChecked()
        self.l21.set_visible(showS); self.l11.set_visible(showS)
        self.d21.set_visible(showR); self.d11.set_visible(showR)
        self.canvas.draw_idle()

    # --- marker handling ---
    def _clear_markers(self):
        for m in self.markers: m.remove()
        self.markers.clear()

    def _on_click(self, event):
        if event.inaxes not in (self.ax21, self.ax11):
            return
        ax = event.inaxes
        if event.button == 1:
            if ax is self.ax21:
                xdata = self.x21
                ydata = self.y21r if self.d21.get_visible() else list(self.l21.get_ydata())
            else:
                xdata = self.x11
                ydata = self.y11r if self.d11.get_visible() else list(self.l11.get_ydata())
            if not xdata:
                return
            idx = int(np.argmin(np.abs(np.array(xdata) - event.xdata)))
            x = xdata[idx]; y = ydata[idx]
            if np.isnan(y):
                return
            mrk = ax.plot(x, y, 'kx', ms=8, mew=2)[0]
            txt = ax.annotate(f"{y:.2f} dB\n{x:.3f} GHz",
                              (x, y),
                              textcoords="offset points",
                              xytext=(5, 5),
                              fontsize=8,
                              color='k',
                              bbox=dict(boxstyle="round,pad=0.2", fc='w', alpha=0.7))
            self.markers.extend([mrk, txt])
            self.canvas.draw_idle()
        elif event.button == 3:
            self._clear_markers()
            self.canvas.draw_idle()

    # --- calibration helpers (unchanged) ---
    def _do_cal(self, msg, ch):
        freqs = np.linspace(MIN_FREQ, MAX_FREQ, CAL_POINTS)
        out = {'freqs': [], 'linear': [], 'db': []}
        dlg = QProgressDialog(msg, "Cancel", 0, len(freqs), self)
        dlg.setWindowModality(Qt.WindowModality.ApplicationModal); dlg.show()
        for i, f in enumerate(freqs):
            if dlg.wasCanceled():
                return None
            dlg.setValue(i); QApplication.processEvents()
            NUM_R = 4 if f < 1e9 else 1
            sdr.tx_lo = sdr.rx_lo = int(f); time.sleep(DWELL)
            for _ in range(CLR_READS):
                sdr.rx()
            acc = np.zeros(NUM_S*NUM_R, np.complex64)
            for j in range(NUM_R):
                r = sdr.rx()
                acc[j*NUM_S:(j+1)*NUM_S] = (r[ch]/2**12)*7
            A = lockin(acc)
            out['freqs'].append(f)
            out['linear'].append(A)
            out['db'].append(to_dB(A))
        dlg.close()
        return {k: np.array(v) for k, v in out.items()}

    # --- S21 helpers ---
    def cal_s21(self):
        self.wk.stop(); self.wk.wait()
        data = self._do_cal("Calibrating S21…", 0)
        if data is not None:
            np.savez("cal_s21.npz", **data)
            self.load_s21()
        self.wk.stop_flag = False
        self.wk.start()

    def load_s21(self):
        d = np.load("cal_s21.npz")
        self.wk.load_cal21({'freqs': d['freqs'], 'db': d['db']})

    # --- S11 helpers ---
    def cal_s11(self):
        self.wk.stop(); self.wk.wait()
        if QMessageBox.information(self, "S11", "Connect 50 Ω LOAD then OK") \
           != QMessageBox.StandardButton.Ok:
            self.wk.stop_flag = False; self.wk.start(); return
        L = self._do_cal("Cal S11 (LOAD)…", 1)
        if L is None:
            self.wk.stop_flag = False; self.wk.start(); return
        if QMessageBox.information(self, "S11", "Connect OPEN then OK") \
           != QMessageBox.StandardButton.Ok:
            self.wk.stop_flag = False; self.wk.start(); return
        O = self._do_cal("Cal S11 (OPEN)…", 1)
        if O is None:
            self.wk.stop_flag = False; self.wk.start(); return

        np.savez("cal_s11_load.npz", **L)
        np.savez("cal_s11_open.npz", **O)
        pd.DataFrame({
            'freq_GHz': L['freqs']/1e9,
            'load_dB':  L['db'],
            'open_dB':  O['db']
        }).to_excel("s11_calibration.xlsx", index=False)
        print("✔️  S11 calibration XLSX written")
        self.load_s11()
        self.wk.stop_flag = False; self.wk.start()

    def load_s11(self):
        L = np.load("cal_s11_load.npz")
        O = np.load("cal_s11_open.npz")
        self.wk.load_cal11(
            {'freqs': L['freqs'], 'linear': L['linear']},
            {'freqs': O['freqs'], 'linear': O['linear']}
        )

    # --- cleanup ---
    def closeEvent(self, e):
        self.wk.stop(); self.wk.wait()
        sdr.tx_destroy_buffer(); sdr.rx_destroy_buffer()
        e.accept()

# ───────────── main ─────────────
if __name__ == "__main__":
    app = QApplication(sys.argv)
    try:
        win = VNA()
    except Exception:
        traceback.print_exc()
        QMessageBox.critical(None, "Init error", "Failed to start:\n" + traceback.format_exc())
        sys.exit(1)
    win.show()
    sys.exit(app.exec())
