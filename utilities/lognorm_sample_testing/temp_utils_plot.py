#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# ---- parameters: must match your C++ run ----
# sigma = 0.2
sigma = 1
# scaleBalls = 1e-5
# a_max = scaleBalls * np.exp(-5.0 * sigma**2 / 2.0)
a_max = 5

# ---- load data ----
data = np.loadtxt("lognorm_samples.txt")

# quick diagnostics
n_total = data.size
n_pos = np.sum(data > 0)
n_zero = np.sum(data == 0)
print(f"Total samples: {n_total}, >0: {n_pos}, ==0: {n_zero}")
if n_pos == 0:
    raise RuntimeError("All samples are <= 0; fix the C++ sampler.")

# keep only positive for log plotting
data = data[data > 0]
dmin, dmax = data.min(), data.max()

# ---- analytical PDF (safe near 0) ----
def lndpdf(a, sigma, a_max):
    # vectorized; return 0 for nonpositive a to avoid log/divide-by-zero
    a = np.asarray(a)
    out = np.zeros_like(a, dtype=float)
    mask = a > 0
    aa = a[mask]
    out[mask] = (2/np.sqrt(np.pi)) / (aa * sigma * 2*np.sqrt(2)) * \
        np.exp(-((np.log(aa/a_max) - sigma**2)**2) / (2*sigma**2))
    return out

# x-range focused on actual data (with a bit of padding)
xmin = max(dmin*0.8, a_max*1e-3)   # avoid going too close to 0
xmin = max(xmin, 1e-20)
# xmax = max(dmax*1.2, 3.0*scaleBalls)
xmax = max(dmax*1.2, 3.0)

# ---- log-spaced bins & PDF curve ----
bins = np.geomspace(xmin, xmax, 80)
a_vals = np.geomspace(xmin, xmax, 800)
pdf_vals = lndpdf(a_vals, sigma, a_max)

# ---- plot ----
plt.figure()
plt.hist(data, bins=bins, density=True, alpha=0.6, label="Sampled radii")
plt.plot(a_vals, pdf_vals, lw=2, label="Analytical lognormal PDF")
plt.xscale("log")
plt.xlabel("Radius a (m)")
plt.ylabel("Probability density")
# plt.title(f"Lognormal sampling (σ={sigma}, scaleBalls={scaleBalls:.1e})")
plt.legend()
plt.tight_layout()
plt.show()
