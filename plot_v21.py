"""Figures for v2.1 post-registration report."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'v19'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'v21_prereg'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

res_all = np.load('v21_results.npy', allow_pickle=True).item()
MASSES = [0.017, 0.05, 0.1, 0.2]
MASS_LABELS = ['17', '50', '100', '200']

# Load published values
r17 = np.load('v17/v17_results.npy', allow_pickle=True).item()
r18 = np.load('v18/v18_results.npy', allow_pickle=True).item()

PUB = {}
for m in MASSES:
    PUB[('C1', m)] = r17[('deg', '3D', m)]['deg']
    PUB[('C2', m)] = r17[('deg', '3D_allpixel', m)]['deg']
    PUB[('C3', m)] = r18['M025'][('deg', '3D_allpixel', m)]['deg']

# ── Figure 15(a): degradation vs mass ──
fig, ax = plt.subplots(1, 1, figsize=(7, 5))
x = np.array(MASSES) * 1000

styles = {
    'C1': dict(color='#2196F3', marker='o', label='f = 1, pixel z (Config 1)'),
    'C2': dict(color='#4CAF50', marker='s', label='f = 1, z every layer (Config 2)'),
    'C3': dict(color='#FF5722', marker='D', label='f = 0.25, z every layer (Config 3)'),
}

for cfg in ['C1', 'C2', 'C3']:
    hv = [PUB[(cfg, m)] for m in MASSES]
    tag = [res_all[cfg][('tag', m)]['deg'] for m in MASSES]
    st = styles[cfg]
    ax.plot(x, hv, st['marker'], color=st['color'], ms=8, mfc='none', mew=1.5, linestyle='--', linewidth=1, alpha=0.6)
    ax.plot(x, tag, st['marker'], color=st['color'], ms=8, mfc=st['color'], linestyle='-', linewidth=1.5, label=st['label'])

# Prediction ranges (shaded) for P1, P2
ax.axhspan(16, 26.7, xmin=0, xmax=0.15, alpha=0.08, color='#FF5722')
# Add a legend entry for open = hard veto, filled = tag
ax.plot([], [], 'ko', ms=7, mfc='none', mew=1.5, linestyle='--', label='Hard veto (published)')
ax.plot([], [], 'ko', ms=7, mfc='k', linestyle='-', label='Likelihood tag (v2.1)')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('ALP mass [MeV]', fontsize=12)
ax.set_ylabel('Limit degradation from conversions', fontsize=12)
ax.set_xticks([17, 50, 100, 200])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.legend(fontsize=8.5, loc='upper right')
ax.set_title('v2.1: Likelihood material tag vs hard material veto', fontsize=11)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig('fig15a_deg_vs_mass.png', dpi=200)
print('Saved fig15a_deg_vs_mass.png')

# ── Figure 15(b): λ distributions for C3, 17 MeV ──
fig, ax = plt.subplots(1, 1, figsize=(7, 4.5))
r3 = res_all['C3']
lam_conv = r3['lam_conv']
lam_sig = r3['lam_sig17']
best_t = r3['best_t_17']

bins = np.linspace(-15, 40, 120)
ax.hist(lam_conv, bins=bins, density=True, alpha=0.5, color='#2196F3', label='Conversions')
ax.hist(lam_sig, bins=bins, density=True, alpha=0.5, color='#FF5722', label='ALP 17 MeV signal')
if best_t is not None:
    ax.axvline(best_t, color='k', linestyle='--', linewidth=1.5, label=f'Threshold t = {best_t}')
ax.set_xlabel('λ (likelihood material tag)', fontsize=12)
ax.set_ylabel('Normalised density', fontsize=12)
ax.set_title('v2.1 Config 3 (f = 0.25, z every layer): λ distributions', fontsize=11)
ax.legend(fontsize=10)
ax.set_xlim(-15, 40)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig('fig15b_lambda_distributions.png', dpi=200)
print('Saved fig15b_lambda_distributions.png')

# ── Figure 15(c): ratio tag/hard-veto vs mass for all configs ──
fig, ax = plt.subplots(1, 1, figsize=(7, 4.5))
for cfg in ['C1', 'C2', 'C3']:
    ratios = [res_all[cfg][('tag', m)]['deg'] / PUB[(cfg, m)] for m in MASSES]
    st = styles[cfg]
    ax.plot(x, ratios, st['marker'] + '-', color=st['color'], ms=8, linewidth=1.5, label=st['label'])
ax.axhline(1.0, color='grey', linestyle=':', linewidth=1)
ax.axhline(1.05, color='red', linestyle=':', linewidth=0.8, alpha=0.5)
ax.axhline(0.95, color='green', linestyle=':', linewidth=0.8, alpha=0.5)
ax.set_xscale('log')
ax.set_xlabel('ALP mass [MeV]', fontsize=12)
ax.set_ylabel('Tag degradation / Hard veto degradation', fontsize=12)
ax.set_xticks([17, 50, 100, 200])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_title('v2.1: Tag / hard veto ratio (< 1 means tag is better)', fontsize=11)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig('fig15c_ratio.png', dpi=200)
print('Saved fig15c_ratio.png')

print('All figures done.')
