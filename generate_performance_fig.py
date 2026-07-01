#!/usr/bin/env python3
"""
Generate wavemaker performance figure with ITTC sea state regions.

Usage:
    python generate_performance_fig.py [--lam 50] [--ocean "N. Atlantic"] [--dpi 300] [--out figure.png]

Arguments:
    --lam    Froude scale factor λ            (default: 50)
    --ocean  Ocean for SS representative Tp   (default: "N. Atlantic"; also "N. Pacific")
    --dpi    Output resolution in DPI         (default: 300)
    --out    Output filename                  (default: auto-generated from lam and ocean)
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.interpolate import UnivariateSpline
from scipy.optimize import brentq

# ── Physical constants ────────────────────────────────────────────────
G = 9.81
H_DEPTH = 3.0  # tank water depth (m)

# ── ITTC (2014) sea state table ───────────────────────────────────────
# Source: ITTC Recommended Procedures and Guidelines, Table 1 (2014).
SEA_STATES = {
    "N. Atlantic": [
        {"ss": 1, "Hs": 0.05,  "Tp":  4.0},
        {"ss": 2, "Hs": 0.30,  "Tp":  6.5},
        {"ss": 3, "Hs": 0.90,  "Tp":  8.5},
        {"ss": 4, "Hs": 1.90,  "Tp": 10.0},
        {"ss": 5, "Hs": 3.25,  "Tp": 12.0},
        {"ss": 6, "Hs": 5.00,  "Tp": 13.5},
        {"ss": 7, "Hs": 7.50,  "Tp": 15.5},
        {"ss": 8, "Hs": 11.50, "Tp": 18.0},
    ],
    "N. Pacific": [
        {"ss": 1, "Hs": 0.05,  "Tp":  4.5},
        {"ss": 2, "Hs": 0.30,  "Tp":  7.5},
        {"ss": 3, "Hs": 0.90,  "Tp": 10.0},
        {"ss": 4, "Hs": 1.90,  "Tp": 12.0},
        {"ss": 5, "Hs": 3.25,  "Tp": 14.0},
        {"ss": 6, "Hs": 5.00,  "Tp": 16.0},
        {"ss": 7, "Hs": 7.50,  "Tp": 17.5},
        {"ss": 8, "Hs": 11.50, "Tp": 20.0},
    ],
}

# Hs ranges from ITTC Table 1; TpHalf is the ITTC modal period half-width.
# Tp box is centered on the selected ocean's representative Tp.
SS_REGIONS = [
    {"ss": 3, "HsMin": 0.5,  "HsMax": 1.25, "TpHalf": 3.0, "color": "#66BB6A"},
    {"ss": 4, "HsMin": 1.25, "HsMax": 2.5,  "TpHalf": 3.0, "color": "#D4E157"},
    {"ss": 5, "HsMin": 2.5,  "HsMax": 4.0,  "TpHalf": 3.0, "color": "#FFA726"},
    {"ss": 6, "HsMin": 4.0,  "HsMax": 6.0,  "TpHalf": 3.5, "color": "#FF7043"},
    {"ss": 7, "HsMin": 6.0,  "HsMax": 9.0,  "TpHalf": 4.0, "color": "#EF5350"},
]


def get_wavelength(T):
    """Wavelength (m) via dispersion relation at depth H_DEPTH."""
    omega = 2 * np.pi / T
    lam0  = G * T**2 / (2 * np.pi)
    return brentq(
        lambda lam: omega**2 - G * (2*np.pi/lam) * np.tanh(2*np.pi*H_DEPTH/lam),
        1e-3, 1.01 * lam0,
    )


def build_figure(lam, ocean, dpi, out):
    # ── Wavemaker spline data (model scale, d = 3 m) ─────────────────
    x1 = np.array([0, 0.17425186, 0.46584599, 0.71508047, 0.95787431,
                   1.22014491, 1.46761441, 1.68111569, 1.85563247,
                   2.02606192, 2.27167354])
    y1 = np.array([0, 0.00999409, 0.03567234, 0.06999511, 0.11503782,
                   0.17535982, 0.24219987, 0.30754080, 0.36451806,
                   0.42013116, 0.49912825])
    x2 = np.array([2.26510903, 2.39265243])
    y2 = np.array([0.50045593, 0.49950394])
    x3 = np.array([2.39280726, 2.56952253, 2.78556290, 3.09260839,
                   3.44535771, 3.87830548, 4.63566990, 5.28558699,
                   5.80331597, 6.29162853])
    y3 = np.array([0.50020426, 0.45109096, 0.39920193, 0.34560953,
                   0.29798440, 0.25608582, 0.20624667, 0.17796777,
                   0.16023373, 0.14659946])

    x_all   = np.concatenate([x1, x2, x3])
    y_all   = np.concatenate([y1, y2, y3])
    idx     = np.argsort(x_all)
    spline  = UnivariateSpline(x_all[idx], y_all[idx], k=3, s=0.001)
    T_plot  = np.linspace(0, 5, 500)
    H_reg   = spline(T_plot)
    H_irr   = H_reg / 1.52

    # ── Limit curves ─────────────────────────────────────────────────
    T_lim    = np.linspace(0.15, 5.0, 500)
    L_lim    = np.array([get_wavelength(T) for T in T_lim])
    H_break  = 0.1  * L_lim * np.tanh(2*np.pi*H_DEPTH/L_lim)
    H_airy   = 0.03 * L_lim
    mask_brk  = (T_lim >= 0.5) & (T_lim <= 2.0) & (H_break <= 0.6)
    mask_airy = (T_lim >= 0.5) & (T_lim <= 3.0) & (H_airy  <= 0.6)

    # ── Depth regime boundary ─────────────────────────────────────────
    T_di = 2*np.pi / np.sqrt(G * (2*np.pi/(2*H_DEPTH)) * np.tanh(np.pi))  # h/L = 1/2

    # ── Figure ───────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 6))

    # depth zone backgrounds
    ax.axvspan(0,    T_di, color='#cce0ff', alpha=0.35)
    ax.axvspan(T_di, 5.0,  color='#d4f5d4', alpha=0.35)
    ax.axvline(T_di, color='steelblue', linewidth=1.2, linestyle=':', alpha=0.8)

    # ── ITTC sea state regions (model scale) ─────────────────────────
    sqrt_lam = np.sqrt(lam)
    for reg in SS_REGIONS:
        rep_tp      = SEA_STATES[ocean][reg["ss"] - 1]["Tp"]
        tp_min      = (rep_tp - reg["TpHalf"]) / sqrt_lam
        tp_max      = (rep_tp + reg["TpHalf"]) / sqrt_lam
        hs_min      = reg["HsMin"] / lam
        hs_max      = reg["HsMax"] / lam
        rect = mpatches.Rectangle(
            (tp_min, hs_min), tp_max - tp_min, hs_max - hs_min,
            linewidth=1.5, edgecolor=reg["color"], facecolor=reg["color"],
            alpha=0.2, zorder=2,
        )
        ax.add_patch(rect)
        ax.text(tp_min + 0.03, hs_max - (hs_max - hs_min) * 0.18,
                f'SS{reg["ss"]}', color=reg["color"],
                fontsize=9, style='italic', clip_on=True, zorder=3)

    # limit curves
    ax.plot(T_lim[mask_brk],  H_break[mask_brk],
            color='red',  linewidth=1.8, linestyle='--',
            label=r'breaking limit ($H/\lambda = 1/10$)')
    ax.plot(T_lim[mask_airy], H_airy[mask_airy],
            color='teal', linewidth=1.5, linestyle='-.',
            label=r'Airy theory limit ($H/\lambda = 0.03$)')

    # wavemaker envelopes
    ax.plot(T_plot, H_reg, label='regular',   color='dodgerblue', linewidth=2.5)
    ax.plot(T_plot, H_irr, label='irregular', color='#9932CC',    linewidth=2.5)

    # depth zone labels
    ax.text(T_di / 2,        0.575, 'deep water',
            ha='center', va='top', fontsize=9, color='steelblue', style='italic')
    ax.text((T_di + 5) / 2, 0.575, 'intermediate water',
            ha='center', va='top', fontsize=9, color='seagreen',  style='italic')

    ax.set_xlabel('period (s)', fontsize=16)
    ax.set_ylabel('height (m)', fontsize=16)
    ax.set_xlim(0, 5)
    ax.set_ylim(0, 0.6)
    ax.tick_params(labelsize=14)
    ax.grid(True, alpha=0.35)
    ax.legend(fontsize=11, loc='upper right', bbox_to_anchor=(1.0, 0.92))
    ax.set_title(f'Wavemaker performance  —  1:{lam:.0f} scale  ·  {ocean}', fontsize=11)

    plt.tight_layout()

    if out is None:
        ocean_tag = "natl" if ocean == "N. Atlantic" else "npac"
        out = f"wavemaker_performance_lam{lam:.0f}_{ocean_tag}.png"

    plt.savefig(out, dpi=dpi)
    print(f"Saved: {out}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Generate wavemaker performance figure with ITTC sea state regions."
    )
    parser.add_argument("--lam",   type=float, default=50,
                        help="Froude scale factor λ (default: 50)")
    parser.add_argument("--ocean", type=str,   default="N. Atlantic",
                        choices=["N. Atlantic", "N. Pacific"],
                        help='Ocean (default: "N. Atlantic")')
    parser.add_argument("--dpi",   type=int,   default=300,
                        help="Output DPI (default: 300)")
    parser.add_argument("--out",   type=str,   default=None,
                        help="Output filename (default: auto-generated)")
    args = parser.parse_args()
    build_figure(args.lam, args.ocean, args.dpi, args.out)


if __name__ == "__main__":
    main()
