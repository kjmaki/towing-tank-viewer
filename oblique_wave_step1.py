#!/usr/bin/env python3
"""
Oblique Wave Basin — Step 1
Original wave (no reflection) with propagation rays and first-reflection locus.

Edit the WAVE PARAMETERS block to change T, a, theta, phi.
Water depth h = 3 m, basin 100 m x 6.5 m.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider
from scipy.optimize import brentq

# =============================================================================
# BASIN (fixed)
# =============================================================================
L = 100.0   # length [m]  x-direction (along-flume)
W = 6.5     # width  [m]  y-direction (cross-flume)
g = 9.81    # m/s²
h = 3.0     # water depth [m]

# =============================================================================
# WAVE PARAMETERS  ← edit here
# =============================================================================
T     = 2.0    # wave period [s]
a     = 0.10   # amplitude [m]
theta = 20.0   # propagation angle from x-axis [deg], positive toward y = W wall
phi   = 0.0    # initial phase offset [rad]

# =============================================================================
# DISPERSION  ω² = g k tanh(kh)
# =============================================================================
omega     = 2.0 * np.pi / T
theta_rad = np.radians(theta)

k   = brentq(lambda kk: omega**2 - g * kk * np.tanh(kk * h), 1e-9, 1e5)
lam = 2.0 * np.pi / k     # wavelength [m]
c   = omega / k            # phase speed [m/s]
cg  = 0.5 * c * (1.0 + 2.0*k*h / np.sinh(2.0*k*h))  # group speed [m/s]
kx  = k * np.cos(theta_rad)
ky  = k * np.sin(theta_rad)

print(f"\nWave parameters:")
print(f"  T = {T} s,  a = {a} m,  theta = {theta} deg,  phi = {phi} rad")
print(f"  omega = {omega:.4f} rad/s")
print(f"  k = {k:.4f} rad/m,  lambda = {lam:.2f} m,  c = {c:.2f} m/s,  cg = {cg:.2f} m/s")
print(f"  kx = {kx:.4f} rad/m,  ky = {ky:.4f} rad/m")
print(f"  h = {h} m,  kh = {k*h:.3f}")

# =============================================================================
# GRID
# =============================================================================
nx, ny = 600, 80
x = np.linspace(0, L, nx)
y = np.linspace(0, W, ny)
X, Y = np.meshgrid(x, y)

def surface(t):
    return a * np.cos(omega * t - kx * X - ky * Y + phi)

# =============================================================================
# GROUP-VELOCITY MASK
# The wave front is the plane  x·cosθ + y·sinθ = cg·t,  perpendicular to the
# propagation direction.  Points behind the front (smaller projection) are
# illuminated; points ahead are masked to zero.
# =============================================================================
def compute_mask(t, smooth=1.0):
    proj = X * np.cos(theta_rad) + Y * np.sin(theta_rad)
    return 0.5 * (1.0 + np.tanh((cg * t - proj) / smooth))

def front_xy(t):
    """
    Endpoints of the group-velocity front line clipped to the basin.
    Line equation: x·cosθ + y·sinθ = cg·t
    Returns ([x0,x1],[y0,y1]) or None if the front is outside the basin.
    """
    d     = cg * t
    cos_t = np.cos(theta_rad)
    sin_t = np.sin(theta_rad)
    pts   = []

    def add(x_, y_):
        if -1e-9 <= x_ <= L + 1e-9 and -1e-9 <= y_ <= W + 1e-9:
            pts.append((float(np.clip(x_, 0, L)), float(np.clip(y_, 0, W))))

    if abs(cos_t) > 1e-9:
        add(d / cos_t,               0.0)   # y = 0 wall
        add((d - W * sin_t) / cos_t, W)     # y = W wall
    if abs(sin_t) > 1e-9:
        add(0.0, d / sin_t)                  # x = 0 wavemaker
        add(L,   (d - L * cos_t) / sin_t)   # x = L far wall

    seen, unique = set(), []
    for p in pts:
        key = (round(p[0], 6), round(p[1], 6))
        if key not in seen:
            seen.add(key); unique.append(p)
    if len(unique) < 2:
        return None
    unique.sort(key=lambda p: p[0] * cos_t + p[1] * sin_t)
    return [unique[0][0], unique[-1][0]], [unique[0][1], unique[-1][1]]

# Time for the group front to sweep fully across the basin
t_max_fill = 1.2 * (L * np.cos(theta_rad) + W * abs(np.sin(theta_rad))) / cg
print(f"  Basin fill time: {t_max_fill:.1f} s  (slider range)\n")

# =============================================================================
# FIRST REFLECTED WAVE
# Reflection from the nearest side wall (y=W for theta>0, y=0 for theta<0).
# Image principle: replace y → 2·wall_y − y in the incident phase.
# The reflected wave front is the image of the incident front through wall_y.
# =============================================================================
wall_y = W if ky > 0.0 else 0.0   # wall hit first

def surface_reflected(t):
    return a * np.cos(omega * t - kx * X + ky * Y + phi - 2.0 * ky * wall_y)

def compute_mask_reflected(t, smooth=1.0):
    """Group-velocity envelope for the first reflected wave."""
    proj_r = X * np.cos(theta_rad) + (2.0 * wall_y - Y) * np.sin(theta_rad)
    return 0.5 * (1.0 + np.tanh((cg * t - proj_r) / smooth))

def front_xy_reflected(t):
    """
    Endpoints of the reflected wave front clipped to the basin.
    Line: x·cosθ + (2·wall_y − y)·sinθ = cg·t
    ↔    x·cosθ − y·sinθ = cg·t − 2·wall_y·sinθ  ≡  dr
    """
    dr    = cg * t - 2.0 * wall_y * np.sin(theta_rad)
    cos_t = np.cos(theta_rad)
    sin_t = np.sin(theta_rad)
    pts   = []

    def add(x_, y_):
        if -1e-9 <= x_ <= L + 1e-9 and -1e-9 <= y_ <= W + 1e-9:
            pts.append((float(np.clip(x_, 0, L)), float(np.clip(y_, 0, W))))

    if abs(cos_t) > 1e-9:
        add(dr / cos_t,               0.0)         # y = 0 wall
        add((dr + W * sin_t) / cos_t, W)           # y = W wall
    if abs(sin_t) > 1e-9:
        add(0.0, -dr / sin_t)                       # x = 0 wavemaker
        add(L,   (L * cos_t - dr) / sin_t)         # x = L far wall

    seen, unique = set(), []
    for p in pts:
        key = (round(p[0], 6), round(p[1], 6))
        if key not in seen:
            seen.add(key); unique.append(p)
    if len(unique) < 2:
        return None
    unique.sort(key=lambda p: p[0])
    return [unique[0][0], unique[-1][0]], [unique[0][1], unique[-1][1]]

# =============================================================================
# WAVE RAYS  (straight paths from wavemaker x=0, no reflection)
# =============================================================================
def make_rays(n=9):
    """Ray paths from n evenly-spaced points on the wavemaker face."""
    rays = []
    for y0 in np.linspace(0, W, n):
        # Parametric: position = (s*cos(theta), y0 + s*sin(theta))
        # Find first boundary crossing
        candidates = []

        # Far end x = L
        if abs(np.cos(theta_rad)) > 1e-12:
            candidates.append(L / np.cos(theta_rad))

        # Side wall y = W (only if wave has positive y-component and not at wall)
        if ky > 1e-12 and (W - y0) > 1e-9:
            candidates.append((W - y0) / np.sin(theta_rad))

        # Side wall y = 0 (only if wave has negative y-component and not at wall)
        if ky < -1e-12 and y0 > 1e-9:
            candidates.append(y0 / abs(np.sin(theta_rad)))

        if not candidates:
            continue

        s_hit = min(s for s in candidates if s > 1e-9)
        xe = np.clip(s_hit * np.cos(theta_rad), 0.0, L)
        ye = np.clip(y0 + s_hit * np.sin(theta_rad), 0.0, W)
        rays.append(([0.0, xe], [y0, ye]))
    return rays

# =============================================================================
# FIRST REFLECTION LOCUS
# The segment on the side wall where wave rays from the full wavemaker first arrive.
#   theta > 0: rays hit y = W,  locus spans x in [0, W/tan(theta)]
#   theta < 0: rays hit y = 0,  locus spans x in [0, W/|tan(theta)|]
# The star marks where the corner ray (from y0 = 0 for theta>0) arrives last.
# =============================================================================
def make_locus():
    if abs(ky) < 1e-10:
        return None

    wall_y = W if ky > 0 else 0.0
    x_far  = W / abs(np.tan(theta_rad))   # ray from the far wavemaker corner
    x_far  = min(x_far, L)                # clip if wave hits far end first

    # For theta>0: corner ray from (0,0) travels the longest before hitting y=W
    # For theta<0: corner ray from (0,W) travels the longest before hitting y=0
    if ky > 0:
        star_y0 = 0.0
    else:
        star_y0 = W

    print(f"\n  First reflection locus: x in [0, {x_far:.2f}] m  on  y = {wall_y} m wall")
    print(f"  Corner ray from (0, {star_y0:.1f}) first hits wall at  ({x_far:.2f}, {wall_y}) m")

    return {
        'lx'    : np.array([0.0, x_far]),
        'ly'    : np.array([wall_y, wall_y]),
        'star'  : (x_far, wall_y),
        'wall_y': wall_y,
        'x_far' : x_far,
    }

rays  = make_rays()
locus = make_locus()
print()

# =============================================================================
# FIGURE LAYOUT
# =============================================================================
levels = np.linspace(-2*a, 2*a, 41)   # doubled range: incident + reflected can reach ±2a

fig = plt.figure(figsize=(18, 3.0))
gs  = gridspec.GridSpec(2, 1, fig, height_ratios=[6, 1], hspace=0.45)
ax    = fig.add_subplot(gs[0])
ax_sl = fig.add_subplot(gs[1])

# --- Initial contourf + contour lines (masked by group-velocity front) ---
t0  = 0.0
Z0  = (surface(t0)           * compute_mask(t0) +
       surface_reflected(t0) * compute_mask_reflected(t0))
cf  = [ax.contourf(X, Y, Z0, levels=levels, cmap='RdBu_r', extend='both', zorder=1)]
cs  = [ax.contour (X, Y, Z0, levels=levels[::5], colors='k',
                   linewidths=0.35, alpha=0.45, zorder=2)]
cb  = fig.colorbar(cf[0], ax=ax, label='η [m]', fraction=0.015, pad=0.02)

# --- Static: wave rays (drawn once, always on top of contourf) ---
for xr, yr in rays:
    ax.plot(xr, yr, 'w--', lw=0.9, alpha=0.65, zorder=3)

# --- Static: first-reflection locus ---
if locus:
    wall_label = f"y = {locus['wall_y']:.1f} m"
    ax.plot(locus['lx'], locus['ly'], color='gold', lw=3.5, zorder=4,
            label=f'1st reflection locus  ({wall_label})')
    ax.plot(locus['star'][0], locus['star'][1],
            marker='*', color='gold', markersize=14, zorder=5,
            label=(f"Corner ray (0, 0) → ({locus['x_far']:.1f} m, {locus['wall_y']:.1f} m)"))

# --- Dynamic: incident group-velocity front ---
_f0 = front_xy(t0)
front_artist, = ax.plot(
    _f0[0] if _f0 else [], _f0[1] if _f0 else [],
    color='lime', lw=2.0, zorder=6, label='incident front')

# --- Dynamic: reflected wave front ---
_fr0 = front_xy_reflected(t0)
front_artist_r, = ax.plot(
    _fr0[0] if _fr0 else [], _fr0[1] if _fr0 else [],
    color='orange', lw=2.0, ls='--', zorder=6, label='reflected front')

ax.set_xlim(0, L)
ax.set_ylim(0, W)
ax.set_aspect('equal')
ax.set_xlabel('x [m]', fontsize=11)
ax.set_ylabel('y [m]', fontsize=11)
ax.legend(loc='upper right', fontsize=9, framealpha=0.8)
title = [ax.set_title(
    f'Oblique wave + 1st reflection  │  '
    f'T = {T} s,  a = {a} m,  θ = {theta}°,  '
    f'λ = {lam:.1f} m,  h = {h} m  │  t = {t0:.2f} s',
    fontsize=10)]

# =============================================================================
# TIME SLIDER
# =============================================================================
slider = Slider(ax_sl, 't  [s]', 0.0, t_max_fill,
               valinit=t0, valstep=t_max_fill / 200)

def update(val):
    t = slider.val
    cf[0].remove()
    cs[0].remove()
    Z = (surface(t)           * compute_mask(t) +
         surface_reflected(t) * compute_mask_reflected(t))
    cf[0] = ax.contourf(X, Y, Z, levels=levels, cmap='RdBu_r', extend='both', zorder=1)
    cs[0] = ax.contour (X, Y, Z, levels=levels[::5], colors='k',
                        linewidths=0.35, alpha=0.45, zorder=2)
    _f  = front_xy(t)
    front_artist.set_data(_f[0] if _f else [], _f[1] if _f else [])
    _fr = front_xy_reflected(t)
    front_artist_r.set_data(_fr[0] if _fr else [], _fr[1] if _fr else [])
    title[0].set_text(
        f'Oblique wave + 1st reflection  │  '
        f'T = {T} s,  a = {a} m,  θ = {theta}°,  '
        f'λ = {lam:.1f} m,  h = {h} m  │  t = {t:.2f} s')
    fig.canvas.draw_idle()

slider.on_changed(update)

fig.subplots_adjust(left=0.04, right=0.91, top=0.88, bottom=0.22, hspace=1.5)
plt.show()
