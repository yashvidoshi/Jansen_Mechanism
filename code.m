"""
theo_jansen_complete.py
=======================
Theo Jansen Walking Mechanism — Complete Simulation
IE410: Introduction to Robotics — Project Part B

All-in-one file combining:
  1. Kinematics engine  (jansen_kinematics.py)
  2. Main simulation    (sim_jansen.py)
  3. Link variation     (link_variation.py)

Produces five output files:
  foot_trajectory.png      — closed shoe-sole foot trajectory
  mechanism_animation.gif  — animated mechanism (dark background, new link colours)
  gait_comparison.png      — simulated vs reference gait curve
  link_variation_m.png     — effect of varying crank length m
  link_variation_h.png     — effect of varying long connector h

Link lengths: Theo Jansen's published "holy numbers" (2007).
Reference: Shin et al. (2018), JMR; Jadav et al., Single-DOF Gait Trainer.

Usage:
    python theo_jansen_complete.py
"""

import os
import copy
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm
import imageio.v2 as imageio

OUT_DIR = os.path.dirname(os.path.abspath(__file__))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — KINEMATICS ENGINE
# ═══════════════════════════════════════════════════════════════════════════════

# Theo Jansen's "holy numbers" (unscaled, Jansen 2007)
DEFAULT_LINKS = {
    'a': 38.0,  'b': 41.5,  'c': 39.3,  'd': 40.1,
    'e': 55.8,  'f': 39.4,  'g': 36.7,  'h': 65.7,
    'i': 49.0,  'j': 50.0,  'k': 61.9,  'l': 7.8,
    'm': 15.0,
}

# Scale factor so plots show convenient mm-like units
SCALE = 5.4

# ── NEW LINK COLOURS ───────────────────────────────────────────────────────────
# Vivid palette chosen for contrast on dark (black) background.
# Format: (start_joint, end_joint, hex_colour, label)
LINK_EDGES = [
    ('O',    'C',    '#FF4444', 'm  (crank)'),           # bright red
    ('O',    'P',    '#AAAAAA', 'a–l (ground)'),          # light grey
    ('C',    'B_up', '#FFD700', 'j  (upper coupler)'),    # gold
    ('P',    'B_up', '#00CFFF', 'b  (upper rocker)'),     # sky blue
    ('B_up', 'D_up', '#39FF14', 'c'),                     # neon green
    ('P',    'D_up', '#BF5FFF', 'd'),                     # violet
    ('C',    'B_lo', '#FF6EC7', 'e  (lower coupler)'),    # hot pink
    ('P',    'B_lo', '#FF9900', 'f'),                     # orange
    ('B_lo', 'D_lo', '#FF4444', 'g'),                     # bright red
    ('D_up', 'D_lo', '#00CFFF', 'h  (long connector)'),   # sky blue
    ('B_lo', 'F',    '#39FF14', 'i  (foot side)'),        # neon green
    ('D_lo', 'F',    '#FFD700', 'k  (foot side)'),        # gold
]

# ── BACKGROUND / FOREGROUND COLOURS ──────────────────────────────────────────
BG_COLOR   = '#0D0D0D'   # near-black background
FG_COLOR   = '#E8E8E8'   # off-white text / axes
GRID_COLOR = '#2A2A2A'   # subtle dark grid
FOOT_COLOR = '#FF00FF'   # magenta foot trace (matches original video)


def _apply_dark_style(ax, fig):
    """Apply dark background to a figure and axes."""
    fig.patch.set_facecolor(BG_COLOR)
    ax.set_facecolor(BG_COLOR)
    for spine in ax.spines.values():
        spine.set_edgecolor(FG_COLOR)
    ax.tick_params(colors=FG_COLOR, which='both')
    ax.xaxis.label.set_color(FG_COLOR)
    ax.yaxis.label.set_color(FG_COLOR)
    ax.title.set_color(FG_COLOR)
    ax.grid(True, color=GRID_COLOR, alpha=0.6, linewidth=0.8)


# ── CIRCLE–CIRCLE INTERSECTION ────────────────────────────────────────────────
def cci(c1, r1, c2, r2, branch='+'):
    """
    Return one intersection point of two circles.

    branch='+' → left  of directed vector c1→c2
    branch='-' → right of directed vector c1→c2
    Returns None if circles do not intersect.
    """
    c1 = np.asarray(c1, dtype=float)
    c2 = np.asarray(c2, dtype=float)
    d_vec = c2 - c1
    d = np.linalg.norm(d_vec)
    if d == 0.0 or d > r1 + r2 or d < abs(r1 - r2):
        return None
    p = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d)
    h_sq = r1 * r1 - p * p
    if h_sq < 0.0:
        return None
    h = np.sqrt(max(h_sq, 0.0))
    midpoint = c1 + p * d_vec / d
    perp = np.array([-d_vec[1] / d, d_vec[0] / d])
    return midpoint + h * perp if branch == '+' else midpoint - h * perp


# ── FORWARD KINEMATICS ────────────────────────────────────────────────────────
def solve_pose(theta, links=None, scale=None):
    """
    Solve full mechanism pose for crank angle theta (radians).
    Returns dict of joint positions, or None if assembly fails.
    """
    L = DEFAULT_LINKS if links is None else links
    s = SCALE        if scale  is None else scale

    a  = L['a'] * s;  b  = L['b'] * s;  c  = L['c'] * s
    d_ = L['d'] * s;  e  = L['e'] * s;  f  = L['f'] * s
    g  = L['g'] * s;  h_ = L['h'] * s;  i_ = L['i'] * s
    j  = L['j'] * s;  k  = L['k'] * s
    l  = L['l'] * s;  m  = L['m'] * s

    O = np.array([0.0,  0.0])
    P = np.array([-a,  -l])

    # Crank tip — mirrored so leg extends left
    C = O + m * np.array([-np.cos(theta), np.sin(theta)])

    B_up = cci(C, j, P, b, branch='-');   
    if B_up is None: return None

    D_up = cci(B_up, c, P, d_, branch='-')
    if D_up is None: return None

    B_lo = cci(C, e, P, f, branch='+')
    if B_lo is None: return None

    D_lo = cci(B_lo, g, D_up, h_, branch='+')
    if D_lo is None: return None

    F = cci(B_lo, i_, D_lo, k, branch='+')
    if F is None: return None

    return {'O': O, 'P': P, 'C': C,
            'B_up': B_up, 'D_up': D_up,
            'B_lo': B_lo, 'D_lo': D_lo,
            'F': F}


# ── SIMULATE ONE CYCLE ────────────────────────────────────────────────────────
def simulate_cycle(n_samples=360, links=None, scale=None,
                   theta_start=0.0, theta_end=2 * np.pi):
    """
    Simulate one full crank revolution.
    Returns (thetas, foot_array, poses_list).
    """
    thetas = np.linspace(theta_start, theta_end, n_samples)
    foot   = np.full((n_samples, 2), np.nan)
    poses  = []
    for i, th in enumerate(thetas):
        pose = solve_pose(th, links, scale)
        poses.append(pose)
        if pose is not None:
            foot[i] = pose['F']
    return thetas, foot, poses


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — PLOTS & ANIMATION
# ═══════════════════════════════════════════════════════════════════════════════

# ── PLOT 1: Foot Trajectory ───────────────────────────────────────────────────
def plot_foot_trajectory(foot, save_path):
    """Static plot of the closed foot-tip trajectory (dark background)."""
    fig, ax = plt.subplots(figsize=(11, 5))
    _apply_dark_style(ax, fig)

    # Main magenta trace
    ax.plot(foot[:, 0], foot[:, 1], color=FOOT_COLOR, lw=2.2,
            label='Foot tip (F) trajectory', zorder=3)

    # Highlight stance phase (lower half of trajectory)
    y_med   = np.median(foot[:, 1])
    stance  = foot[:, 1] < y_med
    ax.plot(foot[stance, 0], foot[stance, 1],
            color='#FF6600', lw=4.0, alpha=0.6,
            label='Stance phase (ground contact)', zorder=2)

    # Start marker
    ax.plot(foot[0, 0], foot[0, 1], 'o',
            color='#00FF88', ms=11, zorder=5,
            label='Start (θ = 0°)')

    ax.set_aspect('equal')
    ax.set_xlabel('X (mm)', fontsize=12)
    ax.set_ylabel('Y (mm)', fontsize=12)
    ax.set_title('Theo Jansen — Foot Tip Trajectory  (one full crank revolution)',
                 fontsize=13)

    leg = ax.legend(loc='lower right', fontsize=10,
                    facecolor='#1A1A1A', edgecolor=FG_COLOR,
                    labelcolor=FG_COLOR)

    stride = np.ptp(foot[:, 0])
    height = np.ptp(foot[:, 1])
    info = (f'Stride length : {stride:.1f} mm\n'
            f'Step height   : {height:.1f} mm\n'
            f'Stride / step : {stride / height:.2f}')
    ax.text(0.02, 0.98, info, transform=ax.transAxes,
            va='top', fontsize=11, color=FG_COLOR,
            bbox=dict(boxstyle='round', facecolor='#1A1A2E', alpha=0.85,
                      edgecolor='#444444'))

    plt.tight_layout()
    plt.savefig(save_path, dpi=140, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"  saved: {save_path}")


# ── ANIMATION FRAME ───────────────────────────────────────────────────────────
def render_frame(pose, foot_history, ax_lim, theta_deg):
    """Render one animation frame — dark background, new link colours."""
    fig, ax = plt.subplots(figsize=(8, 8))
    _apply_dark_style(ax, fig)

    # Built-up foot trace
    if len(foot_history) > 1:
        fh = np.array(foot_history)
        ax.plot(fh[:, 0], fh[:, 1], '-', color=FOOT_COLOR,
                lw=2.0, alpha=0.90, zorder=2)

    # Draw every link with its assigned colour
    for (a, b, colour, _lbl) in LINK_EDGES:
        x = [pose[a][0], pose[b][0]]
        y = [pose[a][1], pose[b][1]]
        ax.plot(x, y, '-', color=colour, lw=5.5,
                solid_capstyle='round', zorder=3)

    # Joint markers
    joint_names = ['O', 'P', 'C', 'B_up', 'D_up', 'B_lo', 'D_lo', 'F']
    for name in joint_names:
        pt = pose[name]
        if name in ('O', 'P'):
            ax.plot(pt[0], pt[1], 's', color='white', ms=9,
                    markeredgecolor='#888888', zorder=5)
        elif name == 'F':
            ax.plot(pt[0], pt[1], 'o', color=FOOT_COLOR, ms=9,
                    markeredgecolor='white', zorder=6)
        else:
            ax.plot(pt[0], pt[1], 'o', color='white', ms=5, zorder=4)

    ax.set_xlim(ax_lim[0], ax_lim[1])
    ax.set_ylim(ax_lim[2], ax_lim[3])
    ax.set_aspect('equal')
    ax.set_xlabel('X (mm)', fontsize=11)
    ax.set_ylabel('Y (mm)', fontsize=11)
    ax.set_title(f'Theo Jansen Mechanism  —  crank angle = {theta_deg:6.1f}°',
                 fontsize=12)

    fig.canvas.draw()
    img = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()
    plt.close(fig)
    return img


def make_animation(thetas, poses, foot, save_path, n_frames=90):
    """Save animated GIF over one full crank cycle."""
    valid_poses = [p for p in poses if p is not None]
    if not valid_poses:
        print("  No valid poses — skipping animation.")
        return

    # Compute axis limits from all joint positions
    all_pts = np.vstack([list(p.values()) for p in valid_poses])
    valid_foot = foot[~np.isnan(foot[:, 0])]
    all_pts = np.vstack([all_pts, valid_foot])
    pad  = 70
    xlim = (all_pts[:, 0].min() - pad, all_pts[:, 0].max() + pad)
    ylim = (all_pts[:, 1].min() - pad, all_pts[:, 1].max() + pad)
    lim  = (xlim[0], xlim[1], ylim[0], ylim[1])

    n    = len(poses)
    step = max(1, n // n_frames)
    frames, foot_history = [], []

    for i in range(0, n, step):
        if poses[i] is None:
            continue
        foot_history.append(foot[i])
        frames.append(
            render_frame(poses[i], foot_history,
                         lim, np.degrees(thetas[i])))

    imageio.mimsave(save_path, frames, duration=0.06, loop=0)
    print(f"  saved: {save_path}  ({len(frames)} frames)")


# ── PLOT 2: Gait Comparison ───────────────────────────────────────────────────
def _reference_gait_curve(stride, height, n=200):
    """Idealised shoe-sole reference curve at given stride/height."""
    t   = np.linspace(0, 2 * np.pi, n)
    rx  = (stride / 2.0) * np.cos(t)
    ry  = np.maximum((height / 2.0) * np.sin(t), 0.05 * np.sin(2 * t))
    return rx, ry


def plot_gait_comparison(foot, save_path):
    """Overlay simulated trajectory vs an idealised reference gait curve."""
    fig, ax = plt.subplots(figsize=(10, 6))
    _apply_dark_style(ax, fig)

    fx = foot[:, 0] - foot[:, 0].mean()
    fy = foot[:, 1] - foot[:, 1].min()

    stride = np.ptp(fx);  height = np.ptp(fy)
    rx, ry = _reference_gait_curve(stride, height)
    ry     = ry - ry.min()

    ax.plot(fx, fy, color=FOOT_COLOR,  lw=2.8, label='Simulated (Jansen)')
    ax.plot(rx, ry, '--', color='#FF6600', lw=2.2, label='Reference gait curve')

    ax.set_aspect('equal')
    ax.set_xlabel('X (mm, centred)',   fontsize=12)
    ax.set_ylabel('Y (mm, ground = 0)', fontsize=12)
    ax.set_title('Simulated Jansen Trajectory  vs  Reference Gait Curve',
                 fontsize=13)
    leg = ax.legend(fontsize=11, facecolor='#1A1A1A',
                    edgecolor=FG_COLOR, labelcolor=FG_COLOR)

    plt.tight_layout()
    plt.savefig(save_path, dpi=140, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"  saved: {save_path}")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — LINK VARIATION (parameter sweep)
# ═══════════════════════════════════════════════════════════════════════════════

def sweep_parameter(param_name, values, n_samples=360):
    """Re-run kinematics for each value of one link length."""
    results = []
    for v in values:
        links = copy.deepcopy(DEFAULT_LINKS)
        links[param_name] = v
        _, foot, _ = simulate_cycle(n_samples=n_samples, links=links)
        valid = ~np.isnan(foot[:, 0])
        results.append((v, foot[valid]))
    return results


def plot_sweep(param_name, results, save_path):
    """Overlay foot trajectories for each sweep value (dark background)."""
    fig, ax = plt.subplots(figsize=(11, 5.5))
    _apply_dark_style(ax, fig)

    # Use a vivid colormap that pops on dark background
    palette = [
        '#FF4444', '#FFD700', '#39FF14',
        '#00CFFF', '#FF6EC7'
    ]

    summary = []
    for idx, ((val, foot), color) in enumerate(zip(results, palette)):
        if len(foot) == 0:
            summary.append(f'  {param_name}={val:.1f}: assembly failed')
            continue
        ax.plot(foot[:, 0], foot[:, 1], '-', color=color, lw=2.4,
                label=f'{param_name} = {val:.1f}')
        stride = np.ptp(foot[:, 0]);  height = np.ptp(foot[:, 1])
        summary.append(
            f'  {param_name}={val:.1f}: stride={stride:.1f} mm, '
            f'step={height:.1f} mm')

    ax.set_aspect('equal')
    ax.set_xlabel('X (mm)', fontsize=12)
    ax.set_ylabel('Y (mm)', fontsize=12)
    nominal = DEFAULT_LINKS[param_name]
    ax.set_title(
        f'Effect of varying  "{param_name}"  on the gait trajectory  '
        f'(nominal = {nominal:.1f})',
        fontsize=13)
    leg = ax.legend(loc='upper right', fontsize=9,
                    facecolor='#1A1A1A', edgecolor=FG_COLOR,
                    labelcolor=FG_COLOR)

    plt.tight_layout()
    plt.savefig(save_path, dpi=140, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"  saved: {save_path}")
    for line in summary:
        print(line)


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — LEGEND FIGURE (link colours reference)
# ═══════════════════════════════════════════════════════════════════════════════

def plot_link_legend(save_path):
    """Small figure showing the colour assigned to each link."""
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    fig.patch.set_facecolor(BG_COLOR)
    ax.set_facecolor(BG_COLOR)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.axis('off')
    ax.set_title('Link Colour Legend', color=FG_COLOR, fontsize=12, pad=8)

    n = len(LINK_EDGES)
    for idx, (_, _, colour, label) in enumerate(LINK_EDGES):
        y = 0.93 - idx * (0.88 / n)
        ax.plot([0.05, 0.22], [y, y], '-', color=colour, lw=5,
                solid_capstyle='round')
        ax.text(0.27, y, label, color=FG_COLOR, va='center', fontsize=9)

    plt.tight_layout()
    plt.savefig(save_path, dpi=130, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"  saved: {save_path}")


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 68)
    print("  Theo Jansen Walking Mechanism — Complete Simulation")
    print("  IE410: Introduction to Robotics")
    print("=" * 68)
    print(f"\nUsing Jansen's holy numbers  (scale × {SCALE}):")
    for k, v in DEFAULT_LINKS.items():
        print(f"   {k:1s} = {v:5.1f}  →  {v * SCALE:7.2f} mm")
    print()

    # ── Step 1: Kinematics ──────────────────────────────────────────────────
    print("─" * 50)
    print("[1/5]  Solving kinematics (360 crank angles) …")
    thetas, foot, poses = simulate_cycle(n_samples=360)
    valid = ~np.isnan(foot[:, 0])
    print(f"       Solved {valid.sum()} / {len(thetas)} poses successfully.")
    stride = np.ptp(foot[valid, 0]);  height = np.ptp(foot[valid, 1])
    print(f"       Stride length (X-span) : {stride:.1f} mm")
    print(f"       Step height   (Y-span)  : {height:.1f} mm")
    print(f"       Stride / step ratio     : {stride / height:.2f}")

    # ── Step 2: Foot trajectory plot ────────────────────────────────────────
    print("\n─" * 50)
    print("[2/5]  Plotting foot trajectory …")
    plot_foot_trajectory(
        foot[valid],
        os.path.join(OUT_DIR, 'foot_trajectory.png'))

    # ── Step 3: Animation ───────────────────────────────────────────────────
    print("\n─" * 50)
    print("[3/5]  Rendering mechanism animation …")
    make_animation(
        thetas, poses, foot,
        os.path.join(OUT_DIR, 'mechanism_animation.gif'),
        n_frames=90)

    # ── Step 4: Gait comparison ─────────────────────────────────────────────
    print("\n─" * 50)
    print("[4/5]  Plotting gait comparison …")
    plot_gait_comparison(
        foot[valid],
        os.path.join(OUT_DIR, 'gait_comparison.png'))

    # ── Step 5: Link variation ──────────────────────────────────────────────
    print("\n─" * 50)
    print("[5/5]  Parameter sweep: varying m (crank length) …")
    res_m = sweep_parameter('m', [12.0, 13.5, 15.0, 16.5, 18.0])
    plot_sweep('m', res_m, os.path.join(OUT_DIR, 'link_variation_m.png'))

    print()
    print("        Parameter sweep: varying h (long connector) …")
    res_h = sweep_parameter('h', [62.0, 64.0, 65.7, 67.5, 69.5])
    plot_sweep('h', res_h, os.path.join(OUT_DIR, 'link_variation_h.png'))

    # ── Bonus: colour legend ────────────────────────────────────────────────
    print("\n─" * 50)
    print("[+]    Generating link colour legend …")
    plot_link_legend(os.path.join(OUT_DIR, 'link_colour_legend.png'))

    print("\n" + "=" * 68)
    print("  All outputs saved to:", OUT_DIR)
    print("=" * 68)


if __name__ == '__main__':
    main()
