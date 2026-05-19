"""
⚡ Generator Design Explorer  v3
==================================
Ring of cylindrical magnets + two positionable coils.

FIXES IN THIS VERSION
---------------------
  • Radial vs Axial actually show different behaviour:
      - Switching snaps coils to their optimal positions/orientations
      - Radial: coils face inward at ring equator
      - Axial:  coils face upward above/below the ring
  • Real physical units on all graphs:
      - Flux in nanoweber  [nWb]   (1 nWb = 1e-9 T·m²)
      - EMF  in millivolt  [mV]    (requires RPM + N_turns)
  • RPM and N_turns sliders so you design for real voltages

PANELS
------
  Top-left   3D animated ring + two coils (glow = current)
  Top-right  B-field XY cross-section at Coil 1 height
  Mid-right  Flux Φ  [nWb] for each coil — one full revolution
  Bot-right  Induced EMF [mV] + combined wired output

CONTROLS
--------
  Coil 1  X / Y / Z / Theta / Phi  — free 3D positioning
  Coil 2  ring-angle / Z offset    — auto-tracks ring surface
  Coil radius  (shared)
  Magnetization toggle  Radial | Axial
  Wiring toggle         Series | Parallel
  RPM slider            rotor speed (physical + visual)
  N turns slider        coil turns (scales EMF linearly)
  Magnets N             4–12 magnets in the ring

REQUIREMENTS
------------
  pip install magpylib matplotlib numpy scipy
"""

import numpy as np
import matplotlib
for _b in ('TkAgg', 'Qt5Agg', 'WXAgg', 'macosx', 'GTK3Agg'):
    try: matplotlib.use(_b); break
    except Exception: pass

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.animation as manim
from matplotlib.widgets import Slider, Button, RadioButtons
from mpl_toolkits.mplot3d import Axes3D                   # noqa: F401
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import magpylib as magpy
from scipy.interpolate import interp1d
import warnings
warnings.filterwarnings('ignore')

# ── palette ────────────────────────────────────────────────────────────────────
BG    = '#0d1117';  PANEL = '#161b22'
CN    = '#58a6ff';  CS    = '#f78166'   # North blue / South red
CC1   = '#e3b341';  CC2   = '#a5d6ff'   # Coil 1 amber / Coil 2 sky
CEMF  = '#3fb950';  CNEG  = '#f78166'   # EMF+ green / EMF- red
CFLX  = '#d2a679'   # flux gold
TEXT  = '#e6edf3';  GRID  = '#21262d';  SLID = '#30363d'

plt.rcParams.update({
    'figure.facecolor': BG, 'axes.facecolor': BG,
    'axes.edgecolor': GRID, 'axes.labelcolor': TEXT,
    'xtick.color': TEXT, 'ytick.color': TEXT,
    'text.color': TEXT, 'font.family': 'monospace',
    'grid.color': GRID, 'grid.alpha': 0.35,
})

# ── unit conversion ────────────────────────────────────────────────────────────
# magpylib getB returns mT; distances in mm
# flux = mean(B·n̂) * area  →  units: mT · mm²  =  1e-3 T · 1e-6 m²  =  1e-9 Wb = 1 nWb
# EMF  = -N · dΦ/dt  =  -N · ω · dΦ/dθ
#      = N · ω [rad/s] · |dΦ/dθ [nWb/rad]| · 1e-9   →  Volts
#      multiply by 1e3 for millivolts
NWB_PER_UNIT = 1.0          # flux already in nWb after mT·mm² → nWb
MV_SCALE     = 1e-6         # nWb/rad → mV:  N · ω · 1e-9 · 1e3 = N · ω · 1e-6


# ══════════════════════════════════════════════════════════════════════════════
# GEOMETRY HELPERS
# ══════════════════════════════════════════════════════════════════════════════

def cylinder_verts(cx, cy, cz, r, h, n=10):
    th = np.linspace(0, 2*np.pi, n, endpoint=False)
    xt, yt = cx + r*np.cos(th), cy + r*np.sin(th)
    zt, zb = cz + h/2, cz - h/2
    v = []
    for i in range(n):
        j = (i+1) % n
        v.append([[xt[i],yt[i],zb],[xt[j],yt[j],zb],
                  [xt[j],yt[j],zt],[xt[i],yt[i],zt]])
    v.append([[xt[i],yt[i],zt] for i in range(n)])
    v.append([[xt[i],yt[i],zb] for i in range(n)])
    return v

def coil_xyz(cx, cy, cz, r, normal, n=80):
    nv = np.array(normal, dtype=float); nv /= np.linalg.norm(nv)
    p1 = np.cross(nv, [0,0,1])
    if np.linalg.norm(p1) < 0.01: p1 = np.cross(nv, [1,0,0])
    p1 /= np.linalg.norm(p1); p2 = np.cross(nv, p1)
    th = np.linspace(0, 2*np.pi, n)
    pts = np.array([cx,cy,cz]) + r*(np.outer(np.cos(th),p1)+np.outer(np.sin(th),p2))
    return pts[:,0], pts[:,1], pts[:,2]

def spherical_normal(theta_deg, phi_deg):
    th, ph = np.deg2rad(theta_deg), np.deg2rad(phi_deg)
    return np.array([np.sin(th)*np.cos(ph), np.sin(th)*np.sin(ph), np.cos(th)])


# ══════════════════════════════════════════════════════════════════════════════
# PHYSICS ENGINE
# ══════════════════════════════════════════════════════════════════════════════

class Engine:
    STRENGTH = 1_200_000   # A/m  (≈ N42 NdFeB)
    RING_R   = 50.0        # mm
    MAG_DIA  = 12.0        # mm diameter
    MAG_H    = 10.0        # mm height

    def __init__(self):
        self.n_magnets = 8
        self.radial    = True      # True = radial (in/out),  False = axial (up/down)

    def build_ring(self, rotor_angle=0.0):
        mags = []
        for i in range(self.n_magnets):
            pa   = rotor_angle + i * 360.0 / self.n_magnets
            rad  = np.deg2rad(pa)
            pos  = np.array([self.RING_R*np.cos(rad), self.RING_R*np.sin(rad), 0.0])
            sign = 1 if i % 2 == 0 else -1
            if self.radial:
                # magnetization points radially outward/inward (alternating)
                mag_dir = sign * np.array([np.cos(rad), np.sin(rad), 0.0])
            else:
                # magnetization points along Z axis (alternating up/down)
                mag_dir = sign * np.array([0.0, 0.0, 1.0])
            m = magpy.magnet.Cylinder(
                magnetization = mag_dir * self.STRENGTH,
                dimension     = (self.MAG_DIA, self.MAG_H),
                position      = pos)
            if self.radial:
                m.rotate_from_angax(pa, 'z')   # align flat face outward
            mags.append(m)
        return magpy.Collection(*mags)

    def mag_info(self, rotor_angle=0.0):
        info = []
        for i in range(self.n_magnets):
            pa  = rotor_angle + i * 360.0 / self.n_magnets
            rad = np.deg2rad(pa)
            info.append((self.RING_R*np.cos(rad), self.RING_R*np.sin(rad), 0.0, i%2==0))
        return info

    def flux_nwb(self, ring, center_mm, normal, coil_r_mm):
        """Return flux in nanoweber [nWb]."""
        nv  = np.array(normal, dtype=float); nv /= np.linalg.norm(nv)
        p1  = np.cross(nv, [0,0,1])
        if np.linalg.norm(p1) < 0.01: p1 = np.cross(nv, [1,0,0])
        p1 /= np.linalg.norm(p1); p2 = np.cross(nv, p1)
        cen = np.array(center_mm, dtype=float)
        pts = [cen]
        for ri in range(1, 6):
            rd = coil_r_mm * ri / 5
            for ci in range(16):
                phi = 2*np.pi*ci/16
                pts.append(cen + rd*(np.cos(phi)*p1 + np.sin(phi)*p2))
        pts = np.array(pts)
        B = magpy.getB(ring, pts, squeeze=True)   # mT
        # mT · mm²  = 1e-3 T · 1e-6 m²  = 1e-9 Wb = 1 nWb
        return float((B @ nv).mean() * np.pi * coil_r_mm**2)

    def field_slice(self, ring, z_mm=0, half=85, res=18):
        xs = np.linspace(-half, half, res)
        gx, gy = np.meshgrid(xs, xs)
        pts = np.column_stack([gx.ravel(), gy.ravel(), np.full(res*res, z_mm)])
        B   = magpy.getB(ring, pts, squeeze=True)
        return gx, gy, B[:,0].reshape(gx.shape), B[:,1].reshape(gy.shape)

    def waveform(self, c1_cen, c1_norm, c2_cen, c2_norm, coil_r, n=120):
        """
        Returns angles [deg], flux1 [nWb], flux2 [nWb],
                emf1_per_rad [nWb/rad], emf2_per_rad [nWb/rad]
        Multiply EMF by  N * omega * 1e-6  to get millivolts.
        """
        print("  Computing waveforms... ", end='', flush=True)
        angles = np.linspace(0, 360, n, endpoint=False)
        f1, f2 = [], []
        for a in angles:
            ring = self.build_ring(a)
            f1.append(self.flux_nwb(ring, c1_cen, c1_norm, coil_r))
            f2.append(self.flux_nwb(ring, c2_cen, c2_norm, coil_r))
        f1, f2 = np.array(f1), np.array(f2)
        da = np.deg2rad(angles[1] - angles[0])
        e1 = -np.gradient(f1, da)
        e2 = -np.gradient(f2, da)
        print("done.")
        return angles, f1, f2, e1, e2


# ══════════════════════════════════════════════════════════════════════════════
# OPTIMAL COIL DEFAULTS PER MAGNETIZATION TYPE
# ══════════════════════════════════════════════════════════════════════════════

def optimal_c1(radial, ring_r=50, coil_r=15):
    """Return (cx,cy,cz,theta,phi) for coil 1 optimal placement."""
    if radial:
        # At ring equator, just outside, normal pointing inward (-x direction)
        d = ring_r + coil_r * 0.5 + 4
        return dict(cx=d, cy=0.0, cz=0.0, theta=90.0, phi=180.0)
    else:
        # Above ring plane, over the edge, normal pointing up
        return dict(cx=ring_r, cy=0.0, cz=ring_r*0.35, theta=0.0, phi=0.0)

def optimal_c2_pos(ring_angle_deg, radial, ring_r=50, coil_r=15, cz=0.0):
    """Return (center_mm, normal) for coil 2 at given ring angle."""
    ang = np.deg2rad(ring_angle_deg)
    if radial:
        d   = ring_r + coil_r * 0.5 + 4
        cen = [d*np.cos(ang), d*np.sin(ang), cz]
        norm = [-np.cos(ang), -np.sin(ang), 0.0]
    else:
        cen  = [ring_r*np.cos(ang), ring_r*np.sin(ang), ring_r*0.35 + cz]
        norm = [0.0, 0.0, 1.0]
    return cen, norm


# ══════════════════════════════════════════════════════════════════════════════
# VISUALISER
# ══════════════════════════════════════════════════════════════════════════════

class Viz:
    MS = 100

    def __init__(self):
        self.eng     = Engine()
        self.angle   = 0.0
        self.rpm     = 300.0
        self.n_turns = 200
        self.series  = True
        self.coil_r  = 15.0

        # Coil 1 free state
        self.c1 = dict(optimal_c1(True, self.eng.RING_R, self.coil_r))
        # Coil 2 parametric state
        self.c2_ring_angle = 0.0
        self.c2_z = 0.0

        print("\n⚡  Generator Design Explorer — initialising...")
        self._recalc()
        self._build_ui()

    # ── coil geometry helpers ─────────────────────────────────────────────────

    def _c1_cen_norm(self):
        cen  = [self.c1['cx'], self.c1['cy'], self.c1['cz']]
        norm = spherical_normal(self.c1['theta'], self.c1['phi'])
        return cen, norm

    def _c2_cen_norm(self):
        return optimal_c2_pos(self.c2_ring_angle, self.eng.radial,
                               self.eng.RING_R, self.coil_r, self.c2_z)

    def _omega(self):
        return self.rpm * 2 * np.pi / 60.0

    def _emf_mv(self, emf_per_rad):
        """Convert nWb/rad → mV using current N_turns and RPM."""
        return emf_per_rad * self.n_turns * self._omega() * MV_SCALE

    # ── recalculate waveforms ─────────────────────────────────────────────────

    def _recalc(self):
        c1c, c1n = self._c1_cen_norm()
        c2c, c2n = self._c2_cen_norm()
        (self._wv_ang,
         self._wv_f1, self._wv_f2,
         self._wv_e1_pr, self._wv_e2_pr) = self.eng.waveform(
            c1c, c1n, c2c, c2n, self.coil_r, n=120)
        self._update_combined()

    def _update_combined(self):
        e1_mv = self._emf_mv(self._wv_e1_pr)
        e2_mv = self._emf_mv(self._wv_e2_pr)
        if self.series:
            self._wv_comb = e1_mv + e2_mv
            self._comb_lbl = 'Series: V_out = EMF1 + EMF2'
        else:
            self._wv_comb = (e1_mv + e2_mv) / 2
            self._comb_lbl = 'Parallel: V_out = avg(EMF1,EMF2)  [×2 current]'
        self._wv_e1_mv = e1_mv
        self._wv_e2_mv = e2_mv
        self._emf_max  = max(abs(self._wv_comb).max(), 1e-12)
        ang = self._wv_ang
        self._e1_itp   = interp1d(ang, e1_mv,            kind='cubic', fill_value='extrapolate')
        self._e2_itp   = interp1d(ang, e2_mv,            kind='cubic', fill_value='extrapolate')
        self._comb_itp = interp1d(ang, self._wv_comb,    kind='cubic', fill_value='extrapolate')
        self._f1_itp   = interp1d(ang, self._wv_f1,      kind='cubic', fill_value='extrapolate')
        self._f2_itp   = interp1d(ang, self._wv_f2,      kind='cubic', fill_value='extrapolate')

    # ── layout ────────────────────────────────────────────────────────────────

    def _build_ui(self):
        self.fig = plt.figure(figsize=(19, 11), facecolor=BG)
        try: self.fig.canvas.manager.set_window_title('Generator Design Explorer')
        except Exception: pass

        gs = gridspec.GridSpec(3, 3, figure=self.fig,
                               left=0.04, right=0.97, top=0.93, bottom=0.30,
                               wspace=0.38, hspace=0.55)

        # 3D view
        self.ax3 = self.fig.add_subplot(gs[:, 0], projection='3d')
        self.ax3.set_facecolor(BG); self.ax3.tick_params(labelsize=6)
        lim = 90
        self.ax3.set_xlim(-lim,lim); self.ax3.set_ylim(-lim,lim); self.ax3.set_zlim(-lim,lim)
        self.ax3.set_xlabel('X (mm)',fontsize=7,labelpad=1)
        self.ax3.set_ylabel('Y (mm)',fontsize=7,labelpad=1)
        self.ax3.set_zlabel('Z (mm)',fontsize=7,labelpad=1)
        self.ax3.grid(True, alpha=0.2)

        # Right panels
        self.ax_b  = self.fig.add_subplot(gs[0, 1:])
        self.ax_f  = self.fig.add_subplot(gs[1, 1:])
        self.ax_e  = self.fig.add_subplot(gs[2, 1:])

        self.ax_b.set_facecolor(PANEL); self.ax_b.set_aspect('equal')
        self.ax_b.set_xlabel('X (mm)',fontsize=8); self.ax_b.set_ylabel('Y (mm)',fontsize=8)
        self.ax_b.tick_params(labelsize=7)

        self.ax_f.set_facecolor(PANEL); self.ax_f.tick_params(labelsize=7)
        self.ax_f.set_xlabel('Rotor angle (°)',fontsize=8)
        self.ax_f.set_ylabel('Flux  Φ  [nWb]', color=CFLX, fontsize=8)
        self.ax_f.grid(True,alpha=0.3); self.ax_f.set_xlim(0,360)

        self.ax_e.set_facecolor(PANEL); self.ax_e.tick_params(labelsize=7)
        self.ax_e.set_xlabel('Rotor angle (°)',fontsize=8)
        self.ax_e.set_ylabel('EMF  [mV]', color=CEMF, fontsize=8)
        self.ax_e.grid(True,alpha=0.3); self.ax_e.set_xlim(0,360)
        self.ax_e.axhline(0, color=TEXT, lw=0.6, alpha=0.4)

        # ── sliders ───────────────────────────────────────────────────────────
        def msl(l, b, label, lo, hi, init, step=None, w=0.13):
            ax = self.fig.add_axes([l, b, w, 0.019]); ax.set_facecolor(PANEL)
            kw = dict(color=SLID, track_color=GRID)
            if step: kw['valstep'] = step
            sl = Slider(ax, label, lo, hi, valinit=init, **kw)
            sl.label.set_color(TEXT); sl.label.set_fontsize(7.5)
            sl.valtext.set_color(CN);  sl.valtext.set_fontsize(7.5)
            return sl

        # Row 1: coil 1 position
        self.sl_c1x  = msl(0.04, 0.265, 'C1 X mm',  -80,  80, self.c1['cx'])
        self.sl_c1y  = msl(0.19, 0.265, 'C1 Y mm',  -80,  80, self.c1['cy'])
        self.sl_c1z  = msl(0.34, 0.265, 'C1 Z mm',  -70,  70, self.c1['cz'])
        # Row 2: coil 1 orientation
        self.sl_c1th = msl(0.04, 0.243, 'C1 θ tilt',  0, 180, self.c1['theta'])
        self.sl_c1ph = msl(0.19, 0.243, 'C1 φ rot',   0, 360, self.c1['phi'])
        self.sl_cr   = msl(0.34, 0.243, 'Coil R mm',  5,  40, self.coil_r)
        # Row 3: coil 2
        self.sl_c2a  = msl(0.04, 0.221, 'C2 angle °', 0, 360, 0.0)
        self.sl_c2z  = msl(0.19, 0.221, 'C2 Z mm',  -70,  70, 0.0)
        # Row 4: rotor + output
        self.sl_rpm  = msl(0.04, 0.199, 'RPM',        10, 3000, self.rpm)
        self.sl_nt   = msl(0.19, 0.199, 'N turns',    10, 2000, self.n_turns, step=10)
        self.sl_nm   = msl(0.34, 0.199, 'Magnets N',   4,   12, 8, step=2)

        for sl, cb in [
            (self.sl_c1x, self._cb_c1),(self.sl_c1y, self._cb_c1),(self.sl_c1z, self._cb_c1),
            (self.sl_c1th,self._cb_c1),(self.sl_c1ph,self._cb_c1),(self.sl_cr,  self._cb_c1),
            (self.sl_c2a, self._cb_c2),(self.sl_c2z, self._cb_c2),
            (self.sl_rpm, self._cb_rpm),(self.sl_nt,  self._cb_turns),
            (self.sl_nm,  self._cb_nm),
        ]:
            sl.on_changed(cb)

        # ── radio buttons ─────────────────────────────────────────────────────
        ax_mag = self.fig.add_axes([0.50, 0.193, 0.11, 0.055])
        ax_mag.set_facecolor(PANEL)
        self.rb_mag = RadioButtons(ax_mag, ('Radial', 'Axial'), activecolor=CN)
        for lbl in self.rb_mag.labels:
            lbl.set_color(TEXT); lbl.set_fontsize(8)
        self.rb_mag.on_clicked(self._cb_mag)

        ax_wir = self.fig.add_axes([0.63, 0.193, 0.12, 0.055])
        ax_wir.set_facecolor(PANEL)
        self.rb_wir = RadioButtons(ax_wir, ('Series','Parallel'), activecolor=CEMF)
        for lbl in self.rb_wir.labels:
            lbl.set_color(TEXT); lbl.set_fontsize(8)
        self.rb_wir.on_clicked(self._cb_wiring)

        # ── buttons ───────────────────────────────────────────────────────────
        ax_pp = self.fig.add_axes([0.77, 0.200, 0.065, 0.038])
        ax_rs = self.fig.add_axes([0.84, 0.200, 0.065, 0.038])
        self.btn_pp = Button(ax_pp,'Pause',color=PANEL,hovercolor=SLID)
        self.btn_rs = Button(ax_rs,'Reset',color=PANEL,hovercolor=SLID)
        self.btn_pp.label.set_color(TEXT); self.btn_rs.label.set_color(TEXT)
        self.btn_pp.on_clicked(self._toggle_play)
        self.btn_rs.on_clicked(self._reset)

        # ── header ────────────────────────────────────────────────────────────
        self.fig.text(0.5, 0.97,
            'GENERATOR DESIGN EXPLORER  —  Radial vs Axial  |  Dual Coil  |  Real Units [nWb / mV]',
            ha='center', fontsize=12, color=CN, fontfamily='monospace', fontweight='bold')

        # Column labels
        self.fig.text(0.04, 0.285, 'COIL 1  (amber) — free position', fontsize=8, color=CC1)
        self.fig.text(0.04, 0.273, 'COIL 2  (blue)  — auto on ring surface', fontsize=8, color=CC2)

        for x, c, lbl in [(0.50,CN,'N pole'),(0.55,CS,'S pole'),
                           (0.61,CC1,'Coil 1'),(0.67,CC2,'Coil 2'),(0.73,CEMF,'I+/EMF+')]:
            self.fig.text(x, 0.285, '●', color=c, fontsize=11)
            self.fig.text(x+0.012, 0.285, lbl, color=TEXT, fontsize=8)

        # ── note about switching ───────────────────────────────────────────────
        self.fig.text(0.50, 0.255,
            'Radial/Axial switch: snaps coils to optimal positions for each type',
            fontsize=7.5, color='#666e79')
        self.fig.text(0.50, 0.243,
            'RPM + N turns scale the mV output without recomputing field',
            fontsize=7.5, color='#666e79')

        # ── initialise background waveforms ───────────────────────────────────
        self._wv_artists = {}
        self._wv_fill    = []
        self._dyn        = {}
        self._refresh_waveform_bg()

    # ── waveform background ────────────────────────────────────────────────────

    def _refresh_waveform_bg(self):
        for k, v in list(self._wv_artists.items()):
            try: v.remove()
            except Exception: pass
        for fc in self._wv_fill:
            try: fc.remove()
            except Exception: pass
        self._wv_artists = {}; self._wv_fill = []

        a   = self._wv_ang
        f1  = self._wv_f1;  f2  = self._wv_f2
        e1  = self._wv_e1_mv; e2 = self._wv_e2_mv; ec = self._wv_comb

        self._wv_artists['f1'], = self.ax_f.plot(a, f1, color=CC1,  lw=1.5, alpha=0.5, zorder=1)
        self._wv_artists['f2'], = self.ax_f.plot(a, f2, color=CC2,  lw=1.5, alpha=0.5, zorder=1)
        self.ax_f.relim(); self.ax_f.autoscale_view()

        self._wv_artists['e1'], = self.ax_e.plot(a, e1, color=CC1,  lw=1.0, alpha=0.5, zorder=1, ls='--', label='EMF1')
        self._wv_artists['e2'], = self.ax_e.plot(a, e2, color=CC2,  lw=1.0, alpha=0.5, zorder=1, ls='--', label='EMF2')
        self._wv_artists['ec'], = self.ax_e.plot(a, ec, color=CEMF, lw=2.5, alpha=0.95, zorder=2, label='Combined')
        self.ax_e.relim(); self.ax_e.autoscale_view()
        self.ax_e.set_title(self._comb_lbl, color=TEXT, fontsize=8.5)

        peak = abs(ec).max()
        mag_str = 'Radial' if self.eng.radial else 'Axial'
        self.ax_f.set_title(f'Flux Φ per coil  [{mag_str} magnetization]', color=TEXT, fontsize=9)

        self._wv_fill = [
            self.ax_e.fill_between(a, 0, ec, where=ec>=0, color=CEMF, alpha=0.10, zorder=0),
            self.ax_e.fill_between(a, 0, ec, where=ec<0,  color=CNEG, alpha=0.10, zorder=0),
        ]
        # Peak annotation
        if 'pk' in self._wv_artists:
            try: self._wv_artists['pk'].remove()
            except Exception: pass
        pk_idx = np.argmax(abs(ec))
        self._wv_artists['pk'] = self.ax_e.annotate(
            f'  pk {peak:.2f} mV',
            xy=(a[pk_idx], ec[pk_idx]),
            color=CEMF, fontsize=8, zorder=6)

    # ── artist lifecycle ───────────────────────────────────────────────────────

    def _rm(self, *keys):
        for k in keys:
            obj = self._dyn.pop(k, None)
            if obj is None: continue
            for o in (obj if isinstance(obj,(list,tuple)) else [obj]):
                try: o.remove()
                except Exception: pass

    # ── 3D scene ──────────────────────────────────────────────────────────────

    def _draw_3d(self, ring_angle):
        self._rm('3d_mags','3d_c1','3d_c2','3d_axis','3d_guide',
                 '3d_c1n','3d_c2n','3d_glow1','3d_glow2','3d_info')
        ax = self.ax3

        # Guide ring
        th = np.linspace(0,2*np.pi,120)
        gl, = ax.plot(self.eng.RING_R*np.cos(th), self.eng.RING_R*np.sin(th),
                      np.zeros(120), '-', color='#2a3038', lw=0.9, alpha=0.5)
        self._dyn['3d_guide'] = gl
        axl,= ax.plot([0,0],[0,0],[-80,80],'--',color='#222a32',lw=1.0,alpha=0.5)
        self._dyn['3d_axis'] = axl

        # Magnets
        polys = []
        for (mx,my,mz,isn) in self.eng.mag_info(ring_angle):
            verts = cylinder_verts(mx,my,mz, self.eng.MAG_DIA/2, self.eng.MAG_H)
            p = Poly3DCollection(verts, alpha=0.88,
                                 facecolor=(CN if isn else CS),
                                 edgecolor='#0d1117', linewidth=0.3, zorder=3)
            ax.add_collection3d(p); polys.append(p)
        self._dyn['3d_mags'] = polys

        # Coil 1
        c1c, c1n = self._c1_cen_norm()
        a1 = []
        for dz in (-1.2, 0, 1.2):
            xs,ys,zs = coil_xyz(c1c[0],c1c[1],c1c[2]+dz, self.coil_r, c1n)
            l, = ax.plot(xs,ys,zs, '-', color=CC1, lw=1.8, alpha=0.9, zorder=4)
            a1.append(l)
        nv1 = np.array(c1n)*self.coil_r*1.5
        ln, = ax.plot([c1c[0],c1c[0]+nv1[0]],[c1c[1],c1c[1]+nv1[1]],[c1c[2],c1c[2]+nv1[2]],
                      '-', color=CC1, lw=2, zorder=5)
        a1.append(ln)
        self._dyn['3d_c1'] = a1

        # Coil 2
        c2c, c2n = self._c2_cen_norm()
        a2 = []
        for dz in (-1.2, 0, 1.2):
            xs,ys,zs = coil_xyz(c2c[0],c2c[1],c2c[2]+dz, self.coil_r, c2n)
            l, = ax.plot(xs,ys,zs, '-', color=CC2, lw=1.8, alpha=0.9, zorder=4)
            a2.append(l)
        nv2 = np.array(c2n)*self.coil_r*1.5
        ln2,= ax.plot([c2c[0],c2c[0]+nv2[0]],[c2c[1],c2c[1]+nv2[1]],[c2c[2],c2c[2]+nv2[2]],
                      '-', color=CC2, lw=2, zorder=5)
        a2.append(ln2)
        self._dyn['3d_c2'] = a2

        mag_lbl = 'Radial' if self.eng.radial else 'Axial'
        ax.set_title(f'{mag_lbl} ring  |  θ={ring_angle:.1f}°  |  {self.eng.n_magnets} magnets  '
                     f'|  {self.rpm:.0f} RPM  |  N={self.n_turns}',
                     color=TEXT, fontsize=8.5, pad=3)

    # ── B-field slice ──────────────────────────────────────────────────────────

    def _draw_bfield(self, ring):
        self._rm('bf_q','bf_c1','bf_c2')
        ax = self.ax_b
        c1c, c1n = self._c1_cen_norm()
        gx, gy, Bx, By = self.eng.field_slice(ring, z_mm=c1c[2], half=85, res=18)
        Bm = np.hypot(Bx,By); Bm[Bm==0]=1e-12
        norm = np.log1p(Bm)/(np.log1p(Bm).max()+1e-12)
        self._dyn['bf_q'] = ax.quiver(gx,gy,Bx/Bm,By/Bm,norm,
                                       cmap='cool',alpha=0.80,scale=30,
                                       width=0.004,pivot='mid',zorder=2)
        # coil footprints
        th = np.linspace(0,2*np.pi,80)
        for (cc, cn, col, key) in [(c1c,c1n,CC1,'bf_c1'),(self._c2_cen_norm()[0],self._c2_cen_norm()[1],CC2,'bf_c2')]:
            p1v = np.cross(cn,[0,0,1])
            if np.linalg.norm(p1v)<0.01: p1v=np.cross(cn,[1,0,0])
            p1v /= np.linalg.norm(p1v); p2v = np.cross(cn,p1v)
            xy = np.array([[cc[0]+self.coil_r*(np.cos(t)*p1v[0]+np.sin(t)*p2v[0]),
                             cc[1]+self.coil_r*(np.cos(t)*p1v[1]+np.sin(t)*p2v[1])] for t in th])
            lc, = ax.plot(xy[:,0],xy[:,1],'-',color=col,lw=2.0,alpha=0.9,zorder=5)
            ax.plot(cc[0],cc[1],'o',color=col,ms=5,zorder=6)
            self._dyn[key] = lc
        ax.set_title(f'B-field XY slice at Z = {c1c[2]:.0f} mm', color=TEXT, fontsize=9)

    # ── waveform cursors ───────────────────────────────────────────────────────

    def _draw_cursors(self, angle):
        self._rm('wc_f1','wc_f2','wd_f1','wd_f2',
                 'wc_e1','wc_e2','wc_ec','wd_e1','wd_e2','wd_ec','w_lbl')
        a    = angle % 360
        f1v  = float(self._f1_itp(a))
        f2v  = float(self._f2_itp(a))
        e1v  = float(self._e1_itp(a))
        e2v  = float(self._e2_itp(a))
        cv   = float(self._comb_itp(a))
        norm = float(np.clip(cv / self._emf_max, -1, 1))

        for ax, lim_fn, cursors_dots in [
            (self.ax_f, lambda: self.ax_f.get_ylim(),
             [('wc_f1',CC1,[a,a],None),('wc_f2',CC2,[a,a],None),
              ('wd_f1',CC1,a,f1v),('wd_f2',CC2,a,f2v)]),
            (self.ax_e, lambda: self.ax_e.get_ylim(),
             [('wc_e1',CC1,[a,a],None),('wc_e2',CC2,[a,a],None),('wc_ec',CEMF,[a,a],None),
              ('wd_e1',CC1,a,e1v),('wd_e2',CC2,a,e2v),('wd_ec',CEMF,a,cv)]),
        ]:
            ylim = lim_fn()
            for key, col, xv, yv in cursors_dots:
                if yv is None:  # vertical cursor
                    self._dyn[key], = ax.plot(xv, list(ylim), color=col, lw=0.9, alpha=0.8, zorder=3)
                else:           # dot
                    self._dyn[key], = ax.plot(xv, yv, 'o', color=col, ms=7, zorder=4)

        col = CEMF if norm >= 0 else CNEG
        lbl = self.ax_e.text(a+4, cv, f' {cv:.2f} mV  {">>>" if norm>=0 else "<<<"}',
                              color=col, fontsize=8, va='center', zorder=5)
        self._dyn['w_lbl'] = lbl
        return norm

    # ── 3D current glow ────────────────────────────────────────────────────────

    def _draw_glow(self, norm):
        self._rm('3d_glow1','3d_glow2','3d_info')
        col   = CEMF if norm >= 0 else CNEG
        alpha = 0.35 + 0.65 * abs(norm)
        lw    = 1.2 + 4.5 * abs(norm)
        c1c,c1n = self._c1_cen_norm()
        c2c,c2n = self._c2_cen_norm()
        for cen, norm_v, key in [(c1c,c1n,'3d_glow1'),(c2c,c2n,'3d_glow2')]:
            xs,ys,zs = coil_xyz(cen[0],cen[1],cen[2], self.coil_r, norm_v)
            g, = self.ax3.plot(xs,ys,zs, '-', color=col, lw=lw, alpha=alpha, zorder=7)
            self._dyn[key] = g
        wir = 'Series' if self.series else 'Parallel'
        mag = 'Radial' if self.eng.radial else 'Axial'
        cv  = float(self._comb_itp(self.angle % 360))
        t = self.ax3.text2D(0.02, 0.04,
            f'V_out = {cv:.2f} mV  {">>>" if norm>=0 else "<<<"}\n{mag} | {wir}',
            transform=self.ax3.transAxes, color=col, fontsize=9, fontweight='bold')
        self._dyn['3d_info'] = t

    # ── animation ─────────────────────────────────────────────────────────────

    def _animate(self, _f):
        if self.playing:
            # advance angle proportional to RPM (visual speed)
            deg_per_frame = self.rpm / 60.0 * (self.MS/1000.0) * 360.0
            self.angle = (self.angle + max(deg_per_frame, 0.5)) % 360
        ring = self.eng.build_ring(self.angle)
        self._draw_3d(self.angle)
        self._draw_bfield(ring)
        norm = self._draw_cursors(self.angle)
        self._draw_glow(norm)
        return []

    # ── callbacks ─────────────────────────────────────────────────────────────

    def _cb_c1(self, _v):
        self.c1 = dict(cx=self.sl_c1x.val, cy=self.sl_c1y.val, cz=self.sl_c1z.val,
                       theta=self.sl_c1th.val, phi=self.sl_c1ph.val)
        self.coil_r = self.sl_cr.val
        self._recalc(); self._refresh_waveform_bg()

    def _cb_c2(self, _v):
        self.c2_ring_angle = self.sl_c2a.val
        self.c2_z          = self.sl_c2z.val
        self._recalc(); self._refresh_waveform_bg()

    def _cb_rpm(self, v):
        self.rpm = v
        self._update_combined(); self._refresh_waveform_bg()

    def _cb_turns(self, v):
        self.n_turns = int(v)
        self._update_combined(); self._refresh_waveform_bg()

    def _cb_nm(self, v):
        self.eng.n_magnets = int(v)
        self._recalc(); self._refresh_waveform_bg()

    def _cb_mag(self, label):
        """Switch magnetization type AND snap coils to optimal positions."""
        self.eng.radial = (label == 'Radial')
        # Snap both coils to their optimal positions for this type
        opt = optimal_c1(self.eng.radial, self.eng.RING_R, self.coil_r)
        self.c1 = dict(opt)
        # Update sliders silently (disconnect first to avoid re-triggering)
        for sl, key in [(self.sl_c1x,'cx'),(self.sl_c1y,'cy'),(self.sl_c1z,'cz'),
                        (self.sl_c1th,'theta'),(self.sl_c1ph,'phi')]:
            sl.set_val(opt[key])
        self._recalc(); self._refresh_waveform_bg()

    def _cb_wiring(self, label):
        self.series = (label == 'Series')
        self._update_combined(); self._refresh_waveform_bg()

    def _toggle_play(self, _e):
        self.playing = not self.playing
        self.btn_pp.label.set_text('Play' if not self.playing else 'Pause')

    def _reset(self, _e):
        self.angle = 0.0

    def run(self):
        self.playing = True
        self.anim = manim.FuncAnimation(
            self.fig, self._animate, interval=self.MS,
            blit=False, cache_frame_data=False)
        print("  Window open. RPM/N_turns sliders update mV instantly (no recompute).\n"
              "  Coil position/magnet changes recompute the field (~10s).\n")
        plt.show()


# ══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    print("""
+================================================================+
|  GENERATOR DESIGN EXPLORER  v3                                 |
+================================================================+
|  WHAT'S FIXED                                                  |
|    Radial vs Axial actually differ — switching snaps coils     |
|    to optimal positions for each magnetization type            |
|    Graphs show real units: Flux [nWb], EMF [mV]               |
|    RPM + N_turns sliders scale mV output without recomputing   |
|                                                                |
|  KEY THINGS TO TRY                                             |
|    1. Hit Radial vs Axial — watch coils snap + waveform change |
|    2. Slide N_turns 10→2000 — see mV scale linearly           |
|    3. Slide RPM 10→3000 — see mV scale linearly               |
|    4. C2 angle at 0°  → in-phase  → series gives 2× peak      |
|       C2 angle at 45° → anti-phase (8 magnets) → they cancel! |
|    5. More magnets → higher frequency, same peak amplitude     |
+================================================================+
""")
    viz = Viz()
    viz.run()