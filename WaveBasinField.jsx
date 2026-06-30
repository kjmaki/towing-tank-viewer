import React, { useRef, useMemo, useEffect, useState, useCallback } from "react";

/**
 * Directional wavemaker field — sidewall-reflection model (Dalrymple, JHR 27:1, 1989).
 *
 * Constant-depth, full-span (a = b) flap wavemaker. Forward "raw" model reconstructs the
 * field that actually results from imposing a snake motion of obliquity theta at x = 0,
 * including the cross-tank mode-beating from the two reflective walls (eqs 18-19, 28-29).
 * "Designer" mode phase-aligns the propagating modes to give a clean oblique crest at x_m
 * (eqs 30-31, flat-bottom limit).
 *
 * Pattern note: the flap geometry / stroke enters only through a single scalar (G, S0) common
 * to every mode, so it cancels under normalization and does NOT affect the contour shapes.
 * The field shape depends solely on (theta, T) -> k, lambda0, and the modal sum.
 *
 * Tank (this build): h = 3 m, width 6.6 m (b = 3.3 m), 12 paddles -> 12 transverse modes/parity.
 */

const G_ACC = 9.81;
const H_DEPTH = 3.0;
const B_HALF = 3.3;          // half-width; 12 x 0.55 m = 6.6 m span
const NM = 12;               // modes per parity = paddle count (transverse resolution cap)
const PX_PER_M = 28;         // square pixels -> crest angles render honestly
const GRID_NX_CAP = 760;
const POLE_TOL = 1e-7;

function solveK(omega) {
  // Newton on omega^2 = g k tanh(k h)
  let k = (omega * omega) / G_ACC;        // deep-water seed
  if (k * H_DEPTH < 1) k = omega / Math.sqrt(G_ACC * H_DEPTH); // shallow seed
  for (let i = 0; i < 80; i++) {
    const th = Math.tanh(k * H_DEPTH);
    const f = G_ACC * k * th - omega * omega;
    const df = G_ACC * th + G_ACC * k * H_DEPTH * (1 - th * th);
    const kn = k - f / df;
    if (!isFinite(kn) || kn <= 0) break;
    if (Math.abs(kn - k) < 1e-12) { k = kn; break; }
    k = kn;
  }
  return k;
}

function buildField(thetaDeg, T, xMax, designer, xm) {
  const omega = (2 * Math.PI) / T;
  const k = solveK(omega);
  const L = (2 * Math.PI) / k;
  const lam0 = k * Math.sin((thetaDeg * Math.PI) / 180);
  const kx0 = Math.sqrt(Math.max(k * k - lam0 * lam0, 0)); // target oblique x-wavenumber (always real here)

  const H = Math.max(2, Math.round(2 * B_HALF * PX_PER_M));
  const width = Math.max(2, Math.round(xMax * PX_PER_M));
  const gridNx = Math.min(width, GRID_NX_CAP);

  // ---- per-mode eigenvalues, x-wavenumbers, complex coefficients ----
  // even: lambda_n = n pi / b (cos), odd: gamma_n = (n+1/2) pi / b (sin)
  const ev = []; // {alpha, p, q, Ar, Ai, prop}
  const od = [];
  const sinLb = Math.sin(lam0 * B_HALF);
  const cosLb = Math.cos(lam0 * B_HALF);

  for (let n = 0; n < NM; n++) {
    const sgn = n % 2 === 0 ? 1 : -1;

    // even mode
    {
      const a = (n * Math.PI) / B_HALF;
      const prop = a < k;
      const p = prop ? Math.sqrt(k * k - a * a) : 0;
      const q = prop ? 0 : Math.sqrt(a * a - k * k);
      // fA = lam0 sin(lam0 b) / (lam0^2 - a^2), with removable-pole limits
      let fA;
      if (Math.abs(lam0 * lam0 - a * a) < POLE_TOL) {
        fA = a < POLE_TOL ? B_HALF : sgn * B_HALF / 2;
      } else {
        fA = (lam0 * sinLb) / (lam0 * lam0 - a * a);
      }
      // A = -2 (-1)^n fA / (kx b);  1/kx for complex kx = p+iq
      const real = (-2 * sgn * fA) / B_HALF;
      const den = p * p + q * q || 1;
      const Ar = (real * p) / den;
      const Ai = (real * -q) / den;
      ev.push({ a, p, q, Ar, Ai, prop });
    }
    // odd mode
    {
      const a = ((n + 0.5) * Math.PI) / B_HALF;
      const prop = a < k;
      const p = prop ? Math.sqrt(k * k - a * a) : 0;
      const q = prop ? 0 : Math.sqrt(a * a - k * k);
      let fB;
      if (Math.abs(lam0 * lam0 - a * a) < POLE_TOL) {
        fB = -sgn * B_HALF / 2;
      } else {
        fB = (lam0 * cosLb) / (lam0 * lam0 - a * a);
      }
      // B = -2 i (-1)^n fB / (kx b);  factor i: (-2 sgn fB / b) * i / kx
      const real = (-2 * sgn * fB) / B_HALF;     // this is the coefficient of i/kx
      const den = p * p + q * q || 1;
      // i / kx = i (p - iq)/den = (q + i p)/den
      const Br = (real * q) / den;
      const Bi = (real * p) / den;
      od.push({ a, p, q, Br, Bi, prop });
    }
  }

  // ---- designer: phase-align propagating modes to x_m; drop evanescent ----
  if (designer) {
    for (const m of ev) {
      if (!m.prop) { m.Ar = 0; m.Ai = 0; continue; }
      const ph = (kx0 - m.p) * xm;           // exp(i ph)
      const c = Math.cos(ph), s = Math.sin(ph);
      const Ar = m.Ar * c - m.Ai * s, Ai = m.Ar * s + m.Ai * c;
      m.Ar = Ar; m.Ai = Ai;
    }
    for (const m of od) {
      if (!m.prop) { m.Br = 0; m.Bi = 0; continue; }
      const ph = (kx0 - m.p) * xm;
      const c = Math.cos(ph), s = Math.sin(ph);
      const Br = m.Br * c - m.Bi * s, Bi = m.Br * s + m.Bi * c;
      m.Br = Br; m.Bi = Bi;
    }
  }

  // ---- precompute A*E(x) per column, and cos/sin(alpha y) per row ----
  const AEr_e = new Float32Array(gridNx * NM), AEi_e = new Float32Array(gridNx * NM);
  const AEr_o = new Float32Array(gridNx * NM), AEi_o = new Float32Array(gridNx * NM);
  for (let cx = 0; cx < gridNx; cx++) {
    const x = (cx / (gridNx - 1 || 1)) * xMax;
    for (let m = 0; m < NM; m++) {
      const e = ev[m];
      const dec = e.q ? Math.exp(-e.q * x) : 1;
      const Er = dec * Math.cos(e.p * x), Ei = dec * Math.sin(e.p * x);
      AEr_e[cx * NM + m] = e.Ar * Er - e.Ai * Ei;
      AEi_e[cx * NM + m] = e.Ar * Ei + e.Ai * Er;
      const o = od[m];
      const deco = o.q ? Math.exp(-o.q * x) : 1;
      const Or = deco * Math.cos(o.p * x), Oi = deco * Math.sin(o.p * x);
      AEr_o[cx * NM + m] = o.Br * Or - o.Bi * Oi;
      AEi_o[cx * NM + m] = o.Br * Oi + o.Bi * Or;
    }
  }
  const cyl = new Float32Array(H * NM), syg = new Float32Array(H * NM);
  for (let iy = 0; iy < H; iy++) {
    const y = -B_HALF + (iy / (H - 1 || 1)) * (2 * B_HALF);
    for (let m = 0; m < NM; m++) {
      cyl[iy * NM + m] = Math.cos(ev[m].a * y);
      syg[iy * NM + m] = Math.sin(od[m].a * y);
    }
  }

  // ---- assemble complex Phi = i * sum ----
  const Phre = new Float32Array(gridNx * H);
  const Phim = new Float32Array(gridNx * H);
  let maxAmp = 1e-9;
  for (let cx = 0; cx < gridNx; cx++) {
    for (let iy = 0; iy < H; iy++) {
      let Sr = 0, Si = 0;
      const cb = cx * NM, yb = iy * NM;
      for (let m = 0; m < NM; m++) {
        const ce = cyl[yb + m], so = syg[yb + m];
        Sr += AEr_e[cb + m] * ce + AEr_o[cb + m] * so;
        Si += AEi_e[cb + m] * ce + AEi_o[cb + m] * so;
      }
      // Phi = i * S
      const pr = -Si, pi = Sr;
      const idx = cx * H + iy;
      Phre[idx] = pr; Phim[idx] = pi;
      const amp = Math.hypot(pr, pi);
      if (amp > maxAmp) maxAmp = amp;
    }
  }

  // ---- diagnostics ----
  const comb = Math.PI / (2 * B_HALF);          // transverse wavenumber spacing
  const nPropEven = Math.max(0, Math.floor((k * B_HALF) / Math.PI - 1e-9) + 1);
  const nPropOdd = Math.max(0, Math.floor((k * B_HALF) / Math.PI - 0.5 + 1e-9) + 1);
  const nSpanned = lam0 / comb;                 // how many transverse modes the oblique wave needs
  let verdict, level;
  if (thetaDeg < 0.5) { verdict = "Plane wave (θ ≈ 0): uniform across the width, no reflection structure."; level = "ok"; }
  else if (lam0 < comb) {
    verdict = "Below transverse cutoff — no oblique mode fits the width. The tank radiates an essentially plane wave with an evanescent fringe near the paddle; true obliquity is not realizable at this period.";
    level = "bad";
  } else if (nSpanned < 2) {
    verdict = "Marginal — only ~1 oblique transverse mode in play, so the crest is coarsely resolved and strongly modulated down-tank.";
    level = "warn";
  } else {
    verdict = `Feasible — ${nSpanned.toFixed(1)} transverse modes resolve the ${thetaDeg}° crest.`;
    level = "ok";
  }

  return { Phre, Phim, gridNx, H, width, maxAmp, k, L, lam0, kx0, omega,
           comb, nPropEven, nPropOdd, nSpanned, verdict, level };
}

function diverging(v, banded) {
  // v in [-1,1] -> blue (trough) / pale (still) / warm-red (crest)
  let t = v;
  if (banded) t = Math.round(t * 7) / 7;
  if (t >= 0) {
    const a = Math.min(t, 1);
    return [248 - a * 42, 244 - a * 175, 236 - a * 196]; // -> ~ (206,69,40)
  } else {
    const a = Math.min(-t, 1);
    return [248 - a * 211, 244 - a * 138, 236 - a * 9];  // -> ~ (37,106,227)
  }
}

export default function WaveBasinField() {
  const [theta, setTheta] = useState(15);
  const [T, setT] = useState(1.0);
  const [xMax, setXMax] = useState(24);
  const [banded, setBanded] = useState(true);
  const [designer, setDesigner] = useState(false);
  const [xm, setXm] = useState(24);
  const [playing, setPlaying] = useState(true);
  const [speed, setSpeed] = useState(1);
  const [phase, setPhase] = useState(0);

  const canvasRef = useRef(null);
  const phaseRef = useRef(0);
  const rafRef = useRef(0);

  // keep x_m inside the tank when xMax shrinks
  useEffect(() => { if (xm > xMax) setXm(xMax); }, [xMax]); // eslint-disable-line

  const field = useMemo(
    () => buildField(theta, T, xMax, designer, Math.min(xm, xMax)),
    [theta, T, xMax, designer, xm]
  );

  const draw = useCallback((ph) => {
    const cv = canvasRef.current;
    if (!cv) return;
    const { Phre, Phim, gridNx, H, width, maxAmp } = field;
    if (cv.width !== width || cv.height !== H) { cv.width = width; cv.height = H; }
    const ctx = cv.getContext("2d");
    const img = ctx.createImageData(width, H);
    const data = img.data;
    const c = Math.cos(ph), s = Math.sin(ph);
    const inv = 1 / maxAmp;
    for (let ix = 0; ix < width; ix++) {
      const cx = gridNx === width ? ix : Math.min(gridNx - 1, Math.round((ix / (width - 1 || 1)) * (gridNx - 1)));
      const base = cx * H;
      for (let iy = 0; iy < H; iy++) {
        const idx = base + iy;
        const v = (Phre[idx] * c + Phim[idx] * s) * inv;
        const [r, g, b] = diverging(v, banded);
        const o = (iy * width + ix) * 4;
        data[o] = r; data[o + 1] = g; data[o + 2] = b; data[o + 3] = 255;
      }
    }
    ctx.putImageData(img, 0, 0);
  }, [field, banded]);

  // animation
  useEffect(() => {
    let last = performance.now();
    const loop = (now) => {
      const dt = (now - last) / 1000; last = now;
      if (playing) {
        phaseRef.current = (phaseRef.current + 2 * Math.PI * speed * dt) % (2 * Math.PI);
        setPhase(phaseRef.current);
      }
      draw(phaseRef.current);
      rafRef.current = requestAnimationFrame(loop);
    };
    rafRef.current = requestAnimationFrame(loop);
    return () => cancelAnimationFrame(rafRef.current);
  }, [playing, speed, draw]);

  // redraw immediately on field/colormap change when paused
  useEffect(() => { draw(phaseRef.current); }, [draw]);

  const onScrub = (val) => {
    phaseRef.current = val; setPhase(val); draw(val);
  };

  const wl = field.L;
  const yTicks = [0, B_HALF, 2 * B_HALF];
  const xTickStep = xMax <= 16 ? 2 : xMax <= 32 ? 4 : 10;
  const xTicks = [];
  for (let xv = 0; xv <= xMax + 1e-6; xv += xTickStep) xTicks.push(+xv.toFixed(2));

  const lvlColor = { ok: "#3ba776", warn: "#c98a1e", bad: "#c0504d" }[field.level];

  return (
    <div style={S.root}>
      <style>{CSS}</style>

      <header style={S.header}>
        <div>
          <div style={S.eyebrow}>DIRECTIONAL WAVEMAKER · SIDEWALL-REFLECTION FIELD</div>
          <h1 style={S.h1}>Wave basin <span style={{ color: "#6aa9ff" }}>η(x, y, t)</span></h1>
        </div>
        <div style={S.tankspec}>
          100 × 6.6 m · h = 3 m<br />
          <span style={{ color: "#7d8aa0" }}>12 paddles · full span</span>
        </div>
      </header>

      {/* readouts */}
      <div style={S.readouts}>
        <Stat label="period T" value={`${T.toFixed(2)} s`} />
        <Stat label="wavelength λ" value={`${wl.toFixed(2)} m`} />
        <Stat label="wavenumber k" value={`${field.k.toFixed(2)} /m`} />
        <Stat label="kh" value={(field.k * H_DEPTH).toFixed(1)} />
        <Stat label="angle θ" value={`${theta}°`} accent />
        <Stat label="propagating modes" value={`${field.nPropEven + field.nPropOdd}`} />
      </div>

      {/* HERO field */}
      <div style={S.stage}>
        <div style={S.yaxis}>
          {yTicks.map((t) => (
            <span key={t} style={{ ...S.ytick, top: `${((2 * B_HALF - t) / (2 * B_HALF)) * 100}%` }}>
              {t.toFixed(1)}
            </span>
          ))}
          <span style={S.ylabel}>y (m)</span>
        </div>

        <div style={S.canvasWrap}>
          <div style={S.scroller}>
            <canvas ref={canvasRef} style={{ width: field.width, height: field.H, display: "block", imageRendering: "auto" }} />
            <div style={S.paddleEdge} title="wavemaker (x = 0)" />
            {designer && (
              <div style={{ ...S.xmMark, left: `${(Math.min(xm, xMax) / xMax) * 100}%` }} title={`clean wave target x_m = ${xm} m`}>
                <span>x<sub>m</sub></span>
              </div>
            )}
            <div style={S.xticks}>
              {xTicks.map((t) => (
                <span key={t} style={{ ...S.xtick, left: `${(t / xMax) * 100}%` }}>{t}</span>
              ))}
            </div>
          </div>
          <div style={S.xlabel}>x — distance from wavemaker (m) →</div>
        </div>

        <div style={S.cbar}>
          <div style={S.cbarGrad} />
          <span style={{ ...S.cbarT, top: 0 }}>+1</span>
          <span style={{ ...S.cbarT, top: "50%" }}>0</span>
          <span style={{ ...S.cbarT, top: "100%" }}>−1</span>
          <span style={S.cbarLabel}>η / η<sub>max</sub></span>
        </div>
      </div>

      {/* feasibility */}
      <div style={{ ...S.verdict, borderColor: lvlColor }}>
        <span style={{ ...S.verdictDot, background: lvlColor }} />
        <span>{field.verdict}</span>
      </div>

      {/* transport */}
      <div style={S.transport}>
        <button style={S.play} onClick={() => setPlaying((p) => !p)}>
          {playing ? "❚❚  pause" : "▶  play"}
        </button>
        <input style={S.scrubInput} type="range" min={0} max={2 * Math.PI} step={0.01}
          value={phase} onChange={(e) => onScrub(parseFloat(e.target.value))} />
        <span style={S.phaseTxt}>phase {(phase / (2 * Math.PI)).toFixed(2)}·2π</span>
        <label style={S.speedWrap}>
          speed
          <input type="range" min={0.1} max={3} step={0.1} value={speed}
            onChange={(e) => setSpeed(parseFloat(e.target.value))} style={S.speedInput} />
          <span style={S.mono}>{speed.toFixed(1)}×</span>
        </label>
      </div>

      {/* controls */}
      <div style={S.controls} className="wbf-ctrl">
        <Slider label="direction θ" unit="°" min={0} max={20} step={1} value={theta} onChange={setTheta}
          note="spread limit for your tank" />
        <Slider label="period T" unit="s" min={0.6} max={2.5} step={0.05} value={T} onChange={setT}
          note={`λ = ${wl.toFixed(2)} m`} />
        <Slider label="view length" unit="m" min={8} max={60} step={1} value={xMax} onChange={setXMax}
          note="x-window down-tank" />

        <div style={S.toggleRow}>
          <Toggle label="contour bands" on={banded} onClick={() => setBanded((b) => !b)} />
          <Toggle label="designer (clean at xₘ)" on={designer} onClick={() => setDesigner((d) => !d)} />
          {designer && (
            <div style={S.xmSlider}>
              <Slider label="test station xₘ" unit="m" min={2} max={xMax} step={1} value={Math.min(xm, xMax)} onChange={setXm} compact />
            </div>
          )}
        </div>
      </div>

      <p style={S.foot}>
        Forward reconstruction of Dalrymple (1989), constant depth, full-span flap (a = b). “Raw” shows the field that
        results from a snake paddle at obliquity θ — including the short-crested seas and sidewall reflections. “Designer”
        phase-aligns the propagating modes for a clean oblique crest at xₘ. Linear theory; the flap stroke is a common
        scalar and does not change the normalized pattern.
      </p>
    </div>
  );
}

function Stat({ label, value, accent }) {
  return (
    <div style={S.stat}>
      <div style={S.statLabel}>{label}</div>
      <div style={{ ...S.statValue, color: accent ? "#6aa9ff" : "#e8edf6" }}>{value}</div>
    </div>
  );
}

function Slider({ label, unit, min, max, step, value, onChange, note, compact }) {
  return (
    <div style={{ ...S.sliderBox, ...(compact ? { padding: "8px 0 0" } : {}) }}>
      <div style={S.sliderHead}>
        <span style={S.sliderLabel}>{label}</span>
        <span style={S.sliderVal}>{typeof value === "number" ? value : ""}<span style={S.unit}>{unit}</span></span>
      </div>
      <input type="range" min={min} max={max} step={step} value={value}
        onChange={(e) => onChange(parseFloat(e.target.value))} style={S.range} />
      {note && <div style={S.sliderNote}>{note}</div>}
    </div>
  );
}

function Toggle({ label, on, onClick }) {
  return (
    <button onClick={onClick} style={{ ...S.toggle, ...(on ? S.toggleOn : {}) }}>
      <span style={{ ...S.toggleDot, ...(on ? S.toggleDotOn : {}) }} />
      {label}
    </button>
  );
}

const CSS = `
  input[type=range]{ -webkit-appearance:none; appearance:none; height:3px; background:#2a3346; border-radius:3px; outline:none; }
  input[type=range]::-webkit-slider-thumb{ -webkit-appearance:none; appearance:none; width:15px; height:15px; border-radius:50%;
    background:#6aa9ff; cursor:pointer; border:2px solid #0d1322; box-shadow:0 0 0 1px #6aa9ff44; }
  input[type=range]::-moz-range-thumb{ width:13px; height:13px; border-radius:50%; background:#6aa9ff; cursor:pointer; border:2px solid #0d1322; }
  @media (max-width:680px){ .wbf-ctrl{ grid-template-columns:1fr !important; } }
`;

const S = {
  root: { fontFamily: "'Inter',system-ui,sans-serif", background: "#0b0f1a", color: "#e8edf6",
    padding: "22px 24px 28px", maxWidth: 1180, margin: "0 auto", borderRadius: 14 },
  header: { display: "flex", justifyContent: "space-between", alignItems: "flex-start", marginBottom: 18 },
  eyebrow: { fontFamily: "'JetBrains Mono',ui-monospace,monospace", fontSize: 10.5, letterSpacing: "0.18em",
    color: "#6e7a93", marginBottom: 6 },
  h1: { margin: 0, fontSize: 26, fontWeight: 600, letterSpacing: "-0.01em" },
  tankspec: { fontFamily: "'JetBrains Mono',ui-monospace,monospace", fontSize: 11.5, textAlign: "right",
    color: "#aeb8cc", lineHeight: 1.5 },
  readouts: { display: "grid", gridTemplateColumns: "repeat(6,1fr)", gap: 1, background: "#1a2233",
    border: "1px solid #1a2233", borderRadius: 8, overflow: "hidden", marginBottom: 16 },
  stat: { background: "#0e1422", padding: "10px 12px" },
  statLabel: { fontSize: 9.5, textTransform: "uppercase", letterSpacing: "0.08em", color: "#6e7a93", marginBottom: 3 },
  statValue: { fontFamily: "'JetBrains Mono',ui-monospace,monospace", fontSize: 16, fontWeight: 500 },

  stage: { display: "grid", gridTemplateColumns: "34px 1fr 56px", gap: 8, alignItems: "stretch", marginBottom: 12 },
  yaxis: { position: "relative", fontFamily: "'JetBrains Mono',monospace", fontSize: 10, color: "#7d8aa0" },
  ytick: { position: "absolute", right: 4, transform: "translateY(-50%)" },
  ylabel: { position: "absolute", left: -6, top: "50%", transform: "rotate(-90deg) translateX(50%)",
    transformOrigin: "left", color: "#5b6680", letterSpacing: "0.1em" },
  canvasWrap: { minWidth: 0 },
  scroller: { position: "relative", overflowX: "auto", border: "1px solid #1d2536",
    borderRadius: 6, background: "#0e1422" },
  paddleEdge: { position: "absolute", left: 0, top: 0, bottom: 0, width: 3,
    background: "linear-gradient(180deg,#6aa9ff,#3ba776)", boxShadow: "0 0 8px #6aa9ff88" },
  xmMark: { position: "absolute", top: 0, bottom: 16, width: 0, borderLeft: "1.5px dashed #f0c674",
    color: "#f0c674", fontSize: 10, fontFamily: "'JetBrains Mono',monospace" },
  xticks: { position: "relative", height: 16 },
  xtick: { position: "absolute", transform: "translateX(-50%)", top: 2, fontSize: 10,
    fontFamily: "'JetBrains Mono',monospace", color: "#7d8aa0" },
  xlabel: { textAlign: "center", fontSize: 11, color: "#8b96ac", marginTop: 4, letterSpacing: "0.04em" },

  cbar: { position: "relative", paddingLeft: 8 },
  cbarGrad: { width: 14, height: "calc(100% - 16px)", borderRadius: 3,
    background: "linear-gradient(180deg,#ce4528,#f8f4ec,#256ae3)", border: "1px solid #1d2536" },
  cbarT: { position: "absolute", left: 26, transform: "translateY(-50%)", fontSize: 10,
    fontFamily: "'JetBrains Mono',monospace", color: "#8b96ac" },
  cbarLabel: { position: "absolute", bottom: -2, left: -2, fontSize: 9.5, color: "#6e7a93", whiteSpace: "nowrap" },

  verdict: { display: "flex", alignItems: "center", gap: 10, padding: "10px 14px", borderRadius: 8,
    background: "#0e1422", border: "1px solid", borderLeft: "3px solid", fontSize: 13, lineHeight: 1.45,
    color: "#cdd6e6", marginBottom: 16 },
  verdictDot: { width: 8, height: 8, borderRadius: "50%", flexShrink: 0 },

  transport: { display: "flex", alignItems: "center", gap: 14, marginBottom: 18, flexWrap: "wrap" },
  play: { fontFamily: "'JetBrains Mono',monospace", fontSize: 12, color: "#0b0f1a", background: "#6aa9ff",
    border: "none", borderRadius: 6, padding: "8px 16px", cursor: "pointer", fontWeight: 600, minWidth: 92 },
  scrubInput: { flex: 1, minWidth: 160 },
  phaseTxt: { fontFamily: "'JetBrains Mono',monospace", fontSize: 11, color: "#8b96ac", minWidth: 116 },
  speedWrap: { display: "flex", alignItems: "center", gap: 8, fontSize: 11, color: "#8b96ac" },
  speedInput: { width: 80 },
  mono: { fontFamily: "'JetBrains Mono',monospace", color: "#aeb8cc", minWidth: 30 },

  controls: { display: "grid", gridTemplateColumns: "repeat(3,1fr)", gap: 20 },
  sliderBox: { padding: "2px 0" },
  sliderHead: { display: "flex", justifyContent: "space-between", alignItems: "baseline", marginBottom: 8 },
  sliderLabel: { fontSize: 12.5, color: "#cdd6e6" },
  sliderVal: { fontFamily: "'JetBrains Mono',monospace", fontSize: 14, color: "#e8edf6" },
  unit: { color: "#6e7a93", fontSize: 11, marginLeft: 2 },
  range: { width: "100%" },
  sliderNote: { fontSize: 10.5, color: "#6e7a93", marginTop: 6, fontFamily: "'JetBrains Mono',monospace" },

  toggleRow: { gridColumn: "1 / -1", display: "flex", gap: 12, flexWrap: "wrap", alignItems: "flex-start",
    borderTop: "1px solid #1a2233", paddingTop: 16, marginTop: 2 },
  toggle: { display: "flex", alignItems: "center", gap: 8, background: "#0e1422", border: "1px solid #232c40",
    color: "#aeb8cc", borderRadius: 7, padding: "8px 14px", fontSize: 12.5, cursor: "pointer" },
  toggleOn: { borderColor: "#6aa9ff66", color: "#e8edf6", background: "#11203a" },
  toggleDot: { width: 8, height: 8, borderRadius: "50%", background: "#3a455f" },
  toggleDotOn: { background: "#6aa9ff", boxShadow: "0 0 6px #6aa9ff" },
  xmSlider: { flex: 1, minWidth: 200 },

  foot: { fontSize: 11, color: "#6e7a93", lineHeight: 1.6, marginTop: 20, borderTop: "1px solid #1a2233", paddingTop: 14 },
};
