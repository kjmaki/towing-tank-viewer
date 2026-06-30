import React, { useRef, useMemo, useEffect, useState, useCallback } from "react";

/**
 * Directional irregular sea in a reflecting towing tank.
 * JONSWAP frequency spectrum x cos^2s directional spread, each (omega, theta) component
 * realized as the Dalrymple (1989) modal tank field (raw paddle or designer), then summed
 * with random phases. Tank: h = 3 m, width 6.6 m (b = 3.3 m), 12 modes/parity.
 *
 *   eta(x,y,t) = sum_i Re{ Phi_i(x,y) e^{-i w_i t} },   Phi_i = sum_j a_ij e^{i eps_ij} Phi_ij
 *
 * Key tank limit: a component at heading theta carries lambda0 = k sin(theta); it only reads
 * as oblique if lambda0 >= pi/2b. Headings within +/- theta_c = asin(pi/2bk) of head-on
 * collapse toward plane. Only the spread's wings are genuinely short-crested.
 */

const G_ACC = 9.81;
const H_DEPTH = 3.0;
const B_HALF = 3.3;
const NM = 12;
const PX_PER_M = 28;
const NY = 120;
const POLE_TOL = 1e-7;

const QUALITY = { low: [14, 7], med: [20, 9], high: [28, 13], vhigh: [60, 15], max: [200, 20] };

function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function solveK(omega) {
  let k = (omega * omega) / G_ACC;
  if (k * H_DEPTH < 1) k = omega / Math.sqrt(G_ACC * H_DEPTH);
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

function jonswap(f, fp, gamma) {
  if (f <= 0) return 0;
  const sig = f <= fp ? 0.07 : 0.09;
  const r = Math.exp(-((f - fp) ** 2) / (2 * sig * sig * fp * fp));
  return Math.pow(f, -5) * Math.exp(-1.25 * Math.pow(fp / f, 4)) * Math.pow(gamma, r);
}

// modal coefficients (complex) for one (k, lambda0); designer aligns propagating modes to xm
function modeCoef(k, lam0, designer, xm) {
  const kx0 = Math.sqrt(Math.max(k * k - lam0 * lam0, 0));
  const sinLb = Math.sin(lam0 * B_HALF), cosLb = Math.cos(lam0 * B_HALF);
  const ev = new Array(NM), od = new Array(NM);
  for (let n = 0; n < NM; n++) {
    const sgn = n % 2 === 0 ? 1 : -1;
    // even
    {
      const a = (n * Math.PI) / B_HALF, prop = a < k;
      const p = prop ? Math.sqrt(k * k - a * a) : 0;
      const q = prop ? 0 : Math.sqrt(a * a - k * k);
      let fA;
      if (Math.abs(lam0 * lam0 - a * a) < POLE_TOL) fA = a < POLE_TOL ? B_HALF : (sgn * B_HALF) / 2;
      else fA = (lam0 * sinLb) / (lam0 * lam0 - a * a);
      const real = (-2 * sgn * fA) / B_HALF, den = p * p + q * q || 1;
      let Ar = (real * p) / den, Ai = (real * -q) / den;
      if (designer) {
        if (!prop) { Ar = 0; Ai = 0; }
        else { const ph = (kx0 - p) * xm, c = Math.cos(ph), s = Math.sin(ph); const r2 = Ar * c - Ai * s; Ai = Ar * s + Ai * c; Ar = r2; }
      }
      ev[n] = { Ar, Ai, p, q };
    }
    // odd
    {
      const a = ((n + 0.5) * Math.PI) / B_HALF, prop = a < k;
      const p = prop ? Math.sqrt(k * k - a * a) : 0;
      const q = prop ? 0 : Math.sqrt(a * a - k * k);
      let fB;
      if (Math.abs(lam0 * lam0 - a * a) < POLE_TOL) fB = (-sgn * B_HALF) / 2;
      else fB = (lam0 * cosLb) / (lam0 * lam0 - a * a);
      const real = (-2 * sgn * fB) / B_HALF, den = p * p + q * q || 1;
      let Br = (real * q) / den, Bi = (real * p) / den; // i/kx folded
      if (designer) {
        if (!prop) { Br = 0; Bi = 0; }
        else { const ph = (kx0 - p) * xm, c = Math.cos(ph), s = Math.sin(ph); const r2 = Br * c - Bi * s; Bi = Br * s + Bi * c; Br = r2; }
      }
      od[n] = { Br, Bi, p, q };
    }
  }
  return { ev, od };
}

function buildSea(p) {
  const { Hs, Tp, gamma, sigmaDeg, theta0Deg, xMax, designer, xm, seed, quality } = p;
  const [Nf, Nd] = QUALITY[quality];
  const rng = mulberry32(seed);
  const fp = 1 / Tp;
  const fLo = 0.55 * fp, fHi = 2.2 * fp, df = (fHi - fLo) / Nf;

  const freqs = [], Sval = [];
  for (let i = 0; i < Nf; i++) { const f = fLo + (i + 0.5) * df; freqs.push(f); Sval.push(jonswap(f, fp, gamma)); }
  const a_f = Sval.map((S) => Math.sqrt(2 * Math.max(S, 0) * df));

  const sigmaR = Math.max((sigmaDeg * Math.PI) / 180, (1.0 * Math.PI) / 180);
  const sExp = Math.max(0.5, 2 / (sigmaR * sigmaR) - 1);
  const spreadHalf = Math.max((0.3 * Math.PI) / 180, Math.min(3 * sigmaR, (78 * Math.PI) / 180));
  const th0 = (theta0Deg * Math.PI) / 180;
  const dirs = [], wj = []; const dth = (2 * spreadHalf) / Nd; let wsum = 0;
  for (let j = 0; j < Nd; j++) {
    const th = th0 - spreadHalf + (j + 0.5) * dth;
    const d = Math.abs(th - th0) < Math.PI / 2 ? Math.pow(Math.cos(th - th0), 2 * sExp) : 0;
    dirs.push(th); wj.push(d * dth); wsum += d * dth;
  }
  for (let j = 0; j < Nd; j++) wj[j] /= wsum || 1;

  // component amplitudes, normalized so total variance -> (Hs/4)^2
  const A = []; let var0 = 0;
  for (let i = 0; i < Nf; i++) { A.push([]); for (let j = 0; j < Nd; j++) { const aij = a_f[i] * Math.sqrt(wj[j]); A[i][j] = aij; var0 += 0.5 * aij * aij; } }
  const scale = Math.sqrt((Hs / 4) ** 2 / Math.max(var0, 1e-12));

  const Ny = NY, Nx = Math.max(60, Math.min(340, Math.round(xMax * 12)));
  // cross-tank basis (geometry only -> shared across all components)
  const cosL = new Float32Array(Ny * NM), sinG = new Float32Array(Ny * NM);
  for (let iy = 0; iy < Ny; iy++) {
    const y = -B_HALF + (iy / (Ny - 1)) * 2 * B_HALF;
    for (let m = 0; m < NM; m++) { cosL[iy * NM + m] = Math.cos((m * Math.PI / B_HALF) * y); sinG[iy * NM + m] = Math.sin(((m + 0.5) * Math.PI / B_HALF) * y); }
  }

  const gRe = [], gIm = [], omega = [];
  for (let i = 0; i < Nf; i++) { gRe.push(new Float32Array(Nx * Ny)); gIm.push(new Float32Array(Nx * Ny)); omega.push(2 * Math.PI * freqs[i]); }

  const comb = Math.PI / (2 * B_HALF);
  let realE = 0, totE = 0;
  const AEer = new Float32Array(NM), AEei = new Float32Array(NM), AEor = new Float32Array(NM), AEoi = new Float32Array(NM);

  for (let i = 0; i < Nf; i++) {
    const w = omega[i], k = solveK(w), gr = gRe[i], gi = gIm[i];
    for (let j = 0; j < Nd; j++) {
      const th = dirs[j], aij = A[i][j] * scale, eps = 2 * Math.PI * rng();
      const lam0 = k * Math.sin(th);
      totE += 0.5 * aij * aij; if (Math.abs(lam0) >= comb) realE += 0.5 * aij * aij;
      const { ev, od } = modeCoef(k, lam0, designer, xm);
      const Wr = aij * Math.cos(eps), Wi = aij * Math.sin(eps);
      const Cr = -Wi, Ci = Wr; // C = i * a e^{i eps}
      for (let cx = 0; cx < Nx; cx++) {
        const x = (cx / (Nx - 1)) * xMax;
        for (let m = 0; m < NM; m++) {
          const e = ev[m]; const dec = e.q ? Math.exp(-e.q * x) : 1; const Er = dec * Math.cos(e.p * x), Ei = dec * Math.sin(e.p * x);
          AEer[m] = e.Ar * Er - e.Ai * Ei; AEei[m] = e.Ar * Ei + e.Ai * Er;
          const o = od[m]; const dc = o.q ? Math.exp(-o.q * x) : 1; const Or = dc * Math.cos(o.p * x), Oi = dc * Math.sin(o.p * x);
          AEor[m] = o.Br * Or - o.Bi * Oi; AEoi[m] = o.Br * Oi + o.Bi * Or;
        }
        const cb = cx * Ny;
        for (let iy = 0; iy < Ny; iy++) {
          let Sr = 0, Si = 0; const yb = iy * NM;
          for (let m = 0; m < NM; m++) { const c = cosL[yb + m], sg = sinG[yb + m]; Sr += AEer[m] * c + AEor[m] * sg; Si += AEei[m] * c + AEoi[m] * sg; }
          const idx = cb + iy; gr[idx] += Cr * Sr - Ci * Si; gi[idx] += Cr * Si + Ci * Sr;
        }
      }
    }
  }

  // calibrate to target Hs over the test region (exclude near-paddle evanescent fringe),
  // mirroring how paddle gain is tuned experimentally to hit a target Hs at the test station
  const cx0 = Math.floor(0.28 * Nx);
  let varSum = 0, cnt = 0;
  for (let cx = cx0; cx < Nx; cx++) for (let iy = 0; iy < Ny; iy++) {
    const idx = cx * Ny + iy; let v = 0; for (let i = 0; i < Nf; i++) v += 0.5 * (gRe[i][idx] ** 2 + gIm[i][idx] ** 2); varSum += v; cnt++;
  }
  const calib = (Hs / 4) / Math.sqrt(Math.max(varSum / cnt, 1e-20));
  for (let i = 0; i < Nf; i++) { const gr = gRe[i], gi = gIm[i]; for (let idx = 0; idx < Nx * Ny; idx++) { gr[idx] *= calib; gi[idx] *= calib; } }

  // spatial rms field
  let maxSig = 1e-9; const sigF = new Float32Array(Nx * Ny);
  for (let idx = 0; idx < Nx * Ny; idx++) { let v = 0; for (let i = 0; i < Nf; i++) v += 0.5 * (gRe[i][idx] ** 2 + gIm[i][idx] ** 2); const s = Math.sqrt(v); sigF[idx] = s; if (s > maxSig) maxSig = s; }

  const kp = solveK(2 * Math.PI * fp);
  const thetaC = (Math.asin(Math.min(comb / kp, 1)) * 180) / Math.PI;
  const Lp = (2 * Math.PI) / kp;

  return { gRe, gIm, omega, Nx, Ny, xMax, maxSig, kp, Lp, thetaC, comb,
    realFrac: realE / Math.max(totE, 1e-12), Nf, Nd, freqs, Sval, fp, fLo, fHi,
    dirs, wj, sExp, spreadHalf, th0 };
}

function diverging(v, banded) {
  let t = Math.max(-1, Math.min(1, v));
  if (banded) t = Math.round(t * 7) / 7;
  if (t >= 0) { const a = t; return [248 - a * 42, 244 - a * 175, 236 - a * 196]; }
  const a = -t; return [248 - a * 211, 244 - a * 138, 236 - a * 9];
}

export default function DirectionalSeaField() {
  const [Hs, setHs] = useState(0.14);
  const [Tp, setTp] = useState(1.2);
  const [gamma, setGamma] = useState(3.3);
  const [sigmaDeg, setSigma] = useState(20);
  const [theta0Deg, setTheta0] = useState(0);
  const [xMax, setXMax] = useState(24);
  const [banded, setBanded] = useState(true);
  const [designer, setDesigner] = useState(false);
  const [xm, setXm] = useState(24);
  const [quality, setQuality] = useState("med");
  const [seed, setSeed] = useState(1234);

  const [playing, setPlaying] = useState(true);
  const [speed, setSpeed] = useState(1);
  const [time, setTime] = useState(0);
  const [sea, setSea] = useState(null);
  const [busy, setBusy] = useState(true);

  const canvasRef = useRef(null), offRef = useRef(null), tRef = useRef(0), rafRef = useRef(0);

  useEffect(() => { if (xm > xMax) setXm(xMax); }, [xMax]); // eslint-disable-line

  // debounced (re)generation — longer delay for heavy quality tiers
  useEffect(() => {
    setBusy(true);
    const delay = quality === "max" ? 900 : quality === "vhigh" ? 350 : 60;
    const id = setTimeout(() => {
      const s = buildSea({ Hs, Tp, gamma, sigmaDeg, theta0Deg, xMax, designer, xm: Math.min(xm, xMax), seed, quality });
      setSea(s); setBusy(false);
    }, delay);
    return () => clearTimeout(id);
  }, [Hs, Tp, gamma, sigmaDeg, theta0Deg, xMax, designer, xm, seed, quality]);

  const displayW = Math.round(xMax * PX_PER_M);
  const displayH = Math.round(2 * B_HALF * PX_PER_M);
  const cMax = 0.8 * Hs;

  const draw = useCallback((t) => {
    const cv = canvasRef.current; if (!cv || !sea) return;
    const { gRe, gIm, omega, Nx, Ny } = sea;
    let off = offRef.current;
    if (!off || off.width !== Nx || off.height !== Ny) { off = document.createElement("canvas"); off.width = Nx; off.height = Ny; offRef.current = off; }
    const octx = off.getContext("2d");
    const img = octx.createImageData(Nx, Ny); const data = img.data;
    const Nf = omega.length; const cs = new Float32Array(Nf), sn = new Float32Array(Nf);
    for (let i = 0; i < Nf; i++) { cs[i] = Math.cos(omega[i] * t); sn[i] = Math.sin(omega[i] * t); }
    const inv = 1 / cMax;
    for (let cx = 0; cx < Nx; cx++) {
      const cb = cx * Ny;
      for (let iy = 0; iy < Ny; iy++) {
        const idx = cb + iy; let eta = 0;
        for (let i = 0; i < Nf; i++) eta += gRe[i][idx] * cs[i] + gIm[i][idx] * sn[i];
        const [r, g, b] = diverging(eta * inv, banded);
        const o = (iy * Nx + cx) * 4; data[o] = r; data[o + 1] = g; data[o + 2] = b; data[o + 3] = 255;
      }
    }
    octx.putImageData(img, 0, 0);
    if (cv.width !== displayW || cv.height !== displayH) { cv.width = displayW; cv.height = displayH; }
    const ctx = cv.getContext("2d");
    ctx.imageSmoothingEnabled = true; ctx.imageSmoothingQuality = "high";
    ctx.drawImage(off, 0, 0, displayW, displayH);
  }, [sea, banded, cMax, displayW, displayH]);

  useEffect(() => {
    let last = performance.now();
    const loop = (now) => {
      const dt = (now - last) / 1000; last = now;
      if (playing && sea) { tRef.current = (tRef.current + speed * dt) % 60; setTime(tRef.current); }
      draw(tRef.current);
      rafRef.current = requestAnimationFrame(loop);
    };
    rafRef.current = requestAnimationFrame(loop);
    return () => cancelAnimationFrame(rafRef.current);
  }, [playing, speed, draw, sea]);

  useEffect(() => { draw(tRef.current); }, [draw]);

  // ---- insets (SVG paths from current params; no regen needed) ----
  const specPath = useMemo(() => {
    if (!sea) return "";
    const { freqs, Sval, Nf } = sea; const mx = Math.max(...Sval, 1e-9);
    const W = 150, Hh = 52;
    return freqs.map((f, i) => {
      const x = ((f - sea.fLo) / (sea.fHi - sea.fLo)) * W;
      const y = Hh - (Sval[i] / mx) * (Hh - 4);
      return `${i === 0 ? "M" : "L"}${x.toFixed(1)},${y.toFixed(1)}`;
    }).join(" ");
  }, [sea]);

  const dirData = useMemo(() => {
    if (!sea) return null;
    const W = 150, Hh = 60, span = 50; // degrees half-range shown
    const th0d = (sea.th0 * 180) / Math.PI;
    const pts = [];
    for (let d = -span; d <= span; d += 2) {
      const rel = (d - th0d) * Math.PI / 180;
      const val = Math.abs(rel) < Math.PI / 2 ? Math.pow(Math.cos(rel), 2 * sea.sExp) : 0;
      const x = ((d + span) / (2 * span)) * W;
      const y = Hh - val * (Hh - 4);
      pts.push(`${d === -span ? "M" : "L"}${x.toFixed(1)},${y.toFixed(1)}`);
    }
    const xc = sea.thetaC, x0 = ((-xc + span) / (2 * span)) * W, x1 = ((xc + span) / (2 * span)) * W;
    const xMean = ((th0d + span) / (2 * span)) * W;
    return { path: pts.join(" "), W, Hh, cutL: x0, cutR: x1, mean: xMean, span };
  }, [sea]);

  const xTickStep = xMax <= 16 ? 2 : xMax <= 32 ? 4 : 10;
  const xTicks = []; for (let xv = 0; xv <= xMax + 1e-6; xv += xTickStep) xTicks.push(+xv.toFixed(2));
  const yTicks = [0, B_HALF, 2 * B_HALF];
  const realPct = sea ? Math.round(sea.realFrac * 100) : 0;
  const lvl = realPct >= 55 ? "ok" : realPct >= 25 ? "warn" : "bad";
  const lvlColor = { ok: "#3ba776", warn: "#c98a1e", bad: "#c0504d" }[lvl];

  return (
    <div style={S.root}>
      <style>{CSS}</style>

      <header style={S.header}>
        <div>
          <div style={S.eyebrow}>DIRECTIONAL IRREGULAR SEA · TANK RESPONSE</div>
          <h1 style={S.h1}>Short-crested sea <span style={{ color: "#6aa9ff" }}>η(x, y, t)</span></h1>
        </div>
        <div style={S.tankspec}>JONSWAP × cos²ˢ<br /><span style={{ color: "#7d8aa0" }}>100 × 6.6 m · h = 3 m</span></div>
      </header>

      <div style={S.readouts}>
        <Stat label="Hs (target)" value={`${(Hs * 100).toFixed(0)} cm`} />
        <Stat label="peak Tp" value={`${Tp.toFixed(2)} s`} />
        <Stat label="spread σ" value={`${sigmaDeg}°`} accent />
        <Stat label="mean dir" value={`${theta0Deg}°`} />
        <Stat label="cutoff θc" value={sea ? `${sea.thetaC.toFixed(1)}°` : "—"} />
        <Stat label="realizable E" value={`${realPct}%`} accent={lvl === "ok"} />
      </div>

      <div style={S.stage}>
        <div style={S.yaxis}>
          {yTicks.map((t) => (<span key={t} style={{ ...S.ytick, top: `${((2 * B_HALF - t) / (2 * B_HALF)) * 100}%` }}>{t.toFixed(1)}</span>))}
          <span style={S.ylabel}>y (m)</span>
        </div>
        <div style={S.canvasWrap}>
          <div style={S.scroller}>
            <canvas ref={canvasRef} style={{ width: displayW, height: displayH, display: "block" }} />
            <div style={S.paddleEdge} title="wavemaker (x = 0)" />
            {designer && (<div style={{ ...S.xmMark, left: `${(Math.min(xm, xMax) / xMax) * 100}%` }} title={`clean target xm = ${xm} m`}><span>xₘ</span></div>)}
            {busy && <div style={S.busy}>generating sea…</div>}
            <div style={S.xticks}>{xTicks.map((t) => (<span key={t} style={{ ...S.xtick, left: `${(t / xMax) * 100}%` }}>{t}</span>))}</div>
          </div>
          <div style={S.xlabel}>x — distance from wavemaker (m) →</div>
        </div>
        <div style={S.cbar}>
          <div style={S.cbarGrad} />
          <span style={{ ...S.cbarT, top: 0 }}>+{(cMax * 100).toFixed(0)}</span>
          <span style={{ ...S.cbarT, top: "50%" }}>0</span>
          <span style={{ ...S.cbarT, top: "100%" }}>−{(cMax * 100).toFixed(0)}</span>
          <span style={S.cbarLabel}>η (cm)</span>
        </div>
      </div>

      <div style={{ ...S.verdict, borderColor: lvlColor }}>
        <span style={{ ...S.verdictDot, background: lvlColor }} />
        <span>
          {sea && (
            <>Headings within <b>±{sea.thetaC.toFixed(1)}°</b> of head-on fall below the transverse cutoff at the peak and render as essentially plane.{" "}
            With σ = {sigmaDeg}° about {theta0Deg}°, <b>{realPct}%</b> of the directional energy lands beyond cutoff and is realized as true short-crestedness
            {realPct < 40 ? " — steer the mean heading off-axis or shorten Tp to recover more spread." : "."}</>
          )}
        </span>
      </div>

      {/* insets */}
      <div style={S.insets}>
        <div style={S.inset}>
          <div style={S.insetTitle}>frequency spectrum S(f)</div>
          <svg width="150" height="52" style={{ overflow: "visible" }}>
            <line x1="0" y1="52" x2="150" y2="52" stroke="#2a3346" />
            {sea && <line x1={((sea.fp - sea.fLo) / (sea.fHi - sea.fLo)) * 150} y1="0" x2={((sea.fp - sea.fLo) / (sea.fHi - sea.fLo)) * 150} y2="52" stroke="#f0c67455" strokeDasharray="2 2" />}
            <path d={specPath} fill="none" stroke="#6aa9ff" strokeWidth="1.6" />
          </svg>
          <div style={S.insetFoot}>fp = {(1 / Tp).toFixed(2)} Hz · γ = {gamma}</div>
        </div>
        <div style={S.inset}>
          <div style={S.insetTitle}>directional spread D(θ)</div>
          {dirData && (
            <svg width="150" height="60" style={{ overflow: "visible" }}>
              <rect x="0" y="0" width={dirData.cutL} height="60" fill="#3ba77645" />
              <rect x={dirData.cutR} y="0" width={150 - dirData.cutR} height="60" fill="#3ba77645" />
              <rect x={dirData.cutL} y="0" width={dirData.cutR - dirData.cutL} height="60" fill="#c0504d55" />
              <line x1={dirData.mean} y1="0" x2={dirData.mean} y2="60" stroke="#f0c674" strokeWidth="1" />
              <line x1="0" y1="60" x2="150" y2="60" stroke="#2a3346" />
              <path d={dirData.path} fill="none" stroke="#6aa9ff" strokeWidth="1.6" />
              <line x1="75" y1="56" x2="75" y2="60" stroke="#4a5568" />
            </svg>
          )}
          <div style={S.insetFoot}><span style={{ color: "#3ba776" }}>▮</span> resolved · <span style={{ color: "#c0504d" }}>▮</span> below cutoff · <span style={{ color: "#f0c674" }}>│</span> mean</div>
          <div style={{...S.insetFoot, marginTop: 3}}>s ≈ 2/σ² − 1 &nbsp;(σ in rad) &nbsp;→&nbsp; s = {sea ? sea.sExp.toFixed(1) : "—"}</div>
        </div>
      </div>

      <div style={S.transport}>
        <button style={S.play} onClick={() => setPlaying((v) => !v)}>{playing ? "❚❚  pause" : "▶  play"}</button>
        <input style={S.scrubInput} type="range" min={0} max={60} step={0.05} value={time}
          onChange={(e) => { const v = parseFloat(e.target.value); tRef.current = v; setTime(v); draw(v); }} />
        <span style={S.phaseTxt}>t = {time.toFixed(1)} s</span>
        <label style={S.speedWrap}>speed<input type="range" min={0.1} max={3} step={0.1} value={speed} onChange={(e) => setSpeed(parseFloat(e.target.value))} style={S.speedInput} /><span style={S.mono}>{speed.toFixed(1)}×</span></label>
        <button style={S.reroll} onClick={() => setSeed((s) => (s * 1103515245 + 12345) & 0x7fffffff)}>↻ new sea</button>
      </div>

      <div style={S.controls} className="wbf-ctrl">
        <Slider label="sig. wave height Hs" unit="cm" min={2} max={40} step={1} value={Math.round(Hs * 100)} onChange={(v) => setHs(v / 100)} note="sets η scale" />
        <Slider label="peak period Tp" unit="s" min={0.7} max={2.2} step={0.05} value={Tp} onChange={setTp} note={sea ? `λp = ${sea.Lp.toFixed(2)} m` : ""} />
        <Slider label="peakedness γ" unit="" min={1} max={7} step={0.1} value={gamma} onChange={setGamma} note="JONSWAP" />
        <Slider label="directional spread σ" unit="°" min={0} max={30} step={1} value={sigmaDeg} onChange={setSigma}
          note={sea ? `s = ${sea.sExp.toFixed(1)}  (exponent in cos²ˢ,  s ≈ 2/σ² − 1)` : "cos²ˢ spread"} />
        <Slider label="mean heading θ₀" unit="°" min={0} max={30} step={1} value={theta0Deg} onChange={setTheta0} note="off tank axis" />
        <Slider label="view length" unit="m" min={8} max={60} step={1} value={xMax} onChange={setXMax} note="x-window" />

        <div style={S.toggleRow}>
          <Toggle label="contour bands" on={banded} onClick={() => setBanded((b) => !b)} />
          <Toggle label="designer (clean at xₘ)" on={designer} onClick={() => setDesigner((d) => !d)} />
          <div style={S.qual}>
            <span style={S.qualLabel}>components</span>
            {Object.keys(QUALITY).map((q) => (
              <button key={q} onClick={() => setQuality(q)}
                title={q === "max" ? "200×20 — slow build (~5 s), best for final renders" : q === "vhigh" ? "60×15 — moderate build (~1 s)" : ""}
                style={{ ...S.qBtn, ...(quality === q ? S.qBtnOn : {}), ...(q === "max" ? { borderColor: "#f0c67444" } : {}) }}>
                {QUALITY[q][0]}×{QUALITY[q][1]}{q === "max" ? " ★" : ""}
              </button>
            ))}
          </div>
          {designer && (<div style={S.xmSlider}><Slider label="test station xₘ" unit="m" min={2} max={xMax} step={1} value={Math.min(xm, xMax)} onChange={setXm} compact /></div>)}
        </div>
      </div>

      <p style={S.foot}>
        Linear superposition of JONSWAP × cos²ˢ components, each realized as the Dalrymple (1989) modal tank field with random phase.
        “Raw” shows the confused short-crested sea the paddles actually produce, with sidewall reflections; “designer” phase-aligns the
        propagating modes for a clean target sea at xₘ. η is the tank response — its local rms departs from the target Hs near reflections
        and where directional energy falls below cutoff. No nonlinear wave–wave interaction.
      </p>
    </div>
  );
}

function Stat({ label, value, accent }) {
  return (<div style={S.stat}><div style={S.statLabel}>{label}</div><div style={{ ...S.statValue, color: accent ? "#6aa9ff" : "#e8edf6" }}>{value}</div></div>);
}
function Slider({ label, unit, min, max, step, value, onChange, note, compact }) {
  return (
    <div style={{ ...S.sliderBox, ...(compact ? { padding: "8px 0 0" } : {}) }}>
      <div style={S.sliderHead}><span style={S.sliderLabel}>{label}</span><span style={S.sliderVal}>{value}<span style={S.unit}>{unit}</span></span></div>
      <input type="range" min={min} max={max} step={step} value={value} onChange={(e) => onChange(parseFloat(e.target.value))} style={S.range} />
      {note && <div style={S.sliderNote}>{note}</div>}
    </div>
  );
}
function Toggle({ label, on, onClick }) {
  return (<button onClick={onClick} style={{ ...S.toggle, ...(on ? S.toggleOn : {}) }}><span style={{ ...S.toggleDot, ...(on ? S.toggleDotOn : {}) }} />{label}</button>);
}

const CSS = `
  input[type=range]{ -webkit-appearance:none; appearance:none; height:3px; background:#2a3346; border-radius:3px; outline:none; }
  input[type=range]::-webkit-slider-thumb{ -webkit-appearance:none; appearance:none; width:15px; height:15px; border-radius:50%; background:#6aa9ff; cursor:pointer; border:2px solid #0d1322; box-shadow:0 0 0 1px #6aa9ff44; }
  input[type=range]::-moz-range-thumb{ width:13px; height:13px; border-radius:50%; background:#6aa9ff; cursor:pointer; border:2px solid #0d1322; }
  @media (max-width:680px){ .wbf-ctrl{ grid-template-columns:1fr 1fr !important; } }
`;

const S = {
  root: { fontFamily: "'Inter',system-ui,sans-serif", background: "#0b0f1a", color: "#e8edf6", padding: "22px 24px 28px", maxWidth: 1180, margin: "0 auto", borderRadius: 14 },
  header: { display: "flex", justifyContent: "space-between", alignItems: "flex-start", marginBottom: 18 },
  eyebrow: { fontFamily: "'JetBrains Mono',ui-monospace,monospace", fontSize: 10.5, letterSpacing: "0.18em", color: "#6e7a93", marginBottom: 6 },
  h1: { margin: 0, fontSize: 26, fontWeight: 600, letterSpacing: "-0.01em" },
  tankspec: { fontFamily: "'JetBrains Mono',ui-monospace,monospace", fontSize: 11.5, textAlign: "right", color: "#aeb8cc", lineHeight: 1.5 },
  readouts: { display: "grid", gridTemplateColumns: "repeat(6,1fr)", gap: 1, background: "#1a2233", border: "1px solid #1a2233", borderRadius: 8, overflow: "hidden", marginBottom: 16 },
  stat: { background: "#0e1422", padding: "10px 12px" },
  statLabel: { fontSize: 9.5, textTransform: "uppercase", letterSpacing: "0.08em", color: "#6e7a93", marginBottom: 3 },
  statValue: { fontFamily: "'JetBrains Mono',ui-monospace,monospace", fontSize: 16, fontWeight: 500 },
  stage: { display: "grid", gridTemplateColumns: "34px 1fr 56px", gap: 8, marginBottom: 12 },
  yaxis: { position: "relative", fontFamily: "'JetBrains Mono',monospace", fontSize: 10, color: "#7d8aa0" },
  ytick: { position: "absolute", right: 4, transform: "translateY(-50%)" },
  ylabel: { position: "absolute", left: -6, top: "50%", transform: "rotate(-90deg) translateX(50%)", transformOrigin: "left", color: "#5b6680", letterSpacing: "0.1em" },
  canvasWrap: { minWidth: 0 },
  scroller: { position: "relative", overflowX: "auto", border: "1px solid #1d2536", borderRadius: 6, background: "#0e1422" },
  paddleEdge: { position: "absolute", left: 0, top: 0, bottom: 0, width: 3, background: "linear-gradient(180deg,#6aa9ff,#3ba776)", boxShadow: "0 0 8px #6aa9ff88" },
  xmMark: { position: "absolute", top: 0, bottom: 16, width: 0, borderLeft: "1.5px dashed #f0c674", color: "#f0c674", fontSize: 10, fontFamily: "'JetBrains Mono',monospace" },
  busy: { position: "absolute", top: 8, right: 10, background: "#11203a", color: "#6aa9ff", fontSize: 11, fontFamily: "'JetBrains Mono',monospace", padding: "3px 8px", borderRadius: 4, border: "1px solid #6aa9ff44" },
  xticks: { position: "relative", height: 16 },
  xtick: { position: "absolute", transform: "translateX(-50%)", top: 2, fontSize: 10, fontFamily: "'JetBrains Mono',monospace", color: "#7d8aa0" },
  xlabel: { textAlign: "center", fontSize: 11, color: "#8b96ac", marginTop: 4, letterSpacing: "0.04em" },
  cbar: { position: "relative", paddingLeft: 8 },
  cbarGrad: { width: 14, height: "calc(100% - 16px)", borderRadius: 3, background: "linear-gradient(180deg,#ce4528,#f8f4ec,#256ae3)", border: "1px solid #1d2536" },
  cbarT: { position: "absolute", left: 26, transform: "translateY(-50%)", fontSize: 10, fontFamily: "'JetBrains Mono',monospace", color: "#8b96ac" },
  cbarLabel: { position: "absolute", bottom: -2, left: -2, fontSize: 9.5, color: "#6e7a93", whiteSpace: "nowrap" },
  verdict: { display: "flex", alignItems: "flex-start", gap: 10, padding: "11px 14px", borderRadius: 8, background: "#0e1422", border: "1px solid", borderLeft: "3px solid", fontSize: 13, lineHeight: 1.5, color: "#cdd6e6", marginBottom: 16 },
  verdictDot: { width: 8, height: 8, borderRadius: "50%", flexShrink: 0, marginTop: 5 },
  insets: { display: "flex", gap: 14, marginBottom: 18, flexWrap: "wrap" },
  inset: { background: "#0e1422", border: "1px solid #1a2233", borderRadius: 8, padding: "10px 14px" },
  insetTitle: { fontSize: 10, textTransform: "uppercase", letterSpacing: "0.08em", color: "#6e7a93", marginBottom: 8 },
  insetFoot: { fontSize: 9.5, color: "#6e7a93", marginTop: 6, fontFamily: "'JetBrains Mono',monospace" },
  transport: { display: "flex", alignItems: "center", gap: 14, marginBottom: 18, flexWrap: "wrap" },
  play: { fontFamily: "'JetBrains Mono',monospace", fontSize: 12, color: "#0b0f1a", background: "#6aa9ff", border: "none", borderRadius: 6, padding: "8px 16px", cursor: "pointer", fontWeight: 600, minWidth: 92 },
  scrubInput: { flex: 1, minWidth: 140 },
  phaseTxt: { fontFamily: "'JetBrains Mono',monospace", fontSize: 11, color: "#8b96ac", minWidth: 78 },
  speedWrap: { display: "flex", alignItems: "center", gap: 8, fontSize: 11, color: "#8b96ac" },
  speedInput: { width: 72 },
  mono: { fontFamily: "'JetBrains Mono',monospace", color: "#aeb8cc", minWidth: 30 },
  reroll: { fontFamily: "'JetBrains Mono',monospace", fontSize: 11.5, color: "#aeb8cc", background: "#0e1422", border: "1px solid #232c40", borderRadius: 6, padding: "8px 12px", cursor: "pointer" },
  controls: { display: "grid", gridTemplateColumns: "repeat(3,1fr)", gap: 20 },
  sliderBox: { padding: "2px 0" },
  sliderHead: { display: "flex", justifyContent: "space-between", alignItems: "baseline", marginBottom: 8 },
  sliderLabel: { fontSize: 12.5, color: "#cdd6e6" },
  sliderVal: { fontFamily: "'JetBrains Mono',monospace", fontSize: 14, color: "#e8edf6" },
  unit: { color: "#6e7a93", fontSize: 11, marginLeft: 2 },
  range: { width: "100%" },
  sliderNote: { fontSize: 10.5, color: "#6e7a93", marginTop: 6, fontFamily: "'JetBrains Mono',monospace" },
  toggleRow: { gridColumn: "1 / -1", display: "flex", gap: 12, flexWrap: "wrap", alignItems: "center", borderTop: "1px solid #1a2233", paddingTop: 16, marginTop: 2 },
  toggle: { display: "flex", alignItems: "center", gap: 8, background: "#0e1422", border: "1px solid #232c40", color: "#aeb8cc", borderRadius: 7, padding: "8px 14px", fontSize: 12.5, cursor: "pointer" },
  toggleOn: { borderColor: "#6aa9ff66", color: "#e8edf6", background: "#11203a" },
  toggleDot: { width: 8, height: 8, borderRadius: "50%", background: "#3a455f" },
  toggleDotOn: { background: "#6aa9ff", boxShadow: "0 0 6px #6aa9ff" },
  qual: { display: "flex", alignItems: "center", gap: 6 },
  qualLabel: { fontSize: 11.5, color: "#8b96ac", marginRight: 2 },
  qBtn: { fontFamily: "'JetBrains Mono',monospace", fontSize: 10.5, color: "#8b96ac", background: "#0e1422", border: "1px solid #232c40", borderRadius: 5, padding: "5px 8px", cursor: "pointer" },
  qBtnOn: { borderColor: "#6aa9ff66", color: "#e8edf6", background: "#11203a" },
  xmSlider: { flex: 1, minWidth: 200 },
  foot: { fontSize: 11, color: "#6e7a93", lineHeight: 1.6, marginTop: 20, borderTop: "1px solid #1a2233", paddingTop: 14 },
};
