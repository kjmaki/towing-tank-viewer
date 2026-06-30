import React, { useState, useMemo } from "react";

const G = 9.81;
const H_DEPTH = 3.0;
const PI = Math.PI;

function solveK(omega) {
  let k = (omega * omega) / G;
  if (k * H_DEPTH < 1) k = omega / Math.sqrt(G * H_DEPTH);
  for (let i = 0; i < 80; i++) {
    const th = Math.tanh(k * H_DEPTH);
    const f = G * k * th - omega * omega;
    const df = G * th + G * k * H_DEPTH * (1 - th * th);
    const kn = k - f / df;
    if (!isFinite(kn) || kn <= 0) break;
    if (Math.abs(kn - k) < 1e-12) { k = kn; break; }
    k = kn;
  }
  return k;
}

// ── Le Méhauté (1976) wave theory diagram ────────────────────────────
// Axes: X = d/gT² (log, 1e-4..1e0), Y = H/gT² (log, 1e-4..1e-1)
// Curves parametrized by kh = u via X(u) = u·tanh(u)/(4π²)

const ML = 68, MR = 16, MT = 22, MB = 54;
const SVG_W = 560, SVG_H = 460;
const PW = SVG_W - ML - MR;
const PH = SVG_H - MT - MB;
const XLO = -4, XHI = 0;
const YLO = -4, YHI = -1;

function px(v) { return ML + (Math.log10(Math.max(v, 1e-15)) - XLO) / (XHI - XLO) * PW; }
function py(v) { return MT + (1 - (Math.log10(Math.max(v, 1e-15)) - YLO) / (YHI - YLO)) * PH; }
function pathOf(arr) { return arr.map((p, i) => `${i ? "L" : "M"}${p[0].toFixed(1)},${p[1].toFixed(1)}`).join(""); }

// log-spaced kh sample: dense at small kh (where curves bend sharply), reaching
// far enough that X = kh·tanh(kh)/(4π²) sweeps past the right edge of the chart.
const N_KH = 500;
const KH_ARR = Array.from({ length: N_KH }, (_, i) => {
  const t = i / (N_KH - 1);
  return Math.pow(10, Math.log10(0.01) + t * (Math.log10(60) - Math.log10(0.01)));
});

function micheXY(u, fraction) {
  const th = Math.tanh(u);
  const X = (u * th) / (4 * PI * PI);
  const Y = (fraction * 0.142 * th * 2 * PI * X) / u; // H/L = fraction·0.142·tanh(kh)
  return [X, Y];
}

function micheCurve(fraction) {
  const pts = [];
  for (const u of KH_ARR) {
    const [X, Y] = micheXY(u, fraction);
    if (X > 9e-5 && X < 1.05 && Y > 9e-5 && Y < 0.13) pts.push([px(X), py(Y)]);
  }
  return pts;
}

const CNOIDAL_PTS = (() => {
  // Ursell number Hλ²/d³ = 26 boundary: Y = 26·u²·X/(4π²).
  // Only meaningful below the breaking limit — stop once it would cross above Miche.
  const pts = [];
  for (const u of KH_ARR) {
    const th = Math.tanh(u);
    const X = (u * th) / (4 * PI * PI);
    const Y = (26 * u * u * X) / (4 * PI * PI);
    const [, Ymiche] = micheXY(u, 1.0);
    if (Y >= Ymiche) break;
    if (X > 9e-5 && X < 1.05 && Y > 9e-5 && Y < 0.13) pts.push([px(X), py(Y)]);
  }
  return pts;
})();

const SHALLOW_PTS = (() => {
  // Solitary-wave / shallow breaking limit: H/d = 0.78
  const pts = [];
  for (let i = 0; i <= 200; i++) {
    const X = Math.pow(10, XLO + (i * 3) / 200);
    const Y = 0.78 * X;
    if (Y > 9e-5 && Y < 0.13) pts.push([px(X), py(Y)]);
  }
  return pts;
})();

const MICHE_PTS = micheCurve(1.0);
const STOKES5_PTS = micheCurve(0.8);
const STOKES3_PTS = micheCurve(0.35);
const STOKES2_PTS = micheCurve(0.13);

const SUP = { "-": "⁻", 0: "⁰", 1: "¹", 2: "²", 3: "³", 4: "⁴" };
function expLabel(e) {
  if (e === 0) return "1";
  return "10" + String(e).split("").map((c) => SUP[c] ?? c).join("");
}

function LeMeHauteDiagram({ xDim, yDim }) {
  const dotX = xDim != null ? px(xDim) : null;
  const dotY = yDim != null ? py(yDim) : null;
  const inPlot = dotX != null && dotX >= ML && dotX <= ML + PW && dotY >= MT && dotY <= MT + PH;

  const xTicks = [1e-4, 1e-3, 1e-2, 1e-1, 1e0];
  const yTicks = [1e-4, 1e-3, 1e-2, 1e-1];

  return (
    <svg
      viewBox={`0 0 ${SVG_W} ${SVG_H}`}
      style={{ width: "100%", maxWidth: SVG_W, display: "block", background: "#0a0e1a", borderRadius: 8 }}
    >
      {/* major grid */}
      {xTicks.map((v) => (
        <line key={`gx${v}`} x1={px(v)} y1={MT} x2={px(v)} y2={MT + PH} stroke="#1a2844" strokeWidth="1" />
      ))}
      {yTicks.map((v) => (
        <line key={`gy${v}`} x1={ML} y1={py(v)} x2={ML + PW} y2={py(v)} stroke="#1a2844" strokeWidth="1" />
      ))}
      {/* minor grid */}
      {[-3.5, -2.5, -1.5, -0.5].map((e) => (
        <line key={`mx${e}`} x1={px(10 ** e)} y1={MT} x2={px(10 ** e)} y2={MT + PH} stroke="#121a2c" strokeWidth="1" />
      ))}
      {[-3.5, -2.5, -1.5].map((e) => (
        <line key={`my${e}`} x1={ML} y1={py(10 ** e)} x2={ML + PW} y2={py(10 ** e)} stroke="#121a2c" strokeWidth="1" />
      ))}

      {/* theory boundary curves */}
      <path d={pathOf(MICHE_PTS)} fill="none" stroke="#c0504d" strokeWidth="2.5" />
      <path d={pathOf(SHALLOW_PTS)} fill="none" stroke="#e0a830" strokeWidth="2" />
      <path d={pathOf(STOKES5_PTS)} fill="none" stroke="#7a9ccf" strokeWidth="1.4" strokeDasharray="7,3" />
      <path d={pathOf(STOKES3_PTS)} fill="none" stroke="#5a80b8" strokeWidth="1.4" strokeDasharray="7,3" />
      <path d={pathOf(STOKES2_PTS)} fill="none" stroke="#4a6aa0" strokeWidth="1.4" strokeDasharray="7,3" />
      <path d={pathOf(CNOIDAL_PTS)} fill="none" stroke="#3ba776" strokeWidth="1.4" strokeDasharray="4,4" />

      {/* region labels — centered in the (constant-ratio) gaps between curves */}
      <text x={px(0.3)} y={py(7e-4)} fill="#3a4d70" fontSize="11.5" textAnchor="middle">Linear wave theory (Airy)</text>
      <text x={px(0.03)} y={py(1.57e-3)} fill="#4a6694" fontSize="9.5" textAnchor="middle">Stokes 2nd</text>
      <text x={px(0.075)} y={py(4.77e-3)} fill="#5a7aac" fontSize="9.5" textAnchor="middle">Stokes 3rd</text>
      <text x={px(0.13)} y={py(1.196e-2)} fill="#88aee0" fontSize="9.5" textAnchor="middle">Stokes 5th</text>
      <text x={px(6e-4)} y={py(1.7e-4)} fill="#2f8a64" fontSize="9.5" textAnchor="start">Cnoidal waves</text>
      <text x={px(0.22)} y={py(2.8e-2) - 5} fill="#c0504d" fontSize="9.5" textAnchor="middle">deep water breaking · H/λ = 0.142</text>
      <text x={px(0.1)} y={py(0.078) - 7} fill="#e0a830" fontSize="9.5" textAnchor="end">
        shallow/solitary breaking · H/d = 0.78
      </text>

      {/* axes box */}
      <rect x={ML} y={MT} width={PW} height={PH} fill="none" stroke="#2a3a5a" strokeWidth="1.5" />

      {xTicks.map((v) => (
        <g key={`xt${v}`}>
          <line x1={px(v)} y1={MT + PH} x2={px(v)} y2={MT + PH + 5} stroke="#3a4a6a" />
          <text x={px(v)} y={MT + PH + 18} fill="#6e7a93" fontSize="11" textAnchor="middle">
            {expLabel(Math.round(Math.log10(v)))}
          </text>
        </g>
      ))}
      {yTicks.map((v) => (
        <g key={`yt${v}`}>
          <line x1={ML - 5} y1={py(v)} x2={ML} y2={py(v)} stroke="#3a4a6a" />
          <text x={ML - 8} y={py(v) + 4} fill="#6e7a93" fontSize="11" textAnchor="end">
            {expLabel(Math.round(Math.log10(v)))}
          </text>
        </g>
      ))}

      <text x={ML + PW / 2} y={SVG_H - 6} fill="#6e7a93" fontSize="12" textAnchor="middle">d / gT²</text>
      <text
        x={13} y={MT + PH / 2}
        fill="#6e7a93" fontSize="12" textAnchor="middle"
        transform={`rotate(-90, 13, ${MT + PH / 2})`}
      >
        H / gT²
      </text>

      {/* current wave condition */}
      {inPlot && (
        <g>
          <line x1={dotX - 9} y1={dotY} x2={dotX + 9} y2={dotY} stroke="#f0c674" strokeWidth="1.5" />
          <line x1={dotX} y1={dotY - 9} x2={dotX} y2={dotY + 9} stroke="#f0c674" strokeWidth="1.5" />
          <circle cx={dotX} cy={dotY} r="5" fill="#f0c674" stroke="#0a0e1a" strokeWidth="1.5" />
        </g>
      )}
      {dotX != null && !inPlot && (
        <text x={ML + PW / 2} y={MT + PH / 2} fill="#f0c674" fontSize="11" textAnchor="middle">
          (T, H) outside chart range
        </text>
      )}
    </svg>
  );
}

// ── Styles ────────────────────────────────────────────────────────────
const S = {
  page: {
    background: "#0b0f1a", minHeight: "100vh", padding: "28px 36px",
    fontFamily: "'JetBrains Mono', 'Fira Mono', monospace",
  },
  section: {
    background: "#0d1220", border: "1px solid #1a2844", borderRadius: 10,
    padding: "22px 28px", marginBottom: 28, maxWidth: 700,
  },
  eyebrow: {
    fontSize: 10.5, letterSpacing: 2.5, color: "#4a6080",
    textTransform: "uppercase", marginBottom: 6,
  },
  h1: { fontSize: 21, fontWeight: 700, color: "#e8edf6", marginBottom: 22 },
  inputRow: { display: "flex", gap: 28, alignItems: "flex-end", flexWrap: "wrap" },
  inputGroup: { display: "flex", flexDirection: "column", gap: 5 },
  label: { fontSize: 10.5, color: "#6e7a93", letterSpacing: 1.5, textTransform: "uppercase" },
  input: {
    background: "#0b1525", border: "1px solid #2a3a5a", borderRadius: 6,
    color: "#e8edf6", padding: "8px 12px", fontSize: 14, width: 130,
    fontFamily: "inherit", outline: "none",
  },
  divider: { border: "none", borderTop: "1px solid #1a2844", margin: "20px 0" },
  grid: { display: "grid", gridTemplateColumns: "1fr 1fr", gap: "14px 40px" },
  statLabel: { fontSize: 10, color: "#6e7a93", letterSpacing: 1.5, textTransform: "uppercase", marginBottom: 2 },
  statValue: { fontSize: 15, color: "#6aa9ff", fontWeight: 600 },
  statUnit: { color: "#4a6080", fontSize: 11.5, marginLeft: 4, fontWeight: 400 },
  statNote: { fontSize: 10, color: "#4a6080", marginTop: 1 },
  gauge: { marginTop: 20 },
  gaugeHeader: { display: "flex", justifyContent: "space-between", alignItems: "baseline", marginBottom: 6 },
  gaugeLabel: { fontSize: 11, color: "#6e7a93" },
  gaugePct: { fontSize: 13, fontWeight: 600 },
  gaugeTrack: { height: 8, background: "#1a2233", borderRadius: 4, overflow: "hidden" },
  gaugeFill: { height: "100%", borderRadius: 4, transition: "width 0.15s, background 0.15s" },
  diagramTitle: { fontSize: 12, color: "#8b96ac", marginBottom: 10 },
  regionTag: {
    display: "inline-block", marginTop: 10, padding: "4px 10px", borderRadius: 5,
    fontSize: 11.5, fontWeight: 600, background: "#0b1525", border: "1px solid #2a3a5a",
  },
};

function Stat({ label, value, unit, note }) {
  return (
    <div>
      <div style={S.statLabel}>{label}</div>
      <div style={S.statValue}>
        {value}
        {unit && <span style={S.statUnit}>{unit}</span>}
      </div>
      {note && <div style={S.statNote}>{note}</div>}
    </div>
  );
}

export default function WaveConditions() {
  const [T, setT] = useState("2.0");
  const [H, setH] = useState("0.10");

  const Tv = parseFloat(T);
  const Hv = parseFloat(H);

  const calc = useMemo(() => {
    if (!(Tv > 0) || !(Hv > 0)) return null;
    const omega = (2 * Math.PI) / Tv;
    const k = solveK(omega);
    const L = (2 * Math.PI) / k;
    const c = L / Tv;
    const kh = k * H_DEPTH;
    const steepness = Hv / L;
    const steepnessBreak = 0.142 * Math.tanh(kh); // Miche (1954)
    const pct = (steepness / steepnessBreak) * 100;

    let depthClass, depthColor;
    if (kh > Math.PI)           { depthClass = "deep water";     depthColor = "#6aa9ff"; }
    else if (kh < Math.PI / 10) { depthClass = "shallow water";  depthColor = "#f0c674"; }
    else                        { depthClass = "intermediate";   depthColor = "#a0c4ff"; }

    let gaugeColor;
    if (pct < 50)       gaugeColor = "#3ba776";
    else if (pct < 75)  gaugeColor = "#f0c674";
    else if (pct < 90)  gaugeColor = "#e07830";
    else                gaugeColor = "#c0504d";

    // Le Méhauté diagram coordinates and applicable-theory classification
    const xDim = H_DEPTH / (G * Tv * Tv);
    const yDim = Hv / (G * Tv * Tv);
    const ursell = (4 * PI * PI * yDim) / (kh * kh * xDim);

    let theory;
    if (pct >= 100) theory = "breaking";
    else if (kh < 1.5 && ursell > 26) theory = "cnoidal / shallow";
    else if (pct >= 80) theory = "Stokes 5th order";
    else if (pct >= 35) theory = "Stokes 3rd order";
    else if (pct >= 13) theory = "Stokes 2nd order";
    else theory = "linear (Airy)";

    return { k, L, c, kh, steepness, steepnessBreak, pct, depthClass, depthColor, gaugeColor, xDim, yDim, theory };
  }, [Tv, Hv]);

  return (
    <div style={S.page}>

      {/* ── Regular Waves ─────────────────────────────────────────── */}
      <div style={S.section}>
        <div style={S.eyebrow}>Wave Conditions · h = {H_DEPTH} m</div>
        <h1 style={S.h1}>Regular waves</h1>

        <div style={S.inputRow}>
          <div style={S.inputGroup}>
            <label style={S.label}>Period T (s)</label>
            <input
              style={S.input} type="number" min="0.1" step="0.1"
              value={T} onChange={e => setT(e.target.value)}
            />
          </div>
          <div style={S.inputGroup}>
            <label style={S.label}>Height H (m)</label>
            <input
              style={S.input} type="number" min="0.001" step="0.01"
              value={H} onChange={e => setH(e.target.value)}
            />
          </div>
        </div>

        {calc && (
          <>
            <hr style={S.divider} />
            <div style={S.grid}>
              <Stat label="Wavenumber k"    value={calc.k.toFixed(4)}          unit="rad/m" />
              <Stat label="Wavelength λ"    value={calc.L.toFixed(2)}           unit="m" />
              <Stat label="Phase speed c"   value={calc.c.toFixed(3)}           unit="m/s" />
              <Stat
                label="Shallowness kh"
                value={calc.kh.toFixed(3)}
                note={calc.depthClass}
              />
              <Stat label="Steepness H/λ"          value={calc.steepness.toFixed(5)} />
              <Stat label="Breaking limit H/λ"      value={calc.steepnessBreak.toFixed(5)} note="Miche (1954)" />
            </div>

            <div style={S.gauge}>
              <div style={S.gaugeHeader}>
                <span style={S.gaugeLabel}>% of breaking limit</span>
                <span style={{ ...S.gaugePct, color: calc.gaugeColor }}>
                  {calc.pct.toFixed(1)}%
                  {calc.pct >= 100 && "  · BREAKING"}
                </span>
              </div>
              <div style={S.gaugeTrack}>
                <div style={{
                  ...S.gaugeFill,
                  width: `${Math.min(calc.pct, 100)}%`,
                  background: calc.gaugeColor,
                }} />
              </div>
            </div>
          </>
        )}
      </div>

      {/* ── Le Méhauté wave theory diagram ────────────────────────── */}
      {calc && (
        <div style={S.section}>
          <div style={S.eyebrow}>Applicable Wave Theory</div>
          <h1 style={S.h1}>Le Méhauté (1976) diagram</h1>
          <div style={S.diagramTitle}>
            Region boundaries computed from linear dispersion; deep-water breaking (Miche 1954),
            shallow/solitary breaking, and Stokes/cnoidal limits (Ursell number ≈ 26).
          </div>
          <LeMeHauteDiagram xDim={calc.xDim} yDim={calc.yDim} />
          <div style={S.regionTag}>
            current condition → <span style={{ color: calc.gaugeColor }}>{calc.theory}</span>
          </div>
        </div>
      )}

      {/* ── Irregular Waves placeholder ───────────────────────────── */}

    </div>
  );
}
