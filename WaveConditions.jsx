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

// ── Wave theory validity diagram ──────────────────────────────────────
// x = d/L (relative depth, log), y = H/L (steepness, log)
// Ur = H·L²/d³ = (H/L)/(d/L)³   →   lines Ur=const are H/L = Ur·(d/L)³
// Breaking (Miche 1944): H/L = 0.142·tanh(2π·d/L)
// Linear valid: π·H/L < 0.1  AND  Ur < 1
// Stokes 5th valid: Ur < 30

const ML = 68, MR = 24, MT = 28, MB = 54;
const SVG_W = 560, SVG_H = 440;
const PW = SVG_W - ML - MR;   // 468
const PH = SVG_H - MT - MB;   // 358
const XLO = -3, XHI = 0;      // d/L: 0.001 → 1
const YLO = -3, YHI = -0.8;   // H/L: 0.001 → ~0.16

function px(v) { return ML + (Math.log10(Math.max(v, 1e-15)) - XLO) / (XHI - XLO) * PW; }
function py(v) { return MT + (1 - (Math.log10(Math.max(v, 1e-15)) - YLO) / (YHI - YLO)) * PH; }
function pathOf(arr) {
  return arr.map((p, i) => `${i ? "L" : "M"}${p[0].toFixed(1)},${p[1].toFixed(1)}`).join("");
}

const N = 600;
const DL_ARR = Array.from({ length: N }, (_, i) =>
  Math.pow(10, -3 + (i / (N - 1)) * 3)
);

// Breaking limit
const MICHE_PTS = (() => {
  const pts = [];
  for (const x of DL_ARR) {
    const y = 0.142 * Math.tanh(2 * PI * x);
    if (x >= 1e-3 && x <= 1.01 && y >= 1e-3 && y <= 0.17) pts.push([px(x), py(y)]);
  }
  return pts;
})();

// Ursell curves: H/L = Ur·(d/L)³, clipped at breaking
function ursellCurve(Ur) {
  const pts = [];
  for (const x of DL_ARR) {
    const y = Ur * x * x * x;
    if (y >= 0.142 * Math.tanh(2 * PI * x)) break;
    if (x >= 1e-3 && x <= 1.01 && y >= 1e-3 && y <= 0.17) pts.push([px(x), py(y)]);
  }
  return pts;
}

// Linear/Airy boundary: min(π·H/L = 0.1, Ur = 1)
//   kH/2 < 0.1  →  H/L < 0.1/π  ≈ 0.0318  (flat in deep water)
//   Ur   < 1    →  H/L < (d/L)³             (steep in shallow water)
const LINEAR_PTS = (() => {
  const pts = [];
  for (const x of DL_ARR) {
    const y = Math.min(0.1 / PI, x * x * x);
    if (x >= 1e-3 && x <= 1.01 && y >= 1e-3 && y <= 0.17) pts.push([px(x), py(y)]);
  }
  return pts;
})();

const UR30_PTS = ursellCurve(30);

// ── Region fills (painter's algorithm: back → front) ──────────────────

// 1. Stokes region = everything below Miche
const STOKES_FILL = (() => {
  if (MICHE_PTS.length < 2) return "";
  const first = MICHE_PTS[0], last = MICHE_PTS[MICHE_PTS.length - 1];
  return [
    ...MICHE_PTS.map((p, i) => `${i ? "L" : "M"}${p[0].toFixed(1)},${p[1].toFixed(1)}`),
    `L${last[0].toFixed(1)},${(MT + PH).toFixed(1)}`,
    `L${first[0].toFixed(1)},${(MT + PH).toFixed(1)}`,
    "Z",
  ].join("");
})();

// 2. Linear region = below LINEAR_PTS
const LINEAR_FILL = (() => {
  if (LINEAR_PTS.length < 2) return "";
  const first = LINEAR_PTS[0], last = LINEAR_PTS[LINEAR_PTS.length - 1];
  return [
    ...LINEAR_PTS.map((p, i) => `${i ? "L" : "M"}${p[0].toFixed(1)},${p[1].toFixed(1)}`),
    `L${last[0].toFixed(1)},${(MT + PH).toFixed(1)}`,
    `L${first[0].toFixed(1)},${(MT + PH).toFixed(1)}`,
    "Z",
  ].join("");
})();

// 3. Cnoidal region = between Ur=30 and Miche (left portion where Ur30 < Miche)
const CNOIDAL_FILL = (() => {
  if (UR30_PTS.length < 2) return "";
  const maxPx = UR30_PTS[UR30_PTS.length - 1][0];
  const micheClip = MICHE_PTS.filter(p => p[0] <= maxPx + 2);
  if (micheClip.length < 2) return "";
  const rev = [...UR30_PTS].reverse();
  return [
    ...micheClip.map((p, i) => `${i ? "L" : "M"}${p[0].toFixed(1)},${p[1].toFixed(1)}`),
    ...rev.map(p => `L${p[0].toFixed(1)},${p[1].toFixed(1)}`),
    "Z",
  ].join("");
})();

// 4. Breaking region = above Miche
const BREAKING_FILL = (() => {
  if (MICHE_PTS.length < 2) return "";
  const first = MICHE_PTS[0], last = MICHE_PTS[MICHE_PTS.length - 1];
  return [
    ...MICHE_PTS.map((p, i) => `${i ? "L" : "M"}${p[0].toFixed(1)},${p[1].toFixed(1)}`),
    `L${last[0].toFixed(1)},${MT.toFixed(1)}`,
    `L${first[0].toFixed(1)},${MT.toFixed(1)}`,
    "Z",
  ].join("");
})();

const SUP = { "-": "⁻", 0: "⁰", 1: "¹", 2: "²", 3: "³", 4: "⁴" };
function expLabel(e) {
  if (e === 0) return "1";
  return "10" + String(e).split("").map((c) => SUP[c] ?? c).join("");
}

function WaveTheoryDiagram({ dOverL, HoverL }) {
  const dotX = dOverL != null ? px(dOverL) : null;
  const dotY = HoverL  != null ? py(HoverL)  : null;
  const inPlot =
    dotX != null && dotX >= ML && dotX <= ML + PW &&
    dotY != null && dotY >= MT && dotY <= MT + PH;

  const xTicks = [1e-3, 1e-2, 1e-1, 1e0];
  const yTicks = [1e-3, 1e-2, 1e-1];

  return (
    <svg
      viewBox={`0 0 ${SVG_W} ${SVG_H}`}
      style={{ width: "100%", maxWidth: SVG_W, display: "block", background: "#0a0e1a", borderRadius: 8 }}
    >
      {/* region fills — back to front */}
      <path d={STOKES_FILL}  fill="#4a7abf" fillOpacity="0.13" />
      <path d={LINEAR_FILL}  fill="#c8a830" fillOpacity="0.13" />
      <path d={CNOIDAL_FILL} fill="#3ba776" fillOpacity="0.18" />
      <path d={BREAKING_FILL} fill="#c0504d" fillOpacity="0.20" />

      {/* major grid */}
      {xTicks.map((v) => (
        <line key={`gx${v}`} x1={px(v)} y1={MT} x2={px(v)} y2={MT + PH} stroke="#1a2844" strokeWidth="1" />
      ))}
      {yTicks.map((v) => (
        <line key={`gy${v}`} x1={ML} y1={py(v)} x2={ML + PW} y2={py(v)} stroke="#1a2844" strokeWidth="1" />
      ))}
      {/* minor grid */}
      {[-2.5, -1.5, -0.5].map((e) => (
        <line key={`mx${e}`} x1={px(10 ** e)} y1={MT} x2={px(10 ** e)} y2={MT + PH} stroke="#121a2c" strokeWidth="1" />
      ))}
      {[-2.5, -1.5].map((e) => (
        <line key={`my${e}`} x1={ML} y1={py(10 ** e)} x2={ML + PW} y2={py(10 ** e)} stroke="#121a2c" strokeWidth="1" />
      ))}

      {/* depth reference lines: d/L = 0.05 (shallow) and d/L = 0.5 (deep) */}
      <line x1={px(0.5)}  y1={MT} x2={px(0.5)}  y2={MT + PH} stroke="#3a4f6a" strokeWidth="1" strokeDasharray="5,4" />
      <line x1={px(0.05)} y1={MT} x2={px(0.05)} y2={MT + PH} stroke="#3a4f6a" strokeWidth="1" strokeDasharray="5,4" />

      {/* theory boundary curves */}
      <path d={pathOf(MICHE_PTS)}  fill="none" stroke="#c0504d" strokeWidth="2.5" />
      <path d={pathOf(UR30_PTS)}   fill="none" stroke="#3ba776" strokeWidth="1.8" strokeDasharray="7,3" />
      <path d={pathOf(LINEAR_PTS)} fill="none" stroke="#c8a830" strokeWidth="1.5" strokeDasharray="5,3" />

      {/* region labels */}
      <text x={px(0.5)}  y={py(6e-3)}  fill="#c8a830" fontSize="10.5" textAnchor="middle">Linear (Airy)</text>
      <text x={px(0.28)} y={py(4.5e-2)} fill="#7ab0e0" fontSize="10.5" textAnchor="middle">Stokes 5th order</text>
      <text x={px(0.065)} y={py(2e-2)} fill="#3ba776" fontSize="9"    textAnchor="middle">cnoidal /</text>
      <text x={px(0.065)} y={py(2e-2) + 12} fill="#3ba776" fontSize="9" textAnchor="middle">shallow</text>
      <text x={px(0.22)} y={MT + 12}  fill="#c0504d" fontSize="10"   textAnchor="middle">Breaking</text>

      {/* curve labels */}
      <text x={px(0.38)} y={py(0.115)} fill="#c0504d" fontSize="8.5" textAnchor="middle">
        H/L = 0.142·tanh(2π·d/L)
      </text>
      <text x={px(0.11)} y={py(1.4e-2)} fill="#3ba776" fontSize="8.5" textAnchor="start">Ur = 30</text>
      <text x={px(0.62)} y={py(0.026)}  fill="#c8a830" fontSize="8.5" textAnchor="end">πH/L = 0.1</text>

      {/* depth zone labels */}
      <text x={(ML + px(0.05)) / 2}         y={MT - 8} fill="#4a6080" fontSize="8.5" textAnchor="middle">shallow</text>
      <text x={(px(0.05) + px(0.5)) / 2}    y={MT - 8} fill="#4a6080" fontSize="8.5" textAnchor="middle">intermediate</text>
      <text x={(px(0.5) + ML + PW) / 2}     y={MT - 8} fill="#4a6080" fontSize="8.5" textAnchor="middle">deep</text>

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

      <text x={ML + PW / 2} y={SVG_H - 6} fill="#6e7a93" fontSize="12" textAnchor="middle">d / L</text>
      <text
        x={13} y={MT + PH / 2}
        fill="#6e7a93" fontSize="12" textAnchor="middle"
        transform={`rotate(-90, 13, ${MT + PH / 2})`}
      >
        H / L
      </text>

      {/* current wave condition dot */}
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
    const omega = (2 * PI) / Tv;
    const k = solveK(omega);
    const L = (2 * PI) / k;
    const c = L / Tv;
    const kh = k * H_DEPTH;

    const dOverL = H_DEPTH / L;
    const HoverL = Hv / L;
    const steepnessBreak = 0.142 * Math.tanh(2 * PI * dOverL); // Miche (1944)
    const pct = (HoverL / steepnessBreak) * 100;
    const urNumber = HoverL / (dOverL * dOverL * dOverL);       // Ur = HL²/d³
    const piHoverL = PI * HoverL;                               // kH/2

    let depthClass, depthColor;
    if (kh > PI)           { depthClass = "deep water";    depthColor = "#6aa9ff"; }
    else if (kh < PI / 10) { depthClass = "shallow water"; depthColor = "#f0c674"; }
    else                   { depthClass = "intermediate";  depthColor = "#a0c4ff"; }

    let gaugeColor;
    if (pct < 50)      gaugeColor = "#3ba776";
    else if (pct < 75) gaugeColor = "#f0c674";
    else if (pct < 90) gaugeColor = "#e07830";
    else               gaugeColor = "#c0504d";

    let theory;
    if (pct >= 100)                           theory = "breaking";
    else if (urNumber > 30)                   theory = "cnoidal / shallow water";
    else if (piHoverL < 0.1 && urNumber < 1) theory = "linear (Airy)";
    else                                      theory = "Stokes (nonlinear)";

    return { k, L, c, kh, dOverL, HoverL, steepnessBreak, pct,
             depthClass, depthColor, gaugeColor, urNumber, piHoverL, theory };
  }, [Tv, Hv]);

  return (
    <div style={S.page}>

      {/* ── Regular Waves ─────────────────────────────────────────── */}
      <div style={S.section}>
        <div style={S.eyebrow}>Wave Conditions · d = {H_DEPTH} m</div>
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
              <Stat label="Wavenumber k"       value={calc.k.toFixed(4)}       unit="rad/m" />
              <Stat label="Wavelength L"       value={calc.L.toFixed(2)}       unit="m" />
              <Stat label="Phase speed c"      value={calc.c.toFixed(3)}       unit="m/s" />
              <Stat label="Shallowness kd"     value={calc.kh.toFixed(3)}      note={calc.depthClass} />
              <Stat label="Steepness H/L"      value={calc.HoverL.toFixed(5)} />
              <Stat label="Breaking limit H/L" value={calc.steepnessBreak.toFixed(5)} note="Miche (1944)" />
              <Stat label="Ursell number Ur"   value={calc.urNumber < 9999 ? calc.urNumber.toFixed(1) : ">9999"} />
              <Stat label="Relative depth d/L" value={calc.dOverL.toFixed(4)} />
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

      {/* ── Wave theory validity diagram ───────────────────────────── */}
      {calc && (
        <div style={S.section}>
          <div style={S.eyebrow}>Applicable Wave Theory</div>
          <h1 style={S.h1}>Wave theory validity diagram</h1>
          <div style={S.diagramTitle}>
            Axes: relative depth d/L and steepness H/L (d = still-water depth = {H_DEPTH} m).
            Breaking limit: Miche (1944). Stokes/cnoidal boundary: Ur = HL²/d³ = 30.
            Linear valid where πH/L &lt; 0.1 and Ur &lt; 1.
          </div>
          <WaveTheoryDiagram dOverL={calc.dOverL} HoverL={calc.HoverL} />
          <div style={S.regionTag}>
            current condition → <span style={{ color: calc.gaugeColor }}>{calc.theory}</span>
          </div>
        </div>
      )}

      {/* ── Irregular Waves placeholder ───────────────────────────── */}

    </div>
  );
}
