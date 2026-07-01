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

// ── Irregular sea-state reference table ───────────────────────────────
// Hs and Tp (peak period) representative values per ITTC sea state number.
// N. Pacific Tp is longer reflecting stronger swell contribution.
const SEA_STATES = {
  "N. Atlantic": [
    { ss: 1, Hs: 0.05,  Tp:  4.0, desc: "Calm (rippled)" },
    { ss: 2, Hs: 0.30,  Tp:  6.5, desc: "Smooth" },
    { ss: 3, Hs: 0.90,  Tp:  8.5, desc: "Slight" },
    { ss: 4, Hs: 1.90,  Tp: 10.0, desc: "Moderate" },
    { ss: 5, Hs: 3.25,  Tp: 12.0, desc: "Rough" },
    { ss: 6, Hs: 5.00,  Tp: 13.5, desc: "Very rough" },
    { ss: 7, Hs: 7.50,  Tp: 15.5, desc: "High" },
    { ss: 8, Hs: 11.50, Tp: 18.0, desc: "Very high" },
  ],
  "N. Pacific": [
    { ss: 1, Hs: 0.05,  Tp:  4.5, desc: "Calm (rippled)" },
    { ss: 2, Hs: 0.30,  Tp:  7.5, desc: "Smooth" },
    { ss: 3, Hs: 0.90,  Tp: 10.0, desc: "Slight" },
    { ss: 4, Hs: 1.90,  Tp: 12.0, desc: "Moderate" },
    { ss: 5, Hs: 3.25,  Tp: 14.0, desc: "Rough" },
    { ss: 6, Hs: 5.00,  Tp: 16.0, desc: "Very rough" },
    { ss: 7, Hs: 7.50,  Tp: 17.5, desc: "High" },
    { ss: 8, Hs: 11.50, Tp: 20.0, desc: "Very high" },
  ],
};

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
  toggleBtn: {
    background: "none", border: "1px solid #2a3a5a", borderRadius: 6,
    color: "#8b96ac", padding: "5px 14px", cursor: "pointer",
    fontSize: 12, fontFamily: "inherit",
  },
  toggleBtnActive: {
    background: "#1a2233", border: "1px solid #6aa9ff66",
    borderRadius: 6, color: "#e8edf6", padding: "5px 14px",
    cursor: "pointer", fontSize: 12, fontFamily: "inherit",
  },
  ssBtn: {
    background: "none", border: "1px solid #1e2d44", borderRadius: 5,
    color: "#6e7a93", padding: "4px 8px", cursor: "pointer",
    fontSize: 11, fontFamily: "inherit",
  },
  ssBtnActive: {
    background: "#1a2a3a", border: "1px solid #4a7abf",
    borderRadius: 5, color: "#a0c4ff", padding: "4px 8px",
    cursor: "pointer", fontSize: 11, fontFamily: "inherit", fontWeight: 600,
  },
  scaleRow: {
    display: "grid", gridTemplateColumns: "1fr 1fr", gap: "8px 32px",
    background: "#090d18", borderRadius: 8, padding: "14px 18px",
    border: "1px solid #1a2844", marginTop: 18,
  },
  scaleHead: { fontSize: 10, color: "#4a6080", letterSpacing: 1.5, textTransform: "uppercase", marginBottom: 8, gridColumn: "span 2", borderBottom: "1px solid #1a2844", paddingBottom: 6 },
  scaleLabel: { fontSize: 10, color: "#6e7a93", letterSpacing: 1.2, textTransform: "uppercase", marginBottom: 2 },
  scaleVal: { fontSize: 16, color: "#6aa9ff", fontWeight: 700 },
  scaleValModel: { fontSize: 16, color: "#f0c674", fontWeight: 700 },
  scaleUnit: { fontSize: 11, color: "#4a6080", marginLeft: 3 },
  formulaNote: { fontSize: 10, color: "#3a4f6a", marginTop: 12, lineHeight: 1.6 },
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

function IrregularSection({ ocean, setOcean, ssIdx, setSsIdx, lambda, setLambda, customHs, setCustomHs, customTp, setCustomTp }) {
  const entry = SEA_STATES[ocean][ssIdx];
  const lam   = parseFloat(lambda);

  const hsVal = customHs && parseFloat(customHs) > 0 ? parseFloat(customHs) : entry.Hs;
  const tpVal = customTp && parseFloat(customTp) > 0 ? parseFloat(customTp) : entry.Tp;
  const valid = lam > 0 && isFinite(lam);

  const modelHs = valid ? hsVal / lam          : null;
  const modelTp = valid ? tpVal / Math.sqrt(lam) : null;

  const fmt = (v, dec) => v != null ? v.toFixed(dec) : "—";

  return (
    <div style={S.section}>
      <div style={S.eyebrow}>Irregular Seas · Froude Scaling</div>
      <h1 style={S.h1}>Sea state scale conversion</h1>

      {/* ocean selector */}
      <div style={{ marginBottom: 18 }}>
        <div style={{ ...S.label, marginBottom: 6 }}>Ocean</div>
        <div style={{ display: "flex", gap: 8 }}>
          {Object.keys(SEA_STATES).map(o => (
            <button key={o} style={ocean === o ? S.toggleBtnActive : S.toggleBtn}
              onClick={() => setOcean(o)}>{o}</button>
          ))}
        </div>
      </div>

      {/* sea state selector */}
      <div style={{ marginBottom: 18 }}>
        <div style={{ ...S.label, marginBottom: 6 }}>Sea state</div>
        <div style={{ display: "flex", gap: 5, flexWrap: "wrap" }}>
          {SEA_STATES[ocean].map((e, i) => (
            <button key={e.ss} style={ssIdx === i ? S.ssBtnActive : S.ssBtn}
              onClick={() => setSsIdx(i)}>
              SS{e.ss}
            </button>
          ))}
        </div>
        <div style={{ fontSize: 10.5, color: "#4a6080", marginTop: 5 }}>
          {entry.desc} · Hs = {entry.Hs} m, Tp = {entry.Tp} s (full scale reference)
        </div>
      </div>

      {/* scale factor */}
      <div style={{ ...S.inputRow, marginBottom: 18 }}>
        <div style={S.inputGroup}>
          <label style={S.label}>Scale factor λ</label>
          <input style={{ ...S.input, width: 110 }} type="number" min="1" step="1"
            value={lambda} onChange={e => setLambda(e.target.value)} />
        </div>
      </div>

      {/* optional full-scale override */}
      <div style={{ marginBottom: 4 }}>
        <div style={{ ...S.label, marginBottom: 6 }}>Full-scale override (optional)</div>
        <div style={{ display: "flex", gap: 20, flexWrap: "wrap" }}>
          <div style={S.inputGroup}>
            <label style={S.label}>Hs (m)</label>
            <input style={{ ...S.input, width: 110 }} type="number" min="0" step="0.1"
              placeholder={String(entry.Hs)}
              value={customHs} onChange={e => setCustomHs(e.target.value)} />
          </div>
          <div style={S.inputGroup}>
            <label style={S.label}>Tp (s)</label>
            <input style={{ ...S.input, width: 110 }} type="number" min="0" step="0.1"
              placeholder={String(entry.Tp)}
              value={customTp} onChange={e => setCustomTp(e.target.value)} />
          </div>
        </div>
      </div>

      {/* results */}
      {valid && (
        <div style={S.scaleRow}>
          <div style={S.scaleHead}>Scaled model conditions  (1 : {lam})</div>

          <div>
            <div style={S.scaleLabel}>Full-scale Hs</div>
            <div style={S.scaleVal}>{fmt(hsVal, 2)}<span style={S.scaleUnit}>m</span></div>
          </div>
          <div>
            <div style={S.scaleLabel}>Model Hs</div>
            <div style={S.scaleValModel}>{fmt(modelHs, 4)}<span style={S.scaleUnit}>m</span></div>
          </div>

          <div>
            <div style={S.scaleLabel}>Full-scale Tp</div>
            <div style={S.scaleVal}>{fmt(tpVal, 2)}<span style={S.scaleUnit}>s</span></div>
          </div>
          <div>
            <div style={S.scaleLabel}>Model Tp</div>
            <div style={S.scaleValModel}>{fmt(modelTp, 3)}<span style={S.scaleUnit}>s</span></div>
          </div>
        </div>
      )}

      <div style={S.formulaNote}>
        Froude scaling: H<sub>s,m</sub> = H<sub>s</sub> / λ &nbsp;·&nbsp; T<sub>p,m</sub> = T<sub>p</sub> / √λ
      </div>
    </div>
  );
}

export default function WaveConditions() {
  const [T, setT] = useState("2.0");
  const [H, setH] = useState("0.10");

  // irregular seas state
  const [ocean,     setOcean]     = useState("N. Atlantic");
  const [ssIdx,     setSsIdx]     = useState(5); // SS6 default
  const [lambda,    setLambda]    = useState("50");
  const [customHs,  setCustomHs]  = useState("");
  const [customTp,  setCustomTp]  = useState("");

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
    const piHoverL = PI * HoverL;

    let depthClass, depthColor;
    if (kh > PI)           { depthClass = "deep water";    depthColor = "#6aa9ff"; }
    else if (kh < PI / 10) { depthClass = "shallow water"; depthColor = "#f0c674"; }
    else                   { depthClass = "intermediate";  depthColor = "#a0c4ff"; }

    let gaugeColor;
    if (pct < 50)      gaugeColor = "#3ba776";
    else if (pct < 75) gaugeColor = "#f0c674";
    else if (pct < 90) gaugeColor = "#e07830";
    else               gaugeColor = "#c0504d";

    // Le Méhauté (1976) classification
    let theory;
    if (pct >= 100)                                                  theory = "breaking";
    else if (urNumber < 26)                                          theory = "linear (Airy)";
    else if (dOverL > 1/25 && pct < 25)                             theory = "Stokes 2nd order";
    else if (dOverL > 1/25)                                          theory = "Stokes 3rd / 4th order";
    else                                                             theory = "cnoidal / stream function";

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

      {/* ── Wavemaker performance figure ──────────────────────────── */}
      <div style={S.section}>
        <div style={S.eyebrow}>Wavemaker Performance · d = {H_DEPTH} m</div>
        <h1 style={S.h1}>Wave height vs period</h1>
        <div style={S.diagramTitle}>
          Regular and irregular performance envelopes with breaking limit and depth regime boundaries.
        </div>
        <img
          src={`${import.meta.env.BASE_URL}comparison-wave-regimes.png`}
          alt="Wavemaker performance envelope"
          style={{ width: "100%", borderRadius: 8, display: "block" }}
        />
        {calc && (
          <div style={S.regionTag}>
            current condition → <span style={{ color: calc.gaugeColor }}>{calc.theory}</span>
          </div>
        )}
      </div>

      {/* ── Irregular Seas ────────────────────────────────────────── */}
      <IrregularSection
        ocean={ocean} setOcean={setOcean}
        ssIdx={ssIdx} setSsIdx={setSsIdx}
        lambda={lambda} setLambda={setLambda}
        customHs={customHs} setCustomHs={setCustomHs}
        customTp={customTp} setCustomTp={setCustomTp}
      />

    </div>
  );
}
