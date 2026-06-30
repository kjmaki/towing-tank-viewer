import React, { useState, useMemo } from "react";

const G = 9.81;
const H_DEPTH = 3.0;

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

const S = {
  page: {
    background: "#0b0f1a", minHeight: "100vh", padding: "28px 36px",
    fontFamily: "'JetBrains Mono', 'Fira Mono', monospace",
  },
  section: {
    background: "#0d1220", border: "1px solid #1a2844", borderRadius: 10,
    padding: "22px 28px", marginBottom: 28, maxWidth: 660,
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

    return { k, L, c, kh, steepness, steepnessBreak, pct, depthClass, depthColor, gaugeColor };
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

      {/* ── Irregular Waves placeholder ───────────────────────────── */}

    </div>
  );
}
