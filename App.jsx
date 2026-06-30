import React, { useState } from 'react'
import WaveBasinField from './WaveBasinField.jsx'
import DirectionalSeaField from './DirectionalSeaField.jsx'
import WaveConditions from './WaveConditions.jsx'

const btn = {
  background: 'none', border: '1px solid #2a3346', color: '#8b96ac',
  padding: '6px 18px', borderRadius: 6, cursor: 'pointer', fontSize: 13,
}
const activeBtn = { ...btn, background: '#1a2233', color: '#e8edf6', borderColor: '#6aa9ff66' }

export default function App() {
  const [tab, setTab] = useState('wave')
  return (
    <div style={{ background: '#0b0f1a', minHeight: '100vh' }}>
      <div style={{ display: 'flex', gap: 8, padding: '12px 24px', borderBottom: '1px solid #1a2233' }}>
        <button style={tab === 'wave' ? activeBtn : btn} onClick={() => setTab('wave')}>
          Regular Waves
        </button>
        <button style={tab === 'sea' ? activeBtn : btn} onClick={() => setTab('sea')}>
          Irregular Waves
        </button>
        <button style={tab === 'conditions' ? activeBtn : btn} onClick={() => setTab('conditions')}>
          Wave Conditions
        </button>
      </div>
      {tab === 'wave' && <WaveBasinField />}
      {tab === 'sea' && <DirectionalSeaField />}
      {tab === 'conditions' && <WaveConditions />}
    </div>
  )
}
