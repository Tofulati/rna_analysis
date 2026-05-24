import React, { useState, useMemo, useCallback, useRef } from 'react';
import { ScatterChart, Scatter, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, ReferenceLine } from 'recharts';
import { Search, Database } from 'lucide-react';
import { useVirtualizer } from '@tanstack/react-virtual';
import { useGeneData } from './useGeneData';

const MAX_SCATTER_POINTS = 2000;

function useDebounce(value, delay) {
  const [debounced, setDebounced] = React.useState(value);
  React.useEffect(() => {
    const t = setTimeout(() => setDebounced(value), delay);
    return () => clearTimeout(t);
  }, [value, delay]);
  return debounced;
}

const calculateLinearRegression = (data, xKey, yKey) => {
  const n = data.length;
  if (n < 2) return { slope: 0, intercept: 0, r2: 0 };
  let sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;
  for (const point of data) {
    const x = point[xKey] * 100;
    const y = point[yKey];
    sumX += x; sumY += y; sumXY += x * y; sumX2 += x * x;
  }
  const slope = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
  const intercept = (sumY - slope * sumX) / n;
  const yMean = sumY / n;
  let ssRes = 0, ssTot = 0;
  for (const point of data) {
    const x = point[xKey] * 100;
    const y = point[yKey];
    ssRes += (y - (slope * x + intercept)) ** 2;
    ssTot += (y - yMean) ** 2;
  }
  return { slope, intercept, r2: ssTot === 0 ? 0 : 1 - ssRes / ssTot };
};

function downsample(arr, maxPoints) {
  if (arr.length <= maxPoints) return arr;
  const step = arr.length / maxPoints;
  return Array.from({ length: maxPoints }, (_, i) => arr[Math.floor(i * step)]);
}

// Neutral palette — no saturated colors
const C = {
  bg:           '#ffffff',
  surface:      '#d0d0d0',
  border:       '#888888',
  borderLight:  '#555555',
  text:         '#000000',
  textMuted:    '#888888',
  textDim:      '#555555',
  accent:       '#d0d0d0',   // selected state — just lighter text/border
  dot:          '#5a5a5a',   // scatter dots
  dotHover:     '#aaaaaa',
  regLine:      '#666666',
};

const btnBase = {
  padding: '0.45rem 1rem',
  border: `1px solid ${C.border}`,
  borderRadius: '4px',
  background: 'transparent',
  color: C.textMuted,
  fontWeight: 400,
  cursor: 'pointer',
  fontSize: '0.85rem',
  fontFamily: 'inherit',
  letterSpacing: '0.02em',
};

const btnActive = {
  ...btnBase,
  border: `1px solid ${C.accent}`,
  color: C.text,
  background: '#888888',
};

export default function RNAModificationAnalyzer() {
  const [selectedSample,  setSelectedSample]  = useState('MR01_1');
  const [selectedRegion,  setSelectedRegion]  = useState('utr5');
  const [selectedModType, setSelectedModType] = useState('ai');
  const [searchTerm,      setSearchTerm]      = useState('');
  const [selectedGene,    setSelectedGene]    = useState(null);

  const debouncedSearch = useDebounce(searchTerm, 200);
  const { geneDataMR01_1, geneDataMR01_2, isLoading, loadError } = useGeneData();
  const geneData = selectedSample === 'MR01_1' ? geneDataMR01_1 : geneDataMR01_2;

  const geneById = useMemo(() => {
    const map = new Map();
    for (const g of geneData) map.set(g.id, g);
    return map;
  }, [geneData]);

  const samples  = [{ id: 'MR01_1', label: 'MR01-1' }, { id: 'MR01_2', label: 'MR01-2' }];
  const regions  = [
    { id: 'utr5',   label: "5′ UTR"     },
    { id: 'utr3',   label: "3′ UTR"     },
    { id: 'exon',   label: 'Exon'       },
    { id: 'intron', label: 'Intron'     },
    { id: 'total',  label: 'Total gene' },
  ];
  const modTypes = [
    { id: 'ai',     label: 'A-to-I'     },
    { id: 'm6a',    label: 'm6A'        },
    { id: 'either', label: 'Either mod' },
  ];

  const scatterDataFull = useMemo(() => {
    const rateKey = `${selectedRegion}_${selectedModType}_rate`;
    const cpmKey  = `${selectedRegion}_cpm`;
    return geneData.map(g => ({ name: g.name, modRate: g[rateKey] || 0, cpm: g[cpmKey] || 0, id: g.id }));
  }, [geneData, selectedRegion, selectedModType]);

  const scatterData = useMemo(() => downsample(scatterDataFull, MAX_SCATTER_POINTS), [scatterDataFull]);

  const regression = useMemo(() => calculateLinearRegression(scatterDataFull, 'modRate', 'cpm'), [scatterDataFull]);

  const regressionPoints = useMemo(() => {
    const xs = scatterDataFull.map(d => d.modRate * 100);
    const minX = Math.min(...xs), maxX = Math.max(...xs);
    return [
      { x: minX / 100, y: regression.slope * minX + regression.intercept },
      { x: maxX / 100, y: regression.slope * maxX + regression.intercept },
    ];
  }, [scatterDataFull, regression]);

  const filteredGenes = useMemo(() => {
    if (!debouncedSearch) return geneData;
    const lower = debouncedSearch.toLowerCase();
    return geneData.filter(g => g.name.toLowerCase().includes(lower));
  }, [geneData, debouncedSearch]);

  const handleScatterClick = useCallback((data) => {
    if (data?.id != null) setSelectedGene(geneById.get(data.id) ?? null);
  }, [geneById]);

  const parentRef = useRef(null);
  const virtualizer = useVirtualizer({
    count: filteredGenes.length,
    getScrollElement: () => parentRef.current,
    estimateSize: () => 36,
    overscan: 5,
  });

  const currentSample  = samples.find(s => s.id === selectedSample);
  const currentRegion  = regions.find(r => r.id === selectedRegion);
  const currentModType = modTypes.find(m => m.id === selectedModType);

  if (isLoading) return (
    <div style={{ minHeight: '100vh', display: 'flex', alignItems: 'center', justifyContent: 'center', background: C.bg, color: C.textMuted, fontSize: '0.9rem', fontFamily: 'system-ui, sans-serif', letterSpacing: '0.05em' }}>
      Pulling data from database...
    </div>
  );

  if (loadError) return (
    <div style={{ minHeight: '100vh', display: 'flex', alignItems: 'center', justifyContent: 'center', background: C.bg, color: '#cc5555', fontSize: '0.9rem', fontFamily: 'monospace' }}>
      Error: {loadError}
    </div>
  );

  const divider = { borderTop: `1px solid ${C.border}`, margin: '0' };
  const sectionPad = { padding: '1.75rem 2rem' };

  return (
    <div style={{ minHeight: '100vh', background: C.bg, fontFamily: '"SF Mono", "JetBrains Mono", "Fira Code", ui-monospace, monospace', color: C.text, fontSize: '0.875rem' }}>
      <div style={{ maxWidth: '1400px', margin: '0 auto', border: `1px solid ${C.border}` }}>

        {/* Header */}
        <div style={{ ...sectionPad, borderBottom: `1px solid ${C.border}` }}>
          <div style={{ display: 'flex', alignItems: 'baseline', gap: '1rem' }}>
            <span style={{ fontSize: '0.8rem', color: C.textDim, letterSpacing: '0.08em', textTransform: 'uppercase' }}>RNA Modification Analysis</span>
            <span style={{ color: C.textDim }}>·</span>
            <span style={{ color: C.textMuted, fontSize: '0.8rem' }}>{geneData.length.toLocaleString()} genes</span>
            <span style={{ color: C.textDim }}>·</span>
            <span style={{ color: C.textMuted, fontSize: '0.8rem' }}>{currentSample.label}</span>
          </div>
        </div>

        {/* Controls */}
        <div style={{ ...sectionPad, borderBottom: `1px solid ${C.border}`, display: 'flex', gap: '2.5rem', flexWrap: 'wrap', alignItems: 'flex-start' }}>
          {[
            { label: 'sample',   items: samples,   selected: selectedSample,  onSelect: (id) => { setSelectedSample(id); setSelectedGene(null); } },
            { label: 'region',   items: regions,   selected: selectedRegion,  onSelect: setSelectedRegion  },
            { label: 'mod type', items: modTypes,  selected: selectedModType, onSelect: setSelectedModType },
          ].map(group => (
            <div key={group.label}>
              <div style={{ fontSize: '0.72rem', color: C.textDim, letterSpacing: '0.1em', textTransform: 'uppercase', marginBottom: '0.6rem' }}>{group.label}</div>
              <div style={{ display: 'flex', gap: '0.5rem', flexWrap: 'wrap' }}>
                {group.items.map(item => (
                  <button key={item.id} onClick={() => group.onSelect(item.id)}
                    style={group.selected === item.id ? btnActive : btnBase}>
                    {item.label}
                  </button>
                ))}
              </div>
            </div>
          ))}
        </div>

        {/* Scatter */}
        <div style={sectionPad}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'baseline', marginBottom: '1.25rem' }}>
            <span style={{ color: C.textMuted, fontSize: '0.8rem' }}>
              CPM vs {currentModType.label} · {currentRegion.label}
              {scatterData.length < scatterDataFull.length && (
                <span style={{ color: C.textDim }}> · {scatterData.length.toLocaleString()} of {scatterDataFull.length.toLocaleString()} shown</span>
              )}
            </span>
            <span style={{ color: C.textDim, fontSize: '0.78rem', letterSpacing: '0.04em' }}>
              R² = {regression.r2.toFixed(4)} &nbsp; y = {regression.slope.toFixed(2)}x + {regression.intercept.toFixed(2)}
            </span>
          </div>

          <ResponsiveContainer width="100%" height={440}>
            <ScatterChart margin={{ top: 10, right: 20, bottom: 50, left: 50 }}>
              <CartesianGrid strokeDasharray="2 4" stroke={C.borderLight} />
              <XAxis
                type="number" dataKey="modRate" name="Modification Rate"
                tickFormatter={v => `${(v * 100).toFixed(1)}%`}
                label={{ value: `${currentModType.label} (%)`, position: 'insideBottom', offset: -12, style: { fontSize: '11px', fill: C.textDim } }}
                tick={{ fontSize: 10, fill: C.textDim }}
                stroke={C.border}
              />
              <YAxis
                type="number" dataKey="cpm" name="CPM"
                label={{ value: 'CPM', angle: -90, position: 'insideLeft', offset: 14, style: { fontSize: '11px', fill: C.textDim } }}
                tick={{ fontSize: 10, fill: C.textDim }}
                stroke={C.border}
              />
              <Tooltip
                isAnimationActive={false}
                content={({ active, payload }) => {
                  if (!active || !payload?.length) return null;
                  const d = payload[0].payload;
                  return (
                    <div style={{ background: C.surface, border: `1px solid ${C.border}`, padding: '0.6rem 0.85rem', fontSize: '0.78rem', fontFamily: 'inherit' }}>
                      <div style={{ color: C.text, marginBottom: '0.3rem' }}>{d.name}</div>
                      <div style={{ color: C.textMuted }}>{currentModType.label}: {(d.modRate * 100).toFixed(2)}%</div>
                      <div style={{ color: C.textMuted }}>CPM: {d.cpm.toFixed(2)}</div>
                    </div>
                  );
                }}
              />
              <Scatter data={scatterData} fill={C.dot} fillOpacity={0.6} isAnimationActive={false} onClick={handleScatterClick} style={{ cursor: 'pointer' }} r={2} />
              <ReferenceLine
                segment={regressionPoints}
                stroke={C.regLine}
                strokeWidth={1.5}
                strokeDasharray="5 4"
                ifOverflow="extendDomain"
              />
            </ScatterChart>
          </ResponsiveContainer>
        </div>

        {/* Gene Browser */}
        <div style={{ ...sectionPad, borderTop: `1px solid ${C.border}` }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1.25rem' }}>
            <span style={{ fontSize: '0.72rem', color: C.textDim, letterSpacing: '0.1em', textTransform: 'uppercase' }}>Gene browser</span>
            <div style={{ position: 'relative' }}>
              <Search size={13} style={{ position: 'absolute', left: '10px', top: '50%', transform: 'translateY(-50%)', color: C.textDim }} />
              <input
                type="text" placeholder="Search…" value={searchTerm}
                onChange={e => setSearchTerm(e.target.value)}
                style={{
                  paddingLeft: '30px', paddingRight: '10px', height: '30px', width: '220px',
                  background: 'transparent', border: `1px solid ${C.border}`, borderRadius: '3px',
                  color: C.text, fontSize: '0.8rem', fontFamily: 'inherit', outline: 'none',
                }}
              />
            </div>
          </div>

          {/* Virtualised list */}
          <div ref={parentRef} style={{ height: '220px', overflowY: 'auto', border: `1px solid ${C.border}` }}>
            <div style={{ height: virtualizer.getTotalSize(), position: 'relative' }}>
              {virtualizer.getVirtualItems().map(row => {
                const gene = filteredGenes[row.index];
                const isSelected = selectedGene?.id === gene.id;
                return (
                  <div key={gene.id} style={{ position: 'absolute', top: 0, left: 0, width: '100%', height: `${row.size}px`, transform: `translateY(${row.start}px)` }}>
                    <button onClick={() => setSelectedGene(gene)} style={{
                      width: '100%', height: '100%',
                      background: isSelected ? '#1e1e1e' : 'transparent',
                      border: 'none',
                      borderBottom: `1px solid ${C.borderLight}`,
                      color: isSelected ? C.text : C.textMuted,
                      textAlign: 'left', padding: '0 1rem',
                      fontSize: '0.82rem', fontFamily: 'inherit', cursor: 'pointer',
                    }}>
                      {gene.name}
                    </button>
                  </div>
                );
              })}
            </div>
          </div>
        </div>

        {/* Raw data table */}
        {selectedGene?.raw_data && (
          <div style={{ borderTop: `1px solid ${C.border}`, ...sectionPad }}>
            <div style={{ fontSize: '0.72rem', color: C.textDim, letterSpacing: '0.1em', textTransform: 'uppercase', marginBottom: '1.25rem' }}>
              {selectedGene.name} · raw data
            </div>

            <div style={{ overflowX: 'auto' }}>
              <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: '0.8rem' }}>
                <thead>
                  <tr style={{ borderBottom: `1px solid ${C.border}` }}>
                    {['Feature', 'Modification', 'Count', 'CPK', 'MR'].map(h => (
                      <th key={h} style={{ padding: '0.5rem 0.75rem', textAlign: h === 'Feature' || h === 'Modification' ? 'left' : 'right', color: C.textDim, fontWeight: 400, letterSpacing: '0.06em', fontSize: '0.72rem', textTransform: 'uppercase' }}>{h}</th>
                    ))}
                  </tr>
                </thead>
                <tbody>
                  {selectedGene.raw_data.map((row, idx) => (
                    <tr key={idx} style={{ borderBottom: `1px solid ${C.borderLight}` }}>
                      <td style={{ padding: '0.45rem 0.75rem', color: C.text }}>{row.feature}</td>
                      <td style={{ padding: '0.45rem 0.75rem', color: C.textMuted }}>{row.modification}</td>
                      <td style={{ padding: '0.45rem 0.75rem', textAlign: 'right', color: C.textMuted }}>{row.count}</td>
                      <td style={{ padding: '0.45rem 0.75rem', textAlign: 'right', color: C.textMuted }}>{row.cpk.toFixed(2)}</td>
                      <td style={{ padding: '0.45rem 0.75rem', textAlign: 'right', color: C.text, fontVariantNumeric: 'tabular-nums' }}>{row.mr.toFixed(6)}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>

            {/* Summary */}
            <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(160px, 1fr))', gap: '1px', marginTop: '1.5rem', border: `1px solid ${C.border}` }}>
              {[
                { label: 'A-to-I rate',   value: `${((selectedGene.total_ai_rate    || 0) * 100).toFixed(2)}%` },
                { label: 'm6A rate',      value: `${((selectedGene.total_m6a_rate   || 0) * 100).toFixed(2)}%` },
                { label: 'Either rate',   value: `${((selectedGene.total_either_rate|| 0) * 100).toFixed(2)}%` },
                { label: 'Total CPM',     value: (selectedGene.total_cpm || 0).toFixed(2) },
              ].map(({ label, value }) => (
                <div key={label} style={{ padding: '1rem 1.25rem', background: C.surface }}>
                  <div style={{ fontSize: '0.7rem', color: C.textDim, letterSpacing: '0.08em', textTransform: 'uppercase', marginBottom: '0.4rem' }}>{label}</div>
                  <div style={{ fontSize: '1.25rem', color: C.text, fontVariantNumeric: 'tabular-nums' }}>{value}</div>
                </div>
              ))}
            </div>
          </div>
        )}

        {!selectedGene && (
          <div style={{ ...sectionPad, borderTop: `1px solid ${C.border}`, color: C.textDim, fontSize: '0.78rem' }}>
            Select a gene to view raw data.
          </div>
        )}

        {/* Footer */}
        <div style={{ borderTop: `1px solid ${C.border}`, padding: '1rem 2rem', display: 'flex', alignItems: 'center', gap: '0.6rem' }}>
          <Database size={12} style={{ color: C.textDim }} />
          <a
            href="https://github.com/Tofulati/rna_analysis/"
            target="_blank"
            rel="noopener noreferrer"
            style={{
              color: C.textDim,
              fontSize: '0.75rem',
              letterSpacing: '0.05em',
              textDecoration: 'none',
              fontFamily: 'inherit',
            }}
            onMouseEnter={e => e.target.style.color = C.text}
            onMouseLeave={e => e.target.style.color = C.textDim}
          >
            github.com/Tofulati/rna_analysis
          </a>
        </div>

      </div>
    </div>
  );
}