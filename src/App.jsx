import React, { useState, useMemo, useEffect } from 'react';
import { ScatterChart, Scatter, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer, Line, LineChart } from 'recharts';
import { TrendingUp, Search, Database } from 'lucide-react';

// Linear regression calculation
const calculateLinearRegression = (data, xKey, yKey) => {
  const n = data.length;
  let sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;
  
  data.forEach(point => {
    const x = point[xKey] * 100;
    const y = point[yKey];
    sumX += x;
    sumY += y;
    sumXY += x * y;
    sumX2 += x * x;
  });
  
  const slope = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
  const intercept = (sumY - slope * sumX) / n;
  
  const yMean = sumY / n;
  let ssRes = 0, ssTot = 0;
  
  data.forEach(point => {
    const x = point[xKey] * 100;
    const y = point[yKey];
    const yPred = slope * x + intercept;
    ssRes += (y - yPred) ** 2;
    ssTot += (y - yMean) ** 2;
  });
  
  const r2 = 1 - (ssRes / ssTot);
  
  return { slope, intercept, r2 };
};

const RNAModificationAnalyzer = () => {
  const [geneDataMR01_1, setGeneDataMR01_1] = useState([]);
  const [geneDataMR01_2, setGeneDataMR01_2] = useState([]);
  const [selectedSample, setSelectedSample] = useState('MR01_1');
  const [selectedRegion, setSelectedRegion] = useState('utr5');
  const [selectedModType, setSelectedModType] = useState('ai');
  const [searchTerm, setSearchTerm] = useState('');
  const [selectedGene, setSelectedGene] = useState(null);
  const [isLoading, setIsLoading] = useState(true);
  const [loadError, setLoadError] = useState(null);
  
  useEffect(() => {
    const loadData = async () => {
      try {
        setIsLoading(true);
        
        const response1 = await fetch('output/gene_data_MR01_1.json');
        if (!response1.ok) throw new Error('Failed to load MR01_1 data');
        const data1 = await response1.json();
        setGeneDataMR01_1(data1);
        
        const response2 = await fetch('output/gene_data_MR01_2.json');
        if (!response2.ok) throw new Error('Failed to load MR01_2 data');
        const data2 = await response2.json();
        setGeneDataMR01_2(data2);
        
        setIsLoading(false);
      } catch (error) {
        console.error('Error loading data:', error);
        setLoadError(error.message);
        setIsLoading(false);
      }
    };
    
    loadData();
  }, []);
  
  const geneData = selectedSample === 'MR01_1' ? geneDataMR01_1 : geneDataMR01_2;
  
  const samples = [
    { id: 'MR01_1', label: 'MR01-1', color: '#9B59B6' },
    { id: 'MR01_2', label: 'MR01-2', color: '#E74C3C' }
  ];
  
  const regions = [
    { id: 'utr5', label: "5' UTR", color: '#FF6B9D', feature: 'UTR_5' },
    { id: 'utr3', label: "3' UTR", color: '#4ECDC4', feature: 'UTR_3' },
    { id: 'exon', label: 'Exonic', color: '#FFE66D', feature: 'Exon' },
    { id: 'intron', label: 'Intronic', color: '#95E1D3', feature: 'Intron' },
    { id: 'total', label: 'Total Gene', color: '#A8E6CF', feature: null }
  ];
  
  const modTypes = [
    { id: 'ai', label: '% A-to-I', color: '#FF6B9D' },
    { id: 'm6a', label: '% m6A', color: '#4ECDC4' },
    { id: 'either', label: '% Either Modification', color: '#C7CEEA' }
  ];
  
  const scatterData = useMemo(() => {
    const rateKey = `${selectedRegion}_${selectedModType}_rate`;
    const cpmKey = `${selectedRegion}_cpm`;
    
    return geneData.map(gene => ({
      name: gene.name,
      modRate: gene[rateKey] || 0,
      cpm: gene[cpmKey] || 0,
      id: gene.id
    }));
  }, [geneData, selectedRegion, selectedModType]);
  
  const regression = useMemo(() => {
    return calculateLinearRegression(scatterData, 'modRate', 'cpm');
  }, [scatterData]);
  
  const regressionLine = useMemo(() => {
    const minX = Math.min(...scatterData.map(d => d.modRate * 100));
    const maxX = Math.max(...scatterData.map(d => d.modRate * 100));
    
    return [
      { x: minX, y: regression.slope * minX + regression.intercept },
      { x: maxX, y: regression.slope * maxX + regression.intercept }
    ];
  }, [scatterData, regression]);
  
  const filteredGenes = useMemo(() => {
    if (!searchTerm) return geneData;
    return geneData.filter(gene => 
      gene.name.toLowerCase().includes(searchTerm.toLowerCase()) ||
      (gene.chromosome && gene.chromosome.toLowerCase().includes(searchTerm.toLowerCase()))
    );
  }, [geneData, searchTerm]);
  
  const currentSample = samples.find(s => s.id === selectedSample);
  const currentRegion = regions.find(r => r.id === selectedRegion);
  const currentModType = modTypes.find(m => m.id === selectedModType);
  
  if (isLoading) {
    return (
      <div style={{ 
        minHeight: '100vh', 
        display: 'flex', 
        alignItems: 'center', 
        justifyContent: 'center',
        background: 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)',
        color: 'white',
        fontSize: '1.3rem',
        fontWeight: 600,
        fontFamily: '"Space Grotesk", system-ui, sans-serif'
      }}>
        Loading RNA modification data...
      </div>
    );
  }
  
  if (loadError) {
    return (
      <div style={{ 
        minHeight: '100vh', 
        display: 'flex', 
        alignItems: 'center', 
        justifyContent: 'center',
        flexDirection: 'column',
        background: '#1a1a2e',
        color: '#eee',
        padding: '2rem',
        fontFamily: '"JetBrains Mono", monospace'
      }}>
        <h1 style={{ fontSize: '1.8rem', marginBottom: '1rem', color: '#ff6b6b' }}>Error Loading Data</h1>
        <p style={{ fontSize: '1rem', marginBottom: '1rem', color: '#ee5a6f' }}>{loadError}</p>
        <p style={{ fontSize: '0.9rem', color: '#c0c0c0' }}>
          Ensure JSON files are in /output directory
        </p>
      </div>
    );
  }
  
  return (
    <div style={{ 
      minHeight: '100vh', 
      background: 'linear-gradient(to bottom, #0f0c29, #302b63, #24243e)', 
      padding: '2.5rem',
      fontFamily: '"Space Grotesk", -apple-system, sans-serif',
      color: '#e0e0e0'
    }}>
      <div style={{ 
        maxWidth: '1500px', 
        margin: '0 auto',
        background: 'rgba(255, 255, 255, 0.03)',
        borderRadius: '16px',
        boxShadow: '0 8px 32px rgba(0, 0, 0, 0.3)',
        overflow: 'hidden',
        border: '1px solid rgba(255, 255, 255, 0.1)',
        backdropFilter: 'blur(10px)'
      }}>
        {/* Header */}
        <div style={{ 
          background: 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)', 
          padding: '2.5rem',
          borderBottom: '2px solid rgba(255, 255, 255, 0.1)'
        }}>
          <h1 style={{ 
            margin: 0, 
            fontSize: '2.2rem', 
            fontWeight: 800,
            color: 'white',
            display: 'flex',
            alignItems: 'center',
            gap: '1rem',
            letterSpacing: '-0.02em'
          }}>
            <Database size={36} strokeWidth={2.5} />
            RNA Modification Analysis
          </h1>
          <p style={{ margin: '0.75rem 0 0 0', fontSize: '1.05rem', color: 'rgba(255, 255, 255, 0.9)' }}>
            Regression analysis of {geneData.length} protein-coding genes • {currentSample.label}
          </p>
        </div>
        
        {/* Controls */}
        <div style={{ padding: '2rem 2.5rem', borderBottom: '1px solid rgba(255, 255, 255, 0.1)', background: 'rgba(0, 0, 0, 0.2)' }}>
          <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(300px, 1fr))', gap: '2.5rem' }}>
            <div>
              <label style={{ 
                display: 'block', 
                fontSize: '0.75rem', 
                fontWeight: 700, 
                color: '#a0a0a0',
                marginBottom: '0.75rem',
                textTransform: 'uppercase',
                letterSpacing: '0.1em'
              }}>
                Sample
              </label>
              <div style={{ display: 'flex', gap: '0.75rem', flexWrap: 'wrap' }}>
                {samples.map(sample => (
                  <button
                    key={sample.id}
                    onClick={() => {
                      setSelectedSample(sample.id);
                      setSelectedGene(null);
                    }}
                    style={{
                      padding: '0.65rem 1.25rem',
                      border: selectedSample === sample.id ? 'none' : '1px solid rgba(255, 255, 255, 0.2)',
                      borderRadius: '8px',
                      background: selectedSample === sample.id ? sample.color : 'rgba(255, 255, 255, 0.05)',
                      color: selectedSample === sample.id ? 'white' : '#d0d0d0',
                      fontWeight: 600,
                      cursor: 'pointer',
                      transition: 'all 0.2s',
                      fontSize: '0.95rem',
                      fontFamily: '"Space Grotesk", sans-serif'
                    }}
                  >
                    {sample.label}
                  </button>
                ))}
              </div>
            </div>
            
            <div>
              <label style={{ 
                display: 'block', 
                fontSize: '0.75rem', 
                fontWeight: 700, 
                color: '#a0a0a0',
                marginBottom: '0.75rem',
                textTransform: 'uppercase',
                letterSpacing: '0.1em'
              }}>
                Genomic Region
              </label>
              <div style={{ display: 'flex', gap: '0.75rem', flexWrap: 'wrap' }}>
                {regions.map(region => (
                  <button
                    key={region.id}
                    onClick={() => setSelectedRegion(region.id)}
                    style={{
                      padding: '0.65rem 1.25rem',
                      border: selectedRegion === region.id ? `2px solid ${region.color}` : '1px solid rgba(255, 255, 255, 0.2)',
                      borderRadius: '8px',
                      background: selectedRegion === region.id ? region.color : 'rgba(255, 255, 255, 0.05)',
                      color: selectedRegion === region.id ? '#1a1a2e' : '#d0d0d0',
                      fontWeight: selectedRegion === region.id ? 700 : 600,
                      cursor: 'pointer',
                      transition: 'all 0.2s',
                      fontSize: '0.95rem',
                      fontFamily: '"Space Grotesk", sans-serif'
                    }}
                  >
                    {region.label}
                  </button>
                ))}
              </div>
            </div>
            
            <div>
              <label style={{ 
                display: 'block', 
                fontSize: '0.75rem', 
                fontWeight: 700, 
                color: '#a0a0a0',
                marginBottom: '0.75rem',
                textTransform: 'uppercase',
                letterSpacing: '0.1em'
              }}>
                Modification Type
              </label>
              <div style={{ display: 'flex', gap: '0.75rem', flexWrap: 'wrap' }}>
                {modTypes.map(modType => (
                  <button
                    key={modType.id}
                    onClick={() => setSelectedModType(modType.id)}
                    style={{
                      padding: '0.65rem 1.25rem',
                      border: selectedModType === modType.id ? `2px solid ${modType.color}` : '1px solid rgba(255, 255, 255, 0.2)',
                      borderRadius: '8px',
                      background: selectedModType === modType.id ? modType.color : 'rgba(255, 255, 255, 0.05)',
                      color: selectedModType === modType.id ? '#1a1a2e' : '#d0d0d0',
                      fontWeight: selectedModType === modType.id ? 700 : 600,
                      cursor: 'pointer',
                      transition: 'all 0.2s',
                      fontSize: '0.95rem',
                      fontFamily: '"Space Grotesk", sans-serif'
                    }}
                  >
                    {modType.label}
                  </button>
                ))}
              </div>
            </div>
          </div>
        </div>
        
        {/* Scatter Plot */}
        <div style={{ padding: '2.5rem' }}>
          <div style={{ 
            background: 'rgba(255, 255, 255, 0.03)',
            borderRadius: '12px',
            padding: '1.5rem',
            border: '1px solid rgba(255, 255, 255, 0.1)'
          }}>
            <div style={{ 
              display: 'flex', 
              justifyContent: 'space-between', 
              alignItems: 'center',
              marginBottom: '2rem',
              paddingBottom: '1.25rem',
              borderBottom: '1px solid rgba(255, 255, 255, 0.1)'
            }}>
              <h2 style={{ margin: 0, fontSize: '1.5rem', fontWeight: 700, color: '#fff' }}>
                CPM vs {currentModType.label} • {currentRegion.label}
              </h2>
              <div style={{ 
                background: 'rgba(102, 126, 234, 0.15)',
                padding: '0.75rem 1.25rem',
                borderRadius: '8px',
                fontSize: '0.9rem',
                border: '1px solid rgba(102, 126, 234, 0.3)',
                fontFamily: '"JetBrains Mono", monospace'
              }}>
                <span style={{ fontWeight: 700, color: '#a9b7ff', marginRight: '1.5rem' }}>
                  R² = {regression.r2.toFixed(4)}
                </span>
                <span style={{ color: '#9ca3af' }}>
                  y = {regression.slope.toFixed(2)}x + {regression.intercept.toFixed(2)}
                </span>
              </div>
            </div>
            
            <ResponsiveContainer width="100%" height={500}>
              <ScatterChart margin={{ top: 20, right: 30, bottom: 60, left: 60 }}>
                <CartesianGrid strokeDasharray="3 3" stroke="rgba(255, 255, 255, 0.1)" />
                <XAxis 
                  type="number"
                  dataKey="modRate"
                  name="Modification Rate"
                  tickFormatter={(value) => `${(value * 100).toFixed(1)}%`}
                  label={{ 
                    value: `${currentModType.label} (%)`,
                    position: 'insideBottom',
                    offset: -10,
                    style: { fontSize: '14px', fontWeight: 700, fill: '#c0c0c0' }
                  }}
                  tick={{ fontSize: 12, fill: '#a0a0a0' }}
                  stroke="rgba(255, 255, 255, 0.3)"
                />
                <YAxis 
                  type="number"
                  dataKey="cpm"
                  name="CPM"
                  label={{ 
                    value: 'Counts Per Million (CPM)',
                    angle: -90,
                    position: 'insideLeft',
                    style: { fontSize: '14px', fontWeight: 700, fill: '#c0c0c0' }
                  }}
                  tick={{ fontSize: 12, fill: '#a0a0a0' }}
                  stroke="rgba(255, 255, 255, 0.3)"
                />
                <Tooltip 
                  content={({ active, payload }) => {
                    if (active && payload && payload.length) {
                      const data = payload[0].payload;
                      return (
                        <div style={{ 
                          background: 'rgba(26, 26, 46, 0.95)',
                          padding: '1rem',
                          border: '1px solid rgba(255, 255, 255, 0.2)',
                          borderRadius: '8px',
                          boxShadow: '0 4px 12px rgba(0, 0, 0, 0.3)'
                        }}>
                          <div style={{ fontWeight: 700, marginBottom: '0.5rem', color: '#fff' }}>
                            {data.name}
                          </div>
                          <div style={{ fontSize: '0.85rem', color: '#c0c0c0' }}>
                            <div>Rate: {(data.modRate * 100).toFixed(2)}%</div>
                            <div>CPM: {data.cpm.toFixed(2)}</div>
                          </div>
                        </div>
                      );
                    }
                    return null;
                  }}
                />
                <Scatter 
                  data={scatterData}
                  fill={currentRegion.color}
                  fillOpacity={0.7}
                  onClick={(data) => {
                    const gene = geneData.find(g => g.id === data.id);
                    setSelectedGene(gene);
                  }}
                  style={{ cursor: 'pointer' }}
                />
                <Line 
                  data={regressionLine}
                  type="monotone"
                  dataKey="y"
                  stroke="#ff6b9d"
                  strokeWidth={3}
                  dot={false}
                  strokeDasharray="5 5"
                />
              </ScatterChart>
            </ResponsiveContainer>
          </div>
        </div>
        
        {/* Gene Browser */}
        <div style={{ padding: '0 2.5rem 2.5rem 2.5rem' }}>
          <div style={{ 
            background: 'rgba(255, 255, 255, 0.03)',
            borderRadius: '12px',
            padding: '2rem',
            border: '1px solid rgba(255, 255, 255, 0.1)'
          }}>
            <div style={{ 
              display: 'flex', 
              justifyContent: 'space-between', 
              alignItems: 'center',
              marginBottom: '2rem'
            }}>
              <h2 style={{ margin: 0, fontSize: '1.5rem', fontWeight: 700, color: '#fff' }}>
                Gene Browser
              </h2>
              <div style={{ position: 'relative', width: '350px' }}>
                <Search 
                  size={20} 
                  style={{ 
                    position: 'absolute',
                    left: '14px',
                    top: '50%',
                    transform: 'translateY(-50%)',
                    color: '#888'
                  }}
                />
                <input
                  type="text"
                  placeholder="Search genes..."
                  value={searchTerm}
                  onChange={(e) => setSearchTerm(e.target.value)}
                  style={{
                    width: '80%',
                    padding: '0.75rem 1rem 0.75rem 2.75rem',
                    border: '1px solid rgba(255, 255, 255, 0.2)',
                    borderRadius: '8px',
                    fontSize: '0.95rem',
                    outline: 'none',
                    background: 'rgba(0, 0, 0, 0.3)',
                    color: '#e0e0e0',
                    fontFamily: '"Space Grotesk", sans-serif'
                  }}
                />
              </div>
            </div>
            
            <div style={{ 
              maxHeight: '280px', 
              overflowY: 'auto',
              border: '1px solid rgba(255, 255, 255, 0.1)',
              borderRadius: '8px',
              marginBottom: '2rem',
              background: 'rgba(0, 0, 0, 0.2)'
            }}>
              <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fill, minmax(150px, 1fr))', gap: '1px', background: 'rgba(255, 255, 255, 0.05)', padding: '1px' }}>
                {filteredGenes.map((gene) => (
                  <button
                    key={gene.id}
                    onClick={() => setSelectedGene(gene)}
                    style={{
                      padding: '0.75rem 1rem',
                      border: 'none',
                      background: selectedGene?.id === gene.id ? 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)' : 'rgba(0, 0, 0, 0.3)',
                      color: selectedGene?.id === gene.id ? 'white' : '#d0d0d0',
                      fontWeight: selectedGene?.id === gene.id ? 700 : 500,
                      cursor: 'pointer',
                      textAlign: 'left',
                      fontSize: '0.9rem',
                      fontFamily: '"Space Grotesk", sans-serif',
                      transition: 'all 0.2s'
                    }}
                  >
                    {gene.name}
                  </button>
                ))}
              </div>
            </div>
            
            {/* Raw Data Table from PKL */}
            {selectedGene && selectedGene.raw_data && (
              <div>
                <h3 style={{ 
                  margin: '0 0 1.5rem 0', 
                  fontSize: '1.3rem', 
                  fontWeight: 700, 
                  color: '#fff',
                  padding: '1rem',
                  background: 'rgba(102, 126, 234, 0.15)',
                  borderRadius: '8px',
                  border: '1px solid rgba(102, 126, 234, 0.3)'
                }}>
                  {selectedGene.name} • Raw Data (from PKL file)
                </h3>
                
                <div style={{ overflowX: 'auto', border: '1px solid rgba(255, 255, 255, 0.1)', borderRadius: '8px' }}>
                  <table style={{ 
                    width: '100%', 
                    borderCollapse: 'collapse',
                    fontSize: '0.9rem',
                    fontFamily: '"JetBrains Mono", monospace'
                  }}>
                    <thead style={{ background: 'rgba(102, 126, 234, 0.2)', borderBottom: '1px solid rgba(255, 255, 255, 0.15)' }}>
                      <tr>
                        <th style={{ padding: '1rem', textAlign: 'left', fontWeight: 700, color: '#a9b7ff' }}>Feature</th>
                        <th style={{ padding: '1rem', textAlign: 'left', fontWeight: 700, color: '#a9b7ff' }}>Modification</th>
                        <th style={{ padding: '1rem', textAlign: 'right', fontWeight: 700, color: '#a9b7ff' }}>
                          Count
                        </th>
                        <th style={{ padding: '1rem', textAlign: 'right', fontWeight: 700, color: '#a9b7ff' }}>
                          CPK
                        </th>
                        <th style={{ padding: '1rem', textAlign: 'right', fontWeight: 700, color: '#a9b7ff' }}>
                          MR (Rate)
                        </th>
                      </tr>
                    </thead>
                    <tbody>
                      {selectedGene.raw_data.map((row, idx) => (
                        <tr key={idx} style={{ 
                          background: idx % 2 === 0 ? 'rgba(0, 0, 0, 0.2)' : 'rgba(0, 0, 0, 0.3)',
                          borderBottom: '1px solid rgba(255, 255, 255, 0.05)'
                        }}>
                          <td style={{ padding: '0.75rem 1rem', fontWeight: 600, color: '#fff' }}>
                            {row.feature}
                          </td>
                          <td style={{ padding: '0.75rem 1rem', color: '#c0c0c0' }}>
                            {row.modification}
                          </td>
                          <td style={{ padding: '0.75rem 1rem', textAlign: 'right', color: '#d0d0d0' }}>
                            {row.count}
                          </td>
                          <td style={{ padding: '0.75rem 1rem', textAlign: 'right', color: '#d0d0d0' }}>
                            {row.cpk.toFixed(2)}
                          </td>
                          <td style={{ padding: '0.75rem 1rem', textAlign: 'right', color: '#ffd700', fontWeight: 600 }}>
                            {row.mr.toFixed(6)}
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
                
                {/* Summary Stats */}
                <div style={{ 
                  marginTop: '2rem',
                  padding: '2rem',
                  background: 'rgba(0, 0, 0, 0.3)',
                  borderRadius: '8px',
                  border: '1px solid rgba(255, 255, 255, 0.1)',
                  display: 'grid',
                  gridTemplateColumns: 'repeat(auto-fit, minmax(220px, 1fr))',
                  gap: '2rem'
                }}>
                  <div>
                    <div style={{ fontSize: '0.75rem', fontWeight: 700, color: '#888', marginBottom: '0.5rem', textTransform: 'uppercase', letterSpacing: '0.1em' }}>
                      Total A-to-I Rate
                    </div>
                    <div style={{ fontSize: '1.6rem', fontWeight: 700, color: '#FF6B9D' }}>
                      {((selectedGene.total_ai_rate || 0) * 100).toFixed(2)}%
                    </div>
                  </div>
                  <div>
                    <div style={{ fontSize: '0.75rem', fontWeight: 700, color: '#888', marginBottom: '0.5rem', textTransform: 'uppercase', letterSpacing: '0.1em' }}>
                      Total m6A Rate
                    </div>
                    <div style={{ fontSize: '1.6rem', fontWeight: 700, color: '#4ECDC4' }}>
                      {((selectedGene.total_m6a_rate || 0) * 100).toFixed(2)}%
                    </div>
                  </div>
                  <div>
                    <div style={{ fontSize: '0.75rem', fontWeight: 700, color: '#888', marginBottom: '0.5rem', textTransform: 'uppercase', letterSpacing: '0.1em' }}>
                      Total Either Rate
                    </div>
                    <div style={{ fontSize: '1.6rem', fontWeight: 700, color: '#C7CEEA' }}>
                      {((selectedGene.total_either_rate || 0) * 100).toFixed(2)}%
                    </div>
                  </div>
                  <div>
                    <div style={{ fontSize: '0.75rem', fontWeight: 700, color: '#888', marginBottom: '0.5rem', textTransform: 'uppercase', letterSpacing: '0.1em' }}>
                      Total CPM
                    </div>
                    <div style={{ fontSize: '1.6rem', fontWeight: 700, color: '#FFE66D' }}>
                      {(selectedGene.total_cpm || 0).toFixed(2)}
                    </div>
                  </div>
                </div>
              </div>
            )}
            
            {!selectedGene && (
              <div style={{
                padding: '4rem',
                textAlign: 'center',
                color: '#666',
                fontSize: '1.1rem',
                fontStyle: 'italic'
              }}>
                Select a gene to view raw data from PKL file
              </div>
            )}
          </div>
        </div>
      </div>
    </div>
  );
};

export default RNAModificationAnalyzer;