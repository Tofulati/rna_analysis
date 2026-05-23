import { useState, useEffect } from 'react'
import { supabase } from './supabaseClient'

// Mirrors what load_analyze.py does locally —
// fetches all rows for both samples and reshapes into
// the same structure the React component expects.

function extractMean(val) {
  if (typeof val === 'number') return val
  if (typeof val === 'string') {
    const match = val.match(/^([0-9.]+)/)
    return match ? parseFloat(match[1]) : 0
  }
  return 0
}

function processGeneRows(rows, sample) {
  // rows = all DB rows for one gene, both samples present in each row
  const cpkCol    = sample === 'MR01_1' ? 'cpk_mr01_1'   : 'cpk_mr01_2'
  const mrCol     = sample === 'MR01_1' ? 'mr01_1'        : 'mr01_2'
  const countCol  = sample === 'MR01_1' ? 'count_mr01_1'  : 'count_mr01_2'

  const regions = ['UTR_5', 'UTR_3', 'Exon', 'Intron']
  const mods    = ['Unmod', 'm6A', 'Inosine']

  const stats = {}
  const rawData = []

  for (const row of rows) {
    const regionKey = row.feature.toLowerCase().replace('_', '')

    // raw data table (mirrors raw_data array in old JSON)
    rawData.push({
      feature:      row.feature,
      modification: row.modification,
      count:        row[countCol] ?? 0,
      cpk:          row[cpkCol]   ?? 0,
      mr:           extractMean(row[mrCol]),
    })

    if (!regions.includes(row.feature)) continue

    // CPM
    const cpkKey = `${regionKey}_cpm`
    if (!stats[cpkKey]) stats[cpkKey] = []
    if ((row[cpkCol] ?? 0) > 0) stats[cpkKey].push(row[cpkCol])

    // Rates
    const rate = extractMean(row[mrCol])
    if (row.modification === 'Inosine') stats[`${regionKey}_ai_rate`]   = rate
    if (row.modification === 'm6A')     stats[`${regionKey}_m6a_rate`]  = rate
    if (row.modification === 'Unmod')   stats[`${regionKey}_unmod_rate`] = rate
  }

  // Collapse CPM arrays to means
  for (const region of regions) {
    const rk     = region.toLowerCase().replace('_', '')
    const cpmKey = `${rk}_cpm`
    stats[cpmKey] = stats[cpmKey]?.length
      ? stats[cpmKey].reduce((a, b) => a + b, 0) / stats[cpmKey].length
      : 0

    // Either rate
    const ai  = stats[`${rk}_ai_rate`]  ?? 0
    const m6a = stats[`${rk}_m6a_rate`] ?? 0
    stats[`${rk}_either_rate`] = Math.min(ai + m6a, 1.0)
  }

  // Totals
  const cpmVals  = regions.map(r => stats[`${r.toLowerCase().replace('_','')} _cpm`]).filter(v => v > 0)
  stats.total_cpm       = cpmVals.length ? cpmVals.reduce((a,b)=>a+b,0)/cpmVals.length : 0
  stats.total_ai_rate   = regions.map(r=>stats[`${r.toLowerCase().replace('_','')} _ai_rate`]??0).reduce((a,b)=>a+b,0)/regions.length
  stats.total_m6a_rate  = regions.map(r=>stats[`${r.toLowerCase().replace('_','')} _m6a_rate`]??0).reduce((a,b)=>a+b,0)/regions.length
  stats.total_either_rate = Math.min(stats.total_ai_rate + stats.total_m6a_rate, 1.0)

  return { stats, rawData }
}

export function useGeneData() {
  const [geneDataMR01_1, setGeneDataMR01_1] = useState([])
  const [geneDataMR01_2, setGeneDataMR01_2] = useState([])
  const [isLoading, setIsLoading]           = useState(true)
  const [loadError, setLoadError]           = useState(null)

  useEffect(() => {
    async function load() {
      try {
        setIsLoading(true)

        // Fetch all rows in batches (Supabase default limit is 1000)
        let allRows = []
        let offset  = 0
        const batchSize = 1000

        while (true) {
          const { data, error } = await supabase
            .from('gene_modifications')
            .select('*')
            .range(offset, offset + batchSize - 1)

          if (error) throw error
          if (!data || data.length === 0) break

          allRows = allRows.concat(data)
          if (data.length < batchSize) break
          offset += batchSize
        }

        // Group rows by gene
        const byGene = {}
        for (const row of allRows) {
          if (!byGene[row.gene]) byGene[row.gene] = []
          byGene[row.gene].push(row)
        }

        // Build per-sample gene arrays (mirrors old JSON structure)
        const genes1 = []
        const genes2 = []
        let id = 1

        for (const [geneName, rows] of Object.entries(byGene).sort()) {
          const { stats: s1, rawData: rd1 } = processGeneRows(rows, 'MR01_1')
          const { stats: s2, rawData: rd2 } = processGeneRows(rows, 'MR01_2')

          genes1.push({ id, name: geneName, raw_data: rd1, ...s1 })
          genes2.push({ id, name: geneName, raw_data: rd2, ...s2 })
          id++
        }

        setGeneDataMR01_1(genes1)
        setGeneDataMR01_2(genes2)
      } catch (err) {
        console.error('Supabase load error:', err)
        setLoadError(err.message)
      } finally {
        setIsLoading(false)
      }
    }

    load()
  }, [])

  return { geneDataMR01_1, geneDataMR01_2, isLoading, loadError }
}