import { useState, useEffect } from 'react'
import { supabase } from './supabaseClient'

function extractMean(val) {
  if (typeof val === 'number') return val
  if (typeof val === 'string') {
    const match = val.match(/^([0-9.]+)/)
    return match ? parseFloat(match[1]) : 0
  }
  return 0
}

function regionKey(feature) {
  return feature.toLowerCase().replace(/_/g, '')
}

function processGeneRows(rows, sample) {
  const cpkCol   = sample === 'MR01_1' ? 'cpk_mr01_1'  : 'cpk_mr01_2'
  const mrCol    = sample === 'MR01_1' ? 'mr01_1'       : 'mr01_2'
  const countCol = sample === 'MR01_1' ? 'count_mr01_1' : 'count_mr01_2'

  const regions = ['UTR_5', 'UTR_3', 'Exon', 'Intron']
  const stats   = {}
  const rawData = []

  // Initialize all keys to 0 so nothing is undefined
  for (const region of regions) {
    const rk = regionKey(region)
    stats[`${rk}_cpm`]         = 0
    stats[`${rk}_ai_rate`]     = 0
    stats[`${rk}_m6a_rate`]    = 0
    stats[`${rk}_unmod_rate`]  = 0
    stats[`${rk}_either_rate`] = 0
  }

  // Accumulate CPK values per region for averaging
  const cpkAccum = {}

  for (const row of rows) {
    const rk = regionKey(row.feature)

    // Build raw data table
    rawData.push({
      feature:      row.feature,
      modification: row.modification,
      count:        row[countCol] ?? 0,
      cpk:          row[cpkCol]   ?? 0,
      mr:           extractMean(row[mrCol]),
    })

    if (!regions.includes(row.feature)) continue

    // Accumulate CPK
    const cpk = row[cpkCol] ?? 0
    if (cpk > 0) {
      if (!cpkAccum[rk]) cpkAccum[rk] = []
      cpkAccum[rk].push(cpk)
    }

    // Rates
    const rate = extractMean(row[mrCol])
    if (row.modification === 'Inosine') stats[`${rk}_ai_rate`]    = rate
    if (row.modification === 'm6A')     stats[`${rk}_m6a_rate`]   = rate
    if (row.modification === 'Unmod')   stats[`${rk}_unmod_rate`] = rate
  }

  // Collapse CPK accumulators and compute either rate
  for (const region of regions) {
    const rk  = regionKey(region)
    const arr = cpkAccum[rk] ?? []
    stats[`${rk}_cpm`] = arr.length
      ? arr.reduce((a, b) => a + b, 0) / arr.length
      : 0

    stats[`${rk}_either_rate`] = Math.min(
      (stats[`${rk}_ai_rate`] ?? 0) + (stats[`${rk}_m6a_rate`] ?? 0),
      1.0
    )
  }

  // Totals — average across all regions
  const cpmVals = regions
    .map(r => stats[`${regionKey(r)}_cpm`])
    .filter(v => v > 0)

  stats.total_cpm = cpmVals.length
    ? cpmVals.reduce((a, b) => a + b, 0) / cpmVals.length
    : 0

  stats.total_ai_rate = regions
    .map(r => stats[`${regionKey(r)}_ai_rate`] ?? 0)
    .reduce((a, b) => a + b, 0) / regions.length

  stats.total_m6a_rate = regions
    .map(r => stats[`${regionKey(r)}_m6a_rate`] ?? 0)
    .reduce((a, b) => a + b, 0) / regions.length

  stats.total_either_rate = Math.min(
    stats.total_ai_rate + stats.total_m6a_rate,
    1.0
  )

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

        // Return cached data immediately if available
        const cached = sessionStorage.getItem('geneData')
        if (cached) {
          const { genes1, genes2 } = JSON.parse(cached)
          setGeneDataMR01_1(genes1)
          setGeneDataMR01_2(genes2)
          setIsLoading(false)
          return
        }

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

        const byGene = {}
        for (const row of allRows) {
          if (!byGene[row.gene]) byGene[row.gene] = []
          byGene[row.gene].push(row)
        }

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

        // Cache for this session so reload is instant
        try {
          sessionStorage.setItem('geneData', JSON.stringify({ genes1, genes2 }))
        } catch (e) {
          // sessionStorage quota exceeded — skip caching, not fatal
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