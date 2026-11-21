# ABC Enhancer-Gene Prediction Pipeline

## Overview
Snakemake-based pipeline computing Activity-by-Contact (ABC) scores to predict enhancer-gene regulatory connections.

**Core formula**:
```
ABC Score = (Activity × Contact) / Σ(Activity × Contact) for all enhancers of that gene
```
- **Activity** = `activity_base_enh` (enhancer activity from accessibility, quantile normalized)
- **Contact** = `hic_contact_pl_scaled_adj` (Hi-C with powerlaw scaling + pseudocount) or `powerlaw_contact` (fallback)

The score represents the fraction of total regulatory input each enhancer contributes to a gene.

## Directory Structure
```
workflow/
├── Snakefile              # Main entry point
├── rules/                 # Snakemake rules (macs2, candidate_regions, neighborhoods, predictions, qc)
├── scripts/               # Python implementation (~17 scripts)
└── envs/                  # Conda environments (abcenv.yml)
```

## Key Files by Function

| Purpose | Files |
|---------|-------|
| Peak calling | `rules/macs2.smk`, external MACS2 |
| Candidate regions | `rules/candidate_regions.smk`, `scripts/peaks.py`, `scripts/makeCandidateRegions.py` |
| Neighborhoods | `rules/neighborhoods.smk`, `scripts/neighborhoods.py` (29KB, most complex), `scripts/run.neighborhoods.py` |
| Predictions | `rules/predictions.smk`, `scripts/predictor.py` (19KB, core engine), `scripts/predict.py` |
| Filtering | `scripts/filter_predictions.py` |
| Hi-C handling | `scripts/hic.py`, `scripts/juicebox_dump.py`, `scripts/makeAverageHiC.py` |
| QC/Metrics | `rules/qc.smk`, `scripts/grabMetrics.py`, `scripts/metrics.py` |
| Utilities | `rules/utils.smk`, `scripts/tools.py` |

## Configuration
- **Main config**: `config/config.yaml` - paths, reference files, rule parameters
- **Samples**: `config/config_biosamples_*.tsv` - per-sample inputs (DHS/ATAC, H3K27ac, HiC)
- **Key params**: `peakExtendFromSummit` (250bp), `nStrongestPeaks` (150k), `hic_resolution` (5000bp only), `use_qnorm` (True critical for thresholds)

## Pipeline Flow
```
MACS2 peaks → Candidate regions → Neighborhoods → Predictions → Filtering → QC
```

## Output Files (per biosample)
```
{biosample}/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed
{biosample}/Peaks/candidateRegions.qc.txt          # QC stats for MACS2/candidate regions
{biosample}/Neighborhoods/EnhancerList.txt, GeneList.txt
{biosample}/Predictions/EnhancerPredictionsAllPutative.tsv.gz
{biosample}/Predictions/EnhancerPredictions_threshold{X}.tsv
{biosample}/Metrics/QCSummary_threshold{X}.tsv
```

## Known Issues
See `TODO.md` for tracked issues.

## Running
```bash
snakemake -s workflow/Snakefile --use-conda
```

## Dependencies
BioConda: snakemake, bedtools, samtools, macs2, pysam
Python: pandas, numpy, scipy, pyranges, hic-straw, matplotlib
