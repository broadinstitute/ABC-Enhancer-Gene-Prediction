# HANDOVER — H3K27ac counting-mode experiment

Testing whether changing how H3K27ac ChIP-seq reads are counted over
candidate enhancer regions (element / element+flank / flank-only, at
multiple flank widths) changes ABC prediction performance against the K562
CRISPR benchmark, for DNase and ATAC accessibility separately.

## Status snapshot

| Area | Status |
|---|---|
| Code changes (`count_bam`, CLI flags, config plumbing) | Done, on branch `h3k27ac_counting_modes` — **not yet committed or merged** (working tree still shows the same uncommitted diff) |
| Unit tests (`tests/test_neighborhoods_counting.py`) | Done, 4/4 passing |
| Regression check (`tests/test_full_abc_run.py`, chr22) | Passing — confirms new params don't change existing behavior |
| ABC run — 12 K562 biosamples (DNase/ATAC × none/element/element_flank@150,200/flank_only@150,200) | Complete — `results/2026_0721_h3k27ac_counting_comparison/` |
| ABC run — 4 additional 500bp-flank biosamples (element_flank@500, flank_only@500, DNase + ATAC) | Complete — added to the same biosamples TSV/results dir |
| CRISPR benchmark — combined DNase+ATAC report (12 predictors, superseded) | Complete but superseded — `CRISPR_comparison/workflow/results/2026_0721_h3k27ac_counting_comparison/` |
| CRISPR benchmark — split DNase-only report (8 predictors incl. 500bp) | Complete — `CRISPR_comparison/workflow/results/2026_0721_h3k27ac_counting_DHS/2026_0721_h3k27ac_counting_DHS_crispr_comparison.html` |
| CRISPR benchmark — split ATAC-only report (8 predictors incl. 500bp) | Complete — `CRISPR_comparison/workflow/results/2026_0721_h3k27ac_counting_ATAC/2026_0721_h3k27ac_counting_ATAC_crispr_comparison.html` |
| **Interpretation of results (which counting mode wins)** | **Not yet done** — reports are generated but not reviewed |

## Key files (this repo)

- `workflow/scripts/neighborhoods.py` — `count_bam`/`run_count_reads`/`count_features_for_bed`/`load_enhancers` gained `count_mode` (`element`/`element_flank`/`flank_only`) and `flank_bp`, applied only to the H3K27ac feature.
- `workflow/scripts/run.neighborhoods.py` — new `--h3k27ac_counting_mode`/`--h3k27ac_flank_bp` CLI flags.
- `workflow/rules/neighborhoods.smk`, `workflow/rules/utils.smk` — new optional per-biosample TSV columns `h3k27ac_counting_mode`/`h3k27ac_flank_bp`, backward-compatible (default to `element`/`0` if absent).
- `config/mine/config_biosamples_h3k27ac_counting.tsv` — 16 K562 biosample rows (2 accessibility assays × 8 H3K27ac variants: none, element, element_flank@150/200/500, flank_only@150/200/500). Real ENCODE file paths (`ENCFF205FNC`/`ENCFF860XAE` DNase, `ENCFF790GFL`/`ENCFF817HMW` H3K27ac, `ENCFF077FBI`/`ENCFF128WZG`/`ENCFF534DCE` ATAC).
- `config/mine/config_h3k27ac_counting.yaml` — `use_qnorm: False` (deliberate — qnorm assumes cross-cell-type comparison, irrelevant here since everything is K562), `results_dir: results/2026_0721_h3k27ac_counting_comparison`.
- `config/mine/submit_h3k27ac_counting.sbatch` + `config/mine/slurm_profile/` — driver submits per-rule jobs sized by each rule's own `resources.mem_mb` rather than one fixed allocation.

## Key files (`CRISPR_comparison_v3/CRISPR_comparison`)

- `config/config_h3k27ac_counting.yml` — **two** comparisons: `2026_0721_h3k27ac_counting_DHS` and `2026_0721_h3k27ac_counting_ATAC` (split per your request; an earlier combined 12-predictor comparison of the same name minus `_DHS`/`_ATAC` also exists but is superseded).
- `resources/pred_config/pred_config_h3k27ac_counting_DHS.tsv` / `..._ATAC.tsv` — DHS uses blue (`element_flank` family) / teal (`flank_only` family); ATAC uses orange / purple. `alpha` (score cutoff) intentionally left blank for all H3K27ac variants since qnorm is off and the pre-calibrated `abc_thresholds.tsv` values assume qnorm'd data.
- `resources/cell_type_mapping/K562_*.txt` — one file per predictor, maps benchmark's `K562` label to the matching ABC biosample name (`predict.py --cellType {biosample}`).
- `workflow/submit_h3k27ac_counting.sbatch` — must be run with cwd = `workflow/` (Snakefile's `configfile:` directive is relative to invocation dir, not repo root).

## Open items / next steps

1. **Review the two HTML reports** (`performance_summary.txt`/`pr_table.txt.gz` alongside each) to see which H3K27ac counting mode performs best, for DNase and ATAC separately.
2. **Commit and push** the `h3k27ac_counting_modes` branch — the code changes are still only in the working tree, not committed. Consider whether to open a PR against `main` (or wherever this branch should land — it was branched off `add_ref_qc`).
3. Known inefficiency (not fixed, out of scope so far): each of the 16 biosample rows reruns MACS2/candidate-regions independently even though only 2 unique accessibility datasets (DNase, ATAC) underlie all of them — 8x redundant compute per assay. Fine for a one-off analysis at this scale; would matter if this pattern is repeated for a larger sweep.
4. Sherlock/Snakemake gotcha discovered while running the CRISPR_comparison benchmark (classic `--cluster` generic mode conflicting with `SLURM_TRES_PER_TASK`, fixed by switching to Snakemake 9 + the native `slurm` executor plugin) is now documented in the global `CLAUDE.md` (`Snakemake on Sherlock (HPC)` section) for reuse in other projects.

*Last updated: 2026-07-22.*
