# ABC Pipeline - Task Tracker

## Critical Issues

### Security
- [x] **Shell injection vulnerability** in `workflow/scripts/getVariantOverlap.py` - FIXED
  - Replaced `os.system()` with `run_piped_commands_safe()` in tools.py
  - Added `tests/test_tools.py` for validation

### Hard-coded Values
- [ ] **Hard-coded thresholds** in `workflow/scripts/getVariantOverlap.py:23-26`
  - Values 0.015 and 0.1 are not configurable
  - No documentation explaining these thresholds
- [ ] **Magic numbers** throughout `predictor.py`, `neighborhoods.py`, `metrics.py`
  - Example: `search_space_bins = math.ceil(5000 / hic_resolution)` (line ~120)
  - Distance percentiles 10th/90th in metrics without explanation

### Memory Management
- [ ] **Memory-intensive operations** in `workflow/scripts/makeAverageHiC.py:91`
  - `pd.concat()` with `join="outer"` can exhaust memory
  - Code comment acknowledges: "TO DO: This is very memory intensive!"
- [ ] **No chunking for large files** in `workflow/scripts/predict.py:274`
  - `pd.concat(all_putative_list)` loads entire dataset into memory

## High Priority Issues

### Missing Validation
- [x] **MACS2 output validation** - FIXED
  - Added `validate_macs_peaks()` in `workflow/scripts/peaks.py`
  - Fails with informative error if MACS2 output is empty
  - Warns if peak count is below threshold (1000)
  - QC stats written to `candidateRegions.qc.txt`
- [x] **Input file and index validation** - FIXED
  - Added `_validate_input_files_exist()` in `workflow/rules/utils.smk`
  - Validates BAM/tagAlign/H3K27ac/HiC files exist during config loading
  - Checks for required index files (.bai for BAM, .tbi for tagAlign.gz)
  - Provides clear error messages with instructions to create indexes
  - Index files added as Snakemake rule inputs for dependency tracking
- [ ] **Chromosome name mismatch not detected**
  - Silent data loss if genes/enhancers/HiC use different naming (chr1 vs 1)
- [ ] **Hi-C file validity check is fragile** in `workflow/scripts/hic.py:49-50`
  - Uses `os.path.getsize() > 100` as validity heuristic

### Error Handling
- [ ] **Incomplete retry mechanism** in `workflow/rules/utils.smk`
  - `enable_retry()` only wraps `load_biosamples_config()`
  - Other file I/O lacks retry logic
- [ ] **Inconsistent exception handling styles**
  - Mix of try/except, assertions, and silent failures across scripts

### Reference Data Dependencies
- [ ] **Implicit qnorm file dependency**
  - `EnhancersQNormRef.K562.txt` required when `use_qnorm: True`
  - No clear error if file missing; may silently skip normalization
- [ ] **abc_thresholds.tsv not validated**
  - Silently defaults to 0.02 if file missing/corrupted

## Medium Priority Issues

### Code Quality
- [ ] **No type hints** across most Python files
  - Dynamic dtype conversions can fail silently
- [ ] **Pandas chained assignment warning** only configured in `tools.py`
  - Inconsistent across other modules
- [ ] **Missing required column validation**
  - Hard-coded column names ("chr", "start", "end", "symbol") without existence checks

### Performance
- [ ] **Inefficient PyRanges join** in `predictor.py:93-95`
  - Creates full cartesian product before distance filtering
- [ ] **No progress indicators** for long-running operations
  - 24+ hour runs provide minimal feedback

### Documentation
- [ ] **Complex parameters undocumented** in config.yaml
  - `--scale_hic_using_powerlaw` flag purpose unclear
  - `hic_pseudocount_distance: 5000` - rationale not explained
- [ ] **No external tool version checking**
  - MACS2, bedtools, samtools versions not validated
- [ ] **`docs/developers/testing.rst` still references CircleCI**
  - CI migrated to GitHub Actions (`.github/workflows/ci.yml`, `nightly-dev.yml`); `.circleci/config.yml` removed
  - Lines 67-78 need updating to describe the GitHub Actions setup instead (workflow file location, badge, required status checks, nightly schedule)

## Low Priority / Future Improvements

- [ ] Add comprehensive test suite (none visible in workflow/)
- [ ] Implement streaming/chunking for large datasets
- [ ] Standardize logging across all scripts
- [ ] Add chromosome name normalization layer
- [ ] Replace `os.system()` calls with subprocess for security
- [ ] Add pre-flight validation rule in Snakefile

---
*Last updated: 2026-08-14*
