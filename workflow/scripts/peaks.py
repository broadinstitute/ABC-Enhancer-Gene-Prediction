import os
import os.path

import pandas as pd
from neighborhoods import run_count_reads
from tools import run_piped_commands

# Minimum peaks threshold - below this we warn users
MIN_PEAKS_WARNING = 1000


def validate_macs_peaks(macs_peaks, min_peaks=MIN_PEAKS_WARNING):
    """
    Validate MACS2 output has sufficient peaks.

    Returns number of peaks found.
    Raises ValueError if file is empty.
    Prints warning if below threshold.
    """
    if not os.path.exists(macs_peaks):
        raise ValueError(f"MACS2 output file not found: {macs_peaks}")

    n_peaks = sum(1 for _ in open(macs_peaks))

    if n_peaks == 0:
        raise ValueError(
            f"MACS2 output is empty: {macs_peaks}\n"
            "Possible causes:\n"
            "  - Input BAM/tagAlign files have no reads\n"
            "  - Input files are corrupted or in wrong format\n"
            "  - MACS2 p-value threshold is too stringent\n"
            "Check your input accessibility files and MACS2 logs."
        )

    if n_peaks < min_peaks:
        print(
            f"WARNING: Only {n_peaks} peaks found in MACS2 output (expected >{min_peaks}).\n"
            f"  File: {macs_peaks}\n"
            "  Candidate regions may be dominated by promoter regions from the includelist.\n"
            "  Consider checking input data quality or adjusting MACS2 parameters."
        )

    return n_peaks


def write_candidate_regions_qc(outdir, n_macs_peaks, n_candidate_regions):
    """Write QC stats for candidate regions to a text file."""
    qc_file = os.path.join(outdir, "candidateRegions.qc.txt")

    with open(qc_file, "w") as f:
        f.write("=== Candidate Regions QC Summary ===\n\n")
        f.write(f"MACS2 peaks (input):          {n_macs_peaks:,}\n")
        f.write(f"Candidate regions (output):   {n_candidate_regions:,}\n")

        # Add warnings if concerning
        if n_candidate_regions == 0:
            f.write("\n*** ERROR: No candidate regions produced! ***\n")
        elif n_macs_peaks < MIN_PEAKS_WARNING:
            f.write(
                f"\n*** WARNING: Low number of MACS2 peaks ({n_macs_peaks:,}). ***\n"
                "Candidate regions may be dominated by promoters from the includelist.\n"
                "Consider checking input data quality.\n"
            )

    print(f"QC summary written to: {qc_file}")
    return qc_file


def count_lines(filepath):
    """Count lines in a file."""
    if not os.path.exists(filepath):
        return 0
    return sum(1 for _ in open(filepath))


def make_candidate_regions_from_summits(
    macs_peaks,
    accessibility_files,
    genome_sizes,
    genome_sizes_bed,
    regions_includelist,
    regions_blocklist,
    n_enhancers,
    peak_extend,
    outdir,
):
    ## Generate enhancer regions from MACS summits: 1. Count reads in DHS peaks 2. Take top N regions, get summits, extend summits, merge
    outfile = os.path.join(
        outdir, os.path.basename(macs_peaks) + ".candidateRegions.bed"
    )

    # Validate MACS2 output before proceeding
    n_macs_peaks = validate_macs_peaks(macs_peaks)

    includelist_command = get_includelist_command(regions_includelist, genome_sizes_bed)
    blocklist_command = get_blocklist_command(regions_blocklist)

    # 1. Count DHS/ATAC reads in candidate regions for all accessibility files provided, and return the filename of the average # reads
    reads_out = get_read_counts(
        accessibility_files, outdir, macs_peaks, genome_sizes, genome_sizes_bed
    )

    # 2. Take top N regions, get summits, extend summits, merge, remove blocklist, add includelist, sort and merge
    # use -sorted in intersect command? Not worth it, both files are small
    piped_cmds = [
        f"bedtools sort -i {reads_out} -faidx {genome_sizes}",
        "bedtools merge -i stdin -c 4 -o max",
        "sort -nr -k 4",
        f"head -n {n_enhancers}",
        f"bedtools intersect -b stdin -a {macs_peaks} -wa",
        'awk \'{{print $1 "\\t" $2 + $10 "\\t" $2 + $10}}\'',
        f"bedtools slop -i stdin -b {peak_extend} -g {genome_sizes}",
        f"bedtools sort -i stdin -faidx {genome_sizes}",
        "bedtools merge -i stdin",
        blocklist_command,
        "cut -f 1-3",
        includelist_command,
        f"bedtools sort -i stdin -faidx {genome_sizes}",
        f"bedtools merge -i stdin > {outfile}",
    ]

    run_piped_commands(piped_cmds)

    # Write QC stats
    n_candidate_regions = count_lines(outfile)
    write_candidate_regions_qc(outdir, n_macs_peaks, n_candidate_regions)


def make_candidate_regions_from_peaks(
    macs_peaks,
    accessibility_files,
    genome_sizes,
    genome_sizes_bed,
    regions_includelist,
    regions_blocklist,
    n_enhancers,
    peak_extend,
    minPeakWidth,
    outdir,
):
    ## Generate enhancer regions from MACS narrowPeak - do not use summits
    outfile = os.path.join(
        outdir, os.path.basename(macs_peaks) + ".candidateRegions.bed"
    )

    # Validate MACS2 output before proceeding
    n_macs_peaks = validate_macs_peaks(macs_peaks)

    includelist_command = get_includelist_command(regions_includelist, genome_sizes_bed)
    blocklist_command = get_blocklist_command(regions_blocklist)

    # 1. Count DHS/ATAC reads in candidate regions
    reads_out = get_read_counts(
        accessibility_files, outdir, macs_peaks, genome_sizes, genome_sizes_bed
    )

    # 2. Take top N regions, extend peaks (min size 500), merge, remove blocklist, add includelist, sort and merge
    # use -sorted in intersect command? Not worth it, both files are small
    piped_cmds = [
        f"bedtools sort -i {reads_out} -faidx {genome_sizes}",
        f"bedtools merge -i stdin -c 4 -o max",
        "sort -nr -k 4",
        f"head -n {n_enhancers}",
        f"bedtools intersect -b stdin -a {macs_peaks} -wa",
        f"bedtools slop -i stdin -b {peak_extend} -g {genome_sizes}",
        f'awk \'{{ l=$3-$2; if (l < {minPeakWidth}) {{ $2 = $2 - int(({minPeakWidth}-l)/2); $3 = $3 + int(({minPeakWidth}-l)/2) }} print $1 "\\t" $2 "\\t" $3}}\'',
        f"bedtools sort -i stdin -faidx {genome_sizes}",
        "bedtools merge -i stdin",
        blocklist_command,
        "cut -f 1-3",
        includelist_command,
        f"bedtools sort -i stdin -faidx {genome_sizes} | bedtools merge -i stdin > {outfile}",
    ]

    run_piped_commands(piped_cmds)

    # Write QC stats
    n_candidate_regions = count_lines(outfile)
    write_candidate_regions_qc(outdir, n_macs_peaks, n_candidate_regions)


def get_includelist_command(regions_includelist, genome_sizes_bed):
    if regions_includelist:
        return f"(bedtools intersect -a {regions_includelist} -b {genome_sizes_bed} -wa | cut -f 1-3 && cat)"
    else:
        return ""


def get_blocklist_command(regions_blocklist):
    if regions_blocklist:
        return f"bedtools intersect -v -wa -a stdin -b {regions_blocklist}"
    else:
        return ""


def get_read_counts(
    accessibility_files, outdir, macs_peaks, genome_sizes, genome_sizes_bed
):
    raw_counts_out = []  # initialize list for output file names
    for access_in in accessibility_files:  # loop through input accessibilty files
        raw_counts_out.append(
            os.path.join(
                outdir,
                os.path.basename(macs_peaks)
                + "."
                + os.path.basename(access_in)
                + ".Counts.bed",
            )
        )

    # 1. Count DHS/ATAC reads in candidate regions for all accessibility files provided, and return the filename of the average # reads
    reads_out = count_reads_over_peaks(
        accessibility_files,
        raw_counts_out,
        macs_peaks,
        genome_sizes,
        genome_sizes_bed,
        outdir,
        use_fast_count=True,
    )
    return reads_out


# count reads over however many DHS files and return average
def count_reads_over_peaks(
    accessibility_files,
    raw_counts_out,
    macs_peaks,
    genome_sizes,
    genome_sizes_bed,
    outdir,
    use_fast_count=True,
):
    for access_in, counts_out in zip(accessibility_files, raw_counts_out):
        run_count_reads(
            access_in,
            counts_out,
            macs_peaks,
            genome_sizes,
            genome_sizes_bed,
            use_fast_count,
        )

    nFiles = len(accessibility_files)
    if nFiles > 1:
        avg_out = os.path.join(
            outdir, os.path.basename(macs_peaks) + ".averageAccessibility.Counts.bed"
        )
        col_names = ["chrom", "start", "end", "count"]
        df1 = pd.read_csv(raw_counts_out[0], sep="\t", names=col_names)
        for i in range(1, nFiles):
            dfx = pd.read_csv(
                raw_counts_out[i], sep="\t", names=col_names, usecols=["count"]
            )
            df1["count"] = df1["count"].add(dfx["count"])
        df1["count"] = df1["count"] / nFiles
        df1.to_csv(avg_out, header=None, index=None, sep="\t")
        return avg_out
    else:
        return raw_counts_out[0]
