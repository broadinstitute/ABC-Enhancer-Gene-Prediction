import os
import os.path

import pandas as pd
from neighborhoods import run_count_reads
from tools import run_piped_commands


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
