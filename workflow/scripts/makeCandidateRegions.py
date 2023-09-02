import argparse
import os

from peaks import make_candidate_regions_from_peaks, make_candidate_regions_from_summits
from tools import write_params


def parseargs(required_args=True):
    class formatter(
        argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter
    ):
        pass

    epilog = ""
    parser = argparse.ArgumentParser(
        description="Make peaks file for a given cell type",
        epilog=epilog,
        formatter_class=formatter,
    )
    readable = argparse.FileType("r")

    parser.add_argument(
        "--narrowPeak",
        required=required_args,
        help="narrowPeak file output by macs2. Must include summits (--call-summits)",
    )
    parser.add_argument(
        "--bam",
        required=required_args,
        nargs="?",
        help="DNAase-Seq or ATAC-Seq bam/tagalign file",
    )
    parser.add_argument(
        "--chrom_sizes",
        required=required_args,
        help="File listing chromosome size annotaions",
    )
    parser.add_argument(
        "--chrom_sizes_bed",
        required=required_args,
        help="File listing chromosome size annotaions",
    )
    parser.add_argument("--outDir", required=required_args)

    parser.add_argument(
        "--nStrongestPeaks",
        default=175000,
        help="Number of peaks to use for defining candidate regions",
    )
    parser.add_argument(
        "--peakExtendFromSummit",
        default=250,
        help="Number of base pairs to extend each preak from its summit (or from both ends of region if using --ignoreSummits)",
    )
    parser.add_argument(
        "--ignoreSummits",
        action="store_true",
        help="Compute peaks using the full peak regions, rather than extending from summit.",
    )
    parser.add_argument(
        "--minPeakWidth",
        default=500,
        help="Candidate regions whose width is below this threshold are expanded to this width. Only used with --ignoreSummits",
    )

    parser.add_argument(
        "--regions_includelist",
        default="",
        help="Bed file of regions to forcibly include in candidate enhancers. Overrides regions_blocklist",
    )
    parser.add_argument(
        "--regions_blocklist",
        default="",
        help="Bed file of regions to forcibly exclude from candidate enhancers",
    )

    args = parser.parse_args()
    return args


def processCellType(args):
    os.makedirs(os.path.join(args.outDir), exist_ok=True)
    write_params(args, os.path.join(args.outDir, "params.txt"))

    # Make candidate regions
    if not args.ignoreSummits:
        make_candidate_regions_from_summits(
            macs_peaks=args.narrowPeak,
            accessibility_file=args.bam,
            genome_sizes=args.chrom_sizes,
            genome_sizes_bed=args.chrom_sizes_bed,
            regions_includelist=args.regions_includelist,
            regions_blocklist=args.regions_blocklist,
            n_enhancers=args.nStrongestPeaks,
            peak_extend=args.peakExtendFromSummit,
            outdir=args.outDir,
        )
    else:
        make_candidate_regions_from_peaks(
            macs_peaks=args.narrowPeak,
            accessibility_file=args.bam,
            genome_sizes=args.chrom_sizes,
            genome_sizes_bed=args.chrom_sizes_bed,
            regions_includelist=args.regions_includelist,
            regions_blocklist=args.regions_blocklist,
            n_enhancers=args.nStrongestPeaks,
            peak_extend=args.peakExtendFromSummit,
            minPeakWidth=args.minPeakWidth,
            outdir=args.outDir,
        )


def main(args):
    processCellType(args)


if __name__ == "__main__":
    args = parseargs()
    main(args)
