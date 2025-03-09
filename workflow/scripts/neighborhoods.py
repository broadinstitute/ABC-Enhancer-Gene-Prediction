import linecache
import os
import os.path
import time
import traceback
from subprocess import PIPE, Popen, check_call, check_output

import numpy as np
import pandas as pd
import pyranges as pr
from pyBigWig import open as open_bigwig
import pysam
from scipy import interpolate
from tools import df_to_pyranges, run_command, run_piped_commands


pd.options.display.max_colwidth = (
    10000  # seems to be necessary for pandas to read long file names... strange
)

BED3_COLS = ["chr", "start", "end"]
BED6_COLS = BED3_COLS + ["name", "score", "strand"]


def read_gene_bed_file(bed_file):
    # BED6 format with extra columns for ensemble info
    columns = BED6_COLS + ["Ensembl_ID", "gene_type"]
    skip = 1 if ("track" in open(bed_file, "r").readline()) else 0
    result = pd.read_table(
        bed_file, names=columns, header=None, skiprows=skip, comment="#"
    )
    ensembl_id_col = result.loc[:, "Ensembl_ID"]
    if ensembl_id_col.isna().all() or not ensembl_id_col.str.contains("EN").all():
        raise Exception("Gene file doesn't follow the correct format with Ensembl info")
    return result


def load_genes(
    file,
    ue_file,
    chrom_sizes,
    outdir,
    expression_table_list,
    gene_id_names,
    primary_id,
    cellType,
    class_gene_file,
):
    # Add
    bed = read_gene_bed_file(file)
    genes = process_gene_bed(bed, gene_id_names, primary_id, chrom_sizes)

    genes[["chr", "start", "end", "name", "score", "strand"]].to_csv(
        os.path.join(outdir, "GeneList.bed"), sep="\t", index=False, header=False
    )

    if len(expression_table_list) > 0:
        # Add expression information
        names_list = []
        print("Using gene expression from files: {} \n".format(expression_table_list))

        for expression_table in expression_table_list:
            try:
                name = os.path.basename(expression_table)
                expr = pd.read_table(
                    expression_table, names=[primary_id, name + ".Expression"]
                )
                expr[name + ".Expression"] = expr[name + ".Expression"].astype(float)
                expr = expr.groupby(primary_id).max()

                genes = genes.merge(
                    expr, how="left", right_index=True, left_on="symbol"
                )
                names_list.append(name + ".Expression")
            except Exception as e:
                print(e)
                traceback.print_exc()
                print("Failed on {}".format(expression_table))

        genes["Expression"] = genes[names_list].mean(axis=1)
        genes["Expression.quantile"] = genes["Expression"].rank(
            method="average", na_option="top", ascending=True, pct=True
        )
    else:
        genes["Expression"] = np.NaN

    # Ubiquitously expressed annotation
    if ue_file is not None:
        ubiq = pd.read_csv(ue_file, sep="\t")
        genes["is_ue"] = genes["name"].isin(ubiq.iloc[:, 0].values.tolist())

    # cell type
    genes["cellType"] = cellType

    # genes for class assignment
    if class_gene_file is None:
        genes_for_class_assignment = genes
    else:
        genes_for_class_assignment = read_bed(class_gene_file)
        genes_for_class_assignment = process_gene_bed(
            genes_for_class_assignment,
            gene_id_names,
            primary_id,
            chrom_sizes,
            fail_on_nonunique=False,
        )

    return genes, genes_for_class_assignment


def annotate_genes_with_features(
    genes,
    genome_sizes,
    genome_sizes_bed,
    chrom_sizes_map,
    features={},
    outdir=".",
    use_fast_count=True,
    default_accessibility_feature="",
):
    # Setup files for counting
    bounds_bed = os.path.join(outdir, "GeneList.bed")
    tss1kb = make_tss_region_file(genes, outdir, genome_sizes, chrom_sizes_map)
    tss1kb_file = os.path.join(outdir, "GeneList.TSS1kb.bed")

    # Count features over genes and promoters
    genes = count_features_for_bed(
        genes,
        bounds_bed,
        genome_sizes,
        genome_sizes_bed,
        features,
        outdir,
        "Genes",
        use_fast_count=use_fast_count,
    )
    tsscounts = count_features_for_bed(
        tss1kb,
        tss1kb_file,
        genome_sizes,
        genome_sizes_bed,
        features,
        outdir,
        "Genes.TSS1kb",
        use_fast_count=use_fast_count,
    )
    tsscounts = tsscounts.drop(["chr", "start", "end", "score", "strand"], axis=1)

    merged = genes.merge(tsscounts, on="name", suffixes=["", ".TSS1Kb"])

    access_col = default_accessibility_feature + ".RPKM.quantile.TSS1Kb"

    if "H3K27ac.RPKM.quantile.TSS1Kb" in merged.columns:
        merged["PromoterActivityQuantile"] = (
            (0.0001 + merged["H3K27ac.RPKM.quantile.TSS1Kb"])
            * (0.0001 + merged[access_col])
        ).rank(method="average", na_option="top", ascending=True, pct=True)
    else:
        merged["PromoterActivityQuantile"] = ((0.0001 + merged[access_col])).rank(
            method="average", na_option="top", ascending=True, pct=True
        )

    merged.to_csv(
        os.path.join(outdir, "GeneList.txt"),
        sep="\t",
        index=False,
        header=True,
        float_format="%.6f",
    )

    return merged


def make_tss_region_file(genes, outdir, sizes, chrom_sizes_map, tss_slop=500):
    # Given a gene file, define 1kb regions around the tss of each gene

    tss1kb = genes.loc[:, ["chr", "start", "end", "name", "score", "strand"]]
    tss1kb["start"] = genes["tss"]
    tss1kb["end"] = genes["tss"]
    tss1kb = df_to_pyranges(tss1kb).slack(tss_slop)
    tss1kb = pr.gf.genome_bounds(tss1kb, chrom_sizes_map).df[
        ["Chromosome", "Start", "End", "name", "score", "strand"]
    ]
    tss1kb.columns = ["chr", "start", "end", "name", "score", "strand"]
    tss1kb_file = os.path.join(outdir, "GeneList.TSS1kb.bed")
    tss1kb.to_csv(tss1kb_file, header=False, index=False, sep="\t")

    # The TSS1kb file should be sorted
    sort_command = "bedtools sort -faidx {sizes} -i {tss1kb_file} > {tss1kb_file}.sorted; mv {tss1kb_file}.sorted {tss1kb_file}".format(
        **locals()
    )
    run_command(sort_command)

    return tss1kb


def process_gene_bed(
    bed, name_cols, main_name, chrom_sizes=None, fail_on_nonunique=True
):
    try:
        bed = bed.drop(
            [
                "thickStart",
                "thickEnd",
                "itemRgb",
                "blockCount",
                "blockSizes",
                "blockStarts",
            ],
            axis=1,
        )
    except Exception as e:
        pass

    assert main_name in name_cols

    names = bed.name.str.split(";", expand=True)
    assert len(names.columns) == len(name_cols.split(","))
    names.columns = name_cols.split(",")
    bed = pd.concat([bed, names], axis=1)

    bed["name"] = bed[main_name]
    # bed = bed.sort_values(by=['chr','start']) #JN Keep original sort order

    bed["tss"] = get_tss_for_bed(bed)

    bed.drop_duplicates(inplace=True)

    # Remove genes that are not defined in chromosomes file
    if chrom_sizes is not None:
        sizes = read_bed(chrom_sizes)
        bed["chr"] = bed["chr"].astype(
            "str"
        )  # JN needed in case chromosomes are all integer
        bed = bed[bed["chr"].isin(set(sizes["chr"].values))]

    # Enforce that gene names should be unique
    if fail_on_nonunique:
        assert len(set(bed["name"])) == len(
            bed["name"]
        ), "Gene IDs are not unique! Failing. Please ensure unique identifiers are passed to --genes"

    return bed


def get_tss_for_bed(bed):
    assert_bed3(bed)
    tss = bed["start"].copy()
    tss.loc[bed.loc[:, "strand"] == "-"] = bed.loc[bed.loc[:, "strand"] == "-", "end"]

    return tss


def assert_bed3(df):
    assert type(df).__name__ == "DataFrame"
    assert "chr" in df.columns
    assert "start" in df.columns
    assert "end" in df.columns
    assert "strand" in df.columns


def load_enhancers(
    outdir=".",
    genome_sizes="",
    genome_sizes_bed="",
    features={},
    genes=None,
    candidate_peaks="",
    skip_rpkm_quantile=False,
    cellType=None,
    tss_slop_for_class_assignment=500,
    use_fast_count=True,
    default_accessibility_feature="",
    qnorm=None,
    class_override_file=None,
    chrom_sizes_map=None,
):
    enhancers = read_bed(candidate_peaks)
    enhancers["chr"] = enhancers["chr"].astype("str")

    enhancers = count_features_for_bed(
        enhancers,
        candidate_peaks,
        genome_sizes,
        genome_sizes_bed,
        features,
        outdir,
        "Enhancers",
        skip_rpkm_quantile,
        use_fast_count,
    )

    # cellType
    if cellType is not None:
        enhancers["cellType"] = cellType

    # Assign categories
    if genes is not None:
        print("Assigning classes to enhancers")
        enhancers = assign_enhancer_classes(
            enhancers, genes, chrom_sizes_map, tss_slop=tss_slop_for_class_assignment
        )

    # TO DO: Should qnorm each bam file separately (before averaging). Currently qnorm being performed on the average
    enhancers = run_qnorm(enhancers, qnorm)
    enhancers = compute_activity(enhancers, default_accessibility_feature)

    enhancers[["chr", "start", "end", "name"]].to_csv(
        os.path.join(outdir, "EnhancerList.bed"),
        sep="\t",
        index=False,
        header=False,
    )
    enhancers.to_csv(
        os.path.join(outdir, "EnhancerList.txt.tmp"),
        sep="\t",
        index=False,
        header=True,
        float_format="%.6f",
    )
    os.rename(
        os.path.join(outdir, "EnhancerList.txt.tmp"),
        os.path.join(outdir, "EnhancerList.txt"),
    )


# Kristy's version
def assign_enhancer_classes(enhancers, genes, chrom_sizes_map, tss_slop=500):
    # build pyranges df
    tss_pyranges = df_to_pyranges(
        genes,
        start_col="tss",
        end_col="tss",
        start_slop=tss_slop,
        end_slop=tss_slop,
        chrom_sizes_map=chrom_sizes_map,
    )
    gene_pyranges = df_to_pyranges(genes)

    def get_class_pyranges(
        enhancers, tss_pyranges=tss_pyranges, gene_pyranges=gene_pyranges
    ):
        """
        Takes in PyRanges objects : Enhancers, tss_pyranges, gene_pyranges
        Returns dataframe with  uid (representing enhancer) and symbol of the gene/promoter that is overlapped
        """

        # genes
        genic_enh = enhancers.join(gene_pyranges, suffix="_genic")
        genic_enh = (
            genic_enh.df[["symbol", "uid"]]
            .groupby("uid", as_index=False)
            .aggregate(lambda x: ",".join(list(set(x))))
        )

        # promoters
        promoter_enh = enhancers.join(tss_pyranges, suffix="_promoter")
        promoter_enh = (
            promoter_enh.df[["symbol", "uid"]]
            .groupby("uid", as_index=False)
            .aggregate(lambda x: ",".join(list(set(x))))
        )

        return genic_enh, promoter_enh

    # label everything as intergenic
    enhancers["class"] = "intergenic"
    enhancers["uid"] = range(enhancers.shape[0])
    enh = df_to_pyranges(enhancers)

    genes, promoters = get_class_pyranges(enh)
    enhancers = enh.df.drop(["Chromosome", "Start", "End"], axis=1)
    enhancers.loc[enhancers["uid"].isin(genes.uid), "class"] = "genic"
    enhancers.loc[enhancers["uid"].isin(promoters.uid), "class"] = "promoter"

    enhancers["isPromoterElement"] = enhancers["class"] == "promoter"
    enhancers["isGenicElement"] = enhancers["class"] == "genic"
    enhancers["isIntergenicElement"] = enhancers["class"] == "intergenic"

    # Output stats
    print("Total enhancers: {}".format(len(enhancers)))
    print("         Promoters: {}".format(sum(enhancers["isPromoterElement"])))
    print("         Genic: {}".format(sum(enhancers["isGenicElement"])))
    print("         Intergenic: {}".format(sum(enhancers["isIntergenicElement"])))

    # Add promoter/genic symbol
    enhancers = enhancers.merge(
        promoters.rename(columns={"symbol": "promoterSymbol"}), on="uid", how="left"
    ).fillna(value={"promoterSymbol": ""})
    enhancers = enhancers.merge(
        genes.rename(columns={"symbol": "genicSymbol"}), on="uid", how="left"
    ).fillna(value={"genicSymbol": ""})
    enhancers.drop(["uid"], axis=1, inplace=True)

    # just to keep things consistent with original code
    enhancers["name"] = enhancers.apply(
        lambda e: "{}|{}:{}-{}".format(e["class"], e.chr, e.start, e.end), axis=1
    )
    return enhancers


def run_count_reads(
    target, output, bed_file, genome_sizes, genome_sizes_bed, use_fast_count
):
    filename = os.path.basename(target)
    if filename.endswith(".bam"):
        count_bam(
            target,
            bed_file,
            output,
            genome_sizes=genome_sizes,
            use_fast_count=use_fast_count,
        )
    elif "tagAlign" in filename:
        count_tagalign(target, bed_file, output, genome_sizes, genome_sizes_bed)
    elif isBigWigFile(filename):
        count_bigwig(target, bed_file, output)
    else:
        raise ValueError(
            "File {} name format doesn't match bam, tagAlign, or bigWig".format(target)
        )
    double_sex_chrom_counts(output)


def double_sex_chrom_counts(output):
    # Double the count values for sex chromosomes to make it seem
    # like they have 2 copies
    awk_command = r"""awk 'BEGIN {FS=OFS="\t"} (substr($1, length($1)) == "X" || substr($1, length($1)) == "Y") { $4 *= 2 } 1' """
    file_creation_command = f"{output} > {output}.tmp && mv {output}.tmp {output}"
    run_command(awk_command + file_creation_command)


def count_bam(
    bamfile, bed_file, output, genome_sizes, use_fast_count=True, verbose=True
):
    reads = pysam.AlignmentFile(bamfile)
    read_chrs = set(reads.references)
    bed_regions = pd.read_table(bed_file, header=None)
    bed_regions = bed_regions[bed_regions.columns[:3]]
    bed_regions.columns = "chr start end".split()
    counts = [
        (reads.count(row.chr, row.start, row.end) if (row.chr in read_chrs) else 0)
        for _, row in bed_regions.iterrows()
    ]
    bed_regions["count"] = counts
    bed_regions.to_csv(output, header=None, index=None, sep="\t")


def count_tagalign(tagalign, bed_file, output, genome_sizes, genome_sizes_bed):
    index_file = tagalign + ".tbi"
    if not os.path.exists(index_file):
        cmd = f"tabix -p bed {tagalign}"
        run_command(cmd)

    remove_alt_chr_cmd = f"bedtools intersect -u -a {tagalign} -b {genome_sizes_bed}"
    coverage_cmd = (
        f"bedtools coverage -counts -sorted -g {genome_sizes} -b - -a {bed_file}"
    )
    awk_cmd = 'awk \'{{print $1 "\\t" $2 "\\t" $3 "\\t" $NF}}\'' + f" > {output}"
    piped_cmds = [remove_alt_chr_cmd, coverage_cmd, awk_cmd]
    run_piped_commands(piped_cmds)


def count_bigwig(target, bed_file, output):
    bw = open_bigwig(target)
    bed = read_bed(bed_file)
    with open(output, "wb") as outfp:
        for chr, start, end, *rest in bed.itertuples(index=False, name=None):
            # if isinstance(name, np.float):
            #     name = ""
            try:
                if chr not in bw.chroms():
                    val = 0
                else:
                    val = (
                        bw.stats(chr, int(start), int(end), type="sum", exact=True)[0]
                        or 0
                    )
            except RuntimeError:
                print("Failed on", chr, start, end)
                raise
            output = ("\t".join([chr, str(start), str(end), str(val)]) + "\n").encode(
                "ascii"
            )
            outfp.write(output)


def isBigWigFile(filename):
    return (
        filename.endswith(".bw")
        or filename.endswith(".bigWig")
        or filename.endswith(".bigwig")
    )


def count_features_for_bed(
    df,
    bed_file,
    genome_sizes,
    genome_sizes_bed,
    features,
    directory,
    filebase,
    skip_rpkm_quantile=False,
    use_fast_count=True,
):
    for feature, feature_bam_list in features.items():
        start_time = time.time()
        if isinstance(feature_bam_list, str):
            feature_bam_list = [feature_bam_list]

        for feature_bam in feature_bam_list:
            df = count_single_feature_for_bed(
                df,
                bed_file,
                genome_sizes,
                genome_sizes_bed,
                feature_bam,
                feature,
                directory,
                filebase,
                skip_rpkm_quantile,
                use_fast_count,
            )

        df = average_features(
            df, feature.replace("feature_", ""), feature_bam_list, skip_rpkm_quantile
        )
        elapsed_time = time.time() - start_time
        print("Feature " + feature + " completed in " + str(elapsed_time))

    return df


def count_single_feature_for_bed(
    df,
    bed_file,
    genome_sizes,
    genome_sizes_bed,
    feature_bam,
    feature,
    directory,
    filebase,
    skip_rpkm_quantile,
    use_fast_count,
):
    orig_df = df.copy()
    orig_shape = df.shape[0]
    feature_name = feature + "." + os.path.basename(feature_bam)
    feature_outfile = os.path.join(
        directory, "{}.{}.CountReads.bedgraph".format(filebase, feature_name)
    )

    print("Generating", feature_outfile)
    print("Counting coverage for {}".format(filebase + "." + feature_name))
    run_count_reads(
        feature_bam,
        feature_outfile,
        bed_file,
        genome_sizes,
        genome_sizes_bed,
        use_fast_count,
    )

    domain_counts = read_bed(feature_outfile)
    score_column = domain_counts.columns[-1]

    total_counts = count_total(feature_bam)

    domain_counts = domain_counts[["chr", "start", "end", score_column]]
    featurecount = feature_name + ".readCount"
    domain_counts.rename(columns={score_column: featurecount}, inplace=True)
    domain_counts["chr"] = domain_counts["chr"].astype("str")

    df = df.merge(domain_counts.drop_duplicates())
    # df = smart_merge(df, domain_counts.drop_duplicates())
    assert df.shape[0] == orig_shape, "Dimension mismatch"

    df[feature_name + ".RPM"] = 1e6 * df[featurecount] / float(total_counts)

    if not skip_rpkm_quantile:
        df[featurecount + ".quantile"] = df[featurecount].rank() / float(len(df))
        df[feature_name + ".RPM.quantile"] = df[feature_name + ".RPM"].rank() / float(
            len(df)
        )
        df[feature_name + ".RPKM"] = (
            1e3 * df[feature_name + ".RPM"] / (df.end - df.start).astype(float)
        )
        df[feature_name + ".RPKM.quantile"] = df[feature_name + ".RPKM"].rank() / float(
            len(df)
        )

    return df[~df.duplicated()]


def average_features(df, feature, feature_bam_list, skip_rpkm_quantile):
    feature_RPM_cols = [
        feature + "." + os.path.basename(feature_bam) + ".RPM"
        for feature_bam in feature_bam_list
    ]

    df[feature + ".RPM"] = df[feature_RPM_cols].mean(axis=1)

    if not skip_rpkm_quantile:
        feature_RPKM_cols = [
            feature + "." + os.path.basename(feature_bam) + ".RPKM"
            for feature_bam in feature_bam_list
        ]
        df[feature + ".RPM.quantile"] = df[feature + ".RPM"].rank() / float(len(df))
        df[feature + ".RPKM"] = df[feature_RPKM_cols].mean(axis=1)
        df[feature + ".RPKM.quantile"] = df[feature + ".RPKM"].rank() / float(len(df))

    return df


# From /seq/lincRNA/Jesse/bin/scripts/JuicerUtilities.R
#
bed_extra_colnames = [
    "name",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "itemRgb",
    "blockCount",
    "blockSizes",
    "blockStarts",
]


# JN: 9/13/19: Don't assume chromosomes start with 'chr'
# chromosomes = ['chr' + str(entry) for entry in list(range(1,23)) + ['M','X','Y']]   # should pass this in as an input file to specify chromosome order
def read_bed(
    filename,
    extra_colnames=bed_extra_colnames,
    chr=None,
    sort=False,
    skip_chr_sorting=True,
):
    skip = 1 if ("track" in open(filename, "r").readline()) else 0
    names = BED3_COLS + extra_colnames
    result = pd.read_table(
        filename, names=names, header=None, skiprows=skip, comment="#"
    )
    result = result.dropna(axis=1, how="all")  # drop empty columns
    assert result.columns[0] == "chr"

    result["chr"] = pd.Categorical(result["chr"], ordered=True)
    if chr is not None:
        result = result[result.chr == chr]
    if not skip_chr_sorting:
        result.sort_values("chr", inplace=True)
    if sort:
        result.sort_values(BED3_COLS, inplace=True)
    return result


def read_bedgraph(filename):
    read_bed(filename, extra_colnames=["score"], skip_chr_sorting=True)


def count_bam_mapped(bam_file):
    # Counts number of reads in a BAM file WITHOUT iterating.  Requires that the BAM is indexed
    # chromosomes = ['chr' + str(x) for x in range(1,23)] + ['chrX'] + ['chrY']
    command = "samtools idxstats " + bam_file
    data = check_output(command, shell=True)
    lines = data.decode("ascii").split("\n")
    # vals = list(int(l.split("\t")[2]) for l in lines[:-1] if l.split("\t")[0] in chromosomes)
    vals = list(int(l.split("\t")[2]) for l in lines[:-1])
    if not sum(vals) > 0:
        raise ValueError("Error counting BAM file: count <= 0")
    return sum(vals)


def count_tagalign_total(tagalign):
    result = int(
        check_output(
            "zcat {} | grep -E 'chr[1-9]|chr1[0-9]|chr2[0-2]|chrX|chrY' | wc -l".format(
                tagalign
            ),
            shell=True,
        )
    )
    assert result > 0
    return result


def count_bigwig_total(bw_file):
    bw = open_bigwig(bw_file)
    result = sum(
        l * bw.stats(ch, 0, l, "mean", exact=True)[0] for ch, l in bw.chroms().items()
    )
    assert (
        abs(result) > 0
    )  ## BigWig could have negative values, e.g. the negative-strand GroCAP bigwigs
    return result


def count_total(infile):
    filename = os.path.basename(infile)
    if "tagAlign" in filename:
        total_counts = count_tagalign_total(infile)
    elif filename.endswith(".bam"):
        total_counts = count_bam_mapped(infile)
    elif isBigWigFile(filename):
        total_counts = count_bigwig_total(infile)
    else:
        raise RuntimeError("Did not recognize file format of: " + infile)

    return total_counts


def parse_params_file(args):
    # Parse parameters file and return params dictionary
    params = {}

    params["default_accessibility_feature"] = determine_accessibility_feature(args)
    params["features"] = get_features(args)

    if args.expression_table:
        params["expression_table"] = args.expression_table.split(",")
    else:
        params["expression_table"] = ""

    return params


def get_features(args):
    features = {}

    if args.H3K27ac:
        features["H3K27ac"] = args.H3K27ac.split(",")

    if args.ATAC:
        features["ATAC"] = args.ATAC.split(",")

    if args.DHS:
        features["DHS"] = args.DHS.split(",")

    if args.supplementary_features is not None:
        supp = pd.read_csv(args.supplementary_features, sep="\t")
        for idx, row in supp.iterrows():
            features[row["feature_name"]] = row["file"].split(",")

    return features


def determine_accessibility_feature(args):
    if args.default_accessibility_feature is not None:
        return args.default_accessibility_feature
    elif (not args.ATAC) and (not args.DHS):
        raise RuntimeError(
            "Both DHS and ATAC have been provided. Must set one file to be the default accessibility feature!"
        )
    elif args.ATAC:
        return "ATAC"
    elif args.DHS:
        return "DHS"
    else:
        raise RuntimeError("At least one of ATAC or DHS must be provided!")


def compute_activity(df, access_col):
    if access_col == "DHS":
        if "H3K27ac.RPM" in df.columns:
            df["activity_base"] = np.sqrt(
                df["normalized_h3k27ac"] * df["normalized_dhs"]
            )
            df["activity_base_no_qnorm"] = np.sqrt(df["H3K27ac.RPM"] * df["DHS.RPM"])
        else:
            df["activity_base"] = df["normalized_dhs"]
            df["activity_base_no_qnorm"] = df["DHS.RPM"]
    elif access_col == "ATAC":
        if "H3K27ac.RPM" in df.columns:
            df["activity_base"] = np.sqrt(
                df["normalized_h3k27ac"] * df["normalized_atac"]
            )
            df["activity_base_no_qnorm"] = np.sqrt(df["H3K27ac.RPM"] * df["ATAC.RPM"])
        else:
            df["activity_base"] = df["normalized_atac"]
            df["activity_base_no_qnorm"] = df["ATAC.RPM"]
    else:
        raise RuntimeError("At least one of ATAC or DHS must be provided!")

    return df


def run_qnorm(df, qnorm, qnorm_method="rank", separate_promoters=True):
    # Quantile normalize epigenetic data to a reference
    #
    # Option to qnorm promoters and nonpromoters separately

    if qnorm is None:
        if "H3K27ac.RPM" in df.columns:
            df["normalized_h3k27ac"] = df["H3K27ac.RPM"]
        if "DHS.RPM" in df.columns:
            df["normalized_dhs"] = df["DHS.RPM"]
        if "ATAC.RPM" in df.columns:
            df["normalized_atac"] = df["ATAC.RPM"]
    else:
        qnorm = pd.read_csv(qnorm, sep="\t")
        nRegions = df.shape[0]
        col_dict = {
            "DHS.RPM": "normalized_dhs",
            "ATAC.RPM": "normalized_atac",
            "H3K27ac.RPM": "normalized_h3k27ac",
        }

        # & operator doesn't work in newer versions https://pastebin.com/8d2Ra9L1
        for col in set(df.columns.intersection(col_dict.keys())):
            # if there is no ATAC.RPM in the qnorm file, but there is ATAC.RPM in enhancers, then qnorm ATAC to DHS
            if col == "ATAC.RPM" and "ATAC.RPM" not in qnorm.columns:
                qnorm["ATAC.RPM"] = qnorm["DHS.RPM"]

            if not separate_promoters:
                qnorm = qnorm.loc[qnorm["enh_class" == "any"]]
                if qnorm_method == "rank":
                    interpfunc = interpolate.interp1d(
                        qnorm["rank"],
                        qnorm[col],
                        kind="linear",
                        fill_value="extrapolate",
                    )
                    df[col_dict[col]] = interpfunc(
                        (1 - df[col + ".quantile"]) * nRegions
                    ).clip(0)
                elif qnorm_method == "quantile":
                    interpfunc = interpolate.interp1d(
                        qnorm["quantile"],
                        qnorm[col],
                        kind="linear",
                        fill_value="extrapolate",
                    )
                    df[col_dict[col]] = interpfunc(df[col + ".quantile"]).clip(0)
            else:
                for enh_class in ["promoter", "nonpromoter"]:
                    this_qnorm = qnorm.loc[qnorm["enh_class"] == enh_class]

                    # Need to recompute quantiles within each class
                    if enh_class == "promoter":
                        this_idx = df.index[
                            np.logical_or(
                                df["class"] == "tss", df["class"] == "promoter"
                            )
                        ]
                    else:
                        this_idx = df.index[
                            np.logical_and(
                                df["class"] != "tss", df["class"] != "promoter"
                            )
                        ]
                    df.loc[this_idx, col + enh_class + ".quantile"] = df.loc[
                        this_idx, col
                    ].rank() / len(this_idx)

                    if qnorm_method == "rank":
                        interpfunc = interpolate.interp1d(
                            this_qnorm["rank"],
                            this_qnorm[col],
                            kind="linear",
                            fill_value="extrapolate",
                        )
                        df.loc[this_idx, col_dict[col]] = interpfunc(
                            (1 - df.loc[this_idx, col + enh_class + ".quantile"])
                            * len(this_idx)
                        ).clip(0)
                    elif qnorm_method == "quantile":
                        interpfunc = interpolate.interp1d(
                            this_qnorm["quantile"],
                            this_qnorm[col],
                            kind="linear",
                            fill_value="extrapolate",
                        )
                        df.loc[this_idx, col_dict[col]] = interpfunc(
                            df.loc[this_idx, col + enh_class + ".quantile"]
                        ).clip(0)

    return df
