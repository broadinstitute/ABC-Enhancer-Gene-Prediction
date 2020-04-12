#! bin/bash 
#################################################################################
## Get list of files from ENCODE Portal for ABC predictions
outdir=$2

## Download metadata files from ENCODE portal
if [ $1 == "hg19" ] | [$1 == "all"]:
	## splitting into pairedends + singleends download 
	wget --quiet -O pairedend_bam_hg19_DHS.tsv "https://www.encodeproject.org/metadata/?type=Experiment&files.assembly=hg19&files.output_type=alignments&files.file_type=bam&files.file_format=bam&files.assay_term_name=DNase-seq&files.run_type=paired-ended" -P $outdir
	wget --quiet -O pairedend_bam_hg19_ATAC.tsv "https://www.encodeproject.org/metadata/?type=Experiment&files.assembly=hg19&files.output_type=alignments&files.file_type=bam&files.file_format=bam&files.assay_term_name=ATAC-seq&files.run_type=paired-ended" -P $outdir
	wget --quiet -O pairedend_bam_hg19_H3K27ac.tsv "https://www.encodeproject.org/metadata/?type=Experiment&files.assembly=hg19&files.output_type=alignments&files.file_type=bam&files.file_format=bam&files.assay_term_name=ChIP-seq&target.label=H3K27ac&files.run_type=paired-ended" -P $outdir

	wget --quiet -O singleend_bam_hg19_DHS.tsv "https://www.encodeproject.org/metadata/?type=Experiment&files.assembly=hg19&files.output_type=alignments&files.file_type=bam&files.file_format=bam&files.assay_term_name=DNase-seq&files.run_type=single-ended" -P $outdir
	wget --quiet -O singleend_bam_hg19_ATAC.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=hg19&assay_title=ATAC-seq&files.run_type=single-ended" -P $outdir
	wget --quiet -O singleend_bam_hg19_H3K27ac.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=hg19&assay_title=Histone+ChIP-seq&target.label=H3K27ac&files.run_type=single-ended" -P $outdir

if [$1 == "GRCh38"] | [$1 == "all"]:
	# hg38 download - split pairedends, singleends
	wget --quiet -O pairedend_bam_GRCh38_DHS.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=GRCh38&assay_title=DNase-seq&files.run_type=paired-ended" -P $outdir
	wget --quiet -O pairedend_bam_GRCh38_ATAC.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=GRCh38&assay_title=ATAC-seq&files.run_type=paired-ended" -P $outdir
	wget --quiet -O pairedend_bam_GRCh38_H3K27ac.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=GRCh38&assay_title=Histone+ChIP-seq&target.label=H3K27ac&files.run_type=paired-ended" -P $outdir

	wget --quiet -O singleend_bam_GRCh38_DHS.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=GRCh38&assay_title=DNase-seq&files.run_type=single-ended" -P $outdir
	wget --quiet -O singleend_bam_GRCh38_ATAC.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=GRCh38&assay_title=ATAC-seq&files.run_type=single-ended" -P $outdir
	wget --quiet -O singleend_bam_GRCh38_H3K27ac.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=GRCh38&assay_title=Histone+ChIP-seq&target.label=H3K27ac&files.run_type=single-ended" -P $outdir

## To do: 
## Get cell types, handling potential treatments, genetic modifications, etc.  Does the collapsed ENCODE portal types help with this?

## Write script to pull more info on each BAM file, and filter based on the following
	## Check read lengths
	## For DNase data — mostly good
	## For Histone data – might often find files from multiple labs (default to Broad)


##################################################################################
## Output the BAM files into the parameters file as input into ABC script



#################################################################################
## Run ABC neighborhoods + predictions



