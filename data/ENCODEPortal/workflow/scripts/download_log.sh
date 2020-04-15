#################################################################################
## Get list of files from ENCODE Portal for ABC predictions

## Download metadata files from ENCODE portal
wget --quiet -O bam_hg19_DHS.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=hg19&assay_title=DNase-seq" 
wget --quiet -O bam_hg19_ATAC.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=hg19&assay_title=ATAC-seq" 
wget --quiet -O bam_hg19_H3K27ac.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=hg19&assay_title=Histone+ChIP-seq&target.label=H3K27ac" 


wget --quiet -O bam_GRCh38_DHS.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=GRCh38&assay_title=DNase-seq" 
wget --quiet -O bam_GRCh38_ATAC.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=GRCh38&assay_title=ATAC-seq" 
wget --quiet -O bam_GRCh38_H3K27ac.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=GRCh38&assay_title=Histone+ChIP-seq&target.label=H3K27ac"

## Download metadata files from ENCODE portal
## seems like there are some entries that are considered BOTH paired-ended and single-ended 
## splitting into pairedends + singleends download 
wget --quiet -O pairedend_bam_hg19_DHS.tsv "https://www.encodeproject.org/metadata/?type=Experiment&files.assembly=hg19&files.output_type=alignments&files.file_type=bam&files.file_format=bam&files.assay_term_name=DNase-seq&files.run_type=paired-ended" 
wget --quiet -O pairedend_bam_hg19_ATAC.tsv "https://www.encodeproject.org/metadata/?type=Experiment&files.assembly=hg19&files.output_type=alignments&files.file_type=bam&files.file_format=bam&files.assay_term_name=ATAC-seq&files.run_type=paired-ended" 
wget --quiet -O pairedend_bam_hg19_H3K27ac.tsv "https://www.encodeproject.org/metadata/?type=Experiment&files.assembly=hg19&files.output_type=alignments&files.file_type=bam&files.file_format=bam&files.assay_term_name=ChIP-seq&target.label=H3K27ac&files.run_type=paired-ended"

wget --quiet -O singleend_bam_hg19_DHS.tsv "https://www.encodeproject.org/metadata/?type=Experiment&files.assembly=hg19&files.output_type=alignments&files.file_type=bam&files.file_format=bam&files.assay_term_name=DNase-seq&files.run_type=single-ended"
wget --quiet -O singleend_bam_hg19_ATAC.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=hg19&assay_title=ATAC-seq&files.run_type=single-ended"
wget --quiet -O singleend_bam_hg19_H3K27ac.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=hg19&assay_title=Histone+ChIP-seq&target.label=H3K27ac&files.run_type=single-ended"


# hg38 download - split pairedends, singleends
wget --quiet -O pairedend_bam_GRCh38_DHS.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=GRCh38&assay_title=DNase-seq&files.run_type=paired-ended"
wget --quiet -O pairedend_bam_GRCh38_ATAC.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=GRCh38&assay_title=ATAC-seq&files.run_type=paired-ended"
wget --quiet -O pairedend_bam_GRCh38_H3K27ac.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=GRCh38&assay_title=Histone+ChIP-seq&target.label=H3K27ac&files.run_type=paired-ended" 

wget --quiet -O singleend_bam_GRCh38_DHS.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=GRCh38&assay_title=DNase-seq&files.run_type=single-ended"
wget --quiet -O singleend_bam_GRCh38_ATAC.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=GRCh38&assay_title=ATAC-seq&files.run_type=single-ended"
wget --quiet -O singleend_bam_GRCh38_H3K27ac.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=GRCh38&assay_title=Histone+ChIP-seq&target.label=H3K27ac&files.run_type=single-ended"

