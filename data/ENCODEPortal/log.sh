#################################################################################
## Get list of files from ENCODE Portal for ABC predictions

## Download metadata files from ENCODE portal
wget --quiet -O bam_hg19_DHS.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=hg19&assay_title=DNase-seq"
wget --quiet -O bam_hg19_ATAC.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=hg19&assay_title=ATAC-seq"
wget --quiet -O bam_hg19_H3K27ac.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=hg19&assay_title=Histone+ChIP-seq&target.label=H3K27ac"

wget --quiet -O bam_GRCh38_DHS.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=GRCh38&assay_title=DNase-seq"
wget --quiet -O bam_GRCh38_ATAC.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=GRCh38&assay_title=ATAC-seq"
wget --quiet -O bam_GRCh38_H3K27ac.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=GRCh38&assay_title=Histone+ChIP-seq&target.label=H3K27ac"

## To do: 
## Get cell types, handling potential treatments, genetic modifications, etc.  Does the collapsed ENCODE portal types help with this?

## Write script to pull more info on each BAM file, and filter based on the following
	## Check read lengths
	## For DNase data — mostly good
	## For Histone data – might often find files from multiple labs (default to Broad)


##################################################################################
# From Kristy:
# grab celltypes
cut -f7 bam_hg19_H3K27ac.tsv | sed 1d > celltypes_hg19_H3K27ac.txt
cut -f7 bam_hg19_DHS.tsv | sed 1d > celltypes_hg19_accessibility.txt
cut -f7 bam_hg19_ATAC.tsv | sed 1d >> celltypes_hg19_accessibility.txt

# remove whitespace
sed -i '1d' celltypes_H3K27ac.tsv
sed -i '1,2d' celltypes_DHS.tsv
# grab common celltypes:
comm -12 <(sort celltypes_H3K27ac.txt) <(sort celltypes_DHS.txt) | sort -u > common_celltypes.txt


##################################################################################
## Output the BAM files into the parameters file as input into ABC script



#################################################################################
## Run ABC neighborhoods + predictions



