#################################################################################
## Get list of files from ENCODE Portal for ABC predictions
# review of fastq and bam 

## Download metadata files from ENCODE portal
wget --quiet -O bam_hg19_DHS.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=hg19&assay_title=DNase-seq&audit.ERROR.category%21=extremely+low+spot+score&aud%5Cit.ERROR.category%21=extremely+low+read+depth&audit.NOT_COMPLIANT.category%21=insufficient+read+depth%5C&files.status=released&files.assembly=hg19&files.file_type=fastq" 
wget --quiet -O bam_hg19_ATAC.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=hg19&assay_title=ATAC-seq&audit.ERROR.category%21=extremely+low+spot+score&aud%5Cit.ERROR.category%21=extremely+low+read+depth&audit.NOT_COMPLIANT.category%21=insufficient+read+depth%5C&files.status=released&files.assembly=hg19&files.file_type=fastq" 
wget --quiet -O bam_hg19_H3K27ac.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=hg19&assay_title=Histone+ChIP-seq&target.label=H3K27ac&audit.ERROR.category%21=extremely+low+spot+score&aud%5Cit.ERROR.category%21=extremely+low+read+depth&audit.NOT_COMPLIANT.category%21=insufficient+read+depth%5C&files.status=released&files.assembly=hg19&files.file_type=fastq" 


wget --quiet -O bam_GRCh38_DHS.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=GRCh38&assay_title=DNase-seq&audit.ERROR.category%21=extremely+low+spot+score&aud%5Cit.ERROR.category%21=extremely+low+read+depth&audit.NOT_COMPLIANT.category%21=insufficient+read+depth%5C&files.status=released&files.assembly=hg19&files.file_type=fastq" 
wget --quiet -O bam_GRCh38_ATAC.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=GRCh38&assay_title=ATAC-seq&audit.ERROR.category%21=extremely+low+spot+score&aud%5Cit.ERROR.category%21=extremely+low+read+depth&audit.NOT_COMPLIANT.category%21=insufficient+read+depth%5C&files.status=released&files.assembly=hg19&files.file_type=fastq" 
wget --quiet -O bam_GRCh38_H3K27ac.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=GRCh38&assay_title=Histone+ChIP-seq&target.label=H3K27ac&audit.ERROR.category%21=extremely+low+spot+score&aud%5Cit.ERROR.category%21=extremely+low+read+depth&audit.NOT_COMPLIANT.category%21=insufficient+read+depth%5C&files.status=released&files.assembly=hg19&files.file_type=fastq"


