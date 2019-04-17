# Activity by Contact Model of Enhancer-Gene Specificity

The Activity-by-Contact (ABC) model predicts which enhancers regulate which genes on a cell type specific basis. The ABC model is motivated by the conceptual notion that the contribution of an enhancer to transcription of a gene is dependent on two quantities: the instrinsic strength of the enhancer ('Activity') and the Contact frequency between the enhancer and gene promoter. The ABC model uses DNase-seq and H3K27ac ChIP-seq as measures of enhancer Activity, and Hi-C as a measure enhancer-promoter Contact frequency (see Description of the ABC Model). The ABC model is also able to make accurate predictions in the absence of cell-type specific Hi-C data (see Hi-C section). Thus the minimal requirements needed in order to generate enhancer-gene predictions using the ABC model are a measure of chromatin accessibility (typically Dnase-seq or ATAC-seq) and a measure of enhancer strength (H3K27ac ChIP-seq) in the cell type of interest.

The ABC model is effectively able to predict the effects of enhancers on gene expression as measured by perturbational experiments. This repository implements the ABC model as described in Fulco et al (BioArxiv 2019). The code in this repository can be used to generate enhancer-gene predictions for any cell type with the required epigenetic data. 

We provide the code needed to run the ABC model as well as small sample data files, example commands, and some general tips and suggestions. 

For each cell-type, the inputs to the ABC model are:

 * Required Inputs
 	* bam file for Dnase-Seq or ATAC-Seq (indexed and sorted)
 	* bam file for H3K27ac ChIP-Seq (indexed and sorted)
 	* bed file containing candidate enhancer regions
 	* bed file containing gene annotations
 	* bed file containing chromosome annotations
 * Optional Inputs
 	* Hi-C data (see the Hi-C section below)
 	* A measure of gene expression (see gene expression section)

In addition the following (non-cell-type specific) genome annotation files are required

 * bed file containing gene annotations
 * bed file containing chromosome annotations

The ABC model produces output in the following directory structure

 * ABC_output
    * Peaks
		* Candidate enhancer regions and files related to MACS2 peak calls. Note this output directory is only applicable if candidate regions are defined using ```curateFeatures.py```
	* Neighborhoods
	  * EnhancerList.txt: Candidate enhancer regions with Dnase-seq and H3K27ac ChIP-seq read counts
	  * GeneList.txt: 
  * Predictions
     * EnhancerPredictions.txt: Enhancer-Gene predictions for highly expressed genes with ABC scores above the provided threshold. This is the main ABC output file. 
     * Predictions.bedpe: Enhancer-Gene predictions in bedpe format, which can be visualized in IGV
     * genes/: Directory containing a separate file for all genes. Each file in this directory contains ABC scores for each candidate enhancer 5mb of the gene. These files contain predicted negatives as well as predicted positives.

## Description of the ABC Model

We designed the Activity by Contact (ABC) score to represent a mechanistic model in which enhancers contact target promoters to activate gene expression. In a simple conception of such a model, the quantitative effect of an enhancer depends on the frequency at which it contacts a promoter multiplied by the strength of the enhancer (i.e., the ability of the enhancer to activate transcription upon contacting a promoter). Moreover, the relative contribution of an element on a gene’s expression (for example as assayed by the proportional decrease in expression upon CRISPR-inhibition) should depend on the element’s effect divided by the total effect of all elements.

To extend this conceptual framework to enable computing the quantitative effects of enhancers on the expression of any gene, we formulated the ABC score:

ABC score for effect of element E on gene G = Activity of E × Contact frequency between E and G /  Sum of (Activity × Contact Frequency) over all candidate elements within 5 Mb.

Operationally, Activity (A) is defined as the geometric mean of the read counts of DNase-seq and H3K27ac ChIP-seq at an element E, and Contact (C) as the KR normalized Hi-C contact frequency between E and the promoter of gene G, and elements are defined as ~500bp regions centered on DHS peaks. 

 
## Running the ABC Model
Running the ABC model consists of the following steps:

 1. Set up cell type configuration files, setup directories
 2. Quantifying the activity level of each candidate enhancer
 3. Making enhancer-gene predictions

The ABC model relies on the following dependancies (tested version provided in parentheses):

```
Python (3.4)
samtools (0.1.19)
bedtools (2.26.0)
Tabix (0.2.5) - Partial dependancy
MACS2 (2.1.1.20160309) - Partial dependancy
Java (1.7) - Partial dependancy

Python packages:
numpy
pandas
intervaltree
psyam
pyBigWig
matplotlib
```


### Step 1. Setting up configuration files and directories

Example configuration files are provided in example/config/. 

**cellTypeParameters.txt**: Add one entry per cell type to the format described in example/config/cellTypeParameters.txt. Replicate experiments should be inserted as comma-delimted entries. 

**genomes.txt**: Add one entry per genome to the format described in example/config/genomes.txt. 'name' corresponds to 'genome' column of cellTypeParameters.txt

Define and make directories

```
OUTDIR=example/ABC_output/
PEAKDIR=$OUTDIR/Peaks/
NBHDDIR=$OUTDIR/Neighborhoods/
PREDDIR=$OUTDIR/Predictions/

mkdir -p $PEAKDIR
mkdir -p $NBHDDIR
mkdir -p $PREDDIR

```

### Step 2. Quantifying Enhancer Activity: 
NOTE: This section assumes candidate enhancer elements have already been defined (See below section on defining candidate elements)

```run.neighborhoods.py``` will count DNase-seq (or ATAC-seq) and H3K27ac ChIP-seq reads in candidate enhancer regions. It also makes GeneList.txt.

Sample Command:

```
python src/run.neighborhoods.py \
--cellType K562 \
--params_file example/config/cellTypeParameters.txt \
--outdir $NBHDDIR \
--genome example/config/genomes.txt \
--candidate_enhancer_regions example/input_data/Chromatin/wgEncodeUwDnaseK562.mergedPeaks.chr22.slop175.bed \
--genes example/config/RefSeqCurated.170308.chr22.small.bed
```
### Step 3. Making predictions

Sample Command:

```
python src/predict.py \
--cellType K562 \
--params_file example/config/cellTypeParameters.txt \
--outdir $PREDDIR \
--HiCdir example/input_data/HiC/bedgraph/ \
--nbhd_directory $NBHDDIR \
--threshold .022
```

The main output file is $PREDDIR/EnhancerPredictions.txt.
It contains all element-gene connections with ABC.score greater than the threshold. The default threshold of 0.022 corresponds to 70% recall and 74% precision in the Fulco et al 2019 dataset.

## Defining Candidate Enhancers
'Candidate elements' are the set of putative enhancers for which ABC scores will be computed. In computing the ABC score, the sum of DNase-seq (or ATAC-seq) and H3K27ac ChIP-seq reads will be counted in the candidate element. Thus the candidate elements should be regions of open (nucleasome depleted) chromatin of sufficient length to capture H3K27ac marks on flanking nucleosomes. In Fulco et al 2019, we defined candidate regions to be 500 bp (150bp of the DHS peak extended 175bp in each direction). 

### Defining candidate elements from a DHS or ATAC bam
A typical way to define candidate elements is by calling peaks from a DNase-seq or ATAC-seq bam file. Below we provide a convenience function for defining candidate regions using the MACS2 peak caller. 

```curateFeatures.py``` is a wrapper around MACS2 which produces candidate regions from a Dnase-seq or ATAC-seq bam file. The script performs the following steps:

 1. Call peaks using MACS2
 2. Resize each peak to be a fixed number of base pairs centered on the peak summit
 3. Count DNase-seq reads in each peak and retain the top N peaks with the most read counts
 4. Remove any blacklisted regions and include any whitelisted regions

Sample command:

```
python src/curateFeatures.py \
--cellType K562 \
--outDir $PEAKDIR \
--params_file example/config/cellTypeParameters.txt \
--genome example/config/genomes.txt \
--regions_blacklist example/config/wgEncodeHg19ConsensusSignalArtifactRegions.bed \
--regions_whitelist example/config/RefSeqCurated.170308.chr22.TSS1KB.bed \
--pval_cutoff .1 \
--nStrongestPeaks 175000 \
--peakExtendFromSummit 250
```
Given that the ABC score uses absolute counts of Dnase-seq reads in each region, ```curateFeatures.py``` attempts to select the strongest peaks as measured by absolute read counts (not read counts relative to some background rate). In order to do this, we first call peaks using a lenient significance threshold (.1 in the above example) and then count reads in each of called peaks. 

We recommend removing elements overlapping regions of the genome that have been observed to accumulate anomalous number of reads in epigenetic sequencing experiments (‘blacklisted regions’). For convenience, we provide the list of blackedlisted regions available from https://sites.google.com/site/anshulkundaje/projects/blacklists.

Different peak calling algorithms will produce varying number of peaks of variable length. Empirically we have noticed that the number of peaks and their width depends on the signal to noise ratio of the DNase-seq dataset. We note that defining candidate elements is an ongoing area of research...

### Defining candidate elements from an ENCODE peak file
```
bedtools slop -b 175 -i example/input_data/Chromatin/wgEncodeUwDnaseK562.mergedPeaks.chr22.bed -g example/config/chr22.bed | bedtools merge -i stdin > example/input_data/Chromatin/wgEncodeUwDnaseK562.mergedPeaks.chr22.slop175.bed
```

## Contact and Hi-C
Given that cell-type specific Hi-C data is more difficult to generate than ATAC-seq or ChIP-seq, we have explored alternatives to using cell-type specific Hi-C data. It is known that Hi-C contact frequencies generally follow a powerlaw relationship (with respect to genomic distance) and that many TADs, loops and other structural features of the 3D genome are **not** cell-type specific. 

Using an average Hi-C profile in the ABC model gives approximately equally good performance as using a cell-type specific profile. To facilitate making ABC predictions in a large panel of cell types, including those without cell type-specific Hi-C data, we have provided an average Hi-C profile (averaged across 8 cell lines) in this repository. 

In the case where cell-type specific Hi-C data is available, we provide a pipeline which takes as input a .hic file, and formats it as the ABC model code expects (see below)

### Description of Average Hi-C data provided
* Generate bedgraphs for ten cell types using pipeline described below
* Powerlaw scale each cell type's bedgraphs to K562
* For each gene, generate an average bedgraph profile by averaging together the bedgraphs from all ten cell types


### Pipeline to Download and Format Hi-C data

When predicting enhancers for a specific gene, the ABC model requires the row of the hic matrix corresponding to the TSS of the gene (given as a begraph). The below pipeline will download a Hi-C matrix from Juicebox (in .hic format) and generate tss-anchored bedgraphs.

Three steps

1. Download raw data using Juicebox
2. Make HiC Bedgraphs
3. Get powerlaw parameters (Optional)

```
#Make directory structure
HICDIR=example/input_data/HiC/
mkdir -p $HICDIR/raw/
mkdir -p $HICDIR/bedgraph/
mkdir -p $HICDIR/powerlaw/
```

```
#Download hic matrix file from juicebox
python src/juicebox_dump.py \
--hic_file https://hicfiles.s3.amazonaws.com/hiseq/k562/in-situ/combined.hic \
--outdir $HICDIR/raw/ \
--chromosomes 22
```

```
#Make a virtual 4C bedgraph anchored at the TSS of each gene
python src/make_bedgraph_from_HiC.py \
--outdir $HICDIR/bedgraph/ \
--genes example/config/RefSeqCurated.170308.chr22.small.bed \
--hic_dir $HICDIR/raw/5kb_resolution_intrachromosomal/
```

```
#Fit HiC data to powerlaw model and extract parameters
python src/compute_powerlaw_fit_from_hic.py \
--bedDir $HICDIR/bedgraph/ \
--outDir $HICDIR/powerlaw/
```

## Gene Expression in ABC
The ABC model is designed to predict the effect of enhancers on expressed genes. If a gene is not expressed in a given cell type (or cell state) then we assume it does not have any activating enhancers (enhancers for which inhibition of the enhancer would lead to decrease in gene expression). Thus we typically only report enhancer-gene connections for expressed genes.

In the absence of expression data, DNase-seq and H3K27ac ChIP-seq at the gene promoter can be used as a proxy for expression. We suggest only considering enhancer-gene connections for genes with sufficiently active promoters (for instance in the top half of gene promoters in the cell type)

## Tips and Suggestions

* Accurate transcription start site annotations are critical.
* We have found that ubiquitously expressed genes appear insensitive to the effects of distal enhancers. For completeness, this code calculates the ABC score for all genes and flags ubiquitously expressed genes.

## Citation

Fulco CP, Nasser J, Jones TR, Munson G, Bergman D, Subramanian V, Grossman SR, Anyoha R, Patwardhan TA, Nguyen TH, Kane M, Doughty B, Perez E, Durand NC, Stamenova EK, Lieberman Aiden E, Lander ES, Engreitz JM. Activity-by-Contact model for enhancer specificity from thousands of CRISPR perturbations. bioRxiv. 2019 Jan 26.
