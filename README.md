# Activity by Contact Model of Enhancer-Gene Specificity

The Activity-by-Contact (ABC) model predicts which enhancers regulate which genes on a cell type specific basis. This repository includes the code needed to run the ABC model as well as small sample data files, example commands, and some general tips and suggestions. We provide a brief description of the model below, see Fulco et al (BioArxiv 2019) for a full description.

## Requirements
For each cell-type, the inputs to the ABC model are:

 * Required Inputs
 	* bam file for Dnase-Seq or ATAC-Seq (indexed and sorted)
 	* bam file for H3K27ac ChIP-Seq (indexed and sorted)
 	* bed file containing candidate enhancer regions
 * Optional Inputs
 	* Hi-C data (see the Hi-C section below)
 	* A measure of gene expression (see gene expression section)

In addition the following (non-cell-type specific) genome annotation files are required

 * bed file containing gene annotations (may change across cell types if using cell-type specific TSS's)
 * bed file containing chromosome annotations

### Dependencies

The codebase relies on the following dependancies (tested version provided in 
parentheses):

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
matplotlib - Partial dependancy
scipy - Partial dependancy
```

## Description of the ABC Model

The Activity by Contact (ABC) model is designed to represent a mechanistic model in which enhancers activate gene expression upon enhancer-promoter contact. In a simple conception of such a model, the quantitative effect of an enhancer depends on the frequency with which it contacts a promoter multiplied by the strength of the enhancer (i.e., the ability of the enhancer to activate transcription upon contacting a promoter). Moreover, the relative contribution of an element on a gene’s expression (for example as assayed by the proportional decrease in expression upon CRISPR-inhibition) should depend on the element’s effect divided by the total effect of all elements.

To convert this conceptual framework into a practical score (which can be applied genome-wide), we formulated the ABC score:

ABC score for effect of element E on gene G = Activity of E × Contact frequency between E and G /  Sum of (Activity × Contact Frequency) over all candidate elements within 5 Mb.

Operationally, Activity (A) is defined as the geometric mean of the read counts of DNase-seq and H3K27ac ChIP-seq at an element E, and Contact (C) as the KR normalized Hi-C contact frequency between E and the promoter of gene G. Elements are defined as ~500bp regions centered on DHS peaks. 

 
## Running the ABC Model
Running the ABC model consists of the following steps:

 1. Define candidate enhancer regions, collect gene annotations and epigenetic data
 2. Quantifying enhancer activity
 3. Computing the ABC Score

### Step 1. Define candidate regions, gather gene annotations and Hi-C data

A typical way to define candidate elements is by calling peaks on a DNase-Seq or ATAC-Seq bam file (see Defining candidate elements from a DNase or ATAC bam). 

Gene annotations should be provided in .bed format.

We recommend either using cell-type specific Hi-C data or average Hi-C data (SEE PROVIDED). If neither of these are available, then the power-law relationship can be used as a proxy for Contact.

### Step 2. Quantifying Enhancer Activity: 

```run.neighborhoods.py``` will count DNase-seq (or ATAC-seq) and H3K27ac ChIP-seq reads in candidate enhancer regions. It also makes GeneList.txt, which includes data about genes and their promoters.

Sample Command:

```
python src/run.neighborhoods.py \
--candidate_enhancer_regions example/input_data/Chromatin/wgEncodeUwDnaseK562.mergedPeaks.slop175.withTSS500bp.chr22.bed \
--H3K27ac example/input_data/Chromatin/wgEncodeBroadHistoneK562H3K27ac_ENCFF000BWZ.q30.chr22.bam \
--DHS example/input_data/Chromatin/wgEncodeUwDnaseK562AlnRep1.chr22.bam \
--expression_table example/input_data/Expression/K562.ENCFF934YBO.TPM.txt \
--genes example/config/RefSeqCurated.170308.bed.CollapsedGeneBounds.chr22.bed \
--chrom_sizes example/config/chr22 \
--ubiquitously_expressed_genes example/config/UbiquitouslyExpressedGenesHG19.txt \
--cellType K562 \
--outdir example/ABC_output/Neighborhoods/ 
```

  * EnhancerList.txt: Candidate enhancer regions with Dnase-seq and H3K27ac ChIP-seq read counts
  * GeneList.txt: Dnase-seq and H3K27ac ChIP-seq read counts on gene bodies and gene promoter regions


### Step 3. Computing the ABC Score

Compute ABC scores by combining Activity as calculated by ```run.neighborhoods.py``` and Hi-C.

Sample Command:

```
python src/predict.py \
--enhancers example/ABC_output/Neighborhoods/EnhancerList.txt \
--genes example/ABC_output/Neighborhoods/GeneList.txt \
--HiCdir example/input_data/HiC/bedgraph/ \
--threshold .02 \
--cellType K562 \
--outdir example/ABC_output/Predictions/ 
```

The main output files are:

* EnhancerPredictions.txt: all element-gene pairs with scores above the provided threshold. Only includes expressed genes
* EnhancerPredictions.bedpe: Same as above in .bedpe format. Can be visualized in IGV.
* EnhancerPredictionsAllPutative.txt.gz: ABC scores for all element-gene pairs. Includes non-expressed genes and pairs with scores below the threshold. (use ```--make_all_putative``` to generate this file)

The default threshold of 0.02 corresponds to 70% recall and 63% precision in the Fulco et al 2019 dataset.
Columns are further defined in https://docs.google.com/spreadsheets/d/1UfoVXoCxUpMNPfGypvIum1-RvS07928grsieiaPX67I/edit?usp=sharing

## Defining Candidate Enhancers
'Candidate elements' are the set of putative enhancers for which ABC scores will be computed. In computing the ABC score, the product of DNase-seq (or ATAC-seq) and H3K27ac ChIP-seq reads will be counted in the candidate element. Thus the candidate elements should be regions of open (nucleasome depleted) chromatin of sufficient length to capture H3K27ac marks on flanking nucleosomes. In Fulco et al 2019, we defined candidate regions to be 500 bp (150bp of the DHS peak extended 175bp in each direction). 

### Defining candidate elements from a DHS or ATAC bam
A typical way to define candidate elements is by calling peaks from a DNase-seq or ATAC-seq bam file. Below we provide a convenience function for defining candidate regions using the MACS2 peak caller. 

```makeCandidateRegions.py``` is a wrapper around MACS2 which produces candidate regions from a Dnase-seq or ATAC-seq bam file. The script performs the following steps:

 1. Call peaks using MACS2
 2. Resize each peak to be a fixed number of base pairs centered on the peak summit
 3. Count DNase-seq reads in each peak and retain the top N peaks with the most read counts
 4. Remove any blacklisted regions and include any whitelisted regions

Sample command:

```
python src/makeCandidateRegions.py \
--bam example/input_data/Chromatin/wgEncodeUwDnaseK562AlnRep1.chr22.bam \
--outDir example/ABC_output/Peaks/ \
--chrom_sizes example/config/chr22 \
--regions_blacklist example/config/wgEncodeHg19ConsensusSignalArtifactRegions.bed \
--regions_whitelist example/config/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS.500bp.chr22.bed \
--pval_cutoff .1 \
--peakExtendFromSummit 250 \
--nStrongestPeaks 3000 
```

We recommend using ```--nStrongestPeaks 150000``` when making genome-wide peak calls. ```3000``` is just used for the small test run on chr22. 

Given that the ABC score uses absolute counts of Dnase-seq reads in each region, ```curateFeatures.py``` attempts to select the strongest peaks as measured by absolute read counts (not read counts relative to some background rate). In order to do this, we first call peaks using a lenient significance threshold (.1 in the above example) and then count reads in each of called peaks. 

We recommend removing elements overlapping regions of the genome that have been observed to accumulate anomalous number of reads in epigenetic sequencing experiments (‘blacklisted regions’). For convenience, we provide the list of blackedlisted regions available from https://sites.google.com/site/anshulkundaje/projects/blacklists.

## Contact and Hi-C
Given that cell-type specific Hi-C data is more difficult to generate than ATAC-seq or ChIP-seq, we have explored alternatives to using cell-type specific Hi-C data. It is known that Hi-C contact frequencies generally follow a powerlaw relationship (with respect to genomic distance) and that many TADs, loops and other structural features of the 3D genome are **not** cell-type specific. 

We have found that, for most genes, using an average Hi-C profile in the ABC model gives approximately equally good performance as using a cell-type specific Hi-C profile. To facilitate making ABC predictions in a large panel of cell types, including those without cell type-specific Hi-C data, we have provided an average Hi-C profile (averaged across 10 cell lines) in this repository. 

In the case where cell-type specific Hi-C data is available, we provide a pipeline which takes as input a .hic file, and formats it as the ABC model code expects (see below)

### Description of Average Hi-C data provided
* Generate bedgraphs for ten cell types using pipeline described below
* Powerlaw scale each cell type's bedgraphs to K562
* For each gene, generate an average bedgraph profile by averaging together the bedgraphs from all ten cell types


### Pipeline to Download and Format Hi-C data

The below pipeline will download a Hi-C matrix from Juicebox (in .hic format) and generate tss-anchored bedgraphs (the format the ABC model expects).

Three steps

1. Download raw data using Juicebox
2. Make HiC Bedgraphs
3. Get powerlaw parameters (Optional)

```
#Download hic matrix file from juicebox
python src/juicebox_dump.py \
--hic_file https://hicfiles.s3.amazonaws.com/hiseq/k562/in-situ/combined.hic \
--juicebox "java -jar juicer_tools.jar" \
--outdir $HICDIR/raw/ \
--chromosomes 22
```

```
#Make a virtual 4C bedgraph anchored at the TSS of each gene
python src/make_bedgraph_from_HiC.py \
--outdir $HICDIR/bedgraph/ \
--genes example/config/RefSeqCurated.170308.bed.CollapsedGeneBounds.chr22.bed \
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

## Quantile Normalization

The ABC Score uses the quantitative signal of Hi-C, ATAC-Seq and H3K27ac ChIP-Seq. As such it is sensitive to differences in these datasets due to experimental or technical procedures. For example, ABC scores computed on an epigenetic dataset with low signal to noise will be lower than ABC scores computed on an epigenetic dataset with high signal to noise. 

In an effort to make ABC scores comparable across cell types, the ABC model code supports quantile normalizing the epigenetic signal in candidate enhancer regions to some reference. The reference we provide in this repository is computed on ~160,000 DHS peaks in K562. 

Empirically, we have found that applying quantile normalization makes ABC predictions more comparable across cell types. However, it may not be applicable to all circumstances.

We recommend using quantile normalization if you are looking to compare ABC scores across cell types. Additionally, the threshold value on the ABC score of .02 (described in Fulco et al) is calculated based on the K562 epigenetic data. 


## Tips and Comments

* Accurate transcription start site annotations are critical.
* We have found that ubiquitously expressed genes appear insensitive to the effects of distal enhancers. For completeness, this code calculates the ABC score for all genes and flags ubiquitously expressed genes.
* The size of candidate enhancer elements is important. For example, if two candidate regions are merged, then the ABC score of the merged region will be approximately the sum of the ABC scores for each individual region.

## Citation

Fulco CP, Nasser J, Jones TR, Munson G, Bergman D, Subramanian V, Grossman SR, Anyoha R, Patwardhan TA, Nguyen TH, Kane M, Doughty B, Perez E, Durand NC, Stamenova EK, Lieberman Aiden E, Lander ES, Engreitz JM. Activity-by-Contact model for enhancer specificity from thousands of CRISPR perturbations. bioRxiv. 2019 Jan 26.
