Reference files to run ABC for the mm10 genome assembly, built using gene annotations from GENCODE vM25.

*Note:* These files have not been thoroughly tested, please provide feedback if you encounter major issues!

- `mm10.chrom.sizes`: mm10 chromosome sizes 
- `mm10.blacklist.ENCFF547MET.bed`: regions of mm10 to exclude from analysis
- `gencode.vM25.protein_coding.TSS500bp.bed`: TSS locations (1 TSS chosen per gene)
- `gencode.vM25.protein_coding.genes.bed`: corresponding gene bounds
- `gencode.vM25.protein_coding.gene_promoter_class.tsv` [for use in downstream E2G models]: classification of genes as ubiquitously-expressed across cell types, derived from mapping annotations from human orthologues 

Example of ABC config reference section using these files
```
### REFERENCE FILES
ref:
        chrom_sizes: "reference/mm10_vM25/mm10.chrom.sizes"
        regions_blocklist: "reference/mm10_vM25/mm10.blacklist.ENCFF547MET.bed"
        ubiquitous_genes: "reference/UbiquitouslyExpressedGenes.txt" # not provided, does not affect predictions
        genes: "reference/mm10_vM25/gencode.vM25.protein_coding.genes.bed"
        genome_tss: "reference/mm10_vM25/gencode.vM25.protein_coding.TSS500bp.bed"
        qnorm: "reference/EnhancersQNormRef.K562.txt" # does not change
        abc_thresholds: "reference/abc_thresholds.tsv" # does not change
```

