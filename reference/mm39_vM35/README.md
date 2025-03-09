Reference files to run ABC for the mm39 genome assembly, built using gene annotations from GENCODE vM35.

*Note:* These files have not been thoroughly tested, please provide feedback if you encounter major issues!

- `mm39.chrom.sizes`: mm39 chromosome sizes 
- `mm39.blacklist.excluderanges.bed`: regions of mm39 to exclude from analysis, sourced from the [excluderanges GitHub](https://github.com/dozmorovlab/excluderanges?tab=readme-ov-file)
- `gencode.vM35.protein_coding.TSS500bp.bed`: TSS locations (1 TSS chosen per gene)
- `gencode.vM35.protein_coding.genes.bed`: corresponding gene bounds
- `gencode.vM35.protein_coding.gene_promoter_class.tsv` [for use in downstream E2G models]: classification of genes as ubiquitously-expressed across cell types, derived from mapping annotations from human orthologues 

Example of ABC config reference section using these files
```
### REFERENCE FILES
ref:
        chrom_sizes: "reference/mm39_vM35/mm10.chrom.sizes"
        regions_blocklist: "reference/mm39_vM35/mm39.blacklist.excluderanges.bed"
        ubiquitous_genes: "reference/UbiquitouslyExpressedGenes.txt" # not provided, does not affect predictions
        genes: "reference/mm39_vM35/gencode.vM35.protein_coding.genes.bed"
        genome_tss: "reference/mm39_vM35/gencode.vM35.protein_coding.TSS500bp.bed"
        qnorm: "reference/EnhancersQNormRef.K562.txt" # does not change
        abc_thresholds: "reference/abc_thresholds.tsv" # does not change
```
