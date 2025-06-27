### hg19 reference files
- _Genes_ ([hg19/CollapsedGeneBounds.hg19.bed](hg19/CollapsedGeneBounds.hg19.bed)) and _promoters_ ([hg19/CollapsedGeneBounds.hg19.TSS500bp.bed](hg19/CollapsedGeneBounds.hg19.TSS500bp.bed))
  - These files are derived from a previous set of RefSeq annotations. We defined one promoter per gene symbol, defined as a 500 base pair region centered around the RefSeq TSS with the largest number of coding isoforms.
  - We curated set of genes by matching the gene symbols with Ensembl IDs using the HUGO database and manual annotation, then used the GENCODE v29 database to annotate each gene with its “gene type”.
  - We retained genes that 1) were included in the GENCODE v29 database AND 2) had a gene type of “protein_coding,” “processed_transcript,” or “lincRNA” for a total of 20,666 genes. 
- _Chromosome sizes:_ [hg19/chrom_sizes.tsv](hg19/chrom_sizes.tsv)
- _Exclusion list:_ [hg19/wgEncodeHg19ConsensusSignalArtifactRegions.bed](hg19/wgEncodeHg19ConsensusSignalArtifactRegions.bed)
  - Derived from the [ENCODE DAC Exclusion List Regions](https://www.encodeproject.org/annotations/ENCSR636HFF/)

### hg38 reference files
- _Genes and promoters (several options):_
   - [hg38/CollapsedGeneBounds.hg38.bed](hg38/CollapsedGeneBounds.hg38.bed)) and [hg38/CollapsedGeneBounds.hg38.TSS500bp.bed](hg38/CollapsedGeneBounds.hg38.TSS500bp.bed)
      - These are the hg19 reference files lifted over into hg38. The gene symbols, gene bounds, and promoter bounds were not otherwise changed.
   - [hg38/CollapsedGeneBound.hg38.GENCODEv43GeneSymbol.bed](hg38/CollapsedGeneBound.hg38.GENCODEv43GeneSymbol.bed)) and [hg38/CollapsedGeneBound.hg38.GENCODEv43GeneSymbol.TSS500bp.bed](hg38/CollapsedGeneBound.hg38.GENCODEv43GeneSymbol.TSS500bp.bed)
      - These files were created to be compatible with GENCODE v43. We filtered the above hg38 reference files to the 20,531 genes with exactly one corresponding gene entry from GENCODE v43.
      - We replaced the RefSeq gene symbols with the corresponding GENCODE v43 gene symbols.
      - The gene bounds and promoter bounds were not otherwise changed.
   - **NOTE**: The [scE2G pipeline, v1.2](https://github.com/EngreitzLab/scE2G/tree/v1.2) by default uses a third version of hg38 [gene](https://github.com/EngreitzLab/scE2G/blob/v1.2/resources/genome_annotations/CollapsedGeneBounds.hg38.intGENCODEv43.bed) and [TSS](https://github.com/EngreitzLab/scE2G/blob/v1.2/resources/genome_annotations/CollapsedGeneBounds.hg38.intGENCODEv43.TSS500bp.bed) annotations. These are equivalent  to the GENCODE v43 compatible files, **except** they retained the original RefSeq gene symbols. 
- _Chromosome sizes:_ [hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv](hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv)
- _Exculsion list:_ [hg38/GRCh38_unified_blacklist.bed](hg38/GRCh38_unified_blacklist.bed)
  - Obtained from the [ENCODE DAC Exclusion List Regions](https://www.encodeproject.org/annotations/ENCSR636HFF/)
 
### mm10 and mm39 reference files
Please reference [mm10_vM25/README.md](mm10_vM25/README.md) and [mm39_vM35/README.md](mm39_vM35/README.md) for details.
