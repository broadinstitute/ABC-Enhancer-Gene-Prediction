# Running ABC pipeline as a WDL

## Initial Setup
1. Install Cromwell and WOMTools locally
   - Download the cromwell and womtools .jar files from [here](https://github.com/broadinstitute/cromwell/releases)
   - For ease of use, put them in the same place somewhere on your machine. In this example, we'll be using `~/jars/`
2. Create your `input.json` based on the one below. If you're not running this in terra, you may want to create a script to automate the creation of this input file.
   ```
   {
     ## Defining candidate elements from a DHS or ATAC bam (makeCandidateRegions.py)
     "ABCpipeline.dnaseqbam": "/path/to/ABC-Enhancer-Gene-Prediction/example/input_data/Chromatin/wgEncodeUwDnaseK562AlnRep1.chr22.bam",
     "ABCpipeline.dnaseqbam_index": "/path/to/ABC-Enhancer-Gene-Prediction/example/input_data/Chromatin/wgEncodeUwDnaseK562AlnRep1.chr22.bam.bai",
     ## outDir
     "ABCpipeline.chrom_sizes": "/path/to/ABC-Enhancer-Gene-Prediction/example/config/chr22",
     "ABCpipeline.regions_blacklist": "/path/to/ABC-Enhancer-Gene-Prediction/example/config/wgEncodeHg19ConsensusSignalArtifactRegions.bed",
     "ABCpipeline.regions_whitelist": "/path/to/ABC-Enhancer-Gene-Prediction/example/config/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS.500bp.chr22.bed",
     "ABCpipeline.makeCandidateRegions.pval_cutoff": 0.1,
     "ABCpipeline.makeCandidateRegions.peakExtendFromSummit": 250,
     "ABCpipeline.makeCandidateRegions.nStrongestPeaks": 3000,

     ## Quantifying Enhancer Activity (run.neighborhoods.py)
     ## candidate_enhancer_regions is output from makeCandidateRegions.py
     "ABCpipeline.genes_bed": "/path/to/ABC-Enhancer-Gene-Prediction/example/config/RefSeqCurated.170308.bed.CollapsedGeneBounds.chr22.bed",
     "ABCpipeline.h3k27ac_bam": "/path/to/ABC-Enhancer-Gene-Prediction/example/input_data/Chromatin/ENCFF384ZZM.chr22.bam",
     "ABCpipeline.h3k27ac_bam_index": "/path/to/ABC-Enhancer-Gene-Prediction/example/input_data/Chromatin/ENCFF384ZZM.chr22.bam.bai",
     ## DHS -- set above
     "ABCpipeline.expression_table": "/path/to/ABC-Enhancer-Gene-Prediction/example/input_data/Expression/K562.ENCFF934YBO.TPM.txt",
     ## chrom_sizes -- set above
     "ABCpipeline.ubiq_genes": "/path/to/ABC-Enhancer-Gene-Prediction/example/config/UbiquitouslyExpressedGenesHG19.txt",
     "ABCpipeline.cellType": "K562",
     ## outDir


     ## Computing the ABC Score (predict.py)
     ## enhancers is output from run.neighborhoods.py
     ## genes is output from run.neighborhoods.py
     "ABCpipeline.HiCdirTar": "/path/to/ABC-Enhancer-Gene-Prediction/example/input_data/HiC/bedgraph.tar.gz", # this is a tarfile of the bedgraph dir
     ## scale_hic_using_powerlaw is a boolean
     "ABCpipeline.makePrediction.threshold": 0.022,
     ## cellType -- set above
     ## outDir
   }

    ```

For flags like "is_paired_end", for a true value, it must be set in the inputs.json file
3. Create the `bedgraph.tar.gz` file used in `ABCpipeline.HiCdirTar` by tarring the `/example/input_data/HiC/bedgraph` directory
4. Download and install docker, we use [Docker Desktop](https://www.docker.com/products/docker-desktop)

## Running the Pipeline
1. Start cromwell by running `java -jar ~/jars/cromwell.jar server` from the command line
2. Run the script by running `java -jar ~/jars/cromwell.jar run /path/to/abc-pipeline.wdl --inputs /path/to/inputs.json`
