# Running ABC pipeline as a WDL

## Initial Setup
1. Install Cromwell and WOMTools locally
   - Download the cromwell and womtools .jar files from [here](https://github.com/broadinstitute/cromwell/releases)
   - For ease of use, put them in the same place somewhere on your machine. In this example, we'll be using `~/jars/`
2. Create your `input.json` based on the one below. If you're not running this in terra, you may want to create a script to automate the creation of this input file.
   ```
   {
        "ABCpipeline.HiCdirTar": "/Users/myessail/broad/ABC-Enhancer-Gene-Prediction/example/input_data/HiC/bedgraph.tar.gz",
        "ABCpipeline.makeCandidateRegions.pval_cutoff": 0.1,
        "ABCpipeline.genes_bed": "/Users/myessail/broad/ABC-Enhancer-Gene-Prediction/example/config/RefSeqCurated.170308.bed.CollapsedGeneBounds.chr22.bed",
        "ABCpipeline.regions_whitelist": "/Users/myessail/broad/ABC-Enhancer-Gene-Prediction/example/config/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS.500bp.chr22.bed",
        "ABCpipeline.h3k27ac_bam": "/Users/myessail/broad/ABC-Enhancer-Gene-Prediction/example/input_data/Chromatin/ENCFF384ZZM.chr22.bam",
        "ABCpipeline.h3k27ac_bam_index": "/Users/myessail/broad/ABC-Enhancer-Gene-Prediction/example/input_data/Chromatin/ENCFF384ZZM.chr22.bam.bai",
        "ABCpipeline.regions_blacklist": "/Users/myessail/broad/ABC-Enhancer-Gene-Prediction/example/config/wgEncodeHg19ConsensusSignalArtifactRegions.bed",
        "ABCpipeline.makeCandidateRegions.nStrongestPeaks": 3000,
        "ABCpipeline.ubiq_genes": "/Users/myessail/broad/ABC-Enhancer-Gene-Prediction/example/config/UbiquitouslyExpressedGenesHG19.txt",
        "ABCpipeline.expression_table": "/Users/myessail/broad/ABC-Enhancer-Gene-Prediction/example/input_data/Expression/K562.ENCFF934YBO.TPM.txt",
        "ABCpipeline.cellType": "K562",
        "ABCpipeline.dnaseqbam": "/Users/myessail/broad/ABC-Enhancer-Gene-Prediction/example/input_data/Chromatin/wgEncodeUwDnaseK562AlnRep1.chr22.bam",
        "ABCpipeline.dnaseqbam_index": "/Users/myessail/broad/ABC-Enhancer-Gene-Prediction/example/input_data/Chromatin/wgEncodeUwDnaseK562AlnRep1.chr22.bam.bai",
        "ABCpipeline.makeCandidateRegions.peakExtendFromSummit": 250,
        "ABCpipeline.makePrediction.threshold": 0.022,
        "ABCpipeline.chrom_sizes": "/Users/myessail/broad/ABC-Enhancer-Gene-Prediction/example/config/chr22"
    }
    ```
3. Create the `bedgraph.tar.gz` file used in `ABCpipeline.HiCdirTar` by tarring the `/example/input_data/HiC/bedgraph` directory
4. Download and install docker, we use [Docker Desktop](https://www.docker.com/products/docker-desktop)

## Running the Pipeline
1. Start cromwell by running `java -jar ~/jars/cromwell.jar server` from the command line
2. Run the script by running `java -jar ~/jars/cromwell.jar run /path/to/abc-pipeline.wdl --inputs /path/toinputs.json`