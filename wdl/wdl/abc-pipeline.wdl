version 1.0

workflow ABCpipeline {
    meta {
        description: "WDL version of the ABC pipeline"
    }

    parameter_meta {
        dnaseqbam: "bam file to call peaks on"
        chrom_sizes: "File listing chromosome size annotations"
        is_paired_end: "Flag denoting whether bam file is paired end. Needed for MACS2. True if present. False if not"
        regions_blacklist: "Bed file of regions to forcibly exclude from candidate enhancers"
        regions_whitelist: "Bed file of regions to forcibly include in candidate enhancers. Overrides regions_blacklist"
        genes_bed: "bed file with gene annotations. Must be in bed-6 format. Will be used to assign TSS to genes."
        #genes_for_class_assignment: "bed gene annotations for assigning elements to promoter/genic/intergenic classes. Will not be used for TSS definition"
        ubiq_genes: "File listing ubiquitously expressed genes. These will be flagged by the model, but this annotation does not affect model predictions"
        gene_name_annotations: "Comma delimited string of names corresponding to the gene identifiers present in the name field of the gene annotation bed file"
        primary_gene_identifier: "Primary identifier used to identify genes. Must be present in gene_name_annotations. The primary identifier must be unique"
        # DHS: "Comma delimited string of DHS .bam files. Either ATAC or DHS must be provided" is this a duplicate?
        # ATAC: "Comma delimited string of ATAC .bam files. Either ATAC or DHS must be provided"
        # default_accessibility_feature: "If both ATAC and DHS are provided, this flag must be set to either 'DHS' or 'ATAC' signifying which datatype to use in computing activity"
        # will want conditional here ^
        expression_table: "Comma delimited string of gene expression files"
        qnorm:"Quantile normalization reference file"
        tss_slop_for_class_assignment: "Consider an element a promoter if it is within this many bp of a tss"
        skip_rpkm_quantile: "Do not compute RPKM and quantiles in EnhancerList.txt"
        use_secondary_counting_method: "Use a slightly slower way to count bam over bed. Also requires more memory. But is more stable"
        chrom_sizes: "Genome file listing chromosome sizes. Also requires associated .bed file"
        enhancer_class_override: "Annotation file to override enhancer class assignment"
        supplementary_features: "Additional features to count over regions"
        cellType: "Name of cell type"
        HiCdirTar: "tarball of the bedgraph directory"
    }

    input {
        File dnaseqbam
        File dnaseqbam_index
        File chrom_sizes
        #Boolean? is_paired_end TODO add a optional for this?
        File? regions_blacklist
        File? regions_whitelist
        File genes_bed
        File? genes_for_class_assignment
        File? ubiq_genes
        String? gene_name_annotations
        String? primary_gene_identifier
        File h3k27ac_bam
        File h3k27ac_bam_index
        File? expression_table
        File? qnorm
        Int? tss_slop_for_class_assignment
        #Boolean? skip_rpkm_quantile
        #Boolean? use_secondary_counting_method
        File? enhancer_class_override
        File? supplementary_features
        String cellType = "defCellType"
        File HiCdirTar
        Float? pval_cutoff
        Int? nStrongestPeaks
        Int? peakExtendFromSummit
    }

    call makeCandidateRegions {
       input:
           bam = dnaseqbam,
           bam_index = dnaseqbam_index,
           chrom_sizes = chrom_sizes,
           #is_paired_end = is_paired_end,
           pval_cutoff = pval_cutoff,
           nStrongestPeaks = nStrongestPeaks,
           peakExtendFromSummit = peakExtendFromSummit,
           regions_blacklist = regions_blacklist,
           regions_whitelist = regions_whitelist,
    }

    call runNeighborhoods {
       input:
           candidate_enhancer_regions = makeCandidateRegions.candidateRegions,
           genes_bed = genes_bed,
           genes_for_class_assignment = genes_for_class_assignment,
           ubiquitously_expressed_genes = ubiq_genes,
           gene_name_annotations = gene_name_annotations,
           primary_gene_identifier = primary_gene_identifier,
           h3k27ac_bam = h3k27ac_bam,
           h3k27ac_bam_index = h3k27ac_bam_index,
           dnase_bam = dnaseqbam,
           dnase_bam_index = dnaseqbam_index,
           expression_table = expression_table,
           qnorm = qnorm,
           tss_slop_for_class_assignment = tss_slop_for_class_assignment,
           #skip_rpkm_quantile = skip_rpkm_quantile,
           #use_secondary_counting_method = use_secondary_counting_method,
           chromosome_sizes = chrom_sizes,
           enhancer_class_override = enhancer_class_override,
           supplementary_features = supplementary_features,
           cellType = cellType
    }

    call makePrediction {
        input:
            enhancerList = runNeighborhoods.enhancerList,
            geneList = runNeighborhoods.geneList,
            HiCdirTar = HiCdirTar,
            cellType = cellType
    }

    output {
       File candidateRegions = makeCandidateRegions.candidateRegions
    }
}


    task makeCandidateRegions {
        input {
            File bam
            File bam_index
            File chrom_sizes
            #Boolean? is_paired_end
            File? regions_blacklist
            File? regions_whitelist
            Float? pval_cutoff
            Int? peakExtendFromSummit
            Int? nStrongestPeaks
        }

        String docker_image = "quay.io/nbarkas/abc-general-container:latest"
        Int num_threads = 1
        String mem_size = "1 GB"


        command {
            set -euo pipefail

            mkdir outputs

            python3 /usr/src/app/src/makeCandidateRegions.py \
                --bam ~{bam} \
                --outDir outputs \
                --chrom_sizes ~{chrom_sizes} \
                ${"--regions_blacklist=" + regions_blacklist} \
                ${"--regions_whitelist=" + regions_whitelist} \
                ${"--pval_cutoff=" + pval_cutoff} \
                ${"--peakExtendFromSummit=" + peakExtendFromSummit} \
                ${"--nStrongestPeaks=" + nStrongestPeaks} \
        }
        output {
            # TODO: Add all the outputs
            File candidateRegions = "outputs/" + basename(bam, ".bam") + ".macs2_peaks.narrowPeak.candidateRegions.bed"
        }
        runtime {
            docker: docker_image
            cpu: num_threads
            memory: mem_size
            disks: "local-disk " + ceil((size(bam, "GiB")) * 1.2) + " HDD"
        }
    }



task runNeighborhoods {
    input {
       File candidate_enhancer_regions
       File genes_bed
       File? genes_for_class_assignment
       File? ubiquitously_expressed_genes
       String? gene_name_annotations
       String? primary_gene_identifier
       File h3k27ac_bam
       File h3k27ac_bam_index
       File dnase_bam # TODO this might not be required if -- ATAC flag is used
       File dnase_bam_index
       File? expression_table
       File? qnorm
       Int? tss_slop_for_class_assignment
       #Boolean? skip_rpkm_quantile
       #Boolean? use_secondary_counting_method
       File chromosome_sizes
       File? enhancer_class_override
       File? supplementary_features
       String cellType = "defCellType"
    }

        String docker_image = "quay.io/nbarkas/abc-general-container:latest"
        Int num_threads = 1
        String mem_size = "1 GB"

    command {
        set -euo pipefail

        python3 /usr/src/app/src/run.neighborhoods.py \
            --candidate_enhancer_regions ~{candidate_enhancer_regions} \
            --genes ~{genes_bed} \
            --outdir outputs/ \
            ${"--genes_for_class_assignment=" + genes_for_class_assignment} \
            ${"--ubiquitously_expressed_genes=" + ubiquitously_expressed_genes} \
            ${"--gene_name_annotations=" + gene_name_annotations} \
            ${"--primary_gene_identifier=" + primary_gene_identifier} \
            --H3K27ac ~{h3k27ac_bam} \
            --DHS ~{dnase_bam} \ ## TODO THIS IS WRONG, what about --ATAC flag
            ${"--expression_table=" + expression_table} \
            ${"--qnorm=" + qnorm} \
            ${"--tss_slop_for_class_assignment=" + tss_slop_for_class_assignment} \
            --chrom_sizes ~{chromosome_sizes} \
            ${"--enhancer_class_override=" + enhancer_class_override} \
            ${"--supplementary_features=" + supplementary_features} \
            --cellType ~{cellType}
    }
    output {
        # TODO: add remain outputs
        File enhancerList = "outputs/EnhancerList.txt"
        File geneList = "outputs/GeneList.txt"
    }
    runtime {
        docker: docker_image
        cpu: num_threads
        memory: mem_size
        disks: "local-disk " + ceil((size(dnase_bam, "GiB") + size(h3k27ac_bam, "GiB")) * 1.2) + " HDD"
    }
}

task makePrediction {
    input {
        File enhancerList
        File geneList
        File HiCdirTar
        Float threshold = "0.022"
        String cellType
    }

    String docker_image = "quay.io/nbarkas/abc-general-container:latest"
    Int num_threads = 1
    String mem_size = "1 GB"

    command {
        set -euo pipefail
        tar -xf ~{HiCdirTar}
        python3 /usr/src/app/src/predict.py \
            --enhancers ~{enhancerList} \
            --genes ~{geneList} \
            --HiCdir "bedgraph" \
            --scale_hic_using_powerlaw \
            --threshold ~{threshold} \
            --cellType ~{cellType} \
            --outdir outputs/
    }
    output {

    }
    runtime {
        docker: docker_image
        cpu: num_threads
        memory: mem_size
        disks: "local-disk " + ceil(size(HiCdirTar, "GiB")) * 3 + " HDD"
    }
}
