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
        dhs_bam: "Comma delimited string of DHS .bam files. Either ATAC or DHS must be provided"
        atac_bam: "Comma delimited string of ATAC .bam files. Either ATAC or DHS must be provided"
        default_accessibility_feature: "If both ATAC and DHS are provided, this flag must be set to either 'DHS' or 'ATAC' signifying which datatype to use in computing activity"
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
        window: "Make predictions for all candidate elements within this distance of the gene's TSS"
        threshold:"Threshold on ABC Score to call a predicted positive"
        HiCdirTar: "tarball of the bedgraph directory with hic bedgraphs"
        tss_hic_contribution: "Weighting of diagonal bin of hic matrix as a percentage of its neighbors"
        hic_pseudocount_distance: "A pseudocount is added equal to the powerlaw fit at this distance"
        scale_hic_using_powerlaw: "Scale Hi-C values using powerlaw relationship"
        hic_gamma: "powerlaw exponent of hic data. Must be positive"
        hic_gamma_reference: "powerlaw exponent to scale to. Must be positive"
        #run_all_genes: "Do not check for gene expression, make predictions for all genes"
        expression_cutoff: "Make predictions for genes with expression higher than this value"
        promoter_activity_quantile_cutoff: "Quantile cutoff on promoter activity. Used to consider a gene 'expressed' in the absence of expression data"
        #skip_gene_files: "Do not make individual gene files"
        #skinny_gene_files: "Use subset of columns for genes files"
        #make_all_putative: "Make big file with concatenation of all genes file"
        tss_slop: "Distance from tss to search for self-promoters"
        #include_chrY: "Include Y chromosome"
    }

    input {
        File dnaseqbam
        File dnaseqbam_index
        File chrom_sizes
        #Boolean? is_paired_end TODO add a optional for this?
        File? regions_blacklist
        File? regions_whitelist
        Float? pval_cutoff = 0.1
        Int? nStrongestPeaks = 3000
        Int? peakExtendFromSummit = 250

        File genes_bed
        File? genes_for_class_assignment
        File? ubiq_genes
        String? gene_name_annotations
        String? primary_gene_identifier
        File h3k27ac_bam
        File h3k27ac_bam_index
        File? dhs_bam
        File? dhs_bam_index
        File? atac_bam
        File? atac_bam_index
        String? default_accessibility_feature
        File? expression_table
        File? qnorm
        Int? tss_slop_for_class_assignment
        #Boolean? skip_rpkm_quantile
        #Boolean? use_secondary_counting_method
        File? enhancer_class_override
        File? supplementary_features
        String cellType = "defCellType"

        Int? window
        Float threshold = 0.022
        File HiCdirTar
        Float? tss_hic_contribution
        Int? hic_pseudocount_distance
        #Boolean? scale_hic_using_powerlaw
        Float? hic_gamma
        Float? hic_gamma_reference
        #Boolean? run_all_genes
        Float? expression_cutoff
        Float? promoter_activity_quantile_cutoff
        #Boolean? skip_gene_files
        #Boolean? skinny_gene_files
        #Boolean? make_all_putative
        Int? tss_slop
        #Boolean? include_chrY
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
           dhs_bam = dhs_bam,
           dhs_bam_index = dhs_bam_index,
           atac_bam = atac_bam,
           atac_bam_index = atac_bam_index,
           default_accessibility_feature = default_accessibility_feature,
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
            window = window,
            threshold = threshold,
            cellType = cellType,
            HiCdirTar = HiCdirTar,
            tss_hic_contribution = tss_hic_contribution,
            hic_pseudocount_distance = hic_pseudocount_distance,
            #scale_hic_using_powerlaw = scale_hic_using_powerlaw,
            #hic_gamma = hic_gamma,
            #hic_gamma_reference = hic_gamma_reference,
            #run_all_genes = run_all_genes,
            expression_cutoff = expression_cutoff,
            promoter_activity_quantile_cutoff = promoter_activity_quantile_cutoff,
            #skip_gene_files = skip_gene_files,
            #skinny_gene_files = skinny_gene_files,
            #make_all_putative = make_all_putative,
            tss_slop = tss_slop,
            #include_chrY = include_chrY,
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
       File? dhs_bam
       File? dhs_bam_index
       File? atac_bam
       File? atac_bam_index
       String? default_accessibility_feature
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

       ## TODO THIS IS WRONG, what about --ATAC flag
    command {
        set -euo pipefail

        python3 /usr/src/app/src/run.neighborhoods.py \
            --candidate_enhancer_regions ~{candidate_enhancer_regions} \
            --genes ~{genes_bed} \
            ${"--genes_for_class_assignment=" + genes_for_class_assignment} \
            ${"--ubiquitously_expressed_genes=" + ubiquitously_expressed_genes} \
            ${"--gene_name_annotations=" + gene_name_annotations} \
            ${"--primary_gene_identifier=" + primary_gene_identifier} \
            --H3K27ac ~{h3k27ac_bam} \
            ${"--DHS=" + dhs_bam} \
            ${"--ATAC=" + atac_bam} \
            ${"--default_accessibility_feature=" + default_accessibility_feature} \
            ${"--expression_table=" + expression_table} \
            ${"--qnorm=" + qnorm} \
            ${"--tss_slop_for_class_assignment=" + tss_slop_for_class_assignment} \
            --chrom_sizes ~{chromosome_sizes} \
            ${"--enhancer_class_override=" + enhancer_class_override} \
            ${"--supplementary_features=" + supplementary_features} \
            --cellType ~{cellType} \
            --outdir outputs/
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
        disks: "local-disk " + ceil((size(dhs_bam, "GiB") + size(h3k27ac_bam, "GiB")) * 1.2) + " HDD"
    }
}

task makePrediction {
    input {
        File enhancerList
        File geneList
        Int? window
        Float threshold
        String cellType
        File HiCdirTar
        Float? tss_hic_contribution
        Int? hic_pseudocount_distance
        #Boolean? scale_hic_using_powerlaw
        Float? hic_gamma
        Float? hic_gamma_reference
        #Boolean? run_all_genes
        Float? expression_cutoff
        Float? promoter_activity_quantile_cutoff
        #Boolean? skip_gene_files
        #Boolean? skinny_gene_files
        #Boolean? make_all_putative
        Int? tss_slop
        #Boolean? include_chrY
    }

    String docker_image = "quay.io/nbarkas/abc-general-container:latest"
    Int num_threads = 1
    String mem_size = "1 GB"

    ## TODO this is wrong-- need to deal with Booleans
    command {
        set -euo pipefail
        tar -xf ~{HiCdirTar}
        python3 /usr/src/app/src/predict.py \
            --enhancers ~{enhancerList} \
            --genes ~{geneList} \
            ${"--window=" + window} \
            --threshold ~{threshold} \
            --cellType ~{cellType} \
            --HiCdir "bedgraph" \
            ${"--tss_hic_contribution=" + tss_hic_contribution} \
            ${"--hic_gamma=" + hic_gamma} \
            ${"--hic_gamma_reference=" + hic_gamma_reference} \
            ${"--expression_cutoff=" + expression_cutoff} \
            ${"--promoter_activity_quantile_cutoff=" + promoter_activity_quantile_cutoff} \
            ${"--tss_slop=" + tss_slop} \
            --scale_hic_using_powerlaw \
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
