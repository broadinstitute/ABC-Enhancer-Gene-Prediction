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
        #gene_name_annotations: "Comma delimited string of names corresponding to the gene identifiers present in the name field of the gene annotation bed file"
        primary_gene_identifier: "Primary identifier used to identify genes. Must be present in gene_name_annotations. The primary identifier must be unique"

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
        #File? genes_for_class_assignment
        File ubiq_genes
        #String gene_name_annotations
        File h3k27ac_bam
        File h3k27ac_bam_index
        File expression_table
        File HiCdirTar
        String cellType = "defCellType"
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
           #genes_for_class_assignment = genes_for_class_assignment,
           h3k27ac_bam = h3k27ac_bam,
           h3k27ac_bam_index = h3k27ac_bam_index,
           dnase_bam = dnaseqbam,
           dnase_bam_index = dnaseqbam_index,
           expression_table = expression_table,
           chromosome_sizes = chrom_sizes,
           ubiq_genes = ubiq_genes,
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
       File h3k27ac_bam
       File h3k27ac_bam_index
       File dnase_bam
       File dnase_bam_index
       File expression_table
       File chromosome_sizes
       File ubiq_genes
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
            --H3K27ac ~{h3k27ac_bam} \
            --DHS ~{dnase_bam} \
            --expression_table ~{expression_table} \
            --chrom_sizes ~{chromosome_sizes} \
            --ubiquitously_expressed_genes ~{ubiq_genes} \
            --cellType ~{cellType} \
            --outdir outputs/
    }
    output {
        # TODO: add remain outpus
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
