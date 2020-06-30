# To Do:
# 1. Support multiple DHS, K27ac files
# 2. Support more input flags
# 3. Python package versions in DockerFile?
# 4. Add more output files

workflow ABCpipeline {
  # If this is defined makeCandidateRegions is not run
  File? precomputedCandidateRegions
  
  # Inputs to makeCandidateRegions
  String? macs_type
  String? genome_string = "hs"
  Boolean? is_paired_end
  File? regions_blacklist
  File? regions_whitelist
  Float? pval_cutoff = 0.1
  Int? nStrongestPeaks = 3000
  Int? peakExtendFromSummit = 250

  # Used in makeCandidateRegions and runNeighborhoods
  File chrom_sizes
  File chrom_sizes_bed = chrom_sizes + ".bed"

  File new_neighborhoods_py


  File genes_bed
  File? genes_for_class_assignment
  File? ubiq_genes
  String? gene_name_annotations
  String? primary_gene_identifier
  File? dhs_bam
  File? dhs_bam_index
  File? atac_bam
  File? atac_bam_index
  String? default_accessibility_feature = "DHS"
  File? expression_table
  File? qnorm
  Int? tss_slop_for_class_assignment
  Boolean? skip_rpkm_quantile
  Boolean? use_secondary_counting_method
  File? enhancer_class_override
  File? supplementary_features
  String cellType = "defCellType"

  Int? window
  Float threshold = 0.022
  Float? tss_hic_contribution
  Int? hic_pseudocount_distance
  Boolean? scale_hic_using_powerlaw
  Float? hic_gamma
  Float? hic_gamma_reference
  Boolean? run_all_genes
  Float? expression_cutoff
  Float? promoter_activity_quantile_cutoff
  Boolean? skip_gene_files
  Boolean? skinny_gene_files
  Boolean? make_all_putative
  Int? tss_slop
  Boolean? include_chrY

  #Runtime
  String docker_image = "quay.io/jnasser/abc-container"
  Int num_threads = 4
  String mem_size = "32G"
  String disks = "local-disk 100 SSD"

  # If candidate regions are not defined run makeCandidateRegions
  # to generate them
  if ( !defined(precomputedCandidateRegions) ) {
      call runMACS {
          input:
              bam = dhs_bam,
              bam_index = dhs_bam_index,
              macs_type = macs_type,
              genome_string = genome_string,
              pval_cutoff = pval_cutoff,
              docker_image = docker_image,
              num_threads = num_threads,
              mem_size = mem_size,
              disks = disks
      }

      if (length(runMACS.narrowPeak) >= 1) {
        File narrowPeak2 = runMACS.narrowPeak[0]
      }

      call sortMACS {
          input:
              narrowPeak = narrowPeak2,
              chrom_sizes = chrom_sizes,
              docker_image = docker_image,
              num_threads = num_threads,
              mem_size = mem_size,
              disks = disks
      }

      if (length(sortMACS.narrowPeakSorted) >= 1) {
        File narrowPeakSorted2 = sortMACS.narrowPeakSorted[0]
      }

      call makeCandidateRegions {
         input:
             narrowPeakSorted = narrowPeakSorted2,
             bam = select_first([dhs_bam]),
             bam_index = select_first([dhs_bam_index]),
             chrom_sizes = chrom_sizes,
             chrom_sizes_bed = chrom_sizes_bed,
             pval_cutoff = pval_cutoff,
             nStrongestPeaks = nStrongestPeaks,
             peakExtendFromSummit = peakExtendFromSummit,
             regions_blacklist = regions_blacklist,
             regions_whitelist = regions_whitelist,
              docker_image = docker_image,
              num_threads = num_threads,
              mem_size = mem_size,
              disks = disks,
              new_neighborhoods_py = new_neighborhoods_py
      }

      if (length(makeCandidateRegions.candidateRegions) >= 1) {
        File candidateRegions2 = makeCandidateRegions.candidateRegions[0]
      }
  }    

    call runNeighborhoods {
       input:
           candidate_enhancer_regions = select_first([precomputedCandidateRegions, candidateRegions2]),
           genes_bed = genes_bed,
           genes_for_class_assignment = genes_for_class_assignment,
           ubiquitously_expressed_genes = ubiq_genes,
           gene_name_annotations = gene_name_annotations,
           primary_gene_identifier = primary_gene_identifier,
           dhs_bam = dhs_bam,
           dhs_bam_index = dhs_bam_index,
           atac_bam = atac_bam,
           atac_bam_index = atac_bam_index,
           default_accessibility_feature = default_accessibility_feature,
           expression_table = expression_table,
           qnorm = qnorm,
           tss_slop_for_class_assignment = tss_slop_for_class_assignment,
           skip_rpkm_quantile = skip_rpkm_quantile,
           use_secondary_counting_method = use_secondary_counting_method,
           chrom_sizes = chrom_sizes,
           chrom_sizes_bed = chrom_sizes_bed,
           enhancer_class_override = enhancer_class_override,
           supplementary_features = supplementary_features,
           cellType = cellType,
              docker_image = docker_image,
              num_threads = num_threads,
              mem_size = mem_size,
              disks = disks,
              new_neighborhoods_py = new_neighborhoods_py
              
    }

    call makePrediction {
        input:
            enhancerList = runNeighborhoods.enhancerList,
            geneList = runNeighborhoods.geneList,
            window = window,
            threshold = threshold,
            cellType = cellType,
            tss_hic_contribution = tss_hic_contribution,
            hic_pseudocount_distance = hic_pseudocount_distance,
            scale_hic_using_powerlaw = scale_hic_using_powerlaw,
            hic_gamma = hic_gamma,
            hic_gamma_reference = hic_gamma_reference,
            run_all_genes = run_all_genes,
            expression_cutoff = expression_cutoff,
            promoter_activity_quantile_cutoff = promoter_activity_quantile_cutoff,
            skip_gene_files = skip_gene_files,
            skinny_gene_files = skinny_gene_files,
            make_all_putative = make_all_putative,
            tss_slop = tss_slop,
            include_chrY = include_chrY,
              docker_image = docker_image,
              num_threads = num_threads,
              mem_size = mem_size,
              disks = disks
    }
}


task runMACS {
        File bam
        File bam_index
        String? macs_type
        String? genome_string
        Float pval_cutoff

        String docker_image
        Int num_threads
        String mem_size
        String disks

    command {
        set -euo pipefail

        mkdir Peaks/

        macs2 callpeak \
          -t ${bam} \
          -n "macs2" \
          -f ${macs_type} \
          -g ${genome_string} \
          -p ${pval_cutoff} \
          --call-summits

    }
    output {
        # TODO: Add all the outputs
        # File candidateRegions = "Peaks/" + basename(bam, ".bam") + ".macs2_peaks.narrowPeak.candidateRegions.bed"
        Array[File] narrowPeak = glob("macs2_peaks.narrowPeak")
    }
    runtime {
        docker: docker_image
        cpu: num_threads
        memory: "32G"
        disks: disks
    }
}

task sortMACS {
        File narrowPeak
        File chrom_sizes

        String docker_image
        Int num_threads
        String mem_size
        String disks

    command {
        set -euo pipefail

        cut -f 1 ${chrom_sizes} > /tmp/known_chroms
        grep -w -f /tmp/known_chroms ${narrowPeak} | \
        bedtools sort \
        -faidx ${chrom_sizes} \
        -i - > "macs2_peaks.narrowPeak.sorted"
    }

    output {
        # TODO: Add all the outputs
        Array[File] narrowPeakSorted = glob("macs2_peaks.narrowPeak.sorted")
    }

    runtime {
        docker: docker_image
        cpu: num_threads
        memory: mem_size
        disks: disks
    }
}

task makeCandidateRegions {
       File narrowPeakSorted
       File bam
       File bam_index
       File chrom_sizes
       File chrom_sizes_bed
       Float pval_cutoff
       Int nStrongestPeaks
       Int peakExtendFromSummit
       File? regions_blacklist
       File? regions_whitelist
       File new_neighborhoods_py

      String docker_image
      Int num_threads
      String mem_size
      String disks

    command {
          set -euo pipefail

          cp ${new_neighborhoods_py} /usr/src/app/src/neighborhoods.py

          pip3 install pysam

          python3 /usr/src/app/src/makeCandidateRegions.py \
              --narrowPeak ${narrowPeakSorted} \
              --bam ${bam} \
              --chrom_sizes ${chrom_sizes} \
              --outDir "." \
              --regions_blacklist ${regions_blacklist} \
              --regions_whitelist ${regions_whitelist} \
              --peakExtendFromSummit ${peakExtendFromSummit} \
              --nStrongestPeaks ${nStrongestPeaks} 
    }

    output {
        Array[File] candidateRegions = glob("macs2_peaks.narrowPeak.sorted.candidateRegions.bed")
        Array[File] rawCounts = glob("*.Counts.bed")
    }

    runtime {
        docker: docker_image
        cpu: num_threads
        memory: mem_size
        disks: disks
    }
}

task runNeighborhoods {
       File candidate_enhancer_regions
       File genes_bed
       File? genes_for_class_assignment
       File? ubiquitously_expressed_genes
       String? gene_name_annotations
       String? primary_gene_identifier
       File? dhs_bam
       File? dhs_bam_index
       File? atac_bam
       File? atac_bam_index
       String? default_accessibility_feature
       File? expression_table
       File? qnorm
       Int? tss_slop_for_class_assignment
       Boolean? skip_rpkm_quantile
       Boolean? use_secondary_counting_method
       File chrom_sizes
       File chrom_sizes_bed
       File? enhancer_class_override
       File? supplementary_features
       String cellType = "defCellType"

       File new_neighborhoods_py

       String docker_image
       Int num_threads
       String mem_size
       String disks

    ## TODO: check about --ATAC flag
    command {
        set -euo pipefail

        mkdir Neighborhoods/

        cp ${new_neighborhoods_py} /usr/src/app/src/neighborhoods.py

        pip3 install pysam

        python3 /usr/src/app/src/run.neighborhoods.py \
            --candidate_enhancer_regions ${candidate_enhancer_regions} \
            --genes ${genes_bed} \
            --ubiquitously_expressed_genes ${ubiquitously_expressed_genes} \
            --DHS ${dhs_bam} \
            --expression_table ${expression_table} \
            --chrom_sizes ${chrom_sizes} \
            --cellType ${cellType} \
            --outdir "."
    }
    output {
        # TODO: add remain outputs
        File enhancerList = "EnhancerList.txt"
        File geneList = "GeneList.txt"
    }
    runtime {
        docker: docker_image
        cpu: num_threads
        memory: mem_size
        disks: disks
    }
}

task makePrediction {
        File enhancerList
        File geneList
        Int? window
        Float threshold
        String cellType
        Float? tss_hic_contribution
        Int? hic_pseudocount_distance
        Boolean? scale_hic_using_powerlaw
        Float? hic_gamma
        Float? hic_gamma_reference
        Boolean? run_all_genes
        Float? expression_cutoff
        Float? promoter_activity_quantile_cutoff
        Boolean? skip_gene_files
        Boolean? skinny_gene_files
        Boolean? make_all_putative
        Int? tss_slop
        Boolean? include_chrY

        String docker_image
        Int num_threads
        String mem_size
        String disks

    command {
        set -euo pipefail
        mkdir Predictions/
        python3 /usr/src/app/src/predict.py \
            --enhancers ${enhancerList} \
            --genes ${geneList} \
            --window ${window} \
            --threshold ${threshold} \
            --cellType ${cellType} \
            --hic_resolution 5000 \
            --hic_gamma ${hic_gamma} \
            --hic_gamma_reference ${hic_gamma_reference} \
            --make_all_putative \
            --score_column "powerlaw.Score" \
            --outdir "."
    }
    output {
        File enhancerPredictions = "EnhancerPredictions.txt"
        File enhancerPredictionsFull = "EnhancerPredictionsFull.txt"
        File params = "parameters.predict.txt"
        File bedpe = "EnhancerPredictions.bedpe"
        File geneStats = "GenePredictionStats.txt"
        Array[File] allPutative = glob("EnhancerPredictionsAllPutative.txt.gz")

    }
    runtime {
        docker: docker_image
        cpu: num_threads
        memory: "32G"
        disks: disks
    }
}