version 1.0

workflow ValidateABC {
    meta {
        description: "Validate ABC pipeline outputs"
    }

    input {
        File candidateRegions 
    }

    call ValidateABCcheck {

    }
}

task ValidateABCcheck {
    input {
        File candidateRegions
    }

    String docker_image = "ubuntu:18.04"
    Int num_threads = 1
    String mem_size = "1 GB"

    command {
        set -euo pipefail
        cat ~{candidateRegions} | md5sum | cut -f 2 > checksum.txt
    }

    output {
        File checksum = "checksum.txt"
    }

    runtime {
        docker: docker_image
        cpu: num_threads
        memory: mem_size
        disks: "local-disk " + ceil((size(bam, "GiB")) * 1.2) + " HDD"
    }
}