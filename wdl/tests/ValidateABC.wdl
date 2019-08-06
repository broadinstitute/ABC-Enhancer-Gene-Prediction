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
        String candidateRegions_targetMD5
    }

    String docker_image = "ubuntu:18.04"
    Int num_threads = 1
    String mem_size = "1 GB"

    command {
        set -euo pipefail
        candidateRegions_hash=$(cat ~{candidateRegions} | md5sum | awk '{print $1}')

        fail=false

        if [ "$candidateRegions_hash" == "~{candidateRegions_targetMD5}" ]; then
            >&2 echo "candidateRegions hash ($candidateRegions_hash) did not match the expected hash (~{candidateRegions_targetMD5})"
            fail=true
        fi

        if [ $fail == "true" ]; then exit 1; fi
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