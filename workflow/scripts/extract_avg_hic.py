import gzip
import os
import subprocess

import click


@click.command()
@click.option("--avg_hic_bed_file", type=str, required=True)
@click.option("--output_dir", type=str, default=".")
def main(avg_hic_bed_file, output_dir):
    output_dir = os.path.join(output_dir, "AvgHiC")
    os.makedirs(output_dir, exist_ok=True)
    file_handles = {}
    with gzip.open(avg_hic_bed_file, "rt") as f:
        for line in f:
            if line.startswith("#"):  # header line
                continue
            chrom = line.split("\t")[0]
            hic_info = "\t".join(line.split("\t")[1:])

            if chrom not in file_handles:
                print(f"Writing lines for {chrom}")
                os.makedirs(os.path.join(output_dir, chrom), exist_ok=True)
                chrom_file = os.path.join(
                    os.path.join(output_dir, chrom), f"{chrom}.bed"
                )
                file_handles[chrom] = open(chrom_file, "w")
            file_handles[chrom].write(hic_info)

    # Close all file handles
    for fh in file_handles.values():
        filename = fh.name
        cmd = f"pigz -f {filename}"
        print(f"Gzipping {filename}")
        subprocess.run(cmd, shell=True, check=True)
        fh.close()


if __name__ == "__main__":
    main()
