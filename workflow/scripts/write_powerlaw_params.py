import os
import subprocess

import click
import pandas as pd


def exec_compute_powerlaw_fit_script(
    hic_file: str, hic_resolution: str, output_file: str
):
    output_dir = os.path.dirname(output_file)
    subprocess.run(
        [
            "python",
            "workflow/scripts/compute_powerlaw_fit_from_hic.py",
            "--hic_file",
            hic_file,
            "--hic_resolution",
            str(hic_resolution),
            "--outDir",
            output_dir,
        ]
    )


@click.command()
@click.option("--hic_file", type=str)
@click.option("--hic_resolution", type=str)
@click.option("--hic_gamma", type=float)
@click.option("--hic_scale", type=float)
@click.option("--output_file", type=str)
def main(hic_file, hic_resolution, hic_gamma, hic_scale, output_file):
    if hic_gamma and hic_scale:
        pd.DataFrame({"hic_gamma": [hic_gamma], "hic_scale": [hic_scale]}).to_csv(
            output_file, sep="\t", index=False
        )
    else:
        if not (hic_file and hic_resolution):
            raise click.MissingParameter(
                "Must provide gamma/scale values or full HiC info (dir,type, and resolution)"
            )
        exec_compute_powerlaw_fit_script(hic_file, hic_resolution, output_file)


if __name__ == "__main__":
    main()
