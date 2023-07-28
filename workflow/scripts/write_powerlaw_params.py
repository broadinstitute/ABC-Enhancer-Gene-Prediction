import click
import os
import pandas as pd
import subprocess


def exec_compute_powerlaw_fit_script(
    hic_dir: str, hic_type: str, hic_resolution: str, output_file: str
):
    output_dir = os.path.dirname(output_file)
    subprocess.run(
        [
            "python",
            "workflow/scripts/compute_powerlaw_fit_from_hic.py",
            "--hic_dir",
            hic_dir,
            "--hic_type",
            hic_type,
            "--hic_resolution",
            str(hic_resolution),
            "--outDir",
            output_dir,
        ]
    )


@click.command()
@click.option("--hic_dir", type=str)
@click.option("--hic_type", type=str)
@click.option("--hic_resolution", type=str)
@click.option("--hic_gamma", type=float)
@click.option("--hic_scale", type=float)
@click.option("--output_file", type=str)
def main(hic_dir, hic_type, hic_resolution, hic_gamma, hic_scale, output_file):
    if hic_gamma and hic_scale:
        pd.DataFrame({"hic_gamma": [hic_gamma], "hic_scale": [hic_scale]}).to_csv(
            output_file, sep="\t", index=False
        )
    else:
        if not (hic_dir and hic_type and hic_resolution):
            raise click.MissingParameter(
                "Must provide gamma/scale values or full HiC info (dir,type, and resolution)"
            )
        exec_compute_powerlaw_fit_script(hic_dir, hic_type, hic_resolution, output_file)


if __name__ == "__main__":
    main()
