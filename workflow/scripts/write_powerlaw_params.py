import click
import os
import pandas as pd
import subprocess


def exec_compute_powerlaw_fit_script(hic_dir: str, hic_type: str, output_file: str):
    output_dir = os.path.dirname(output_file)
    subprocess.run(
        [
            "python",
            "workflow/scripts/compute_powerlaw_fit_from_hic.py",
            "--hic_dir",
            hic_dir,
            "--hic_type",
            hic_type,
            "--outDir",
            output_dir,
        ]
    )


def is_float(string):
    try:
        float_value = float(string)
        return True
    except ValueError:
        return False


@click.command()
@click.option("--hic_dir", type=str)
@click.option("--hic_type", type=str)
@click.option("--hic_gamma", type=str)
@click.option("--hic_scale", type=str)
@click.option("--output_file", type=str)
def main(hic_dir, hic_type, hic_gamma, hic_scale, output_file):
    if is_float(hic_gamma) and is_float(hic_scale):
        pd.DataFrame({"hic_gamma": [hic_gamma], "hic_scale": [hic_scale]}).to_csv(
            output_file, sep="\t", index=False
        )
    else:
        if not hic_dir or not hic_type:
            raise click.MissingParameter("Must provide gamma/scale values or HiC dir")
        exec_compute_powerlaw_fit_script(hic_dir, hic_type, output_file)


if __name__ == "__main__":
    main()
