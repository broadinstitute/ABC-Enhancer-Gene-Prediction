import logging
import os
import time
import unittest
from typing import Callable, List

import pandas as pd
import yaml
from test_utils import (
    compare_files,
    get_biosample_names,
    get_filtered_dataframe,
    run_cmd,
)

logging.basicConfig(level=logging.INFO)

CONFIG_FILE = "tests/config/generic_config.yml"
with open(CONFIG_FILE, "r") as file:
    CONFIG = yaml.safe_load(file)
COLUMNS_TO_COMPARE = [
    "chr",
    "start",
    "end",
    "name",
    "class",
    "TargetGene",
    "ABC.Score.Numerator",
    "ABC.Score",
]
TEST_OUTPUT_DIR = CONFIG["predictions_results_dir"]
EXPECTED_OUTPUT_DIR = f"tests/expected_output/{CONFIG['TEST_CONFIG_NAME']}"
PREDICTION_FILE = "Predictions/EnhancerPredictionsAllPutative.txt.gz"
INTERMEDIATE_FILES = [
    # EnhancerList is not good to compare b/c column names change randomly
    # "Neighborhoods/EnhancerList.txt",
    "Neighborhoods/GeneList.txt",
    "Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed",
]


class TestFullABCRun(unittest.TestCase):
    def compare_intermediate_files(self, biosample: str) -> None:
        for file in INTERMEDIATE_FILES:
            test_file = os.path.join(TEST_OUTPUT_DIR, biosample, file)
            expected_file = os.path.join(EXPECTED_OUTPUT_DIR, biosample, file)
            if not compare_files(test_file, expected_file):
                logging.error(f"Intermediate file comparison failed for {file}")

    def compare_prediction_file(self, biosample: str) -> None:
        test_file = os.path.join(TEST_OUTPUT_DIR, biosample, PREDICTION_FILE)
        expected_file = os.path.join(EXPECTED_OUTPUT_DIR, biosample, PREDICTION_FILE)
        pd.testing.assert_frame_equal(
            get_filtered_dataframe(test_file, COLUMNS_TO_COMPARE),
            get_filtered_dataframe(expected_file, COLUMNS_TO_COMPARE),
        )

    def run_test(self, config_file: str) -> None:
        start = time.time()
        cmd = f"snakemake -j4 -F --configfile {config_file}"
        run_cmd(cmd)
        time_taken = time.time() - start

        biosample_names = get_biosample_names(CONFIG["biosamplesTable"])
        for biosample in biosample_names:
            self.compare_intermediate_files(biosample)
            self.compare_prediction_file(biosample)

        # Make sure the test doesn't take too long
        # May need to adjust as more biosamples are added, but we should keep
        # tests quick (so don't run ABC on all chromosomes)
        max_time = 60 * 3  # 3 min
        self.assertLessEqual(
            time_taken,
            max_time,
            msg=f"Running ABC took too long: {int(time_taken)} seconds",
        )

    def test_full_abc_run(self) -> None:
        self.run_test(CONFIG_FILE)


if __name__ == "__main__":
    unittest.main()
