import logging
import os
import time
import unittest
from typing import Dict

import numpy as np
import pandas as pd
import glob
import yaml
from utils import get_biosample_names, get_filtered_dataframe, read_file, run_cmd

logging.basicConfig(level=logging.INFO)

CONFIG_FILE = "tests/config/generic_config.yml"
with open(CONFIG_FILE, "r") as file:
    CONFIG = yaml.safe_load(file)
COLUMNS_TO_COMPARE: Dict[str, type] = {
    "chr": str,
    "start": np.int64,
    "end": np.int64,
    "name": str,
    "class": str,
    "TargetGene": str,
    "ABC.Score.Numerator": np.float64,
    "ABC.Score": np.float64,
    "powerlaw.Score": np.float64,
}
TEST_OUTPUT_DIR = CONFIG["results_dir"]
EXPECTED_OUTPUT_DIR = f"tests/expected_output/{CONFIG['TEST_CONFIG_NAME']}"
ALL_PUTATIVE_PRED_FILE = "Predictions/EnhancerPredictionsAllPutative.tsv.gz"
THRESHOLDED_PRED_FILE = (
    "Predictions/EnhancerPredictionsFull_threshold*_self_promoter.tsv"
)
INTERMEDIATE_FILES = [
    # EnhancerList is not good to compare b/c genicSymbol column isn't deterministic
    # "Neighborhoods/EnhancerList.txt",
    "Neighborhoods/GeneList.txt",
    "Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed",
]


class TestFullABCRun(unittest.TestCase):
    def compare_intermediate_files(self, biosample: str) -> None:
        for file in INTERMEDIATE_FILES:
            test_file = os.path.join(TEST_OUTPUT_DIR, biosample, file)
            expected_file = os.path.join(EXPECTED_OUTPUT_DIR, biosample, file)
            msg = f"Intermediate file comparison failed for {file}"
            test_contents = read_file(test_file)
            expected_contents = read_file(expected_file)
            self.assertEqual(test_contents, expected_contents, msg)

    def compare_all_prediction_file(self, biosample: str, pred_file) -> None:
        test_file = os.path.join(TEST_OUTPUT_DIR, biosample, pred_file)
        expected_file = os.path.join(EXPECTED_OUTPUT_DIR, biosample, pred_file)
        print(f"Comparing biosample: {biosample} for pred_file: {pred_file}")
        pd.testing.assert_frame_equal(
            get_filtered_dataframe(test_file, COLUMNS_TO_COMPARE),
            get_filtered_dataframe(expected_file, COLUMNS_TO_COMPARE),
        )

    def compare_thresholded_prediction_file(self, biosample: str) -> None:
        test_files = glob.glob(
            os.path.join(TEST_OUTPUT_DIR, biosample, THRESHOLDED_PRED_FILE)
        )
        expected_files = glob.glob(
            os.path.join(EXPECTED_OUTPUT_DIR, biosample, THRESHOLDED_PRED_FILE)
        )
        if len(test_files) != 1:
            raise Exception(
                f"Multiple or no test thresholded files found. Please clean up. {test_files}"
            )
        if len(expected_files) != 1:
            raise Exception(
                f"Multiple or no expected thresholded files found. Please clean up. {expected_files}"
            )
        test_file = test_files[0]
        expected_file = expected_files[0]
        print(
            f"Comparing biosample: {biosample} for pred_file: {os.path.basename(test_file)}"
        )
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
            self.compare_all_prediction_file(biosample, ALL_PUTATIVE_PRED_FILE)
            self.compare_thresholded_prediction_file(biosample)

        # Make sure the test doesn't take too long
        # May need to adjust as more biosamples are added, but we should keep
        # tests quick, so don't run ABC on all chromosomes
        max_time = 60 * 8  # 8 min
        self.assertLessEqual(
            time_taken,
            max_time,
            msg=f"Running ABC took too long: {int(time_taken)} seconds",
        )

    def test_full_abc_run(self) -> None:
        self.run_test(CONFIG_FILE)


if __name__ == "__main__":
    unittest.main()
