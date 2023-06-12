import gzip
import logging
import os
import subprocess
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Callable, List

import pandas as pd

logging.basicConfig(level=logging.INFO)

CONFIG_FILE = "tests/config.yaml"
TEST_OUTPUT_DIR = "tests/test_output"
EXPECTED_OUTPUT_DIR = "tests/expected_output"
PREDICTION_FILE = "K562/Predictions/EnhancerPredictionsAllPutative.txt.gz"
COLUMNS_TO_COMPARE = ['chr', 'start', 'end', 'name', 'class', 'TargetGene', "ABC.Score.Numerator", "ABC.Score"]
INTERMEDIATE_FILES=[
    # EnhancerList is not good to compare b/c column names change randomly
    # "K562/Neighborhoods/EnhancerList.txt", 
    "K562/Neighborhoods/GeneList.txt", 
    "K562/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed"
]

def read_file(file_name: str) -> List[str]:
    open_fn: Callable = open
    if file_name.endswith(".gz"):
        open_fn = gzip.open
        
    with open_fn(file_name, 'r') as f:
        return sorted(f.readlines())
    

def compare_files(test_file: str, expected_file: str) -> bool:
    test_contents = read_file(test_file)
    expected_contents = read_file(expected_file)
    return test_contents == expected_contents

def get_filtered_dataframe(file: str) -> pd.DataFrame:
    return pd.read_csv(file, sep='\t', compression='gzip')[COLUMNS_TO_COMPARE]

def run_cmd(cmd: str, raise_ex: bool = True) -> bool:
    try:
        subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error: {e.output}")
        if raise_ex:
            raise
        return False
    return True

class TestFullABCRun(unittest.TestCase):
    def test_full_abc_run(self):
        cmd = f"snakemake -j1 -F --configfile {CONFIG_FILE}"
        run_cmd(cmd)
        
        # compare intermediate files       
        for file in INTERMEDIATE_FILES:
            test_file = os.path.join(TEST_OUTPUT_DIR, file)
            expected_file = os.path.join(EXPECTED_OUTPUT_DIR, file)
            if not compare_files(test_file, expected_file):
                logging.error(f"Intermediate file comparison failed for {file}")
        
        test_file = os.path.join(TEST_OUTPUT_DIR, PREDICTION_FILE)
        expected_file = os.path.join(EXPECTED_OUTPUT_DIR, PREDICTION_FILE)
        pd.testing.assert_frame_equal(get_filtered_dataframe(test_file), get_filtered_dataframe(expected_file))
        
        
if __name__ == '__main__':
    unittest.main()