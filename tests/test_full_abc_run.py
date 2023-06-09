import unittest
import subprocess
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import gzip
from typing import Callable, List
import logging

logging.basicConfig(level=logging.INFO)

CONFIG_FILE = "tests/config.yaml"
TEST_OUTPUT_DIR = "tests/test_output"
EXPECTED_OUTPUT_DIR = "tests/expected_output"
MAIN_OUTPUT_FILE = "K562/Predictions/EnhancerPredictionsAllPutative.txt.gz"

OTHER_FILES_TO_COMPARE=[
    # EnhancerList is not good to compare b/c column names change randomly
    # "K562/Neighborhoods/EnhancerList.txt", 
    "K562/Neighborhoods/GeneList.txt", 
    "K562/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed"
]

def read_file_sorted(file_name: str) -> List[str]:
    open_fn: Callable = open
    if file_name.endswith(".gz"):
        open_fn = gzip.open
        
    with open_fn(file_name, 'r') as f:
        return sorted(f.readlines())
    

def compare_files(test_file: str, expected_file: str) -> bool:
    test_contents = read_file_sorted(test_file)
    expected_contents = read_file_sorted(expected_file)
    return test_contents == expected_contents

class TestFullABCRun(unittest.TestCase):

    def run_cmd(self, cmd: str, raise_ex: bool = True) -> bool:
        try:
            subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Error: {e.output}")
            if raise_ex:
                raise
            return False
        return True

    def test_full_abc_run(self):
        cmd = f"snakemake --use-conda -j1 -F --configfile {CONFIG_FILE}"
        self.run_cmd(cmd)
        
        # compare files       
        for file in OTHER_FILES_TO_COMPARE:
            test_file = os.path.join(TEST_OUTPUT_DIR, file)
            expected_file = os.path.join(EXPECTED_OUTPUT_DIR, file)
            if not compare_files(test_file, expected_file):
                logging.error(f"File comparison failed for {file}")
        
        test_file = os.path.join(TEST_OUTPUT_DIR, MAIN_OUTPUT_FILE)
        expected_file = os.path.join(EXPECTED_OUTPUT_DIR, MAIN_OUTPUT_FILE)
        msg = f"Prediction mismatch. \nTest file: {test_file} \nExpected file: {expected_file}"
        self.assertTrue(compare_files(test_file, expected_file))
        

if __name__ == '__main__':
    unittest.main()