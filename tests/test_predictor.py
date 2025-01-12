import os
import sys
import unittest

SCRIPTS_DIR = os.path.abspath("workflow/scripts")
sys.path.insert(0, SCRIPTS_DIR)
from predictor import add_hic_from_hic_file
import pandas as pd

HIC_FILE = "https://encode-public.s3.amazonaws.com/2022/05/15/0571c671-3645-4f92-beae-51dfd3f42c36/ENCFF621AIY.hic"

# this file has 3k rows of E-G pairs with valid contact values
# contact values were generated from the original doubly stochastic method
# w/ juicebox data. prediction df was saved in add_hic_from_directory, prior
# to the columns getting dropped
# pred = pred[pred['hic_contact'] > 0].sample(n=5000)
# pred.to_csv('tests/test_data/test_pred_df.tsv', sep='\t', index=False)
PRED_DF_FILE = "tests/test_data/test_pred_df.tsv"


class TestPredictor(unittest.TestCase):
    def test_adding_hic(self) -> None:
        """
        Verifies hic streaming fills in hic contact the same way as juicebox
        - normalizing hic matrix to be doubly stochastic
        - imputing diagonals of matrix based on neighboring bins
        """
        pred = pd.read_csv(PRED_DF_FILE, sep="\t")
        expected_hic_contact_values = pred["hic_contact"]
        pred.drop("hic_contact", axis=1, inplace=True)
        pred = add_hic_from_hic_file(pred, HIC_FILE, "chr22", 5000)
        test_hic_contact_values = pred["hic_contact"]
        pd.testing.assert_series_equal(
            test_hic_contact_values, expected_hic_contact_values, atol=1e-9, rtol=0
        )


if __name__ == "__main__":
    unittest.main()
