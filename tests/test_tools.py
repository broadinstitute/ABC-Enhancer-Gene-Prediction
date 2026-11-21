import os
import sys
import tempfile
import unittest

SCRIPTS_DIR = os.path.abspath("workflow/scripts")
sys.path.insert(0, SCRIPTS_DIR)
from tools import run_piped_commands_safe


class TestTools(unittest.TestCase):
    def test_run_piped_commands_safe(self) -> None:
        """Test piped commands with gzip (mimics getVariantOverlap usage)."""
        with tempfile.NamedTemporaryFile(delete=False, suffix=".gz") as f:
            outfile = f.name

        try:
            run_piped_commands_safe([["echo", "header"], ["gzip"]], outfile=outfile)
            run_piped_commands_safe([["echo", "data"], ["gzip"]], outfile=outfile, append=True)

            result = run_piped_commands_safe([["zcat", outfile]])
            self.assertEqual(result.strip(), b"header\ndata")
        finally:
            os.unlink(outfile)


if __name__ == "__main__":
    unittest.main()
