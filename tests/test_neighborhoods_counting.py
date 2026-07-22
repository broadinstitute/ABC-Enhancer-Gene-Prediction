import os
import sys
import tempfile
import unittest

import pysam

SCRIPTS_DIR = os.path.abspath("workflow/scripts")
sys.path.insert(0, SCRIPTS_DIR)
from neighborhoods import count_bam

CHROM = "chr1"
CHROM_LEN_IN_BAM = 2000  # header contig length, deliberately larger than chrom_sizes_map's declared length
CHROM_SIZES_MAP = {CHROM: 1000}  # used to test right-edge clipping independent of the bam contig length
FLANK_BP = 50

# (position, length) for each synthetic single-end read
READS = {
    "R1": (420, 20),  # inside element1 [400,500)
    "R2": (460, 20),  # inside element1
    "R3": (360, 20),  # inside element1's left flank [350,400)
    "R4": (510, 20),  # inside element1's right flank [500,550)
    "R5": (700, 20),  # unrelated control read, should never be counted
    "R_boundary": (5, 10),  # inside element2's clipped left flank [0,20)
    "R6": (1010, 20),  # inside element3's *unclipped* right flank, must be excluded once clipped to chrom len 1000
    "R7": (950, 15),  # inside element3's left flank [920,970)
    "R8": (992, 5),  # inside element3's clipped right flank [990,1000)
}

# element rows: (chr, start, end)
ELEMENT1 = (CHROM, 400, 500)
ELEMENT2 = (CHROM, 20, 40)  # near chromosome start, tests left-edge clipping at 0
ELEMENT3 = (CHROM, 970, 990)  # near chrom_sizes_map's declared end, tests right-edge clipping


def _make_bam(path):
    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": CHROM_LEN_IN_BAM, "SN": CHROM}]}
    unsorted_path = path + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted_path, "wb", header=header) as outf:
        for name, (pos, length) in READS.items():
            a = pysam.AlignedSegment()
            a.query_name = name
            a.query_sequence = "A" * length
            a.flag = 0
            a.reference_id = 0
            a.reference_start = pos
            a.mapping_quality = 60
            a.cigar = [(0, length)]
            outf.write(a)
    pysam.sort("-o", path, unsorted_path)
    pysam.index(path)
    os.unlink(unsorted_path)


def _make_bed(path, rows):
    with open(path, "w") as f:
        for chrom, start, end in rows:
            f.write(f"{chrom}\t{start}\t{end}\n")


def _read_counts(path):
    with open(path) as f:
        return [int(line.strip().split("\t")[3]) for line in f]


class TestH3K27acCountingModes(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.bam_path = os.path.join(self.tmpdir, "test.bam")
        self.bed_path = os.path.join(self.tmpdir, "elements.bed")
        _make_bam(self.bam_path)
        _make_bed(self.bed_path, [ELEMENT1, ELEMENT2, ELEMENT3])

    def _count(self, count_mode):
        output = os.path.join(self.tmpdir, f"counts.{count_mode}.bed")
        count_bam(
            self.bam_path,
            self.bed_path,
            output,
            genome_sizes=None,
            count_mode=count_mode,
            flank_bp=FLANK_BP,
            chrom_sizes_map=CHROM_SIZES_MAP,
        )
        return _read_counts(output)

    def test_element_mode_matches_current_behavior(self) -> None:
        """element mode counts only within [start,end), unaffected by flank_bp."""
        self.assertEqual(self._count("element"), [2, 0, 0])

    def test_element_flank_mode_extends_symmetrically(self) -> None:
        """element_flank counts within [start-N,end+N), clipped to chrom bounds."""
        self.assertEqual(self._count("element_flank"), [4, 1, 2])

    def test_flank_only_mode_excludes_element(self) -> None:
        """flank_only sums the two flanking windows, excluding the element itself."""
        self.assertEqual(self._count("flank_only"), [2, 1, 2])

    def test_right_edge_clipping_uses_chrom_sizes_map(self) -> None:
        """
        Element3's right flank is clipped to chrom_sizes_map's declared length (1000),
        even though the bam's own contig is longer (2000). R6 at pos 1010 sits inside
        the *unclipped* window but must be excluded once clipping is applied.
        """
        element_flank_counts = self._count("element_flank")
        flank_only_counts = self._count("flank_only")
        # element3 is the 3rd row (index 2): expected count is 2 (R7 + R8), not 3 (R6 too)
        self.assertEqual(element_flank_counts[2], 2)
        self.assertEqual(flank_only_counts[2], 2)


if __name__ == "__main__":
    unittest.main()
