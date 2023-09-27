from hic import load_hic_avg

data = load_hic_avg(
    "/oak/stanford/groups/engreitz/Projects/ABC/HiC/avg_track2/chr1/chr1.avg.gz", 5000
)
import time

time.sleep(10)
