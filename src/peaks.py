import pandas as pd
import numpy as np
import os
import os.path
from tools import *
from neighborhoods import *

def run_macs2_callpeak(infile, outdir, outfile_prefix, MACS_format, pval_cutoff, genome_sizes, force=False):

	if force or (not os.path.exists(os.path.join(outdir, outfile_prefix + '_peaks.narrowPeak'))) or (os.path.getsize(os.path.join(outdir, outfile_prefix + '_peaks.narrowPeak')) == 0):
		#Macs does not support bgzipped files. Make a temp file
		remove_tmp = False
		if infile.endswith(".bgz"):
			remove_tmp = True
			tmpfile = infile.replace(".bgz", ".tmp.gz")
			run_command("bgzip -cd {} | gzip > {}".format(infile, tmpfile))
			infile = tmpfile

		macs_command = "macs2 callpeak -t {} -n {} -f {} -g hs -p {} --call-summits --outdir {}".format(infile, outfile_prefix, MACS_format, pval_cutoff, outdir)
		print("Running: {}".format(macs_command))
		p = Popen(macs_command, stdout=PIPE, stderr=PIPE, shell=True)
		(stdoutdata, stderrdata) = p.communicate()

		outfile = os.path.join(outdir, outfile_prefix  + '_peaks.narrowPeak')

		sort_command = "bedtools sort -faidx {} -i {} > {} ; mv {} {}; ".format(genome_sizes, outfile, outfile + '.tmp', outfile + '.tmp', outfile)
		print("Running: {}".format(sort_command))
		p = Popen(sort_command, stdout=PIPE, stderr=PIPE, shell=True)
		(stdoutdata, stderrdata) = p.communicate()

		if remove_tmp:
			run_command("rm {}".format(tmpfile))
	else:
		print('{}_peaks.narrowPeak already exists. Not recreating'.format(outfile_prefix))

# def make_v_plot(bamfile, bedfile, out, extension):
# 	print("Making v plot for file {}...".format(bamfile))
# 	command = reuse("Python-2.7") + "python /seq/lincRNA/Jesse/bin/scripts/pyMakeVplot_JN.py -a {} -b {} -c 4 -o {} -e {} -p ends -v -u".format(bamfile, bedfile, out, extension)
# 	print("Running: {}".format(command))

# 	p = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
# 	(stdoutdata, stderrdata) = p.communicate()

# 	err = str(stderrdata, 'utf-8')
# 	print(err)

# 	if "Error" not in err:
# 		data = pd.read_table(out, header=None)

# 		window = 20
# 		avg = np.convolve(data.ix[:,0], np.ones(window), 'same')/window/np.mean(data.ix[1:200,0]) #20bp moving average
# 		#tss_enrichment = avg[1000]
# 		tss_enrichment = np.amax(avg)

# 		print("Done computing v plot. TSS enrichment = {}".format(tss_enrichment))
# 	else:
# 		tss_enrichment = np.nan
# 		print("Making vplot failed. Error is: ".format(err))

# 	return(tss_enrichment)

def make_candidate_regions_from_summits(macs_peaks, accessibility_file, genome_sizes, regions_whitelist, regions_blacklist, n_enhancers, peak_extend, outdir):
    ## Generate enhancer regions from MACS summits
    # 1. Count reads in DHS peaks
    # 2. Take top N regions, get summits, extend summits, merge

    outfile = os.path.join(outdir, os.path.basename(macs_peaks) + ".candidateRegions.bed")
    raw_counts_outfile = os.path.join(outdir, os.path.basename(macs_peaks) + os.path.basename(accessibility_file) + ".Counts.bed")

    if regions_whitelist:
    	whitelist_command = "| (bedtools intersect -a {regions_whitelist} -b {genome_sizes}.bed -wa | bedtools slop -i stdin -b {peak_extend} -g {genome_sizes} | cut -f 1-3 && cat) |"
    else:
    	whitelist_command = ""

    if regions_blacklist:
    	blacklist_command = "bedtools intersect -v -wa -a stdin -b {regions_blacklist} | "
    else:
    	blacklist_command = ""

    #1. Count DHS/ATAC reads in candidate regions
    run_count_reads(accessibility_file, raw_counts_outfile, macs_peaks, genome_sizes)

    #2. Take top N regions, get summits, extend summits, merge, remove blacklist, add whitelist, sort and merge
    #use -sorted in intersect command? Not worth it, both files are small
    command = "bedtools sort -i {raw_counts_outfile} -faidx {genome_sizes} | bedtools merge -i stdin -c 4 -o max | sort -nr -k 4 | head -n {n_enhancers} |" + \
        "bedtools intersect -b stdin -a {macs_peaks} -wa |" + \
        "awk '{{print $1 "\t" $2 + $10 "\t" $2 + $10}}' |" + \
        "bedtools slop -i stdin -b {peak_extend} -g {genome_sizes} |" + \
        "bedtools sort -i stdin -faidx {genome_sizes} |" + \
        "bedtools merge -i stdin | " + \
        blacklist_command + \
        "cut -f 1-3 " + whitelist_command + \
        "bedtools sort -i stdin -faidx {genome_sizes} | bedtools merge -i stdin > {outfile}"

    command = command.format(**locals())
    print("Running: " + command)
    out = getoutput(command)
    print(out)

    return out

def get_macs_format(feature_type):
	if feature_type == "ATAC":
		return "BAMPE"
	elif feature_type == "DHS":
		return "AUTO"
	else:
		print("Unknown feature type")
		return None

def compute_macs_stats(peak_file):
	npeaks = int(check_output("cut -f1-3 {} | uniq | wc -l".format(peak_file), shell=True))

	return npeaks
