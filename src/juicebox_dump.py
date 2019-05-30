## Download and dump data from Aiden Lab juicebox, 
## matching the file formats and directory structure of the Rao et al. 2014 data from GEO
#
# use Java-1.8

import argparse
import subprocess

def parseargs():
    parser = argparse.ArgumentParser(description='Download and dump HiC data')
    parser.add_argument('--hic_file', required=True, help="Path or url to .hic file.")
    parser.add_argument('--juicebox', required=True, default="", help="path to juicebox executable or java command invoking juicer_tools.jar. eg: 'java -jar juicer_tools.jar'")
    parser.add_argument('--resolution', default=5000, help="Resolution of HiC to download. In units of bp.")
    parser.add_argument('--outdir', default=".")
    parser.add_argument('--obskr', action="store_true", help="Only download the KR observed matrix (as opposed to the Raw matrix and the KR norm vector separately")
    parser.add_argument('--chromosomes', default="all", help="comma delimited list of chromosomes to download")

    return parser.parse_args()

def main(args):

    if args.chromosomes == "all":
        chromosomes = list(range(1,23)) + ['X']
    else:
        chromosomes = args.chromosomes.split(",")

    for chromosome in chromosomes:
        print("Starting chr" + str(chromosome) + " ... ")
        outdir = "{0}/{2}kb_resolution_intrachromosomal/chr{1}/".format(args.outdir, chromosome, int(args.resolution/1000))
        command = "mkdir -p " + outdir
        out = subprocess.getoutput(command)

        if args.obskr:
	        ## Download observed matrix with KR normalization
	        command = args.juicebox + " dump observed KR {0} {1} {1} BP {3} {2}/chr{1}_5kb.KRobserved".format(args.hic_file, chromosome, outdir, args.resolution)
	        print(command)
	        out = subprocess.getoutput(command)
        else:
	        ## Download raw observed matrix
	        command = args.juicebox + " dump observed NONE {0} {1} {1} BP {3} {2}/chr{1}_5kb.RAWobserved".format(args.hic_file, chromosome, outdir, args.resolution)
	        print(command)
	        out = subprocess.getoutput(command)
	        
	        ## Download KR norm file
	        command = args.juicebox + " dump norm KR {0} {1} BP {3} {2}/chr{1}_5kb.KRnorm".format(args.hic_file, chromosome, outdir, args.resolution)
	        out = subprocess.getoutput(command)
	        print(out)

if __name__ == '__main__':
    args = parseargs()
    main(args)
