#!/usr/bin/env python
# -*- coding: utf-8 -*-

# CPCantalapiedra 2017

import sys, gzip
from optparse import OptionParser

################ MAIN
################

################ GLOBALS AND PARAMETERS

COUNTS_FILE_FIELD_TARGET = 0
COUNTS_FILE_FIELD_POSITION = 1
COUNTS_FILE_FIELD_COUNT = 2

SAMPLES_FILE_FIELD_SAMPLE = 0
SAMPLES_FILE_FIELD_PATH = 1

DEFAULT_DIR_PARAM = "./"

## Usage
__usage = "usage: kmeleon_table.py [OPTIONS] SAMPLES_COUNTS_FILE\n"+\
          "The SAMPLES_COUNTS_FILE has a row for each sample, with tab-separated columns:"+\
          "1st column is the sample name, 2nd column is the path to the counts file,"+\
          "which has the format of the output of kmeleon_count.py\n"+\
          "Note that this software outputs to stderr and stdout.\n\n"+\
          "typical command: kmeleon_table.py -D demo_data demo_data/demo_samples.list > demo_data/demo_counts.table"
optParser = OptionParser(__usage)

optParser.add_option('-D', '--DIR', action='store', dest='dir_param', type='string', \
                    help='The directory where the files with counts for the samples are located.'+\
                    'Note that this is unnecessary if the directory has been included in every path within the SAMPLES_COUNTS_FILES'
                    '(default: '+str(DEFAULT_DIR_PARAM)+')')

########### Read parameters
###########

(options, arguments) = optParser.parse_args()

sys.stderr.write("Running kmeleon table...\n")

# THIS IS MANDATORY
if not arguments or len(arguments)==0:
    optParser.exit(0, "You may wish to run '-help' option.\n")
# INPUT FILE
samples_counts_fn = arguments[0]

## depth
if options.dir_param:
    dir_param = options.dir_param
else:
    dir_param = DEFAULT_DIR_PARAM

#if (len(sys.argv)==4):
#    samples_list_fn = sys.argv[1]
#    chrom = sys.argv[2]
#    counts_dir = sys.argv[3]
#else:
#    sys.stderr.write("\n")
#    sys.stderr.write("Wrong parameter number for kmeleon table.\n")
#    sys.stderr.write("A typical kmeleon table command would be:\n")
#    sys.stderr.write("\tkmeleon_table.py samples_list_file chr1 counts\n")
#    sys.stderr.write("\n")
#    sys.stderr.write("We need 3 parameters:\n")
#    sys.stderr.write("\t- samples_list_file: a file with a list of files whose counts will be joined in a single table.\n")
#    sys.stderr.write("\t- target: chromosome number, for example.\n")
#    sys.stderr.write("\t- counts directory: a path to the directory where the files with counts (from kmeleon count) can be found.\n")
#    sys.stderr.write("\n")
#    sys.exit(-1)

samples_list_file = open(samples_counts_fn, 'r')

samples_list = []
samples_fn_dict = {}

sys.stderr.write("Reading sample names file...\n")
for sample in samples_list_file:
    sample_data = sample.strip().split("\t")
    
    sample_id = sample_data[SAMPLES_FILE_FIELD_SAMPLE]
    sample_file = sample_data[SAMPLES_FILE_FIELD_PATH]
    samples_list.append(sample_id)
    sample_counts_fn = dir_param+"/"+sample_file
    samples_fn_dict[sample_id] = sample_counts_fn
    #sys.stderr.write("Added: "+sample_counts_fn+"\n")

samples_list_file.close()

sys.stderr.write("Number of files to process: "+str(len(samples_list))+"\n")

refs_dict = {}
refs_list = []
pos_dict = {}
pos_list = []
for curr_sample in samples_list:
    sample_counts_fn = samples_fn_dict[curr_sample]
    sys.stderr.write("Reading "+sample_counts_fn+"\n")
    
    # Open .gz or regular text file
    ext_pos = sample_counts_fn.rfind(".")
    ext = sample_counts_fn[ext_pos+1:]
    if ext == "gz" or ext == "gzip":
        sample_counts_file = gzip.open(sample_counts_fn, 'r')
    else:
        sample_counts_file = open(sample_counts_fn, 'r')
    
    # Read every row in the file of counts of the current sample
    for i, line in enumerate(sample_counts_file):
        
        if i==0: continue # skip header
        
        line_data = line.strip().split("\t")
        
        curr_target = line_data[COUNTS_FILE_FIELD_TARGET]
        curr_pos = long(line_data[COUNTS_FILE_FIELD_POSITION])
        kmer_count = int(line_data[COUNTS_FILE_FIELD_COUNT])
        
        if curr_target in refs_dict:
            pos_dict = refs_dict[curr_target]["pos_dict"]
            pos_list = refs_dict[curr_target]["pos_list"]
        else:
            pos_dict = {}
            pos_list = []
            refs_dict[curr_target] = {"pos_dict":pos_dict, "pos_list":pos_list}
            refs_list.append(curr_target)
        
        if curr_pos in pos_dict:
            pos_sample_dict = pos_dict[curr_pos]
            if curr_sample in pos_sample_dict:
                raise Exception("Already in pos_sample_dict")
            else:
                pos_sample_dict[curr_sample] = kmer_count
        else:
            pos_sample_dict = {}
            pos_dict[curr_pos] = pos_sample_dict
            pos_sample_dict[curr_sample] = kmer_count
            pos_list.append(curr_pos)
            
        #if i>=10000: break
    
    
    sample_counts_file.close()

sys.stderr.write("All sample counts files processed.\n")

#print pos_dict

# Create output
# header

sys.stdout.write("@Target\tPosition\t")
for sample in samples_list:
    sys.stdout.write(sample+"\t")
sys.stdout.write("\n")

# rows
for target in refs_list:
    pos_dict = refs_dict[target]["pos_dict"]
    pos_list = refs_dict[target]["pos_list"]
    
    for pos in sorted(pos_list):
        sys.stdout.write(str(target)+"\t"+str(pos)+"\t")
        for sample in samples_list:
            if sample in pos_dict[pos]:
                kmer_count = pos_dict[pos][sample]
                sys.stdout.write(str(kmer_count)+"\t")
            else:
                sys.stdout.write("0\t")
            
        sys.stdout.write("\n")  

sys.stderr.write("Finished kmeleon table.\n")

## END
