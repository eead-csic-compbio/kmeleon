#!/usr/bin/env python
# -*- coding: utf-8 -*-

# CPCantalapiedra 2017

import sys, gzip
from optparse import OptionParser

from src.counts.Count import *
from src.intervals.Interval import *

################ MAIN
################

################ GLOBALS AND PARAMETERS

#COUNTS_FILE_FIELD_TARGET = 0
#COUNTS_FILE_FIELD_POSITION = 1
#COUNTS_FILE_FIELD_COUNT = 2

SAMPLES_FILE_FIELD_SAMPLE = 0
SAMPLES_FILE_FIELD_PATH = 1

DEFAULT_DIR_PARAM = "./"

## Usage
__usage = "Usage: kmeleon_table.py [OPTIONS] SAMPLES_FILE\n\n"+\
          "The SAMPLES_FILE has a row for each sample, with tab-separated columns: "+\
          "1st column is the sample name, 2nd column is the path to the counts or intervals file, "+\
          "which has the format of the output of \"kmeleon_count.py\" or \"kmeleon_intervals -m windows\", respectively.\n"+\
          "Note that this software outputs to stderr and stdout.\n\n"+\
          "typical command: kmeleon_table.py -D demo_data demo_data/demo_counts.list > demo_data/demo_counts.table\n"+\
          "typical command: kmeleon_table.py -w -D demo_data demo_data/demo_windows.list > demo_data/demo_windows.table"
optParser = OptionParser(__usage)

optParser.add_option('-D', '--DIR', action='store', dest='dir_param', type='string', \
                    help='The directory where the files with counts or intervals are located. '+\
                    'Note that this is unnecessary if the directory has been included in every path within the SAMPLES_FILE '
                    '(default: '+str(DEFAULT_DIR_PARAM)+')')

optParser.add_option('-w', action='store_true', dest="windows_mode", \
                     help="If this flag is set, the input data should be intervals of fixed length and positions for all samples. "+\
                      'i.e. as those produced by "kmeleon_intervals.py -m windows". '+\
                      'If -w is not set, the input data expected will be kmer counts as those produced with "kmeleon_count.py".')

###########
########### Functions

############### Functions for
############### COUNTS
def f_process_counts_files(samples_list, samples_fn_dict):
    refs_dict = {}
    refs_list = []
    pos_dict = {}
    pos_list = []
    for curr_sample in samples_list:
        sample_fn = samples_fn_dict[curr_sample]
        sys.stderr.write("Reading "+sample_fn+"\n")
        
        # Open .gz or regular text file
        ext_pos = sample_fn.rfind(".")
        ext = sample_fn[ext_pos+1:]
        if ext == "gz" or ext == "gzip":
            sample_file = gzip.open(sample_fn, 'r')
        else:
            sample_file = open(sample_fn, 'r')
            
        # Read every row in the file of counts of the current sample
        for i, line in enumerate(sample_file):
            
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
        
        sample_file.close()
    
    sys.stderr.write("All sample counts files processed.\n")
    
    return (refs_list, refs_dict, pos_list, pos_dict)

def f_print_counts_header(samples_list):
    # Create output header
    
    sys.stdout.write("@Target\tPosition\t")
    for sample in samples_list:
        sys.stdout.write(sample+"\t")
    sys.stdout.write("\n")
    
    return
    
def f_print_counts_rows(refs_list, refs_dict, pos_list, pos_dict):
    # Create output: rows
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
        
    return


############### Functions for
############### WINDOWS
def f_process_windows_files(samples_list, samples_fn_dict):
    refs_dict = {}
    refs_list = []
    pos_dict = {}
    pos_list = []
    for curr_sample in samples_list:
        sample_fn = samples_fn_dict[curr_sample]
        sys.stderr.write("Reading "+sample_fn+"\n")
        
        # Open .gz or regular text file
        ext_pos = sample_fn.rfind(".")
        ext = sample_fn[ext_pos+1:]
        if ext == "gz" or ext == "gzip":
            sample_file = gzip.open(sample_fn, 'r')
        else:
            sample_file = open(sample_fn, 'r')
            
        # Read every row in the file of counts of the current sample
        for i, line in enumerate(sample_file):
            
            if i==0: continue # skip header
            
            line_data = line.strip().split("\t")
            
            curr_target = line_data[INTERVALS_FILE_FIELD_TARGET]
            start_pos = long(line_data[INTERVALS_FILE_FIELD_START])
            end_pos = long(line_data[INTERVALS_FILE_FIELD_END])
            kmer_count = float(line_data[INTERVALS_FILE_FIELD_COUNT])
            
            pos_key=str(start_pos)+":"+str(end_pos)
            
            if curr_target in refs_dict:
                pos_dict = refs_dict[curr_target]["pos_dict"]
                pos_list = refs_dict[curr_target]["pos_list"]
            else:
                pos_dict = {}
                pos_list = []
                refs_dict[curr_target] = {"pos_dict":pos_dict, "pos_list":pos_list}
                refs_list.append(curr_target)
            
            if pos_key in pos_dict:
                pos_sample_dict = pos_dict[pos_key]
                if curr_sample in pos_sample_dict:
                    raise Exception("Already in pos_sample_dict")
                else:
                    pos_sample_dict[curr_sample] = kmer_count
            else:
                pos_sample_dict = {}
                pos_dict[pos_key] = pos_sample_dict
                pos_sample_dict[curr_sample] = kmer_count
                pos_list.append(pos_key)
                
            #if i>=10000: break
        
        sample_file.close()
    
    sys.stderr.write("All sample counts files processed.\n")
    
    return (refs_list, refs_dict, pos_list, pos_dict)

def f_print_windows_header(samples_list):
    # Create output header
    
    sys.stdout.write("@Target\tStart\tEnd\t")
    for sample in samples_list:
        sys.stdout.write(sample+"\t")
    sys.stdout.write("\n")
    
    return

def f_print_windows_rows(refs_list, refs_dict, pos_list, pos_dict):
    # Create output: rows
    for target in refs_list:
        pos_dict = refs_dict[target]["pos_dict"]
        pos_list = refs_dict[target]["pos_list"]
        
        for pos in sorted(pos_list):
            splitpos = pos.split(":")
            start_pos = splitpos[0]
            end_pos = splitpos[1]
            
            sys.stdout.write(str(target)+"\t"+str(start_pos)+"\t"+str(end_pos)+"\t")
            
            for sample in samples_list:
                if sample in pos_dict[pos]:
                    kmer_count = pos_dict[pos][sample]
                    sys.stdout.write(str(kmer_count)+"\t")
                else:
                    sys.stdout.write("0\t")
                
            sys.stdout.write("\n")
        
    return

###########
########### MAIN

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
    
## counts or intervals
if options.windows_mode:
    windows_mode = True
else:
    windows_mode = False

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
    sample_fn = dir_param+"/"+sample_file
    samples_fn_dict[sample_id] = sample_fn
    #sys.stderr.write("Added: "+sample_fn+"\n")

samples_list_file.close()

sys.stderr.write("Number of files to process: "+str(len(samples_list))+"\n")

# WINDOWS
if windows_mode:
    (refs_list, refs_dict, pos_list, pos_dict) = f_process_windows_files(samples_list, samples_fn_dict)
    
    f_print_windows_header(samples_list)
    f_print_windows_rows(refs_list, refs_dict, pos_list, pos_dict)
    
# COUNTS
else:
    (refs_list, refs_dict, pos_list, pos_dict) = f_process_counts_files(samples_list, samples_fn_dict)
    
    f_print_counts_header(samples_list)
    f_print_counts_rows(refs_list, refs_dict, pos_list, pos_dict)

sys.stderr.write("Finished kmeleon table.\n")

## END