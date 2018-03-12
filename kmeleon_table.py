#!/usr/bin/env python
# -*- coding: utf-8 -*-

# CPCantalapiedra 2017

import sys, gzip

sys.stderr.write("Running kmeleon table...\n")

if (len(sys.argv)==4):
    samples_list_fn = sys.argv[1]
    chrom = sys.argv[2]
    counts_dir = sys.argv[3]
else:
    sys.stderr.write("\n")
    sys.stderr.write("Wrong parameter number for kmeleon table.\n")
    sys.stderr.write("A typical kmeleon table command would be:\n")
    sys.stderr.write("\tkmeleon_table.py samples_list_file chr1 counts\n")
    sys.stderr.write("\n")
    sys.stderr.write("We need 3 parameters:\n")
    sys.stderr.write("\t- samples_list_file: a file with a list of files whose counts will be joined in a single table.\n")
    sys.stderr.write("\t- target: chromosome number, for example.\n")
    sys.stderr.write("\t- counts directory: a path to the directory where the files with counts (from kmeleon count) can be found.\n")
    sys.stderr.write("\n")
    sys.exit(-1)

samples_list_file = open(samples_list_fn, 'r')

samples_list = []
samples_counts_fn_dict = {}

sys.stderr.write("Reading sample names file...\n")
for sample in samples_list_file:
    sample_id = sample.strip()
    samples_list.append(sample_id)
    sample_counts_fn = counts_dir+"/"+sample_id+"."+chrom+".counts.out.gz"
    samples_counts_fn_dict[sample_id] = sample_counts_fn
    #sys.stderr.write("Added: "+sample_counts_fn+"\n")

samples_list_file.close()

sys.stderr.write("Number of files to process: "+str(len(samples_list))+"\n")

pos_dict = {}
pos_list = []
for curr_sample in samples_list:
    sample_counts_fn = samples_counts_fn_dict[curr_sample]
    sys.stderr.write("Reading "+sample_counts_fn+"\n")
    sample_counts_file = gzip.open(sample_counts_fn, 'r')
    
    for i, line in enumerate(sample_counts_file):
        line_data = line.strip().split("\t")
        curr_pos = long(line_data[0])
        kmer_count = int(line_data[1])
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

sys.stdout.write("Pos\t")
for sample in samples_list:
    sys.stdout.write(sample+"\t")
sys.stdout.write("\n")

# rows

for pos in sorted(pos_list):#pos_dict:
    sys.stdout.write(str(pos)+"\t")
    for sample in samples_list:
        if sample in pos_dict[pos]:
            kmer_count = pos_dict[pos][sample]
            sys.stdout.write(str(kmer_count)+"\t")
        else:
            sys.stdout.write("0\t")
    sys.stdout.write("\n")

sys.stderr.write("Finished kmers_counts_join.py\n")

## END
