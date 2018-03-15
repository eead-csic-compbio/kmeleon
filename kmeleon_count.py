#!/usr/bin/env python
# -*- coding: utf-8 -*-

# CPCantalapiedra 2017

import sys, gzip
from optparse import OptionParser

################ MAIN
################

sys.stderr.write("Running kmeleon count...\n")

################ GLOBALS AND PARAMETERS

DEFAULT_DEPTH_PARAM = 0

## Usage
__usage = "usage: kmeleon_count.py [OPTIONS] KMERS_FILE\n"+\
          "Note that this software outputs to stderr and stdout.\n\n"+\
          "typical command: kmeleon_count.py -d 4 demo_data/demo.2_19.kmers > demo_data/demo.2_19.counts"
optParser = OptionParser(__usage)

optParser.add_option('-d', '--depth', action='store', dest='depth_param', type='int', \
                    help='The minimum times a k-mer is found to be reported in the output.'+\
                    '(default: '+str(DEFAULT_DEPTH_PARAM)+')')

(options, arguments) = optParser.parse_args()

# THIS IS MANDATORY
if not arguments or len(arguments)==0:
    optParser.exit(0, "You may wish to run '-help' option.\n")
# INPUT FILE
kmers_file = arguments[0]

## depth
if options.depth_param:
    depth_param = options.depth_param
else:
    depth_param = DEFAULT_DEPTH_PARAM

#########################
######################### BEGIN

pos_dict = {}
pos_list = []

ext_pos = kmers_file.rfind(".")
ext = kmers_file[ext_pos+1:]

if ext == "gz" or ext == "gzip":
    kmers_fileobj = gzip.open(kmers_file, 'r')
else:
    kmers_fileobj = open(kmers_file, 'r')
    
for i, line in enumerate(kmers_fileobj):
    
    if i==0: continue #header
    
    #if i>1000000: break
    line_data = line.strip().split("\t")
    
    #sys.stderr.write(str(line_data)+"\n")
    
    md_z_count = int(line_data[2])
    if md_z_count < depth_param: continue
    
    pos = int(line_data[0])
    md_z = line_data[1]
    
    if pos in pos_dict:
        #pos_kmer_count = pos_dict[pos]
        #pos_kmer_count+=1#.append(md_z)
        pos_dict[pos]+=1
    else:
        pos_dict[pos] = 1#[md_z]
        pos_list.append(pos)

kmers_fileobj.close()

#sys.stdout.write("@ Position\tNumber_diff_kmers\n")
for pos in sorted(pos_list):
    if pos in pos_dict:
        pos_kmer_count = pos_dict[pos]
        sys.stdout.write(str(pos)+"\t"+str(pos_kmer_count)+"\n")
    else:
        raise Exception("Position "+str(pos)+" is not in dict")

sys.stderr.write("Finished.\n")

## END
