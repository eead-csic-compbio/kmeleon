#!/usr/bin/env python
# -*- coding: utf-8 -*-

# CPCantalapiedra 2017

import sys, gzip
from optparse import OptionParser

################ MAIN
################

################ GLOBALS AND PARAMETERS

KMERS_FILE_FIELD_TARGET = 0
KMERS_FILE_FIELD_POSITION = 1
KMERS_FILE_FIELD_KMER = 2
KMERS_FILE_FIELD_DEPTH = 3

DEFAULT_DEPTH_PARAM = 0

## Usage
__usage = "usage: kmeleon_count.py [OPTIONS] KMERS_FILE\n"+\
          "The KMERS_FILE has the format of the output of kmeleon_extract.py\n"+\
          "Note that this software outputs to stderr and stdout.\n\n"+\
          "typical command: kmeleon_count.py -d 4 demo_data/demo.2_19.kmers > demo_data/demo.2_19.counts"
optParser = OptionParser(__usage)

optParser.add_option('-d', '--depth', action='store', dest='depth_param', type='int', \
                    help='The minimum times a k-mer is found to be reported in the output and counted. '+\
                    '(default: '+str(DEFAULT_DEPTH_PARAM)+'). '+\
                    ' This is indeed the same parameter as in kmeleon_extract.py.')

########### Read parameters
###########

(options, arguments) = optParser.parse_args()

sys.stderr.write("Running kmeleon count...\n")

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

refs_dict = {}
refs_list = []

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
    
    md_z_count = int(line_data[KMERS_FILE_FIELD_DEPTH])
    if md_z_count < depth_param: continue
    
    target = line_data[KMERS_FILE_FIELD_TARGET]
    pos = int(line_data[KMERS_FILE_FIELD_POSITION])
    md_z = line_data[KMERS_FILE_FIELD_KMER]
    
    if target in refs_dict:
        pos_dict = refs_dict[target]["pos_dict"]
        pos_list = refs_dict[target]["pos_list"]
    else:
        pos_dict = {}
        pos_list = []
        refs_dict[target] = {"pos_dict":pos_dict, "pos_list":pos_list}
        refs_list.append(target)
    
    if pos in pos_dict:
        #pos_kmer_count = pos_dict[pos]
        #pos_kmer_count+=1#.append(md_z)
        pos_dict[pos]["kc"]+=1
        pos_dict[pos]["dp"]+=md_z_count
    else:
        pos_dict[pos] = {"dp":md_z_count, "kc":1} #kc = k-mer count; dp = depth
        pos_list.append(pos)

kmers_fileobj.close()

sys.stdout.write("@Target\tPosition\tkc\tdp\n")
for target in refs_list:
    pos_list = refs_dict[target]["pos_list"]
    pos_dict = refs_dict[target]["pos_dict"]
    
    for pos in sorted(pos_list):
        if pos in pos_dict:
            pos_kmer_count = pos_dict[pos]["kc"]
            pos_depth = pos_dict[pos]["dp"]
            sys.stdout.write(str(target)+"\t"+str(pos)+"\t"+str(pos_kmer_count)+"\t"+str(pos_depth)+"\n")
        else:
            raise Exception("Position "+str(pos)+" is not in dict")

sys.stderr.write("Finished.\n")

## END
