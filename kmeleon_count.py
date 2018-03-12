#!/usr/bin/env python
# -*- coding: utf-8 -*-

# CPCantalapiedra 2017

import sys, gzip

sys.stderr.write("Running kmeleon count...\n")

if (len(sys.argv)==3):
    kmers_file = sys.argv[1]
    min_count = int(sys.argv[2])#4
else:
    sys.stderr.write("\n")
    sys.stderr.write("Wrong parameter number for kmeleon count.\n")
    sys.stderr.write("A typical kmeleon count command would be:\n")
    sys.stderr.write("\tkmeleon_count.py kmers_file 4\n")
    sys.stderr.write("\n")
    sys.stderr.write("We need 2 parameters:\n")
    sys.stderr.write("\t- kmers_file: a file with kmers by position \
(the one reported by kmeleon extract with the 'depths' option, for example). \
It could be in '.gz' extension or not.\n")
    sys.stderr.write("\t- min_depth: a number which is the minimum depth of a kmer to be considered and counted. For instance 4.\n")
    sys.stderr.write("\n")
    sys.exit(-1)


pos_dict = {}
pos_list = []

ext = kmers_file[-2:]
if ext == "gz":
    kmers_fileobj = gzip.open(kmers_file, 'r')
else:
    kmers_fileobj = open(kmers_file, 'r')
    
for i, line in enumerate(kmers_fileobj):
    
    if i==0: continue #header
    
    #if i>1000000: break
    line_data = line.strip().split("\t")
    
    #sys.stderr.write(str(line_data)+"\n")
    
    md_z_count = int(line_data[2])
    if md_z_count < min_count: continue
    
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
