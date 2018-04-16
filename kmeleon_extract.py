#!/usr/bin/env python
# -*- coding: utf-8 -*-

## CPCantalapiedra 2017 from a code by
# CPCantalapiedra 2015
######################################
# This is a script which takes a BAM file as input
# and the user has to specify also the contig/chr
# to be processed and a window size.
# From those parameters, the script reports, for each
# mapped position to the contig/chr of the BAM file,
# all the different kmers (MD:Z formatted in fact)
# of size equal to the window size specified.
######################################

import sys
import pysam
from optparse import OptionParser

############## FUNCTION DEFINITIONS
##############

# Open the file to be read (either as SAM or BAM)
def get_input_file(input_file_param, file_type):
    ret_file = None
    if file_type == FILE_TYPE_SAM:
        ret_file = pysam.AlignmentFile(input_file_param, "r")
    elif file_type == FILE_TYPE_BAM:
        ret_file = pysam.AlignmentFile(input_file_param, "rb")
    else:
        raise Exception("Unrecognized file type "+str(file_type)+".")
    
    return ret_file

# Count reads in a mappings file
def get_samfile_count(input_file, target, start, end):
    if target == DEFAULT_TARGET_PARAM:
        ret_count = input_file.count()
    elif start == DEFAULT_START_PARAM or end == DEFAULT_END_PARAM:
        ret_count = input_file.count(target)
    else:
        ret_count = input_file.count(target, start, end)
    
    return ret_count

# Obtain the iterator depending on target or not
def get_samfile_iter(input_file, target, start, end):
    if target == DEFAULT_TARGET_PARAM:
        ret_iter = input_file.fetch() #(until_eof=True)
    elif start == DEFAULT_START_PARAM or end == DEFAULT_END_PARAM:
        ret_iter = input_file.fetch(target)
    else:
        ret_iter = input_file.fetch(target, start, end)#samfile.fetch('chr1', 100, 120):
    
    return ret_iter

# Ask whether a value is a number or not
def is_number(s):
    try:
        int(s)
        return True
    except Exception:
        return False

# A function which expands the MD:Z field
# to single values. E.g.:
# 3G2A --> [1, 1, 1, G, 1, 1, A]
def f_expand_md_z(md_z_field):
    expanded_md_z = []
    
    prev_is_number = False
    curr_number = -1
    for md_z_data in md_z_field:
        if is_number(md_z_data):
            if prev_is_number:
                curr_number = curr_number*10 + int(md_z_data)
            else:
                curr_number = int(md_z_data)
            
            prev_is_number = True
        else:
            if prev_is_number:
                expanded_md_z.extend([1]*int(curr_number))
            
            expanded_md_z.append(md_z_data)
            
            prev_is_number = False
            curr_number = -1
    
    if prev_is_number:
        expanded_md_z.extend([1]*int(curr_number))
    
    #sys.stderr.write("Expanded MD:Z "+str(expanded_md_z)+"\n")
    
    return expanded_md_z

# A function which collapses an array of single values,
# obtained previously with f_expand_md_z function,
# to a MD:Z formatted field. E.g.:
# [1, 1, 1, G, 1, 1, A] --> 3G2A
def f_collapse_md_z(md_z):
    # NOTE: md_z must be an expanded_md_z
    
    collapsed_md_z = []
    sum_md_z_number = 0
    for md_z_data in md_z:
        if is_number(md_z_data):
            md_z_number = int(md_z_data)
            if md_z_number != 1:
                raise Exception("Is really an expanded MD:Z "+str(md_z)+"?")
            else:
                sum_md_z_number += 1
        else:
            # Add the previous numeric field if existed
            if sum_md_z_number > 0:
                collapsed_md_z.append(sum_md_z_number)
                sum_md_z_number = 0
            
            collapsed_md_z.append(md_z_data)
    
    # Add the last numeric field if necessary
    if sum_md_z_number > 0:
        collapsed_md_z.append(sum_md_z_number)
    
    #sys.stdout.write("Collapsed\n")
    #sys.stdout.write(str(collapsed_md_z)+"\n")
    
    return collapsed_md_z

# A function which clips an MD:Z field (md_z parameter),
# previously expanded with the f_expand_md_z function,
# based on a central position (pos parameter)
# and a half-window size (window_side parameter). E.g.:
# md_z = [1, 1, 1, G, 1, 1, A]
# pos = 3
# window_side = 1
# returns [1, 1, G]
def f_clip_md_z(md_z, pos, window_side):
    # NOTE: md_z is expected to be an expanded_md_z
    clipped_md_z = []
    
    start_clip = pos - window_side - 1
    end_clip = pos + window_side
    
    clipped_md_z = md_z[start_clip:end_clip]
    
    #sys.stdout.write("Clipped in "+str(pos)+", interval ["+str(start_clip)+", "+str(end_clip)+")\n")
    #sys.stdout.write(str(clipped_md_z)+"\n")
    
    return clipped_md_z

# A function which takes a read mapping position (read_map_position),
# a relative position of a kmer within the read (internal_pos parameter)
# and a left clipping based on soft clipping of the CIGAR (read_clipping parameter)
# to return the absolute position of a kmer
def f_get_kmer_contig_position(internal_pos, read_map_pos, read_clipping):
    kmer_pos = None
    
    # consider only left clipping (position 0 of the tuple)
    #kmer_pos = read_map_pos + (internal_pos - 1) + read_clipping[0]
    kmer_pos = read_map_pos + (internal_pos - 1) + read_clipping
    
    #sys.stdout.write("\tkmer_pos: "+str(kmer_pos)+"\n")
    
    return kmer_pos

# Adds a kmer to an absolute position,
# and counts how many times such kmer has been
# added to such position (depth of kmer at that position)
def f_add_kmer_depths(kmer_ref_id, kmer_ref_name, kmer_position, md_z):
    
    #sys.stdout.write(str(kmer_position)+" - "+str(md_z)+"\n")
    
    md_z_string = "".join([str(x) for x in md_z])
    
    # Obtain positions of the given reference_id
    if kmer_ref_id in refs_dict:
        pos_dict = refs_dict[kmer_ref_id]["pos_dict"]
        pos_list = refs_dict[kmer_ref_id]["pos_list"]
    else:
        pos_dict = {}
        pos_list = []
        refs_dict[kmer_ref_id] = {"pos_dict":pos_dict, "pos_list":pos_list, "ref_name":kmer_ref_name}
        # add new reference to list
        refs_list.append(kmer_ref_id)
    
    # add the kmer or increase its count
    if kmer_position in pos_dict:
        pos_kmers_dict = pos_dict[kmer_position]
        if md_z_string in pos_kmers_dict:
            pos_kmers_dict[md_z_string] += 1
        else:
            pos_kmers_dict[md_z_string] = 1
        
    else:
        pos_kmers_dict = {}
        pos_kmers_dict[md_z_string] = 1
        pos_dict[kmer_position] = pos_kmers_dict
        # additionally, add position to list
        pos_list.append(kmer_position)
    
    return

# A function which adds the insertions to the MD:Z field
# NOTE: the MD:Z field specifies deletions (insertions in the reference)
# but no insertions (deletions in the reference).
# However, the insertions are specified in the CIGAR as the length of the insertion and "I"
# Thus, we add the insertions to the expanded md z so that the kmers with differences due
# to insertions could be identified.
# This is not the case of deletions (insertions in the reference), which are indicated
# in the MD:Z field as "^" followed by the nucleotides missing in the read. Note that
# this string could be left untouched to differenciate kmers with different deletions.
# However, processing deletions further than that should be done with care, since
# a deletion followed by SNPs could not be unambiguously resolved from MD:Z of some mappers.
# E.g. "^AAAGG" is the last G the last nt of the deletion or a G SNP...
# Some mapper add a 0 between the indel and a SNP, but this is not always the case and it is
# not done when the sequence following the read is not an SNP (is a match), and therefore is
# difficult to process
def f_add_insertions(expanded_md_z, insertions):
    res_md = expanded_md_z # final res_md result of adding insertions
    
    for position in insertions:
        res_md = res_md[:position]+["I"]+res_md[position+1:]
    
    return res_md

## Prints all the kmers in a specific position
def f_print_pos_kmers(kmer_ref_name, kmer_position, pos_kmers_dict, min_depth):
    
    for kmer in pos_kmers_dict:
        kmer_count = pos_kmers_dict[kmer]
        if kmer_count >= min_depth:
            sys.stdout.write(str(kmer_ref_name)+"\t"+str(kmer_position)+"\t"+str(kmer)+"\t"+str(kmer_count)+"\n")
        
    return

##
def f_print_previous_kmers(curr_ref_id, curr_pos, min_depth, flush_interval):
    new_refs_list = []
    new_refs_dict = {}
    
    ################ Generate output
    # Rows
    for read_ref_id in refs_list:
        pos_list = refs_dict[read_ref_id]["pos_list"]
        pos_dict = refs_dict[read_ref_id]["pos_dict"]
        read_ref_name = refs_dict[read_ref_id]["ref_name"]
        
        if read_ref_id in new_refs_dict:
            new_pos_list = new_refs_dict[read_ref_id]["pos_list"]
            new_pos_dict = new_refs_dict[read_ref_id]["pos_dict"]
        else:
            new_pos_list = []
            new_pos_dict = {}
            new_refs_dict[read_ref_id] = {"pos_dict":new_pos_dict, "pos_list":new_pos_list, "ref_name":read_ref_name}
            new_refs_list.append(read_ref_id)
        
        for kmer_position in pos_list:
            if kmer_position in pos_dict:
                pos_kmers_dict = pos_dict[kmer_position]
                
                #sys.stderr.write("Position "+str(kmer_position)+"\n")
                if curr_ref_id != read_ref_id:
                    f_print_pos_kmers(read_ref_name, kmer_position, pos_kmers_dict, min_depth)
                    
                elif curr_pos - flush_interval > kmer_position:
                    f_print_pos_kmers(read_ref_name, kmer_position, pos_kmers_dict, min_depth)
                else:
                    
                    # Obtain positions of the given reference_id
                    new_pos_list.append(kmer_position)
                    new_pos_dict[kmer_position] = pos_kmers_dict
                
            else:
                raise Exception("Position "+str(kmer_position)+" is not in dict")
    
    #print "AFTER FLUSH"
    #print curr_pos
    #print len(new_list)
    #print len(new_dict)
    #print pos_list[-1]
    #print pos_dict[pos_list[-1]]
    
    return (new_refs_list, new_refs_dict)

################ MAIN
################

################ GLOBALS AND PARAMETERS

VERBOSE_ALL = 5
verbosity = 0

FLUSH_INTERVAL = 10000

CIGAR_INSERTION_FIELD = 1
CIGAR_SOFTCLIP_FIELD = 4

FILE_TYPE_BAM = "BAM"
FILE_TYPE_SAM = "SAM"

DEFAULT_TARGET_PARAM = "all"
DEFAULT_START_PARAM = -1
DEFAULT_END_PARAM = -1
DEFAULT_KMER_PARAM = 50
DEFAULT_DEPTH_PARAM = 0

## Usage
__usage = "usage: kmeleon_extract.py [OPTIONS] -b BAM_FILE|-s SAM_FILE\n"+\
          "Note that this software outputs to stderr and stdout.\n\n"+\
          "typical command: kmeleon_extract.py -d 4 -b demo_data/demo.2_19.bam > demo_data/demo.2_19.kmers"
optParser = OptionParser(__usage)

optParser.add_option('-t', '--target', action='store', dest='target_param', type='string', \
                    help='A chromosome number or name, or a specific contig, or "all" to process all the mappings. '+\
                    '(default: "all".)')

optParser.add_option('--start', action='store', dest='start_param', type='int', \
                    help='The -t target parameter is required. Starting basepairs '+\
                    'position to process within the given target. '+\
                    '(default: '+str(DEFAULT_START_PARAM)+')')

optParser.add_option('--end', action='store', dest='end_param', type='int', \
                    help='The -t target parameter is required. '+\
                    'Ending basepairs position to process within the given target. '+\
                    '(default: '+str(DEFAULT_END_PARAM)+')')

optParser.add_option('-k', '--kmer', action='store', dest='kmer_param', type='int', \
                    help='The length of fragments (k-mers) to parse. '+\
                    '(default: '+str(DEFAULT_KMER_PARAM)+')')

optParser.add_option('-d', '--depth', action='store', dest='depth_param', type='int', \
                    help='The minimum times a k-mer is found to be reported in the output. '+\
                    '(default: '+str(DEFAULT_DEPTH_PARAM)+')')

optParser.add_option('-b', '--bam', action='store', dest='bam_param', type='string', \
                    help='A BAM file to process. '+\
                    'Either the -b or the -s option is required.')

optParser.add_option('-s', '--sam', action='store', dest='sam_param', type='string', \
                    help='A SAM file to process. '+\
                    'Either the -s or the -b option is required.')

########### Read parameters
###########

(options, arguments) = optParser.parse_args()

sys.stderr.write("Running kmeleon extract...\n")

## Target
if options.target_param:
    target_param = options.target_param
else:
    target_param = DEFAULT_TARGET_PARAM

## Start
if options.start_param:
    if target_param == DEFAULT_TARGET_PARAM:
        optParser.exit(0, "The --start and --end options require a specific --target (-t) option.\n")
    else:
        start_param = options.start_param
else:
    start_param = DEFAULT_START_PARAM

## End
if options.end_param:
    if target_param == DEFAULT_TARGET_PARAM:
        optParser.exit(0, "The --start and --end options require a specific --target (-t) option.\n")
    else:
        end_param = options.end_param
else:
    end_param = DEFAULT_END_PARAM

## kmer
if options.kmer_param:
    kmer_param = options.kmer_param
else:
    kmer_param = DEFAULT_KMER_PARAM

## depth
if options.depth_param:
    depth_param = options.depth_param
else:
    depth_param = DEFAULT_DEPTH_PARAM
    
## BAM file
file_type = ""
if options.bam_param:
    input_file_param = options.bam_param
    file_type = FILE_TYPE_BAM
elif options.sam_param:
        input_file_param = options.sam_param
        file_type = FILE_TYPE_SAM
else:
    optParser.exit(0, "A BAM file or SAM file must be specified (with the -b or -s options, respectively).\n")

#########################
######################### BEGIN

#store_depths = True

sys.stderr.write("Reading file "+input_file_param+" to process target "+target_param+"\n")

#kmer_param = 50
window_side = kmer_param / 2

refs_dict = {} # chromosome, contig, ...
refs_list = [] # a list of chromosomes, contigs, ...
pos_dict = {} # a base pairs position
pos_list = [] # a list of base pairs positions

## TODO: Add exception handling to this
input_file = get_input_file(input_file_param, file_type)

numreads = get_samfile_count(input_file, target_param, start_param, end_param)

sys.stderr.write("Num of reads to process "+str(numreads)+"\n")

numreads_processed = 0
header_printed = False
first_to_flush_pos = -1

#for line in open(samfile, 'r'):
samfile_iter = get_samfile_iter(input_file, target_param, start_param, end_param)

for read in samfile_iter:
    
    #print read
    
    if read.is_unmapped: continue
    
    read_id = read.query_name
    # this is the chromosome or contig of the current mapping
    # however this is an id, to obtain the name we need to do:
    # input_file.getrname(reference_id)
    read_ref_id = read.reference_id
    read_ref_name = input_file.getrname(read_ref_id)
    
    try:
        read_md_z = read.get_tag("MD")
    except KeyError:
        if verbosity == VERBOSE_ALL:
            sys.stderr.write("WARNING: No MD:Z field in read "+read_id+"\n")
        #sys.stderr.write(str(read)+"\n")
        #sys.stderr.write("\n")
        #break
        continue
    
    expanded_md_z = f_expand_md_z(read_md_z)
    
    cigartuples = read.cigartuples
    
    if CIGAR_INSERTION_FIELD in [x[0] for x in cigartuples]:
        #print expanded_md_z
        #print read.cigarstring
        #print cigartuples
        position = 0
        insertions = []
        for cigartuple in cigartuples:
            if cigartuple[0]==CIGAR_SOFTCLIP_FIELD:
                continue
            if cigartuple[0]==CIGAR_INSERTION_FIELD:
                for i in range(1, cigartuple[1]+1):
                    insertions.append(position+i)
                
            position+=cigartuple[1]
            
        #print insertions
        
        expanded_md_z = f_add_insertions(expanded_md_z, insertions)
        #print expanded_md_z
    
    # obtain starting and ending position of read
    read_pos_start = 1
    read_pos_end = read.query_alignment_end - read.query_alignment_start
    ref_pos_start = read.reference_start+1
    if first_to_flush_pos == -1:
        first_to_flush_pos = ref_pos_start # a flag to trigger output process
    ref_pos_end = read.reference_end
    
    # we need to consider clipping when obtaining the final position of the
    # kmer in the contig
    read_left_clip = read.query_alignment_start
    
    #sys.stdout.write("\tRel Position: "+str(read_pos_start)+"-"+str(read_pos_end)+"\n")
    #sys.stdout.write("\tLeft clip (soft mask): "+str(read_left_clip)+"\n")
    #sys.stdout.write("\tAbs Position: "+str(ref_pos_start)+"-"+str(ref_pos_end)+"\n")
    
    pos_range = range(read_pos_start + window_side, read_pos_end - window_side + 1)
    if len(pos_range) > 0:
        for i in pos_range:
            
            clipped_md_z = f_clip_md_z(expanded_md_z, i, window_side)
            
            collapsed_md_z = f_collapse_md_z(clipped_md_z)
            
            #kmer_contig_position = f_get_kmer_contig_position(i, read_map_pos, read_clipping)
            kmer_contig_position = f_get_kmer_contig_position(i, ref_pos_start, read_left_clip)
            
            #if store_depths:
            f_add_kmer_depths(read_ref_id, read_ref_name, kmer_contig_position, collapsed_md_z)
            #else:
            #    f_add_kmer(kmer_contig_position, collapsed_md_z)
        
        ################ Generate output
        # Header
        if not header_printed:
            #if store_depths:
            sys.stdout.write("@Target\tPosition\tkmer(md_z)\tdp\n")
            #else:
            #    sys.stdout.write("@ Position\tMD_Z\n")
            header_printed = True
        
        if ref_pos_start - FLUSH_INTERVAL*2 > first_to_flush_pos:
            (refs_list, refs_dict) = f_print_previous_kmers(read_ref_id, ref_pos_start, depth_param, FLUSH_INTERVAL)
            #print "POS RANGE"
            #print read
            #print pos_range
            #print read_pos_start
            #print read_pos_end
            #print window_side
            first_to_flush_pos = refs_dict[refs_list[0]]["pos_list"][0]
    
    numreads_processed+=1
    if (numreads_processed % 1000 == 0):
        sys.stderr.write(str(numreads_processed)+" reads processed (out of "+str(numreads)+").\n")

input_file.close()

################ Generate output for the remaining data
# Rows
for read_ref_id in refs_list:
    pos_list = refs_dict[read_ref_id]["pos_list"]
    pos_dict = refs_dict[read_ref_id]["pos_dict"]
    read_ref_name = refs_dict[read_ref_id]["ref_name"]
    
    for kmer_position in pos_list:
        if kmer_position in pos_dict:
            #sys.stderr.write("Position "+str(kmer_position)+"\n")
            pos_kmers_dict = pos_dict[kmer_position]
            f_print_pos_kmers(read_ref_name, kmer_position, pos_kmers_dict, depth_param)
            
        else:
            raise Exception("Position "+str(kmer_position)+" is not in dict")
    
sys.stderr.write("Finished.\n")

## END
