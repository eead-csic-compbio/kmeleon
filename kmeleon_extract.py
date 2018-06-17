#!/usr/bin/env python
# -*- coding: utf-8 -*-

## CPCantalapiedra 2018 from a code by
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

# A function which takes a read mapping position (read_map_position),
# a relative position of a kmer within the read (internal_pos parameter)
# and a left clipping based on soft clipping of the CIGAR (read_clipping parameter)
# to return the absolute position of a kmer
def f_get_kmer_contig_position(internal_pos, read_map_pos, read_clipping):
    kmer_pos = None
    
    # consider only left clipping (position 0 of the tuple)
    kmer_pos = read_map_pos + (internal_pos - 1)
    #kmer_pos = read_map_pos + (internal_pos - 1) + read_clipping
    
    #sys.stdout.write("\tkmer_pos: "+str(kmer_pos)+"\n")
    
    return kmer_pos

# Adds a kmer to an absolute position,
# and counts how many times such kmer has been
# added to such position (depth of kmer at that position)
def f_add_kmer_depths(refs_dict, refs_list, kmer_ref_id, kmer_ref_name, kmer_position, md_z, kmer_seq):
    
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
            pos_kmers_dict[md_z_string]["count"] += 1
        else:
            pos_kmers_dict[md_z_string] = {"count":1, "kmer_seq":kmer_seq}
        
    else:
        pos_kmers_dict = {}
        pos_kmers_dict[md_z_string] = {"count":1, "kmer_seq":kmer_seq}
        pos_dict[kmer_position] = pos_kmers_dict
        # additionally, add position to list
        pos_list.append(kmer_position)
    
    return

# A function which processes the insertions and deletions to the MD:Z field
def f_add_insertions(md_z, cigartuples, read_seq, read_left_clip):
        ret_read_seq = read_seq
        ret_md_z = md_z
        
        insertions = []
        if CIGAR_FIELD_INSERTION in [x[0] for x in cigartuples] or \
            CIGAR_FIELD_DELETION in [x[0] for x in cigartuples]:
            position = 0
            insertions = []
            
            #M BAM_CMATCH 0
            #I BAM_CINS 1
            #D BAM_CDEL 2
            #N BAM_CREF_SKIP 3
            #S BAM_CSOFT_CLIP 4
            #H BAM_CHARD_CLIP 5
            #P BAM_CPAD 6
            #= BAM_CEQUAL 7
            #X BAM_CDIFF 8
            #B BAM_CBACK 9
            
            for cigartuple in cigartuples:
                cigarfield = cigartuple[0]
                cigarspan = cigartuple[1]
                
                if cigarfield == CIGAR_FIELD_SOFT_CLIP or \
                    cigarfield == CIGAR_FIELD_HARD_CLIP:
                    continue
                
                elif cigarfield == CIGAR_FIELD_INSERTION:
                    #for i in range(1, cigarspan+1):
                    #    insertions.append(position+i)
                    ret_md_z = ret_md_z[:position]+[SYMBOL_INSERTION]*cigarspan+ret_md_z[position:]
                    
                    position+=cigarspan
                    
                elif cigarfield == CIGAR_FIELD_DELETION:
                    ret_md_z = ret_md_z[:position]+ret_md_z[position+cigarspan+1:]
                    
                    ret_md_z = ret_md_z[:position]+[SYMBOL_DELETION]*cigarspan+ret_md_z[position:]
                    
                    read_pos = position+read_left_clip
                    ret_read_seq = ret_read_seq[:read_pos]+"".join([SYMBOL_DELETION]*cigarspan)+ret_read_seq[read_pos:]
                    
                    position+=cigarspan
                    
                    #ret_md_z = ret_md_z[:position]+ret_md_z[position+cigarspan+1:]
                    #position+=1
                    
                elif cigarfield == CIGAR_FIELD_MATCH:
                    position+=cigarspan
                    
                else:
                    sys.stderr.write("**WARNING: Skipping cigarfield "+str(cigarfield)+"\n")
                    position+=cigarspan
                
            #for position in insertions:
            #    ret_md_z = ret_md_z[:position-1]+["I"]+ret_md_z[position-1:]
            #print expanded_md_z
        
        return (ret_md_z, ret_read_seq)

## Prints all the kmers in a specific position
def f_print_pos_kmers(kmer_ref_name, kmer_position, pos_kmers_dict, min_depth):
    
    for kmer in pos_kmers_dict:
        kmer_count = pos_kmers_dict[kmer]["count"]
        kmer_seq = pos_kmers_dict[kmer]["kmer_seq"]
        if kmer_count >= min_depth:
            if read_param:
                sys.stdout.write(str(kmer_ref_name)+"\t"+str(kmer_position)+"\t"+
                                 str(kmer)+"\t"+kmer_seq+"\t"+str(kmer_count)+"\n")
            else:
                sys.stdout.write(str(kmer_ref_name)+"\t"+str(kmer_position)+"\t"+
                                 str(kmer)+"\t"+str(kmer_count)+"\n")
        
    return

## Prints all the kmers from already-processed reads (by position, according to flush_interval)
def f_print_previous_kmers(refs_dict, refs_list, curr_ref_id, curr_pos, min_depth, flush_interval):
    new_refs_list = []
    new_refs_dict = {}
    
    ################ Generate output
    # Rows
    for read_ref_id in refs_list:
        
        pos_list = refs_dict[read_ref_id]["pos_list"]
        pos_dict = refs_dict[read_ref_id]["pos_dict"]
        read_ref_name = refs_dict[read_ref_id]["ref_name"]
        
        # Obtain or create the info for this reference,
        # which is common for all the kmers in such reference
        if read_ref_id in new_refs_dict:
            new_pos_list = new_refs_dict[read_ref_id]["pos_list"]
            new_pos_dict = new_refs_dict[read_ref_id]["pos_dict"]
        else:
            new_pos_list = []
            new_pos_dict = {}
            new_refs_dict[read_ref_id] = {"pos_dict":new_pos_dict,
                                          "pos_list":new_pos_list,
                                          "ref_name":read_ref_name}
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
                    
                    #if read_ref_id in new_refs_dict:
                    #    new_pos_list = new_refs_dict[read_ref_id]["pos_list"]
                    #    new_pos_dict = new_refs_dict[read_ref_id]["pos_dict"]
                    #else:
                    #    new_pos_list = []
                    #    new_pos_dict = {}
                    #    new_refs_dict[read_ref_id] = {"pos_dict":new_pos_dict,
                    #                                  "pos_list":new_pos_list,
                    #                                  "ref_name":read_ref_name}
                    #    new_refs_list.append(read_ref_id)
                    #    
                    # Obtain positions of the given reference_id
                    new_pos_list.append(kmer_position)
                    new_pos_dict[kmer_position] = pos_kmers_dict
                    
                    #if read_ref_id not in new_refs_dict:
                    #    new_refs_dict[read_ref_id] = {"pos_dict":new_pos_dict,
                    #                                  "pos_list":new_pos_list,
                    #                                  "ref_name":read_ref_name}
                    #    new_refs_list.append(read_ref_id)
                
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

VERBOSE_ALL = 5 # maximum verbosity
verbosity = 0 # verbosity parameter

# FLUSH_INTERVAL controls the distance to
# already-processed reads whose kmers will
# be output while processing resumes with
# the reads in the current position.
FLUSH_INTERVAL = 10000

# The next are the fields found in the
# cigar tuples from pysam
CIGAR_FIELD_MATCH = 0
CIGAR_FIELD_INSERTION = 1
CIGAR_FIELD_DELETION = 2
CIGAR_FIELD_SOFT_CLIP = 4
CIGAR_FIELD_HARD_CLIP = 5

# The complete set of codes for cigar tuples
# in pysam is:
    #M BAM_CMATCH 0
    #I BAM_CINS 1
    #D BAM_CDEL 2
    #N BAM_CREF_SKIP 3
    #S BAM_CSOFT_CLIP 4
    #H BAM_CHARD_CLIP 5
    #P BAM_CPAD 6
    #= BAM_CEQUAL 7
    #X BAM_CDIFF 8
    #B BAM_CBACK 9

# Symbols use to fill the
# MD field with insertions and deletions
SYMBOL_INSERTION = "I"
SYMBOL_DELETION = "^"

FILE_TYPE_BAM = "BAM"
FILE_TYPE_SAM = "SAM"

DEFAULT_TARGET_PARAM = "all"
DEFAULT_START_PARAM = -1
DEFAULT_END_PARAM = -1
DEFAULT_KMER_PARAM = 50
DEFAULT_DEPTH_PARAM = 0
DEFAULT_QUALITY_PARAM = -1

## Usage
##
__usage = "usage: kmeleon_extract.py [OPTIONS] -b BAM_FILE|-s SAM_FILE\n"+\
          "Note that this software outputs to stderr and stdout.\n\n"+\
          "typical command: kmeleon_extract.py -q 60 -r -d 4 -b demo_data/demo.2_19.bam > demo_data/demo.2_19.kmers\n\n"+\
          "See options running kmeleon_extract.py --help\n"
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

optParser.add_option('-q', '--quality', action='store', dest='quality_param', type='int', \
                    help='The minimum mapping quality (MAPQ) of a read to be considered. '+\
                    '(default: '+str(DEFAULT_QUALITY_PARAM)+')')

optParser.add_option('-d', '--depth', action='store', dest='depth_param', type='int', \
                    help='The minimum times a k-mer is found to be reported in the output. '+\
                    '(default: '+str(DEFAULT_DEPTH_PARAM)+')')

optParser.add_option('-b', '--bam', action='store', dest='bam_param', type='string', \
                    help='A BAM file to process. '+\
                    'Either the -b or the -s option is required.')

optParser.add_option('-s', '--sam', action='store', dest='sam_param', type='string', \
                    help='A SAM file to process. '+\
                    'Either the -s or the -b option is required.')

optParser.add_option('-r', '--read', action='store_true', dest='read_param', \
                    help='When -r is set, the kmer sequence is also output.')

########### Read and prepare parameters
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
        optParser.exit(0, __usage+"\n"+"The --start and --end options require a specific --target (-t) option.\n")
    else:
        start_param = options.start_param
else:
    start_param = DEFAULT_START_PARAM

## End
if options.end_param:
    if target_param == DEFAULT_TARGET_PARAM:
        optParser.exit(0, __usage+"\n"+"The --start and --end options require a specific --target (-t) option.\n")
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

## minimum read mapping quality
if options.quality_param:
    quality_param = options.quality_param
else:
    quality_param = DEFAULT_QUALITY_PARAM
    
## include DNA sequence of the kmer in the output
if options.read_param:
    read_param = True
else:
    read_param = False
    
## BAM file
file_type = ""
if options.bam_param:
    input_file_param = options.bam_param
    file_type = FILE_TYPE_BAM
elif options.sam_param:
        input_file_param = options.sam_param
        file_type = FILE_TYPE_SAM
else:
    optParser.exit(0, __usage+"\n"+"A BAM file or SAM file must be specified (with the -b or -s options, respectively).\n")

#########################
######################### BEGIN
#########################

sys.stderr.write("Reading file "+input_file_param+" to process target "+target_param+"\n")

# each kmer is extracted from the read
# from the current position +- window_side
# <--------------i-------------->
window_side = kmer_param / 2

# The global structures which keep the positions and kmers
# for each chromosome/contig in the reference
refs_dict = {} # chromosome, contig, ... --> positions --> kmers
refs_list = [] # the list of chromosomes, contigs, ...

## TODO: Add exception handling to this
# Open the input (BAM/SAM) file
input_file = get_input_file(input_file_param, file_type)

# Count the number of reads to process
numreads = get_samfile_count(input_file, target_param, start_param, end_param)

sys.stderr.write("Num of reads to process "+str(numreads)+"\n")

numreads_processed = 0
numreads_lowqual = 0
numreads_unmapped = 0
header_printed = False
first_to_flush_pos = -1

# Obtain the iterator to read the input file
samfile_iter = get_samfile_iter(input_file, target_param, start_param, end_param)

# For each read in the input file
#
for read in samfile_iter:
    
    # CPCantalapiedra 20180615
    # skip those reads under the quality (MAPQ) threshold
    if read.mapping_quality < quality_param:
        numreads_lowqual+=1
    
    # skip unmapped reads
    elif read.is_unmapped:
        numreads_unmapped+=1
        
    else:
        read_id = read.query_name
        # this is the chromosome or contig of the current mapping
        # however this is an id, to obtain the name we need to do:
        # input_file.getrname(reference_id)
        read_ref_id = read.reference_id
        read_ref_name = input_file.getrname(read_ref_id)
        
        if read.is_reverse: read_strand = "-"
        else: read_strand = "+"
        
        try:
            read_md_z = read.get_tag("MD")
            # not sure whether the next is really necessary
            if read_md_z.startswith("0"): read_md_z = read_md_z[1:]
        except KeyError:
            sys.stderr.write("WARNING: No MD:Z field in read "+read_id+"\n")
            if verbosity == VERBOSE_ALL:
                sys.stderr.write(str(read)+"\n")
            sys.stderr.write("\t you could wish to run samtools calmd to obtain MD:Z field.\n")
            #sys.stderr.write(str(read)+"\n")
            #sys.stderr.write("\n")
            #break
            continue
        
        # Add insertions to expanded MD:Z
        cigartuples = read.cigartuples
        read_left_clip = read.query_alignment_start
        
        read_seq = read.seq
        
        ########### The next are very important steps
        ########### in which the MD:Z and CIGAR fields are
        ########### processed to obtain a final list
        ########### of symbols representing the alignment
        ########### of the read with the reference
        ########### kmers will be retrieved afterwards
        ########### from this list of symbols
        # 3G^T2A --> [1, 1, 1, G, ^, T, 1, 1, A]
        expanded_md_z = f_expand_md_z(read_md_z)
        
        (final_md_z, read_seq) = f_add_insertions(expanded_md_z, cigartuples, read_seq, read_left_clip)
        
        # obtain starting and ending position of read
        
        read_pos_start = 0 #read_pos_start = 1
        #read_pos_end = read.query_alignment_end - read.query_alignment_start
        read_pos_end = len(final_md_z)
        
        # positions of the read alignment in the reference
        ref_pos_start = read.reference_start
        ref_pos_end = read.reference_end
        
        pos_range = range(read_pos_start + window_side, read_pos_end - window_side)
        
        #if read_id == "HWI-ST1450_131:6:2308:3608:17497#5@0" or \
        #    read_id == "HWI-ST1450_131:6:1314:15403:49296#5@0":
        #    print "ID: "+read_id
        #    print "SEQ: "+read_seq
        #    print "REF: "+read_ref_name
        #    print "strand: "+read_strand
        #    print "CIGAR: "+read.cigarstring
        #    print "MD:Z: "+read_md_z
        #    print "cigar tuples: "+" ".join([str(x) for x in cigartuples])
        #    print "expanded MD:Z: "+str(expanded_md_z)
        #    print "MD:Z with insertions: "+str(final_md_z)
        #    print "\tlen: "+str(len(final_md_z))
        #    print "Ref start: "+str(read.reference_start)
        #    print "Ref end: "+str(read.reference_end)
        #    print "Read align start: "+str(read.query_alignment_start)
        #    print "Read align end: "+str(read.query_alignment_end)
        #    print "Read left clip: "+str(read_left_clip)
        #    print "Read pos end: "+str(read_pos_end)
        #    print "Pos range: "+str(pos_range)
        #    print ""
        
        # Initialize the flush variable to control when to output
        # kmers from already processed reads
        if first_to_flush_pos == -1:
            first_to_flush_pos = ref_pos_start # a flag to trigger output process
        
        ################# From this point, the read is traversed
        ################# to extract and record each single kmer
        if len(pos_range) > 0:
            prev_insertions = 0
            for i in pos_range:
                
                # this is the symbol in the center of the kmer
                central_nt = final_md_z[i]
                
                # if the symbol is an insertion, the kmer is skipped
                # since there is no actual kmer mapped to the reference
                # with the nucleotides which are in the read insertion
                if central_nt == SYMBOL_INSERTION: continue
                
                clip_start = i - window_side
                clip_end = i + window_side + 1
                
                # Insertions have to be included in the kmer
                # so that this kmer can be differentiated from
                # those kmers without the insertions.
                # Therefore, the kmer size must be increased to make
                # room for the inserted nucleotides.
                #
                # TODO: this probably should be a while() loop
                # checking whether new insertions have been added after increasing
                # the clipping interval.
                # Check whether the loop works, and remove the previous code
                clipped_md_z = final_md_z[clip_start:clip_end]
                num_insertions = clipped_md_z.count(SYMBOL_INSERTION)
                prev_insertions = -1
                while (num_insertions != prev_insertions):
                    clip_end = i + window_side + 1 + num_insertions
                    clipped_md_z = final_md_z[clip_start:clip_end]
                    prev_insertions = num_insertions
                    num_insertions = clipped_md_z.count(SYMBOL_INSERTION)
                
                # if the clipping goes beyond the length of the kmer
                # then the kmer is not output, since will cause shorter
                # kmers than expected
                if clip_end > len(final_md_z): continue
                
                #clip_end = i + window_side + 1 + num_insertions
                #clipped_md_z = final_md_z[clip_start:clip_end]
                
                # clip the kmer sequence.
                # the soft clipping is not included in the MD:Z field
                # BUT it is included in the read sequence, and thus
                # the soft left clipping has to be considered when
                # extracting the kmer sequence from the read
                kmer_seq = read_seq[clip_start+read_left_clip:clip_end+read_left_clip]
                
                ###################################################
                # the next represents the final kmer to be recorded
                collapsed_md_z = f_collapse_md_z(clipped_md_z)
                
                # obtain the kmer position
                #kmer_contig_position = f_get_kmer_contig_position(i, ref_pos_start, read_left_clip)
                kmer_contig_position = ref_pos_start + (i - 1)
                
                # Insertions before the central nt of this kmer
                # have to be counted, so that the actual position
                # of the kmer is computed correctly
                prev_insertions = final_md_z[:i].count(SYMBOL_INSERTION)
                kmer_contig_position -= prev_insertions
                
                #if read_id == "HWI-ST1450_131:6:1213:21225:9778#5@0" or \
                #    read_id == "HWI-ST1450_131:6:2211:4897:90286#5@0":
                # if kmer_contig_position == 308378329:
                #     print "ID: "+read_id
                #     print "Pos in range: "+str(i)
                #     print "Central nt: "+str(central_nt)
                #     print "Prev insertions: "+str(prev_insertions)
                #     print "Num insertions: "+str(num_insertions)
                #     print "Read MD:Z "+str(final_md_z)
                #     print "\tlen: "+str(len(final_md_z))
                #     print "Clip start: "+str(clip_start)
                #     print "Clip left: "+str(read_left_clip)
                #     print "Clip end: "+str(clip_end)
                #     print "Clipped md:z: "+str(clipped_md_z)
                #     print "\tlen: "+str(len(clipped_md_z))
                #     print "Collapsed md:z"+str(collapsed_md_z)
                #     print "kmer position: "+str(kmer_contig_position)
                #     print "kmer seq: "+kmer_seq
                #     print "\tlen: "+str(len(kmer_seq))
                #     print ""
                
                # record the kmer in its position
                #
                f_add_kmer_depths(refs_dict, refs_list,
                                  read_ref_id, read_ref_name,
                                  kmer_contig_position, collapsed_md_z, kmer_seq)
            
            ################ Generate output
            #
            # Header (only the first time)
            if not header_printed:
                #if store_depths:
                if read_param:
                    sys.stdout.write("@Target\tPosition\tkmer(md_z)\tdp\tkmer(DNA)\n")
                else:
                    sys.stdout.write("@Target\tPosition\tkmer(md_z)\tdp\n")
                #else:
                #    sys.stdout.write("@ Position\tMD_Z\n")
                header_printed = True
            
            #
            # kmers (for those reads which have been processed already)
            if ref_pos_start - FLUSH_INTERVAL*2 > first_to_flush_pos:
                (refs_list, refs_dict) = f_print_previous_kmers(refs_dict, refs_list,
                                                                read_ref_id, ref_pos_start,
                                                                depth_param, FLUSH_INTERVAL)
                
                if len(refs_list) > 0 and len(refs_dict[refs_list[0]]["pos_list"]) > 0:
                    first_to_flush_pos = refs_dict[refs_list[0]]["pos_list"][0]
                else:
                    first_to_flush_pos = -1 # will be updated to ref_pos_start in the next iteration
        
        # Here, the read has been processed completely
        # and its kmers have been extracted and recorded.
        numreads_processed+=1
    
    numreads_total = numreads_processed+numreads_lowqual+numreads_unmapped
    if (numreads_total % 1000 == 0):
        sys.stderr.write(str(numreads_total)+" reads parsed (out of "+str(numreads)+").\n")

input_file.close()

################ Generate output for the remaining data
#
# kmers which have not been "flushed"
for read_ref_id in refs_list:
    pos_list = refs_dict[read_ref_id]["pos_list"]
    pos_dict = refs_dict[read_ref_id]["pos_dict"]
    read_ref_name = refs_dict[read_ref_id]["ref_name"]
    
    for kmer_position in pos_list:
        if kmer_position in pos_dict:
            pos_kmers_dict = pos_dict[kmer_position]
            f_print_pos_kmers(read_ref_name, kmer_position, pos_kmers_dict, depth_param)
            
        else:
            raise Exception("Position "+str(kmer_position)+" is not in dict")

sys.stderr.write("Num of reads processed "+str(numreads_processed)+"\n")
sys.stderr.write("Num of low quality reads "+str(numreads_lowqual)+"\n")
sys.stderr.write("Num of unmapped reads "+str(numreads_unmapped)+"\n")
numreads_total = numreads_processed+numreads_lowqual+numreads_unmapped
sys.stderr.write("Total reads "+str(numreads_total)+"\n")

sys.stderr.write("Finished.\n")

## END