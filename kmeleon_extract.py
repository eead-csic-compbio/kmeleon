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

sys.stderr.write("Running kmeleon extract...\n")

############## FUNCTION DEFINITIONS
##############

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

# Adds a kmer to an absolute position
def f_add_kmer(kmer_position, md_z):
    
    #sys.stdout.write(str(kmer_position)+" - "+str(md_z)+"\n")
    
    md_z_string = "".join([str(x) for x in md_z])
    
    if kmer_position in pos_dict:
        pos_kmers_list = pos_dict[kmer_position]
        if md_z_string not in set(pos_kmers_list):
            pos_kmers_list.append(md_z_string)
        
    else:
        pos_dict[kmer_position] = [md_z_string]
        
        # additionally, add position to list
        pos_list.append(kmer_position)
    
    return

# Adds a kmer to an absolute position,
# and counts how many times such kmer has been
# added to such position (depth of kmer at that position)
def f_add_kmer_depths(kmer_position, md_z):
    
    #sys.stdout.write(str(kmer_position)+" - "+str(md_z)+"\n")
    
    md_z_string = "".join([str(x) for x in md_z])
    
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

##
def f_print_previous_kmers(curr_pos, flush_interval):
    new_list = []
    new_dict = {}
    
    ################ Generate output
    # Rows
    for kmer_position in pos_list:
        if kmer_position in pos_dict:
            pos_kmers_dict = pos_dict[kmer_position]
            
            #sys.stderr.write("Position "+str(kmer_position)+"\n")
            if curr_pos - flush_interval > kmer_position:
                for kmer in pos_kmers_dict:
                    if store_depths:
                        kmer_count = pos_kmers_dict[kmer]
                        sys.stdout.write(str(kmer_position)+"\t"+str(kmer)+"\t"+str(kmer_count)+"\n")
                    else:
                        sys.stdout.write(str(kmer_position)+"\t"+str(kmer)+"\n")
            else:
                new_list.append(kmer_position)
                new_dict[kmer_position] = pos_kmers_dict
            
        else:
            raise Exception("Position "+str(kmer_position)+" is not in dict")
    
    #print "AFTER FLUSH"
    #print curr_pos
    #print len(new_list)
    #print len(new_dict)
    #print pos_list[-1]
    #print pos_dict[pos_list[-1]]
    
    return (new_list, new_dict)


################ MAIN
################

VERBOSE_ALL = 5
verbosity = 0

CIGAR_INSERTION_FIELD = 1
CIGAR_SOFTCLIP_FIELD = 4

if ((len(sys.argv)!=5) and (len(sys.argv)!=6)):
    sys.stderr.write("\n")
    sys.stderr.write("Wrong parameter number for kmeleon extract.\n")
    sys.stderr.write("A typical kmeleon extract command would be:\n")
    sys.stderr.write("\tkmeleon_extract.py chr1 mappings.sam 50 10000\n")
    sys.stderr.write("or, requesting depths:\n")
    sys.stderr.write("\tkmeleon_extract.py chr1 mappings.sam 50 10000 depths\n")
    sys.stderr.write("\n")
    sys.stderr.write("We need 5 parameters and another optional one:\n")
    sys.stderr.write("\t- Target: for example the chromosome number.\n")
    sys.stderr.write("\t- Mappings file: the mappings to process. A SAM or BAM file would be great.\n")
    sys.stderr.write("\t- kmer size: a number which indicates the size of sequences to search for. 50 for example.\n")
    sys.stderr.write("\t- flush interval: a number which allows controlling how much memory is used by kmeleon. 10000 or 20000 is ok in general.\n")
    sys.stderr.write("\t- Optional argument: if the word 'depths' is given, the depth of each kmer found will be also reported.\n")
    sys.stderr.write("\n")
    sys.exit(-1)

else:
    target_seq = sys.argv[1]
    sam_file = sys.argv[2]
    window_size = int(sys.argv[3])
    FLUSH_INTERVAL = int(sys.argv[4]) #20000
    store_depths = False
    if len(sys.argv)==6:
        if sys.argv[5] == "depths":
            store_depths = True

sys.stderr.write("Reading file "+sam_file+" to process target "+target_seq+"\n")

#window_size = 50
window_side = window_size / 2

pos_dict = {}
pos_list = []

## TODO: Add exception handling to this
samfile = pysam.AlignmentFile(sam_file, "rb")

numreads = samfile.count(target_seq)

numreads_processed = 0
header_printed = False
first_to_flush_pos = -1

#for line in open(samfile, 'r'):
for read in samfile.fetch(target_seq):#samfile.fetch('chr1', 100, 120):
    
    #print read
    
    read_id = read.query_name
    
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
            
            if store_depths:
                f_add_kmer_depths(kmer_contig_position, collapsed_md_z)
            else:
                f_add_kmer(kmer_contig_position, collapsed_md_z)
        
        ################ Generate output
        # Header
        if not header_printed:
            if store_depths:
                sys.stdout.write("@ Position\tMD_Z\tcount\n")
            else:
                sys.stdout.write("@ Position\tMD_Z\n")
            header_printed = True
        
        if ref_pos_start - FLUSH_INTERVAL*2 > first_to_flush_pos:
            (pos_list, pos_dict) = f_print_previous_kmers(ref_pos_start, FLUSH_INTERVAL)
            #print "POS RANGE"
            #print read
            #print pos_range
            #print read_pos_start
            #print read_pos_end
            #print window_side
            first_to_flush_pos = pos_list[0]
    
    numreads_processed+=1
    if (numreads_processed % 1000 == 0):
        sys.stderr.write(str(numreads_processed)+" reads processed (out of "+str(numreads)+").\n")

samfile.close()

################ Generate output for the remaining data
# Rows
for kmer_position in pos_list:
    if kmer_position in pos_dict:
        #sys.stderr.write("Position "+str(kmer_position)+"\n")
        pos_kmers_dict = pos_dict[kmer_position]
        for kmer in pos_kmers_dict:
            if store_depths:
                kmer_count = pos_kmers_dict[kmer]
                sys.stdout.write(str(kmer_position)+"\t"+str(kmer)+"\t"+str(kmer_count)+"\n")
            else:
                sys.stdout.write(str(kmer_position)+"\t"+str(kmer)+"\n")
        
    else:
        raise Exception("Position "+str(kmer_position)+" is not in dict")
    
sys.stderr.write("Finished.\n")

## END
