#!/usr/bin/env python
# -*- coding: utf-8 -*-

# CPCantalapiedra 2018

import sys
from util import *
from ..counts.Count import *

def f_is_new_interval(curr_target, curr_pos, curr_count, prev_target, prev_pos, prev_count):
    is_new_interval = False
    
    if curr_target == prev_target and \
            curr_pos == prev_pos + 1 and \
            curr_count == prev_count:
                is_new_interval = False
    else:
        is_new_interval = True
    
    return is_new_interval

def f_new_interval(target, start, count):
    new_interval = {}
    
    new_interval["target"] = target
    new_interval["start"] = start-1 # BED format is 0-based
    new_interval["end"] = start
    new_interval["count"] = count
    
    return new_interval

def f_add_to_interval(interval, target, pos, count):
    
    interval["end"] = pos
    
    return interval

### This is the main function for MODE_CONSTANT
# This is the simplest mode. Just, the counts are parsed
# and every new interval identified as separated from the next,
# either due to a basepair span or a change in kmer count,
# is printed. Therefore, if span_param==0 the raw intervals
# representing the kmers are printed.
def f_intervals(counts_fileobj, span_param, min_kc_param):
    
    prev_data = None
    prev_target = ""
    prev_pos = -1
    prev_count = -1
    curr_interval = None
    for i, line in enumerate(counts_fileobj):
        
        if i==0: continue
        
        line_data = line.strip().split("\t")
        curr_target = line_data[COUNTS_FILE_FIELD_TARGET]
        curr_pos = long(line_data[COUNTS_FILE_FIELD_POSITION])
        curr_count = int(line_data[COUNTS_FILE_FIELD_COUNT])
        
        #print prev_data
        #print line_data
        
        is_new_interval = f_is_new_interval(curr_target, curr_pos, curr_count, prev_target, prev_pos, prev_count)
        
        #print is_new_interval
        
        if is_new_interval:
            if curr_interval:
                prev_interval = f_add_to_interval(curr_interval, prev_target, prev_pos, prev_count)
                if (f_get_interval_length(prev_interval) >= span_param) and\
                    (f_get_interval_count(prev_interval) >= min_kc_param):
                    
                    f_print_interval(prev_interval)
                
            curr_interval = f_new_interval(curr_target, curr_pos, curr_count)
        
        prev_target = curr_target
        prev_pos = curr_pos
        prev_count = curr_count
        prev_data = line_data
    
    # Last interval
    if curr_interval:
        prev_interval = f_add_to_interval(curr_interval, prev_target, prev_pos, prev_count)
        if (f_get_interval_length(prev_interval) >= span_param) and\
            (f_get_interval_count(prev_interval) >= min_kc_param):
            
            f_print_interval(prev_interval)
    
    return

## END