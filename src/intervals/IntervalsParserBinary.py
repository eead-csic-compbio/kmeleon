#!/usr/bin/env python
# -*- coding: utf-8 -*-

# CPCantalapiedra 2018

import sys
from util import *
from ..counts.Count import *

def f_is_new_interval(curr_target, curr_pos, curr_count, prev_target, prev_pos, prev_count):
    is_new_interval = False
    
    if curr_target == prev_target and curr_pos == prev_pos + 1:
        
        if curr_count == 1 and prev_count == 1:
            is_new_interval = False
        elif curr_count > 1 and prev_count > 1:
            is_new_interval = False
        else: # curr_count == 1 and prev_count > 1 OR curr_count > 1 and prev_count == 1
            is_new_interval = True
    else:
        is_new_interval = True
    
    return is_new_interval

def f_new_interval(target, start, count):
    new_interval = {}
    
    new_interval["target"] = target
    new_interval["start"] = start-1
    new_interval["end"] = start
    new_interval["count"] = count
    new_interval["counts"] = [count]
    
    return new_interval

def f_add_to_interval(interval, target, pos, count):
    
    interval["end"] = pos
    interval["counts"].append(count)
    interval["count"] = sum(interval["counts"])*1.0/len(interval["counts"])
    
    return interval
    
# This is the main function of the MODE_BINARY
# This is like the constant mode, but here only intervals
# with k-mer count == 1 or k-mer count > 1 are differentiated to each other.
# Therefore, it changes the identification of a new interval, specially for
# those basepairs with k-mer count > 1, since it is considered the same
# whether is k-mer count == 2 or ==3, for example, and longer intervals will
# be created in general. This way, the k-mer count given for an interval of
# k-mer count > 1 is an average of the counts found for all the basepairs
# in the interval, instead of a direct translation of the k-mer count of
# every position.
def f_intervals(counts_fileobj, span_param):
    
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
                if f_get_interval_length(prev_interval) >= span_param:
                    f_print_interval(prev_interval)
                
            curr_interval = f_new_interval(curr_target, curr_pos, curr_count)
        else:
            curr_interval = f_add_to_interval(curr_interval, prev_target, prev_pos, prev_count)
        
        prev_target = curr_target
        prev_pos = curr_pos
        prev_count = curr_count
        prev_data = line_data
    
    # Last interval
    if curr_interval:
        prev_interval = f_add_to_interval(curr_interval, prev_target, prev_pos, prev_count)
        if f_get_interval_length(prev_interval) >= span_param:
            f_print_interval(prev_interval)
    
    return

## END