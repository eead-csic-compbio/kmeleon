#!/usr/bin/env python
# -*- coding: utf-8 -*-

# CPCantalapiedra 2018

import sys
from util import *
from ..counts.Count import *

from IntervalsParserRaw import f_is_new_interval
from IntervalsParserBinary import f_add_to_interval

## MODE SMOOTH

def f_new_interval(prev_interval, target, start, count):
    new_interval = {}
    
    new_interval["target"] = target
    new_interval["start"] = start-1
    new_interval["end"] = start
    new_interval["count"] = count
    new_interval["counts"] = [count]
    
    new_interval["prev_interval"] = prev_interval
    if prev_interval:
        prev_interval["next_interval"] = new_interval
    
    new_interval["next_interval"] = None
    
    #new_interval["smoothed"] = False # this is for later use, when joining intervals together
    new_interval["dropped"] = False
    
    return new_interval

def f_add_next_interval(prev_interval, curr_interval):
    
    if prev_interval: prev_interval["next_interval"] = curr_interval
    
    if curr_interval: curr_interval["prev_interval"] = prev_interval
    
    return

def f_is_adjacent(interval, next_interval):
    is_adjacent = False
    
    if interval and next_interval:
        if interval["target"] == next_interval["target"] and \
            interval["end"] == next_interval["start"]:
            is_adjacent = True
        else:
            is_adjacent = False
    else:
        is_adjacent = False
    
    return is_adjacent

def f_drop_interval(interval):
    
    prev_interval = interval["prev_interval"]
    next_interval = interval["next_interval"]
    
    # update the links of the next and prev intervals
    if prev_interval:
        prev_interval["next_interval"] = next_interval
    if next_interval:
        next_interval["prev_interval"] = prev_interval
    
    interval["dropped"] = True
    
    return

def f_smooth_prev_interval(prev_interval, interval):
    
    prev_interval["end"] = interval["end"]
    
    prev_interval["counts"].extend(len(interval["counts"])*[prev_interval["count"]])
    prev_interval["count"] = sum(prev_interval["counts"])*1.0/len(prev_interval["counts"])
    
    f_drop_interval(interval)
    
    return prev_interval

def f_smooth_next_interval(interval, next_interval):
    
    next_interval["start"] = interval["start"]
    
    next_interval["counts"].extend(len(interval["counts"])*[next_interval["count"]])
    next_interval["count"] = sum(next_interval["counts"])*1.0/len(next_interval["counts"])
    
    f_drop_interval(interval)
    
    return interval

def f_smooth_intervals(intervals_list, span_param):
    
    did_change = True
    loops = 0
    
    while (did_change):
        did_change = False
        
        sorted_list = sorted(intervals_list, key=lambda x:x["count"], reverse=True)
        
        loops += 1
        
        for interval in sorted_list:
            
            # If the interval has already been dropped
            # just continue with the next one
            if interval["dropped"]: continue
            
            interval_length = f_get_interval_length(interval)
            interval_count = f_get_interval_count(interval)
            
            prev_interval = interval["prev_interval"]
            prev_is_adjacent = f_is_adjacent(prev_interval, interval)
            prev_interval_count = f_get_interval_count(prev_interval)
            next_interval = interval["next_interval"]
            next_is_adjacent = f_is_adjacent(interval, next_interval)
            next_interval_count = f_get_interval_count(next_interval)
            
            # If the interval is large enough, only fuse it with
            # the previous-next intervals if they have the same "count"
            if interval_length >= span_param:
                
                if prev_is_adjacent and prev_interval_count == interval_count:
                    new_interval = f_smooth_prev_interval(prev_interval, interval)
                    did_change = True
                elif next_is_adjacent and next_interval_count == interval_count:
                    new_interval = f_smooth_next_interval(interval, next_interval)
                    did_change = True
                # else: pass
            
            # If the interval is short, we have to smooth it if possible
            #
            else:
                if prev_is_adjacent and next_is_adjacent:
                    if prev_interval_count >= next_interval_count:
                        new_interval = f_smooth_prev_interval(prev_interval, interval)
                        did_change = True
                        
                    else: # prev_interval["count"] < next_interval["count"]:
                        new_interval = f_smooth_next_interval(interval, next_interval)
                        did_change = True
                    
                elif prev_is_adjacent:
                    new_interval = f_smooth_prev_interval(prev_interval, interval)
                    did_change = True
                    
                elif next_is_adjacent:
                    new_interval = f_smooth_next_interval(interval, next_interval)
                    did_change = True
                    
                else: # Too short interval and no other interval is adjacent: drop the interval
                    f_drop_interval(interval)

        intervals_list = [interval for interval in intervals_list if interval["dropped"] == False]
    
    sys.stderr.write("LOOPS "+str(loops)+"\n")
    
    return intervals_list

def f_store_intervals(counts_fileobj):
    intervals_list = []
    
    prev_data = None
    prev_target = ""
    prev_pos = -1
    prev_count = -1
    curr_interval = None
    prev_interval = None
    for i, line in enumerate(counts_fileobj):
        
        if i==0: continue
        
        line_data = line.strip().split("\t")
        curr_target = line_data[COUNTS_FILE_FIELD_TARGET]
        curr_pos = long(line_data[COUNTS_FILE_FIELD_POSITION])
        curr_count = int(line_data[COUNTS_FILE_FIELD_COUNT])
        
        #print prev_data
        #print line_data
        
        is_new_interval = f_is_new_interval(curr_target, curr_pos, curr_count, prev_target, prev_pos, prev_count)
        
        if is_new_interval:
            if curr_interval:
                prev_interval = f_add_to_interval(curr_interval, prev_target, prev_pos, prev_count)
            
            curr_interval = f_new_interval(prev_interval, curr_target, curr_pos, curr_count)
            intervals_list.append(curr_interval)
            
            f_add_next_interval(prev_interval, curr_interval)
        
        prev_target = curr_target
        prev_pos = curr_pos
        prev_count = curr_count
        prev_data = line_data
        
        #print ""
    
    # Last interval
    if curr_interval:
        prev_interval = f_add_to_interval(curr_interval, prev_target, prev_pos, prev_count)
    
    curr_interval = f_new_interval(prev_interval, curr_target, curr_pos, curr_count)
    intervals_list.append(curr_interval)
    
    f_add_next_interval(prev_interval, curr_interval)
    
    return intervals_list
    
# This is the main function of MODE_SMOOTH
# In this mode, the intervals which are shorter than span_param
# are not filtered out directly. Instead, the intervals are smoothed
# to the adjacent intervals.
# Therefore, if the adjacent intervals are of smaller k-mer count
# the interval count will be changed to that smaller k-mer count.
# However, if the adjacent intervals are of greater k-mer count
# the interval count will be changed to that greater k-mer count.
# Note that first of carrying over the smoothing, the intervals
# are sorted, so that those intervals with greater k-mer count are
# processed first.
def f_intervals(counts_fileobj, span_param):
    
    # First, create the intervals from the kmer_counts
    
    intervals_list = f_store_intervals(counts_fileobj)
    
    # Second, once that the intervals are stored
    # sort them, those with largest kmer_count first
    # For each interval, check the length, and the previous and next intervals
    # Fuse the interval with a previous or next one if needed
    
    intervals_list = f_smooth_intervals(intervals_list, span_param)
    
    # Finally, print the smoothed intervals
    
    sorted_list = sorted(intervals_list, key=lambda x:x["start"])
    for interval in sorted_list:
        f_print_interval(interval)
    
    return

## END