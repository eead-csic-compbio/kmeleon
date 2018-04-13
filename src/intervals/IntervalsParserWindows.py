#!/usr/bin/env python
# -*- coding: utf-8 -*-

# CPCantalapiedra 2018

import sys
from util import *
from ..counts.Count import *

from numpy import mean

## MODE WINDOWS

def f_is_new_interval(curr_target, curr_pos, curr_count, prev_target, prev_pos, prev_count):
    is_new_interval = False
    
    if curr_target == prev_target and \
            curr_pos == prev_pos + 1:
                is_new_interval = False
    else:
        is_new_interval = True
    
    return is_new_interval

def f_new_interval(target, start, window_param, count=0):
    new_interval = {}
    
    new_interval["target"] = target
    
    # the new interval previously used start as start and end position
    # (see f_new_interval function commented below)
    # however, to make all the samples comparable
    # all the intervals start always in the closest factor
    # of the window size
    # e.g. a start at position 1615 will create a window starting
    # at 1500 if the window_param = 500,
    # or at 1600 if the window_param = 100 or 50
    remaind = start % window_param
    intvstart = start - remaind + 1
    new_interval["start"] = intvstart-1
    new_interval["end"] = intvstart-1 # an empty window
    new_interval["count"] = count
    new_interval["counts"] = []
    
    return new_interval

#def f_new_interval(target, start, window_param, count=0):
#    new_interval = {}
#    
#    new_interval["target"] = target
#    new_interval["start"] = start-1
#    new_interval["end"] = start-1 # an empty window
#    new_interval["count"] = count
#    new_interval["counts"] = []
#    
#    return new_interval

def f_add_to_interval(interval, target, pos, count):
    interval["end"] = pos
    interval["counts"].append(count)
    return

def f_compute_count(window):
    window["count"] = mean([count for count in window["counts"] if count>0])*1.0
    #window["count"] = mean(window["counts"])*1.0
    return

def f_get_remaining(window, window_param):
    remaining = 0
    
    if window:
        remaining = window_param - (window["end"] - window["start"])
    
    return remaining

def f_fill_window(window, fill):
    for i in range(1, fill+1):
        window["end"] += 1
        window["counts"].append(0)
        
    return

def f_windows(counts_fileobj, window_param):
    intervals_list = []
    
    current_window = None
    prev_target = ""
    prev_pos = -1
    prev_count = -1
    for i, line in enumerate(counts_fileobj):
        
        if i==0: continue
        
        line_data = line.strip().split("\t")
        curr_target = line_data[COUNTS_FILE_FIELD_TARGET]
        curr_pos = long(line_data[COUNTS_FILE_FIELD_POSITION])
        curr_count = int(line_data[COUNTS_FILE_FIELD_COUNT])
        
        if not current_window:
            current_window = f_new_interval(curr_target, curr_pos, window_param)
            f_add_to_interval(current_window, curr_target, curr_pos, curr_count)
        else: # there is current window
            is_new_interval = f_is_new_interval(curr_target, curr_pos, curr_count, prev_target, prev_pos, prev_count)
            
            if is_new_interval:
                remaining = f_get_remaining(current_window, window_param)
                distance = curr_pos - prev_pos
                
                # check it is within range of this window
                if curr_target == prev_target and remaining >= distance:
                    f_fill_window(current_window, distance - 1)
                else:
                    # finish the previous window
                    if current_window:
                        f_fill_window(current_window, remaining)
                        f_compute_count(current_window)
                        intervals_list.append(current_window)
                    # create a new window
                    current_window = f_new_interval(curr_target, curr_pos, window_param)
            # else: fill the current window
            
            f_add_to_interval(current_window, curr_target, curr_pos, curr_count)
            # if the current window has enough size, store it
            window_length = f_get_interval_length(current_window)
            if window_length >= window_param:
                f_compute_count(current_window)
                intervals_list.append(current_window)
                current_window = None
        
        prev_target = curr_target
        prev_pos = curr_pos
        prev_count = curr_count
    
    # last unfinished window
    if current_window:
        remaining = f_get_remaining(current_window, window_param)
        if remaining > 0:
            f_fill_window(current_window, remaining)
        f_compute_count(current_window)
        intervals_list.append(current_window)
    
    return intervals_list

# This is the main function of MODE_WINDOWS
def f_intervals(counts_fileobj, window_param):
    
    # Calculate median k-mer count over moving windows
    
    intervals_list = f_windows(counts_fileobj, window_param)
    
    # Finally, print the intervals
    
    sorted_list = sorted(intervals_list, key=lambda x:x["start"])
    for interval in sorted_list:
        f_print_interval(interval)
    
    return

## END