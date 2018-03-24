#!/usr/bin/env python
# -*- coding: utf-8 -*-

# CPCantalapiedra 2018

import sys

# Base functions

def f_print_interval(interval):
    if interval:
        sys.stdout.write(interval["target"]+"\t"+
                     str(interval["start"])+"\t"+
                     str(interval["end"])+"\t"+
                     str(interval["count"])+"\n")
    else:
        raise Exception("Trying to print non-existing interval.")
    
    return

def f_get_interval_length(interval):
    interval_length = 0
    if interval:
        interval_length = interval["end"] - interval["start"]
    
    return interval_length

def f_get_interval_count(interval):
    interval_count = 0
    if interval:
        interval_count = interval["count"]
        
    return interval_count

## END
