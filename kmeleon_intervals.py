#!/usr/bin/env python
# -*- coding: utf-8 -*-

# CPCantalapiedra 2018

import sys, gzip
from optparse import OptionParser

import src.intervals.IntervalsParserRaw
import src.intervals.IntervalsParserBinary
import src.intervals.IntervalsParserSmooth
import src.intervals.IntervalsParserWindows
    
#########
######### BEGIN

################ GLOBALS AND PARAMETERS

#COUNTS_FILE_FIELD_TARGET = 0
#COUNTS_FILE_FIELD_POSITION = 1
#COUNTS_FILE_FIELD_COUNT = 2

# MODES
# These are the different modes to compute the intervals
# Check the definitions of the functions f_intervals_MODE
# for each respective MODE for further explanations.
MODE_CONSTANT = "raw"
MODE_BINARY = "binary"
MODE_SMOOTH = "smooth"
MODE_WINDOWS = "windows"

DEFAULT_MODE = MODE_BINARY
DEFAULT_SPAN = 50
DEFAULT_MIN_KC = 0
DEFAULT_WINDOW = 500

## Usage
__usage = "usage: kmeleon_intervals.py [OPTIONS] COUNTS_FILE\n"+\
          "The COUNTS file has the format of the output of kmeleon_count.py\n"+\
          "Note that this software outputs to stderr and stdout.\n\n"+\
          "typical command:\nkmeleon_intervals.py -m binary -s 50 demo_data/demo.2_19.counts > demo_data/demo.2_19.intervals"
optParser = OptionParser(__usage)

optParser.add_option('-m', '--mode', action='store', dest='mode_param', type='string', \
                    help='How intervals are computed. Either "'+MODE_WINDOWS+'", "'+MODE_CONSTANT+\
                    '", "'+MODE_BINARY+'" or "'+MODE_SMOOTH+'". '+\
                    '(default: '+str(DEFAULT_MODE)+'). These modes are explained in the README file.')

optParser.add_option('-s', '--span', action='store', dest='span_param', type='int', \
                    help='The minimum span of bases to be considered a whole interval. '+\
                    'Not used with option -w. '+\
                    '(default: '+str(DEFAULT_SPAN)+')')

optParser.add_option('--min-kc', action='store', dest='min_kc_param', type='int', \
                    help='Minimum kmer_count in a interval to be reported. '+\
                    'Set to 0 to report all the intervals. '+\
                    '(default: '+str(DEFAULT_MIN_KC)+')')

optParser.add_option('-d', '--diploid', action='store_true', dest='diploid_param', \
                    help='Requires -m binary. When -d is set, positions with kmer_count=2 '+\
                    'will be treated as regular heterozygous variants and the interval value '+\
                    ' will be set to 0 as with kmer_count==1. If -d is not set, intervals with kmer_count==2 will be '+\
                    'set to 1, as with kmer_count>2.')

optParser.add_option('-w', '--window', action='store', dest='window_param', type='int', \
                    help='Requires -m windows. The size of the window for which mean k-mer counts will be processed. '+\
                    '(default: '+str(DEFAULT_WINDOW)+')')

########### Read parameters
###########

(options, arguments) = optParser.parse_args()

sys.stderr.write("Running kmeleon intervals...\n")

# THIS IS MANDATORY
if not arguments or len(arguments)==0:
    optParser.exit(0, "You may wish to run '-help' option.\n")
# INPUT FILE
counts_file = arguments[0]

## mode
if options.mode_param:
    mode_param = options.mode_param
else:
    mode_param = DEFAULT_MODE
    
## span
if options.span_param != None:
    span_param = options.span_param
else:
    span_param = DEFAULT_SPAN

## min_kc
if options.min_kc_param != None:
    min_kc_param = options.min_kc_param
else:
    min_kc_param = DEFAULT_MIN_KC
    
## how to manage heterozygous variants in binary mode
if options.diploid_param:
    diploid_param = True
else:
    diploid_param = False

## window mode
if options.window_param != None:
    window_param = options.window_param
else:
    window_param = DEFAULT_WINDOW

refs_dict = {}
refs_list = []

ext_pos = counts_file.rfind(".")
ext = counts_file[ext_pos+1:]

if ext == "gz" or ext == "gzip":
    counts_fileobj = gzip.open(counts_file, 'r')
else:
    counts_fileobj = open(counts_file, 'r')

if mode_param == MODE_CONSTANT:
    src.intervals.IntervalsParserRaw.f_intervals(counts_fileobj, span_param, min_kc_param)
elif mode_param == MODE_BINARY:
    src.intervals.IntervalsParserBinary.f_intervals(counts_fileobj, span_param, min_kc_param, diploid_param)
elif mode_param == MODE_SMOOTH:
    src.intervals.IntervalsParserSmooth.f_intervals(counts_fileobj, span_param, min_kc_param)
elif mode_param == MODE_WINDOWS:
    src.intervals.IntervalsParserWindows.f_intervals(counts_fileobj, window_param, min_kc_param)
else:
    raise Exception("Unknown mode "+str(mode_param)+".")

sys.stderr.write("Finished.\n")

## END
