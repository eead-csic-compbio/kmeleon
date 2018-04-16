# kmeleon

A pipeline for k-mer count analysis of genomic diversity.

#### Dependencies

One of the kmeleon scripts, kmeleon_extract.py uses the Python API [pysam](http://pysam.readthedocs.io/en/latest/api.html) to read the SAM/BAM mappings file.
Another one, kmeleon_intersection_matrix.sh uses [bedops](https://bedops.readthedocs.io/en/latest/) to perform comparison of interval files in BED format.

#### Running kmeleon

There are several steps to run a complete kmeleon analysis:
* [Extract the kmers from the mappings: kmeleon_extract.py](https://github.com/eead-csic-compbio/kmeleon#1st-extract-the-kmers-from-the-mappings)
* [Count the number of kmers at each genomic position: kmeleon_count.py](https://github.com/eead-csic-compbio/kmeleon#2nd-count-the-number-of-kmers-at-each-genomic-position)

There are also other tools to manage the results:
* [Translate the position-based kmer counts to intervals of constant kmer count](https://github.com/eead-csic-compbio/kmeleon#Translate-the-position-based-kmer-counts-to-intervals-of-constant-kmer-count)
* [Join the kmer counts of all samples in a single table: kmeleon_table.py](https://github.com/eead-csic-compbio/kmeleon#Join-the-kmer-counts-of-all-samples-in-a-single-table)
* [Obtain a matrix of pairwise distances or similarities based on intersection of intervals between samples](https://github.com/eead-csic-compbio/kmeleon#Obtain-a-matrix-of-pairwise-distances-or-similarities-based-on-intersection-of-intervals-between-samples)

## 1st) Extract the kmers from the mappings

This is done by running the script **kmeleon_extract.py**, which has the next parameters:

```
Usage: kmeleon_extract.py [OPTIONS] -b BAM_FILE|-s SAM_FILE
Note that this software outputs to stderr and stdout.

typical command: kmeleon_extract.py -d 4 -b demo_data/demo.2_19.bam > demo_data/demo.2_19.kmers

Options:
  -h, --help            show this help message and exit
  -t TARGET_PARAM, --target=TARGET_PARAM
                        A chromosome number or name, or a specific contig, or
                        "all" to process all the mappings. (default: "all".)
  --start=START_PARAM   The -t target parameter is required. Starting
                        basepairs position to process within the given target.
                        (default: -1)
  --end=END_PARAM       The -t target parameter is required. Ending basepairs
                        position to process within the given target. (default:
                        -1)
  -k KMER_PARAM, --kmer=KMER_PARAM
                        The length of fragments (k-mers) to parse. (default:
                        50)
  -d DEPTH_PARAM, --depth=DEPTH_PARAM
                        The minimum times a k-mer is found to be reported in
                        the output. (default: 0)
  -b BAM_PARAM, --bam=BAM_PARAM
                        A BAM file to process. Either the -b or the -s option
                        is required.
  -s SAM_PARAM, --sam=SAM_PARAM
                        A SAM file to process. Either the -s or the -b option
                        is required.
```

It generates a file with 4 columns, as shown in its header:

```
@Target Position	kmer(md_z)	dp
chr1   2917	51	1
chr1   2918	51	1
chr1   2919	51	1
chr1    2920	51	1
chr1    2921	51	1
chr1    2922	51	1
```

* @ is a symbol used to differentiate the header from the other rows
* Target: the chromosome or contig of the current position.
* Position: the position within the target (chromosome)
* kmer(md_z): a kmer found in such position, in SAM/BAM file MD_Z field format
* dp: depth, i.e. times this kmer has been observed at this position.

## 2nd) Count the number of kmers at each genomic position

This is done by running the script **kmeleon_count.py**, which has the next parameters:

```
Usage: kmeleon_count.py [OPTIONS] KMERS_FILE
The KMERS_FILE has the format of the output of kmeleon_extract.py
Note that this software outputs to stderr and stdout.

typical command: kmeleon_count.py -d 4 demo_data/demo.2_19.kmers > demo_data/demo.2_19.counts

Options:
  -h, --help            show this help message and exit
  -d DEPTH_PARAM, --depth=DEPTH_PARAM
                        The minimum times a k-mer is found to be reported in
                        the output.(default: 0)
```

It generates a file with 4 columns:

```
@Target Position    kc	dp
chr1    11865	1	14
chr1    11866	1	14
chr1    11867	1	17
chr1    11868	1	21
chr1    11869	1	20
chr1    11870	1	16
```

* @ is a symbol used to differentiate the header from the other rows
* Target: the chromosome or contig of the current position.
* Position: the current position within the target.
* kc: the number of different kmers found at the current position.
* dp: depth, total number of kmers found at the current position.

## Translate the position-based kmer counts to intervals of constant kmer count

This is done by running the script **kmeleon_intervals.py** and the next parameters:

```
Usage: kmeleon_intervals.py [OPTIONS] COUNTS_FILE
The COUNTS file has the format of the output of kmeleon_count.py
Note that this software outputs to stderr and stdout.

typical command:
kmeleon_intervals.py -m binary -s 50 demo_data/demo.2_19.counts > demo_data/demo.2_19.intervals

Options:
  -h, --help            show this help message and exit
  -m MODE_PARAM, --mode=MODE_PARAM
                        How intervals are computed. Either "windows", "raw",
                        "binary" or "smooth". (default: binary). These modes
                        are explained in the README file.
  -s SPAN_PARAM, --span=SPAN_PARAM
                        The minimum span of bases to be considered a whole
                        interval. Not used with option -w. (default: 50)
  --min-kc=MIN_KC_PARAM
                        Minimum kmer_count in a interval to be reported. Set
                        to 0 to report all the intervals. (default: 0)
  -d, --diploid         Requires -m binary. When -d is set, positions with
                        kmer_count=2 will be treated as regular heterozygous
                        variants and the interval value  will be set to 0 as
                        with kmer_count==1. If -d is not set, intervals with
                        kmer_count==2 will be set to 1, as with kmer_count>2.
  -w WINDOW_PARAM, --window=WINDOW_PARAM
                        Requires -m windows. The size of the window for which
                        mean k-mer counts will be processed. (default: 500)
```
kmeleon_intervals.py reads a file with kmer counts to create intervals of constant or homogeneous kmer count.
It generates a file in BED-like format without header and with 4 columns:

```
chr1	41868	41871	1
chr1	41911	41925	1
chr1	41977	41982	1
chr1	41985	42023	1
chr1	42054	42056	1
chr1	42169	42231	1
chr1	42274	42348	1
chr1	42349	42375	2
chr1	42376	42390	1
chr1	42391	42420	2
```

* 1st column is the target name
* 2nd column is the starting position of the interval. Note that this field is 0-based, as in BED format
* 3rd column is the ending position of the interval
* 4th column shows the kmer count, depending on the mode used to compute the intervals (see below).

**What do the different modes do?**

In **"raw" mode** (-m raw) the intervals and their associated values correspond directly to the kmer counts found.

For example, the kmer counts in:

```
@Target Position    kc	dp
chr1    11865	1	14
chr1    11866	1	14
chr1    11867	2	17
chr1    11868	2	21
chr1    11869	1	20
chr1    11870	1	16
chr1    11880	3	20
chr1    11881	3	16
```

would translate into:

```
chr1	11864	11866	1
chr1	11866	11868	2
chr1	11868	11870	1
chr1	11880	11881	3

```

In **"binary" mode** (-m binary) the values of the intervals reflect whether there is kmer count > 1 ("1") or not ("0").
Therefore, sucessive bases with kmer count > 1 will be considered of the same interval even if their kmer counts are different.
Note that if the -d (--diploid) parameter is set, the values of the intervals will be "1" for kmer_count > 2 and
"0" for kmer_count==1 or kmer_count==2. The -d option should be set when we expect regular heterozygous variants in our data,
and it should not be set when we do not expect heterozygous variants (haploid species, pure lines, etc).
Thus, "binary" mode will generate longer intervals than "raw" mode, in general.
For example, the previous kmer counts would be:

```
chr1	11864	11866	0
chr1	11866	11868	1
chr1	11868	11870	0
chr1	11880	11881	1

```

Note that, for both "raw" and "binary" modes, the --span (-s) param will set the minimum length
of any interval to be reported in the output. No filtering will be applied when it is set to 0 (-s 0).

In **"smooth" mode** (-m smooth) the intervals are usually handled as in "raw" mode.
However, when an interval is shorter than the --span parameter it is not filtered out directly.
Instead, the value of the interval is adjusted according to adjacent intervals:
* If the adjacent intervals are of smaller k-mer count
the interval count will be changed to that smaller k-mer count.
* If the adjacent intervals are of greater k-mer count
the interval count will be changed to that greater k-mer count.

Note that first of all, the intervals are internally sorted,
and thus those intervals with largest k-mer count are processed first.

For example, with --span 50, the kmer counts in:

```
@Target Position    kc	dp
chr1    11865	1	14
chr1    11866	1	14
chr1    11867	2	17
...
chr1    11887	2	21
chr1    11888	1	20
...
chr1    11988	1	16
chr1    11989	3	20
...
chr1    12989	3	16
```

would translate into: (the positions with kc = 2 have been smoothed to kc = 1)

```
chr1	11864	11988	1
chr1	11988	12989	3

```

In **"windows" mode** (-m windows) each reported interval corresponds to a sliding window of fixed length.
In this mode the --span parameter is ignored. Instead, the --window (-w) parameter can be used to determine
the size of sliding windows to use. The value returned for each window is the mean of kmer counts for all the
bases within that window.

For example, with --window 500, the kmer counts in:

```
@Target Position    kc	dp
chr1    11865	1	14
chr1    11866	1	14
chr1    11867	2	17
...
chr1    11887	2	21
chr1    11888	1	20
...
chr1    11988	1	16
chr1    11989	3	20
...
chr1    12989	3	16
```

would translate into:

```
chr1	11864	12364	1.7
chr1	12364	12864	3
chr1    12864   13364   1.4

```

## Join the kmer counts of all samples in a single table

This is done by running the script **kmeleon_table.py**, which has the next parameters:

```
Usage: kmeleon_table.py [OPTIONS] SAMPLES_FILE

The SAMPLES_FILE has a row for each sample, with tab-separated columns: 1st column is the sample name, 2nd column is the path to the counts or intervals file, which has the format of the output of "kmeleon_count.py" or "kmeleon_intervals -m windows", respectively.
Note that this software outputs to stderr and stdout.

typical command: kmeleon_table.py -D demo_data demo_data/demo_counts.list > demo_data/demo_counts.table
typical command: kmeleon_table.py -w -D demo_data demo_data/demo_windows.list > demo_data/demo_windows.table

Options:
  -h, --help            show this help message and exit
  -D DIR_PARAM, --DIR=DIR_PARAM
                        The directory where the files with counts or intervals
                        are located. Note that this is unnecessary if the
                        directory has been included in every path within the
                        SAMPLES_FILE (default: ./)
  -w                    If this flag is set, the input data should be
                        intervals of fixed length and positions for all
                        samples. i.e. as those produced by
                        "kmeleon_intervals.py -m windows". If -w is not set,
                        the input data expected will be kmer counts as those
                        produced with "kmeleon_count.py".
```

An example of a SAMPLES_FILE can be found at [demo_data/demo_counts.list](demo_data/demo_counts.list) for counts data
or at [demo_data/demo_counts.list](demo_data/demo_counts.list) for intervals data.

An example of the output when joining counts data:

```
@Target Position	sample1	sample2	...	sampleN
chr1    11865   1	1	...	1
chr1    11866   1	2	...	1
chr1    11867   1	2	...	1
chr1    11868   3	2	...	1
chr1    11869   3	1	...	0
chr1    11870   3	1	...	0
```

* @ is a symbol used to differentiate the header from the other rows
* Target: the chromosome or contig of the current position.
* Position: the current position within the target.
* sample1...sampleN: the kmers_count (number of different kmers found) in the current position for each sample

An example of the output when joining intervals (-w option) data:

```
@Target Start   END	sample1	sample2	...	sampleN
chr1    11500   12000	1.6	...	1.0
chr1    12000   12500   1.0	2.4	...	1.0
chr1    12500   13000   1.0	2.6	...	1.1
```

## Obtain a matrix of pairwise distances or similarities based on intersection of intervals between samples

This is done by running the script kmeleon_intersection_matrix.sh and the next parameters:

`kmeleon_intersection_matrix.sh samplesfiles_list_file distances`

- samplesfiles_list_file: a file with a list of duples of sample (1st column) and the file with kmercount intervals of such sample (2nd column).
- distances: whether to report similarities (no) or distances (yes).

For example:

`kmeleon_intersection_matrix.sh my_samples_files yes`

it generates a file with a header and an n*n matrix, where n = number of samples.:

```
sample1 sample2 ...     sampleN
0   0.6       ...     0.4
0.5   0	...     0.5
0.3	0.4	...     0
```
NOTE that the values of the matrix will be distances, if distances=yes, as shown above, or similarities, if distances=no.
