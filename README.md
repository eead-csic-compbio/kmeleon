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
* [Join the kmer counts of all samples in a single table: kmeleon_table.py](https://github.com/eead-csic-compbio/kmeleon#3rd-join-the-kmer-counts-of-all-samples-in-a-single-table)
* [Translate the position-based kmer counts to intervals of constant kmer count](https://github.com/eead-csic-compbio/kmeleon#4th-translate-the-position-based-kmer-counts-to-intervals-of-constant-kmer-count)
* [Obtain a matrix of pairwise distances or similarities based on intersection of intervals between samples](https://github.com/eead-csic-compbio/kmeleon#5th-obtain-a-matrix-of-pairwise-distances-or-similarities-based-on-intersection-of-intervals-between-samples)

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
                        "all" to process all the mappings.(default: "all".
  --start=START_PARAM   The -t target parameter is required. Starting
                        basepairs position to process within the given
                        target.(default: -1)
  --end=END_PARAM       The -t target parameter is required. Ending basepairs
                        position to process within the given target.(default:
                        -1)
  -k KMER_PARAM, --kmer=KMER_PARAM
                        A number, which will translate to "k"=number+1, to
                        parse fragments (k-mers) of length "k".(default: 50)
  -d DEPTH_PARAM, --depth=DEPTH_PARAM
                        The minimum times a k-mer is found to be reported in
                        the output.(default: 0)
  -b BAM_PARAM, --bam=BAM_PARAM
                        A BAM file to process.Either the -b or the -s option
                        is required.
  -s SAM_PARAM, --sam=SAM_PARAM
                        A SAM file to process.Either the -s or the -b option
                        is required.
```

It generates a file with 4 columns, as shown in its header:

```
@Target Position	kmer(MD_Z)	depth
chr1   2917	51	1
chr1   2918	51	1
chr1   2919	51	1
chr1    2920	51	1
chr1    2921	51	1
chr1    2922	51	1
```

* @ is a symbol used to differentiate the header from the other rows
* Position: the position within the target (chromosome)
* kmer(MD_Z): a kmer found in such position, in SAM/BAM file MD_Z field format
* depth: the times that this kmer has been observed at this position.

## 2nd) Count the number of kmers at each genomic position

This is done by running the script **kmeleon_count.py**, which has the next parameters:

```
Usage: kmeleon_count.py [OPTIONS] KMERS_FILE
Note that this software outputs to stderr and stdout.

typical command: kmeleon_count.py -d 4 demo_data/demo.2_19.kmers > demo_data/demo.2_19.counts

Options:
  -h, --help            show this help message and exit
  -d DEPTH_PARAM, --depth=DEPTH_PARAM
                        The minimum times a k-mer is found to be reported in
                        the output.(default: 0)
```

It generates a file with 3 columns:

```
@Target Position    kmers_count
chr1    11865	1
chr1    11866	1
chr1    11867	1
chr1    11868	1
chr1    11869	1
chr1    11870	1
```

* @ is a symbol used to differentiate the header from the other rows
* Target: the chromosome or contig of the current position.
* Position: the current position within the target.
* kmers_count(MD_Z): the number of different kmers found a the current position.

## 3rd) Join the kmer counts of all samples in a single table

This is done by running the script kmeleon_table.py, which has the next parameters:

```
Usage: kmeleon_table.py [OPTIONS] SAMPLES_COUNTS_FILE
The SAMPLES_COUNTS_FILE has a row for each sample, with tab-separated columns:1st column is the sample name, 2nd column is the path to the counts file,which has the format of the output of kmeleon_count.py
Note that this software outputs to stderr and stdout.

typical command: kmeleon_table.py -D demo_data demo_data/demo_samples_list > demo_data/demo_counts.table

Options:
  -h, --help            show this help message and exit
  -D DIR_PARAM, --DIR=DIR_PARAM
                        The directory where the files with counts for the
                        samples are located.Note that this is unnecessary if
                        the directory has been included in every path within
                        the SAMPLES_COUNTS_FILES(default: ./)
```

An example of a SAMPLES_COUNTS_FILE can be found at [demo_data/demo_samples_list](demo_data/demo_samples_list)
It generates a file with a header and n+1 columns, where n = number of samples.:

```
@Target Position	sample1	sample2	...	sampleN
11865   1	1	...	1
11866   1	2	...	1
11867   1	2	...	1
11868   3	2	...	1
11869   3	1	...	0
11870   3	1	...	0
```

* @ is a symbol used to differentiate the header from the other rows
* Target: the chromosome or contig of the current position.
* Position: the current position within the target.
* sample1...sampleN: the kmers_count (number of different kmers found) in the current position for each sample

## 4th) Translate the position-based kmer counts to intervals of constant kmer count

This is done by running the script kmeleon_intervals.sh and the next parameters:

`kmeleon_intervals.sh sample target counts_DIR minspan kmercount binary`

- sample: the sample to process.
- target: for example the chromosome name or number.
- counts_DIR: a path to the directory  where the files with counts (from kmeleon count) can be found. Note that the files in this folder must be in '.gz' format, \
and will have the next filename format:
`counts_DIR/sample_name.target.counts.out.gz`
- minspan: minimum length (consecutive genome nucleotides) of the resulting interval to be reported as such. E.g. 50
- kmercount: whether report all intervals (all), only those with kmercount=1 (uniq), or those with kmercount>1 (multiple)
- binary: whether the raw count value is used to compute intervals (no), or just differentiating kmercount=1 of kmercount>1 (yes)

For example:

`kmeleon_intervals.sh sample_name chr1H_part1 counts 50 multiple yes`

will report all the intervals (kmercount=all), joining all positions with kmercount>1 in intervals (binary=yes), and at least a length of 50 nts (minspan=50)

it generates a file in BED-like format without header and with 4 columns:

```
chr1H_part1	41868	41871	0
chr1H_part1	41911	41925	0
chr1H_part1	41977	41982	0
chr1H_part1	41985	42023	0
chr1H_part1	42054	42056	0
chr1H_part1	42169	42231	0
chr1H_part1	42274	42348	0
chr1H_part1	42349	42375	1
chr1H_part1	42376	42390	0
chr1H_part1	42391	42420	1
```

* 1st column is the target name@ is a symbol used to differentiate the header from the other rows
* 2nd column is the starting position of the interval. Note that this field is 0-based, as it follows the BED format
* 3rd column is the ending position of the interval
* 4th column shows whether kmercount=1 (0) or kmercount>1 (1)

NOTE that the values of the 4th column would show the actual kmercount with the option binary=no

## 5th) Obtain a matrix of pairwise distances or similarities based on intersection of intervals between samples

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
