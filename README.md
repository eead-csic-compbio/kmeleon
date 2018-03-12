# kmeleon

A pipeline for k-mer count analysis of genomic diversity.

There are several steps to run a complete kmeleon analysis:
* [Extract the kmers from the mappings: kmeleon_extract.py](https://github.com/eead-csic-compbio/kmeleon#1st-extract-the-kmers-from-the-mappings)
* Count the number of kmers at each genomic position: kmeleon_count.py

## 1st) Extract the kmers from the mappings

This is done by running the script kmeleon_extract.py and the next parameters:


`kmeleon_extract.py target mappings kmer_size flush depths`

- target: for example the chromosome name or number.
- mappings: the file with mappings to process. A SAM or BAM file would be great.
- kmer_size: a number which indicates the size of sequences to search for. 50 for example.
- flush: a number which allows controlling how much memory is used by kmeleon while reading the SAM/BAM file. 10000 or 20000 is ok in general.
- depths: (Optional argument) if the word 'depths' is given, the depth of each kmer found will be also reported.

For example:

`kmeleon_extract.py chr1 mappings.sam 50 10000 depths`

it generates a file with 3 columns, as shown in its header:

```
@ Position	MD_Z	count
2917	51	1
2918	51	1
2919	51	1
2920	51	1
2921	51	1
2922	51	1
```

or

`kmeleon_extract.py chr1 mappings.sam 50 10000`

to obtain only the first 2 columns (and a lighter file):

```
@ Position	MD_Z
2917	51
2918	51
2919	51
2920	51
2921	51
2922	51
```

* @ is a symbol used to differentiate the header from the other rows
* Position: the position within the target (chromosome)
* MD_Z: a kmer found in such position, in SAM/BAM file MD_Z field format
