## CPCantalapiedra 2018
## 
## A script which converts a file with kmer counts per position
## to kmer counts per interval
##
## If "binary" is set to 1 the intervals are just reported as
## having 1 kmer ("0") or more than 1 kmer ("1").
##
## Else, intervals are restricted to same kmer count, and thus
## a segmentation algorithm like those from CNV would improve
## this script.

sample="$1";
target="$2";
countsDIR="$3";
minspan="$4"; # minimum difference between start and end position to be considered a whole interval
kmercount="$5"; # all, uniq, multiple
binary="$6"; # no: the raw kmer count is used for intervals; yes: 0 if kmer_count=1, 1 if kmer_count>1

zcat "$countsDIR"/"$sample"."$target".counts.out.gz |
awk '
BEGIN{
binary='"$binary"';
pos=-1;count=-1;inipos=-1;
minspan='"$minspan"';
kmercount='"$kmercount"';
}
{
pos=$1;
if (binary=="no") {
	count=$2;
} else { # binary=="yes"
	if ($2==1) count=0;
	else count=1;
}
if (pos!=prevpos+1 || count!=prevcount) {
	if (inipos!=-1 && inipos<prevpos) {
		span=prevpos-inipos+1;
		if (kmercount=="all" || 
			(kmercount=="uniq" && binary=="no" && prevcount==1) ||
			(kmercount=="uniq" && binary=="yes" && prevcount==0) ||
			(kmercount=="multiple" && binary=="no" && prevcount>1) ||
			(kmercount=="multiple" && binary=="yes" && prevcount==1)
			) {
			if (span>=minspan){
				print "'"$full"'""\t"inipos"\t"prevpos"\t"prevcount; 
			}
		}
	}
	inipos=pos;
}
prevpos=pos; prevcount=count;
}
END{
if (inipos!=-1 && inipos<prevpos) {
	span=prevpos-inipos+1;
	if (kmercount=="all" ||	
		(kmercount=="uniq" && binary=="no" && prevcount==1) ||
		(kmercount=="uniq" && binary=="yes" && prevcount==0) ||
		(kmercount=="multiple" && binary=="no" && prevcount>1) ||
		(kmercount=="multiple" && binary=="yes" && prevcount==1)		
		)	{
		if (span>=minspan){
			print "'"$full"'""\t"inipos"\t"prevpos"\t"prevcount;
		}
	}
}
}
'

## END
