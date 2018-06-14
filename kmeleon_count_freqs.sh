## CPCantalapiedra 2018
##

samplename="$1";
samplefile="$2";
outdir="$3";
removetmp="$4"; # This option makes sense only when "$samplefile" is a .gz file
# The uncompressed file will be removed if $removetmp=="1"

tempdir="$outdir";

printf "kmeleon_count_freqs.sh $samplename $samplefile $outdir ...\n";

ext=$(echo "$samplefile" | sed 's#.*\.##g');

# If the file name ends in "gz", uncompress it
#
if [ "$ext" == "gz" ]; then
	uncompfile=$(echo "$samplefile" | sed 's#\.gz##');
	if [ ! -f "$uncompfile" ]; then
		gunzip -c "$samplefile" > "$uncompfile";
	fi;
	targetfile="$uncompfile";
else
	targetfile="$samplefile";
fi;

####### Obtain times each kc value is observed
#######

# header
printf "kc\tcount\n" > "$outdir"/"$samplename".kc.stats.out;

# rows

cat "$targetfile" | awk '{if (NR==1) next; print $0}' | 
cut -f 3 | sort -T "$tempdir" -n | uniq -c | awk '{print $2"\t"$1}' | 
sort -k1,1n \
>> "$outdir"/"$samplename".kc.stats.out

# If an uncompressed file was created, remove it
#
if [ "$removetmp" == "1" ]; then
if [ "$ext" == "gz" ]; then
	if [ -f "$uncompfile" ]; then
		rm "$uncompfile";
	fi;
fi;
fi;

printf "END kmeleon_count_freqs.sh $samplename $samplefile $outdir .\n";

## END
