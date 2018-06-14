## CPCantalapiedra 2018
##

indir="$1";
sampleslist="$2";
outfile="$3";

#bamsdir="$1";
#sampleslist="$2";
#statsdir="$3";

(printf "sample\tgenomesize\tmapped\tfracmapped\tkc1\tfrackc1\t";
printf "kc_gt_1\tfrackc_gt_1\tkc2\tfrackc2\tkc3\tfrackc3\t";
printf "kc4\tfrackc4\tkc5\tfrackc5\tkc6_10\tfrackc6_10\t";
printf "kc11_100\tfrackc11_100\tkc_gt_100\tfrackc_gt_100\n";

cat "$sampleslist" | while read line; do

sample=$(echo "$line" | awk '{print $1}');
genomesize=$(echo "$line" | awk '{print $2}');

printf "$sample\t";

mappedsize=$(cat "$indir"/"$sample".kc.stats.out | awk '{if (NR==1) next; print $2}' | awk '{sum+=$0}END{print sum}')

if [ ! "$genomesize" ]; then
printf "NA\t$mappedsize\tNA\t";
else
mappedfrac=$(echo "" | awk '{print '"$mappedsize"'/'"$genomesize"'}');
printf "$genomesize\t$mappedsize\t$mappedfrac\t";
fi;

# kc==1
kcsize=$(cat "$indir"/"$sample".kc.stats.out | awk '{if (NR==1) next; if ($1==1) print $2}' | awk '{sum+=$0}END{print sum}');

if [ "$kcsize" ]; then
printf "$kcsize\t"
kcfrac=$(echo "" | awk '{print '"$kcsize"'/'"$mappedsize"'}');
printf "$kcfrac\t";
else
printf  "0\t0\t";
fi;

# kc>1
kcsize=$(cat "$indir"/"$sample".kc.stats.out | awk '{if (NR==1) next; if ($1>1) print $2}' | awk '{sum+=$0}END{print sum}');
if [ "$kcsize" ]; then
printf "$kcsize\t"
kcfrac=$(echo "" | awk '{print '"$kcsize"'/'"$mappedsize"'}');
printf "$kcfrac\t";
else
printf 	"0\t0\t";
fi;

# kc==2
kcsize=$(cat "$indir"/"$sample".kc.stats.out | awk '{if (NR==1) next; if ($1==2) print $2}' | awk '{sum+=$0}END{print sum}');
if [ "$kcsize" ]; then
printf "$kcsize\t"
kcfrac=$(echo "" | awk '{print '"$kcsize"'/'"$mappedsize"'}');
printf "$kcfrac\t";
else
printf 	"0\t0\t";
fi;

# kc==3
kcsize=$(cat "$indir"/"$sample".kc.stats.out | awk '{if (NR==1) next; if ($1==3) print $2}' | awk '{sum+=$0}END{print sum}');
if [ "$kcsize" ]; then
printf "$kcsize\t"
kcfrac=$(echo "" | awk '{print '"$kcsize"'/'"$mappedsize"'}');
printf "$kcfrac\t";
else
printf 	"0\t0\t";
fi;

# kc==4
kcsize=$(cat "$indir"/"$sample".kc.stats.out | awk '{if (NR==1) next; if ($1==4) print $2}' | awk '{sum+=$0}END{print sum}');
if [ "$kcsize" ]; then
printf "$kcsize\t"
kcfrac=$(echo "" | awk '{print '"$kcsize"'/'"$mappedsize"'}');
printf "$kcfrac\t";
else
printf 	"0\t0\t";
fi;

# kc==5
kcsize=$(cat "$indir"/"$sample".kc.stats.out | awk '{if (NR==1) next; if ($1==5) print $2}' | awk '{sum+=$0}END{print sum}');
if [ "$kcsize" ]; then
printf "$kcsize\t"
kcfrac=$(echo "" | awk '{print '"$kcsize"'/'"$mappedsize"'}');
printf "$kcfrac\t";
else
printf 	"0\t0\t";
fi;

# 6<=kc<=10
kcsize=$(cat "$indir"/"$sample".kc.stats.out | awk '{if (NR==1) next; if ($1>=6 && $1<=10) print $2}' | awk '{sum+=$0}END{print sum}');
if [ "$kcsize" ]; then
printf "$kcsize\t"
kcfrac=$(echo "" | awk '{print '"$kcsize"'/'"$mappedsize"'}');
printf "$kcfrac\t";
else
printf 	"0\t0\t";
fi;

# 11<=kc<=100
kcsize=$(cat "$indir"/"$sample".kc.stats.out | awk '{if (NR==1) next; if ($1>=11 && $1<=100) print $2}' | awk '{sum+=$0}END{print sum}');
if [ "$kcsize" ]; then
printf "$kcsize\t"
kcfrac=$(echo "" | awk '{print '"$kcsize"'/'"$mappedsize"'}');
printf "$kcfrac\t";
else
printf 	"0\t0\t";
fi;

# kc > 100
kcsize=$(cat "$indir"/"$sample".kc.stats.out | awk '{if (NR==1) next; if ($1>100) print $2}' | awk '{sum+=$0}END{print sum}');
if [ "$kcsize" ]; then
printf "$kcsize\t"
kcfrac=$(echo "" | awk '{print '"$kcsize"'/'"$mappedsize"'}');
printf "$kcfrac";
else
printf 	"0\t0";
fi;

printf "\n";
done);
#done) > "$outfile";
#kc_stats."$bamsdir".kc.freqs

## END
