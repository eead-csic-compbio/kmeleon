## CPCantalapiedra 2018
##

indir="$1";
sampleslist="$2";
#outfile="$3";

numsamples=$(cat "$sampleslist" | wc -l);

cat "$sampleslist" | cut -f 1 | while read sample; do
cat "$indir"/"$sample".kc.stats.out | awk '{if (NR==1) next; print $0}';
done | awk '{
kc=$1;
numsamples='"$numsamples"';
if (kc in kcs) kcs[kc]++; else kcs[kc]=1;
}END{
for (a in kcs) printf a"\t"kcs[a]"\t"kcs[a]/numsamples"\n"
}' | sort -k1,1n;
#}' | sort -k1,1n \
#> "$outfile"

## END
