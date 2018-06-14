## CPCantalapiedra 2018
##

statsdir="$1";
sampleslist="$2";
outprefix="$3";

printf "Generating #samples=f(kc) table...\n";

outhist="$outprefix".kc.hist;

kmeleon_count_stats_samples.sh "$statsdir" "$sampleslist" > "$outhist".tab;
### it generates
# - A file "$outprefix".kc.hist.tab

printf "\t$outhist was created.\n"

printf "Generating #samples=f(kc) plot...\n";

kmeleon_count_stats_samples_plot.R "$outhist".tab "$outhist".png "percen"
### It generates
# - An R plot "$outprefix".kc.hist.png

printf "\t$outhist.png was created.\n"

####################

printf "Generating table with percent of kc>1 per sample...\n";

outfreqs="$outprefix".kc.freqs;

kmeleon_count_stats_freqs.sh "$statsdir" "$sampleslist" > "$outfreqs".tab;
### it generates
# - A file "$outprefix".kc.freqs.tab

printf "\t$outfreqs was created.\n"

printf "Generating plot with percent of kc>1 per sample...\n";

kmeleon_count_stats_freqs_plot.A.R "$outfreqs".tab "$outfreqs".A.png

printf "\t$outfreqs.A.png was created.\n"

printf "Generating percent of kc per sample distribution plot...\n";

kmeleon_count_stats_freqs_plot.B.R "$outfreqs".tab "$outfreqs".B.png

printf "\t$outfreqs.B.png was created.\n"

#./kc_stats_freqs_plot.R "$outfreqs" "$outfreqs".png
### it generates
# - An R plot "$outprefix".kc.freqs.A.png
# - An R plot "$outprefix".kc.freqs.B.png

printf "Finished.\n"

####################

## END
