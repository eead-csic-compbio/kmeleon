## CPCantalapiedra 2018
## 
## A script to perform pairwise comparison of intervals of kmercount

samplesfileslist="$1"
distances="$2"; # no, pairwise intersection; yes, pairwise distances

# sample1 sample1.intervals.file
# sample2 sample2.intervals.file
# ...

## Header

(printf "sample\t";
cat "$samplesfileslist" | cut -f 1 |
while read sample; do
printf "$sample\t";
done;
printf "\n";

## Pairwise comparisons
cat "$samplesfileslist" | while read sampleline1; do
        sample1=$(echo "$sampleline1" | awk '{print $1}');
        printf "$sample1\t";

        intfile1=$(echo "$sampleline1" | awk '{print $2}');
        numint1=$(cat "$intfile1" | wc -l);

        cat "$samplesfileslist" | while read sampleline2; do
                sample2=$(echo "$sampleline2" | awk '{print $1}');
                intfile2=$(echo "$sampleline2" | awk '{print $2}');

                numint=$(bedops -i "$intfile1" "$intfile2" | wc -l);
                freq1=$(echo "" | awk '{print '"$numint"'/'"$numint1"'}');

                printf "$freq1\t";
        done;
printf "\n";
done) |
awk '
BEGIN{
distances='"$distances"';
}
{
if (NR==1) {
	for (i=2;i<=NF;i++) {
                if (i==NF) printf $i;
                else printf $i"\t";
        }
	printf "\n"; next;
}
for (i=2;i<=NF;i++) {x[NR-1,i-1]=$i;}
}END{
# sqrt to obtain the length of a row of	the matrix, since the length of	the matrix is rows*cols
for (i=1;i<=sqrt(length(x));i++) {
        for (j=1;j<=sqrt(length(x));j++) {
		if (distances=="no"){
			if (j==sqrt(length(x))) printf x[i,j];
			else printf x[i,j]"\t";
		} else { # distances=="yes"
	                simil=x[i,j]*x[j,i];
	                dist=1-simil;
	                if (j==sqrt(length(x))) printf dist;
	                else printf dist"\t";
		}
        }
	printf "\n";
};
}'

## END
