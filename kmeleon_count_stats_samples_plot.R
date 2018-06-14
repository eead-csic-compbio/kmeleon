#! /usr/bin/Rscript

## CPCantalapiedra 2018

args <- commandArgs(trailingOnly = TRUE);

datafile=args[1];
outfile=args[2];
column=args[3];

datatable <- read.csv(datafile, sep="\t", header=F, col.names=c("kc", "count", "percen"));

png(outfile);
par(mar=c(5,5,1,1)+.1)
if (column == "count"){
    plot(count~kc, datatable, pch = 19, cex = 0.4, cex.lab = 1.5, cex.axis = 1.4,  
    xlab="kmer count (kc)", ylab = "# of accessions", las=1);
    lines(stats::lowess(datatable$kc, datatable$count), type = "l", col = "red");
} else { # column == "percen"
    plot(percen~kc, datatable, pch = 19, cex = 0.4, cex.lab = 1.5, cex.axis = 1.4,  
    xlab="kmer count (kc)", ylab = "% of accessions", las=1);
    lines(stats::lowess(datatable$kc, datatable$percen), type = "l", col = "red");
}



tsssst <- dev.off()

## END