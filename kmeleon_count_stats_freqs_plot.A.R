#! /usr/bin/Rscript

## CPCantalapiedra 2018

args <- commandArgs(trailingOnly = TRUE);

datafile=args[1];
outfile=args[2];

datatable <- read.csv(datafile, sep="\t", header=T);

datatable <- datatable[order(datatable$frackc_gt_1),]

numgenotypes=nrow(datatable);

#head(datatable,3)

png(outfile, width=800, height=400);
par(mar=c(7,4,1,1)+.1)
par(mgp = c(3,1,0))

pointcex = 0.5;

plot(1, type="n", axes=FALSE, yaxt="n", xaxt = "n", ylab = "freq kc > 1", xlab = "",
#            xaxs="i", yaxs="i",
            xlim=c(1,numgenotypes), ylim=c(0,max(datatable$frackc_gt_1)), las=2)

points(x=1:numgenotypes, y=datatable$frackc_gt_1, pch = 19, cex = pointcex, col = "black")
lines(datatable$frackc_gt_1, type = "l", col = "black");

points(x=1:numgenotypes, y=datatable$frackc2, pch = 19, cex = pointcex - 0.1, col = "red")
lines(datatable$frackc2, type = "l", col = "red");

# lines connecting both series
show <- function(sample, a, b){
    lines(c(sample,sample), c(a, b), type = "l", col = "gray10", lty=3)
}
tssst <- sapply(1:nrow(datatable), function(i) show(i, datatable[i,"frackc2"], datatable[i,"frackc_gt_1"]))

axis(1, at=1:numgenotypes, labels=datatable$sample, las=3, cex.axis = 0.8)
ticks = seq(0,round(max(datatable$frackc_gt_1),2), 0.01);
axis(2, at=ticks, labels=format(ticks, nsmall=2), las=1)

# grid lines
yaxp <- par("yaxp")
abline(h=ticks, lty=6, col = "gray")

# medians
abline(h=median(datatable$frackc_gt_1), lty=1, col = "black")
abline(h=median(datatable$frackc2), lty=1, col = "red")

tsssst <- dev.off()

## END