#! /usr/bin/Rscript

## CPCantalapiedra 2018

args <- commandArgs(trailingOnly = TRUE);

datafile=args[1];
outfile=args[2];

datatable <- read.csv(datafile, sep="\t", header=T);

datatable <- datatable[order(datatable$frackc_gt_1),]

numgenotypes=nrow(datatable);

#head(datatable,3)

################################

origtable = datatable
meds = c()

# kc1
mim = min(datatable$frackc1)
mam = max(datatable$frackc1)
datatable$frackc1 = (datatable$frackc1 - mim) / (mam - mim)
meds = c(meds, median(datatable$frackc1))

# kc_gt_1
mim = min(datatable[!is.na(datatable$frackc_gt_100),]$frackc_gt_1)
mam = max(datatable[!is.na(datatable$frackc_gt_100),]$frackc_gt_1)
datatable$frackc_gt_1 = (datatable$frackc_gt_1 - mim) / (mam - mim)
meds = c(meds, median(datatable$frackc_gt_1))

# kc2
mim = min(datatable[!is.na(datatable$frackc_gt_100),]$frackc2)
mam = max(datatable[!is.na(datatable$frackc_gt_100),]$frackc2)
datatable$frackc2 = (datatable$frackc2 - mim) / (mam - mim)
meds = c(meds, median(datatable$frackc2))

# kc3
mim = min(datatable[!is.na(datatable$frackc_gt_100),]$frackc3)
mam = max(datatable[!is.na(datatable$frackc_gt_100),]$frackc3)
datatable$frackc3 = (datatable$frackc3 - mim) / (mam - mim)
meds = c(meds, median(datatable$frackc3))

# kc4
mim = min(datatable[!is.na(datatable$frackc_gt_100),]$frackc4)
mam = max(datatable[!is.na(datatable$frackc_gt_100),]$frackc4)
datatable$frackc4 = (datatable$frackc4 - mim) / (mam - mim)
meds = c(meds, median(datatable$frackc4))

# kc5
mim = min(datatable[!is.na(datatable$frackc_gt_100),]$frackc5)
mam = max(datatable[!is.na(datatable$frackc_gt_100),]$frackc5)
datatable$frackc5 = (datatable$frackc5 - mim) / (mam - mim)
meds = c(meds, median(datatable$frackc5))

# kc6_10
mim = min(datatable[!is.na(datatable$frackc_gt_100),]$frackc6_10)
mam = max(datatable[!is.na(datatable$frackc_gt_100),]$frackc6_10)
datatable$frackc6_10 = (datatable$frackc6_10 - mim) / (mam - mim)
meds = c(meds, median(datatable$frackc6_10))

# kc11_100
mim = min(datatable[!is.na(datatable$frackc_gt_100),]$frackc11_100)
mam = max(datatable[!is.na(datatable$frackc_gt_100),]$frackc11_100)
datatable$frackc11_100 = (datatable$frackc11_100 - mim) / (mam - mim)
meds = c(meds, median(datatable$frackc11_100))

# kc_gt_100
mim = min(datatable[!is.na(datatable$frackc_gt_100),]$frackc_gt_100)
mam = max(datatable[!is.na(datatable$frackc_gt_100),]$frackc_gt_100)
datatable$frackc_gt_100 = (datatable$frackc_gt_100 - mim) / (mam - mim)
meds = c(meds, median(datatable$frackc_gt_100))

pointcex = 2;
num_decs = 2;

png(outfile, width=1100, height=1200)
par(mar=c(8,4,1,1)+.1)
par(mgp = c(3,1,0))
plot(1, type="n", yaxt = "n", xaxt = "n", ylab = "", xlab = "", 
            ylim=c(0,max(datatable$frackc1)), xlim=c(1,10), las=2)

title(ylab="% of mapped", line=1.5, cex.lab=pointcex + 1.0)
title(xlab="kmer count", line=5, cex.lab=pointcex + 1.0)

f_add_kc <- function(kcdata, origdata, order){
    points(rep(order, numgenotypes), kcdata, pch = 19, cex = pointcex)
    mim = min(kcdata)
    mam = max(kcdata)
    med = median(kcdata)
    if (round(min(origdata)*100, num_decs) < 10) {
        origmim = format(round(min(origdata)*100, num_decs+1), nsmall = num_decs+1)
    } else {
        origmim = format(round(min(origdata)*100, num_decs), nsmall = num_decs)
    }
    if (round(max(origdata)*100, num_decs) < 10) {
        origmam = format(round(max(origdata)*100, num_decs+1), nsmall = num_decs+1)
    } else {
        origmam = format(round(max(origdata)*100, num_decs), nsmall = num_decs)
    }
    if (round(median(origdata)*100, num_decs) < 10){
        origmed = format(round(median(origdata)*100, num_decs+1), nsmall = num_decs+1)
    } else {
        origmed = format(round(median(origdata)*100, num_decs), nsmall = num_decs)
    }
    
    text(x=rep(order, 3)+0.05, y=c(mim, med, mam), labels=c(origmim, origmed, origmam), cex= pointcex+0.2, pos = 4)
    return();
}

tsssst <- f_add_kc(datatable$frackc1, origtable$frackc1, 1)
tsssst <- f_add_kc(datatable$frackc_gt_1, origtable$frackc_gt_1, 2)
tsssst <- f_add_kc(datatable$frackc2, origtable$frackc2, 3)
tsssst <- f_add_kc(datatable$frackc3, origtable$frackc3, 4)
tsssst <- f_add_kc(datatable$frackc4, origtable$frackc4, 5)
tsssst <- f_add_kc(datatable$frackc5, origtable$frackc5, 6)
tsssst <- f_add_kc(datatable$frackc6_10, origtable$frackc6_10, 7)
tsssst <- f_add_kc(datatable$frackc11_100, origtable$frackc11_100, 8)
tsssst <- f_add_kc(datatable$frackc_gt_100, origtable$frackc_gt_100, 9)

lines(meds, type = "l", col = "red", cex = pointcex, lwd=pointcex);
points(c(1:9), meds, pch = "-", cex = pointcex+4.5, col = "red")

lablist = c("1", ">1", "2", "3", "4", "5", "6-10", "11-100", ">100");
axis(1, at=1:9, cex.axis = pointcex+1.0, labels = FALSE) # labels=lablist)
text(x=1:9, y=par()$usr[3]-0.02,
labels = lablist, srt = 45, adj = 1, xpd = TRUE, cex = pointcex + 0.8)
tsssst <- dev.off()

## END