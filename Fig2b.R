library(tximport)
library(readr)
library(DESeq2)
library(ensembldb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(ggbeeswarm)
library(RColorBrewer)
library(pheatmap)
library(gplots)
dir =  "~/"

dat = read.csv("TPM_BATF2.csv",  header = T)

gene.name = c( "Arg1", "Ifih1", "Ifna1", "Ifna11", "Ifna12", "Ifna13", "Ifna14", "Ifna15", "Ifna16", "Ifna2", "Ifna4",
                            "Ifna5", "Ifna6", "Ifna7", "Ifna9", "Ifnb1" , "Mrc1", "Mx1", "Mx2", "Oasl1", "Tnf")
id = match(gene.name, dat[, 1])
sub.dat = dat[id, 2:13]
sub.mean = apply(sub.dat, 1, mean)
dat1 = sub.dat[, c(1, 2, 7, 8, 3, 4, 9, 10, 5, 6, 11, 12)]

a = (dat1 - sub.mean)/sub.sd
a = as.matrix(a)
rownames(a) = gene.name

rc <- rainbow(nrow(a), start = 0, end = .3)
cc <- rainbow(ncol(a), start = 0, end = .3)

a1 = a

id = a1 > 3
a1[id] = 3
id = a1 < -3
a1[id] = -3


dat2 = as.matrix(dat1)

logdat = read.csv("BATF2_KO_counts.csv",  header = T)

gene.name = c('Batf2', 'Tmem173', 'Mavs', 'Ifnb1', 'Ifna1',  'Ifna2', 'Ifna4', 'Ifna5', 'Ifna6', 'Ifna7', 'Ifna9', 'Ifna11', 'Ifna12',
                   'Ifna13',  'Ifna14',  'Ifna15', 'Ifna16', 'Mx1' , 'Mx2', 'Oasl1', 'Tnf', 'Arg1', 
                            'Mrc1')

id = match(gene.name, logdat[, 1])
sub.logdat = log2(logdat[id, 2:13] + 1)
sublog.mean = apply(sub.logdat, 1, mean)
logdat1 = sub.logdat[, c(1, 2, 7, 8, 3, 4, 9, 10, 5, 6, 11, 12)]

sub.logsd = apply(sub.logdat, 1, sd)

loga = (logdat1 - sublog.mean)/sub.logsd
loga = as.matrix(loga)
rownames(loga) = gene.name

pdf("heatmap_Zscore_logBATF2_V2.pdf", height = 8, width = 8)
colnames(loga) = c("WT NT", "WT NT", "KO NT", "KO NT", "WT ISD", "WT ISD", "KO ISD", "KO ISD", "WT poly IC", 
        "WT poly IC", "KO poly IC", "KO poly IC")
heatmap.2(loga, dendrogram='none', col = colorRampPalette(c('blue', 'green', 'white', 'yellow', 'red'))(256), 
		Colv = FALSE, density.info="none", Rowv = FALSE, trace="none"
		, lmat=rbind(c(4, 2), c(1, 3)), lhei=c(2, 8), lwid=c(4, 1))
dev.off()


