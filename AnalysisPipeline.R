

# __     ___      _ _ _     __  __       _      ____  _   _    _                   
# \ \   / (_)_ __(_) (_)___|  \/  | __ _| | ___|  _ \| \ | |  / \   ___  ___  __ _ 
#  \ \ / /| | '__| | | / __| |\/| |/ _` | |/ _ \ |_) |  \| | / _ \ / __|/ _ \/ _` |           
#   \ V / | | |  | | | \__ \ |  | | (_| | |  __/  _ <| |\  |/ ___ \\__ \  __/ (_| |           
#    \_/  |_|_|  |_|_|_|___/_|  |_|\__,_|_|\___|_| \_\_| \_/_/   \_\___/\___|\__, |
#                                                                               |_|               
#  


## Load annotations
grpTrinotate = read.csv("Annotations/Trinotate_report_dvir1.06_subset.txt", header = T, sep = "\t", na.strings = ".", stringsAsFactors=FALSE)
amrTrinotate = read.csv("Annotations/Trinotate_report_amrTrin3_subset.txt", header = T, sep = "\t", na.strings = ".", stringsAsFactors=FALSE)
lumTrinotate = read.csv("Annotations/Trinotate_report_lumTrin3_subset.txt", header = T, sep = "\t", na.strings = ".", stringsAsFactors=FALSE)
novTrinotate = read.csv("Annotations/Trinotate_report_novTrin3_subset.txt", header = T, sep = "\t", na.strings = ".", stringsAsFactors=FALSE)
virTrinotate = read.csv("Annotations/Trinotate_report_virTrin3_subset.txt", header = T, sep = "\t", na.strings = ".", stringsAsFactors=FALSE)

gffRecord = read.table("Annotations/FBgn_ID_name_coordinates.txt", header = T)
TrinOrths = read.table("Annotations/Trin3.orthology.txt", header = T, sep = "\t")

## Load mel ortholog info
melOrths = read.table(file = "Other.Drosophilas/Dmel/mel_orths.txt", header = T)
melSFPs = read.table(file = "Other.Drosophilas/Dmel/ACPlist.Findlay.20130301.txt", header = T, sep = "\t")
melRPKM = read.table(file = "Other.Drosophilas/Dmel/mel.modEncode.RPKM.matrix", header = T, sep = "\t")


## Combine mel info with dvir1.06 Trinotate annotation
melOrthsAll = aggregate(mel_GeneSymbol~FBgn_ID, data = melOrths, toString)
Annots = merge(merge(melOrthsAll, grpTrinotate, all=TRUE), gffRecord, all=TRUE)

## Read in sample info:
grpSamples_data = read.table("Annotations/samples.txt", header=F, check.names=F, fill=T)
grpSamples_data = grpSamples_data[grpSamples_data[,2] != '',]
amrSamples_data = subset(grpSamples_data, grepl("amr", grpSamples_data$V1))
lumSamples_data = subset(grpSamples_data, grepl("lum", grpSamples_data$V1))
novSamples_data = subset(grpSamples_data, grepl("nov", grpSamples_data$V1))
virSamples_data = subset(grpSamples_data, grepl("vir", grpSamples_data$V1))

## Read in PAML and KaKs data 
tmp.FB.names = unique(subset(Annots, select=c("FBgn_ID", "FBtr_ID")))
paml.data = read.csv(file = "PAML.Files/PAML.branchSite.ALL.results.txt", header = T, sep = "\t")
paml.data = merge(tmp.FB.names, paml.data, all=T)
paml.data = merge(gffRecord, paml.data, all=T)
KaKs.data = read.csv(file = "PAML.Files/KaKs.ALL.results.txt", header = T, sep = "\t", check.names = F)
KaKs.data = merge(tmp.FB.names, KaKs.data, all=T)
KaKs.data = merge(gffRecord, KaKs.data, all=T)

#### Calaculate LRT, pValues and FDR for PAML data
paml.data$Damr_LRT = 2*(paml.data$Damr_brSt_H1 - paml.data$Damr_brSt_H0)
paml.data$Damr_pValue = pchisq(q = paml.data$Damr_LRT, df = 1, lower.tail = F)
paml.data$Damr_FDR = p.adjust(p = paml.data$Damr_pValue, method = "fdr")

paml.data$Dlum_LRT = 2*(paml.data$Dlum_brSt_H1 - paml.data$Dlum_brSt_H0)
paml.data$Dlum_pValue = pchisq(q = paml.data$Dlum_LRT, df = 1, lower.tail = F)
paml.data$Dlum_FDR = p.adjust(p = paml.data$Dlum_pValue, method = "fdr")

paml.data$Dnov_LRT = 2*(paml.data$Dnov_brSt_H1 - paml.data$Dnov_brSt_H0)
paml.data$Dnov_pValue = pchisq(q = paml.data$Dnov_LRT, df = 1, lower.tail = F)
paml.data$Dnov_FDR = p.adjust(p = paml.data$Dnov_pValue, method = "fdr")

paml.data$Dvir_LRT = 2*(paml.data$Dvir_brSt_H1 - paml.data$Dvir_brSt_H0)
paml.data$Dvir_pValue = pchisq(q = paml.data$Dvir_LRT, df = 1, lower.tail = F)
paml.data$Dvir_FDR = p.adjust(p = paml.data$Dvir_pValue, method = "fdr")

paml.data$DamrNov_LRT = 2*(paml.data$DamrNov_brSt_H1 - paml.data$DamrNov_brSt_H0)
paml.data$DamrNov_pValue = pchisq(q = paml.data$DamrNov_LRT, df = 1, lower.tail = F)
paml.data$DamrNov_FDR = p.adjust(p = paml.data$DamrNov_pValue, method = "fdr")

## Read in count and TPM matrices:
# all samples mapped to dvir1.06 transcriptome:
grpCountsMatrix = read.table("ExpressionData/genes_dvir1.06_allSamples.counts.matrix", header=T, row.names=1, com='', check.names=F)
grpTpmMatrix.notCrossNorm = read.table("ExpressionData/genes_dvir1.06_allSamples.TPM.not_cross_norm.counts_by_min_TPM", header = T)
grpTmmMatrix = read.table("ExpressionData/genes_dvir1.06_allSamples.TMM.EXPR.matrix", header=T, row.names=1, com='', check.names=F)
# D. americana mapped to Trinity transcripts
amrCountsMatrix = read.table("ExpressionData/amrTrin3.genes.counts.matrix", header=T, row.names=1, com='', check.names=F)
amrTpmMatrix.notCrossNorm = read.table("ExpressionData/amrTrin3.genes.TPM.not_cross_norm.counts_by_min_TPM", header = T)
amrTmmMatrix = read.table("ExpressionData/amrTrin3.genes.TMM.EXPR.matrix", header=T, row.names=1, com='', check.names=F)
# D. lummei mapped to Trinity transcripts
lumCountsMatrix = read.table("ExpressionData/lumTrin3.genes.counts.matrix", header=T, row.names=1, com='', check.names=F)
lumTpmMatrix.notCrossNorm = read.table("ExpressionData/lumTrin3.genes.TPM.not_cross_norm.counts_by_min_TPM", header = T)
lumTmmMatrix = read.table("ExpressionData/lumTrin3.genes.TMM.EXPR.matrix", header=T, row.names=1, com='', check.names=F)
# D. novamexicana mapped to Trinity transcripts
novCountsMatrix = read.table("ExpressionData/novTrin3.genes.counts.matrix", header=T, row.names=1, com='', check.names=F)
novTpmMatrix.notCrossNorm = read.table("ExpressionData/novTrin3.genes.TPM.not_cross_norm.counts_by_min_TPM", header = T)
novTmmMatrix = read.table("ExpressionData/novTrin3.genes.TMM.EXPR.matrix", header=T, row.names=1, com='', check.names=F)
# D. virilis mapped to Trinity transcripts
virCountsMatrix = read.table("ExpressionData/virTrin3.genes.counts.matrix", header=T, row.names=1, com='', check.names=F)
virTpmMatrix.notCrossNorm = read.table("ExpressionData/virTrin3.genes.TPM.not_cross_norm.counts_by_min_TPM", header = T)
virTmmMatrix = read.table("ExpressionData/virTrin3.genes.TMM.EXPR.matrix", header=T, row.names=1, com='', check.names=F)


########################################## QC ###########################################

## Barplot of gene counts by sample
libSizes = as.data.frame(colSums(grpCountsMatrix))
libSizes = cbind(sample = row.names(libSizes), libSizes)
row.names(libSizes)= NULL
colnames(libSizes) = c("sample", "Total_reads")
pdf("ManuscripPlots/Figure.S1.libSizes.pdf", width = 5.5, height = 3.2)
ggplot(libSizes, aes(sample, Total_reads)) + 
  geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0)) + 
  geom_hline(yintercept = 20000000)
dev.off()

## Boxplot of log10(TPM) across all samples
m.expData=melt(as.matrix(grpTmmMatrix))
colnames(m.expData) = c("gene_id", "replicate", "TPM")
m.expData.exp= cSplit(as.data.frame(m.expData), "replicate", "_")
m.expData=data.frame(m.expData, m.expData.exp$replicate_1, m.expData.exp$replicate_2, m.expData.exp$replicate_3)
colnames(m.expData) = c("gene_id", "replicate", "TPM", "species", "tissue", "rep_num")
m.expData$TPM = m.expData$TPM + 1
pdf("ManuscripPlots/Figure.S2.expData_by_sample.pdf", width = 5.5, height = 3.7)
ggplot(m.expData) + 
  geom_boxplot(aes(x = replicate, y = log10(TPM), fill = tissue, colour = species), size = 0.3, alpha = I(1/3)) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0)) + 
  scale_fill_hue(l = 50, h.start = 200)
dev.off()

 
# ## Estimate of the number of expressed genes (Brian Haas' method)
# extract the data between 10 TPM and 100 TPM
amr_filt_data = amrTpmMatrix.notCrossNorm[amrTpmMatrix.notCrossNorm[,1] > -100 & amrTpmMatrix.notCrossNorm[,1] < -10,]
lum_filt_data = lumTpmMatrix.notCrossNorm[lumTpmMatrix.notCrossNorm[,1] > -100 & lumTpmMatrix.notCrossNorm[,1] < -10,]
nov_filt_data = novTpmMatrix.notCrossNorm[novTpmMatrix.notCrossNorm[,1] > -100 & novTpmMatrix.notCrossNorm[,1] < -10,]
vir_filt_data = virTpmMatrix.notCrossNorm[virTpmMatrix.notCrossNorm[,1] > -100 & virTpmMatrix.notCrossNorm[,1] < -10,]
grp_filt_data = grpTpmMatrix.notCrossNorm[grpTpmMatrix.notCrossNorm[,1] > -100 & grpTpmMatrix.notCrossNorm[,1] < -10,]
# perform a linear regression on this filtered subset of the data
amr_fit = lm(amr_filt_data[,2] ~ amr_filt_data[,1])
lum_fit = lm(lum_filt_data[,2] ~ lum_filt_data[,1])
nov_fit = lm(nov_filt_data[,2] ~ nov_filt_data[,1])
vir_fit = lm(vir_filt_data[,2] ~ vir_filt_data[,1])
grp_fit = lm(grp_filt_data[,2] ~ grp_filt_data[,1])
print(amr_fit)
print(lum_fit)
print(nov_fit)
print(vir_fit)
print(grp_fit)
# plot it
amr_fit_plot=ggplot(amrTpmMatrix.notCrossNorm, aes(neg_min_tpm,num_features)) + geom_point() + scale_x_continuous(limits=c(-100,0)) + scale_y_continuous(limits=c(0,20000)) + geom_smooth(data=amr_filt_data, method = "lm") + geom_hline(yintercept = 13358, colour = "green") + ggtitle("amrTrinity")
lum_fit_plot=ggplot(lumTpmMatrix.notCrossNorm, aes(neg_min_tpm,num_features)) + geom_point() + scale_x_continuous(limits=c(-100,0)) + scale_y_continuous(limits=c(0,20000)) + geom_smooth(data=lum_filt_data, method = "lm") + geom_hline(yintercept = 14805, colour = "green") + ggtitle("lumTrinity")
nov_fit_plot=ggplot(novTpmMatrix.notCrossNorm, aes(neg_min_tpm,num_features)) + geom_point() + scale_x_continuous(limits=c(-100,0)) + scale_y_continuous(limits=c(0,20000)) + geom_smooth(data=nov_filt_data, method = "lm") + geom_hline(yintercept = 13646, colour = "green") + ggtitle("novTrinity")
vir_fit_plot=ggplot(virTpmMatrix.notCrossNorm, aes(neg_min_tpm,num_features)) + geom_point() + scale_x_continuous(limits=c(-100,0)) + scale_y_continuous(limits=c(0,20000)) + geom_smooth(data=vir_filt_data, method = "lm") + geom_hline(yintercept = 14616, colour = "green") + ggtitle("virTrinity")
grp_fit_plot=ggplot(grpTpmMatrix.notCrossNorm, aes(neg_min_tpm,num_features)) + geom_point() + scale_x_continuous(limits=c(-100,0)) + scale_y_continuous(limits=c(0,20000)) + geom_smooth(data=grp_filt_data, method = "lm") + geom_hline(yintercept = 9474, colour = "green") + ggtitle("dvir_1.06")
pdf("ManuscripPlots/Figure.S3.linearReg_of_expGenes.pdf", width = 5.7, height = 3.7)
plot_grid(amr_fit_plot, lum_fit_plot, nov_fit_plot, vir_fit_plot, grp_fit_plot,nrow = 2)
dev.off()


## calculate dispersion
d = DGEList(counts = grpCountsMatrix, group = grpSamples_data$V1)
d = calcNormFactors(d)
d = estimateCommonDisp(d)
d = estimateTagwiseDisp(d)
summary(d$tagwise.dispersion)
## Plot biological coefficient of variation
pdf("ManuscripPlots/Figure.S4.BCV.pdf", width = 5.5, height = 3.2)
plotBCV(d)
dev.off()
## Plot grouping of samples
pdf("ManuscripPlots/Figure.S5.MDS.pdf", width = 6.5, height = 5.6)
plotMDS(d, method = "bcv", col=as.numeric(d$samples$group))
dev.off()

## Plot sample correlation
data = log2(grpCountsMatrix+1)
data = as.matrix(data)
sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
sample_dist = as.dist(1-sample_cor)
hc_samples = hclust(sample_dist, method='complete')
pdf("ManuscripPlots/Figure.S6.sampleCorr.pdf", width = 8.6, height = 7.8)
heatmap.3(sample_cor, dendrogram='both', Rowv=as.dendrogram(hc_samples), Colv=as.dendrogram(hc_samples), col = colorpanel(75, '#dd70cd','black','#afc64f'), scale='none', symm=TRUE, key=TRUE,density.info='none', trace='none', symkey=FALSE, symbreaks=F, margins=c(10,10), cexCol=0.8, cexRow=0.8, cex.main=1.8, main=paste("sample correlation matrix\n(all replicates)"))
dev.off()

# normalize by DESeq method:
meta = data.frame(row.names=colnames(grpCountsMatrix), condition=grpSamples_data$V1)
grpCountData=round(grpCountsMatrix)
grpCountData_normByDESeq = newCountDataSet(grpCountData, meta)
grpCountData_normByDESeq = estimateSizeFactors(grpCountData_normByDESeq)
grpCountData_normByDESeq = data.frame(counts(grpCountData_normByDESeq, normalized=T))

# Cross-replicate comparison plots
amr.AG.1_vs_amr.AG.2 = MA_BPlot(grpCountData_normByDESeq, "amr_AG_1", "amr_AG_2")
amr.AG.1_vs_amr.AG.3 = MA_BPlot(grpCountData_normByDESeq, "amr_AG_1", "amr_AG_3")
amr.AG.2_vs_amr.AG.3 = MA_BPlot(grpCountData_normByDESeq, "amr_AG_2", "amr_AG_3")
amr.EB.1_vs_amr.EB.2 = MA_BPlot(grpCountData_normByDESeq, "amr_EB_1", "amr_EB_2")
amr.TS.1_vs_amr.TS.2 = MA_BPlot(grpCountData_normByDESeq, "amr_TS_1", "amr_TS_2")
amr.TS.1_vs_amr.TS.3 = MA_BPlot(grpCountData_normByDESeq, "amr_TS_1", "amr_TS_3")
amr.TS.2_vs_amr.TS.3 = MA_BPlot(grpCountData_normByDESeq, "amr_TS_2", "amr_TS_3")
pdf("ManuscripPlots/Figure.S7a.crossRep.Damr.pdf", width = 16.8, height = 12)
plot_grid(amr.AG.1_vs_amr.AG.2, amr.AG.1_vs_amr.AG.3, amr.AG.2_vs_amr.AG.3, amr.EB.1_vs_amr.EB.2, amr.TS.1_vs_amr.TS.2, amr.TS.1_vs_amr.TS.3, amr.TS.2_vs_amr.TS.3, ncol = 2, nrow = 4)
dev.off()

lum.AG.1_vs_lum.AG.2 = MA_BPlot(grpCountData_normByDESeq, "lum_AG_1", "lum_AG_2")
lum.AG.1_vs_lum.AG.3 = MA_BPlot(grpCountData_normByDESeq, "lum_AG_1", "lum_AG_3")
lum.AG.2_vs_lum.AG.3 = MA_BPlot(grpCountData_normByDESeq, "lum_AG_2", "lum_AG_3")
lum.EB.1_vs_lum.EB.2 = MA_BPlot(grpCountData_normByDESeq, "lum_EB_1", "lum_EB_2")
lum.TS.1_vs_lum.TS.2 = MA_BPlot(grpCountData_normByDESeq, "lum_TS_1", "lum_TS_2")
lum.TS.1_vs_lum.TS.3 = MA_BPlot(grpCountData_normByDESeq, "lum_TS_1", "lum_TS_3")
lum.TS.2_vs_lum.TS.3 = MA_BPlot(grpCountData_normByDESeq, "lum_TS_2", "lum_TS_3")
pdf("ManuscripPlots/Figure.S7b.crossRep.Dlum.pdf", width = 16.8, height = 12)
plot_grid(lum.AG.1_vs_lum.AG.2, lum.AG.1_vs_lum.AG.3, lum.AG.2_vs_lum.AG.3, lum.EB.1_vs_lum.EB.2, lum.TS.1_vs_lum.TS.2, lum.TS.1_vs_lum.TS.3, lum.TS.2_vs_lum.TS.3, ncol = 2, nrow = 4)
dev.off()

nov.AG.1_vs_nov.AG.2 = MA_BPlot(grpCountData_normByDESeq, "nov_AG_1", "nov_AG_2")
nov.AG.1_vs_nov.AG.3 = MA_BPlot(grpCountData_normByDESeq, "nov_AG_1", "nov_AG_3")
nov.AG.2_vs_nov.AG.3 = MA_BPlot(grpCountData_normByDESeq, "nov_AG_2", "nov_AG_3")
nov.EB.1_vs_nov.EB.2 = MA_BPlot(grpCountData_normByDESeq, "nov_EB_1", "nov_EB_2")
nov.TS.1_vs_nov.TS.2 = MA_BPlot(grpCountData_normByDESeq, "nov_TS_1", "nov_TS_2")
nov.TS.1_vs_nov.TS.3 = MA_BPlot(grpCountData_normByDESeq, "nov_TS_1", "nov_TS_3")
nov.TS.2_vs_nov.TS.3 = MA_BPlot(grpCountData_normByDESeq, "nov_TS_2", "nov_TS_3")
pdf("ManuscripPlots/Figure.S7c.crossRep.Dnov.pdf", width = 16.8, height = 12)
plot_grid(nov.AG.1_vs_nov.AG.2, nov.AG.1_vs_nov.AG.3, nov.AG.2_vs_nov.AG.3, nov.EB.1_vs_nov.EB.2, nov.TS.1_vs_nov.TS.2, nov.TS.1_vs_nov.TS.3, nov.TS.2_vs_nov.TS.3, ncol = 2, nrow = 4)
dev.off()

vir.AG.1_vs_vir.AG.2 = MA_BPlot(grpCountData_normByDESeq, "vir_AG_1", "vir_AG_2")
vir.AG.1_vs_vir.AG.3 = MA_BPlot(grpCountData_normByDESeq, "vir_AG_1", "vir_AG_3")
vir.AG.2_vs_vir.AG.3 = MA_BPlot(grpCountData_normByDESeq, "vir_AG_2", "vir_AG_3")
vir.EB.1_vs_vir.EB.2 = MA_BPlot(grpCountData_normByDESeq, "vir_EB_1", "vir_EB_2")
vir.TS.1_vs_vir.TS.2 = MA_BPlot(grpCountData_normByDESeq, "vir_TS_1", "vir_TS_2")
vir.TS.1_vs_vir.TS.3 = MA_BPlot(grpCountData_normByDESeq, "vir_TS_1", "vir_TS_3")
vir.TS.2_vs_vir.TS.3 = MA_BPlot(grpCountData_normByDESeq, "vir_TS_2", "vir_TS_3")
pdf("ManuscripPlots/Figure.S7d.crossRep.Dvir.pdf", width = 16.8, height = 12)
plot_grid(vir.AG.1_vs_vir.AG.2, vir.AG.1_vs_vir.AG.3, vir.AG.2_vs_vir.AG.3, vir.EB.1_vs_vir.EB.2, vir.TS.1_vs_vir.TS.2, vir.TS.1_vs_vir.TS.3, vir.TS.2_vs_vir.TS.3, ncol = 2, nrow = 4)
dev.off()


## Filter count data by minimum count across ANY sample (200 in this case)
amr_max_gene_expr_per_row = apply(amrCountsMatrix, 1, max)
lum_max_gene_expr_per_row = apply(lumCountsMatrix, 1, max)
nov_max_gene_expr_per_row = apply(novCountsMatrix, 1, max)
vir_max_gene_expr_per_row = apply(virCountsMatrix, 1, max)
grp_max_gene_expr_per_row = apply(grpCountsMatrix, 1, max)
amrCountsMatrix.min200count = amrCountsMatrix[amr_max_gene_expr_per_row >= 200,,drop=F ]
lumCountsMatrix.min200count = lumCountsMatrix[lum_max_gene_expr_per_row >= 200,,drop=F ]
novCountsMatrix.min200count = novCountsMatrix[nov_max_gene_expr_per_row >= 200,,drop=F ]
virCountsMatrix.min200count = virCountsMatrix[vir_max_gene_expr_per_row >= 200,,drop=F ]
grpCountsMatrix.min200count = grpCountsMatrix[grp_max_gene_expr_per_row >= 200,,drop=F ]

## Filter data by minimum CPM across ANY sample
# d.byCPM = d
# head(cpm(d.byCPM))
# apply(d.byCPM$counts, 2, sum)
# keep = rowSums(cpm(d.byCPM)>10) >= 2
# countsMatrix.min10CPM = d.byCPM[keep,]
# dim(countsMatrix.min10CPM)


#########################################################################################
### Summary TPM table and matrix for gene level plots (includes all replicates) ######### 
# dvir1.06
grpTPM.tmp=grpTmmMatrix
colnames(grpTPM.tmp) = grpSamples_data$V1
m.grpTPM.tmp = as.data.frame(melt(as.matrix(grpTPM.tmp)))
m.grpTPM.tmp = cSplit(as.data.frame(m.grpTPM.tmp), "X2", "_")
m.grpTPM.tmp = data.frame(m.grpTPM.tmp$X1, m.grpTPM.tmp$X2_1, m.grpTPM.tmp$X2_2, m.grpTPM.tmp$value)
colnames(m.grpTPM.tmp) = c("FBgn_ID", "species", "tissue", "TPM")
m.grpTPM.tmp.c = summarySE(m.grpTPM.tmp, measurevar = "TPM", groupvars = c("FBgn_ID", "species", "tissue"))
# fbgn_to_geneName=subset(gffRecord, select=c(FBgn_ID, gene_name))
# TPMse = merge(fbgn_to_geneName, m.grpTPM.tmp.c, all=TRUE)
# write.table(x = TPMse, file = "ExpressionData/dvir1.06_TPMse.txt", quote = F, sep = "\t", row.names = F)
TPMse = read.table(file = "ExpressionData/dvir1.06_TPMse.txt", header = T, sep = "\t")
grpMeanTPMmatrix=cast(m.grpTPM.tmp.c, FBgn_ID~species+tissue, value ="TPM")
# # amrTrinity
# amrTPM.tmp=amrTmmMatrix
# colnames(amrTPM.tmp) = amrSamples_data$V1
# m.amrTPM.tmp = as.data.frame(melt(as.matrix(amrTPM.tmp)))
# m.amrTPM.tmp = cSplit(as.data.frame(m.amrTPM.tmp), "X2", "_")
# m.amrTPM.tmp=data.frame(m.amrTPM.tmp$X1, m.amrTPM.tmp$X2_1, m.amrTPM.tmp$X2_2, m.amrTPM.tmp$value)
# colnames(m.amrTPM.tmp) = c("trinity_id", "species", "tissue", "TPM")
# TPMse_amrTrin = summarySE(m.amrTPM.tmp, measurevar = "TPM", groupvars = c("trinity_id", "species", "tissue"))
TPMse_amrTrin = read.table("ExpressionData/amrTrin3_TPMse.txt", header = T, sep = "\t")
amrMeanTPMmatrix=cast(TPMse_amrTrin, trinity_id~species+tissue, value ="TPM")
# # lumTrinity
# lumTPM.tmp=lumTmmMatrix
# colnames(lumTPM.tmp) = lumSamples_data$V1
# m.lumTPM.tmp = as.data.frame(melt(as.matrix(lumTPM.tmp)))
# m.lumTPM.tmp = cSplit(as.data.frame(m.lumTPM.tmp), "X2", "_")
# m.lumTPM.tmp=data.frame(m.lumTPM.tmp$X1, m.lumTPM.tmp$X2_1, m.lumTPM.tmp$X2_2, m.lumTPM.tmp$value)
# colnames(m.lumTPM.tmp) = c("trinity_id", "species", "tissue", "TPM")
# TPMse_lumTrin = summarySE(m.lumTPM.tmp, measurevar = "TPM", groupvars = c("trinity_id", "species", "tissue"))
TPMse_lumTrin = read.table("ExpressionData/lumTrin3_TPMse.txt", header = T, sep = "\t")
lumMeanTPMmatrix=cast(TPMse_lumTrin, trinity_id~species+tissue, value ="TPM")
# # novTrinity
# novTPM.tmp=novTmmMatrix
# colnames(novTPM.tmp) = novSamples_data$V1
# m.novTPM.tmp = as.data.frame(melt(as.matrix(novTPM.tmp)))
# m.novTPM.tmp = cSplit(as.data.frame(m.novTPM.tmp), "X2", "_")
# m.novTPM.tmp=data.frame(m.novTPM.tmp$X1, m.novTPM.tmp$X2_1, m.novTPM.tmp$X2_2, m.novTPM.tmp$value)
# colnames(m.novTPM.tmp) = c("trinity_id", "species", "tissue", "TPM")
# TPMse_novTrin = summarySE(m.novTPM.tmp, measurevar = "TPM", groupvars = c("trinity_id", "species", "tissue"))
TPMse_novTrin = read.table("ExpressionData/novTrin3_TPMse.txt", header = T, sep = "\t")
novMeanTPMmatrix=cast(TPMse_novTrin, trinity_id~species+tissue, value ="TPM")
# # virTrinity
# virTPM.tmp=virTmmMatrix
# colnames(virTPM.tmp) = virSamples_data$V1
# m.virTPM.tmp = as.data.frame(melt(as.matrix(virTPM.tmp)))
# m.virTPM.tmp = cSplit(as.data.frame(m.virTPM.tmp), "X2", "_")
# m.virTPM.tmp=data.frame(m.virTPM.tmp$X1, m.virTPM.tmp$X2_1, m.virTPM.tmp$X2_2, m.virTPM.tmp$value)
# colnames(m.virTPM.tmp) = c("trinity_id", "species", "tissue", "TPM")
# TPMse_virTrin = summarySE(m.virTPM.tmp, measurevar = "TPM", groupvars = c("trinity_id", "species", "tissue"))
TPMse_virTrin = read.table("ExpressionData/virTrin3_TPMse.txt", header = T, sep = "\t")
virMeanTPMmatrix=cast(TPMse_virTrin, trinity_id~species+tissue, value ="TPM")


##########################################################################################
################# Remove "bad" replicates  for DE analysis ###############################
### Based on the QC analysis above, some replicates show inconsistencies that are likely due to cross tissue contamination
### during dissections. The DE analysis will exclude these replicates
## Define good replicates (propper replicate grouping and correlation)
amrGoodReps = c("amr_AG_1","amr_AG_3","amr_CR_2","amr_CR_3","amr_EB_1","amr_TS_2","amr_TS_3")
lumGoodReps = c("lum_AG_1","lum_AG_2","lum_AG_3","lum_CR_1","lum_CR_3","lum_EB_1","lum_EB_2","lum_TS_1","lum_TS_2","lum_TS_3")
novGoodReps = c("nov_AG_2","nov_AG_3","nov_CR_1","nov_CR_3","nov_EB_1","nov_EB_2","nov_TS_1","nov_TS_2","nov_TS_3")
virGoodReps = c("vir_AG_1","vir_AG_2","vir_AG_3","vir_CR_1","vir_CR_3","vir_EB_1","vir_TS_1","vir_TS_2","vir_TS_3")
grpGoodReps = c(amrGoodReps, lumGoodReps, novGoodReps, virGoodReps)
## Create counts matrix with good replicates only
amrCountsMatrix.min200count.BRR=subset(amrCountsMatrix.min200count, select=amrGoodReps)
lumCountsMatrix.min200count.BRR=subset(lumCountsMatrix.min200count, select=lumGoodReps)
novCountsMatrix.min200count.BRR=subset(novCountsMatrix.min200count, select=novGoodReps)
virCountsMatrix.min200count.BRR=subset(virCountsMatrix.min200count, select=virGoodReps)
grpCountsMatrix.min200count.BRR=subset(grpCountsMatrix.min200count, select=grpGoodReps)
## Create normalized TPM matrix with good replicates only
amrTmmMatrix.BRR=subset(amrTmmMatrix, select=amrGoodReps)
lumTmmMatrix.BRR=subset(lumTmmMatrix, select=lumGoodReps)
novTmmMatrix.BRR=subset(novTmmMatrix, select=novGoodReps)
virTmmMatrix.BRR=subset(virTmmMatrix, select=virGoodReps)
grpTmmMatrix.BRR=subset(grpTmmMatrix, select=grpGoodReps)
## Rename columns to keep replicate order
# count matrices
colnames(amrCountsMatrix.min200count.BRR) = c("amr_AG_1","amr_AG_2","amr_CR_1","amr_CR_2","amr_EB_1","amr_TS_1","amr_TS_2")
colnames(lumCountsMatrix.min200count.BRR) = c("lum_AG_1","lum_AG_2","lum_AG_3","lum_CR_1","lum_CR_2","lum_EB_1","lum_EB_2","lum_TS_1","lum_TS_2","lum_TS_3")
colnames(novCountsMatrix.min200count.BRR) = c("nov_AG_1","nov_AG_2","nov_CR_1","nov_CR_2","nov_EB_1","nov_EB_2","nov_TS_1","nov_TS_2","nov_TS_3")
colnames(virCountsMatrix.min200count.BRR) = c("vir_AG_1","vir_AG_2","vir_AG_3","vir_CR_1","vir_CR_2","vir_EB_1","vir_TS_1","vir_TS_2","vir_TS_3")
colnames(grpCountsMatrix.min200count.BRR) = c(colnames(amrCountsMatrix.min200count.BRR), colnames(lumCountsMatrix.min200count.BRR), colnames(novCountsMatrix.min200count.BRR), colnames(virCountsMatrix.min200count.BRR))
# TPM matrices
colnames(amrTmmMatrix.BRR) = colnames(amrCountsMatrix.min200count.BRR)
colnames(lumTmmMatrix.BRR) = colnames(lumCountsMatrix.min200count.BRR)
colnames(novTmmMatrix.BRR) = colnames(novCountsMatrix.min200count.BRR)
colnames(virTmmMatrix.BRR) = colnames(virCountsMatrix.min200count.BRR)
colnames(grpTmmMatrix.BRR) = colnames(grpCountsMatrix.min200count.BRR)

##############################################################################################
### Summary TPM table and matrix for gene level plots (includes good replicates only) ########
# dvir1.06
grpTPM2.tmp = grpTmmMatrix.BRR
colnames(grpTPM2.tmp) = c(amrGoodReps, lumGoodReps, novGoodReps, virGoodReps)
m.grpTPM2.tmp = as.data.frame(melt(as.matrix(grpTPM2.tmp)))
m.grpTPM2.tmp = cSplit(as.data.frame(m.grpTPM2.tmp), "X2", "_")
m.grpTPM2.tmp = data.frame(m.grpTPM2.tmp$X1, m.grpTPM2.tmp$X2_1,m.grpTPM2.tmp$X2_2, m.grpTPM2.tmp$value)
colnames(m.grpTPM2.tmp) = c("FBgn_ID", "species", "tissue", "TPM")
m.grpTPM2.tmp.c = summarySE(m.grpTPM2.tmp, measurevar = "TPM", groupvars = c("FBgn_ID", "species", "tissue"))
fbgn_to_geneName=subset(gffRecord, select=c(FBgn_ID, gene_name))
TPMseBRR = merge(fbgn_to_geneName, m.grpTPM2.tmp.c, all=TRUE)
grpMeanTPMmatrix.BRR=cast(m.grpTPM2.tmp.c, FBgn_ID~species+tissue, value ="TPM")
# # amrTrinity
# amrTPM2.tmp = amrTmmMatrix.BRR
# colnames(amrTPM2.tmp) = amrGoodReps
# m.amrTPM2.tmp = as.data.frame(melt(as.matrix(amrTPM2.tmp)))
# m.amrTPM2.tmp = cSplit(as.data.frame(m.amrTPM2.tmp), "X2", "_")
# m.amrTPM2.tmp = data.frame(m.amrTPM2.tmp$X1, m.amrTPM2.tmp$X2_1,m.amrTPM2.tmp$X2_2, m.amrTPM2.tmp$value)
# colnames(m.amrTPM2.tmp) = c("trinity_id", "species", "tissue", "TPM")
# TPMseBRR_amrTrin = summarySE(m.amrTPM2.tmp, measurevar = "TPM", groupvars = c("trinity_id", "species", "tissue"))
# write.table(x = TPMseBRR_amrTrin, file = "ExpressionData/amrTrin3_TPMseBRR.txt", quote = F, sep = "\t", row.names = F)
TPMseBRR_amrTrin=read.table(file = "ExpressionData/amrTrin3_TPMseBRR.txt", header = T, sep = "\t")
amrMeanTPMmatrix.BRR=cast(TPMseBRR_amrTrin, trinity_id~species+tissue, value ="TPM")
# # lumTrinity
# lumTPM2.tmp = lumTmmMatrix.BRR
# colnames(lumTPM2.tmp) = lumGoodReps
# m.lumTPM2.tmp = as.data.frame(melt(as.matrix(lumTPM2.tmp)))
# m.lumTPM2.tmp = cSplit(as.data.frame(m.lumTPM2.tmp), "X2", "_")
# m.lumTPM2.tmp = data.frame(m.lumTPM2.tmp$X1, m.lumTPM2.tmp$X2_1,m.lumTPM2.tmp$X2_2, m.lumTPM2.tmp$value)
# colnames(m.lumTPM2.tmp) = c("trinity_id", "species", "tissue", "TPM")
# TPMseBRR_lumTrin = summarySE(m.lumTPM2.tmp, measurevar = "TPM", groupvars = c("trinity_id", "species", "tissue"))
# write.table(x = TPMseBRR_lumTrin, file = "ExpressionData/lumTrin3_TPMseBRR.txt", quote = F, sep = "\t", row.names = F)
TPMseBRR_lumTrin=read.table(file = "ExpressionData/lumTrin3_TPMseBRR.txt", header = T, sep = "\t")
lumMeanTPMmatrix.BRR=cast(TPMseBRR_lumTrin, trinity_id~species+tissue, value ="TPM")
# # novTrinity
# novTPM2.tmp = novTmmMatrix.BRR
# colnames(novTPM2.tmp) = novGoodReps
# m.novTPM2.tmp = as.data.frame(melt(as.matrix(novTPM2.tmp)))
# m.novTPM2.tmp = cSplit(as.data.frame(m.novTPM2.tmp), "X2", "_")
# m.novTPM2.tmp = data.frame(m.novTPM2.tmp$X1, m.novTPM2.tmp$X2_1,m.novTPM2.tmp$X2_2, m.novTPM2.tmp$value)
# colnames(m.novTPM2.tmp) = c("trinity_id", "species", "tissue", "TPM")
# TPMseBRR_novTrin = summarySE(m.novTPM2.tmp, measurevar = "TPM", groupvars = c("trinity_id", "species", "tissue"))
# write.table(x = TPMseBRR_novTrin, file = "ExpressionData/novTrin3_TPMseBRR.txt", quote = F, sep = "\t", row.names = F)
TPMseBRR_novTrin=read.table(file = "ExpressionData/novTrin3_TPMseBRR.txt", header = T, sep = "\t")
novMeanTPMmatrix.BRR=cast(TPMseBRR_novTrin, trinity_id~species+tissue, value ="TPM")
# # virTrinity
# virTPM2.tmp = virTmmMatrix.BRR
# colnames(virTPM2.tmp) = virGoodReps
# m.virTPM2.tmp = as.data.frame(melt(as.matrix(virTPM2.tmp)))
# m.virTPM2.tmp = cSplit(as.data.frame(m.virTPM2.tmp), "X2", "_")
# m.virTPM2.tmp = data.frame(m.virTPM2.tmp$X1, m.virTPM2.tmp$X2_1,m.virTPM2.tmp$X2_2, m.virTPM2.tmp$value)
# colnames(m.virTPM2.tmp) = c("trinity_id", "species", "tissue", "TPM")
# TPMseBRR_virTrin = summarySE(m.virTPM2.tmp, measurevar = "TPM", groupvars = c("trinity_id", "species", "tissue"))
# write.table(x = TPMseBRR_virTrin, file = "ExpressionData/virTrin3_TPMseBRR.txt", quote = F, sep = "\t", row.names = F)
TPMseBRR_virTrin=read.table(file = "ExpressionData/virTrin3_TPMseBRR.txt", header = T, sep = "\t")
virMeanTPMmatrix.BRR=cast(TPMseBRR_virTrin, trinity_id~species+tissue, value ="TPM")

##########################################################
##########################################################
## Create specificity matrices (dvir1.06)

# 1. within species
#D.amr
amr.dvir1.06.MeanTPMmatrix = subset(grpMeanTPMmatrix.BRR, select=c("FBgn_ID", "amr_AG", "amr_CR", "amr_EB", "amr_TS"))
tmp.amr.dvir1.06.MeanTPMmatrix = amr.dvir1.06.MeanTPMmatrix
rownames(tmp.amr.dvir1.06.MeanTPMmatrix) = tmp.amr.dvir1.06.MeanTPMmatrix[,1]
tmp.amr.dvir1.06.MeanTPMmatrix[,1] = NULL
amr.dvir1.06_Specificity_table = YazSpecificity(tmp.amr.dvir1.06.MeanTPMmatrix)
amr.dvir1.06_Specificity_table = as.data.frame(amr.dvir1.06_Specificity_table)
rm(tmp.amr.dvir1.06.MeanTPMmatrix)
#D.lum
lum.dvir1.06.MeanTPMmatrix = subset(grpMeanTPMmatrix.BRR, select=c("FBgn_ID", "lum_AG", "lum_CR", "lum_EB", "lum_TS"))
tmp.lum.dvir1.06.MeanTPMmatrix = lum.dvir1.06.MeanTPMmatrix
rownames(tmp.lum.dvir1.06.MeanTPMmatrix) = tmp.lum.dvir1.06.MeanTPMmatrix[,1]
tmp.lum.dvir1.06.MeanTPMmatrix[,1] = NULL
lum.dvir1.06_Specificity_table = YazSpecificity(tmp.lum.dvir1.06.MeanTPMmatrix)
lum.dvir1.06_Specificity_table = as.data.frame(lum.dvir1.06_Specificity_table)
rm(tmp.lum.dvir1.06.MeanTPMmatrix)
#D.nov
nov.dvir1.06.MeanTPMmatrix = subset(grpMeanTPMmatrix.BRR, select=c("FBgn_ID", "nov_AG", "nov_CR", "nov_EB", "nov_TS"))
tmp.nov.dvir1.06.MeanTPMmatrix = nov.dvir1.06.MeanTPMmatrix
rownames(tmp.nov.dvir1.06.MeanTPMmatrix) = tmp.nov.dvir1.06.MeanTPMmatrix[,1]
tmp.nov.dvir1.06.MeanTPMmatrix[,1] = NULL
nov.dvir1.06_Specificity_table = YazSpecificity(tmp.nov.dvir1.06.MeanTPMmatrix)
nov.dvir1.06_Specificity_table = as.data.frame(nov.dvir1.06_Specificity_table)
rm(tmp.nov.dvir1.06.MeanTPMmatrix)
#D.vir
vir.dvir1.06.MeanTPMmatrix = subset(grpMeanTPMmatrix.BRR, select=c("FBgn_ID", "vir_AG", "vir_CR", "vir_EB", "vir_TS"))
tmp.vir.dvir1.06.MeanTPMmatrix = vir.dvir1.06.MeanTPMmatrix
rownames(tmp.vir.dvir1.06.MeanTPMmatrix) = tmp.vir.dvir1.06.MeanTPMmatrix[,1]
tmp.vir.dvir1.06.MeanTPMmatrix[,1] = NULL
vir.dvir1.06_Specificity_table = YazSpecificity(tmp.vir.dvir1.06.MeanTPMmatrix)
vir.dvir1.06_Specificity_table = as.data.frame(vir.dvir1.06_Specificity_table)
rm(tmp.vir.dvir1.06.MeanTPMmatrix)

# 2. across species and tissues:
#grp
tmp.grpMeanTPMmatrix = grpMeanTPMmatrix.BRR
rownames(tmp.grpMeanTPMmatrix) = tmp.grpMeanTPMmatrix[,1]
tmp.grpMeanTPMmatrix[,1] = NULL
dvir1.06_Specificity_table = YazSpecificity(tmp.grpMeanTPMmatrix)
dvir1.06_Specificity_table = as.data.frame(dvir1.06_Specificity_table)
rm(tmp.grpMeanTPMmatrix)

# 3. within tissues, across species:
#AG
tmp.AG.grpMeanTPMmatrix = subset(grpMeanTPMmatrix.BRR, select=c("FBgn_ID", "amr_AG", "lum_AG", "nov_AG", "vir_AG"))
rownames(tmp.AG.grpMeanTPMmatrix) = tmp.AG.grpMeanTPMmatrix[,1]
tmp.AG.grpMeanTPMmatrix[,1] = NULL
AG.dvir1.06_Specificity_table = YazSpecificity(tmp.AG.grpMeanTPMmatrix)
AG.dvir1.06_Specificity_table = as.data.frame(AG.dvir1.06_Specificity_table)
rm(tmp.AG.grpMeanTPMmatrix)
#EB
tmp.EB.grpMeanTPMmatrix = subset(grpMeanTPMmatrix.BRR, select=c("FBgn_ID", "amr_EB", "lum_EB", "nov_EB", "vir_EB"))
rownames(tmp.EB.grpMeanTPMmatrix) = tmp.EB.grpMeanTPMmatrix[,1]
tmp.EB.grpMeanTPMmatrix[,1] = NULL
EB.dvir1.06_Specificity_table = YazSpecificity(tmp.EB.grpMeanTPMmatrix)
EB.dvir1.06_Specificity_table = as.data.frame(EB.dvir1.06_Specificity_table)
rm(tmp.EB.grpMeanTPMmatrix)
#TS
tmp.TS.grpMeanTPMmatrix = subset(grpMeanTPMmatrix.BRR, select=c("FBgn_ID", "amr_TS", "lum_TS", "nov_TS", "vir_TS"))
rownames(tmp.TS.grpMeanTPMmatrix) = tmp.TS.grpMeanTPMmatrix[,1]
tmp.TS.grpMeanTPMmatrix[,1] = NULL
TS.dvir1.06_Specificity_table = YazSpecificity(tmp.TS.grpMeanTPMmatrix)
TS.dvir1.06_Specificity_table = as.data.frame(TS.dvir1.06_Specificity_table)
rm(tmp.TS.grpMeanTPMmatrix)

#############################################################################
#############################################################################
############  edgeR DE analysis I: Trinity ind. species #####################

### Identify tissue-biased genes for each species (>4-fold, < 0.001 FDR)

## D. americana
# Basic GLM set-up
amr.group = factor(c(1,1,2,2,3,4,4))
amr.design = model.matrix(~0+amr.group)
colnames(amr.design)=levels(factor(amrSamples_data$V1))
amrExpStd=DGEList(counts = amrCountsMatrix.min200count.BRR, group = amr.group)
amrExpStd=calcNormFactors(amrExpStd)
amrExpStd=estimateDisp(amrExpStd, amr.design, robust = T)
amr.fit = glmFit(amrExpStd, amr.design)
# AG-biased genes
amr.AG.Contrasts=makeContrasts(AG.v.CR=amr_AG-amr_CR, AG.v.EB=amr_AG-amr_EB, AG.v.TS=amr_AG-amr_TS, levels = amr.design)
amr.lrt.AG.v.rest = glmLRT(amr.fit, contrast = amr.AG.Contrasts)
amr.lrt.AG.v.rest.tTags = topTags(amr.lrt.AG.v.rest, n = NULL)
amr.lrt.AG.v.rest.tTags.table = amr.lrt.AG.v.rest.tTags$table
amr.AG.list=rownames(subset(amr.lrt.AG.v.rest.tTags.table, logFC.AG.v.CR > 2 & logFC.AG.v.EB > 2 & logFC.AG.v.TS > 2 & FDR<0.001))
# extract candidate SFPs from AG list
amr.SFP.list=unique(subset(amrTrinotate, gene_id %in% amr.AG.list & SignalP != "NA")$gene_id)
# EB-biased genes
amr.EB.Contrasts=makeContrasts(EB.v.AG=amr_EB-amr_AG, EB.v.CR=amr_EB-amr_CR, EB.v.TS=amr_EB-amr_TS, levels = amr.design)
amr.lrt.EB.v.rest = glmLRT(amr.fit, contrast = amr.EB.Contrasts)
amr.lrt.EB.v.rest.tTags = topTags(amr.lrt.EB.v.rest, n = NULL)
amr.lrt.EB.v.rest.tTags.table = amr.lrt.EB.v.rest.tTags$table
amr.EB.list=rownames(subset(amr.lrt.EB.v.rest.tTags.table, logFC.EB.v.AG > 2 & logFC.EB.v.CR > 2 & logFC.EB.v.TS > 2 & FDR<0.001))
# TS-biased genes
amr.TS.Contrasts=makeContrasts(TS.v.AG=amr_TS-amr_AG, TS.v.CR=amr_TS-amr_CR, TS.v.EB=amr_TS-amr_EB, levels = amr.design)
amr.lrt.TS.v.rest = glmLRT(amr.fit, contrast = amr.TS.Contrasts)
amr.lrt.TS.v.rest.tTags = topTags(amr.lrt.TS.v.rest, n = NULL)
amr.lrt.TS.v.rest.tTags.table = amr.lrt.TS.v.rest.tTags$table
amr.TS.list=rownames(subset(amr.lrt.TS.v.rest.tTags.table, logFC.TS.v.AG > 2 & logFC.TS.v.CR > 2 & logFC.TS.v.EB > 2 & FDR<0.001))

## D. lummei
# Basic GLM set-up
lum.group = factor(c(1,1,1,2,2,3,3,4,4,4))
lum.design = model.matrix(~0+lum.group)
colnames(lum.design)=levels(factor(lumSamples_data$V1))
lumExpStd=DGEList(counts = lumCountsMatrix.min200count.BRR, group = lum.group)
lumExpStd=calcNormFactors(lumExpStd)
lumExpStd=estimateDisp(lumExpStd, lum.design, robust = T)
lum.fit = glmFit(lumExpStd, lum.design)
# AG-biased genes
lum.AG.Contrasts=makeContrasts(AG.v.CR=lum_AG-lum_CR, AG.v.EB=lum_AG-lum_EB, AG.v.TS=lum_AG-lum_TS, levels = lum.design)
lum.lrt.AG.v.rest = glmLRT(lum.fit, contrast = lum.AG.Contrasts)
lum.lrt.AG.v.rest.tTags = topTags(lum.lrt.AG.v.rest, n = NULL)
lum.lrt.AG.v.rest.tTags.table = lum.lrt.AG.v.rest.tTags$table
lum.AG.list=rownames(subset(lum.lrt.AG.v.rest.tTags.table, logFC.AG.v.CR > 2 & logFC.AG.v.EB > 2 & logFC.AG.v.TS > 2 & FDR<0.001))
# extract candidate SFPs from AG list
lum.SFP.list=unique(subset(lumTrinotate, gene_id %in% lum.AG.list & SignalP != "NA")$gene_id)
# EB-biased genes
lum.EB.Contrasts=makeContrasts(EB.v.AG=lum_EB-lum_AG, EB.v.CR=lum_EB-lum_CR, EB.v.TS=lum_EB-lum_TS, levels = lum.design)
lum.lrt.EB.v.rest = glmLRT(lum.fit, contrast = lum.EB.Contrasts)
lum.lrt.EB.v.rest.tTags = topTags(lum.lrt.EB.v.rest, n = NULL)
lum.lrt.EB.v.rest.tTags.table = lum.lrt.EB.v.rest.tTags$table
lum.EB.list=rownames(subset(lum.lrt.EB.v.rest.tTags.table, logFC.EB.v.AG > 2 & logFC.EB.v.CR > 2 & logFC.EB.v.TS > 2 & FDR<0.001))
# TS-biased genes
lum.TS.Contrasts=makeContrasts(TS.v.AG=lum_TS-lum_AG, TS.v.CR=lum_TS-lum_CR, TS.v.EB=lum_TS-lum_EB, levels = lum.design)
lum.lrt.TS.v.rest = glmLRT(lum.fit, contrast = lum.TS.Contrasts)
lum.lrt.TS.v.rest.tTags = topTags(lum.lrt.TS.v.rest, n = NULL)
lum.lrt.TS.v.rest.tTags.table = lum.lrt.TS.v.rest.tTags$table
lum.TS.list=rownames(subset(lum.lrt.TS.v.rest.tTags.table, logFC.TS.v.AG > 2 & logFC.TS.v.CR > 2 & logFC.TS.v.EB > 2 & FDR<0.001))

## D. novamexicana
# Basic GLM set-up
nov.group = factor(c(1,1,2,2,3,3,4,4,4))
nov.design = model.matrix(~0+nov.group)
colnames(nov.design)=levels(factor(novSamples_data$V1))
novExpStd=DGEList(counts = novCountsMatrix.min200count.BRR, group = nov.group)
novExpStd=calcNormFactors(novExpStd)
novExpStd=estimateDisp(novExpStd, nov.design, robust = T)
nov.fit = glmFit(novExpStd, nov.design)
# AG-biased genes
nov.AG.Contrasts=makeContrasts(AG.v.CR=nov_AG-nov_CR, AG.v.EB=nov_AG-nov_EB, AG.v.TS=nov_AG-nov_TS, levels = nov.design)
nov.lrt.AG.v.rest = glmLRT(nov.fit, contrast = nov.AG.Contrasts)
nov.lrt.AG.v.rest.tTags = topTags(nov.lrt.AG.v.rest, n = NULL)
nov.lrt.AG.v.rest.tTags.table = nov.lrt.AG.v.rest.tTags$table
nov.AG.list=rownames(subset(nov.lrt.AG.v.rest.tTags.table, logFC.AG.v.CR > 2 & logFC.AG.v.EB > 2 & logFC.AG.v.TS > 2 & FDR<0.001))
# extract candidate SFPs from AG list
nov.SFP.list=unique(subset(novTrinotate, gene_id %in% nov.AG.list & SignalP != "NA")$gene_id)
# EB-biased genes
nov.EB.Contrasts=makeContrasts(EB.v.AG=nov_EB-nov_AG, EB.v.CR=nov_EB-nov_CR, EB.v.TS=nov_EB-nov_TS, levels = nov.design)
nov.lrt.EB.v.rest = glmLRT(nov.fit, contrast = nov.EB.Contrasts)
nov.lrt.EB.v.rest.tTags = topTags(nov.lrt.EB.v.rest, n = NULL)
nov.lrt.EB.v.rest.tTags.table = nov.lrt.EB.v.rest.tTags$table
nov.EB.list=rownames(subset(nov.lrt.EB.v.rest.tTags.table, logFC.EB.v.AG > 2 & logFC.EB.v.CR > 2 & logFC.EB.v.TS > 2 & FDR<0.001))
# TS-biased genes
nov.TS.Contrasts=makeContrasts(TS.v.AG=nov_TS-nov_AG, TS.v.CR=nov_TS-nov_CR, TS.v.EB=nov_TS-nov_EB, levels = nov.design)
nov.lrt.TS.v.rest = glmLRT(nov.fit, contrast = nov.TS.Contrasts)
nov.lrt.TS.v.rest.tTags = topTags(nov.lrt.TS.v.rest, n = NULL)
nov.lrt.TS.v.rest.tTags.table = nov.lrt.TS.v.rest.tTags$table
nov.TS.list=rownames(subset(nov.lrt.TS.v.rest.tTags.table, logFC.TS.v.AG > 2 & logFC.TS.v.CR > 2 & logFC.TS.v.EB > 2 & FDR<0.001))

## D. virilis
# Basic GLM set-up
vir.group = factor(c(1,1,1,2,2,3,4,4,4))
vir.design = model.matrix(~0+vir.group)
colnames(vir.design)=levels(factor(virSamples_data$V1))
virExpStd=DGEList(counts = virCountsMatrix.min200count.BRR, group = vir.group)
virExpStd=calcNormFactors(virExpStd)
virExpStd=estimateDisp(virExpStd, vir.design, robust = T)
vir.fit = glmFit(virExpStd, vir.design)
# AG-biased genes
vir.AG.Contrasts=makeContrasts(AG.v.CR=vir_AG-vir_CR, AG.v.EB=vir_AG-vir_EB, AG.v.TS=vir_AG-vir_TS, levels = vir.design)
vir.lrt.AG.v.rest = glmLRT(vir.fit, contrast = vir.AG.Contrasts)
vir.lrt.AG.v.rest.tTags = topTags(vir.lrt.AG.v.rest, n = NULL)
vir.lrt.AG.v.rest.tTags.table = vir.lrt.AG.v.rest.tTags$table
vir.AG.list=rownames(subset(vir.lrt.AG.v.rest.tTags.table, logFC.AG.v.CR > 2 & logFC.AG.v.EB > 2 & logFC.AG.v.TS > 2 & FDR<0.001))
# extract candidate SFPs from AG list
vir.SFP.list=unique(subset(virTrinotate, gene_id %in% vir.AG.list & SignalP != "NA")$gene_id)
# EB-biased genes
vir.EB.Contrasts=makeContrasts(EB.v.AG=vir_EB-vir_AG, EB.v.CR=vir_EB-vir_CR, EB.v.TS=vir_EB-vir_TS, levels = vir.design)
vir.lrt.EB.v.rest = glmLRT(vir.fit, contrast = vir.EB.Contrasts)
vir.lrt.EB.v.rest.tTags = topTags(vir.lrt.EB.v.rest, n = NULL)
vir.lrt.EB.v.rest.tTags.table = vir.lrt.EB.v.rest.tTags$table
vir.EB.list=rownames(subset(vir.lrt.EB.v.rest.tTags.table, logFC.EB.v.AG > 2 & logFC.EB.v.CR > 2 & logFC.EB.v.TS > 2 & FDR<0.001))
# TS-biased genes
vir.TS.Contrasts=makeContrasts(TS.v.AG=vir_TS-vir_AG, TS.v.CR=vir_TS-vir_CR, TS.v.EB=vir_TS-vir_EB, levels = vir.design)
vir.lrt.TS.v.rest = glmLRT(vir.fit, contrast = vir.TS.Contrasts)
vir.lrt.TS.v.rest.tTags = topTags(vir.lrt.TS.v.rest, n = NULL)
vir.lrt.TS.v.rest.tTags.table = vir.lrt.TS.v.rest.tTags$table
vir.TS.list=rownames(subset(vir.lrt.TS.v.rest.tTags.table, logFC.TS.v.AG > 2 & logFC.TS.v.CR > 2 & logFC.TS.v.EB > 2 & FDR<0.001))

#############################################################################
#############################################################################
###############  edgeR DE analysis I: group dvir1.06 ########################

grp.group = factor(c(1,1,2,2,3,4,4,5,5,5,6,6,7,7,8,8,8,9,9,10,10,11,11,12,12,12,13,13,13,14,14,15,16,16,16))
grp.design = model.matrix(~0+grp.group)
colnames(grp.design)=levels(factor(grpSamples_data$V1))
grpExpStd=DGEList(counts = grpCountsMatrix.min200count.BRR, group = grp.group)
grpExpStd=calcNormFactors(grpExpStd)
grpExpStd=estimateDisp(grpExpStd, grp.design, robust = T)
grp.fit = glmFit(grpExpStd, grp.design)

### 1. Identify tissue-biased genes within species (>4-fold, < 0.001 FDR)

## D. americana
# AG-biased genes
amr.dvir1.06.AG.Contrasts=makeContrasts(AG.v.CR=amr_AG-amr_CR, AG.v.EB=amr_AG-amr_EB, AG.v.TS=amr_AG-amr_TS, levels = grp.design)
amr.dvir1.06.lrt.AG.v.rest = glmLRT(grp.fit, contrast = amr.dvir1.06.AG.Contrasts)
amr.dvir1.06.lrt.AG.v.rest.tTags = topTags(amr.dvir1.06.lrt.AG.v.rest, n = NULL)
amr.dvir1.06.lrt.AG.v.rest.tTags.table = amr.dvir1.06.lrt.AG.v.rest.tTags$table
amr.dvir1.06.AG.list=rownames(subset(amr.dvir1.06.lrt.AG.v.rest.tTags.table, logFC.AG.v.CR > 2 & logFC.AG.v.EB > 2 & logFC.AG.v.TS > 2 & FDR<0.001))
# extract candidate SFPs from AG list
amr.dvir1.06.SFP.list=unique(subset(Annots, FBgn_ID %in% amr.dvir1.06.AG.list & SignalP == "YES")$FBgn_ID)
# EB-biased genes
amr.dvir1.06.EB.Contrasts=makeContrasts(EB.v.CR=amr_EB-amr_CR, EB.v.AG=amr_EB-amr_AG, EB.v.TS=amr_EB-amr_TS, levels = grp.design)
amr.dvir1.06.lrt.EB.v.rest = glmLRT(grp.fit, contrast = amr.dvir1.06.EB.Contrasts)
amr.dvir1.06.lrt.EB.v.rest.tTags = topTags(amr.dvir1.06.lrt.EB.v.rest, n = NULL)
amr.dvir1.06.lrt.EB.v.rest.tTags.table = amr.dvir1.06.lrt.EB.v.rest.tTags$table
amr.dvir1.06.EB.list=rownames(subset(amr.dvir1.06.lrt.EB.v.rest.tTags.table, logFC.EB.v.CR > 2 & logFC.EB.v.AG > 2 & logFC.EB.v.TS > 2 & FDR<0.001))
# TS-biased genes
amr.dvir1.06.TS.Contrasts=makeContrasts(TS.v.CR=amr_TS-amr_CR, TS.v.EB=amr_TS-amr_EB, TS.v.AG=amr_TS-amr_AG, levels = grp.design)
amr.dvir1.06.lrt.TS.v.rest = glmLRT(grp.fit, contrast = amr.dvir1.06.TS.Contrasts)
amr.dvir1.06.lrt.TS.v.rest.tTags = topTags(amr.dvir1.06.lrt.TS.v.rest, n = NULL)
amr.dvir1.06.lrt.TS.v.rest.tTags.table = amr.dvir1.06.lrt.TS.v.rest.tTags$table
amr.dvir1.06.TS.list=rownames(subset(amr.dvir1.06.lrt.TS.v.rest.tTags.table, logFC.TS.v.CR > 2 & logFC.TS.v.EB > 2 & logFC.TS.v.AG > 2 & FDR<0.001))

## D. lummei
# AG-biased genes
lum.dvir1.06.AG.Contrasts=makeContrasts(AG.v.CR=lum_AG-lum_CR, AG.v.EB=lum_AG-lum_EB, AG.v.TS=lum_AG-lum_TS, levels = grp.design)
lum.dvir1.06.lrt.AG.v.rest = glmLRT(grp.fit, contrast = lum.dvir1.06.AG.Contrasts)
lum.dvir1.06.lrt.AG.v.rest.tTags = topTags(lum.dvir1.06.lrt.AG.v.rest, n = NULL)
lum.dvir1.06.lrt.AG.v.rest.tTags.table = lum.dvir1.06.lrt.AG.v.rest.tTags$table
lum.dvir1.06.AG.list=rownames(subset(lum.dvir1.06.lrt.AG.v.rest.tTags.table, logFC.AG.v.CR > 2 & logFC.AG.v.EB > 2 & logFC.AG.v.TS > 2 & FDR<0.001))
# extract candidate SFPs from AG list
lum.dvir1.06.SFP.list=unique(subset(Annots, FBgn_ID %in% lum.dvir1.06.AG.list & SignalP == "YES")$FBgn_ID)
# EB-biased genes
lum.dvir1.06.EB.Contrasts=makeContrasts(EB.v.CR=lum_EB-lum_CR, EB.v.AG=lum_EB-lum_AG, EB.v.TS=lum_EB-lum_TS, levels = grp.design)
lum.dvir1.06.lrt.EB.v.rest = glmLRT(grp.fit, contrast = lum.dvir1.06.EB.Contrasts)
lum.dvir1.06.lrt.EB.v.rest.tTags = topTags(lum.dvir1.06.lrt.EB.v.rest, n = NULL)
lum.dvir1.06.lrt.EB.v.rest.tTags.table = lum.dvir1.06.lrt.EB.v.rest.tTags$table
lum.dvir1.06.EB.list=rownames(subset(lum.dvir1.06.lrt.EB.v.rest.tTags.table, logFC.EB.v.CR > 2 & logFC.EB.v.AG > 2 & logFC.EB.v.TS > 2 & FDR<0.001))
# TS-biased genes
lum.dvir1.06.TS.Contrasts=makeContrasts(TS.v.CR=lum_TS-lum_CR, TS.v.EB=lum_TS-lum_EB, TS.v.AG=lum_TS-lum_AG, levels = grp.design)
lum.dvir1.06.lrt.TS.v.rest = glmLRT(grp.fit, contrast = lum.dvir1.06.TS.Contrasts)
lum.dvir1.06.lrt.TS.v.rest.tTags = topTags(lum.dvir1.06.lrt.TS.v.rest, n = NULL)
lum.dvir1.06.lrt.TS.v.rest.tTags.table = lum.dvir1.06.lrt.TS.v.rest.tTags$table
lum.dvir1.06.TS.list=rownames(subset(lum.dvir1.06.lrt.TS.v.rest.tTags.table, logFC.TS.v.CR > 2 & logFC.TS.v.EB > 2 & logFC.TS.v.AG > 2 & FDR<0.001))

## D. novamexicana
# AG-biased genes
nov.dvir1.06.AG.Contrasts=makeContrasts(AG.v.CR=nov_AG-nov_CR, AG.v.EB=nov_AG-nov_EB, AG.v.TS=nov_AG-nov_TS, levels = grp.design)
nov.dvir1.06.lrt.AG.v.rest = glmLRT(grp.fit, contrast = nov.dvir1.06.AG.Contrasts)
nov.dvir1.06.lrt.AG.v.rest.tTags = topTags(nov.dvir1.06.lrt.AG.v.rest, n = NULL)
nov.dvir1.06.lrt.AG.v.rest.tTags.table = nov.dvir1.06.lrt.AG.v.rest.tTags$table
nov.dvir1.06.AG.list=rownames(subset(nov.dvir1.06.lrt.AG.v.rest.tTags.table, logFC.AG.v.CR > 2 & logFC.AG.v.EB > 2 & logFC.AG.v.TS > 2 & FDR<0.001))
# extract candidate SFPs from AG list
nov.dvir1.06.SFP.list=unique(subset(Annots, FBgn_ID %in% nov.dvir1.06.AG.list & SignalP == "YES")$FBgn_ID)
# EB-biased genes
nov.dvir1.06.EB.Contrasts=makeContrasts(EB.v.CR=nov_EB-nov_CR, EB.v.AG=nov_EB-nov_AG, EB.v.TS=nov_EB-nov_TS, levels = grp.design)
nov.dvir1.06.lrt.EB.v.rest = glmLRT(grp.fit, contrast = nov.dvir1.06.EB.Contrasts)
nov.dvir1.06.lrt.EB.v.rest.tTags = topTags(nov.dvir1.06.lrt.EB.v.rest, n = NULL)
nov.dvir1.06.lrt.EB.v.rest.tTags.table = nov.dvir1.06.lrt.EB.v.rest.tTags$table
nov.dvir1.06.EB.list=rownames(subset(nov.dvir1.06.lrt.EB.v.rest.tTags.table, logFC.EB.v.CR > 2 & logFC.EB.v.AG > 2 & logFC.EB.v.TS > 2 & FDR<0.001))
# TS-biased genes
nov.dvir1.06.TS.Contrasts=makeContrasts(TS.v.CR=nov_TS-nov_CR, TS.v.EB=nov_TS-nov_EB, TS.v.AG=nov_TS-nov_AG, levels = grp.design)
nov.dvir1.06.lrt.TS.v.rest = glmLRT(grp.fit, contrast = nov.dvir1.06.TS.Contrasts)
nov.dvir1.06.lrt.TS.v.rest.tTags = topTags(nov.dvir1.06.lrt.TS.v.rest, n = NULL)
nov.dvir1.06.lrt.TS.v.rest.tTags.table = nov.dvir1.06.lrt.TS.v.rest.tTags$table
nov.dvir1.06.TS.list=rownames(subset(nov.dvir1.06.lrt.TS.v.rest.tTags.table, logFC.TS.v.CR > 2 & logFC.TS.v.EB > 2 & logFC.TS.v.AG > 2 & FDR<0.001))

## D. virilis
# AG-biased genes
vir.dvir1.06.AG.Contrasts=makeContrasts(AG.v.CR=vir_AG-vir_CR, AG.v.EB=vir_AG-vir_EB, AG.v.TS=vir_AG-vir_TS, levels = grp.design)
vir.dvir1.06.lrt.AG.v.rest = glmLRT(grp.fit, contrast = vir.dvir1.06.AG.Contrasts)
vir.dvir1.06.lrt.AG.v.rest.tTags = topTags(vir.dvir1.06.lrt.AG.v.rest, n = NULL)
vir.dvir1.06.lrt.AG.v.rest.tTags.table = vir.dvir1.06.lrt.AG.v.rest.tTags$table
vir.dvir1.06.AG.list=rownames(subset(vir.dvir1.06.lrt.AG.v.rest.tTags.table, logFC.AG.v.CR > 2 & logFC.AG.v.EB > 2 & logFC.AG.v.TS > 2 & FDR<0.001))
# extract candidate SFPs from AG list
vir.dvir1.06.SFP.list=unique(subset(Annots, FBgn_ID %in% vir.dvir1.06.AG.list & SignalP == "YES")$FBgn_ID)
# EB-biased genes
vir.dvir1.06.EB.Contrasts=makeContrasts(EB.v.CR=vir_EB-vir_CR, EB.v.AG=vir_EB-vir_AG, EB.v.TS=vir_EB-vir_TS, levels = grp.design)
vir.dvir1.06.lrt.EB.v.rest = glmLRT(grp.fit, contrast = vir.dvir1.06.EB.Contrasts)
vir.dvir1.06.lrt.EB.v.rest.tTags = topTags(vir.dvir1.06.lrt.EB.v.rest, n = NULL)
vir.dvir1.06.lrt.EB.v.rest.tTags.table = vir.dvir1.06.lrt.EB.v.rest.tTags$table
vir.dvir1.06.EB.list=rownames(subset(vir.dvir1.06.lrt.EB.v.rest.tTags.table, logFC.EB.v.CR > 2 & logFC.EB.v.AG > 2 & logFC.EB.v.TS > 2 & FDR<0.001))
# TS-biased genes
vir.dvir1.06.TS.Contrasts=makeContrasts(TS.v.CR=vir_TS-vir_CR, TS.v.EB=vir_TS-vir_EB, TS.v.AG=vir_TS-vir_AG, levels = grp.design)
vir.dvir1.06.lrt.TS.v.rest = glmLRT(grp.fit, contrast = vir.dvir1.06.TS.Contrasts)
vir.dvir1.06.lrt.TS.v.rest.tTags = topTags(vir.dvir1.06.lrt.TS.v.rest, n = NULL)
vir.dvir1.06.lrt.TS.v.rest.tTags.table = vir.dvir1.06.lrt.TS.v.rest.tTags$table
vir.dvir1.06.TS.list=rownames(subset(vir.dvir1.06.lrt.TS.v.rest.tTags.table, logFC.TS.v.CR > 2 & logFC.TS.v.EB > 2 & logFC.TS.v.AG > 2 & FDR<0.001))

# Create lists containing tissue-biased candidates by species
AG_candidates = list(D.ame = amr.dvir1.06.AG.list, 
                      D.lum = lum.dvir1.06.AG.list, 
                      D.nov = nov.dvir1.06.AG.list, 
                      D.vir = vir.dvir1.06.AG.list)
SFP_candidates = list(D.ame = amr.dvir1.06.SFP.list, 
                       D.lum = lum.dvir1.06.SFP.list, 
                       D.nov = nov.dvir1.06.SFP.list, 
                       D.vir = vir.dvir1.06.SFP.list)
EB_candidates = list(D.ame = amr.dvir1.06.EB.list, 
                      D.lum = lum.dvir1.06.EB.list, 
                      D.nov = nov.dvir1.06.EB.list, 
                      D.vir = vir.dvir1.06.EB.list)
TS_candidates = list(D.ame = amr.dvir1.06.TS.list, 
                      D.lum = lum.dvir1.06.TS.list, 
                      D.nov = nov.dvir1.06.TS.list, 
                      D.vir = vir.dvir1.06.TS.list)

# Rearrange into lists of lists, and partition by species
AG_combs = unlist(lapply(1:length(AG_candidates), function(j) combn(names(AG_candidates), j, simplify = FALSE)), recursive = FALSE)
SFP_combs = unlist(lapply(1:length(SFP_candidates), function(j) combn(names(SFP_candidates), j, simplify = FALSE)), recursive = FALSE)
EB_combs = unlist(lapply(1:length(EB_candidates), function(j) combn(names(EB_candidates), j, simplify = FALSE)), recursive = FALSE)
TS_combs = unlist(lapply(1:length(TS_candidates), function(j) combn(names(TS_candidates), j, simplify = FALSE)), recursive = FALSE)
names(AG_combs) = sapply(AG_combs, function(i) paste0(i, collapse = ","))
names(SFP_combs) = sapply(SFP_combs, function(i) paste0(i, collapse = ","))
names(EB_combs) = sapply(EB_combs, function(i) paste0(i, collapse = ","))
names(TS_combs) = sapply(TS_combs, function(i) paste0(i, collapse = ","))
AG_elements = lapply(AG_combs, function(i) Setdiff(AG_candidates[i], AG_candidates[setdiff(names(AG_candidates), i)]))
SFP_elements = lapply(SFP_combs, function(i) Setdiff(SFP_candidates[i], SFP_candidates[setdiff(names(SFP_candidates), i)]))
EB_elements = lapply(EB_combs, function(i) Setdiff(EB_candidates[i], EB_candidates[setdiff(names(EB_candidates), i)]))
TS_elements = lapply(TS_combs, function(i) Setdiff(TS_candidates[i], TS_candidates[setdiff(names(TS_candidates), i)]))

# example to show the partitioning of genes within each tissue element
sapply(AG_elements, length)

### Draw a VennDiagram of each element
AG_candidates_Vdiag=venn.diagram(AG_candidates, NULL, fill=c("#670066", "#86bd38", "#29b5ff", "#717e36"), alpha=c(0.5,0.5,0.5,0.5), cex = 0.8, cat.fontface= 4, cat.cex = 1, resolution = 1000, main = "acc. glands", main.cex = 1.6)
SFP_candidates_Vdiag=venn.diagram(SFP_candidates, NULL, fill=c("#670066", "#86bd38", "#29b5ff", "#717e36"), alpha=c(0.5,0.5,0.5,0.5), cex = 0.8, cat.fontface= 4, cat.cex = 1, resolution = 1000, main = "SFPs", main.cex = 1.6)
EB_candidates_Vdiag=venn.diagram(EB_candidates, NULL, fill=c("#670066", "#86bd38", "#29b5ff", "#717e36"), alpha=c(0.5,0.5,0.5,0.5), cex = 0.8, cat.fontface= 4, cat.cex = 1, resolution = 1000, main = "ejac. bulb", main.cex = 1.6)
TS_candidates_Vdiag=venn.diagram(TS_candidates, NULL, fill=c("#670066", "#86bd38", "#29b5ff", "#717e36"), alpha=c(0.5,0.5,0.5,0.5), cex = 0.8, cat.fontface= 4, cat.cex = 1, resolution = 1000, main = "testes", main.cex = 1.6)

# PLOT for 
pdf("ManuscripPlots/Figure.S8.VennDiagram.pdf", width = 6.9, height = 4.9)
grid.arrange(gTree(children=AG_candidates_Vdiag), gTree(children =SFP_candidates_Vdiag), gTree(children=EB_candidates_Vdiag), gTree(children=TS_candidates_Vdiag), ncol=2)
dev.off()


############################################################################
### 2. Asses differential expression of tissue-biased genes between species.

## Accessory Glands 
crossSpecies.AG.Contrasts=makeContrasts(amr.v.lum=amr_AG-lum_AG, amr.v.nov=amr_AG-nov_AG, amr.v.vir=amr_AG-vir_AG, lum.v.nov=lum_AG-nov_AG, lum.v.vir=lum_AG-vir_AG, nov.v.vir=nov_AG-vir_AG, levels = grp.design)

amr.v.lum.lrt.crossSpecies.AG.Contrasts = glmLRT(grp.fit, contrast = crossSpecies.AG.Contrasts[,"amr.v.lum"])
amr.v.nov.lrt.crossSpecies.AG.Contrasts = glmLRT(grp.fit, contrast = crossSpecies.AG.Contrasts[,"amr.v.nov"])
amr.v.vir.lrt.crossSpecies.AG.Contrasts = glmLRT(grp.fit, contrast = crossSpecies.AG.Contrasts[,"amr.v.vir"])
lum.v.nov.lrt.crossSpecies.AG.Contrasts = glmLRT(grp.fit, contrast = crossSpecies.AG.Contrasts[,"lum.v.nov"])
lum.v.vir.lrt.crossSpecies.AG.Contrasts = glmLRT(grp.fit, contrast = crossSpecies.AG.Contrasts[,"lum.v.vir"])
nov.v.vir.lrt.crossSpecies.AG.Contrasts = glmLRT(grp.fit, contrast = crossSpecies.AG.Contrasts[,"nov.v.vir"])

amr.v.lum.lrt.crossSpecies.AG.Contrasts.tTags = topTags(amr.v.lum.lrt.crossSpecies.AG.Contrasts, n = NULL)
amr.v.nov.lrt.crossSpecies.AG.Contrasts.tTags = topTags(amr.v.nov.lrt.crossSpecies.AG.Contrasts, n = NULL)
amr.v.vir.lrt.crossSpecies.AG.Contrasts.tTags = topTags(amr.v.vir.lrt.crossSpecies.AG.Contrasts, n = NULL)
lum.v.nov.lrt.crossSpecies.AG.Contrasts.tTags = topTags(lum.v.nov.lrt.crossSpecies.AG.Contrasts, n = NULL)
lum.v.vir.lrt.crossSpecies.AG.Contrasts.tTags = topTags(lum.v.vir.lrt.crossSpecies.AG.Contrasts, n = NULL)
nov.v.vir.lrt.crossSpecies.AG.Contrasts.tTags = topTags(nov.v.vir.lrt.crossSpecies.AG.Contrasts, n = NULL)

amr.v.lum.lrt.crossSpecies.AG.Contrasts.tTags.table = amr.v.lum.lrt.crossSpecies.AG.Contrasts.tTags$table
amr.v.nov.lrt.crossSpecies.AG.Contrasts.tTags.table = amr.v.nov.lrt.crossSpecies.AG.Contrasts.tTags$table
amr.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table = amr.v.vir.lrt.crossSpecies.AG.Contrasts.tTags$table
lum.v.nov.lrt.crossSpecies.AG.Contrasts.tTags.table = lum.v.nov.lrt.crossSpecies.AG.Contrasts.tTags$table
lum.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table = lum.v.vir.lrt.crossSpecies.AG.Contrasts.tTags$table
nov.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table = nov.v.vir.lrt.crossSpecies.AG.Contrasts.tTags$table

# combine specifity data with tTag tables:
AG.dvir1.06_Specificity_table$FBgn_ID = rownames(AG.dvir1.06_Specificity_table )
rownames(AG.dvir1.06_Specificity_table ) = NULL

# subset AG samples only
dvir1.06_Specificity_table.amr.v.lum.samples = subset(AG.dvir1.06_Specificity_table, select=c("FBgn_ID", "amr_AG", "lum_AG"))
amr.v.lum.lrt.crossSpecies.AG.Contrasts.tTags.table$FBgn_ID = rownames(amr.v.lum.lrt.crossSpecies.AG.Contrasts.tTags.table)
rownames(amr.v.lum.lrt.crossSpecies.AG.Contrasts.tTags.table) = NULL
amr.v.lum.lrt.crossSpecies.AG.Contrasts.tTags.table = merge(amr.v.lum.lrt.crossSpecies.AG.Contrasts.tTags.table, dvir1.06_Specificity_table.amr.v.lum.samples)
amr.v.lum.lrt.crossSpecies.AG.Contrasts.tTags.table$comparison = "D.lum vs D.amr"
amr.v.lum.lrt.crossSpecies.AG.Contrasts.tTags.table$tissue = "Accessory glands"
colnames(amr.v.lum.lrt.crossSpecies.AG.Contrasts.tTags.table)[7:8] = c("S.right", "S.left")

dvir1.06_Specificity_table.amr.v.nov.samples = subset(AG.dvir1.06_Specificity_table, select=c("FBgn_ID", "amr_AG", "nov_AG"))
amr.v.nov.lrt.crossSpecies.AG.Contrasts.tTags.table$FBgn_ID = rownames(amr.v.nov.lrt.crossSpecies.AG.Contrasts.tTags.table)
rownames(amr.v.nov.lrt.crossSpecies.AG.Contrasts.tTags.table) = NULL
amr.v.nov.lrt.crossSpecies.AG.Contrasts.tTags.table = merge(amr.v.nov.lrt.crossSpecies.AG.Contrasts.tTags.table, dvir1.06_Specificity_table.amr.v.nov.samples)
amr.v.nov.lrt.crossSpecies.AG.Contrasts.tTags.table$comparison = "D.nov vs D.amr"
amr.v.nov.lrt.crossSpecies.AG.Contrasts.tTags.table$tissue = "Accessory glands"
colnames(amr.v.nov.lrt.crossSpecies.AG.Contrasts.tTags.table)[7:8] = c("S.right", "S.left")

dvir1.06_Specificity_table.amr.v.vir.samples = subset(AG.dvir1.06_Specificity_table, select=c("FBgn_ID", "amr_AG", "vir_AG"))
amr.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table$FBgn_ID = rownames(amr.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table)
rownames(amr.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table) = NULL
amr.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table = merge(amr.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table, dvir1.06_Specificity_table.amr.v.vir.samples)
amr.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table$comparison = "D.vir vs D.amr"
amr.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table$tissue = "Accessory glands"
colnames(amr.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table)[7:8] = c("S.right", "S.left")

dvir1.06_Specificity_table.lum.v.nov.samples = subset(AG.dvir1.06_Specificity_table, select=c("FBgn_ID", "lum_AG", "nov_AG"))
lum.v.nov.lrt.crossSpecies.AG.Contrasts.tTags.table$FBgn_ID = rownames(lum.v.nov.lrt.crossSpecies.AG.Contrasts.tTags.table)
rownames(lum.v.nov.lrt.crossSpecies.AG.Contrasts.tTags.table) = NULL
lum.v.nov.lrt.crossSpecies.AG.Contrasts.tTags.table = merge(lum.v.nov.lrt.crossSpecies.AG.Contrasts.tTags.table, dvir1.06_Specificity_table.lum.v.nov.samples)
lum.v.nov.lrt.crossSpecies.AG.Contrasts.tTags.table$comparison = "D.nov vs. D.lum"
lum.v.nov.lrt.crossSpecies.AG.Contrasts.tTags.table$tissue = "Accessory glands"
colnames(lum.v.nov.lrt.crossSpecies.AG.Contrasts.tTags.table)[7:8] = c("S.right", "S.left")

dvir1.06_Specificity_table.lum.v.vir.samples = subset(AG.dvir1.06_Specificity_table, select=c("FBgn_ID", "lum_AG", "vir_AG"))
lum.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table$FBgn_ID = rownames(lum.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table)
rownames(lum.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table) = NULL
lum.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table = merge(lum.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table, dvir1.06_Specificity_table.lum.v.vir.samples)
lum.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table$comparison = "D.vir vs D.lum"
lum.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table$tissue = "Accessory glands"
colnames(lum.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table)[7:8] = c("S.right", "S.left")

dvir1.06_Specificity_table.nov.v.vir.samples = subset(AG.dvir1.06_Specificity_table, select=c("FBgn_ID", "nov_AG", "vir_AG"))
nov.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table$FBgn_ID = rownames(nov.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table)
rownames(nov.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table) = NULL
nov.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table = merge(nov.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table, dvir1.06_Specificity_table.nov.v.vir.samples)
nov.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table$comparison = "D.vir vs. D.nov"
nov.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table$tissue = "Accessory glands"
colnames(nov.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table)[7:8] = c("S.right", "S.left")

allAG.tmp = rbind(amr.v.lum.lrt.crossSpecies.AG.Contrasts.tTags.table, amr.v.nov.lrt.crossSpecies.AG.Contrasts.tTags.table, amr.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table, lum.v.nov.lrt.crossSpecies.AG.Contrasts.tTags.table, lum.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table, nov.v.vir.lrt.crossSpecies.AG.Contrasts.tTags.table)
allAG.tmp = subset(allAG.tmp, FBgn_ID %in% unlist(AG_candidates))
allAG.tmp = subset(allAG.tmp, FBgn_ID %!in% unlist(SFP_candidates))
# extract SFP genes from above tmp file
allSFP.tmp = subset(allAG.tmp, FBgn_ID %in% unlist(SFP_candidates))
allSFP.tmp$tissue = "SFPs"

## EJaculatory Bulb
crossSpecies.EB.Contrasts=makeContrasts(amr.v.lum=amr_EB-lum_EB, amr.v.nov=amr_EB-nov_EB, amr.v.vir=amr_EB-vir_EB, lum.v.nov=lum_EB-nov_EB, lum.v.vir=lum_EB-vir_EB, nov.v.vir=nov_EB-vir_EB, levels = grp.design)

amr.v.lum.lrt.crossSpecies.EB.Contrasts = glmLRT(grp.fit, contrast = crossSpecies.EB.Contrasts[,"amr.v.lum"])
amr.v.nov.lrt.crossSpecies.EB.Contrasts = glmLRT(grp.fit, contrast = crossSpecies.EB.Contrasts[,"amr.v.nov"])
amr.v.vir.lrt.crossSpecies.EB.Contrasts = glmLRT(grp.fit, contrast = crossSpecies.EB.Contrasts[,"amr.v.vir"])
lum.v.nov.lrt.crossSpecies.EB.Contrasts = glmLRT(grp.fit, contrast = crossSpecies.EB.Contrasts[,"lum.v.nov"])
lum.v.vir.lrt.crossSpecies.EB.Contrasts = glmLRT(grp.fit, contrast = crossSpecies.EB.Contrasts[,"lum.v.vir"])
nov.v.vir.lrt.crossSpecies.EB.Contrasts = glmLRT(grp.fit, contrast = crossSpecies.EB.Contrasts[,"nov.v.vir"])

amr.v.lum.lrt.crossSpecies.EB.Contrasts.tTags = topTags(amr.v.lum.lrt.crossSpecies.EB.Contrasts, n = NULL)
amr.v.nov.lrt.crossSpecies.EB.Contrasts.tTags = topTags(amr.v.nov.lrt.crossSpecies.EB.Contrasts, n = NULL)
amr.v.vir.lrt.crossSpecies.EB.Contrasts.tTags = topTags(amr.v.vir.lrt.crossSpecies.EB.Contrasts, n = NULL)
lum.v.nov.lrt.crossSpecies.EB.Contrasts.tTags = topTags(lum.v.nov.lrt.crossSpecies.EB.Contrasts, n = NULL)
lum.v.vir.lrt.crossSpecies.EB.Contrasts.tTags = topTags(lum.v.vir.lrt.crossSpecies.EB.Contrasts, n = NULL)
nov.v.vir.lrt.crossSpecies.EB.Contrasts.tTags = topTags(nov.v.vir.lrt.crossSpecies.EB.Contrasts, n = NULL)

amr.v.lum.lrt.crossSpecies.EB.Contrasts.tTags.table = amr.v.lum.lrt.crossSpecies.EB.Contrasts.tTags$table
amr.v.nov.lrt.crossSpecies.EB.Contrasts.tTags.table = amr.v.nov.lrt.crossSpecies.EB.Contrasts.tTags$table
amr.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table = amr.v.vir.lrt.crossSpecies.EB.Contrasts.tTags$table
lum.v.nov.lrt.crossSpecies.EB.Contrasts.tTags.table = lum.v.nov.lrt.crossSpecies.EB.Contrasts.tTags$table
lum.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table = lum.v.vir.lrt.crossSpecies.EB.Contrasts.tTags$table
nov.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table = nov.v.vir.lrt.crossSpecies.EB.Contrasts.tTags$table

# combine specifity data with tTag tables:
EB.dvir1.06_Specificity_table$FBgn_ID = rownames(EB.dvir1.06_Specificity_table )
rownames(EB.dvir1.06_Specificity_table ) = NULL

# subset EB samples only
dvir1.06_Specificity_table.amr.v.lum.samples = subset(EB.dvir1.06_Specificity_table, select=c("FBgn_ID", "amr_EB", "lum_EB"))
amr.v.lum.lrt.crossSpecies.EB.Contrasts.tTags.table$FBgn_ID = rownames(amr.v.lum.lrt.crossSpecies.EB.Contrasts.tTags.table)
rownames(amr.v.lum.lrt.crossSpecies.EB.Contrasts.tTags.table) = NULL
amr.v.lum.lrt.crossSpecies.EB.Contrasts.tTags.table = merge(amr.v.lum.lrt.crossSpecies.EB.Contrasts.tTags.table, dvir1.06_Specificity_table.amr.v.lum.samples)
amr.v.lum.lrt.crossSpecies.EB.Contrasts.tTags.table$comparison = "D.lum vs D.amr"
amr.v.lum.lrt.crossSpecies.EB.Contrasts.tTags.table$tissue = "Ejaculatory bulb"
colnames(amr.v.lum.lrt.crossSpecies.EB.Contrasts.tTags.table)[7:8] = c("S.right", "S.left")

dvir1.06_Specificity_table.amr.v.nov.samples = subset(EB.dvir1.06_Specificity_table, select=c("FBgn_ID", "amr_EB", "nov_EB"))
amr.v.nov.lrt.crossSpecies.EB.Contrasts.tTags.table$FBgn_ID = rownames(amr.v.nov.lrt.crossSpecies.EB.Contrasts.tTags.table)
rownames(amr.v.nov.lrt.crossSpecies.EB.Contrasts.tTags.table) = NULL
amr.v.nov.lrt.crossSpecies.EB.Contrasts.tTags.table = merge(amr.v.nov.lrt.crossSpecies.EB.Contrasts.tTags.table, dvir1.06_Specificity_table.amr.v.nov.samples)
amr.v.nov.lrt.crossSpecies.EB.Contrasts.tTags.table$comparison = "D.nov vs D.amr"
amr.v.nov.lrt.crossSpecies.EB.Contrasts.tTags.table$tissue = "Ejaculatory bulb"
colnames(amr.v.nov.lrt.crossSpecies.EB.Contrasts.tTags.table)[7:8] = c("S.right", "S.left")

dvir1.06_Specificity_table.amr.v.vir.samples = subset(EB.dvir1.06_Specificity_table, select=c("FBgn_ID", "amr_EB", "vir_EB"))
amr.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table$FBgn_ID = rownames(amr.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table)
rownames(amr.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table) = NULL
amr.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table = merge(amr.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table, dvir1.06_Specificity_table.amr.v.vir.samples)
amr.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table$comparison = "D.vir vs D.amr"
amr.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table$tissue = "Ejaculatory bulb"
colnames(amr.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table)[7:8] = c("S.right", "S.left")

dvir1.06_Specificity_table.lum.v.nov.samples = subset(EB.dvir1.06_Specificity_table, select=c("FBgn_ID", "lum_EB", "nov_EB"))
lum.v.nov.lrt.crossSpecies.EB.Contrasts.tTags.table$FBgn_ID = rownames(lum.v.nov.lrt.crossSpecies.EB.Contrasts.tTags.table)
rownames(lum.v.nov.lrt.crossSpecies.EB.Contrasts.tTags.table) = NULL
lum.v.nov.lrt.crossSpecies.EB.Contrasts.tTags.table = merge(lum.v.nov.lrt.crossSpecies.EB.Contrasts.tTags.table, dvir1.06_Specificity_table.lum.v.nov.samples)
lum.v.nov.lrt.crossSpecies.EB.Contrasts.tTags.table$comparison = "D.nov vs. D.lum"
lum.v.nov.lrt.crossSpecies.EB.Contrasts.tTags.table$tissue = "Ejaculatory bulb"
colnames(lum.v.nov.lrt.crossSpecies.EB.Contrasts.tTags.table)[7:8] = c("S.right", "S.left")

dvir1.06_Specificity_table.lum.v.vir.samples = subset(EB.dvir1.06_Specificity_table, select=c("FBgn_ID", "lum_EB", "vir_EB"))
lum.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table$FBgn_ID = rownames(lum.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table)
rownames(lum.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table) = NULL
lum.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table = merge(lum.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table, dvir1.06_Specificity_table.lum.v.vir.samples)
lum.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table$comparison = "D.vir vs D.lum"
lum.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table$tissue = "Ejaculatory bulb"
colnames(lum.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table)[7:8] = c("S.right", "S.left")

dvir1.06_Specificity_table.nov.v.vir.samples = subset(EB.dvir1.06_Specificity_table, select=c("FBgn_ID", "nov_EB", "vir_EB"))
nov.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table$FBgn_ID = rownames(nov.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table)
rownames(nov.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table) = NULL
nov.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table = merge(nov.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table, dvir1.06_Specificity_table.nov.v.vir.samples)
nov.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table$comparison = "D.vir vs. D.nov"
nov.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table$tissue = "Ejaculatory bulb"
colnames(nov.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table)[7:8] = c("S.right", "S.left")

allEB.tmp = rbind(amr.v.lum.lrt.crossSpecies.EB.Contrasts.tTags.table, amr.v.nov.lrt.crossSpecies.EB.Contrasts.tTags.table, amr.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table, lum.v.nov.lrt.crossSpecies.EB.Contrasts.tTags.table, lum.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table, nov.v.vir.lrt.crossSpecies.EB.Contrasts.tTags.table)
allEB.tmp = subset(allEB.tmp, FBgn_ID %in% unlist(EB_candidates))



## Tetes 
crossSpecies.TS.Contrasts=makeContrasts(amr.v.lum=amr_TS-lum_TS, amr.v.nov=amr_TS-nov_TS, amr.v.vir=amr_TS-vir_TS, lum.v.nov=lum_TS-nov_TS, lum.v.vir=lum_TS-vir_TS, nov.v.vir=nov_TS-vir_TS, levels = grp.design)

amr.v.lum.lrt.crossSpecies.TS.Contrasts = glmLRT(grp.fit, contrast = crossSpecies.TS.Contrasts[,"amr.v.lum"])
amr.v.nov.lrt.crossSpecies.TS.Contrasts = glmLRT(grp.fit, contrast = crossSpecies.TS.Contrasts[,"amr.v.nov"])
amr.v.vir.lrt.crossSpecies.TS.Contrasts = glmLRT(grp.fit, contrast = crossSpecies.TS.Contrasts[,"amr.v.vir"])
lum.v.nov.lrt.crossSpecies.TS.Contrasts = glmLRT(grp.fit, contrast = crossSpecies.TS.Contrasts[,"lum.v.nov"])
lum.v.vir.lrt.crossSpecies.TS.Contrasts = glmLRT(grp.fit, contrast = crossSpecies.TS.Contrasts[,"lum.v.vir"])
nov.v.vir.lrt.crossSpecies.TS.Contrasts = glmLRT(grp.fit, contrast = crossSpecies.TS.Contrasts[,"nov.v.vir"])

amr.v.lum.lrt.crossSpecies.TS.Contrasts.tTags = topTags(amr.v.lum.lrt.crossSpecies.TS.Contrasts, n = NULL)
amr.v.nov.lrt.crossSpecies.TS.Contrasts.tTags = topTags(amr.v.nov.lrt.crossSpecies.TS.Contrasts, n = NULL)
amr.v.vir.lrt.crossSpecies.TS.Contrasts.tTags = topTags(amr.v.vir.lrt.crossSpecies.TS.Contrasts, n = NULL)
lum.v.nov.lrt.crossSpecies.TS.Contrasts.tTags = topTags(lum.v.nov.lrt.crossSpecies.TS.Contrasts, n = NULL)
lum.v.vir.lrt.crossSpecies.TS.Contrasts.tTags = topTags(lum.v.vir.lrt.crossSpecies.TS.Contrasts, n = NULL)
nov.v.vir.lrt.crossSpecies.TS.Contrasts.tTags = topTags(nov.v.vir.lrt.crossSpecies.TS.Contrasts, n = NULL)

amr.v.lum.lrt.crossSpecies.TS.Contrasts.tTags.table = amr.v.lum.lrt.crossSpecies.TS.Contrasts.tTags$table
amr.v.nov.lrt.crossSpecies.TS.Contrasts.tTags.table = amr.v.nov.lrt.crossSpecies.TS.Contrasts.tTags$table
amr.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table = amr.v.vir.lrt.crossSpecies.TS.Contrasts.tTags$table
lum.v.nov.lrt.crossSpecies.TS.Contrasts.tTags.table = lum.v.nov.lrt.crossSpecies.TS.Contrasts.tTags$table
lum.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table = lum.v.vir.lrt.crossSpecies.TS.Contrasts.tTags$table
nov.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table = nov.v.vir.lrt.crossSpecies.TS.Contrasts.tTags$table

# combine specifity data with tTag tables:
TS.dvir1.06_Specificity_table$FBgn_ID = rownames(TS.dvir1.06_Specificity_table )
rownames(TS.dvir1.06_Specificity_table ) = NULL

# subset TS samples only
dvir1.06_Specificity_table.amr.v.lum.samples = subset(TS.dvir1.06_Specificity_table, select=c("FBgn_ID", "amr_TS", "lum_TS"))
amr.v.lum.lrt.crossSpecies.TS.Contrasts.tTags.table$FBgn_ID = rownames(amr.v.lum.lrt.crossSpecies.TS.Contrasts.tTags.table)
rownames(amr.v.lum.lrt.crossSpecies.TS.Contrasts.tTags.table) = NULL
amr.v.lum.lrt.crossSpecies.TS.Contrasts.tTags.table = merge(amr.v.lum.lrt.crossSpecies.TS.Contrasts.tTags.table, dvir1.06_Specificity_table.amr.v.lum.samples)
amr.v.lum.lrt.crossSpecies.TS.Contrasts.tTags.table$comparison = "D.lum vs D.amr"
amr.v.lum.lrt.crossSpecies.TS.Contrasts.tTags.table$tissue = "Testes"
colnames(amr.v.lum.lrt.crossSpecies.TS.Contrasts.tTags.table)[7:8] = c("S.right", "S.left")

dvir1.06_Specificity_table.amr.v.nov.samples = subset(TS.dvir1.06_Specificity_table, select=c("FBgn_ID", "amr_TS", "nov_TS"))
amr.v.nov.lrt.crossSpecies.TS.Contrasts.tTags.table$FBgn_ID = rownames(amr.v.nov.lrt.crossSpecies.TS.Contrasts.tTags.table)
rownames(amr.v.nov.lrt.crossSpecies.TS.Contrasts.tTags.table) = NULL
amr.v.nov.lrt.crossSpecies.TS.Contrasts.tTags.table = merge(amr.v.nov.lrt.crossSpecies.TS.Contrasts.tTags.table, dvir1.06_Specificity_table.amr.v.nov.samples)
amr.v.nov.lrt.crossSpecies.TS.Contrasts.tTags.table$comparison = "D.nov vs D.amr"
amr.v.nov.lrt.crossSpecies.TS.Contrasts.tTags.table$tissue = "Testes"
colnames(amr.v.nov.lrt.crossSpecies.TS.Contrasts.tTags.table)[7:8] = c("S.right", "S.left")

dvir1.06_Specificity_table.amr.v.vir.samples = subset(TS.dvir1.06_Specificity_table, select=c("FBgn_ID", "amr_TS", "vir_TS"))
amr.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table$FBgn_ID = rownames(amr.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table)
rownames(amr.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table) = NULL
amr.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table = merge(amr.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table, dvir1.06_Specificity_table.amr.v.vir.samples)
amr.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table$comparison = "D.vir vs D.amr"
amr.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table$tissue = "Testes"
colnames(amr.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table)[7:8] = c("S.right", "S.left")

dvir1.06_Specificity_table.lum.v.nov.samples = subset(TS.dvir1.06_Specificity_table, select=c("FBgn_ID", "lum_TS", "nov_TS"))
lum.v.nov.lrt.crossSpecies.TS.Contrasts.tTags.table$FBgn_ID = rownames(lum.v.nov.lrt.crossSpecies.TS.Contrasts.tTags.table)
rownames(lum.v.nov.lrt.crossSpecies.TS.Contrasts.tTags.table) = NULL
lum.v.nov.lrt.crossSpecies.TS.Contrasts.tTags.table = merge(lum.v.nov.lrt.crossSpecies.TS.Contrasts.tTags.table, dvir1.06_Specificity_table.lum.v.nov.samples)
lum.v.nov.lrt.crossSpecies.TS.Contrasts.tTags.table$comparison = "D.nov vs. D.lum"
lum.v.nov.lrt.crossSpecies.TS.Contrasts.tTags.table$tissue = "Testes"
colnames(lum.v.nov.lrt.crossSpecies.TS.Contrasts.tTags.table)[7:8] = c("S.right", "S.left")

dvir1.06_Specificity_table.lum.v.vir.samples = subset(TS.dvir1.06_Specificity_table, select=c("FBgn_ID", "lum_TS", "vir_TS"))
lum.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table$FBgn_ID = rownames(lum.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table)
rownames(lum.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table) = NULL
lum.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table = merge(lum.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table, dvir1.06_Specificity_table.lum.v.vir.samples)
lum.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table$comparison = "D.vir vs D.lum"
lum.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table$tissue = "Testes"
colnames(lum.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table)[7:8] = c("S.right", "S.left")

dvir1.06_Specificity_table.nov.v.vir.samples = subset(TS.dvir1.06_Specificity_table, select=c("FBgn_ID", "nov_TS", "vir_TS"))
nov.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table$FBgn_ID = rownames(nov.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table)
rownames(nov.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table) = NULL
nov.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table = merge(nov.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table, dvir1.06_Specificity_table.nov.v.vir.samples)
nov.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table$comparison = "D.vir vs. D.nov"
nov.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table$tissue = "Testes"
colnames(nov.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table)[7:8] = c("S.right", "S.left")


allTS.tmp = rbind(amr.v.lum.lrt.crossSpecies.TS.Contrasts.tTags.table, amr.v.nov.lrt.crossSpecies.TS.Contrasts.tTags.table, amr.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table, lum.v.nov.lrt.crossSpecies.TS.Contrasts.tTags.table, lum.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table, nov.v.vir.lrt.crossSpecies.TS.Contrasts.tTags.table)
allTS.tmp = subset(allTS.tmp, FBgn_ID %in% unlist(TS_candidates))

#### Put it all together
# tissue-biased transcripts:
crossSpecies.ALL.df = rbind(allAG.tmp, allEB.tmp, allTS.tmp)
crossSpecies.ALL.df$Sig = ifelse(crossSpecies.ALL.df$FDR < 0.01, "YES", "NO")
crossSpecies.ALL.df = merge(crossSpecies.ALL.df, gffRecord)

pdf("ManuscripPlots/Figure.S9.tiss-biased.volcanoPlots.pdf", width = 12, height = 6)
ggplot() + 
  geom_point(data = subset(crossSpecies.ALL.df, logFC > 0), aes (logFC, -log10(PValue), colour = Sig, size = S.right), alpha = 0.6) + 
  geom_point(data = subset(crossSpecies.ALL.df, logFC < 0), aes (logFC, -log10(PValue), colour = Sig, size = S.left), alpha = 0.6) + 
  facet_grid(tissue~comparison, scales = "free") + 
  geom_text_repel(data=subset(crossSpecies.ALL.df, S.left > 0.75 & logFC > 2 | S.right > 0.75 & logFC > 2), aes(logFC, -log10(PValue), label = gene_name), size =2.5, force = 10, point.padding = unit(0.5, "lines"), fontface = "bold.italic") + 
  scale_colour_manual(values = c("#88d542", "#0eacc2")) + 
  labs(size="cross-species\nspecificity", colour="FDR < 0.01?") +
  scale_size(range = c(-2, 2)) + 
  ylab("-log10(p-value)") +
  theme(axis.title.x = element_text(face = "bold", size = 12, vjust=0.1), axis.text.x=element_text(face = "bold", size = 8),axis.text.y = element_text(face = "bold", size = 12), axis.title.y = element_text(face = "bold.italic", size = 12, vjust=0.1), strip.text=element_text(face="bold.italic", size = 11))
dev.off()

# SFP candidates
allSFP.tmp$Sig = ifelse(allSFP.tmp$FDR < 0.01, "YES", "NO")
allSFPs = merge(allSFP.tmp, gffRecord)

pdf("ManuscripPlots/Figure.NA.SFPs.volcanoPlots.pdf", width = 12.5, height = 2.73)
ggplot() + 
  geom_point(data = subset(allSFPs, logFC > 0), aes (logFC, -log10(PValue), colour = Sig, size = S.right), alpha = 0.7) + 
  geom_point(data = subset(allSFPs, logFC < 0), aes (logFC, -log10(PValue), colour = Sig, size = S.left), alpha = 0.7) + 
  facet_grid(~comparison, scales = "free") + 
  geom_text_repel(data=subset(allSFPs, S.right > 0.75 & logFC > 2 | S.left > 0.75 & logFC < -2), aes(logFC, -log10(PValue), label = gene_name), size =3.5, force = 30, fontface = "bold.italic", colour = "black") + 
  scale_colour_manual(values = c("#88d542", "#0eacc2")) + 
  labs(size="cross-species\nspecificity", colour="FDR < 0.01?") + 
  scale_size(range = c(-3, 5)) + 
  ylab("-log10(p-value)") +
  theme(axis.title.x = element_text(face = "bold", size = 12, vjust=0.1), axis.text.x=element_text(face = "bold", size = 10), axis.text.y = element_text(face = "bold", size = 10), axis.title.y = element_text(face = "bold.italic", size = 12, vjust=0.1), strip.text=element_text(face="bold.italic", size = 12)) 
dev.off()

#### 3d Plot (not very useful)
#plot_ly(subset(amr.lrt.crossSpecies.TS.Contrasts.tTags.table, rownames(amr.lrt.crossSpecies.TS.Contrasts.tTags.table) %in% amr.crossSpecies.TS.list), x=~logFC.amr.v.lum, y=~logFC.amr.v.nov, z=~logFC.amr.v.vir, color = ~-log10(FDR), type = "scatter3d", colors = )



##### FIND De Novo Transcripts ############################################################
###########################################################################################                                       
##### Here's a method to find Trinity transcript that are not found in dvir1.06 annotation,
##### then checking whether orthologues exist in the other Trinity assemblies

# D. americana SFPs
selectionCols = c("dvir1.06_BlastX_topHit", "dvir1.06_BlastP_topHit", "gene_id", "prot_id")
amr.SFP.dvir1.06.orths = subset(amrTrinotate, gene_id %in% amr.SFP.list)[selectionCols]
amr.SFP.dvir1.06.orths = droplevels(amr.SFP.dvir1.06.orths)
amr.SFP.dvir1.06.orths = amr.SFP.dvir1.06.orths[order(amr.SFP.dvir1.06.orths$dvir1.06_BlastX_topHit), ]
amr.SFP.dvir1.06.orths[is.na(amr.SFP.dvir1.06.orths)] = "NoHit"
amr.SFP.dvir1.06.orths = subset(amr.SFP.dvir1.06.orths, prot_id != "NoHit")
amr.SFP.dvir1.06.orths = aggregate(gene_id~dvir1.06_BlastX_topHit, data = amr.SFP.dvir1.06.orths, toString)
amr.SFP_no_dvir1.06_hits = subset(amr.SFP.dvir1.06.orths, dvir1.06_BlastX_topHit == "NoHit")$gene_id
amr.SFP_no_dvir1.06_hits = unique(strsplit(amr.SFP_no_dvir1.06_hits, ", ")[[1]])
amr.SFP_no_dvir1.06_hits_Trin_hits = as.character(unique(subset(TrinOrths, Gene %in% amr.SFP_no_dvir1.06_hits)$Gene))
setdiff(amr.SFP_no_dvir1.06_hits, amr.SFP_no_dvir1.06_hits_Trin_hits)
# pdf("Plots/genePlots_amr.SFP_no_dvir1.06_hits.pdf", height = 3)
# lapply(amr.SFP_no_dvir1.06_hits, plotGeneT, object=TPMseBRR_amrTrin)
# dev.off()

# D. lummei SFPs
lum.SFP.dvir1.06.orths = subset(lumTrinotate, gene_id %in% lum.SFP.list)[selectionCols]
lum.SFP.dvir1.06.orths = droplevels(lum.SFP.dvir1.06.orths)
lum.SFP.dvir1.06.orths = lum.SFP.dvir1.06.orths[order(lum.SFP.dvir1.06.orths$dvir1.06_BlastX_topHit), ]
lum.SFP.dvir1.06.orths[is.na(lum.SFP.dvir1.06.orths)] = "NoHit"
lum.SFP.dvir1.06.orths = subset(lum.SFP.dvir1.06.orths, prot_id != "NoHit")
lum.SFP.dvir1.06.orths = aggregate(gene_id~dvir1.06_BlastX_topHit, data = lum.SFP.dvir1.06.orths, toString)
lum.SFP_no_dvir1.06_hits = subset(lum.SFP.dvir1.06.orths, dvir1.06_BlastX_topHit == "NoHit")$gene_id
lum.SFP_no_dvir1.06_hits = unique(strsplit(lum.SFP_no_dvir1.06_hits, ", ")[[1]])
lum.SFP_no_dvir1.06_hits_Trin_hits = as.character(unique(subset(TrinOrths, Gene %in% lum.SFP_no_dvir1.06_hits)$Gene))
setdiff(lum.SFP_no_dvir1.06_hits, lum.SFP_no_dvir1.06_hits_Trin_hits)

# D. novamexicana SFPs
nov.SFP.dvir1.06.orths = subset(novTrinotate, gene_id %in% nov.SFP.list)[selectionCols]
nov.SFP.dvir1.06.orths = droplevels(nov.SFP.dvir1.06.orths)
nov.SFP.dvir1.06.orths = nov.SFP.dvir1.06.orths[order(nov.SFP.dvir1.06.orths$dvir1.06_BlastX_topHit), ]
nov.SFP.dvir1.06.orths[is.na(nov.SFP.dvir1.06.orths)] = "NoHit"
nov.SFP.dvir1.06.orths = subset(novSFP.dvir1.06.orths, prot_id != "NoHit")
nov.SFP.dvir1.06.orths = aggregate(gene_id~dvir1.06_BlastX_topHit, data = nov.SFP.dvir1.06.orths, toString)
nov.SFP_no_dvir1.06_hits = subset(nov.SFP.dvir1.06.orths, dvir1.06_BlastX_topHit == "NoHit")$gene_id
nov.SFP_no_dvir1.06_hits = unique(strsplit(nov.SFP_no_dvir1.06_hits, ", ")[[1]])
nov.SFP_no_dvir1.06_hits_Trin_hits = as.character(unique(subset(TrinOrths, Gene %in% nov.SFP_no_dvir1.06_hits)$Gene))
setdiff(nov.SFP_no_dvir1.06_hits, nov.SFP_no_dvir1.06_hits_Trin_hits)

# D. virilis SFPs
vir.SFP.dvir1.06.orths = subset(virTrinotate, gene_id %in% vir.SFP.list)[selectionCols]
vir.SFP.dvir1.06.orths = droplevels(vir.SFP.dvir1.06.orths)
vir.SFP.dvir1.06.orths = vir.SFP.dvir1.06.orths[order(vir.SFP.dvir1.06.orths$dvir1.06_BlastX_topHit), ]
vir.SFP.dvir1.06.orths[is.na(vir.SFP.dvir1.06.orths)] = "NoHit"
vir.SFP.dvir1.06.orths = subset(vir.SFP.dvir1.06.orths, prot_id != "NoHit")
vir.SFP.dvir1.06.orths = aggregate(gene_id~dvir1.06_BlastX_topHit, data = vir.SFP.dvir1.06.orths, toString)
vir.SFP_no_dvir1.06_hits = subset(vir.SFP.dvir1.06.orths, dvir1.06_BlastX_topHit == "NoHit")$gene_id
vir.SFP_no_dvir1.06_hits = unique(strsplit(vir.SFP_no_dvir1.06_hits, ", ")[[1]])
vir.SFP_no_dvir1.06_hits_Trin_hits = as.character(unique(subset(TrinOrths, Gene %in% vir.SFP_no_dvir1.06_hits)$Gene))
setdiff(vir.SFP_no_dvir1.06_hits, vir.SFP_no_dvir1.06_hits_Trin_hits)

########################################################################

# D. americana AGs
amr.AG.dvir1.06.orths = subset(amrTrinotate, gene_id %in% amr.AG.list)[selectionCols]
amr.AG.dvir1.06.orths = droplevels(amr.AG.dvir1.06.orths)
amr.AG.dvir1.06.orths = amr.AG.dvir1.06.orths[order(amr.AG.dvir1.06.orths$dvir1.06_BlastX_topHit), ]
amr.AG.dvir1.06.orths[is.na(amr.AG.dvir1.06.orths)] = "NoHit"
amr.AG.dvir1.06.orths = subset(amr.AG.dvir1.06.orths, prot_id != "NoHit")
amr.AG.dvir1.06.orths = aggregate(gene_id~dvir1.06_BlastX_topHit, data = amr.AG.dvir1.06.orths, toString)
amr.AG_no_dvir1.06_hits = subset(amr.AG.dvir1.06.orths, dvir1.06_BlastX_topHit == "NoHit")$gene_id
amr.AG_no_dvir1.06_hits = unique(strsplit(amr.AG_no_dvir1.06_hits, ", ")[[1]])
amr.AG_no_dvir1.06_hits_Trin_hits = as.character(unique(subset(TrinOrths, Gene %in% amr.AG_no_dvir1.06_hits)$Gene))
tmpAGamrList = setdiff(amr.AG_no_dvir1.06_hits, amr.AG_no_dvir1.06_hits_Trin_hits)


# D. lummei AGs
lum.AG.dvir1.06.orths = subset(lumTrinotate, gene_id %in% lum.AG.list)[selectionCols]
lum.AG.dvir1.06.orths = droplevels(lum.AG.dvir1.06.orths)
lum.AG.dvir1.06.orths = lum.AG.dvir1.06.orths[order(lum.AG.dvir1.06.orths$dvir1.06_BlastX_topHit), ]
lum.AG.dvir1.06.orths[is.na(lum.AG.dvir1.06.orths)] = "NoHit"
lum.AG.dvir1.06.orths = aggregate(gene_id~dvir1.06_BlastX_topHit, data = lum.AG.dvir1.06.orths, toString)
lum.AG_no_dvir1.06_hits = subset(lum.AG.dvir1.06.orths, dvir1.06_BlastX_topHit == "NoHit")$gene_id
lum.AG_no_dvir1.06_hits = unique(strsplit(lum.AG_no_dvir1.06_hits, ", ")[[1]])
lum.AG_no_dvir1.06_hits_Trin_hits = as.character(unique(subset(TrinOrths, Gene %in% lum.AG_no_dvir1.06_hits)$Gene))
setdiff(lum.AG_no_dvir1.06_hits, lum.AG_no_dvir1.06_hits_Trin_hits)

# D. novamexicana AGs
nov.AG.dvir1.06.orths = subset(novTrinotate, gene_id %in% nov.AG.list)[selectionCols]
nov.AG.dvir1.06.orths = droplevels(nov.AG.dvir1.06.orths)
nov.AG.dvir1.06.orths = nov.AG.dvir1.06.orths[order(nov.AG.dvir1.06.orths$dvir1.06_BlastX_topHit), ]
nov.AG.dvir1.06.orths[is.na(nov.AG.dvir1.06.orths)] = "NoHit"
nov.AG.dvir1.06.orths = aggregate(gene_id~dvir1.06_BlastX_topHit, data = nov.AG.dvir1.06.orths, toString)
nov.AG_no_dvir1.06_hits = subset(nov.AG.dvir1.06.orths, dvir1.06_BlastX_topHit == "NoHit")$gene_id
nov.AG_no_dvir1.06_hits = unique(strsplit(nov.AG_no_dvir1.06_hits, ", ")[[1]])
nov.AG_no_dvir1.06_hits_Trin_hits = as.character(unique(subset(TrinOrths, Gene %in% nov.AG_no_dvir1.06_hits)$Gene))
setdiff(nov.AG_no_dvir1.06_hits, nov.AG_no_dvir1.06_hits_Trin_hits)

# D. virilis AGs
vir.AG.dvir1.06.orths = subset(virTrinotate, gene_id %in% vir.AG.list)[selectionCols]
vir.AG.dvir1.06.orths = droplevels(vir.AG.dvir1.06.orths)
vir.AG.dvir1.06.orths = vir.AG.dvir1.06.orths[order(vir.AG.dvir1.06.orths$dvir1.06_BlastX_topHit), ]
vir.AG.dvir1.06.orths[is.na(vir.AG.dvir1.06.orths)] = "NoHit"
vir.AG.dvir1.06.orths = aggregate(gene_id~dvir1.06_BlastX_topHit, data = vir.AG.dvir1.06.orths, toString)
vir.AG_no_dvir1.06_hits = subset(vir.AG.dvir1.06.orths, dvir1.06_BlastX_topHit == "NoHit")$gene_id
vir.AG_no_dvir1.06_hits = unique(strsplit(vir.AG_no_dvir1.06_hits, ", ")[[1]])
vir.AG_no_dvir1.06_hits_Trin_hits = as.character(unique(subset(TrinOrths, Gene %in% vir.AG_no_dvir1.06_hits)$Gene))
setdiff(vir.AG_no_dvir1.06_hits, vir.AG_no_dvir1.06_hits_Trin_hits)




########################################################################
### Distribution of tissue-biased transcripts on chromosomes ###########
TotalGeneNumber = as.data.frame(table(factor(subset(gffRecord, grepl("Chr", chromosome))$chromosome)))
colnames(TotalGeneNumber) = c("chromosome", "All genes")
total_genes = nrow(gffRecord)
TotalGeneNumber$proportion = (TotalGeneNumber$`All genes`/total_genes)

## for chromosome distribution plots
# americana
genomeNumber.amrAG = length(amr.dvir1.06.AG.list)
chromNumber.amrAG=as.data.frame(table(factor(subset(gffRecord, FBgn_ID %in% amr.dvir1.06.AG.list & grepl("Chr", chromosome) & chromosome != "Chr_6")$chromosome)))
colnames(chromNumber.amrAG) = c("chromosome", "Observed_biased_genes")
chromNumber.amrAG$species = "D. americana"
chromNumber.amrAG$tissue = "Accessory glands"
chromNumber.amrAG = merge(TotalGeneNumber, chromNumber.amrAG)
chromNumber.amrAG$`Expected genes` = genomeNumber.amrAG*chromNumber.amrAG$proportion
chromNumber.amrAG$`obs.exp` = chromNumber.amrAG$Observed_biased_genes/chromNumber.amrAG$`Expected genes`

genomeNumber.amrSFP = length(amr.dvir1.06.SFP.list)
chromNumber.amrSFP=as.data.frame(table(factor(subset(gffRecord, FBgn_ID %in% amr.dvir1.06.SFP.list & grepl("Chr", chromosome) & chromosome != "Chr_6")$chromosome)))
colnames(chromNumber.amrSFP) = c("chromosome", "Observed_biased_genes")
chromNumber.amrSFP$species = "D. americana"
chromNumber.amrSFP$tissue = "SFPs"
chromNumber.amrSFP = merge(TotalGeneNumber, chromNumber.amrSFP)
chromNumber.amrSFP$`Expected genes` = genomeNumber.amrSFP*chromNumber.amrSFP$proportion
chromNumber.amrSFP$`obs.exp` = chromNumber.amrSFP$Observed_biased_genes/chromNumber.amrSFP$`Expected genes`

genomeNumber.amrEB = length(amr.dvir1.06.EB.list)
chromNumber.amrEB=as.data.frame(table(factor(subset(gffRecord, FBgn_ID %in% amr.dvir1.06.EB.list & grepl("Chr", chromosome) & chromosome != "Chr_6")$chromosome)))
colnames(chromNumber.amrEB) = c("chromosome", "Observed_biased_genes")
chromNumber.amrEB$species = "D. americana"
chromNumber.amrEB$tissue = "Ejaculatory Bulb"
chromNumber.amrEB = merge(TotalGeneNumber, chromNumber.amrEB)
chromNumber.amrEB$`Expected genes` = genomeNumber.amrEB*chromNumber.amrEB$proportion
chromNumber.amrEB$`obs.exp` = chromNumber.amrEB$Observed_biased_genes/chromNumber.amrEB$`Expected genes`
genomeNumber.amrTS = length(amr.dvir1.06.TS.list)
chromNumber.amrTS=as.data.frame(table(factor(subset(gffRecord, FBgn_ID %in% amr.dvir1.06.TS.list & grepl("Chr", chromosome) & chromosome != "Chr_6")$chromosome)))
colnames(chromNumber.amrTS) = c("chromosome", "Observed_biased_genes")
chromNumber.amrTS$species = "D. americana"
chromNumber.amrTS$tissue = "Testes"
chromNumber.amrTS = merge(TotalGeneNumber, chromNumber.amrTS)
chromNumber.amrTS$`Expected genes` = genomeNumber.amrTS*chromNumber.amrTS$proportion
chromNumber.amrTS$`obs.exp` = chromNumber.amrTS$Observed_biased_genes/chromNumber.amrTS$`Expected genes`
# lummei
genomeNumber.lumAG = length(lum.dvir1.06.AG.list)
chromNumber.lumAG=as.data.frame(table(factor(subset(gffRecord, FBgn_ID %in% lum.dvir1.06.AG.list & grepl("Chr", chromosome) & chromosome != "Chr_6")$chromosome)))
colnames(chromNumber.lumAG) = c("chromosome", "Observed_biased_genes")
chromNumber.lumAG$species = "D. lummei"
chromNumber.lumAG$tissue = "Accessory glands"
chromNumber.lumAG = merge(TotalGeneNumber, chromNumber.lumAG)
chromNumber.lumAG$`Expected genes` = genomeNumber.lumAG*chromNumber.lumAG$proportion
chromNumber.lumAG$`obs.exp` = chromNumber.lumAG$Observed_biased_genes/chromNumber.lumAG$`Expected genes`

genomeNumber.lumSFP = length(lum.dvir1.06.SFP.list)
chromNumber.lumSFP=as.data.frame(table(factor(subset(gffRecord, FBgn_ID %in% lum.dvir1.06.SFP.list & grepl("Chr", chromosome) & chromosome != "Chr_6")$chromosome)))
colnames(chromNumber.lumSFP) = c("chromosome", "Observed_biased_genes")
chromNumber.lumSFP$species = "D. lummei"
chromNumber.lumSFP$tissue = "SFPs"
chromNumber.lumSFP = merge(TotalGeneNumber, chromNumber.lumSFP)
chromNumber.lumSFP$`Expected genes` = genomeNumber.lumSFP*chromNumber.lumSFP$proportion
chromNumber.lumSFP$`obs.exp` = chromNumber.lumSFP$Observed_biased_genes/chromNumber.lumSFP$`Expected genes`

genomeNumber.lumEB = length(lum.dvir1.06.EB.list)
chromNumber.lumEB=as.data.frame(table(factor(subset(gffRecord, FBgn_ID %in% lum.dvir1.06.EB.list & grepl("Chr", chromosome) & chromosome != "Chr_6")$chromosome)))
colnames(chromNumber.lumEB) = c("chromosome", "Observed_biased_genes")
chromNumber.lumEB$species = "D. lummei"
chromNumber.lumEB$tissue = "Ejaculatory Bulb"
chromNumber.lumEB = merge(TotalGeneNumber, chromNumber.lumEB)
chromNumber.lumEB$`Expected genes` = genomeNumber.lumEB*chromNumber.lumEB$proportion
chromNumber.lumEB$`obs.exp` = chromNumber.lumEB$Observed_biased_genes/chromNumber.lumEB$`Expected genes`
genomeNumber.lumTS = length(lum.dvir1.06.TS.list)
chromNumber.lumTS=as.data.frame(table(factor(subset(gffRecord, FBgn_ID %in% lum.dvir1.06.TS.list & grepl("Chr", chromosome) & chromosome != "Chr_6")$chromosome)))
colnames(chromNumber.lumTS) = c("chromosome", "Observed_biased_genes")
chromNumber.lumTS$species = "D. lummei"
chromNumber.lumTS$tissue = "Testes"
chromNumber.lumTS = merge(TotalGeneNumber, chromNumber.lumTS)
chromNumber.lumTS$`Expected genes` = genomeNumber.lumTS*chromNumber.lumTS$proportion
chromNumber.lumTS$`obs.exp` = chromNumber.lumTS$Observed_biased_genes/chromNumber.lumTS$`Expected genes`
# novamexicana
genomeNumber.novAG = length(nov.dvir1.06.AG.list)
chromNumber.novAG=as.data.frame(table(factor(subset(gffRecord, FBgn_ID %in% nov.dvir1.06.AG.list & grepl("Chr", chromosome) & chromosome != "Chr_6")$chromosome)))
colnames(chromNumber.novAG) = c("chromosome", "Observed_biased_genes")
chromNumber.novAG$species = "D. novamexicana"
chromNumber.novAG$tissue = "Accessory glands"
chromNumber.novAG = merge(TotalGeneNumber, chromNumber.novAG)
chromNumber.novAG$`Expected genes` = genomeNumber.novAG*chromNumber.novAG$proportion
chromNumber.novAG$`obs.exp` = chromNumber.novAG$Observed_biased_genes/chromNumber.novAG$`Expected genes`

genomeNumber.novSFP = length(nov.dvir1.06.SFP.list)
chromNumber.novSFP=as.data.frame(table(factor(subset(gffRecord, FBgn_ID %in% nov.dvir1.06.SFP.list & grepl("Chr", chromosome) & chromosome != "Chr_6")$chromosome)))
colnames(chromNumber.novSFP) = c("chromosome", "Observed_biased_genes")
chromNumber.novSFP$species = "D. novamexicana"
chromNumber.novSFP$tissue = "SFPs"
chromNumber.novSFP = merge(TotalGeneNumber, chromNumber.novSFP)
chromNumber.novSFP$`Expected genes` = genomeNumber.novSFP*chromNumber.novSFP$proportion
chromNumber.novSFP$`obs.exp` = chromNumber.novSFP$Observed_biased_genes/chromNumber.novSFP$`Expected genes`

genomeNumber.novEB = length(nov.dvir1.06.EB.list)
chromNumber.novEB=as.data.frame(table(factor(subset(gffRecord, FBgn_ID %in% nov.dvir1.06.EB.list & grepl("Chr", chromosome) & chromosome != "Chr_6")$chromosome)))
colnames(chromNumber.novEB) = c("chromosome", "Observed_biased_genes")
chromNumber.novEB$species = "D. novamexicana"
chromNumber.novEB$tissue = "Ejaculatory Bulb"
chromNumber.novEB = merge(TotalGeneNumber, chromNumber.novEB)
chromNumber.novEB$`Expected genes` = genomeNumber.novEB*chromNumber.novEB$proportion
chromNumber.novEB$`obs.exp` = chromNumber.novEB$Observed_biased_genes/chromNumber.novEB$`Expected genes`
genomeNumber.novTS = length(nov.dvir1.06.TS.list)
chromNumber.novTS=as.data.frame(table(factor(subset(gffRecord, FBgn_ID %in% nov.dvir1.06.TS.list & grepl("Chr", chromosome) & chromosome != "Chr_6")$chromosome)))
colnames(chromNumber.novTS) = c("chromosome", "Observed_biased_genes")
chromNumber.novTS$species = "D. novamexicana"
chromNumber.novTS$tissue = "Testes"
chromNumber.novTS = merge(TotalGeneNumber, chromNumber.novTS)
chromNumber.novTS$`Expected genes` = genomeNumber.novTS*chromNumber.novTS$proportion
chromNumber.novTS$`obs.exp` = chromNumber.novTS$Observed_biased_genes/chromNumber.novTS$`Expected genes`
# virilis
genomeNumber.virAG = length(vir.dvir1.06.AG.list)
chromNumber.virAG=as.data.frame(table(factor(subset(gffRecord, FBgn_ID %in% vir.dvir1.06.AG.list & grepl("Chr", chromosome) & chromosome != "Chr_6")$chromosome)))
colnames(chromNumber.virAG) = c("chromosome", "Observed_biased_genes")
chromNumber.virAG$species = "D. virilis"
chromNumber.virAG$tissue = "Accessory glands"
chromNumber.virAG = merge(TotalGeneNumber, chromNumber.virAG)
chromNumber.virAG$`Expected genes` = genomeNumber.virAG*chromNumber.virAG$proportion
chromNumber.virAG$`obs.exp` = chromNumber.virAG$Observed_biased_genes/chromNumber.virAG$`Expected genes`

genomeNumber.virSFP = length(vir.dvir1.06.SFP.list)
chromNumber.virSFP=as.data.frame(table(factor(subset(gffRecord, FBgn_ID %in% vir.dvir1.06.SFP.list & grepl("Chr", chromosome) & chromosome != "Chr_6")$chromosome)))
colnames(chromNumber.virSFP) = c("chromosome", "Observed_biased_genes")
chromNumber.virSFP$species = "D. virilis"
chromNumber.virSFP$tissue = "SFPs"
chromNumber.virSFP = merge(TotalGeneNumber, chromNumber.virSFP)
chromNumber.virSFP$`Expected genes` = genomeNumber.virSFP*chromNumber.virSFP$proportion
chromNumber.virSFP$`obs.exp` = chromNumber.virSFP$Observed_biased_genes/chromNumber.virSFP$`Expected genes`

genomeNumber.virEB = length(vir.dvir1.06.EB.list)
chromNumber.virEB=as.data.frame(table(factor(subset(gffRecord, FBgn_ID %in% vir.dvir1.06.EB.list & grepl("Chr", chromosome) & chromosome != "Chr_6")$chromosome)))
colnames(chromNumber.virEB) = c("chromosome", "Observed_biased_genes")
chromNumber.virEB$species = "D. virilis"
chromNumber.virEB$tissue = "Ejaculatory Bulb"
chromNumber.virEB = merge(TotalGeneNumber, chromNumber.virEB)
chromNumber.virEB$`Expected genes` = genomeNumber.virEB*chromNumber.virEB$proportion
chromNumber.virEB$`obs.exp` = chromNumber.virEB$Observed_biased_genes/chromNumber.virEB$`Expected genes`
genomeNumber.virTS = length(vir.dvir1.06.TS.list)
chromNumber.virTS=as.data.frame(table(factor(subset(gffRecord, FBgn_ID %in% vir.dvir1.06.TS.list & grepl("Chr", chromosome) & chromosome != "Chr_6")$chromosome)))
colnames(chromNumber.virTS) = c("chromosome", "Observed_biased_genes")
chromNumber.virTS$species = "D. virilis"
chromNumber.virTS$tissue = "Testes"
chromNumber.virTS = merge(TotalGeneNumber, chromNumber.virTS)
chromNumber.virTS$`Expected genes` = genomeNumber.virTS*chromNumber.virTS$proportion
chromNumber.virTS$`obs.exp` = chromNumber.virTS$Observed_biased_genes/chromNumber.virTS$`Expected genes`
## put it all togetehr
tissue_biased.numbers = rbind(chromNumber.amrSFP, chromNumber.amrAG, chromNumber.amrEB, chromNumber.amrTS, chromNumber.lumSFP, chromNumber.lumAG, chromNumber.lumEB, chromNumber.lumTS, chromNumber.novSFP, chromNumber.novAG, chromNumber.novEB, chromNumber.novTS, chromNumber.virSFP, chromNumber.virAG, chromNumber.virEB, chromNumber.virTS)

tissue_biased.numbers.c = summarySE(tissue_biased.numbers, measurevar = "obs.exp", groupvars = c("chromosome", "tissue"))
tissue_biased.numbers.c$chromosome = factor (tissue_biased.numbers.c$chromosome, levels = c("Chr_X", "Chr_2", "Chr_3", "Chr_4", "Chr_5", "Chr_6"))
tissue_biased.numbers.c2 = summarySE(tissue_biased.numbers, measurevar = "Observed_biased_genes", groupvars = c("chromosome", "tissue"))
tissue_biased.numbers.c3 = summarySE(tissue_biased.numbers, measurevar = "Expected genes", groupvars = c("chromosome", "tissue"))
dataChiSq = data.frame(tissue_biased.numbers.c2$chromosome, tissue_biased.numbers.c2$tissue, tissue_biased.numbers.c2$Observed_biased_genes, tissue_biased.numbers.c3$`Expected genes`, tissue_biased.numbers.c$obs.exp+0.2)
colnames(dataChiSq) = c("chromosome", "tissue", "obs", "exp", "obs.exp")
dataChiSq$ChiSq = ((dataChiSq$obs-dataChiSq$exp)^2)/dataChiSq$exp
dataChiSq$pval = 1-(pchisq(dataChiSq$ChiSq, df = 1))

label.05 = subset(dataChiSq, pval<0.05 & pval > 0.01)
label.01 = subset(dataChiSq, pval<0.01 & pval > 0.001)
label.001 = subset(dataChiSq, pval<0.001)

tissue_biased.numbers.c$tissue = factor(tissue_biased.numbers.c$tissue, levels = c("SFPs", "Accessory glands", "Ejaculatory Bulb", "Testes"))

pdf(file = "ManuscripPlots/Figure.NA.chromDist.pdf", width = 8.10, height = 2.4)
ggplot(tissue_biased.numbers.c, aes(x=chromosome, y=obs.exp, fill=chromosome)) +
  geom_bar(position=position_dodge(), stat="identity") +
  scale_fill_manual(values=c("#CC79A7", "#56B4E9", "#56B4E9", "#56B4E9", "#56B4E9")) +
  facet_wrap(~tissue, nrow = 1) +
  geom_errorbar(aes(ymin=obs.exp-ci, ymax=obs.exp+ci), width=.3, position=position_dodge(.9), size = 0.75) +
  ylab("Observed/expected number\n of genes per chromosome") +
  geom_hline(yintercept = 1, colour = "black", linetype = "dashed") +
  scale_x_discrete(labels=c("Chr_X" = "X", "Chr_2" = "2", "Chr_3" = "3", "Chr_4" = "4", "Chr_5" = "5")) +
  theme(panel.background = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(face = "bold", size = 14, vjust=0.1), axis.text.y = element_text(face = "bold", size = 12), legend.position="none", strip.text=element_text(face="bold", size = 12), axis.title=element_text(face="bold", size = 12)) +
  geom_text(data = label.05, label = "*", size = 6, colour = "red", fontface=2) +
  geom_text(data = label.01, label = "**", size = 6, colour = "red", fontface=2) +
  geom_text(data = label.001, label = "***", size = 6, colour = "red", fontface=2)
dev.off()


##########################################################################################
##########################################################################################
##########################################################################################

#### Subset PAML data by tissue and significant brSt genes
AG.paml.data = subset(paml.data, FBgn_ID %in% AG_elements$`D.ame,D.lum,D.nov,D.vir`)
AG.paml.data.sig = subset(AG.paml.data, Damr_FDR < 0.05 | Dlum_FDR < 0.05 | Dnov_FDR < 0.05 | Dvir_FDR < 0.05)
AG.paml.data.sig.list = as.character(unique(AG.paml.data.sig$FBgn_ID))

SFP.paml.data = subset(paml.data, FBgn_ID %in% SFP_elements$`D.ame,D.lum,D.nov,D.vir`)
SFP.paml.data.sig = subset(SFP.paml.data, Damr_FDR < 0.05 | Dlum_FDR < 0.05 | Dnov_FDR < 0.05 | Dvir_FDR < 0.05)
SFP.paml.data.sig.list = as.character(unique(SFP.paml.data.sig$FBgn_ID))

EB.paml.data = subset(paml.data, FBgn_ID %in% EB_elements$`D.ame,D.lum,D.nov,D.vir`)
EB.paml.data.sig = subset(EB.paml.data, Damr_FDR < 0.05 | Dlum_FDR < 0.05 | Dnov_FDR < 0.05 | Dvir_FDR < 0.05)
EB.paml.data.sig.list = as.character(unique(EB.paml.data.sig$FBgn_ID))

TS.paml.data = subset(paml.data, FBgn_ID %in% TS_elements$`D.ame,D.lum,D.nov,D.vir`)
TS.paml.data.sig = subset(TS.paml.data, Damr_FDR < 0.05 | Dlum_FDR < 0.05 | Dnov_FDR < 0.05 | Dvir_FDR < 0.05)
TS.paml.data.sig.list = as.character(unique(TS.paml.data.sig$FBgn_ID))

genome.omega = subset(paml.data, omega != "NA" & omega < 100)
genome.omega = subset(genome.omega, select=c(FBtr_ID, gene_name, chromosome, omega))
genome.omega$Class = "All genes"

AG.omega = subset(paml.data, FBgn_ID %in% AG_elements$`D.ame,D.lum,D.nov,D.vir` & omega != "NA" & omega < 100)
AG.omega = subset(AG.omega, select=c(FBtr_ID, gene_name, chromosome, omega))
AG.omega$Class = "AG biased"

EB.omega = subset(paml.data, FBgn_ID %in% EB_elements$`D.ame,D.lum,D.nov,D.vir` & omega != "NA" & omega < 100)
EB.omega = subset(EB.omega, select=c(FBtr_ID, gene_name, chromosome, omega))
EB.omega$Class = "EB biased"

TS.omega = subset(paml.data, FBgn_ID %in% TS_elements$`D.ame,D.lum,D.nov,D.vir` & omega != "NA" & omega < 100)
TS.omega = subset(TS.omega, select=c(FBtr_ID, gene_name, chromosome, omega))
TS.omega$Class = "testes biased"

SFP.omega = subset(paml.data, FBgn_ID %in% SFP_elements$`D.ame,D.lum,D.nov,D.vir` & omega != "NA" & omega < 100)
SFP.omega = subset(SFP.omega, select=c(FBtr_ID, gene_name, chromosome, omega))
SFP.omega$Class = "SFPs"

omegaData.df = rbind(genome.omega, AG.omega, EB.omega, TS.omega, SFP.omega)
omegaData.df$Class = factor (omegaData.df$Class, levels = c("All genes", "EB biased", "testes biased", "AG biased", "SFPs"))

genome.avg.omega = mean(subset(omegaData.df, Class == "All genes" & omega != "NA" & omega < 100)$omega)
genome.se.omega = sd(subset(omegaData.df, Class == "All genes" & omega != "NA" & omega < 100)$omega)/sqrt(length(subset(omegaData.df, Class == "All genes" & omega != "NA" & omega < 100)$omega))

EB.avg.omega = mean(subset(omegaData.df, Class == "EB biased" & omega != "NA" & omega < 100)$omega)
EB.se.omega = sd(subset(omegaData.df, Class == "EB biased" & omega != "NA" & omega < 100)$omega)/sqrt(length(subset(omegaData.df, Class == "EB biased" & omega != "NA" & omega < 100)$omega))

AG.avg.omega = mean(subset(omegaData.df, Class == "AG biased" & omega != "NA" & omega < 100)$omega)
AG.se.omega = sd(subset(omegaData.df, Class == "AG biased" & omega != "NA" & omega < 100)$omega)/sqrt(length(subset(omegaData.df, Class == "AG biased" & omega != "NA" & omega < 100)$omega))

TS.avg.omega = mean(subset(omegaData.df, Class == "testes biased" & omega != "NA" & omega < 100)$omega)
TS.se.omega = sd(subset(omegaData.df, Class == "testes biased" & omega != "NA" & omega < 100)$omega)/sqrt(length(subset(omegaData.df, Class == "testes biased" & omega != "NA" & omega < 100)$omega))

SFP.avg.omega = mean(subset(omegaData.df, Class == "SFPs" & omega != "NA" & omega < 100)$omega)
SFP.se.omega = sd(subset(omegaData.df, Class == "SFPs" & omega != "NA" & omega < 100)$omega)/sqrt(length(subset(omegaData.df, Class == "SFPs" & omega != "NA" & omega < 100)$omega))

meanOmega.df=data.frame(Class = c("All", "EB", "testes", "AG", "SFPs"), omega=c(genome.avg.omega, EB.avg.omega, TS.avg.omega, AG.avg.omega, SFP.avg.omega), se = c(genome.se.omega, EB.se.omega, TS.se.omega, AG.se.omega, SFP.se.omega))

meanOmega.df$Class = factor (meanOmega.df$Class, levels = c("All", "EB", "testes", "AG", "SFPs"))

pdf("ManuscripPlots/Figure.NA.meanOmega.pdf", width = 4.25, height = 2.25)
ggplot(meanOmega.df, aes(Class, omega, colour=Class)) + 
  geom_point(size = 2) + 
  geom_errorbar(aes(ymin=omega-se, ymax=omega+se), width=.2, position=position_dodge(.9)) + 
  scale_colour_manual(values=c("black", "#b04fb0", "#b85516", "#f0cc35", "#647700")) + 
  theme_bw() + 
  theme(legend.position="none") +
  ylab(expression(omega)) + 
  xlab("Gene class")
dev.off()

paml.data$chromosome = factor (paml.data$chromosome, levels = c("Chr_X", "Chr_2", "Chr_3", "Chr_4", "Chr_5"))

## Define PMPZ QTL regions (from Ahmed-Braimah 2016, G3)
C2inv.qtl = data.frame(xmin=17747413.5, xmax=34500000, ymin=0, ymax = 2.5, chromosome = "Chr_2")
C5.qtl = data.frame(xmin=c(13800000, 16300000, 22800000), xmax=c(14750000, 21700000, 25000000), ymin=0, ymax = 2.5, chromosome = "Chr_5")

(gg.SFP_and_EB=ggplot(subset(paml.data, FBgn_ID %in% SFP_elements$`D.ame,D.lum,D.nov,D.vir` ), aes(max, omega, colour = I("#647700"))) + 
  geom_hline(yintercept = genome.avg.omega, linetype="dashed", colour = "red") + 
  geom_point(size=2, alpha=0.75) + 
  facet_grid(~chromosome, scales = "free_x") + 
  geom_point(data = subset(paml.data, FBgn_ID %in% EB_elements$`D.ame,D.lum,D.nov,D.vir`& chromosome != "scaffold_12481"), aes(colour = I("#b04fb0"))) + 
  geom_text_repel(data=subset(paml.data, FBgn_ID %in% SFP_elements$`D.ame,D.lum,D.nov,D.vir` & omega > 0.8), aes(label = gene_name), size =3, force = 30, colour = "black") + 
  scale_colour_manual(name = "", values =c("#647700"="#647700","#b04fb0"="#b04fb0"), labels = c("SFPs","EB biased")) + 
  scale_x_continuous(breaks=c(5000000, 10000000, 15000000, 20000000, 25000000, 30000000), labels=expression("5", "10", "15", "20", "25", "30")) + 
  xlab ("Chromosome coordinates (Mb)") +
  ylab(expression(omega)) + 
  geom_rect(data=C2inv.qtl, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill="red"), alpha=0.1, inherit.aes = FALSE) + 
  geom_rect(data=C5.qtl, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="red", alpha=0.1, inherit.aes = FALSE) + 
  scale_fill_manual(name = '', values = "red",labels = "PMPZ\nQTL region") +
  theme_bw())

(gg.AG=ggplot(subset(paml.data, FBgn_ID %in% AG_elements$`D.ame,D.lum,D.nov,D.vir` & omega < 900 & grepl("Chr", chromosome) & FBgn_ID %!in% SFP_elements$`D.ame,D.lum,D.nov,D.vir`), aes(max, omega, colour = I("#f0cc35"))) + 
  geom_point(alpha=0.55) + 
  facet_grid(~chromosome, scale = "free_x") + 
  scale_x_continuous(breaks=c(5000000, 10000000, 15000000, 20000000, 25000000, 30000000), labels=expression("5", "10", "15", "20", "25", "30")) + 
  xlab ("Chromosome coordinates (Mb)") + 
  scale_colour_manual(name = "", values ="#f0cc35", labels = "Acc. Gl.\nbiased") + 
  scale_y_continuous(breaks=c(0.0,0.5,1.0,1.5,2.0,2.5), labels=expression("0.0", "0.5", "1.0", "1.5", "2.0", "2.5")) + 
  geom_hline(yintercept = genome.avg.omega, linetype="dashed", colour = "red") +
  ylab(expression(omega)) +
  geom_rect(data=C2inv.qtl, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill="red"), alpha=0.1, inherit.aes = FALSE) + 
  geom_rect(data=C5.qtl, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="red", alpha=0.1, inherit.aes = FALSE) + 
  scale_fill_manual(name = '', values = "red",labels = "PMPZ\nQTL region") +
  theme_bw())

(gg.TS=ggplot(subset(paml.data, FBgn_ID %in% TS_elements$`D.ame,D.lum,D.nov,D.vir` & omega < 900 & grepl("Chr", chromosome)), aes(max, omega, colour = I("#b85516"))) + 
  geom_point(alpha=0.55) + 
  facet_grid(~chromosome, scale = "free_x") + 
  scale_x_continuous(breaks=c(5000000, 10000000, 15000000, 20000000, 25000000, 30000000), labels=expression("5", "10", "15", "20", "25", "30")) + 
  xlab ("Chromosome coordinates (Mb)") + 
  scale_colour_manual(name = "", values ="#cb6751", labels = "Testes\nbiased") + 
  scale_y_continuous(breaks=c(0.0,0.5,1.0,1.5,2.0,2.5), labels=expression("0.0", "0.5", "1.0", "1.5", "2.0", "2.5")) + 
  geom_hline(yintercept = genome.avg.omega, linetype="dashed", colour = "red") +
  ylab(expression(omega)) +
  geom_rect(data=C2inv.qtl, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill="red"), alpha=0.1, inherit.aes = FALSE) + 
  geom_rect(data=C5.qtl, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="red", alpha=0.1, inherit.aes = FALSE) + 
  scale_fill_manual(name = '', values = "red",labels = "PMPZ\nQTL region") +
  theme_bw())

pdf("ManuscripPlots/Figure.NA.omegaAcrossChrs.pdf", width = 8.5, height = 5.6)
plot_grid(gg.SFP_and_EB, gg.AG, gg.TS, ncol = 1)
dev.off()

###################
##################

### PAML FDR across chromosomes

SFP.subset.paml = subset(paml.data, FBgn_ID %in% SFP_elements$`D.ame,D.lum,D.nov,D.vir`)
SFP.subset.paml$tissue = "SFPs"
AG.subset.paml = subset(paml.data, FBgn_ID %in% AG_elements$`D.ame,D.lum,D.nov,D.vir`)
AG.subset.paml$tissue = "AG-biased"
EB.subset.paml = subset(paml.data, FBgn_ID %in% EB_elements$`D.ame,D.lum,D.nov,D.vir`)
EB.subset.paml$tissue = "EB-biased"
TS.subset.paml = subset(paml.data, FBgn_ID %in% TS_elements$`D.ame,D.lum,D.nov,D.vir`)
TS.subset.paml$tissue = "TS-biased"
subset.paml = rbind(SFP.subset.paml, AG.subset.paml, EB.subset.paml, TS.subset.paml)
subset.paml = subset(subset.paml, select=c("FBgn_ID","gene_name","chromosome","min","Damr_FDR","Dlum_FDR","Dnov_FDR","Dvir_FDR","tissue"))
subset.paml.m = melt(subset.paml, id.vars = c("FBgn_ID","gene_name","chromosome","min","tissue"))  


pdf("ManuscripPlots/Figure.NA.PAMLfdr.chrom.pdf", width = 9, height = 5)
ggplot(subset(subset.paml.m, grepl("Chr", subset.paml.m$chromosome)), aes(min, -log10(value), colour = variable)) + 
  geom_point(size =2, alpha=0.5) +
  facet_grid(tissue~chromosome, scales = "free")+
  geom_hline(yintercept = 1.3, linetype="dashed", colour = "purple") +
  geom_text_repel(data=subset(subset.paml.m, value < 0.05& grepl("Chr", subset.paml.m$chromosome)), 
                  aes(label = gene_name, colour = variable), size =3, force = 30) +
  scale_x_continuous(breaks=c(5000000, 10000000, 15000000, 20000000, 25000000, 30000000), 
                     labels=expression("5", "10", "15", "20", "25", "30")) + 
  xlab ("Chromosome coordinates (Mb)") +
  ylab("-log10(FDR)") +
  scale_colour_manual(values=c("#b38c3a", "#8f73c9", "#5ea46d", "#ca587a")) + 
  theme_bw()
dev.off()


### PAML PValues across chromosomes
SFP.subset.paml = subset(paml.data, FBgn_ID %in% SFP_elements$`D.ame,D.lum,D.nov,D.vir`)
SFP.subset.paml$tissue = "SFPs"
AG.subset.paml = subset(paml.data, FBgn_ID %in% AG_elements$`D.ame,D.lum,D.nov,D.vir`)
AG.subset.paml$tissue = "AG-biased"
EB.subset.paml = subset(paml.data, FBgn_ID %in% EB_elements$`D.ame,D.lum,D.nov,D.vir`)
EB.subset.paml$tissue = "EB-biased"
TS.subset.paml = subset(paml.data, FBgn_ID %in% TS_elements$`D.ame,D.lum,D.nov,D.vir`)
TS.subset.paml$tissue = "TS-biased"
subset.paml = rbind(SFP.subset.paml, AG.subset.paml, EB.subset.paml, TS.subset.paml)
subset.paml = subset(subset.paml, select=c("FBgn_ID","gene_name","chromosome","min","Damr_pValue","Dlum_pValue","Dnov_pValue","Dvir_pValue","tissue"))
subset.paml.m = melt(subset.paml, id.vars = c("FBgn_ID","gene_name","chromosome","min","tissue"))  

ggplot(subset(subset.paml.m, grepl("Chr", subset.paml.m$chromosome)), aes(min, -log10(value), colour = variable)) + 
  geom_point(size =2, alpha=0.5) +
  facet_grid(tissue~chromosome, scales = "free")+
  geom_hline(yintercept = 1.3, linetype="dashed", colour = "purple") +
  geom_text_repel(data=subset(subset.paml.m, tissue == "AG-biased" & value < 0.0001 & grepl("Chr", subset.paml.m$chromosome)),
                  aes(label = gene_name, colour = variable), size =3, force = 30) +
  scale_x_continuous(breaks=c(5000000, 10000000, 15000000, 20000000, 25000000, 30000000), 
                     labels=expression("5", "10", "15", "20", "25", "30")) + 
  xlab ("Chromosome coordinates (Mb)") +
  ylab("-log10(FDR)") +
  scale_colour_manual(values=c("#b38c3a", "#8f73c9", "#5ea46d", "#ca587a")) + 
  theme_bw()


###################
##################

### pairwise Ka/Ks across chromosomes
KaKs.data$chromosome = factor (KaKs.data$chromosome, levels = c("Chr_X", "Chr_2", "Chr_3", "Chr_4", "Chr_5"))

pdf("ManuscripPlots/Figure.NA.pwKaKs.SFPs.pdf", width = 9, height = 7)
ggplot(subset(KaKs.data, FBgn_ID %in% SFP_elements$`D.ame,D.lum,D.nov,D.vir` & `Ka/Ks` < 80), aes(max, `Ka/Ks`, colour = I("#647700"))) + 
  geom_hline(yintercept = genome.avg.omega, linetype="dashed", colour = "red") + 
  geom_point(aes(size = `Ka/Ks`), alpha=0.75) + 
  facet_grid(COMPARISON~chromosome, scales = "free_x") + 
  scale_size(range=c(-1,4)) +
  geom_text_repel(data=subset(KaKs.data, FBgn_ID %in% SFP_elements$`D.ame,D.lum,D.nov,D.vir` & `Ka/Ks` > 1 & `Ka/Ks` < 80), aes(label = gene_name), size =3, force = 30, colour = "black") + 
  scale_x_continuous(breaks=seq(from=0, to=30000000, by = 5000000), labels=as.character(seq(0,30,5))) + 
  xlab ("Chromosome coordinates (Mb)") +
  ylab("Ka/Ks") + 
  geom_rect(data=C2inv.qtl, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill="red"), alpha=0.1, inherit.aes = FALSE) + 
  geom_rect(data=C5.qtl, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="red", alpha=0.1, inherit.aes = FALSE) + 
  scale_fill_manual(name = '', values = "red",labels = "PMPZ\nQTL region") +
  theme_bw()
dev.off()

###############################################################################
###############################################################################


###### Set up mel.Encode TPM summary
mel.FBgn_ID_to_GeneSymbol= read.csv("Other.Drosophilas/Dmel/mel.FBgn_ID-to-GeneSymbol.txt", header = T, sep = "\t")
melRPKM.tmp=melRPKM
m.melRPKM.tmp = as.data.frame(melt(melRPKM.tmp))
m.melRPKM.tmp = cSplit(as.data.frame(m.melRPKM.tmp), "variable", "_")
m.melRPKM.tmp=data.frame(m.melRPKM.tmp$mel_FBgn_ID, m.melRPKM.tmp$variable_1, m.melRPKM.tmp$variable_2, m.melRPKM.tmp$variable_3, m.melRPKM.tmp$variable_4, m.melRPKM.tmp$variable_5, m.melRPKM.tmp$value)
colnames(m.melRPKM.tmp) = c("mel_FBgn_ID", "stage", "sex", "status", "age", "tissue", "RPKM")
melRPKM.data = merge(m.melRPKM.tmp, mel.FBgn_ID_to_GeneSymbol, all=TRUE)
melRPKM.data$age =factor(melRPKM.data$age, levels=c("1days","4days","5days","20days","30days"))


###############################################################
##################### GO analysis #############################

## GO analysis

# General set up
gene_lengths = read.table("GO.analysis/FBgn_lengths.txt", header=T, row.names=1)
gene_lengths = as.matrix(gene_lengths[,1,drop=F])
GO_info = read.table("GO.analysis/Trinotate_report_dvir1.06_gene_ontology.txt", header=F, row.names=1,stringsAsFactors=F)
GO_info_listed = apply(GO_info, 1, function(x) unlist(strsplit(x,',')))
names(GO_info_listed) = rownames(GO_info)
features_with_GO = rownames(GO_info)
lengths_features_with_GO = gene_lengths[features_with_GO,]
get_GO_term_descr =  function(x) {
  d = 'none';
  go_info = GOTERM[[x]];
  if (length(go_info) >0) { d = paste(Ontology(go_info), Term(go_info), sep=' ');}
  return(d);
}

# create gene lists and factor labeling
SFP_factors = as.data.frame(SFP_elements$`D.ame,D.lum,D.nov,D.vir`)
SFP_factors$V1 = "SFP-biased"
rownames(SFP_factors) = SFP_elements$`D.ame,D.lum,D.nov,D.vir`
SFP_factors = subset(SFP_factors, select = "V1")

AG_factors = as.data.frame(AG_elements$`D.ame,D.lum,D.nov,D.vir`)
AG_factors$V1 = "AG-biased"
rownames(AG_factors) = AG_elements$`D.ame,D.lum,D.nov,D.vir`
AG_factors = subset(AG_factors, select = "V1")

EB_factors = as.data.frame(EB_elements$`D.ame,D.lum,D.nov,D.vir`)
EB_factors$V1 = "EB-biased"
rownames(EB_factors) = EB_elements$`D.ame,D.lum,D.nov,D.vir`
EB_factors = subset(EB_factors, select = "V1")

TS_factors = as.data.frame(TS_elements$`D.ame,D.lum,D.nov,D.vir`)
TS_factors$V1 = "TS-biased"
rownames(TS_factors) = TS_elements$`D.ame,D.lum,D.nov,D.vir`
TS_factors = subset(TS_factors, select = "V1")

factor_labeling = rbind(SFP_factors, AG_factors, EB_factors, TS_factors)
colnames(factor_labeling) = c('type')
factor_list = unique(factor_labeling[,1])

# build pwf based on ALL DE features
cat_genes_vec = as.integer(features_with_GO %in% rownames(factor_labeling))
pwf=nullp(cat_genes_vec,bias.data=lengths_features_with_GO)
rownames(pwf) = names(GO_info_listed)

GO_enriched_list = list()

for (feature_cat in factor_list) {
  message('Processing category: ', feature_cat)
  cat_genes_vec = as.integer(features_with_GO %in% rownames(factor_labeling)[factor_labeling$type == feature_cat])
  pwf$DEgenes = cat_genes_vec
  res = goseq(pwf,gene2cat=GO_info_listed)
  ## over-represented categories:
  pvals = res$over_represented_pvalue
  pvals[pvals > 1 -1e-10] = 1-1e-10
  q = qvalue(pvals)
  res$over_represented_FDR = q$qvalues
  go_enrich_factor = feature_cat
  enrich_result_table = res[res$over_represented_pvalue<=0.05,]
  descr = unlist(lapply(enrich_result_table$category, get_GO_term_descr))
  enrich_result_table$go_term = descr
  enrich_result_table$factor = go_enrich_factor
  GO_enriched_list[[feature_cat]] = enrich_result_table
}

GO_enrichment_data = rbindlist(GO_enriched_list)

## plot it

pdf("ManuscripPlots/Figure.NA.SFP_GO.pdf", width = 7.92, height = 5.92)
ggplot(subset(GO_enrichment_data, over_represented_FDR < 0.05 & factor == "SFP-biased"), 
       aes(category, -log10(over_represented_pvalue), size = numDEInCat, colour = ontology)) + 
  geom_point()  + 
  xlab(NULL) + 
  geom_text_repel(data = subset(GO_enrichment_data, over_represented_FDR < 0.05 & factor == "SFP-biased" & category != "GO:0016702" & category != "GO:0016706"), 
                  aes(category, -log10(over_represented_pvalue),label=term), 
                  force = 8, 
                  inherit.aes = F, 
                  box.padding = unit(0.35, "lines"), 
                  point.padding = unit(0.5, "lines"), 
                  fontface = "bold", 
                  size = 3) + 
  theme(axis.text.x = element_text(angle = 45, face = "bold", vjust = 1, hjust = 1)) + 
  scale_size(range = c(5,12)) + 
  scale_colour_manual(values=c("#c27d92", "#a8a34b", "#8e61bf")) + 
  scale_y_continuous(limits=c(3, 10))
dev.off()

pdf("ManuscripPlots/Figure_AG_GO.pdf", width = 5.8, height = 4.6)
(gg_AG_GO = ggplot(subset(GO_enrichment_data, over_represented_FDR < 0.05 & factor == "AG-biased"), 
       aes(category, -log10(over_represented_pvalue), 
           size = numDEInCat, colour = ontology)) + 
  geom_point()  + 
  xlab(NULL) + 
  geom_text_repel(aes(category, -log10(over_represented_pvalue),label=term), force = 100, inherit.aes = F) + 
  scale_colour_manual(values=c("#c27d92", "#a8a34b", "#8e61bf")) +
  ggtitle("AG-biased genes") + 
  theme(axis.text.x = element_text(angle = 45, face = "bold", vjust = 0.5)))
dev.off()

pdf("ManuscripPlots/Figure_EB_GO.pdf", width = 9.7, height = 4.6)
(gg_EB_GO = ggplot(subset(GO_enrichment_data, over_represented_FDR < 0.05 & factor == "EB-biased"), 
       aes(category, -log10(over_represented_pvalue), size = numDEInCat, colour = ontology)) + 
  geom_point()  + 
  xlab(NULL) + 
  geom_text_repel(data=subset(GO_enrichment_data, over_represented_FDR < 0.05 & factor == "EB-biased" & numDEInCat > 10),
                  aes(category, -log10(over_represented_pvalue),label=term), force = 10, inherit.aes = F) + 
  scale_colour_manual(values=c("#c27d92", "#a8a34b", "#8e61bf")) +
  ggtitle("EB-biased genes") + 
  theme(axis.text.x = element_text(angle = 45, face = "bold", vjust = 1, hjust = 1, size =7.5)))
dev.off()

# pdf("ManuscripPlots/Figure_AG_and_EB_GO.pdf", width = 15.5, height = 4.6)
topRow = plot_grid(gg_AG_GO, gg_EB_GO, labels = c('A', 'B'), align = 'h', rel_widths = c(1, 1.73))
# dev.off()

# pdf("ManuscripPlots/Figure_TS_GO.pdf", width = 15.5, height = 4.6)
(gg_TS_GO = ggplot(subset(GO_enrichment_data, over_represented_FDR < 0.001 & factor == "TS-biased"), 
       aes(category, -log10(over_represented_pvalue), size = numDEInCat, colour = ontology)) + 
  geom_point()  + 
  xlab(NULL) + 
  geom_text_repel(data = subset(GO_enrichment_data, over_represented_FDR < 0.001 & factor == "TS-biased" & numDEInCat > 100), aes(category, -log10(over_represented_pvalue),label=term), force = 20, inherit.aes = F) +
  scale_colour_manual(values=c("#c27d92", "#a8a34b", "#8e61bf")) +
  scale_size(range=c(1,10)) +
  ggtitle("Testes-biased genes") + 
  theme(axis.text.x = element_text(angle = 45, face = "bold", vjust = 1, hjust = 1, size =7.5)))

pdf("ManuscripPlots/Figure_ALL_biased_GO.pdf", width = 15.5, height = 9.4)
plot_grid(topRow, gg_TS_GO, labels = c('', 'C'), ncol = 1)
dev.off()



###############################
###############################
## output gene info

## Extract gene info for each tissue-biased element
AG.geneinfo = sapply(AG_elements$`D.ame,D.lum,D.nov,D.vir`, geneLookupG, complete =T)
SFP.geneinfo = sapply(SFP_elements$`D.ame,D.lum,D.nov,D.vir`, geneLookupG, complete =T)
EB.geneinfo = sapply(EB_elements$`D.ame,D.lum,D.nov,D.vir`, geneLookupG, complete =T)
TS.geneinfo = sapply(TS_elements$`D.ame,D.lum,D.nov,D.vir`, geneLookupG, complete =T)

## Output individual gene info files
for (i in names(AG.geneinfo)){
  write.table(AG.geneinfo[[i]], paste("GeneInfo/AG-biased/AG.", i, ".txt", sep=""), quote = F, col.names = F, sep = "\t")
}
for (i in names(SFP.geneinfo)){
  write.table(SFP.geneinfo[[i]], paste("GeneInfo/SFPS/SFP.", i, ".txt", sep=""), quote = F, col.names = F, sep = "\t")
}
for (i in names(EB.geneinfo)){
  write.table(EB.geneinfo[[i]], paste("GeneInfo/EB-biased/EB.", i, ".txt", sep=""), quote = F, col.names = F, sep = "\t")
}
for (i in names(TS.geneinfo)){
  write.table(TS.geneinfo[[i]], paste("GeneInfo/TS-biased/TS.", i, ".txt", sep=""), quote = F, col.names = F, sep = "\t")
}
