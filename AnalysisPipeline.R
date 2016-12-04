#################################################################################################################################################
#################################################################################################################################################
#
# __     ___      _ _ _     __  __       _      ____  _   _    _                   
# \ \   / (_)_ __(_) (_)___|  \/  | __ _| | ___|  _ \| \ | |  / \   ___  ___  __ _ 
#  \ \ / /| | '__| | | / __| |\/| |/ _` | |/ _ \ |_) |  \| | / _ \ / __|/ _ \/ _` |           
#   \ V / | | |  | | | \__ \ |  | | (_| | |  __/  _ <| |\  |/ ___ \\__ \  __/ (_| |           
#    \_/  |_|_|  |_|_|_|___/_|  |_|\__,_|_|\___|_| \_\_| \_/_/   \_\___/\___|\__, |
#                                                                               |_|               Yasir H. Ahmed-Braimah
# 
#################################################################################################################################################
#################################################################################################################################################


## Load annotations information
grpTrinotate = read.csv(file = "Annotations/Trinotate_report_dvir1.06_subset.txt", sep = "\t", header = T, na.strings = ".", stringsAsFactors=FALSE)
gffRecord = read.table(file = "Annotations/FBgn_ID_name_coordinates.txt", header = T)
melOrths = read.table(file = "Annotations/mel_orths.txt", header = T)
melSFPs = read.table(file = "Annotations/ACPlist.Findlay.20130301.txt", header = T, sep = "\t")
melRPKM = read.table(file = "Annotations/mel.modEncode.RPKM.matrix", header = T, sep = "\t")
melOrthsAll = aggregate(mel_GeneSymbol~FBgn_ID, data = melOrths, toString)
tmp.merged = merge(melOrthsAll, grpTrinotate, all=TRUE)
Annots = merge(tmp.merged, gffRecord, all=TRUE)
rm(grpTrinotate, melOrthsAll, tmp.merged)
TrinOrths = read.table("Annotations/Trin3.orthology.txt", header = T, sep = "\t")
amrTrinotate = read.csv("Annotations/Trinotate_report_amrTrin3_subset.txt", header = T, sep = "\t", na.strings = ".", stringsAsFactors=FALSE)
lumTrinotate = read.csv("Annotations/Trinotate_report_lumTrin3_subset.txt", header = T, sep = "\t", na.strings = ".")
novTrinotate = read.csv("Annotations/Trinotate_report_novTrin3_subset.txt", header = T, sep = "\t", na.strings = ".")
virTrinotate = read.csv("Annotations/Trinotate_report_virTrin3_subset.txt", header = T, sep = "\t", na.strings = ".")
## Read in sample info:
grpSamples_data = read.table("Annotations/samples.txt", header=F, check.names=F, fill=T)
grpSamples_data = grpSamples_data[grpSamples_data[,2] != '',]
amrSamples_data = subset(grpSamples_data, grepl("amr", grpSamples_data$V1))
lumSamples_data = subset(grpSamples_data, grepl("lum", grpSamples_data$V1))
novSamples_data = subset(grpSamples_data, grepl("nov", grpSamples_data$V1))
virSamples_data = subset(grpSamples_data, grepl("vir", grpSamples_data$V1))


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


# ## Barplot of gene counts by sample
# libSizes = as.data.frame(colSums(grpCountsMatrix))
# libSizes = cbind(sample = row.names(libSizes), libSizes)
# row.names(libSizes)= NULL
# colnames(libSizes) = c("sample", "Total_reads")
# ggplot(libSizes, aes(sample, Total_reads)) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = -90, hjust = 0)) + geom_hline(yintercept = 20000000)
# 
# ## Boxplot of log10(TPM) across all samples
# m.expData=melt(as.matrix(grpTmmMatrix))
# colnames(m.expData) = c("gene_id", "replicate", "TPM")
# m.expData.exp= within(m.expData, replicate=data.frame(do.call('rbind', strsplit(as.character(replicate),'_',fixed=TRUE))))
# m.expData=data.frame(m.expData, m.expData.exp$replicate$X1, m.expData.exp$replicate$X2, m.expData.exp$replicate$X3)
# colnames(m.expData) = c("gene_id", "replicate", "TPM", "species", "tissue", "rep_num")
# m.expData$TPM = m.expData$TPM + 1
# p = ggplot(m.expData)
# p = p + geom_boxplot(aes(x = replicate, y = log10(TPM), fill = tissue, colour = species), size = 0.3, alpha = I(1/3))
# p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
# p = p + scale_fill_hue(l = 50, h.start = 200)
# p
# 
# ## Estimate of the number of expressed genes (Brian Haas' method)
# # extract the data between 10 TPM and 100 TPM
# amr_filt_data = amrTpmMatrix.notCrossNorm[amrTpmMatrix.notCrossNorm[,1] > -100 & amrTpmMatrix.notCrossNorm[,1] < -10,]
# lum_filt_data = lumTpmMatrix.notCrossNorm[lumTpmMatrix.notCrossNorm[,1] > -100 & lumTpmMatrix.notCrossNorm[,1] < -10,]
# nov_filt_data = novTpmMatrix.notCrossNorm[novTpmMatrix.notCrossNorm[,1] > -100 & novTpmMatrix.notCrossNorm[,1] < -10,]
# vir_filt_data = virTpmMatrix.notCrossNorm[virTpmMatrix.notCrossNorm[,1] > -100 & virTpmMatrix.notCrossNorm[,1] < -10,]
# grp_filt_data = grpTpmMatrix.notCrossNorm[grpTpmMatrix.notCrossNorm[,1] > -100 & grpTpmMatrix.notCrossNorm[,1] < -10,]
# # perform a linear regression on this filtered subset of the data
# amr_fit = lm(amr_filt_data[,2] ~ amr_filt_data[,1])
# lum_fit = lm(lum_filt_data[,2] ~ lum_filt_data[,1])
# nov_fit = lm(nov_filt_data[,2] ~ nov_filt_data[,1])
# vir_fit = lm(vir_filt_data[,2] ~ vir_filt_data[,1])
# grp_fit = lm(grp_filt_data[,2] ~ grp_filt_data[,1])
# print(amr_fit)
# print(lum_fit)
# print(nov_fit)
# print(vir_fit)
# print(grp_fit)
# # plot it
# amr_fit_plot=ggplot(amrTpmMatrix.notCrossNorm, aes(neg_min_tpm,num_features)) + geom_point() + scale_x_continuous(limits=c(-100,0)) + scale_y_continuous(limits=c(0,20000)) + geom_smooth(data=amr_filt_data, method = "lm") + geom_hline(yintercept = 13358, colour = "green") + ggtitle("amrTrinity")
# lum_fit_plot=ggplot(lumTpmMatrix.notCrossNorm, aes(neg_min_tpm,num_features)) + geom_point() + scale_x_continuous(limits=c(-100,0)) + scale_y_continuous(limits=c(0,20000)) + geom_smooth(data=lum_filt_data, method = "lm") + geom_hline(yintercept = 14805, colour = "green") + ggtitle("lumTrinity")
# nov_fit_plot=ggplot(novTpmMatrix.notCrossNorm, aes(neg_min_tpm,num_features)) + geom_point() + scale_x_continuous(limits=c(-100,0)) + scale_y_continuous(limits=c(0,20000)) + geom_smooth(data=nov_filt_data, method = "lm") + geom_hline(yintercept = 13646, colour = "green") + ggtitle("novTrinity")
# vir_fit_plot=ggplot(virTpmMatrix.notCrossNorm, aes(neg_min_tpm,num_features)) + geom_point() + scale_x_continuous(limits=c(-100,0)) + scale_y_continuous(limits=c(0,20000)) + geom_smooth(data=vir_filt_data, method = "lm") + geom_hline(yintercept = 14616, colour = "green") + ggtitle("virTrinity")
# grp_fit_plot=ggplot(grpTpmMatrix.notCrossNorm, aes(neg_min_tpm,num_features)) + geom_point() + scale_x_continuous(limits=c(-100,0)) + scale_y_continuous(limits=c(0,20000)) + geom_smooth(data=grp_filt_data, method = "lm") + geom_hline(yintercept = 9474, colour = "green") + ggtitle("dvir_1.06")
# plot_grid(amr_fit_plot, lum_fit_plot, nov_fit_plot, vir_fit_plot, grp_fit_plot,nrow = 2)
# 
# 
# ## calculate dispersion
# d = DGEList(counts = grpCountsMatrix, group = grpSamples_data$V1)
# d = calcNormFactors(d)
# d = estimateCommonDisp(d)
# d = estimateTagwiseDisp(d)
# summary(d$tagwise.dispersion)
# ## Plot biological coefficient of variation
# plotBCV(d)
# ## Plot grouping of samples
# plotMDS(d, method = "bcv", col=as.numeric(d$samples$group))
# 
# ## Plot sample correlation
# data = log2(grpCountsMatrix+1)
# data = as.matrix(data)
# sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
# sample_dist = as.dist(1-sample_cor)
# hc_samples = hclust(sample_dist, method='complete')
# heatmap.3(sample_cor, dendrogram='both', Rowv=as.dendrogram(hc_samples), Colv=as.dendrogram(hc_samples), col = greenred(75), scale='none', symm=TRUE, key=TRUE,density.info='none', trace='none', symkey=FALSE, symbreaks=F, margins=c(10,10), cexCol=1, cexRow=1, cex.main=0.75, main=paste("sample correlation matrix"))

## Filter count data by minimum count across ANY sample
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

## normalize by DESeq method:
# meta = data.frame(row.names=colnames(grpCountsMatrix), condition=grpSamples_data$V1)
# grpCountData=round(grpCountsMatrix)
# grpCountData_normByDESeq = newCountDataSet(grpCountData, meta)
# grpCountData_normByDESeq = estimateSizeFactors(grpCountData_normByDESeq)
# grpCountData_normByDESeq = data.frame(counts(grpCountData_normByDESeq, normalized=T))

# Example MA plot between replicates
# MA_BPlot(amrCountsMatrix.min400count, "amr_AG_1", "amr_AG_2")

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

#############################################################################
################# Filter by "good" replicates ###############################
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

##############################
## Create specificity matrices (dvir1.06)
amr.dvir1.06.MeanTPMmatrix = subset(grpMeanTPMmatrix.BRR, select=c("FBgn_ID", "amr_AG", "amr_CR", "amr_EB", "amr_TS"))
tmp.amr.dvir1.06.MeanTPMmatrix = amr.dvir1.06.MeanTPMmatrix
rownames(tmp.amr.dvir1.06.MeanTPMmatrix) = tmp.amr.dvir1.06.MeanTPMmatrix[,1]
tmp.amr.dvir1.06.MeanTPMmatrix[,1] = NULL
amr.dvir1.06_Specificity_table = YazSpecificity(tmp.amr.dvir1.06.MeanTPMmatrix)
amr.dvir1.06_Specificity_table = as.data.frame(amr.dvir1.06_Specificity_table)
rm(tmp.amr.dvir1.06.MeanTPMmatrix)

lum.dvir1.06.MeanTPMmatrix = subset(grpMeanTPMmatrix.BRR, select=c("FBgn_ID", "lum_AG", "lum_CR", "lum_EB", "lum_TS"))
tmp.lum.dvir1.06.MeanTPMmatrix = lum.dvir1.06.MeanTPMmatrix
rownames(tmp.lum.dvir1.06.MeanTPMmatrix) = tmp.lum.dvir1.06.MeanTPMmatrix[,1]
tmp.lum.dvir1.06.MeanTPMmatrix[,1] = NULL
lum.dvir1.06_Specificity_table = YazSpecificity(tmp.lum.dvir1.06.MeanTPMmatrix)
lum.dvir1.06_Specificity_table = as.data.frame(lum.dvir1.06_Specificity_table)
rm(tmp.lum.dvir1.06.MeanTPMmatrix)

nov.dvir1.06.MeanTPMmatrix = subset(grpMeanTPMmatrix.BRR, select=c("FBgn_ID", "nov_AG", "nov_CR", "nov_EB", "nov_TS"))
tmp.nov.dvir1.06.MeanTPMmatrix = nov.dvir1.06.MeanTPMmatrix
rownames(tmp.nov.dvir1.06.MeanTPMmatrix) = tmp.nov.dvir1.06.MeanTPMmatrix[,1]
tmp.nov.dvir1.06.MeanTPMmatrix[,1] = NULL
nov.dvir1.06_Specificity_table = YazSpecificity(tmp.nov.dvir1.06.MeanTPMmatrix)
nov.dvir1.06_Specificity_table = as.data.frame(nov.dvir1.06_Specificity_table)
rm(tmp.nov.dvir1.06.MeanTPMmatrix)

vir.dvir1.06.MeanTPMmatrix = subset(grpMeanTPMmatrix.BRR, select=c("FBgn_ID", "vir_AG", "vir_CR", "vir_EB", "vir_TS"))
tmp.vir.dvir1.06.MeanTPMmatrix = vir.dvir1.06.MeanTPMmatrix
rownames(tmp.vir.dvir1.06.MeanTPMmatrix) = tmp.vir.dvir1.06.MeanTPMmatrix[,1]
tmp.vir.dvir1.06.MeanTPMmatrix[,1] = NULL
vir.dvir1.06_Specificity_table = YazSpecificity(tmp.vir.dvir1.06.MeanTPMmatrix)
vir.dvir1.06_Specificity_table = as.data.frame(vir.dvir1.06_Specificity_table)
rm(tmp.vir.dvir1.06.MeanTPMmatrix)

tmp.grpMeanTPMmatrix = grpMeanTPMmatrix.BRR
rownames(tmp.grpMeanTPMmatrix) = tmp.grpMeanTPMmatrix[,1]
tmp.grpMeanTPMmatrix[,1] = NULL
dvir1.06_Specificity_table = YazSpecificity(tmp.grpMeanTPMmatrix)
dvir1.06_Specificity_table = as.data.frame(dvir1.06_Specificity_table)
rm(tmp.grpMeanTPMmatrix)

## Create specificity matrices (Trinity)
tmp.amrMeanTPMmatrix.BRR = amrMeanTPMmatrix.BRR
rownames(tmp.amrMeanTPMmatrix.BRR) = tmp.amrMeanTPMmatrix.BRR[,1]
tmp.amrMeanTPMmatrix.BRR[,1] = NULL
amr_Specificity_table = YazSpecificity(tmp.amrMeanTPMmatrix.BRR)
amr_Specificity_table = as.data.frame(amr_Specificity_table)
rm(tmp.amrMeanTPMmatrix.BRR)

tmp.lumMeanTPMmatrix.BRR = lumMeanTPMmatrix.BRR
rownames(tmp.lumMeanTPMmatrix.BRR) = tmp.lumMeanTPMmatrix.BRR[,1]
tmp.lumMeanTPMmatrix.BRR[,1] = NULL
lum_Specificity_table = YazSpecificity(tmp.lumMeanTPMmatrix.BRR)
lum_Specificity_table = as.data.frame(lum_Specificity_table)
rm(tmp.lumMeanTPMmatrix.BRR)

tmp.novMeanTPMmatrix.BRR = novMeanTPMmatrix.BRR
rownames(tmp.novMeanTPMmatrix.BRR) = tmp.novMeanTPMmatrix.BRR[,1]
tmp.novMeanTPMmatrix.BRR[,1] = NULL
nov_Specificity_table = YazSpecificity(tmp.novMeanTPMmatrix.BRR)
nov_Specificity_table = as.data.frame(nov_Specificity_table)
rm(tmp.novMeanTPMmatrix.BRR)

tmp.virMeanTPMmatrix.BRR = virMeanTPMmatrix.BRR
rownames(tmp.virMeanTPMmatrix.BRR) = tmp.virMeanTPMmatrix.BRR[,1]
tmp.virMeanTPMmatrix.BRR[,1] = NULL
vir_Specificity_table = YazSpecificity(tmp.virMeanTPMmatrix.BRR)
vir_Specificity_table = as.data.frame(vir_Specificity_table)
rm(tmp.virMeanTPMmatrix.BRR)

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
amr.SFP.list=unique(subset(amrTrinotate, gene_id %in% amr.AG.list & SignalP != "NA")$gene_id)
amr.SFP.list = droplevels(amr.SFP.list)
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

### Identify tissue-biased genes within species (>4-fold, < 0.001 FDR)

## D. americana
# AG-biased genes
amr.dvir1.06.AG.Contrasts=makeContrasts(AG.v.CR=amr_AG-amr_CR, AG.v.EB=amr_AG-amr_EB, AG.v.TS=amr_AG-amr_TS, levels = grp.design)
amr.dvir1.06.lrt.AG.v.rest = glmLRT(grp.fit, contrast = amr.dvir1.06.AG.Contrasts)
amr.dvir1.06.lrt.AG.v.rest.tTags = topTags(amr.dvir1.06.lrt.AG.v.rest, n = NULL)
amr.dvir1.06.lrt.AG.v.rest.tTags.table = amr.dvir1.06.lrt.AG.v.rest.tTags$table
amr.dvir1.06.AG.list=rownames(subset(amr.dvir1.06.lrt.AG.v.rest.tTags.table, logFC.AG.v.CR > 2 & logFC.AG.v.EB > 2 & logFC.AG.v.TS > 2 & FDR<0.001))
amr.dvir1.06.AG.list.byTPM=subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% amr.dvir1.06.AG.list & amr_AG > 10 & amr_CR < 5 & amr_EB < 5 & amr_TS < 5)$FBgn_ID
amr.dvir1.06.SFP.list=unique(subset(Annots, FBgn_ID %in% amr.dvir1.06.AG.list & SignalP == "YES")$FBgn_ID)
# EB-biased genes
amr.dvir1.06.EB.Contrasts=makeContrasts(EB.v.CR=amr_EB-amr_CR, EB.v.AG=amr_EB-amr_AG, EB.v.TS=amr_EB-amr_TS, levels = grp.design)
amr.dvir1.06.lrt.EB.v.rest = glmLRT(grp.fit, contrast = amr.dvir1.06.EB.Contrasts)
amr.dvir1.06.lrt.EB.v.rest.tTags = topTags(amr.dvir1.06.lrt.EB.v.rest, n = NULL)
amr.dvir1.06.lrt.EB.v.rest.tTags.table = amr.dvir1.06.lrt.EB.v.rest.tTags$table
amr.dvir1.06.EB.list=rownames(subset(amr.dvir1.06.lrt.EB.v.rest.tTags.table, logFC.EB.v.CR > 2 & logFC.EB.v.AG > 2 & logFC.EB.v.TS > 2 & FDR<0.001))
amr.dvir1.06.EB.list.byTPM=subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% amr.dvir1.06.EB.list & amr_AG < 5 & amr_CR < 5 & amr_EB > 10 & amr_TS < 5)$FBgn_ID
# TS-biased genes
amr.dvir1.06.TS.Contrasts=makeContrasts(TS.v.CR=amr_TS-amr_CR, TS.v.EB=amr_TS-amr_EB, TS.v.AG=amr_TS-amr_AG, levels = grp.design)
amr.dvir1.06.lrt.TS.v.rest = glmLRT(grp.fit, contrast = amr.dvir1.06.TS.Contrasts)
amr.dvir1.06.lrt.TS.v.rest.tTags = topTags(amr.dvir1.06.lrt.TS.v.rest, n = NULL)
amr.dvir1.06.lrt.TS.v.rest.tTags.table = amr.dvir1.06.lrt.TS.v.rest.tTags$table
amr.dvir1.06.TS.list=rownames(subset(amr.dvir1.06.lrt.TS.v.rest.tTags.table, logFC.TS.v.CR > 2 & logFC.TS.v.EB > 2 & logFC.TS.v.AG > 2 & FDR<0.001))
amr.dvir1.06.TS.list.byTPM=subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% amr.dvir1.06.TS.list & amr_AG < 5 & amr_CR < 5 & amr_EB < 5 & amr_TS > 10)$FBgn_ID

## D. lummei
# AG-biased genes
lum.dvir1.06.AG.Contrasts=makeContrasts(AG.v.CR=lum_AG-lum_CR, AG.v.EB=lum_AG-lum_EB, AG.v.TS=lum_AG-lum_TS, levels = grp.design)
lum.dvir1.06.lrt.AG.v.rest = glmLRT(grp.fit, contrast = lum.dvir1.06.AG.Contrasts)
lum.dvir1.06.lrt.AG.v.rest.tTags = topTags(lum.dvir1.06.lrt.AG.v.rest, n = NULL)
lum.dvir1.06.lrt.AG.v.rest.tTags.table = lum.dvir1.06.lrt.AG.v.rest.tTags$table
lum.dvir1.06.AG.list=rownames(subset(lum.dvir1.06.lrt.AG.v.rest.tTags.table, logFC.AG.v.CR > 2 & logFC.AG.v.EB > 2 & logFC.AG.v.TS > 2 & FDR<0.001))
lum.dvir1.06.AG.list.byTPM=subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% lum.dvir1.06.AG.list & lum_AG > 10 & lum_CR < 5 & lum_EB < 5 & lum_TS < 5)$FBgn_ID
lum.dvir1.06.SFP.list=unique(subset(Annots, FBgn_ID %in% lum.dvir1.06.AG.list & SignalP == "YES")$FBgn_ID)
# EB-biased genes
lum.dvir1.06.EB.Contrasts=makeContrasts(EB.v.CR=lum_EB-lum_CR, EB.v.AG=lum_EB-lum_AG, EB.v.TS=lum_EB-lum_TS, levels = grp.design)
lum.dvir1.06.lrt.EB.v.rest = glmLRT(grp.fit, contrast = lum.dvir1.06.EB.Contrasts)
lum.dvir1.06.lrt.EB.v.rest.tTags = topTags(lum.dvir1.06.lrt.EB.v.rest, n = NULL)
lum.dvir1.06.lrt.EB.v.rest.tTags.table = lum.dvir1.06.lrt.EB.v.rest.tTags$table
lum.dvir1.06.EB.list=rownames(subset(lum.dvir1.06.lrt.EB.v.rest.tTags.table, logFC.EB.v.CR > 2 & logFC.EB.v.AG > 2 & logFC.EB.v.TS > 2 & FDR<0.001))
lum.dvir1.06.EB.list.byTPM=subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% lum.dvir1.06.EB.list & lum_AG < 5 & lum_CR < 5 & lum_EB > 10 & lum_TS < 5)$FBgn_ID
# TS-biased genes
lum.dvir1.06.TS.Contrasts=makeContrasts(TS.v.CR=lum_TS-lum_CR, TS.v.EB=lum_TS-lum_EB, TS.v.AG=lum_TS-lum_AG, levels = grp.design)
lum.dvir1.06.lrt.TS.v.rest = glmLRT(grp.fit, contrast = lum.dvir1.06.TS.Contrasts)
lum.dvir1.06.lrt.TS.v.rest.tTags = topTags(lum.dvir1.06.lrt.TS.v.rest, n = NULL)
lum.dvir1.06.lrt.TS.v.rest.tTags.table = lum.dvir1.06.lrt.TS.v.rest.tTags$table
lum.dvir1.06.TS.list=rownames(subset(lum.dvir1.06.lrt.TS.v.rest.tTags.table, logFC.TS.v.CR > 2 & logFC.TS.v.EB > 2 & logFC.TS.v.AG > 2 & FDR<0.001))
lum.dvir1.06.TS.list.byTPM=subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% lum.dvir1.06.TS.list & lum_AG < 5 & lum_CR < 5 & lum_EB < 5 & lum_TS > 10)$FBgn_ID

## D. novamexicana
# AG-biased genes
nov.dvir1.06.AG.Contrasts=makeContrasts(AG.v.CR=nov_AG-nov_CR, AG.v.EB=nov_AG-nov_EB, AG.v.TS=nov_AG-nov_TS, levels = grp.design)
nov.dvir1.06.lrt.AG.v.rest = glmLRT(grp.fit, contrast = nov.dvir1.06.AG.Contrasts)
nov.dvir1.06.lrt.AG.v.rest.tTags = topTags(nov.dvir1.06.lrt.AG.v.rest, n = NULL)
nov.dvir1.06.lrt.AG.v.rest.tTags.table = nov.dvir1.06.lrt.AG.v.rest.tTags$table
nov.dvir1.06.AG.list=rownames(subset(nov.dvir1.06.lrt.AG.v.rest.tTags.table, logFC.AG.v.CR > 2 & logFC.AG.v.EB > 2 & logFC.AG.v.TS > 2 & FDR<0.001))
nov.dvir1.06.AG.list.byTPM=subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% nov.dvir1.06.AG.list & nov_AG > 10 & nov_CR < 5 & nov_EB < 5 & nov_TS < 5)$FBgn_ID
nov.dvir1.06.SFP.list=unique(subset(Annots, FBgn_ID %in% nov.dvir1.06.AG.list & SignalP == "YES")$FBgn_ID)
# EB-biased genes
nov.dvir1.06.EB.Contrasts=makeContrasts(EB.v.CR=nov_EB-nov_CR, EB.v.AG=nov_EB-nov_AG, EB.v.TS=nov_EB-nov_TS, levels = grp.design)
nov.dvir1.06.lrt.EB.v.rest = glmLRT(grp.fit, contrast = nov.dvir1.06.EB.Contrasts)
nov.dvir1.06.lrt.EB.v.rest.tTags = topTags(nov.dvir1.06.lrt.EB.v.rest, n = NULL)
nov.dvir1.06.lrt.EB.v.rest.tTags.table = nov.dvir1.06.lrt.EB.v.rest.tTags$table
nov.dvir1.06.EB.list=rownames(subset(nov.dvir1.06.lrt.EB.v.rest.tTags.table, logFC.EB.v.CR > 2 & logFC.EB.v.AG > 2 & logFC.EB.v.TS > 2 & FDR<0.001))
nov.dvir1.06.EB.list.byTPM=subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% nov.dvir1.06.EB.list & nov_AG < 5 & nov_CR < 5 & nov_EB > 10 & nov_TS < 5)$FBgn_ID
# TS-biased genes
nov.dvir1.06.TS.Contrasts=makeContrasts(TS.v.CR=nov_TS-nov_CR, TS.v.EB=nov_TS-nov_EB, TS.v.AG=nov_TS-nov_AG, levels = grp.design)
nov.dvir1.06.lrt.TS.v.rest = glmLRT(grp.fit, contrast = nov.dvir1.06.TS.Contrasts)
nov.dvir1.06.lrt.TS.v.rest.tTags = topTags(nov.dvir1.06.lrt.TS.v.rest, n = NULL)
nov.dvir1.06.lrt.TS.v.rest.tTags.table = nov.dvir1.06.lrt.TS.v.rest.tTags$table
nov.dvir1.06.TS.list=rownames(subset(nov.dvir1.06.lrt.TS.v.rest.tTags.table, logFC.TS.v.CR > 2 & logFC.TS.v.EB > 2 & logFC.TS.v.AG > 2 & FDR<0.001))
nov.dvir1.06.TS.list.byTPM=subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% nov.dvir1.06.TS.list & nov_AG < 5 & nov_CR < 5 & nov_EB < 5 & nov_TS > 10)$FBgn_ID

## D. virilis
# AG-biased genes
vir.dvir1.06.AG.Contrasts=makeContrasts(AG.v.CR=vir_AG-vir_CR, AG.v.EB=vir_AG-vir_EB, AG.v.TS=vir_AG-vir_TS, levels = grp.design)
vir.dvir1.06.lrt.AG.v.rest = glmLRT(grp.fit, contrast = vir.dvir1.06.AG.Contrasts)
vir.dvir1.06.lrt.AG.v.rest.tTags = topTags(vir.dvir1.06.lrt.AG.v.rest, n = NULL)
vir.dvir1.06.lrt.AG.v.rest.tTags.table = vir.dvir1.06.lrt.AG.v.rest.tTags$table
vir.dvir1.06.AG.list=rownames(subset(vir.dvir1.06.lrt.AG.v.rest.tTags.table, logFC.AG.v.CR > 2 & logFC.AG.v.EB > 2 & logFC.AG.v.TS > 2 & FDR<0.001))
vir.dvir1.06.AG.list.byTPM=subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% vir.dvir1.06.AG.list & vir_AG > 10 & vir_CR < 5 & vir_EB < 5 & vir_TS < 5)$FBgn_ID
vir.dvir1.06.SFP.list=unique(subset(Annots, FBgn_ID %in% vir.dvir1.06.AG.list & SignalP == "YES")$FBgn_ID)
# EB-biased genes
vir.dvir1.06.EB.Contrasts=makeContrasts(EB.v.CR=vir_EB-vir_CR, EB.v.AG=vir_EB-vir_AG, EB.v.TS=vir_EB-vir_TS, levels = grp.design)
vir.dvir1.06.lrt.EB.v.rest = glmLRT(grp.fit, contrast = vir.dvir1.06.EB.Contrasts)
vir.dvir1.06.lrt.EB.v.rest.tTags = topTags(vir.dvir1.06.lrt.EB.v.rest, n = NULL)
vir.dvir1.06.lrt.EB.v.rest.tTags.table = vir.dvir1.06.lrt.EB.v.rest.tTags$table
vir.dvir1.06.EB.list=rownames(subset(vir.dvir1.06.lrt.EB.v.rest.tTags.table, logFC.EB.v.CR > 2 & logFC.EB.v.AG > 2 & logFC.EB.v.TS > 2 & FDR<0.001))
vir.dvir1.06.EB.list.byTPM=subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% vir.dvir1.06.EB.list & vir_AG < 5 & vir_CR < 5 & vir_EB > 10 & vir_TS < 5)$FBgn_ID
# TS-biased genes
vir.dvir1.06.TS.Contrasts=makeContrasts(TS.v.CR=vir_TS-vir_CR, TS.v.EB=vir_TS-vir_EB, TS.v.AG=vir_TS-vir_AG, levels = grp.design)
vir.dvir1.06.lrt.TS.v.rest = glmLRT(grp.fit, contrast = vir.dvir1.06.TS.Contrasts)
vir.dvir1.06.lrt.TS.v.rest.tTags = topTags(vir.dvir1.06.lrt.TS.v.rest, n = NULL)
vir.dvir1.06.lrt.TS.v.rest.tTags.table = vir.dvir1.06.lrt.TS.v.rest.tTags$table
vir.dvir1.06.TS.list=rownames(subset(vir.dvir1.06.lrt.TS.v.rest.tTags.table, logFC.TS.v.CR > 2 & logFC.TS.v.EB > 2 & logFC.TS.v.AG > 2 & FDR<0.001))
vir.dvir1.06.TS.list.byTPM=subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% vir.dvir1.06.TS.list & vir_AG < 5 & vir_CR < 5 & vir_EB < 5 & vir_TS > 10)$FBgn_ID

# Create list objects containins tissue-biased candidates by species
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
#sapply(AG_vetted_elements, length)



### VennDiagram: This is how to draw a VennDiagram using several lists of genes
AG_candidates_Vdiag=venn.diagram(AG_candidates, NULL, fill=c("#670066", "#86bd38", "#29b5ff", "#717e36"), alpha=c(0.5,0.5,0.5,0.5), cex = 1.5, cat.fontface= 4, cat.cex = 1.25, resolution = 1000, main = "AG-biased")
SFP_candidates_Vdiag=venn.diagram(SFP_candidates, NULL, fill=c("#670066", "#86bd38", "#29b5ff", "#717e36"), alpha=c(0.5,0.5,0.5,0.5), cex = 1.5, cat.fontface= 4, cat.cex = 1.25, resolution = 1000, main = "SFP-biased")
EB_candidates_Vdiag=venn.diagram(EB_candidates, NULL, fill=c("#670066", "#86bd38", "#29b5ff", "#717e36"), alpha=c(0.5,0.5,0.5,0.5), cex = 1.5, cat.fontface= 4, cat.cex = 1.25, resolution = 1000, main = "EB-biased")
TS_candidates_Vdiag=venn.diagram(TS_candidates, NULL, fill=c("#670066", "#86bd38", "#29b5ff", "#717e36"), alpha=c(0.5,0.5,0.5,0.5), cex = 1.5, cat.fontface= 4, cat.cex = 1.25, resolution = 1000, main = "TS-biased")
grid.arrange(gTree(children=AG_candidates_Vdiag), gTree(children =SFP_candidates_Vdiag), gTree(children=EB_candidates_Vdiag), gTree(children=TS_candidates_Vdiag), ncol=2)

## write out tissue-biased transcripts for GO analysis
# write.table(AG_elements$`D.ame,D.lum,D.nov,D.vir`, file = "GO.analysis/AG-biased_genes.list", quote = F, row.names = F, col.names = F)
# write.table(SFP_elements$`D.ame,D.lum,D.nov,D.vir`, file = "GO.analysis/SFP-biased_genes.list", quote = F, row.names = F, col.names = F)
# write.table(EB_elements$`D.ame,D.lum,D.nov,D.vir`, file = "GO.analysis/EB-biased_genes.list", quote = F, row.names = F, col.names = F)
# write.table(TS_elements$`D.ame,D.lum,D.nov,D.vir`, file = "GO.analysis/TS-biased_genes.list", quote = F, row.names = F, col.names = F)

## refine the candidates:
## Vetted AG candidates
AG_vetted_elements = NULL
AG_vetted_elements$D.ame = as.character(subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% AG_elements$D.ame & lum_AG < 5 & nov_AG < 5 & vir_AG < 5 & amr_AG > 10)$FBgn_ID)
AG_vetted_elements$D.lum = as.character(subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% AG_elements$D.lum & lum_AG > 10 & nov_AG < 5 & vir_AG < 5 & amr_AG < 5 )$FBgn_ID)
AG_vetted_elements$D.nov = as.character(subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% AG_elements$D.nov & lum_AG < 5 & nov_AG > 10 & vir_AG < 5 & amr_AG < 5 )$FBgn_ID)
AG_vetted_elements$D.vir = as.character(subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% AG_elements$D.vir & lum_AG < 5 & nov_AG < 5 & vir_AG > 10 & amr_AG < 5 )$FBgn_ID)
AG_vetted_elements$`D.ame,D.lum` = as.character(subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% AG_elements$`D.ame,D.lum` & lum_AG > 10 & nov_AG < 5 & vir_AG < 5 & amr_AG > 10 )$FBgn_ID)
AG_vetted_elements$`D.ame,D.nov` = as.character(subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% AG_elements$`D.ame,D.nov` & lum_AG < 5 & nov_AG > 10 & vir_AG < 5 & amr_AG > 10 )$FBgn_ID)
AG_vetted_elements$`D.ame,D.vir` = as.character(subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% AG_elements$`D.ame,D.vir` & lum_AG < 5 & nov_AG < 5 & vir_AG > 10 & amr_AG > 10 )$FBgn_ID)
AG_vetted_elements$`D.lum,D.nov` = as.character(subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% AG_elements$`D.lum,D.nov` & lum_AG > 10 & nov_AG > 10 & vir_AG < 5 & amr_AG < 5 )$FBgn_ID)
AG_vetted_elements$`D.lum,D.vir` = as.character(subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% AG_elements$`D.lum,D.vir` & lum_AG > 10 & nov_AG < 5 & vir_AG > 10 & amr_AG < 5 )$FBgn_ID)
AG_vetted_elements$`D.nov,D.vir` = as.character(subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% AG_elements$`D.nov,D.vir` & lum_AG < 5 & nov_AG > 10 & vir_AG > 10 & amr_AG < 5 )$FBgn_ID)
AG_vetted_elements$`D.ame,D,lum,D.nov` = as.character(subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% AG_elements$`D.ame,D,lum,D.nov` & lum_AG > 10 & nov_AG > 10 & vir_AG < 5 & amr_AG > 10 )$FBgn_ID)
AG_vetted_elements$`D.ame,D.lum,D.vir` = as.character(subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% AG_elements$`D.ame,D.lum,D.vir` & lum_AG > 10 & nov_AG < 5 & vir_AG > 10 & amr_AG > 10 )$FBgn_ID)
AG_vetted_elements$`D.ame,D.nov,D.vir` = as.character(subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% AG_elements$`D.ame,D.nov,D.vir` & lum_AG < 5 & nov_AG > 10 & vir_AG > 10 & amr_AG > 10 )$FBgn_ID)
AG_vetted_elements$`D.lum,D.nov,D.vir` = as.character(subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% AG_elements$`D.lum,D.nov,D.vir` & lum_AG > 10 & nov_AG > 10 & vir_AG > 10 & amr_AG < 5 )$FBgn_ID)
AG_vetted_elements$`D.ame,D.lum,D.nov,D.vir` = as.character(AG_elements$`D.ame,D.lum,D.nov,D.vir`)

amr_vetted_AG_candidates = c(AG_vetted_elements$D.ame, AG_vetted_elements$`D.ame,D.lum`, AG_vetted_elements$`D.ame,D.nov`, AG_vetted_elements$`D.ame,D.vir`, AG_vetted_elements$`D.ame,D,lum,D.nov`, AG_vetted_elements$`D.ame,D.lum,D.vir`, AG_vetted_elements$`D.ame,D.nov,D.vir`, AG_vetted_elements$`D.ame,D.lum,D.nov,D.vir`)
lum_vetted_AG_candidates = c(AG_vetted_elements$D.lum, AG_vetted_elements$`D.ame,D.lum`, AG_vetted_elements$`D.lum,D.nov`, AG_vetted_elements$`D.lum,D.vir`, AG_vetted_elements$`D.ame,D,lum,D.nov`, AG_vetted_elements$`D.ame,D.lum,D.vir`, AG_vetted_elements$`D.lum,D.nov,D.vir`, AG_vetted_elements$`D.ame,D.lum,D.nov,D.vir`)
nov_vetted_AG_candidates = c(AG_vetted_elements$D.nov, AG_vetted_elements$`D.ame,D.nov`, AG_vetted_elements$`D.lum,D.nov`, AG_vetted_elements$`D.nov,D.vir`, AG_vetted_elements$`D.ame,D,lum,D.nov`, AG_vetted_elements$`D.ame,D.nov,D.vir`, AG_vetted_elements$`D.lum,D.nov,D.vir`, AG_vetted_elements$`D.ame,D.lum,D.nov,D.vir`)
vir_vetted_AG_candidates = c(AG_vetted_elements$D.vir, AG_vetted_elements$`D.ame,D.vir`, AG_vetted_elements$`D.lum,D.vir`, AG_vetted_elements$`D.nov,D.vir`, AG_vetted_elements$`D.ame,D,lum,D.vir`, AG_vetted_elements$`D.ame,D.nov,D.vir`, AG_vetted_elements$`D.lum,D.nov,D.vir`, AG_vetted_elements$`D.ame,D.lum,D.nov,D.vir`)
unique(amr_vetted_AG_candidates)

vetted_AG_candidates = list(D.ame = amr_vetted_AG_candidates, 
                             D.lum = lum_vetted_AG_candidates, 
                             D.nov = nov_vetted_AG_candidates, 
                             D.vir = vir_vetted_AG_candidates)
vetted_AG_combs = unlist(lapply(1:length(vetted_AG_candidates), function(j) combn(names(vetted_AG_candidates), j, simplify = FALSE)), recursive = FALSE)
names(vetted_AG_combs) = sapply(vetted_AG_combs, function(i) paste0(i, collapse = ","))
vetted_AG_elements = lapply(vetted_AG_combs, function(i) Setdiff(vetted_AG_candidates[i], vetted_AG_candidates[setdiff(names(vetted_AG_candidates), i)]))
vetted_AG_candidates_Vdiag=venn.diagram(vetted_AG_candidates, NULL, fill=c("#670066", "#86bd38", "#29b5ff", "#717e36"), alpha=c(0.5,0.5,0.5,0.5), cex = 1.5, cat.fontface= 4, cat.cex = 1.25, resolution = 1000, main = "AG-biased")
grid.arrange(gTree(children=vetted_AG_candidates_Vdiag))

###### Output plots
## Dame
# pdf("Plots/dvir1.06_AG-candidates_Dame.pdf", height = 3)
# lapply(AG_vetted_elements$D.ame, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_AG-candidates_Dlum.pdf", height = 3)
# lapply(AG_vetted_elements$D.lum, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_AG-candidates_Dnov.pdf", height = 3)
# lapply(AG_vetted_elements$D.nov, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_vetted_AG-candidates_Dvir.pdf", height = 3)
# lapply(AG_vetted_elements$D.vir, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_vetted_AG-candidates_D.ame-D.lum-D.nov-D.vir.pdf", height = 3)
# lapply(vetted_AG_elements$`D.ame,D.lum,D.nov,D.vir`, plotGeneG, object=TPMse)
# dev.off()
# 
## SFP candidates
pdf("Plots/dvir1.06_SFP-candidates_Dame.pdf", height = 3)
lapply(SFP_elements$D.ame, plotGeneG, object=TPMse)
dev.off()
pdf("Plots/dvir1.06_SFP-candidates_Dlum.pdf", height = 3)
lapply(SFP_elements$D.lum, plotGeneG, object=TPMse)
dev.off()
pdf("Plots/dvir1.06_SFP-candidates_Dnov.pdf", height = 3)
lapply(SFP_elements$D.nov, plotGeneG, object=TPMse)
dev.off()
pdf("Plots/dvir1.06_SFP-candidates_Dvir.pdf", height = 3)
lapply(SFP_elements$D.vir, plotGeneG, object=TPMse)
dev.off()
pdf("Plots/dvir1.06_SFP-candidates_D.ame-D.lum.pdf", height = 3)
lapply(SFP_elements$`D.ame,D.lum`, plotGeneG, object=TPMse)
dev.off()
pdf("Plots/dvir1.06_SFP-candidates_D.ame-D.nov.pdf", height = 3)
lapply(SFP_elements$`D.ame,D.nov`, plotGeneG, object=TPMse)
dev.off()
pdf("Plots/dvir1.06_SFP-candidates_D.ame-D.vir.pdf", height = 3)
lapply(SFP_elements$`D.ame,D.vir`, plotGeneG, object=TPMse)
dev.off()
pdf("Plots/dvir1.06_SFP-candidates_D.lum-D.nov.pdf", height = 3)
lapply(SFP_elements$`D.lum,D.nov`, plotGeneG, object=TPMse)
dev.off()
pdf("Plots/dvir1.06_SFP-candidates_D.lum-D.vir.pdf", height = 3)
lapply(SFP_elements$`D.lum,D.vir`, plotGeneG, object=TPMse)
dev.off()
pdf("Plots/dvir1.06_SFP-candidates_D.nov-D.vir.pdf", height = 3)
lapply(SFP_elements$`D.nov,D.vir`, plotGeneG, object=TPMse)
dev.off()
pdf("Plots/dvir1.06_SFP-candidates_D.ame-D.lum-D.nov.pdf", height = 3)
lapply(SFP_elements$`D.ame,D.lum,D.nov`, plotGeneG, object=TPMse)
dev.off()
pdf("Plots/dvir1.06_SFP-candidates_D.ame-D.lum-D.vir.pdf", height = 3)
lapply(SFP_elements$`D.ame,D.lum,D.vir`, plotGeneG, object=TPMse)
dev.off()
pdf("Plots/dvir1.06_SFP-candidates_D.ame-D.nov-D.vir.pdf", height = 3)
lapply(SFP_elements$`D.ame,D.nov,D.vir`, plotGeneG, object=TPMse)
dev.off()
pdf("Plots/dvir1.06_SFP-candidates_D.lum-D.nov-D.vir.pdf", height = 3)
lapply(SFP_elements$`D.lum,D.nov,D.vir`, plotGeneG, object=TPMse)
dev.off()
pdf("Plots/dvir1.06_SFP-candidates_D.ame-D.lum-D.nov-D.vir.pdf", height = 3)
lapply(SFP_elements$`D.ame,D.lum,D.nov,D.vir`, plotGeneG, object=TPMse)
dev.off()

# ## EB candidates
# pdf("Plots/dvir1.06_EB-candidates_Dame.pdf", height = 3)
# lapply(EB_elements$D.ame, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_EB-candidates_Dlum.pdf", height = 3)
# lapply(EB_elements$D.lum, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_EB-candidates_Dnov.pdf", height = 3)
# lapply(EB_elements$D.nov, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_EB-candidates_Dvir.pdf", height = 3)
# lapply(EB_elements$D.vir, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_EB-candidates_D.ame-D.lum.pdf", height = 3)
# lapply(EB_elements$`D.ame,D.lum`, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_EB-candidates_D.ame-D.nov.pdf", height = 3)
# lapply(EB_elements$`D.ame,D.nov`, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_EB-candidates_D.ame-D.vir.pdf", height = 3)
# lapply(EB_elements$`D.ame,D.vir`, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_EB-candidates_D.lum-D.nov.pdf", height = 3)
# lapply(EB_elements$`D.lum,D.nov`, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_EB-candidates_D.lum-D.vir.pdf", height = 3)
# lapply(EB_elements$`D.lum,D.vir`, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_EB-candidates_D.nov-D.vir.pdf", height = 3)
# lapply(EB_elements$`D.nov,D.vir`, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_EB-candidates_D.ame-D.lum-D.nov.pdf", height = 3)
# lapply(EB_elements$`D.ame,D.lum,D.nov`, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_EB-candidates_D.ame-D.lum-D.vir.pdf", height = 3)
# lapply(EB_elements$`D.ame,D.lum,D.vir`, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_EB-candidates_D.ame-D.nov-D.vir.pdf", height = 3)
# lapply(EB_elements$`D.ame,D.nov,D.vir`, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_EB-candidates_D.lum-D.nov-D.vir.pdf", height = 3)
# lapply(EB_elements$`D.lum,D.nov,D.vir`, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_EB-candidates_D.ame-D.lum-D.nov-D.vir.pdf", height = 3)
# lapply(EB_elements$`D.ame,D.lum,D.nov,D.vir`, plotGeneG, object=TPMse)
# dev.off()
# 
# 
# ## TS candidates
# pdf("Plots/dvir1.06_TS-candidates_Dame.pdf", height = 3)
# lapply(TS_elements$D.ame, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_TS-candidates_Dlum.pdf", height = 3)
# lapply(TS_elements$D.lum, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_TS-candidates_Dnov.pdf", height = 3)
# lapply(TS_elements$D.nov, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_TS-candidates_Dvir.pdf", height = 3)
# lapply(TS_elements$D.vir, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_TS-candidates_D.ame-D.lum.pdf", height = 3)
# lapply(TS_elements$`D.ame,D.lum`, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_TS-candidates_D.ame-D.nov.pdf", height = 3)
# lapply(TS_elements$`D.ame,D.nov`, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_TS-candidates_D.ame-D.vir.pdf", height = 3)
# lapply(TS_elements$`D.ame,D.vir`, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_TS-candidates_D.lum-D.nov.pdf", height = 3)
# lapply(TS_elements$`D.lum,D.nov`, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_TS-candidates_D.lum-D.vir.pdf", height = 3)
# lapply(TS_elements$`D.lum,D.vir`, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_TS-candidates_D.nov-D.vir.pdf", height = 3)
# lapply(TS_elements$`D.nov,D.vir`, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_TS-candidates_D.ame-D.lum-D.nov.pdf", height = 3)
# lapply(TS_elements$`D.ame,D.lum,D.nov`, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_TS-candidates_D.ame-D.lum-D.vir.pdf", height = 3)
# lapply(TS_elements$`D.ame,D.lum,D.vir`, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_TS-candidates_D.ame-D.nov-D.vir.pdf", height = 3)
# lapply(TS_elements$`D.ame,D.nov,D.vir`, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_TS-candidates_D.lum-D.nov-D.vir.pdf", height = 3)
# lapply(TS_elements$`D.lum,D.nov,D.vir`, plotGeneG, object=TPMse)
# dev.off()
# pdf("Plots/dvir1.06_TS-candidates_D.ame-D.lum-D.nov-D.vir.pdf", height = 3)
# lapply(TS_elements$`D.ame,D.lum,D.nov,D.vir`, plotGeneG, object=TPMse)
# dev.off()

# 
########################################################################



##### Produce plots for SFP candidates from Trinity analysis
# pdf("Plots/amr.SFP.list.pdf", height = 3)
# lapply(amr.SFP.list, plotGeneT, object = TPMse_amrTrin)
# dev.off()
# upon inspection, the following are weak candidates: amr_c17060_g4, amr_c20147_g6, amr_c21372_g8
amr.SFP.list = amr.SFP.list[amr.SFP.list != "amr_c17060_g4"] 
amr.SFP.list = amr.SFP.list[amr.SFP.list != "amr_c20147_g6"] 
amr.SFP.list = amr.SFP.list[amr.SFP.list != "amr_c21372_g8"] 

##### Here's a method to find Trinity transcript that are not found in dvir1.06 annotation,
##### then checking whether orthologues exist in the other Trinity assemblies
selectionCols = c("dvir1.06_BlastX_topHit", "gene_id")
amr.SFP.dvir1.06.orths = subset(amrTrinotate, gene_id %in% amr.SFP.list)[selectionCols]
amr.SFP.dvir1.06.orths = droplevels(amr.SFP.dvir1.06.orths)
amr.SFP.dvir1.06.orths = amr.SFP.dvir1.06.orths[order(amr.SFP.dvir1.06.orths$dvir1.06_BlastX_topHit), ]
amr.SFP.dvir1.06.orths[is.na(amr.SFP.dvir1.06.orths)] = "NoHit"
amr.SFP.dvir1.06.orths = aggregate(gene_id~dvir1.06_BlastX_topHit, data = amr.SFP.dvir1.06.orths, toString)
amr.SFP_no_dvir1.06_hits = subset(amr.SFP.dvir1.06.orths, dvir1.06_BlastX_topHit == "NoHit")$gene_id
amr.SFP_no_dvir1.06_hits = unique(strsplit(amr.SFP_no_dvir1.06_hits, ",")[[1]])


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

#pdf(file = "chromDist.figure.pdf", width = 6, height = 2.8)
(chrGG = ggplot(tissue_biased.numbers.c, aes(x=chromosome, y=obs.exp, fill=chromosome)))
(chrGG = chrGG + geom_bar(position=position_dodge(), stat="identity"))
(chrGG = chrGG + scale_fill_manual(values=c("#CC79A7", "#56B4E9", "#56B4E9", "#56B4E9", "#56B4E9")))
(chrGG = chrGG + facet_wrap(~tissue, nrow = 1))
(chrGG = chrGG + geom_errorbar(aes(ymin=obs.exp-ci, ymax=obs.exp+ci), width=.3, position=position_dodge(.9), size = 0.75))
(chrGG = chrGG + ylab("Observed/expected number\n of genes per chromosome"))
(chrGG = chrGG + geom_hline(yintercept = 1, colour = "black", linetype = "dashed"))
(chrGG = chrGG + scale_x_discrete(labels=c("Chr_X" = "X", "Chr_2" = "2", "Chr_3" = "3", "Chr_4" = "4", "Chr_5" = "5")))
(chrGG = chrGG + theme(panel.background = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(face = "bold", size = 14, vjust=0.1), axis.text.y = element_text(face = "bold", size = 12), legend.position="none", strip.text=element_text(face="bold", size = 12), axis.title=element_text(face="bold", size = 12)))
(chrGG = chrGG + geom_text(data = label.05, label = "*", size = 6, colour = "red", fontface=2))
(chrGG = chrGG + geom_text(data = label.01, label = "**", size = 6, colour = "red", fontface=2))
(chrGG = chrGG + geom_text(data = label.001, label = "***", size = 6, colour = "red", fontface=2))
#dev.off()


# example heatmap for gene list
amr.tissueBiased.matrix = subset(grpMeanTPMmatrix, FBgn_ID %in% SFP_elements$`D.ame,D.lum,D.nov,D.vir`)
rownames(amr.tissueBiased.matrix) = amr.tissueBiased.matrix[,1]
amr.tissueBiased.matrix[,1] = NULL
YazHeatmap(amr.tissueBiased.matrix, clustering = "both")

####### Y chromosome analyses #######
Y.genes=as.list(read.table("Annotations/Y_Chromosome_genes.txt"))

################## PAML and KaKs ##################
tmp.FB.names = unique(subset(Annots, select=c("FBgn_ID", "FBtr_ID")))
paml.data = read.csv(file = "PAML.Files/PAML.branchSite.ALL.results.txt", header = T, sep = "\t")
paml.data = merge(tmp.FB.names, paml.data, all=T)
paml.data = merge(gffRecord, paml.data, all=T)
KaKs.data = read.csv(file = "PAML.Files/KaKs.ALL.results.txt", header = T, sep = "\t", check.names = F)
KaKs.data$COMPARISON = paste(KaKs.data$SEQ1, KaKs.data$SEQ2, sep="-")
KaKs.data = merge(tmp.FB.names, KaKs.data, all=T)
KaKs.data = merge(gffRecord, KaKs.data, all=T)

ggplot(subset(KaKs.data, FBgn_ID %in% TS_elements$`D.ame,D.lum,D.nov,D.vir` & Ka/Ks < 40 & grepl("Chr", chromosome)), aes(chromosome, Ka/Ks)) + geom_boxplot()

TS.paml.data = subset(paml.data, FBgn_ID %in% TS_elements$`D.ame,D.lum,D.nov,D.vir`)
TS.paml.data.sig = subset(TS.paml.data, Damr_FDR < 0.05 | Dlum_FDR < 0.05 | Dnov_FDR < 0.05 | Dvir_FDR < 0.05)
TS.paml.data.sig.list = as.character(TS.paml.data.sig$FBgn_ID)

AG.paml.data = subset(paml.data, FBgn_ID %in% AG_elements$`D.ame,D.lum,D.nov,D.vir`)
AG.paml.data.sig = subset(AG.paml.data, Damr_FDR < 0.05 | Dlum_FDR < 0.05 | Dnov_FDR < 0.05 | Dvir_FDR < 0.05)
AG.paml.data.sig.list = as.character(AG.paml.data.sig$FBgn_ID)

# pdf("Plots/SigBranchSite.testes_biased_genes.barPlots.pdf", height = 3)
# lapply(TS.paml.data.sig.list, plotGeneG, object = TPMse)
# dev.off()


ggplot(TS.paml.data.sig, aes(min, omega)) + geom_point() + facet_grid(~chromosome)

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

ggplot(meanOmega.df, aes(Class, omega, colour=Class)) + geom_point(size = 1) + geom_errorbar(aes(ymin=omega-se, ymax=omega+se), width=.1, position=position_dodge(.9)) + theme(legend.position="none")

paml.data$chromosome = factor (paml.data$chromosome, levels = c("Chr_X", "Chr_2", "Chr_3", "Chr_4", "Chr_5"))

gg1=ggplot(subset(paml.data, FBgn_ID %in% SFP_elements$`D.ame,D.lum,D.nov,D.vir` ), aes(max, omega, colour = I("#7aa457"))) + geom_hline(yintercept = genome.avg.omega, linetype="dashed", colour = "red") + geom_point(size=2, alpha=0.75) + facet_grid(~chromosome, scales = "free_x") + geom_point(data = subset(paml.data, FBgn_ID %in% EB_elements$`D.ame,D.lum,D.nov,D.vir`& chromosome != "scaffold_12481"), aes(colour = I("#9e6ebd"))) + geom_text_repel(data=subset(paml.data, FBgn_ID %in% SFP_elements$`D.ame,D.lum,D.nov,D.vir` & omega > 0.8), aes(label = gene_name), size =3, force = 30, colour = "black") + scale_colour_manual(name = "", values =c("#7aa457"="#7aa457","#9e6ebd"="#9e6ebd"), labels = c("SFPs","EB biased")) + scale_x_continuous(breaks=c(5000000, 10000000, 15000000, 20000000, 25000000, 30000000), labels=expression("5", "10", "15", "20", "25", "30")) + xlab ("Chromosome coordinates (Mb)")

gg2=ggplot(subset(paml.data, FBgn_ID %in% TS_elements$`D.ame,D.lum,D.nov,D.vir` & omega < 900 & grepl("Chr", chromosome)), aes(max, omega, colour = I("#7aa457"))) + geom_point(alpha=0.55) + facet_grid(~chromosome, scale = "free_x") + scale_x_continuous(breaks=c(5000000, 10000000, 15000000, 20000000, 25000000, 30000000), labels=expression("5", "10", "15", "20", "25", "30")) + xlab ("Chromosome coordinates (Mb)") + scale_colour_manual(name = "", values ="#cb6751", labels = "testes biased")

plot_grid(gg1, gg2, ncol = 1)

###############################################################################
###############################################################################


## get dvir1.06 matrix of mel_SFP orthologues
virOrths_melSFPs_unique=as.character(unique(subset(melOrths, mel_FBgn_ID %in% melSFPs$mel_FBgn_ID)$FBgn_ID))
virOrths_melSFPs_unique.matrix = subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% virOrths_melSFPs_unique)

jointNames=subset(melOrths, mel_FBgn_ID %in% melSFPs$mel_FBgn_ID)
melSFPs.with.virOrth.and.aggtd.vir.FBgn_ID = aggregate(mel_FBgn_ID~FBgn_ID, data = jointNames, toString)
single_orth_table=as.data.frame(gsub(".*, ", "",melSFPs.with.virOrth.and.aggtd.vir.FBgn_ID$mel_FBgn_ID))
jointNames = cbind(melSFPs.with.virOrth.and.aggtd.vir.FBgn_ID$FBgn_ID, single_orth_table)
colnames(jointNames) = c("FBgn_ID", "mel_FBgn_ID")
#jointNames$FBgn_ID=levels(droplevels(jointNames$FBgn_ID))

melSFPs.RPKM.matrix = subset(melRPKM.reduced.matrix, mel_FBgn_ID %in% jointNames$mel_FBgn_ID)

somethingTmp = merge(jointNames, melSFPs.RPKM.matrix, all=TRUE)
somethingTmp = subset(somethingTmp, Adult_Male_mated_4days_AccGlnd != "NA" & FBgn_ID != "NA")

#unique(subset(somethingTmp, mel_FBgn_ID %in% melSFPs.with.virOrth.and.aggtd.vir.FBgn_ID$mel_FBgn_ID))
somethingTmp2 = merge(somethingTmp, Gene_order, all=TRUE)
melData.to.plot.heatmap=somethingTmp2[order(somethingTmp2$number),]
melData.to.plot.heatmap = subset (melData.to.plot.heatmap, mel_FBgn_ID != "NA" & Adult_Male_mated_4days_AccGlnd != "NA")
melData.to.plot.heatmap = subset (melData.to.plot.heatmap, select = c("FBgn_ID", "Adult_Male_mated_4days_AccGlnd", "Adult_Male_mated_4days_head", "Adult_Male_mated_4days_testis"))
colnames(melData.to.plot.heatmap) = c("FBgn_ID", "AG", "Head", "testis")
rownames(melData.to.plot.heatmap) = melData.to.plot.heatmap$FBgn_ID
melData.to.plot.heatmap[,1] = NULL

virOrths_melSFPs_unique.matrix = subset(grpMeanTPMmatrix.BRR, FBgn_ID %in% jointNames$FBgn_ID)
rownames(virOrths_melSFPs_unique.matrix) = virOrths_melSFPs_unique.matrix$FBgn_ID
virOrths_melSFPs_unique.matrix[,1] = NULL
YazHeatmap(virOrths_melSFPs_unique.matrix, clustering = "both", labRow = T)

#write.table(melSFPs.RPKM.matrix, file = "melSFPs.RPKM.matrix", quote = F, sep = "\t", row.names = T, col.names = T)


colnames(melSFPs.RPKM.matrix)
nrow(melSFPs.RPKM.matrix)
length(melSFPs.with.virOrth.list)
nrow(subset(melSFPs.RPKM.matrix, rownames(melSFPs.RPKM.matrix) %in% melSFPs.with.virOrth.list))
nrow(heatmap_data)

nrow(virOrths_melSFPs_unique.matrix)

aggregate(mel_FBgn_ID~FBgn_ID, data = jointNames, toString)


###### Set up mel.Encode TPM summary
melRPKM.reduced.matrix = subset(melRPKM, select =c("mel_FBgn_ID", "Adult_Male_mated_4days_AccGlnd", "Adult_Male_mated_4days_head", "Adult_Male_mated_4days_testis"))
mel.FBgn_ID_to_GeneSymbol= read.csv("Annotations/mel.FBgn_ID-to-GeneSymbol.txt", header = T, sep = "\t")
melRPKM.tmp=melRPKM
m.melRPKM.tmp = as.data.frame(melt(as.matrix(melRPKM.tmp)))
m.melRPKM.tmp = within(m.melRPKM.tmp, X2=data.frame(do.call('rbind', strsplit(as.character(X2),'_',fixed=TRUE))))
m.melRPKM.tmp=data.frame(m.melRPKM.tmp$X1, m.melRPKM.tmp$X2$X1, m.melRPKM.tmp$X2$X2, m.melRPKM.tmp$X2$X3, m.melRPKM.tmp$X2$X4, m.melRPKM.tmp$X2$X5, m.melRPKM.tmp$value)
colnames(m.melRPKM.tmp) = c("mel_FBgn_ID", "stage", "sex", "status", "age", "tissue", "RPKM")
melRPKM.data = merge(m.melRPKM.tmp, mel.FBgn_ID_to_GeneSymbol, all=TRUE)
melRPKM.data$age =factor(melRPKM.data$age, levels=c("1days","4days","5days","20days","30days"))
#to plot on of the genes:
ggplot(subset(melRPKM.data, mel_FBgn_ID == "FBgn0037039"), aes(tissue, as.numeric(RPKM), fill = status)) + geom_bar(stat = "identity", position = "dodge") + facet_grid(age~sex, scales = "free") + theme(axis.text.x=element_text(angle=45, vjust = 0.1))



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

GO_enrichment_data = do.call(rbind, GO_enriched_list)
