## orthologue heatmaps

## Plot for 
data = virOrths_melSFPs_unique.matrix
sample_types = colnames(data)
nsamples = length(sample_types)
sample_type_list = list()
for (i in 1:nsamples) {
  sample_type_list[[sample_types[i]]] = sample_types[i]
}
sample_factoring = colnames(data)
for (i in 1:nsamples) {
  sample_type = sample_types[i]
  replicates_want = sample_type_list[[sample_type]]
  sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}

# parse gene factors
gene_factor_table = read.table("GO.analysis/tissue-biased-gene-factors.txt", header=F, row.names=2)
gene_factors = unique(gene_factor_table[,1])
gene_factor_colors = rainbow(length(gene_factors))
names(gene_factor_colors) = gene_factors
# data = data[,colSums(data)>=10]
# data = data[rowSums(data)>=10,]
# initial_matrix = data # store before doing various data transformations
data = log2(data+1)
sample_factoring = colnames(data)
for (i in 1:nsamples) {
  sample_type = sample_types[i]
  replicates_want = sample_type_list[[sample_type]]
  sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}
sampleAnnotations = matrix(ncol=ncol(data),nrow=nsamples)
for (i in 1:nsamples) {
  sampleAnnotations[i,] = colnames(data) %in% sample_type_list[[sample_types[i]]]
}
sampleAnnotations = apply(sampleAnnotations, 1:2, function(x) as.logical(x))
sampleAnnotations = sample_matrix_to_color_assignments(sampleAnnotations, col=sample_colors)
rownames(sampleAnnotations) = as.vector(sample_types)
colnames(sampleAnnotations) = colnames(data)
data = as.matrix(data) # convert to matrix
#write.table(data, file="tissue-biased-genes-meanTPMmatrix.txt.minCol10.minRow10.log2.dat", quote=F, sep='	');
sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
#write.table(sample_cor, file="tissue-biased-genes-meanTPMmatrix.txt.minCol10.minRow10.log2.sample_cor.dat", quote=F, sep='	')
sample_dist = as.dist(1-sample_cor)
hc_samples = hclust(sample_dist, method='complete')
gene_cor = NULL
if (is.null(gene_cor)) { gene_cor = cor(t(data), method='pearson', use='pairwise.complete.obs') }
gene_dist = as.dist(1-gene_cor)
if (nrow(data) <= 1) { message('Too few genes to generate heatmap'); quit(status=0); }
hc_genes = hclust(gene_dist, method='complete')
myheatcol = colorpanel(75, 'blue','black','yellow')
 data = t(scale(t(data), scale=F)) # center rows, mean substracted
#write.table(data, file="tissue-biased-genes-meanTPMmatrix.txt.minCol10.minRow10.log2.centered.dat", quote=F, sep='	');
heatmap_data = data
#pdf("tissue-biased-genes-meanTPMmatrix.txt.minCol10.minRow10.log2.centered.genes_vs_samples_heatmap.pdf")
gene_factor_row_vals = as.factor(gene_factor_table[rownames(heatmap_data),])
names(gene_factor_row_vals) = rownames(heatmap_data)
gene_factors_here = unique(gene_factor_row_vals)
names(gene_factors_here) = gene_factors_here
num_gene_factors_here = length(gene_factors_here)
geneFactorColors = c("gray", "purple", "blue", "green", "orange")
if (sum(gene_factors_here %in% sample_types) == num_gene_factors_here) {
  geneFactorColors = sample_colors[names(gene_factors_here)]
}
geneFactorAnnotations = matrix(nrow=nrow(heatmap_data), ncol=num_gene_factors_here)
for (i in 1:num_gene_factors_here) {
  geneFactorAnnotations[,i] = gene_factor_row_vals %in% gene_factors_here[i]
}
geneFactorAnnotations = apply(geneFactorAnnotations, 1:2, function(x) as.logical(x))
geneFactorAnnotations = t(sample_matrix_to_color_assignments(t(geneFactorAnnotations), col=geneFactorColors))
rownames(geneFactorAnnotations) = rownames(heatmap_data)
colnames(geneFactorAnnotations) = gene_factors_here

## set heatmap scale limits here
heatmap_data[heatmap_data < -4] = -4
heatmap_data[heatmap_data > 4] = 4

## plot it
heatmap.3(heatmap_data, dendrogram='both', Rowv=as.dendrogram(hc_genes), Colv=as.dendrogram(hc_samples), col=myheatcol, scale="none", density.info="none", trace="none", key=TRUE, keysize=1.5, cexCol=1, margins=c(10,10), cex.main=0.75, main=paste("samples vs. features", "tissue-biased-genes-meanTPMmatrix.txt.minCol10.minRow10.log2.centered" ) , RowSideColors = geneFactorAnnotations)
#dev.off()


Gene_order <- c("FBgn0198297", "FBgn0198070", "FBgn0211506", "FBgn0211265", "FBgn0211264", "FBgn0210455", "FBgn0208251", "FBgn0211572", "FBgn0201870", "FBgn0208483", "FBgn0207310", "FBgn0199826", "FBgn0208636", "FBgn0204528", "FBgn0210972", "FBgn0208258", "FBgn0208634", "FBgn0197757", "FBgn0211681", "FBgn0208635", "FBgn0210106", "FBgn0207309", "FBgn0198723", "FBgn0204606", "FBgn0201151", "FBgn0198354", "FBgn0206811", "FBgn0211178", "FBgn0210785", "FBgn0211353", "FBgn0209429", "FBgn0198486", "FBgn0209622", "FBgn0198727", "FBgn0203651", "FBgn0197681", "FBgn0203960", "FBgn0208385", "FBgn0201114", "FBgn0199796", "FBgn0204797", "FBgn0197781", "FBgn0204112", "FBgn0208460", "FBgn0203780", "FBgn0208863", "FBgn0202385", "FBgn0211239", "FBgn0014844", "FBgn0207525", "FBgn0200743", "FBgn0208533", "FBgn0209422", "FBgn0205391", "FBgn0208806", "FBgn0202095", "FBgn0208531", "FBgn0200130", "FBgn0207752", "FBgn0207077", "FBgn0206035", "FBgn0205490", "FBgn0201667", "FBgn0211705", "FBgn0209081", "FBgn0198061", "FBgn0210119", "FBgn0205696", "FBgn0203001", "FBgn0208144", "FBgn0197362", "FBgn0211267", "FBgn0209928", "FBgn0208484", "FBgn0204667", "FBgn0199055", "FBgn0211028", "FBgn0199487", "FBgn0205283", "FBgn0204128", "FBgn0211755", "FBgn0207151", "FBgn0207835", "FBgn0208133", "FBgn0200191", "FBgn0210936", "FBgn0201487")
Gene_order = as.data.frame(Gene_order)
Gene_order$number = seq(1,87)
colnames(Gene_order) = c("FBgn_ID", "number")
######################################################################################
######################################################################################

melData.to.plot.heatmap = merge(melMatrix.reduced.SFPs, Gene_order, all=TRUE)
melData.to.plot.heatmap=melData.to.plot.heatmap[order(melData.to.plot.heatmap$number),]
melData.to.plot.heatmap = subset (melData.to.plot.heatmap, mel_FBgn_ID != "NA" & Adult_Male_mated_4days_AccGlnd != "NA")
melData.to.plot.heatmap = subset (melData.to.plot.heatmap, select = c("FBgn_ID", "Adult_Male_mated_4days_AccGlnd", "Adult_Male_mated_4days_head", "Adult_Male_mated_4days_testis"))
colnames(melData.to.plot.heatmap) = c("FBgn_ID", "AG", "Head", "testis")
rownames(melData.to.plot.heatmap) = melData.to.plot.heatmap$FBgn_ID
melData.to.plot.heatmap[,1] = NULL


######################################################################################

data = melData.to.plot.heatmap
sample_types = colnames(data)
nsamples = length(sample_types)
sample_colors = rainbow(nsamples)
sample_type_list = list()
for (i in 1:nsamples) {
  sample_type_list[[sample_types[i]]] = sample_types[i]
}
sample_factoring = colnames(data)
for (i in 1:nsamples) {
  sample_type = sample_types[i]
  replicates_want = sample_type_list[[sample_type]]
  sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}
#data = data[,colSums(data)>=10]
#data = data[rowSums(data)>=1,]
initial_matrix = data # store before doing various data transformations
data = log2(data+1)
sample_factoring = colnames(data)
for (i in 1:nsamples) {
  sample_type = sample_types[i]
  replicates_want = sample_type_list[[sample_type]]
  sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}
sampleAnnotations = matrix(ncol=ncol(data),nrow=nsamples)
for (i in 1:nsamples) {
  sampleAnnotations[i,] = colnames(data) %in% sample_type_list[[sample_types[i]]]
}
sampleAnnotations = apply(sampleAnnotations, 1:2, function(x) as.logical(x))
sampleAnnotations = sample_matrix_to_color_assignments(sampleAnnotations, col=sample_colors)
rownames(sampleAnnotations) = as.vector(sample_types)
colnames(sampleAnnotations) = colnames(data)
data = as.matrix(data) # convert to matrix
# write.table(data, file="ordered.mel.matrix.minCol10.minRow1.log2.dat", quote=F, sep='	');
sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
# write.table(sample_cor, file="ordered.mel.matrix.minCol10.minRow1.log2.sample_cor.dat", quote=F, sep='	')
sample_dist = as.dist(1-sample_cor)
hc_samples = hclust(sample_dist, method='complete')
gene_cor = NULL
if (is.null(gene_cor)) { gene_cor = cor(t(data), method='pearson', use='pairwise.complete.obs') }
gene_dist = as.dist(1-gene_cor)
if (nrow(data) <= 1) { message('Too few genes to generate heatmap'); quit(status=0); }
myheatcol = colorpanel(75, 'blue','black','yellow')
data = t(scale(t(data), scale=F)) # center rows, mean substracted
# write.table(data, file="ordered.mel.matrix.minCol10.minRow1.log2.centered.dat", quote=F, sep='	');
heatmap_data = data
# pdf("ordered.mel.matrix.minCol10.minRow1.log2.centered.genes_vs_samples_heatmap.pdf")
heatmap_data[heatmap_data < -4] = -4
heatmap_data[heatmap_data > 4] = 4
heatmap.3(heatmap_data, dendrogram='col', Rowv=F, Colv=as.dendrogram(hc_samples), col=myheatcol, scale="none", density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=1, margins=c(10,10), cex.main=0.75, main=paste("samples vs. features", "ordered.mel.matrix.minCol10.minRow1.log2.centered" ) )
# dev.off()
