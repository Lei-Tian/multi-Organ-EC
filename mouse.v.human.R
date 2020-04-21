###Comparing bulk RNA -eq of human fetal tissue-specific ECs to scRNA-Seq of corresponding
###adult mouse ECs

##Bulk RNA Seq data analysis
#load packages
require(edgeR)
require(Seurat)
require(AnnotationDbi)
require(org.Hs.eg.db)
require(ggplot2)
require(pheatmap)
require(ggpubr)

#read in count data
seqdata = read.delim("Marcu_Raw_counts.txt", stringsAsFactors = FALSE, header = TRUE)
counts = seqdata[,-1]
rownames(counts) = seqdata[,1]

#get rid of non EC samples
counts = counts[,-(25:30)] 

#convert to counts to CPM
cpm = cpm(counts)

#only keep genes with >0.5 CPM in at least 2 samples
thresh = cpm > 0.5
keep = rowSums(thresh) >= 2
counts.keep = counts[keep,]
logcpm = cpm(counts.keep, log = TRUE)

#change Entrez IDs to gene names and remove NAs
entrez = rownames(logcpm)
symbols = mapIds(org.Hs.eg.db, entrez, 'SYMBOL', 'ENTREZID')
rownames(logcpm) = symbols
a = !is.na(rownames(logcpm))
expr = logcpm[a,]

#split into fresh and cultured
expr.cultured = expr[,1:12]
expr.fresh = expr[,13:24]

#split into organ specific
kidney.cultured = expr.cultured[,c(1,5,9)]
lung.cultured = expr.cultured[,c(2,6,10)]
liver.cultured = expr.cultured[,c(3,7,11)]
heart.cultured = expr.cultured[,c(4,8,12)]

kidney.fresh = expr.fresh[,c(1,5,9)]
lung.fresh = expr.fresh[,c(2,6,10)]
liver.fresh = expr.fresh[,c(3,7,11)]
heart.fresh = expr.fresh[,c(4,8,12)]

##Load in single cell data
EC.kidney = readRDS("EC.kidney.RDS")
EC.liver = readRDS("EC.liver.RDS")
EC.lung = readRDS("EC.lung.RDS")
EC.heart = readRDS("EC.heart.RDS")

#get matrices of single cell expression data
kidney.sc.expr = as.matrix(EC.kidney@data)
liver.sc.expr = as.matrix(EC.liver@data)
lung.sc.expr = as.matrix(EC.lung@data)
heart.sc.expr = as.matrix(EC.heart@data)

#average across all cells / samples 
kidney.sc = data.frame("sc.exprs"=rowMeans(kidney.sc.expr))
liver.sc = data.frame("sc.exprs"=rowMeans(liver.sc.expr))
lung.sc = data.frame("sc.exprs"=rowMeans(lung.sc.expr))
heart.sc = data.frame("sc.exprs"=rowMeans(heart.sc.expr))

kidney.b = data.frame("b.exprs"=rowMeans(kidney.fresh))
liver.b = data.frame("b.exprs"=rowMeans(liver.fresh))
lung.b = data.frame("b.exprs"=rowMeans(lung.fresh))
heart.b = data.frame("b.exprs"=rowMeans(heart.fresh))

#capitalize rownames of sc data to allow for subsequent merge
rownames(kidney.sc) = toupper(rownames(kidney.sc))
rownames(liver.sc) = toupper(rownames(liver.sc))
rownames(lung.sc) = toupper(rownames(lung.sc))
rownames(heart.sc) = toupper(rownames(heart.sc))

#merge single cell and bulk datasets
kidney.sc.b = merge(kidney.sc, kidney.b, by.x=0, by.y=0)
liver.sc.b = merge(liver.sc, liver.b, by.x=0, by.y=0)
lung.sc.b = merge(lung.sc, lung.b, by.x=0, by.y=0)
heart.sc.b = merge(heart.sc, heart.b, by.x=0, by.y=0)

#re-format dataframe
rownames(kidney.sc.b) = kidney.sc.b[,1]
kidney.sc.b=kidney.sc.b[,-1]

rownames(liver.sc.b) = liver.sc.b[,1]
liver.sc.b=liver.sc.b[,-1]

rownames(lung.sc.b) = lung.sc.b[,1]
lung.sc.b=lung.sc.b[,-1]

rownames(heart.sc.b) = heart.sc.b[,1]
heart.sc.b=heart.sc.b[,-1]

#plot
kidney.cor.plot = ggscatter(kidney.sc.b,  x = "sc.exprs", y = "b.exprs", shape = 19, size =0.5, ylab = expression(paste("Human bulk RNA-Seq ( log CPM )" )), xlab = expression(paste("Mouse sc RNA Seq (log" ["2"], " counts)")))
kidney.cor.plot + stat_smooth()

liver.cor.plot = ggscatter(liver.sc.b,  x = "sc.exprs", y = "b.exprs", shape = 19, size =0.5, ylab = expression(paste("Human bulk RNA-Seq ( log CPM )" )), xlab = expression(paste("Mouse sc RNA Seq (log" ["2"], " counts)")))
liver.cor.plot + stat_smooth()

lung.cor.plot = ggscatter(lung.sc.b,  x = "sc.exprs", y = "b.exprs", shape = 19, size =0.5, ylab = expression(paste("Human bulk RNA-Seq ( log CPM )" )), xlab = expression(paste("Mouse sc RNA Seq (log" ["2"], " counts)")))
lung.cor.plot + stat_smooth()

heart.cor.plot = ggscatter(heart.sc.b,  x = "sc.exprs", y = "b.exprs", shape = 19, size =0.5, ylab = expression(paste("Human bulk RNA-Seq ( log CPM )" )), xlab = expression(paste("Mouse sc RNA Seq (log" ["2"], " counts)")))
heart.cor.plot + stat_smooth()

#calculate correlation coefficients
kidney.cor = cor.test(kidney.sc.b$sc.exprs, kidney.sc.b$b.exprs, method = 'spearman')
liver.cor = cor.test(liver.sc.b$sc.exprs, liver.sc.b$b.exprs, method = 'spearman')
lung.cor = cor.test(lung.sc.b$sc.exprs, lung.sc.b$b.exprs, method = 'spearman')
heart.cor = cor.test(heart.sc.b$sc.exprs, heart.sc.b$b.exprs, method = 'spearman')

##Make correlation heatmap
#merge bulk and single cell by all pairwise comparisons
ks.kb = merge(kidney.sc, kidney.b, by.x=0, by.y=0)
rownames(ks.kb) = ks.kb[,1]
ks.kb=ks.kb[,-1]

ks.lib = merge(kidney.sc, liver.b, by.x=0, by.y=0)
rownames(ks.lib) = ks.lib[,1]
ks.lib=ks.lib[,-1]

ks.lub= merge(kidney.sc, lung.b, by.x=0, by.y=0)
rownames(ks.lub) = ks.lub[,1]
ks.lub=ks.lub[,-1]

ks.hb = merge(kidney.sc, heart.b, by.x=0, by.y=0)
rownames(ks.hb) = ks.hb[,1]
ks.hb=ks.hb[,-1]

lis.kb = merge(liver.sc, kidney.b, by.x=0, by.y=0)
rownames(lis.kb) = lis.kb[,1]
lis.kb=lis.kb[,-1]

lis.lib = merge(liver.sc, liver.b, by.x=0, by.y=0)
rownames(lis.lib) = lis.lib[,1]
lis.lib=lis.lib[,-1]

lis.lub= merge(liver.sc, lung.b, by.x=0, by.y=0)
rownames(lis.lub) = lis.lub[,1]
lis.lub=lis.lub[,-1]

lis.hb = merge(liver.sc, heart.b, by.x=0, by.y=0)
rownames(lis.hb) = lis.hb[,1]
lis.hb=lis.hb[,-1]

lus.kb = merge(lung.sc, kidney.b, by.x=0, by.y=0)
rownames(lus.kb) = lus.kb[,1]
lus.kb=lus.kb[,-1]

lus.lib = merge(lung.sc, liver.b, by.x=0, by.y=0)
rownames(lus.lib) = lus.lib[,1]
lus.lib=lus.lib[,-1]

lus.lub = merge(lung.sc, lung.b, by.x=0, by.y=0)
rownames(lus.lub) = lus.lub[,1]
lus.lub=lus.lub[,-1]

lus.hb = merge(lung.sc, heart.b, by.x=0, by.y=0)
rownames(lus.hb) = lus.hb[,1]
lus.hb=lus.hb[,-1]

hs.kb = merge(heart.sc, kidney.b, by.x=0, by.y=0)
rownames(hs.kb) = hs.kb[,1]
hs.kb=hs.kb[,-1]

hs.lib = merge(heart.sc, liver.b, by.x=0, by.y=0)
rownames(hs.lib) = hs.lib[,1]
hs.lib=hs.lib[,-1]

hs.lub = merge(heart.sc, lung.b, by.x=0, by.y=0)
rownames(hs.lub) = hs.lub[,1]
hs.lub=hs.lub[,-1]

hs.hb = merge(heart.sc, heart.b, by.x=0, by.y=0)
rownames(hs.hb) = hs.hb[,1]
hs.hb=hs.hb[,-1]

#Make correlation matrix
comparison.mat = matrix(c("hs.hb","hs.kb", "hs.lib", "hs.lub",  "ks.hb", "ks.kb",
                          "ks.lib", "ks.lub", "lis.hb", "lis.kb", "lis.lib", "lis.lub",
                          "lus.hb", "lus.kb", "lus.lib", "lus.lub"), nrow =4,
                        ncol = 4)

cor.mat.spearman = matrix(data = NA, nrow = nrow(comparison.mat), 
                          ncol = ncol(comparison.mat))

for (i in 1:ncol(comparison.mat)){
  for(j in 1:nrow(comparison.mat)){
    cor.mat.spearman[i,j] = cor(get(comparison.mat[i,j]), method = "spearman")[2,1]
  }
}

#format correlation matrix for ggplot2
rownames(cor.mat.spearman) = c("Heart", "Kidney", "Liver", "Lung")
colnames(cor.mat.spearman) = c("Heart", "Kidney", "Liver", "Lung")


#make heatmap for global correlation matrix in Figure 7A
pheatmap(
  mat = cor.mat.spearman, 
  cluster_rows = FALSE, 
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE)

##Make heatmap showing expression of tissue-specific EC markers, identified in the Tabula
##Muris dataset, in human fetal ECs

#list of tissue-specific EC DEGs from Figure 2A
heart.degs = c("WT1", "SLC28A2", "EEPD1", "KCNA5", "CA8", "FBLN1", "MEOX2", "RFTN1",
               "LAMB1", "MYADM")
kidney.degs = c("DRAM1", "DKK2", "ESM1", "IGFBP5", "PBX1", "BOC", "IGFBP3", "IRX3",
                "TNFAIP2", "PTPRU")
liver.degs = c("DNASE1L3", "CLEC4G", "FCGR2B", "STAB2", "OIT3", "BMP2", "AASS", 
               "MRC1", "PLXNC1", "WNT2")
lung.degs = c("GRTP1", "ADRB1", "SCN7A", "TMEM100", "HPGD", "FOXF1", "NCKAP5",
              "RASGEF1A", "FENDRR", "PRX")

degs.all = c(heart.degs, kidney.degs, liver.degs, lung.degs)

#arrange bulk RNA seq data by tissue
expr.ordered = cbind(heart.fresh, kidney.fresh, liver.fresh,lung.fresh)

#pull out gene expression data for tissue-specific DEGs identified in 
b.expr.s.degs = expr.ordered[(rownames(expr.ordered) %in% degs.all),,drop=FALSE]
b.expr.s.degs = b.expr.s.degs[match(degs.all,rownames(b.expr.s.degs)),]

colnames(b.expr.s.degs) = c("Heart 1", "Heart 2", "Heart 3", "Kidney 1", "Kidney 2",
                            "Kidney 3", "Liver 1", "Liver 2", "Liver 3", "Lung 1",
                            "Lung 2", "Lung 3")

#draw heatmap
pheatmap(
  mat = b.expr.s.degs, 
  cluster_rows = FALSE, 
  cluster_cols = FALSE,
  scale="row",
  gaps_row=c(10,20,30),
  gaps_col=c(3,6,9))



