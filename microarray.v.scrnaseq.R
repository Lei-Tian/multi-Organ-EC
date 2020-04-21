###COMPARE GENE EXPRESSION DATA BETWEEN MICROARRAY (Nolan et al. Dev Cell 2013) 
###AND SC-RNA-SEQ (Tabula Muris)
#load packages
require("AnnotationDbi")
require("BiocManager")
require("dplyr")
require("oligo")
require("pd.mogene.1.0.st.v1")
require("rstudioapi")
require("mogene10sttranscriptcluster.db")
require("data.table")
require("Seurat")
require("tcltk")
require("ggpubr")
require("ggplot2")
require("Biobase")
require("dendextend")
require("dendsort")
require("RColorBrewer")
require("stringr")
require("oligoClasses")
require("limma")
require("gplots")
require("geneplotter")


##SCATTERPLOTS OF GLOBAL GENE EXPRESSION IN TISSUE-SPECIFIC ECS AS DETERMINED BY MICROARRAY AND
##SCRNA-SEQ - SUPPLEMENTAL FIGURE 3A-F

#ask user for directory containing organ-specific folders
path = tk_choose.dir(caption = "Select Directory Containing Organ-Specific Cel Files")

#get probeset ID info from annotation package
annot = as.data.frame(mogene10sttranscriptclusterSYMBOL)

#loop through each folder, read in raw cel data and average
folder.names = dir(path)

for (i in 1:length(folder.names)){
  setwd(file.path(path, paste(folder.names[i])))
  cels  = list.celfiles(full.names=TRUE)
  gfs = read.celfiles(cels)
  gfs.rma = rma(gfs, target = "core")
  gfs.rma.avg = data.frame(exprs = rowMeans(gfs.rma@assayData$exprs))
  gfs.rma.avg.merge = data.frame(exprs = merge(annot, gfs.rma.avg, by.x=1, by.y=0))
  assign(paste(folder.names[i],".ma", sep=""), data.frame(ma.exprs = gfs.rma.avg.merge[,3], symbols = gfs.rma.avg.merge[,2], row.names = NULL)) 
  
}

rm(gfs, gfs.rma, gfs.rma.avg, gfs.rma.avg.merge)

#Get Average Raw scRNA-Seq Data For Each Organ
#read in Seurat object, pull out raw data from male mice and convert to a data frame
file = tk_choose.files(caption = "Select Seurat Object")
EC=readRDS(file)
sc.data = as.matrix(GetAssayData(EC, slot ="data"))
mouse.sex=EC@meta.data$mouse.sex
mouse.male=mouse.sex %in% "M"
sc.data = sc.data[,mouse.male]

#get raw data from each organ
brain.sc = as.matrix(sc.data[,grep("Brain", colnames(sc.data))])
brain.sc=data.frame("sc.exprs"=rowMeans(brain.sc))

heart.sc = as.matrix(sc.data[,grep("Heart", colnames(sc.data))])
heart.sc=data.frame("sc.exprs"=rowMeans(heart.sc))

kidney.sc = as.matrix(sc.data[,grep("Kidney", colnames(sc.data))])
kidney.sc = data.frame("sc.exprs"=rowMeans(kidney.sc))

liver.sc = as.matrix(sc.data[,grep("Liver", colnames(sc.data))])
liver.sc = data.frame("sc.exprs"=rowMeans(liver.sc))

lung.sc = as.matrix(sc.data[,grep("Lung", colnames(sc.data))])
lung.sc = data.frame("sc.exprs"=rowMeans(lung.sc))

muscle.sc = as.matrix(sc.data[,grep("Muscle", colnames(sc.data))])
muscle.sc = data.frame("sc.exprs"=rowMeans(muscle.sc))

#merge microarray and sc RNA Seq data sets
brain.ma.sc = merge(brain.ma, brain.sc, by.x = "symbols", by.y=0)
colnames(brain.ma.sc)=c("symbol", "ma.expression", "sc.expression")
brain.df = as.data.frame(brain.ma.sc[,2:3])

heart.ma.sc = merge(heart.ma, heart.sc, by.x = "symbols", by.y=0)
colnames(heart.ma.sc)=c("symbol", "ma.expression", "sc.expression")
heart.df = as.data.frame(heart.ma.sc[,2:3])

kidney.ma.sc = merge(kidney.ma, kidney.sc, by.x = "symbols", by.y=0)
colnames(kidney.ma.sc)=c("symbol", "ma.expression", "sc.expression")
kidney.df = as.data.frame(kidney.ma.sc[,2:3])

liver.ma.sc = merge(liver.ma, liver.sc, by.x = "symbols", by.y=0)
colnames(liver.ma.sc)=c("symbol", "ma.expression", "sc.expression")
liver.df = as.data.frame(liver.ma.sc[,2:3])

lung.ma.sc = merge(lung.ma, lung.sc, by.x = "symbols", by.y=0)
colnames(lung.ma.sc)=c("symbol", "ma.expression", "sc.expression")
lung.df = as.data.frame(lung.ma.sc[,2:3])

muscle.ma.sc = merge(muscle.ma, muscle.sc, by.x = "symbols", by.y=0)
colnames(muscle.ma.sc)=c("symbol", "ma.expression", "sc.expression")
muscle.df = as.data.frame(muscle.ma.sc[,2:3])

#plot global correlations - Supplemental Figure 3A-F
brain.cor.plot = ggscatter(brain.df,  x = "sc.expression", y = "ma.expression", shape = 19, size =0.5, ylab = expression(paste("Microarray ( log" ["2"], " intensity)")), xlab = expression(paste("sc RNA Seq (log" ["2"], " counts)")))
brain.cor.plot + stat_smooth()

heart.cor.plot = ggscatter(heart.df,  x = "sc.expression", y = "ma.expression", shape = 19, size =0.5, ylab = expression(paste("Microarray ( log" ["2"], " intensity)")), xlab = expression(paste("sc RNA Seq (log" ["2"], " counts)")))
heart.cor.plot + stat_smooth()

kidney.cor.plot = ggscatter(kidney.df,  x = "sc.expression", y = "ma.expression", shape = 19, size =0.5, ylab = expression(paste("Microarray ( log" ["2"], " intensity)")), xlab = expression(paste("sc RNA Seq (log" ["2"], " counts)")))
kidney.cor.plot + stat_smooth()

liver.cor.plot = ggscatter(liver.df,  x = "sc.expression", y = "ma.expression", shape = 19, size =0.5, ylab = expression(paste("Microarray ( log" ["2"], " intensity)")), xlab = expression(paste("sc RNA Seq (log" ["2"], " counts)")))
liver.cor.plot + stat_smooth()

lung.cor.plot = ggscatter(lung.df,  x = "sc.expression", y = "ma.expression", shape = 19, size =0.5, ylab = expression(paste("Microarray ( log" ["2"], " intensity)")), xlab = expression(paste("sc RNA Seq (log" ["2"], " counts)")))
lung.cor.plot + stat_smooth()

muscle.cor.plot = ggscatter(muscle.df,  x = "sc.expression", y = "ma.expression", shape = 19, size =0.5, ylab = expression(paste("Microarray ( log" ["2"], " intensity)")), xlab = expression(paste("sc RNA Seq (log" ["2"], " counts)")))
muscle.cor.plot + stat_smooth()


##MAKE HEATMAP OF TISSUE-SPECIFIC EC CORRELATION COEFFICIENTS BETWEEN MICROARRAY AND SCRNA-SEQ

#merge microarray and sc RNA Seq data sets for all pairwise comparisons
bm.bs = merge(brain.ma, brain.sc, by.x = "symbols", by.y=0)
colnames(bm.bs)=c("symbol", "ma.expression", "sc.expression")
bm.bs = as.data.frame(bm.bs[,2:3])

bm.hs = merge(brain.ma, heart.sc, by.x = "symbols", by.y=0)
colnames(bm.hs)=c("symbol", "ma.expression", "sc.expression")
bm.hs = as.data.frame(bm.hs[,2:3])

bm.ks = merge(brain.ma, kidney.sc, by.x = "symbols", by.y=0)
colnames(bm.ks)=c("symbol", "ma.expression", "sc.expression")
bm.ks = as.data.frame(bm.ks[,2:3])

bm.lis = merge(brain.ma, liver.sc, by.x = "symbols", by.y=0)
colnames(bm.lis)=c("symbol", "ma.expression", "sc.expression")
bm.lis = as.data.frame(bm.lis[,2:3])

bm.lus = merge(brain.ma, lung.sc, by.x = "symbols", by.y=0)
colnames(bm.lus)=c("symbol", "ma.expression", "sc.expression")
bm.lus = as.data.frame(bm.lus[,2:3])

bm.ms = merge(brain.ma, muscle.sc, by.x = "symbols", by.y=0)
colnames(bm.ms)=c("symbol", "ma.expression", "sc.expression")
bm.ms = as.data.frame(bm.ms[,2:3])

hm.bs = merge(heart.ma, brain.sc, by.x = "symbols", by.y=0)
colnames(hm.bs)=c("symbol", "ma.expression", "sc.expression")
hm.bs = as.data.frame(hm.bs[,2:3])

hm.hs = merge(heart.ma, heart.sc, by.x = "symbols", by.y=0)
colnames(hm.hs)=c("symbol", "ma.expression", "sc.expression")
hm.hs = as.data.frame(hm.hs[,2:3])

hm.ks = merge(heart.ma, kidney.sc, by.x = "symbols", by.y=0)
colnames(hm.ks)=c("symbol", "ma.expression", "sc.expression")
hm.ks = as.data.frame(hm.ks[,2:3])

hm.lis = merge(heart.ma, liver.sc, by.x = "symbols", by.y=0)
colnames(hm.lis)=c("symbol", "ma.expression", "sc.expression")
hm.lis = as.data.frame(hm.lis[,2:3])

hm.lus = merge(heart.ma, lung.sc, by.x = "symbols", by.y=0)
colnames(hm.lus)=c("symbol", "ma.expression", "sc.expression")
hm.lus = as.data.frame(hm.lus[,2:3])

hm.ms = merge(heart.ma, muscle.sc, by.x = "symbols", by.y=0)
colnames(hm.ms)=c("symbol", "ma.expression", "sc.expression")
hm.ms = as.data.frame(hm.ms[,2:3])

km.bs = merge(kidney.ma, brain.sc, by.x = "symbols", by.y=0)
colnames(km.bs)=c("symbol", "ma.expression", "sc.expression")
km.bs = as.data.frame(km.bs[,2:3])

km.hs = merge(kidney.ma, heart.sc, by.x = "symbols", by.y=0)
colnames(km.hs)=c("symbol", "ma.expression", "sc.expression")
km.hs = as.data.frame(km.hs[,2:3])

km.ks = merge(kidney.ma, kidney.sc, by.x = "symbols", by.y=0)
colnames(km.ks)=c("symbol", "ma.expression", "sc.expression")
km.ks = as.data.frame(km.ks[,2:3])

km.lis = merge(kidney.ma, liver.sc, by.x = "symbols", by.y=0)
colnames(km.lis)=c("symbol", "ma.expression", "sc.expression")
km.lis = as.data.frame(km.lis[,2:3])

km.lus = merge(kidney.ma, lung.sc, by.x = "symbols", by.y=0)
colnames(km.lus)=c("symbol", "ma.expression", "sc.expression")
km.lus = as.data.frame(km.lus[,2:3])

km.ms = merge(kidney.ma, muscle.sc, by.x = "symbols", by.y=0)
colnames(km.ms)=c("symbol", "ma.expression", "sc.expression")
km.ms = as.data.frame(km.ms[,2:3])

lim.bs = merge(liver.ma, brain.sc, by.x = "symbols", by.y=0)
colnames(lim.bs)=c("symbol", "ma.expression", "sc.expression")
lim.bs = as.data.frame(lim.bs[,2:3])

lim.hs = merge(liver.ma, heart.sc, by.x = "symbols", by.y=0)
colnames(lim.hs)=c("symbol", "ma.expression", "sc.expression")
lim.hs = as.data.frame(lim.hs[,2:3])

lim.ks = merge(liver.ma, kidney.sc, by.x = "symbols", by.y=0)
colnames(lim.ks)=c("symbol", "ma.expression", "sc.expression")
lim.ks = as.data.frame(lim.ks[,2:3])

lim.lis = merge(liver.ma, liver.sc, by.x = "symbols", by.y=0)
colnames(lim.lis)=c("symbol", "ma.expression", "sc.expression")
lim.lis = as.data.frame(lim.lis[,2:3])

lim.lus = merge(liver.ma, lung.sc, by.x = "symbols", by.y=0)
colnames(lim.lus)=c("symbol", "ma.expression", "sc.expression")
lim.lus = as.data.frame(lim.lus[,2:3])

lim.ms = merge(liver.ma, muscle.sc, by.x = "symbols", by.y=0)
colnames(lim.ms)=c("symbol", "ma.expression", "sc.expression")
lim.ms = as.data.frame(lim.ms[,2:3])

lum.bs = merge(lung.ma, brain.sc, by.x = "symbols", by.y=0)
colnames(lum.bs)=c("symbol", "ma.expression", "sc.expression")
lum.bs = as.data.frame(lum.bs[,2:3])

lum.hs = merge(lung.ma, heart.sc, by.x = "symbols", by.y=0)
colnames(lum.hs)=c("symbol", "ma.expression", "sc.expression")
lum.hs = as.data.frame(lum.hs[,2:3])

lum.ks = merge(lung.ma, kidney.sc, by.x = "symbols", by.y=0)
colnames(lum.ks)=c("symbol", "ma.expression", "sc.expression")
lum.ks = as.data.frame(lum.ks[,2:3])

lum.lis = merge(lung.ma, liver.sc, by.x = "symbols", by.y=0)
colnames(lum.lis)=c("symbol", "ma.expression", "sc.expression")
lum.lis = as.data.frame(lum.lis[,2:3])

lum.lus = merge(lung.ma, lung.sc, by.x = "symbols", by.y=0)
colnames(lum.lus)=c("symbol", "ma.expression", "sc.expression")
lum.lus = as.data.frame(lum.lus[,2:3])

lum.ms = merge(lung.ma, muscle.sc, by.x = "symbols", by.y=0)
colnames(lum.ms)=c("symbol", "ma.expression", "sc.expression")
lum.ms = as.data.frame(lum.ms[,2:3])

mm.bs = merge(muscle.ma, brain.sc, by.x = "symbols", by.y=0)
colnames(mm.bs)=c("symbol", "ma.expression", "sc.expression")
mm.bs = as.data.frame(mm.bs[,2:3])

mm.hs = merge(muscle.ma, heart.sc, by.x = "symbols", by.y=0)
colnames(mm.hs)=c("symbol", "ma.expression", "sc.expression")
mm.hs = as.data.frame(mm.hs[,2:3])

mm.ks = merge(muscle.ma, kidney.sc, by.x = "symbols", by.y=0)
colnames(mm.ks)=c("symbol", "ma.expression", "sc.expression")
mm.ks = as.data.frame(mm.ks[,2:3])

mm.lis = merge(muscle.ma, liver.sc, by.x = "symbols", by.y=0)
colnames(mm.lis)=c("symbol", "ma.expression", "sc.expression")
mm.lis = as.data.frame(mm.lis[,2:3])

mm.lus = merge(muscle.ma, lung.sc, by.x = "symbols", by.y=0)
colnames(mm.lus)=c("symbol", "ma.expression", "sc.expression")
mm.lus = as.data.frame(mm.lus[,2:3])

mm.ms = merge(muscle.ma, muscle.sc, by.x = "symbols", by.y=0)
colnames(mm.ms)=c("symbol", "ma.expression", "sc.expression")
mm.ms = as.data.frame(mm.ms[,2:3])

#make correlation matrix
comparison.mat = t(matrix(c("bm.bs", "bm.hs", "bm.ks", "bm.lis", "bm.lus", "bm.ms", "hm.bs", "hm.hs", "hm.ks", "hm.lis", "hm.lus", "hm.ms", "km.bs", "km.hs", "km.ks", "km.lis", "km.lus", "km.ms", "lim.bs", "lim.hs", "lim.ks", "lim.lis", "lim.lus", "lim.ms", "lum.bs", "lum.hs", "lum.ks", "lum.lis", "lum.lus", "lum.ms", "mm.bs", "mm.hs", "mm.ks", "mm.lis", "mm.lus", "mm.ms"), nrow = 6, ncol = 6))
cor.mat.spearman = matrix(data = NA, nrow = nrow(comparison.mat), ncol = ncol(comparison.mat))

for (i in 1:ncol(comparison.mat)){
  for(j in 1:nrow(comparison.mat)){
    cor.mat.spearman[i,j] = cor(get(comparison.mat[i,j]), method = "spearman")[2,1]
  }
}

cor.df.spearman = as.data.frame(cor.mat.spearman)

#generate heatmap seen in Supplemental Figure 3G
pheatmap(
  mat = cor.df.spearman, 
  cluster_rows = FALSE, 
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE)

##SCATTERPLOT OF EXPRESSION OF TOP 100 TISSUE-SPECIFIC EC MARKERS, DETERMINED FROM OVERLAP
##BETWEEN MICROARRAY AND SCRNASEQ DATA, IN TISSUE-SPECIFIC ECS ANALYZED BY MICROARRAY AND 
##SCRNA-SEQ

#Rread in CEL files
path = tk_choose.dir(caption = "Select Directory Containing All Cel Files")
cels = list.celfiles(path, full.names=TRUE)
raw.data = read.celfiles(cels)

#run RMA algorithm
data.ma = oligo::rma(raw.data, target = "core")

#Filter out transcripts that do not have intensities >4 in at least 3 arrays.
#This intensity value is based on visual inspection of the histogram of genewise 
#median intensities
man.threshold = 4
samples.cutoff = 3
idx.man.threshold = apply(Biobase::exprs(data.ma), 1, function(x){
  +     sum(x>man.threshold) >= samples.cutoff})
data.manfiltered = subset(data.ma, idx.man.threshold)

#annotate transcript clusters
annot = AnnotationDbi::select(mogene10sttranscriptcluster.db, keys = featureNames(data.manfiltered), columns = c("SYMBOL", "GENENAME"), keytype = "PROBEID")
annot = subset(annot, !is.na(SYMBOL))

#group annot by PROBEID
annot.grouped = group_by(annot, PROBEID)
#summarize the groups and find # genes that PROBEIDs are mapped to
annot.summarized = dplyr::summarize(annot.grouped, no_of_matches = n_distinct(SYMBOL))
#remove transcripts with multiple mappings
annot.filtered = filter(annot.summarized, no_of_matches > 1)
probe.stats = annot.filtered
ids.to.exclude = (featureNames(data.manfiltered) %in% probe.stats$PROBEID)
data.ma.final = subset(data.manfiltered, !ids.to.exclude)

#add the SYMBOL and GENENAME columns to the feature data
fData(data.ma.final)$PROBEID = rownames(fData(data.ma.final))
fData(data.ma.final) = left_join(fData(data.ma.final), annot)
rownames(fData(data.ma.final)) = fData(data.ma.final)$PROBEID

#Find markers of tissue-specific ECs from microarray data
#create experimental design matrix
design.matrix = model.matrix(~0 + factor(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6)))
colnames(design.matrix) = c("Brain", "Heart", "Kidney", "Liver", "Lung", "Muscle")

#fit data to linear model specififed by the design matrix
fit = lmFit(data.ma.final, design.matrix)

#make contrast matrix to setup pairwise comparisons
contrast.matrix.brain = makeContrasts(Brain-(Heart+Kidney+Liver+Lung+Muscle)/5, levels = design.matrix)
contrast.matrix.heart = makeContrasts(Heart-(Brain+Kidney+Liver+Lung+Muscle)/5, levels = design.matrix)
contrast.matrix.kidney = makeContrasts(Kidney-(Brain+Heart+Liver+Lung+Muscle)/5, levels = design.matrix) 
contrast.matrix.liver = makeContrasts(Liver-(Brain+Heart+Kidney+Lung+Muscle)/5, levels = design.matrix)
contrast.matrix.lung = makeContrasts(Lung-(Brain+Heart+Kidney+Liver+Muscle)/5, levels = design.matrix)
contrast.matrix.muscle = makeContrasts(Muscle-(Brain+Heart+Kidney+Liver+Lung)/5, levels = design.matrix)

#given a linear model fit to microarray data, compute estimated coefficients and SEs for 
#a given contrast matrix
fit2.brain = contrasts.fit(fit, contrast.matrix.brain)
fit2.heart = contrasts.fit(fit, contrast.matrix.heart)
fit2.kidney = contrasts.fit(fit, contrast.matrix.kidney)
fit2.liver = contrasts.fit(fit, contrast.matrix.liver)
fit2.lung = contrasts.fit(fit, contrast.matrix.lung)
fit2.muscle = contrasts.fit(fit, contrast.matrix.muscle)

#calculate moderated t-statistic and log-odds of diff expression
fit2.brain = eBayes(fit2.brain)
fit2.heart = eBayes(fit2.heart)
fit2.kidney = eBayes(fit2.kidney)
fit2.liver = eBayes(fit2.liver)
fit2.lung = eBayes(fit2.lung)
fit2.muscle = eBayes(fit2.muscle)

##Scatterplots
#get tables of top 150 DEGs
brain.top = as.data.frame(topTable(fit2.brain, number = 150, adjust = "BH"))
heart.top = as.data.frame(topTable(fit2.heart, number = 150, adjust = "BH"))
kidney.top = as.data.frame(topTable(fit2.kidney, number = 150, adjust = "BH"))
liver.top = as.data.frame(topTable(fit2.liver, number = 150, adjust = "BH"))
lung.top = as.data.frame(topTable(fit2.lung, number = 150, adjust = "BH"))
muscle.top = as.data.frame(topTable(fit2.muscle, number = 150, adjust = "BH"))

#get lists of top 150 DEGs
brain.top.genes = brain.top$SYMBOL
heart.top.genes = heart.top$SYMBOL
kidney.top.genes = kidney.top$SYMBOL
liver.top.genes = liver.top$SYMBOL
lung.top.genes = lung.top$SYMBOL
muscle.top.genes = muscle.top$SYMBOL

##Get Average scRNA-Seq Data For Each Organ
#read in Seurat object, pull out raw data from male mice and convert to a data frame
file = tk_choose.files(caption = "Select Seurat Object")
EC=readRDS(file)
sc.data = as.matrix(GetAssayData(EC, slot = "data"))
mouse.sex=EC@meta.data$mouse.sex
mouse.male=mouse.sex %in% "M"
sc.data = sc.data[,mouse.male]

#get raw data from each organ
brain.sc = as.matrix(sc.data[,grep("Brain", colnames(sc.data))])
brain.sc=data.frame("sc.exprs"=rowMeans(brain.sc))
brain.sc = cbind("symbols" = rownames(brain.sc), brain.sc)
row.names(brain.sc)=NULL

heart.sc = as.matrix(sc.data[,grep("Heart", colnames(sc.data))])
heart.sc=data.frame("sc.exprs"=rowMeans(heart.sc))
heart.sc = cbind("symbols" = rownames(heart.sc), heart.sc)
row.names(heart.sc)=NULL

kidney.sc = as.matrix(sc.data[,grep("Kidney", colnames(sc.data))])
kidney.sc = data.frame("sc.exprs"=rowMeans(kidney.sc))
kidney.sc = cbind("symbols" = rownames(kidney.sc), kidney.sc)
row.names(kidney.sc)=NULL

liver.sc = as.matrix(sc.data[,grep("Liver", colnames(sc.data))])
liver.sc = data.frame("sc.exprs"=rowMeans(liver.sc))
liver.sc = cbind("symbols" = rownames(liver.sc), liver.sc)
row.names(liver.sc)=NULL

lung.sc = as.matrix(sc.data[,grep("Lung", colnames(sc.data))])
lung.sc = data.frame("sc.exprs"=rowMeans(lung.sc))
lung.sc = cbind("symbols" = rownames(lung.sc), lung.sc)
row.names(lung.sc)=NULL

muscle.sc = as.matrix(sc.data[,grep("Muscle", colnames(sc.data))])
muscle.sc = data.frame("sc.exprs"=rowMeans(muscle.sc))
muscle.sc = cbind("symbols" = rownames(muscle.sc), muscle.sc)
row.names(muscle.sc)=NULL

##Get average MA data for each organ
path = selectDirectory(caption = "Select Directory Containing Organ-Specific Folders with Cel Files", label = "Select", path = NULL)
folder.names = dir(path)
annot = as.data.frame(mogene10sttranscriptclusterSYMBOL)
for (i in 1:length(folder.names)){
  setwd(file.path(path, paste(folder.names[i])))
  cels  = list.celfiles(full.names=TRUE)
  gfs = read.celfiles(cels)
  gfs.rma = rma(gfs, target = "core")
  gfs.rma.avg = data.frame(exprs = rowMeans(gfs.rma@assayData$exprs))
  gfs.rma.avg.merge = data.frame(exprs = merge(annot, gfs.rma.avg, by.x=1, by.y=0))
  assign(paste(folder.names[i],".ma", sep=""), data.frame(ma.exprs = gfs.rma.avg.merge[,3], symbols = gfs.rma.avg.merge[,2], row.names = NULL)) 
  
}

#get rid of duplicates
brain.ma = brain.ma[!duplicated(brain.ma$symbols),]
heart.ma = heart.ma[!duplicated(heart.ma$symbols),]
kidney.ma = kidney.ma[!duplicated(kidney.ma$symbols),]
liver.ma = liver.ma[!duplicated(liver.ma$symbols),]
lung.ma = lung.ma[!duplicated(lung.ma$symbols),]
muscle.ma = muscle.ma[!duplicated(muscle.ma$symbols),]

##Subset data by top DEGs
brain.ma.top=brain.ma[brain.ma$symbols %in% brain.top.genes,]
heart.ma.top=heart.ma[heart.ma$symbols %in% heart.top.genes,]
kidney.ma.top=kidney.ma[kidney.ma$symbols %in% kidney.top.genes,]
liver.ma.top=liver.ma[liver.ma$symbols %in% liver.top.genes,]
lung.ma.top=lung.ma[lung.ma$symbols %in% lung.top.genes,]
muscle.ma.top=muscle.ma[muscle.ma$symbols %in% muscle.top.genes,]

#create data frames with sc and ma expression of top DEGs
brain.top = merge(brain.ma[brain.ma$symbols %in% brain.top.genes,], brain.sc[brain.sc$symbols %in% brain.top.genes,], by = "symbols")
heart.top = merge(heart.ma[heart.ma$symbols %in% heart.top.genes,], heart.sc[heart.sc$symbols %in% heart.top.genes,], by = "symbols")
kidney.top = merge(kidney.ma[kidney.ma$symbols %in% kidney.top.genes,], kidney.sc[kidney.sc$symbols %in% kidney.top.genes,], by = "symbols")
liver.top = merge(liver.ma[liver.ma$symbols %in% liver.top.genes,], liver.sc[liver.sc$symbols %in% liver.top.genes,], by = "symbols")
lung.top = merge(lung.ma[lung.ma$symbols %in% lung.top.genes,], lung.sc[lung.sc$symbols %in% lung.top.genes,], by = "symbols")
muscle.top = merge(muscle.ma[muscle.ma$symbols %in% muscle.top.genes,], muscle.sc[muscle.sc$symbols %in% muscle.top.genes,], by = "symbols")

#sort by MA expression 
brain.top = brain.top[order(-brain.top$ma.exprs),]
heart.top = heart.top[order(-heart.top$ma.exprs),]
kidney.top = kidney.top[order(-kidney.top$ma.exprs),]
liver.top = liver.top[order(-liver.top$ma.exprs),]
lung.top = lung.top[order(-lung.top$ma.exprs),]
muscle.top = muscle.top[order(-muscle.top$ma.exprs),]

#only take top 100 DEGs
brain.top = brain.top[1:100,]
heart.top = heart.top[1:100,]
kidney.top = kidney.top[1:100,]
liver.top = liver.top[1:100,]
lung.top = lung.top[1:100,]
muscle.top = muscle.top[1:100,]

#Generate scatterplots seen in Supplemental Figure 4A-F
brain.deg.cor = ggscatter(brain.top,  x = "sc.exprs", y = "ma.exprs", shape = 19, size =0.5, ylab = expression(paste("Microarray ( log" ["2"], " intensity)")), xlab = expression(paste("sc RNA Seq (log" ["2"], " counts)")))
brain.deg.cor + stat_smooth()

heart.deg.cor = ggscatter(heart.top,  x = "sc.exprs", y = "ma.exprs", shape = 19, size =0.5, ylab = expression(paste("Microarray ( log" ["2"], " intensity)")), xlab = expression(paste("sc RNA Seq (log" ["2"], " counts)")))
heart.deg.cor + stat_smooth()

kidney.deg.cor = ggscatter(kidney.top,  x = "sc.exprs", y = "ma.exprs", shape = 19, size =0.5, ylab = expression(paste("Microarray ( log" ["2"], " intensity)")), xlab = expression(paste("sc RNA Seq (log" ["2"], " counts)")))
kidney.deg.cor + stat_smooth()

liver.deg.cor = ggscatter(liver.top,  x = "sc.exprs", y = "ma.exprs", shape = 19, size =0.5, ylab = expression(paste("Microarray ( log" ["2"], " intensity)")), xlab = expression(paste("sc RNA Seq (log" ["2"], " counts)")))
liver.deg.cor + stat_smooth()

lung.deg.cor = ggscatter(lung.top,  x = "sc.exprs", y = "ma.exprs", shape = 19, size =0.5, ylab = expression(paste("Microarray ( log" ["2"], " intensity)")), xlab = expression(paste("sc RNA Seq (log" ["2"], " counts)")))
lung.deg.cor + stat_smooth()

muscle.deg.cor = ggscatter(muscle.top,  x = "sc.exprs", y = "ma.exprs", shape = 19, size =0.5, ylab = expression(paste("Microarray ( log" ["2"], " intensity)")), xlab = expression(paste("sc RNA Seq (log" ["2"], " counts)")))
muscle.deg.cor + stat_smooth()

##Clustering and heatmap
#put top DEGs into table
brain.top = as.data.frame(topTable(fit2.brain, number = 200, adjust = "bonferroni"))
heart.top = as.data.frame(topTable(fit2.heart, number = 200, adjust = "bonferroni"))
kidney.top = as.data.frame(topTable(fit2.kidney, number = 200, adjust = "bonferroni"))
liver.top = as.data.frame(topTable(fit2.liver, number = 200, adjust = "bonferroni"))
lung.top = as.data.frame(topTable(fit2.lung, number = 200, adjust = "bonferroni"))
muscle.top = as.data.frame(topTable(fit2.muscle, number = 200, adjust = "bonferroni"))

#get list of top positive DEGs 
brain.ma.degs = brain.top[which(brain.top$logFC > 0),2]
heart.ma.degs = heart.top[which(heart.top$logFC > 0),2]
kidney.ma.degs = kidney.top[which(kidney.top$logFC > 0),2]
liver.ma.degs = liver.top[which(liver.top$logFC > 0),2]
lung.ma.degs = lung.top[which(lung.top$logFC > 0),2]
muscle.ma.degs = muscle.top[which(muscle.top$logFC > 0),2]

#Get list of tissue-specific DEGs from sc RNA-Seq data
#load seurat object
file = tk_choose.files(caption = "Select Seurat Object")
EC = readRDS(file)

#subset Seurat object based on the tissues from the microarray data and male mice 
#(only males were used in the Rafii study)
EC.deg = subset(EC, idents = c("Brain", "Heart", "Kidney", "Liver", "Lung",
                                      "Skeletal Muscle"))
Idents(EC.deg) = "mouse.sex"
EC.deg = subset(EC.deg, idents = "M")
Idents(EC.deg) = "orig.ident"

#get DEGs from each organ
brain.sc.degs.all = FindMarkers(EC.deg, ident.1 = "Brain", only.pos = TRUE, test.use = "wilcox", logfc.threshold = log(2))
heart.sc.degs.all = FindMarkers(EC.deg, ident.1 = "Heart", only.pos = TRUE, test.use = "wilcox", logfc.threshold = log(2))
kidney.sc.degs.all = FindMarkers(EC.deg, ident.1 = "Kidney", only.pos = TRUE, test.use = "wilcox", logfc.threshold = log(2))
liver.sc.degs.all = FindMarkers(EC.deg, ident.1 = "Liver", only.pos = TRUE, test.use = "wilcox", logfc.threshold = log(2))
lung.sc.degs.all = FindMarkers(EC.deg, ident.1 = "Lung", only.pos = TRUE, test.use = "wilcox", logfc.threshold = log(2))
muscle.sc.degs.all = FindMarkers(EC.deg, ident.1 = "Limb_Muscle", only.pos = TRUE, test.use = "wilcox", logfc.threshold = log(2))

#put tissue-specific DEGs from sc RNA Seq into an ordered list
brain.sc.degs = data.frame(order.sc = 1:nrow(brain.sc.degs.all), symbols = rownames(brain.sc.degs.all))
heart.sc.degs = data.frame(order.sc = 1:nrow(heart.sc.degs.all), symbols = rownames(heart.sc.degs.all))
kidney.sc.degs = data.frame(order.sc = 1:nrow(kidney.sc.degs.all), symbols = rownames(kidney.sc.degs.all))
liver.sc.degs = data.frame(order.sc = 1:nrow(liver.sc.degs.all), symbols = rownames(liver.sc.degs.all))
lung.sc.degs = data.frame(order.sc = 1:nrow(lung.sc.degs.all), symbols = rownames(lung.sc.degs.all))
muscle.sc.degs = data.frame(order.sc = 1:nrow(muscle.sc.degs.all), symbols = rownames(muscle.sc.degs.all))

#Order the MA DEG lists
brain.ma.degs = data.frame(order.ma = 1:length(brain.ma.degs), symbols = brain.ma.degs)
heart.ma.degs = data.frame(order.ma = 1:length(heart.ma.degs), symbols = heart.ma.degs)
kidney.ma.degs = data.frame(order.ma = 1:length(kidney.ma.degs), symbols = kidney.ma.degs)
liver.ma.degs = data.frame(order.ma = 1:length(liver.ma.degs), symbols = liver.ma.degs)
lung.ma.degs = data.frame(order.ma = 1:length(lung.ma.degs), symbols = lung.ma.degs)
muscle.ma.degs = data.frame(order.ma = 1:length(muscle.ma.degs), symbols = muscle.ma.degs)

#Get intersection of MA and SC DEG lists. Order the list by summing the 
#SC and MA ranks and re-ordering
brain.degs = merge(x=brain.sc.degs, y=brain.ma.degs, by = "symbols")
brain.degs$order.total = brain.degs$order.ma + brain.degs$order.sc
brain.degs[,2:3] = NULL
brain.degs = brain.degs[order(brain.degs$order.total),]
brain.degs.top10 = as.character(brain.degs[1:10,1])

heart.degs = merge(x=heart.sc.degs, y=heart.ma.degs, by = "symbols")
heart.degs$order.total = heart.degs$order.ma + heart.degs$order.sc
heart.degs[,2:3] = NULL
heart.degs = heart.degs[order(heart.degs$order.total),]
heart.degs.top10 = as.character(heart.degs[1:10,1])

kidney.degs = merge(x=kidney.sc.degs, y=kidney.ma.degs, by = "symbols")
kidney.degs$order.total = kidney.degs$order.ma + kidney.degs$order.sc
kidney.degs[,2:3] = NULL
kidney.degs = kidney.degs[order(kidney.degs$order.total),]
kidney.degs.top10 = as.character(kidney.degs[1:10,1])

liver.degs = merge(x=liver.sc.degs, y=liver.ma.degs, by = "symbols")
liver.degs$order.total = liver.degs$order.ma + liver.degs$order.sc
liver.degs[,2:3] = NULL
liver.degs = liver.degs[order(liver.degs$order.total),]
liver.degs.top10 = as.character(liver.degs[1:10,1])

lung.degs = merge(x=lung.sc.degs, y=lung.ma.degs, by = "symbols")
lung.degs$order.total = lung.degs$order.ma + lung.degs$order.sc
lung.degs[,2:3] = NULL
lung.degs = lung.degs[order(lung.degs$order.total),]
lung.degs.top10 = as.character(lung.degs[1:10,1])

muscle.degs = merge(x=muscle.sc.degs, y=muscle.ma.degs, by = "symbols")
muscle.degs$order.total = muscle.degs$order.ma + muscle.degs$order.sc
muscle.degs[,2:3] = NULL
muscle.degs = muscle.degs[order(muscle.degs$order.total),]
muscle.degs.top10 = as.character(muscle.degs[1:10,1])

##Create matrix of ranked expression data from MA and SC
#get microarray expression data from ExpressionSet object
exprs.ma = Biobase::exprs(data.ma.final)
#convert the PROBEID rownames to a column and add to the expression matrix
exprs.ma.r2c = cbind(PROBEID = as.numeric(rownames(exprs.ma)), exprs.ma)
colnames(annot)[1] = "PROBEID"
exprs.ma.annot = merge(exprs.ma.r2c, annot, by = "PROBEID")
exprs.ma = exprs.ma.annot[c(20,2:19)]
exprs.ma = distinct(exprs.ma)

#get sc RNA expression data from Seurat object
exprs.sc = as.matrix(GetAssayData(EC.deg, slot = "data"))
exprs.sc = as.data.frame(exprs.sc)
exprs.sc = cbind(SYMBOL = rownames(exprs.sc), exprs.sc)

#merge ma and sc data frames by symbols and rank expression
colnames(exprs.ma)[1] = "SYMBOL"
exprs.all = merge(exprs.ma, exprs.sc, by = "SYMBOL")
exprs.all.ranked = cbind(exprs.all$SYMBOL, apply(exprs.all[,2:1481], 2, rank, ties.method = "average"))
colnames(exprs.all.ranked)[1] = "symbols"
exprs.all.ranked = data.frame(exprs.all.ranked, stringsAsFactors = FALSE)
exprs.all.ranked = exprs.all.ranked[!duplicated(exprs.all.ranked$symbols),]

#Cluster and heatmap
#get color palette
col.pal = RColorBrewer::brewer.pal(9, "Reds")

#set up empty vector to populate with tissue names
tissue.ids = vector(mode = "character", length = 1480)

#Find location of each tissue-specific sample. Add the tissue name to the corresponding location 
#in the empty vector.
samples = colnames(exprs.all.ranked)[2:1481]
brain.locs = grep("Brain", samples)
tissue.ids[brain.locs] = "Brain"
heart.locs = grep("Heart", samples)
tissue.ids[heart.locs] = "Heart"
kidney.locs = grep("Kidney", samples)
tissue.ids[kidney.locs] = "Kidney"
liver.locs = grep("Liver", samples)
tissue.ids[liver.locs] = "Liver"
lung.locs = grep("Lung", samples)
tissue.ids[lung.locs] = "Lung"
muscle.locs = grep("Muscle", samples)
tissue.ids[muscle.locs] = "Muscle"

#add tissue names to sample IDs
sample.tissue.ids = data.frame(Tissue = tissue.ids)
rownames(sample.tissue.ids) = samples

#form single list of DEGs from each tissue and create data frame to indicate the tissue they are from
degs.all = c(brain.degs.top10, heart.degs.top10, kidney.degs.top10, liver.degs.top10, lung.degs.top10, muscle.degs.top10)
degs.tissues = c(rep("Brain", times = 10), rep("Heart", times = 10), rep("Kidney", times = 10), rep("Liver", times = 10), rep("Lung", times = 10), rep("Muscle", times = 10))
degs.tissue.ids = data.frame(degs = degs.all, Tissue = degs.tissues)

#get rid of two genes because they are markers for 2 different tissues
degs.tissue.ids = degs.tissue.ids[-c(17,52),]
rownames(degs.tissue.ids) = degs.tissue.ids[,1]
degs.tissue.ids[,1] = NULL

#pull out expression for DEGs and convert to data frame
exprs.all.ranked.degs = data.frame(exprs.all.ranked[(exprs.all.ranked[,1] %in% degs.all),1:1481], stringsAsFactors = FALSE)
rnames = exprs.all.ranked.degs[,1]
rownames(exprs.all.ranked.degs) = rnames
exprs.all.ranked.degs[,1] = NULL

exprs.all.ranked.degs.mat = data.matrix(exprs.all.ranked[(exprs.all.ranked[,1] %in% degs.all),1:1481])
rownames(exprs.all.ranked.degs.mat) = rnames
exprs.all.ranked.degs.mat = exprs.all.ranked.degs.mat[,-1]

#sort dendrogram by the order that you inputted the data (as close as possible)
order.desired = c("Heart", "Muscle", "Brain", "Kidney", "Lung", "Liver")
sample.order = rownames(sample.tissue.ids[order(match(sample.tissue.ids$Tissue, order.desired)), , drop = FALSE])
sample.order=gsub(" ", ".", sample.order)
mat.cluster.cols = hclust(dist(t(exprs.all.ranked.degs.mat)))
mat.cluster.cols = rotate(mat.cluster.cols, sample.order)

degs.tissue.ids.ordered = degs.tissue.ids[order(match(degs.tissue.ids$Tissue, order.desired)), , drop = FALSE]
mat.cluster.rows = dendextend::rotate(hclust(dist(exprs.all.ranked.degs.mat)), rownames(degs.tissue.ids.ordered))

#change annotation order to get same order in legend
sample.tissue.ids.ordered = sample.tissue.ids[order(match(sample.tissue.ids$Tissue, order.desired)), , drop = FALSE]
degs.tissue.ids.ordered = degs.tissue.ids[order(match(degs.tissue.ids$Tissue, order.desired)), , drop = FALSE]

#set up vector of column names to show position of MA data
loc.mas = grep("CEL", mat.cluster.cols$labels)
mas = vector(mode = "character", length = length(mat.cluster.cols$labels))
mas[loc.mas] = mat.cluster.cols$labels[loc.mas]

#generate heatmap seen in Supplemental Figure 4B
col.pal = RColorBrewer::brewer.pal(6, "Set1")
annot.colors = list(Tissue = c(Heart = "#27b24a", Muscle = "#cc71ad", Brain = "#b59f31", Kidney = "#2cb789", Lung = "#24b1e6", Liver = "#19bcc2"))
pheatmap::pheatmap(
  mat = exprs.all.ranked.degs.mat, 
  cluster_rows = mat.cluster.rows, 
  cluster_cols = mat.cluster.cols, 
  annotation_col = sample.tissue.ids.ordered, 
  annotation_row = degs.tissue.ids.ordered, 
  annotation_legend = FALSE,
  annotation_colors = annot.colors, 
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  show_rownames = TRUE, 
  show_colnames = FALSE,
  fontsize_row = 5,
  fontsize_col = 4
)