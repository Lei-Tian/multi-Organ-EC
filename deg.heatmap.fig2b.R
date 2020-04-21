###Generate heatmap showing expression of tissue-specific EC markers in ECs
###from different tissues

#Load packages
require(Seurat)
require(pheatmap)

#Load seurat object
ec = readRDS(file)

#Find markers (differentially expressed genes) in tissue-specific ECs 
ec.markers = FindAllMarkers(ec, logfc.threshold = 0.25, test.use = "wilcox",
                            min.pct = 0.1, only.pos = TRUE)

#Get DEGs from each group of tissue-specific ECs, sort by logFC, and only
#keep DEGs with adjusted p-value <0.05
adipose.top = ec.markers[(ec.markers$cluster %in% "Adipose Tissue"),]
adipose.top = adipose.top[order(adipose.top$avg_logFC, decreasing = TRUE),]
adipose.top= adipose.top[which(adipose.top$p_val_adj < 0.05),]

aorta.top = ec.markers[(ec.markers$cluster %in% "Aorta"),]
aorta.top = aorta.top[order(aorta.top$avg_logFC, decreasing = TRUE),]
aorta.top= aorta.top[which(aorta.top$p_val_adj < 0.05),]

brain.top = ec.markers[(ec.markers$cluster %in% "Brain"),]
brain.top = brain.top[order(brain.top$avg_logFC, decreasing = TRUE),]
brain.top= brain.top[which(brain.top$p_val_adj < 0.05),]

diaphragm.top = ec.markers[(ec.markers$cluster %in% "Diaphragm"),]
diaphragm.top = diaphragm.top[order(diaphragm.top$avg_logFC, decreasing = TRUE),]
diaphragm.top= diaphragm.top[which(diaphragm.top$p_val_adj < 0.05),]

heart.top = ec.markers[(ec.markers$cluster %in% "Heart"),]
heart.top = heart.top[order(heart.top$avg_logFC, decreasing = TRUE),]
heart.top= heart.top[which(heart.top$p_val_adj < 0.05),]

kidney.top = ec.markers[(ec.markers$cluster %in% "Kidney"),]
kidney.top = kidney.top[order(kidney.top$avg_logFC, decreasing = TRUE),]
kidney.top= kidney.top[which(kidney.top$p_val_adj < 0.05),]

liver.top = ec.markers[(ec.markers$cluster %in% "Liver"),]
liver.top = liver.top[order(liver.top$avg_logFC, decreasing = TRUE),]
liver.top= liver.top[which(liver.top$p_val_adj < 0.05),]

lung.top = ec.markers[(ec.markers$cluster %in% "Lung"),]
lung.top = lung.top[order(lung.top$avg_logFC, decreasing = TRUE),]
lung.top= lung.top[which(lung.top$p_val_adj < 0.05),]

mammary.gland.top = ec.markers[(ec.markers$cluster %in% "Mammary Gland"),]
mammary.gland.top = mammary.gland.top[order(mammary.gland.top$avg_logFC, decreasing = TRUE),]
mammary.gland.top= mammary.gland.top[which(mammary.gland.top$p_val_adj < 0.05),]

pancreas.top = ec.markers[(ec.markers$cluster %in% "Pancreas"),]
pancreas.top = pancreas.top[order(pancreas.top$avg_logFC, decreasing = TRUE),]
pancreas.top= pancreas.top[which(pancreas.top$p_val_adj < 0.05),]

muscle.top = ec.markers[(ec.markers$cluster %in% "Skeletal Muscle"),]
muscle.top = muscle.top[order(muscle.top$avg_logFC, decreasing = TRUE),]
muscle.top= muscle.top[which(muscle.top$p_val_adj < 0.05),]

trachea.top = ec.markers[(ec.markers$cluster %in% "Trachea"),]
trachea.top = trachea.top[order(trachea.top$avg_logFC, decreasing = TRUE),]
trachea.top= trachea.top[which(trachea.top$p_val_adj < 0.05),]

#Get the top 10 DEGs from each group of tissue-specific ECs
top10 = c(adipose.top$gene[1:10],aorta.top$gene[1:10],brain.top$gene[1:10],
          diaphragm.top$gene[1:10],heart.top$gene[1:10],kidney.top$gene[1:10], 
          liver.top$gene[1:10],lung.top$gene[1:10],mammary.gland.top$gene[1:10], 
          pancreas.top$gene[1:10],muscle.top$gene[1:10],trachea.top$gene[1:10])

exprs.mat = GetAssayData(ec, slot = "data")

#re-order expression matrix
loc.adipose = grep("Fat", colnames(exprs.mat))
loc.aorta = grep("Aorta", colnames(exprs.mat))
loc.brain = grep("Brain", colnames(exprs.mat))
loc.diaphragm = grep("Diaphragm", colnames(exprs.mat))
loc.heart = grep("Heart", colnames(exprs.mat))
loc.kidney = grep("Kidney", colnames(exprs.mat))
loc.liver = grep("Liver", colnames(exprs.mat))
loc.lung = grep("Lung", colnames(exprs.mat))
loc.mammary.gland = grep("Mammary", colnames(exprs.mat))
loc.pancreas = grep("Pancreas", colnames(exprs.mat))
loc.muscle = grep("Muscle", colnames(exprs.mat))
loc.trachea = grep("Trachea", colnames(exprs.mat))

 exprs.mat.ordered = exprs.mat[,c(loc.adipose, loc.aorta, loc.brain, loc.diaphragm, loc.heart,
                                   loc.kidney, loc.liver, loc.lung, loc.mammary.gland, loc.pancreas,
                                   loc.muscle,loc.trachea)]

#only get expression from DEGs to generate expression matrix for heatmap
exprs.adipose.degs = exprs.mat.ordered[(rownames(exprs.mat) %in% adipose.top$gene[1:10]),,drop=FALSE]
exprs.adipose.degs = exprs.adipose.degs[order(match(rownames(exprs.adipose.degs), 
                                                    adipose.top$gene[1:10])),]

exprs.aorta.degs = exprs.mat.ordered[(rownames(exprs.mat) %in% aorta.top$gene[1:10]),,drop=FALSE]
exprs.aorta.degs = exprs.aorta.degs[order(match(rownames(exprs.aorta.degs), 
                                                    aorta.top$gene[1:10])),]

exprs.brain.degs = exprs.mat.ordered[(rownames(exprs.mat) %in% brain.top$gene[1:10]),,drop=FALSE]
exprs.brain.degs = exprs.brain.degs[order(match(rownames(exprs.brain.degs), 
                                                    brain.top$gene[1:10])),]

exprs.diaphragm.degs = exprs.mat.ordered[(rownames(exprs.mat) %in% diaphragm.top$gene[1:10]),,drop=FALSE]
exprs.diaphragm.degs = exprs.diaphragm.degs[order(match(rownames(exprs.diaphragm.degs), 
                                                    diaphragm.top$gene[1:10])),]

exprs.heart.degs = exprs.mat.ordered[(rownames(exprs.mat) %in% heart.top$gene[1:10]),,drop=FALSE]
exprs.heart.degs = exprs.heart.degs[order(match(rownames(exprs.heart.degs), 
                                                    heart.top$gene[1:10])),]

exprs.kidney.degs = exprs.mat.ordered[(rownames(exprs.mat) %in% kidney.top$gene[1:10]),,drop=FALSE]
exprs.kidney.degs = exprs.kidney.degs[order(match(rownames(exprs.kidney.degs), 
                                                    kidney.top$gene[1:10])),]

exprs.liver.degs = exprs.mat.ordered[(rownames(exprs.mat) %in% liver.top$gene[1:10]),,drop=FALSE]
exprs.liver.degs = exprs.liver.degs[order(match(rownames(exprs.liver.degs), 
                                                    liver.top$gene[1:10])),]

exprs.lung.degs = exprs.mat.ordered[(rownames(exprs.mat) %in% lung.top$gene[1:10]),,drop=FALSE]
exprs.lung.degs = exprs.lung.degs[order(match(rownames(exprs.lung.degs), 
                                                    lung.top$gene[1:10])),]

exprs.mammary.gland.degs = exprs.mat.ordered[(rownames(exprs.mat) %in% mammary.gland.top$gene[1:10]),,drop=FALSE]
exprs.mammary.gland.degs = exprs.mammary.gland.degs[order(match(rownames(exprs.mammary.gland.degs), 
                                                    mammary.gland.top$gene[1:10])),]

exprs.pancreas.degs = exprs.mat.ordered[(rownames(exprs.mat) %in% pancreas.top$gene[1:10]),,drop=FALSE]
exprs.pancreas.degs = exprs.pancreas.degs[order(match(rownames(exprs.pancreas.degs), 
                                                    pancreas.top$gene[1:10])),]

exprs.muscle.degs = exprs.mat.ordered[(rownames(exprs.mat) %in% muscle.top$gene[1:10]),,drop=FALSE]
exprs.muscle.degs = exprs.muscle.degs[order(match(rownames(exprs.muscle.degs), 
                                                    muscle.top$gene[1:10])),]

exprs.trachea.degs = exprs.mat.ordered[(rownames(exprs.mat) %in% trachea.top$gene[1:10]),,drop=FALSE]
exprs.trachea.degs = exprs.trachea.degs[order(match(rownames(exprs.trachea.degs), 
                                                    trachea.top$gene[1:10])),]

exprs.degs = rbind(exprs.adipose.degs, exprs.aorta.degs, exprs.brain.degs,
                   exprs.diaphragm.degs, exprs.heart.degs, exprs.kidney.degs,
                   exprs.liver.degs, exprs.lung.degs, exprs.mammary.gland.degs,
                   exprs.pancreas.degs, exprs.muscle.degs, exprs.trachea.degs)

#generate heatmap
heatmap = pheatmap(
  mat = exprs.degs, 
  cluster_rows = FALSE, 
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 5,
  gaps_row=c(10,20,30,40,50,60,70,80,90,100,110),
  gaps_col = c(630,805,1454,1532,2694,2820,3002,3692,3739,3789,3928))

  
