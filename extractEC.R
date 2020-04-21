a=read.csv("cells.csv")

t="Aorta"
load(paste("../data/facs_",t,"_seurat_tiss.Robj",sep=""))
tiss@meta.data$orig.ident=rep(t,length(tiss@meta.data$orig.ident))
b=as.vector(a[which(a$tissue==t),"cell"])
EC=SubsetData(tiss, cells.use = b) 
rm(tiss)
EC@meta.data$subannotation=rep("EC",nrow(EC@meta.data))

t="Brain_Non-Myeloid"
load(paste("../data/facs_",t,"_seurat_tiss.Robj",sep=""))
tiss@meta.data$orig.ident=rep("Brain",length(tiss@meta.data$orig.ident))
b=as.vector(a[which(a$tissue=="Brain_Non-Myeloid"),"cell"])
tmp=SubsetData(tiss, cells.use = b) 
EC=MergeSeurat(object1=EC,object2=tmp,add.cell.id1="Aorta",add.cell.id2="Brain",do.normalize = FALSE)
rm(tiss)


t="Diaphragm"
load(paste("../data/facs_",t,"_seurat_tiss.Robj",sep=""))
tiss@meta.data$orig.ident=rep(t,length(tiss@meta.data$orig.ident))
b=as.vector(a[which(a$tissue==t),"cell"])
tmp=SubsetData(tiss, cells.use = b) 
tmp@meta.data$subannotation=rep("EC",nrow(tmp@meta.data))
EC=MergeSeurat(object1=EC,object2=tmp,add.cell.id2=t,do.normalize = FALSE)
rm(tiss)

t="Fat"
load(paste("../data/facs_",t,"_seurat_tiss.Robj",sep=""))
tiss@meta.data$orig.ident=rep(t,length(tiss@meta.data$orig.ident))
b=as.vector(a[which(a$tissue==t),"cell"])
tmp=SubsetData(tiss, cells.use = b) 
tmp@meta.data$subannotation=rep("EC",nrow(tmp@meta.data))
EC=MergeSeurat(object1=EC,object2=tmp,add.cell.id2=t,do.normalize = FALSE)
rm(tiss)

t="Heart"
load(paste("../data/facs_",t,"_seurat_tiss.Robj",sep=""))
tiss@meta.data$orig.ident=rep(t,length(tiss@meta.data$orig.ident))
b=as.vector(a[which(a$tissue==t),"cell"])
tmp=SubsetData(tiss, cells.use = b) 
tmp@meta.data$subannotation=rep("EC",nrow(tmp@meta.data))
EC=MergeSeurat(object1=EC,object2=tmp,add.cell.id2=t,do.normalize = FALSE)
rm(tiss)

t="Kidney"
load(paste("../data/facs_",t,"_seurat_tiss.Robj",sep=""))
tiss@meta.data$orig.ident=rep(t,length(tiss@meta.data$orig.ident))
b=as.vector(a[which(a$tissue==t),"cell"])
tmp=SubsetData(tiss, cells.use = b) 
tmp@meta.data$subannotation=rep("EC",nrow(tmp@meta.data))
EC=MergeSeurat(object1=EC,object2=tmp,add.cell.id2=t,do.normalize = FALSE)
rm(tiss)

t="Limb_Muscle"
load(paste("../data/facs_",t,"_seurat_tiss.Robj",sep=""))
tiss@meta.data$orig.ident=rep(t,length(tiss@meta.data$orig.ident))
b=as.vector(a[which(a$tissue==t),"cell"])
tmp=SubsetData(tiss, cells.use = b) 
tmp@meta.data$subannotation=rep("EC",nrow(tmp@meta.data))
EC=MergeSeurat(object1=EC,object2=tmp,add.cell.id2=t,do.normalize = FALSE)
rm(tiss)

t="Liver"
load(paste("../data/facs_",t,"_seurat_tiss.Robj",sep=""))
tiss@meta.data$orig.ident=rep(t,length(tiss@meta.data$orig.ident))
b=as.vector(a[which(a$tissue==t),"cell"])
tmp=SubsetData(tiss, cells.use = b) 
tmp@meta.data$subannotation=rep("EC",nrow(tmp@meta.data))
EC=MergeSeurat(object1=EC,object2=tmp,add.cell.id2=t,do.normalize = FALSE)
rm(tiss)

t="Lung"
load(paste("../data/facs_",t,"_seurat_tiss.Robj",sep=""))
tissue@meta.data$orig.ident=rep(t,length(tissue@meta.data$orig.ident))
b=as.vector(a[which(a$tissue==t),"cell"])
tmp=SubsetData(tissue, cells.use = b) 
tmp@meta.data$subannotation=rep("EC",nrow(tmp@meta.data))
EC=MergeSeurat(object1=EC,object2=tmp,add.cell.id2=t,do.normalize = FALSE)
rm(tissue)

t="Mammary_Gland"
load(paste("../data/facs_",t,"_seurat_tiss.Robj",sep=""))
tiss@meta.data$orig.ident=rep(t,length(tiss@meta.data$orig.ident))
b=as.vector(a[which(a$tissue==t),"cell"])
tmp=SubsetData(tiss, cells.use = b) 
tmp@meta.data$subannotation=rep("EC",nrow(tmp@meta.data))
EC=MergeSeurat(object1=EC,object2=tmp,add.cell.id2=t,do.normalize = FALSE)
rm(tiss)

t="Pancreas"
load(paste("../data/facs_",t,"_seurat_tiss.Robj",sep=""))
tiss@meta.data$orig.ident=rep(t,length(tiss@meta.data$orig.ident))
b=as.vector(a[which(a$tissue==t),"cell"])
tmp=SubsetData(tiss, cells.use = b) 
tmp@meta.data$subannotation=rep("EC",nrow(tmp@meta.data))
EC=MergeSeurat(object1=EC,object2=tmp,add.cell.id2=t,do.normalize = FALSE)
rm(tiss)

t="Trachea"
load(paste("../data/facs_",t,"_seurat_tiss.Robj",sep=""))
tiss@meta.data$orig.ident=rep(t,length(tiss@meta.data$orig.ident))
b=as.vector(a[which(a$tissue==t),"cell"])
tmp=SubsetData(tiss, cells.use = b) 
tmp@meta.data$subannotation=rep("EC",nrow(tmp@meta.data))
EC=MergeSeurat(object1=EC,object2=tmp,add.cell.id2=t,do.normalize = FALSE)
rm(tiss)

table(EC@meta.data$orig.ident)
save(EC,file="EC.Robj")

EC=FilterCells(EC,subset.names=c("nGene","nReads"),low.thresholds = c(500,50000))
EC=NormalizeData(EC,scale.factor=1e6)
EC=ScaleData(EC, vars.to.regress = c("nReads", "percent.ribo","percent.ercc"))
EC=FindVariableGenes(EC,do.plot=TRUE,x.high.cutoff = Inf,y.cutoff = 0.5)
EC=RunPCA(EC,pc.genes=EC@var.genes)
EC=FindClusters(EC,reduction.type = "pca",dims.use=1:20,force.recalc=TRUE)
EC=RunTSNE(EC,dims.use=1:20,check_duplicates=F)

pdf("EC.TSNE.pdf")
TSNEPlot(EC,pt.size=1.2)
TSNEPlot(EC,pt.size=1.2, group.by = "orig.ident")
TSNEPlot(EC,pt.size=1.2, group.by = "mouse.id")
TSNEPlot(EC,pt.size=1.2, group.by = "mouse.sex")
for(i in c(1:12)){
mycols=rep("gray",12)
mycols[i]="red"
t=names(table(EC@meta.data$orig.ident))[i]
mycells=EC@cell.names[which(EC@meta.data$orig.ident==t)]
TSNEPlot(EC,pt.size=1.2, group.by = "orig.ident",colors.use = mycols)
TSNEPlot(EC,pt.size=1.2, group.by = "mouse.id",cells.use = mycells)
TSNEPlot(EC,pt.size=1.2, group.by = "mouse.sex",cells.use = mycells)
}
dev.off()

lung=FilterCells(lung,subset.names=c("nGene","nReads"),low.thresholds = c(500,50000))
lung=NormalizeData(lung,scale.factor=1e6)
lung=ScaleData(lung, vars.to.regress = c("nReads", "percent.ribo","percent.ercc"))
lung=FindVariableGenes(lung,do.plot=TRUE,x.high.cutoff = Inf,y.cutoff = 0.5)
lung=RunPCA(lung,pc.genes=lung@var.genes)
lung=FindClusters(lung,reduction.type = "pca",dims.use=1:20,force.recalc=TRUE)
lung=RunTSNE(lung,dims.use=1:20,check_duplicates=F)
pdf("Lung.EC.TSNE.pdf")
TSNEPlot(lung,pt.size=1.2)
TSNEPlot(lung,pt.size=1.2, group.by = "mouse.id")
TSNEPlot(lung,pt.size=1.2, group.by = "mouse.sex")
dev.off()
a=FindAllMarkers(lung,logfc.threshold = log(2),min.pct=0.25,only.pos = TRUE)
write.table(a,"Lung.clusters.DEGs.txt",quote=F,sep="\t")

lung <- SetAllIdent(lung, id = "mouse.sex")
a=FindAllMarkers(lung,logfc.threshold = log(2),min.pct=0.25,only.pos = TRUE)
write.table(a,"Lung.genders.DEGs.txt",quote=F,sep="\t")



tsne=EC@dr$tsne@cell.embeddings
Cluster=EC@ident
tissue=EC@meta.data$orig.ident
results=cbind(tsne,Cluster,tissue)
write.table(results,"tsne.txt",quote=F,sep="\t")

save(EC,file="EC.TSNE.Robj")

pdf("cluster.tree.pdf")
BuildClusterTree(EC,do.reorder = TRUE,reorder.numeric = TRUE,pcs.use = 1:20)
dev.off()

a=FindAllMarkers(EC,logfc.threshold = log(2),min.pct=0.25,only.pos = TRUE)
write.table(a,"DEGs.txt",quote=F,sep="\t")

tissues=c("Aorta","Brain","Diaphragm","Fat","Heart","Kidney","Lung","Liver","Mammary_Gland","Limb_Muscle","Pancreas","Trachea")
tmp=tissue
for(i in 1:12){
  tmp[which(tmp==tissues[i])]=i
}
names(tmp)=names(EC@ident)
EC@ident=factor(tmp,ordered = TRUE, levels = c(1:12))
a=FindAllMarkers(EC,logfc.threshold = log(2),min.pct=0.25,only.pos = TRUE)
a$cluster=tissues[a$cluster]
write.table(a,"Tissues.DEGs.txt",quote=F,sep="\t")

library(dplyr)
a=a[which(a$p_val_adj<0.05),]
top10 <- a %>% group_by(cluster) %>% top_n(10, avg_logFC)
EC <- SetAllIdent(EC, id = "orig.ident")
pdf("DEG.heatmap.pdf",height=15,width=12)
DoHeatmap(object = EC, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE,rotate.key=TRUE,group.label.rot=TRUE)
dev.off()
write.table(top10,"top10.txt",quote=F,sep="\t",row.names = F)


load("../../../SHARED/Curated_Bioinformatic_data/00020_David_Paik_scRNA-seq_project_2016/leit/Processed data files/combine_alignment.Rda")
tmp=all@ident[which(all@ident==4)]
ecs=names(tmp)
pdf("DEG.heatmap.David.scECs.noScale.pdf",height=15,width=12)
DoHeatmap(object = all, genes.use = toupper(top10$gene),cells.use =ecs, slim.col.label = TRUE, remove.key = TRUE,use.scaled = FALSE)
dev.off()


library(dplyr)
load("../../../combine_alignment.Rda")
david.ec=SubsetData(all,ident.use=4)
david.ec <- SetAllIdent(david.ec, id = "orig.ident")
a=read.table("../Tissues.DEGs.txt")
a=a[which(a$p_val_adj<0.05),]
top10 <- a %>% group_by(cluster) %>% top_n(10, avg_logFC)
genes=toupper(top10$gene)[toupper(top10$gene) %in% row.names(david.ec@data)]
library(RColorBrewer)
RdBu = rev(brewer.pal(11, name="BrBG"))
pdf("david.heatmap.pdf",height = 15,width=9)
DoHeatmap(david.ec,genes.use=genes ,slim.col.label = TRUE, remove.key = TRUE, col.low =RdBu[1],col.mid =RdBu[6], col.high = RdBu[10])
dev.off()

load("G:/scRNA-Seq/Tabula Muris/data/facs_Lung_seurat_tiss.Robj")
load(paste0("G:/scRNA-Seq/Tabula Muris/data/facs_Heart_seurat_tiss.Robj"))
all=MergeSeurat(object1=tissue,object2=tiss,add.cell.id1="Lung",add.cell.id2="Heart",do.normalize = FALSE)
all@meta.data$subannotation=NULL
all@meta.data$annotation=NULL
for(t in c("Aorta","Bladder","Brain_Myeloid","Brain_Non-Myeloid","Diaphragm","Fat","Kidney",
           "Limb_Muscle","Liver","Mammary_Gland","Pancreas","Skin","Thymus","Spleen","Marrow",
           "Large_Intestine","Tongue","Trachea")){
  print(t)
  load(paste0("G:/scRNA-Seq/Tabula Muris/data/facs_",t,"_seurat_tiss.Robj"))
  tiss@meta.data$subannotation=NULL
  tiss@meta.data$annotation=NULL
  all=MergeSeurat(object1=all,object2=tiss,add.cell.id2=t,do.normalize = FALSE)
}

a=read.csv("../data/tabula-muris-master/00_data_ingest/18_global_annotation_csv/update/annotations_facs.csv")
a$cell=paste0(a$tissue,"_",a$cell)
row.names(a)=a$cell
a$cell_ontology_class=paste0(a$tissue,"_",a$cell_ontology_class)

all_qc=SubsetData(all, cells.use = a$cell)
all_qc@meta.data$cell_type=a[row.names(all_qc@meta.data),"cell_ontology_class"]
all_qc@ident=as.factor(all_qc@meta.data$cell_type)
names(all_qc@ident)=row.names(all_qc@meta.data)

all_qc=FilterCells(all_qc,subset.names=c("nGene","nReads"),low.thresholds = c(500,50000))
all_qc=NormalizeData(all_qc,scale.factor=1e6)
all_qc=ScaleData(all_qc, vars.to.regress = c("nReads", "percent.ribo","percent.ercc"))

pdf("Kcna5.Tmem100.violine.pdf",width=120,heigh=30)
VlnPlot(object = all_qc, features.plot = c(genes[1]),group.by="cell_type",point.size.use = 0,do.sort = TRUE,x.lab.rot=T,remove.legend = T)+
        theme(axis.title.x = element_blank(),axis.text.y  = element_text( vjust=0.5, size=34,face="bold"),
              axis.text.x  = element_text( vjust=0.5, size=34,face="bold",angle = 90),
              plot.title = element_text(lineheight=.8, size=40, face="bold.italic"))
VlnPlot(object = all_qc, features.plot = c(genes[2]),group.by="cell_type",point.size.use = 0,do.sort = TRUE,x.lab.rot=T,remove.legend = T)+
        theme(axis.title.x = element_blank(),axis.text.y  = element_text( vjust=0.5, size=34,face="bold"),
              axis.text.x  = element_text( vjust=0.5, size=34,face="bold",angle = 90),
              plot.title = element_text(lineheight=.8, size=40, face="bold.italic"))
dev.off()

results=c()
for(t in c("Aorta","Brain_Non-Myeloid","Diaphragm","Fat","Heart","Kidney","Limb_Muscle","Mammary_Gland","Pancreas","Trachea")){
  print(t)
  ec=paste0(t,"_endothelial cell")
  ec.markers <- FindMarkers(object = all_qc, ident.1 = ec,logfc.threshold =log(2), min.pct = 0.25, only.pos = TRUE)
  ec.markers=ec.markers[which(ec.markers$p_val_adj<0.05),]
  ec.markers$tissue=t
  ec.markers$gene=row.names(ec.markers)
  results=rbind(results,ec.markers)
}

t="Liver"
ec="Liver_endothelial cell of hepatic sinusoid"
ec.markers <- FindMarkers(object = all_qc, ident.1 = ec,logfc.threshold =log(2), min.pct = 0.25, only.pos = TRUE)
ec.markers=ec.markers[which(ec.markers$p_val_adj<0.05),]
ec.markers$tissue=t
ec.markers$gene=row.names(ec.markers)
results=rbind(results,ec.markers)

t="Lung"
ec="Lung_lung endothelial cell"
ec.markers <- FindMarkers(object = all_qc, ident.1 = ec,logfc.threshold =log(2), min.pct = 0.25, only.pos = TRUE)
ec.markers=ec.markers[which(ec.markers$p_val_adj<0.05),]
ec.markers$tissue=t
ec.markers$gene=row.names(ec.markers)
results=rbind(results,ec.markers)

write.table(results,"EC.vs.others.txt",quote=F,sep="\t")
