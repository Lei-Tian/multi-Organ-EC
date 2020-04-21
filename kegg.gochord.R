###KEGG pathway analysis and GOChord plot generation 

#load packages
require(Seurat)
require(AnnotationDbi)
require(org.Mm.eg.db)
require(GeneAnswers)
require(GOplot)
require(KEGG.db)

#Load seurat object from file
seurat.object = readRDS(file)

#Find differentially expressed genes
markers = FindAllMarkers(seurat.object, logfc.threshold = 0.25, test.use = "wilcox",
                            min.pct = 0.1, only.pos = TRUE)

##Run KEGG enrichment analysis
#Get ENTREZ IDs of DEGs
genes.symbols = as.character(markers$gene)
genes.entrez = mapIds(org.Mm.eg.db, gene.symbols, 'ENTREZID', 'SYMBOL')
genes.entrez = genes.entrez[!is.na(genes.entrez)]

#run KEGG enrichment
kegg = geneAnswersBuilder(genes.entrez, 'org.Mm.eg.db',
                                categoryType = 'KEGG', testType = 'hyperG',
                                pvalueT=0.1, FDR.correction = TRUE)
kegg = geneAnswersReadable(genes.kegg, verbose = FALSE)

##GO chord plot
#set up matrix with rows as DEGs, first 5 columns as top 5 enriched KEGG pathways
#and 6th column as logFC
a = markers$avg_logFC
b = matrix(0, nrow = nrow(markers), ncol = 5)
c = cbind(b,a)
colnames(c)[1:5] = names(kegg@genesInCategory[1:5]) 
colnames(c)[6] = "logFC"
rownames(c) = markers$gene

#if gene i appears in the list of a KEGG pathway, put 1 in the column for
#that KEGG pathway
for (i in 1:nrow(c)){
  for (j in 1:5){
    if (rownames(c)[i] %in% kegg@genesInCategory[[j]]){
      c[i,names(kegg@genesInCategory)[[j]]]=1
      
    }
  }
}
c = c[rowSums(c[,(1:5)])>0,]

#generate GO chord plots
GOChord(c, space = 0.02, gene.order = 'logFC', gene.space = 0.2, gene.size = 5)



