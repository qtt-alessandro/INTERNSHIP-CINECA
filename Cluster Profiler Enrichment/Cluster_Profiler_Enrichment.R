lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)

library("rjson")


json_file_0 <- "Input_data/convol_filter0.json"
conv_filt_0 <- fromJSON(file=json_file_0)


json_file_2 <- "Input_data/dense.json"
dense <- fromJSON(file=json_file_2)


filter_0     <- conv_filt_0$`TCGA-HNSC`
filter_dense <- dense$`TCGA-HNSC`
selection = "HNSC"


library(org.Hs.eg.db)


gene_to_ID <- function(my.symbols){ 

  hs <- org.Hs.eg.db
  symbol_to_id <- select(hs, 
                         keys = my.symbols,
                         columns = c("ENTREZID", "SYMBOL"),
                         keytype = "SYMBOL")
  
  x2 <- na.omit(symbol_to_id$ENTREZID)
  return(x2)
}



filter_0_ID    <- gene_to_ID(filter_0)

dense_ID        <- gene_to_ID(filter_dense)

  
library(clusterProfiler)
library(ggplot2)

cellular_components <- function(input_genes,organism,class){
out_folder <- paste("CC_",class, sep="")
dir.create(out_folder)
setwd(out_folder)
ggo <- groupGO(gene = input_genes ,OrgDb=organism,ont="CC", level=3, readable=TRUE)
ego <- enrichGO(gene = input_genes,OrgDb =organism, ont="CC", pAdjustMethod="BH", pvalueCutoff  = 0.01,qvalueCutoff  = 0.05,readable=TRUE)
#gse <- gseGO(geneList = unstim,keyType="SYMBOL",OrgDb = 'org.Hs.eg.db', ont="CC", nPerm=1000,minGSSize=80,maxGSSize=500,pvalueCutoff = 0.05,verbose = FALSE)
write.csv(as.data.frame(ggo),"ggo_CC.csv")
write.csv(as.data.frame(ego),"ego_CC.csv")
#write.csv(as.data.frame(gse),"CC/gse_CC.csv")
png("barplot_CC_ggo.png",width = 780, height = 480)
plt1 <- barplot(ggo, drop=TRUE, showCategory=12)
print(plt1)
dev.off()
png("barplot_CC_ego.png",width = 780, height = 480)
plt2<-barplot(ego, showCategory=8)
print(plt2)
dev.off()
png("dotplot_CC.png",width = 780, height = 480)
plt3<-dotplot(ego)
print(plt3)
dev.off()
#png("CC/enrich_map_CC.png",width = 780, height = 480)
#emapplot(gse, showCategory = 10)
#dev.off()
png("cnet_plot_CC.png",width = 780, height = 480)
plt4<-cnetplot(ego, categorySize="pvalue")
print(plt4)
dev.off()
#png("CC/ridgeplot_CC.png",width = 780, height = 480)
#ridgeplot(gse) + labs(x = "enrichment distribution")
#dev.off()
}


#Biological process
biological_process <- function(input_genes,organism,class){
out_folder <- paste("BP_",class, sep="")
dir.create(out_folder)
setwd(out_folder)
ggo <- groupGO(gene = input_genes,OrgDb = organism,ont="BP", level=3, readable=TRUE)
ego <- enrichGO(gene = input_genes,OrgDb = organism, ont="BP", pAdjustMethod="BH", pvalueCutoff  = 0.01,qvalueCutoff  = 0.05,readable=TRUE)
#gse <- gseGO(geneList = input_genes,keyType="SYMBOL",OrgDb='org.Hs.eg.db', ont="BP", nPerm=1000,minGSSize=80,maxGSSize=500,pvalueCutoff = 0.05,verbose = FALSE)
write.csv(as.data.frame(ggo),"ggo_BP.csv")
write.csv(as.data.frame(ego),"ego_BP.csv")
#write.csv(as.data.frame(gse),"BP/gse_BP.csv")
png("barplot_BP_ggo.png",width = 780, height = 480)
plt1<-barplot(ggo, drop=TRUE, showCategory=12)
print(plt1)
dev.off()
png("barplot_BP_ego.png",width = 780, height = 480)
plt2<-barplot(ego, showCategory=8)
print(plt2)
dev.off()
png("dotplot_BP.png",width = 780, height = 480)
plt3<-dotplot(ego)
print(plt3)
dev.off()
#png("BP/enrich_map_BP.png",width = 780, height = 480)
#emapplot(gse, showCategory = 10)
#dev.off()
png("cnet_plot_BP.png",width = 780, height = 480)
plt4<-cnetplot(ego, categorySize="pvalue" ,circular = TRUE, colorEdge = TRUE)

print(plt4)
dev.off()
#png("BP/ridgeplot_BP.png",width = 780, height = 480)
#ridgeplot(gse) + labs(x = "enrichment distribution")
#dev.off()
}


#Molecular function
molecolar_function <-function(input_genes,organism,class){
out_folder <- paste("MF_",class, sep="")
dir.create(out_folder)
setwd(out_folder)
ggo <- groupGO(gene = input_genes,OrgDb = organism,ont="MF", level=3, readable=TRUE)
ego <- enrichGO(gene = input_genes,OrgDb = organism, ont="MF", pAdjustMethod="BH", pvalueCutoff  = 0.01,qvalueCutoff  = 0.05,readable=TRUE)
#gse <- gseGO(geneList = gene_list,keyType="SYMBOL",OrgDb=org.Dm.eg.db, ont="MF", nPerm=1000,minGSSize=80,maxGSSize=500,pvalueCutoff = 0.05,verbose = FALSE)
write.csv(as.data.frame(ggo),"ggo_MF.csv")
write.csv(as.data.frame(ego),"ego_MF.csv")
#write.csv(as.data.frame(gse),"MF/gse_MF.csv")
png("barplot_MF_ggo.png",width = 780, height = 480)
plt1<-barplot(ggo, drop=TRUE, showCategory=12)
print(plt1)
dev.off()
png("barplot_MF_ego.png",width = 780, height = 480)
plt2<-barplot(ego, showCategory=8)
print(plt2)
dev.off()
png("dotplot_MF.png",width = 780, height = 480)
plt3<-dotplot(ego)
print(plt3)
dev.off()
#png("MF/enrich_map_MF.png",width = 780, height = 480)
#emapplot(gse, showCategory = 10)
#dev.off()
png("cnnet_plot_MF.png",width = 780, height = 480)
plt4<-cnetplot(ego, categorySize="pvalue", circular = TRUE, colorEdge = TRUE)
print(plt4)
dev.off()
#png("MF/ridgeplot_MF.png",width = 780, height = 480)
#ridgeplot(gse) + labs(x = "enrichment distribution")
#dev.off()
}


out_name_filter0 <- paste("filter_0_",selection, sep="")
biological_process(input_genes=filter_0_ID,organism='org.Hs.eg.db',class=out_name_filter0)
molecolar_function(input_genes=filter_0_ID,organism='org.Hs.eg.db',class=out_name_filter0)
cellular_components(input_genes=filter_0_ID,organism='org.Hs.eg.db',class=out_name_filter0)


out_name_filter_dense <- paste("dense_",selection, sep="")
biological_process(input_genes=dense_ID,organism='org.Hs.eg.db',class=out_name_filter_dense)
molecolar_function(input_genes=dense_ID,organism='org.Hs.eg.db',class=out_name_filter_dense)
cellular_components(input_genes=dense_ID,organism='org.Hs.eg.db',class=out_name_filter_dense)

