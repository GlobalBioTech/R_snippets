#安装相关包：
if ( FALSE) {
  install.packages("digest")
  install.packages("stringi")
  source("https://bioconductor.org/biocLite.R")
  biocLite("clusterProfiler")  #富集
  biocLite("topGO")          #画图
  biocLite("Rgraphviz")       #供topGO调用
  biocLite("pathview")         #kegg pathway
  biocLite("DOSE")
  biocLite("stringi")
  biocLite("pathview")
  biocLite("org.Hs.eg.db")  #人的数据库
  biocLite("org.Rn.eg.db")  #大鼠的数据库
  biocLite("org.Mm.eg.db")  #小鼠的数据库
}

#第一步：加载相关包
{
  library(DOSE)
  library(GO.db)
  library(topGO)
  library(GSEABase)
  library(clusterProfiler)
  library(pathview)
  
}
##########################################################################################
#第二步：加载物种信息
#GO富集的已有19个物种查看网址：http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
#GO富集其他物种需要自建物种资料；
#KEGG富集物种表可在http://www.genome.jp/kegg/catalog/org_list.html  上查看
#物种；人
{
  library(org.Hs.eg.db)
  species_db_name<- org.Hs.eg.db
  species_kegg_name<- "hsa"
  
  }
#物种；小鼠
{
  library(org.Mm.eg.db)
  species_db_name<- org.Mm.eg.db
  species_kegg_name<- "mmu"
  
  }
#物种；大鼠
{
  library(org.Rn.eg.db)
  species_db_name<- org.Rn.eg.db
  species_kegg_name<- "rno" 
  }

#针对不常见物种，以水稻为例：
{
  library(AnnotationHub)
  hub<-AnnotationHub()
  query(hub,"Oryza sativa")    #以水稻为例
  将得到以下数据：
  #AnnotationHub with 2 records
  # snapshotDate(): 2017-10-27 
  # $dataprovider: Inparanoid8, ftp://ftp.ncbi.nlm...
  # $species: Oryza sativa, Oryza sativa_Japonica_...
  # $rdataclass: Inparanoid8Db, OrgDb
  # additional mcols(): taxonomyid, genome,
  #   description, coordinate_1_based,
  #   maintainer, rdatadateadded, preparerclass,
  #   tags, rdatapath, sourceurl, sourcetype 
  # retrieve records with, e.g.,
  #   'object[["AH10561"]]' 
              title                                  
    AH10561 | hom.Oryza_sativa.inp8.sqlite           
    AH59059 | org.Oryza_sativa_Japonica_Group.eg.s...
  # 其中AH59059就是我们需要得到水稻数据库：
    species_db_name<-hub[["AH59059"]]
}
##########################################################################################
#第三步：数据预处理
{
  first_name<- "ENSEMBL"
}
{
  first_name<- "ENTREZID"
}
{
  first_name<- "SYMBOL"
}

#针对deseq2的标准结果文件，有倍数信息：
{
data<- read.csv(file="diff-mRNA-deseq2_result.csv",header = T)
data<- data[match(ENSMBL,data$X),]
gene_symbol <- as.character(data[,1])
keytypes(species_db_name)
gene_entrez <- mapIds(x=species_db_name,
                      keys=gene_symbol,
                      keytype = first_name,
                      column="ENTREZID")
gene_fold<- as.character(data$log2FoldChange)
genelist<- cbind(gene_entrez,gene_fold)
genelist<- as.data.frame(genelist)    
genelist<- genelist[!is.na(genelist[,1]),]
genelist<- genelist[-2,]
genelist<-genelist[!duplicated(genelist$gene_entrez),]
genelist<- na.omit(genelist)
genename <- genelist$gene_entrez
gene_entrez<- as.character(genename)
genefoldchange<-genelist$gene_fold
genefoldchange<-as.numeric((as.character(genefoldchange)))
genelist2<-as.data.frame(genefoldchange)
gene_444<- na.omit(gene_entrez)
rownames(genelist2)<- gene_444
gene_symbol<- mapIds(x=species_db_name,keys=gene_entrez,keytype = "ENTREZID",column="SYMBOL")  #GO分析时，用symbol更方便
}

#如果是非deseq2标准输入，并且不去匹配deseq2的差异倍数信息；
{
data<- read.csv(file="mus-deseq2.csv",header=T)
gene_symbol <- as.character(data$mRNA)
gene_entrez<- gene_symbol
gene_symbol<- mapIds(x=species_db_name,keys=gene_entrez,keytype = "ENTREZID",column="SYMBOL")  #GO分析时，用symbol更方便
gene_entrez <- mapIds(x=species_db_name,
                      keys=gene_symbol,
                      keytype = first_name,
                      column="ENTREZID")
gene_entrez<-na.omit(gene_entrez)
}
#如果是非deseq2标准输入，但是去匹配deseq2的差异倍数信息；
{
data1<- read.csv(file="mRNA-deseq2_result.csv",header=T)
ENSMBL<- data1$mRNA
ENSMBL<- ENSMBL[!duplicated(ENSMBL)]
data<- read.csv(file="mRNA-deseq2_result.csv",header = T)
data<- data[data$padj<0.05,]
data<- data[match(ENSMBL,data$id),]
gene_symbol <- as.character(data$id)
keytypes(species_db_name)
gene_entrez <- mapIds(x=species_db_name,
                      keys=gene_symbol,
                      keytype = first_name,
                      column="ENTREZID")
gene_fold<- as.character(data$log2FoldChange)
genelist<- cbind(gene_entrez,gene_fold)
genelist<- as.data.frame(genelist)    
genelist<- genelist[!is.na(genelist[,1]),]
genelist<- genelist[-2,]
genelist<-genelist[!duplicated(genelist$gene_entrez),]
genelist<- na.omit(genelist)
genename <- genelist$gene_entrez
gene_entrez<- as.character(genename)
genefoldchange<-genelist$gene_fold
genefoldchange<-as.numeric((as.character(genefoldchange)))
genelist2<-as.data.frame(genefoldchange)
gene_444<- na.omit(gene_entrez)
rownames(genelist2)<- gene_444
gene_symbol<- mapIds(x=species_db_name,keys=gene_entrez,keytype = "ENTREZID",column="SYMBOL")  #GO分析时，用symbol更方
}
#针对ANOVA分析的结果，进行分析：
{
  data<- read.csv(file="anova_test_A-B-filter.csv",header = T)
  gene_name <- as.character(data[,1])  #默认第一列为基因名
  gene_entrez <- mapIds(x=species_db_name,
                        keys=gene_name,
                        keytype = first_name,
                        column="ENTREZID")
  gene_fold<- as.character(data[,4])    #第4列为差异倍数
  genelist<- cbind(gene_entrez,gene_fold)
  genelist<- as.data.frame(genelist)    
  genelist<- genelist[!is.na(genelist[,1]),]
  genelist<- genelist[-2,]
  genelist<-genelist[!duplicated(genelist$gene_entrez),]
  genelist<- na.omit(genelist)
  genename <- genelist$gene_entrez
  gene_entrez<- as.character(genename)
  genefoldchange<-genelist$gene_fold
  genefoldchange<-as.numeric((as.character(genefoldchange)))
  genelist2<-as.data.frame(genefoldchange)
  gene_444<- na.omit(gene_entrez)
  rownames(genelist2)<- gene_444
  gene_symbol<- mapIds(x=species_db_name,keys=gene_entrez,keytype = "ENTREZID",column="SYMBOL")  #GO分析时，用symbol更方便
}

##########################################################################################
#第四步,GO分析：
#ALL GO

erich.go.ALL<- enrichGO(gene=gene_symbol,OrgDb = species_db_name,keyType ="SYMBOL",ont="ALL" ,qvalueCutoff = 0.05 )
write.csv(as.data.frame(erich.go.ALL),"GO_ALL.csv",row.names=F )           	#输出表格，只需要输出ALL GO的即可

#GO BP

erich.go.BP<- enrichGO(gene=gene_symbol,OrgDb = species_db_name,keyType ="SYMBOL",ont="BP" ,qvalueCutoff = 0.05 )
pdf("GO-BP-barplot.pdf",width = 12)
barplot(erich.go.BP,color="qvalue")    #条形图
dev.off()
pdf("GO-BP-dotplot.pdf",width = 12)
dotplot(erich.go.BP,color="qvalue")   #点状图
dev.off()
pdf("GOgraph-BP.pdf")
plotGOgraph(erich.go.BP)  #有向无环图
dev.off()

#GO CC

erich.go.CC<- enrichGO(gene=gene_symbol,OrgDb = species_db_name,keyType ="SYMBOL",ont="CC" ,qvalueCutoff = 0.05 )
pdf("GO-CC-barplot.pdf",width = 12)
barplot(erich.go.CC,color="qvalue")    #条形图
dev.off()
pdf("GO-CC-dotplot.pdf",width = 12)
dotplot(erich.go.CC,color="qvalue")   #点状图
dev.off()
pdf("GOgraph-CC.pdf")
plotGOgraph(erich.go.CC)  #有向无环图
dev.off()

#GO MF

erich.go.MF<- enrichGO(gene=gene_symbol,OrgDb = species_db_name,keyType ="SYMBOL",ont="MF" ,qvalueCutoff = 0.05 )
pdf("GO-MF-barplot.pdf",width = 12)
barplot(erich.go.MF)    #条形图
dev.off()
pdf("GO-MF-dotplot.pdf",width = 12)
dotplot(erich.go.MF)   #点状图
dev.off()
pdf("GOgraph-MF.pdf")
plotGOgraph(erich.go.MF)  #有向无环图
dev.off()

---------------------------------
#综合数据画三个GO类型的柱状图
#erich.go.BP<- enrichGO(gene=gene_symbol,OrgDb = species_db_name,keyType ="SYMBOL",ont="BP" ,qvalueCutoff = 0.05 )
#收集前10行条目：
if (  nrow(as.data.frame(erich.go.BP)) > 10 ){
ego_result_BP <- na.omit(as.data.frame(erich.go.BP))[1:10, ]
BP_number<- 10
} else {
  ego_result_BP <- na.omit(as.data.frame(erich.go.BP)[1:nrow(as.data.frame(erich.go.BP)), ])
  BP_number<- nrow(as.data.frame(erich.go.BP))
  }

#erich.go.CC<- enrichGO(gene=gene_symbol,OrgDb = species_db_name,keyType ="SYMBOL",ont="CC" ,qvalueCutoff = 0.05 )
#收集前10行条目：
if (  nrow(as.data.frame(erich.go.CC)) > 10 ){
  ego_result_CC <- na.omit(as.data.frame(erich.go.CC))[1:10, ]
  CC_number<- 10
} else {
  ego_result_CC <- na.omit(as.data.frame(erich.go.CC)[1:nrow(as.data.frame(erich.go.CC)), ])
  CC_number<- nrow(as.data.frame(erich.go.CC))
  }

#erich.go.MF<- enrichGO(gene=gene_symbol,OrgDb = species_db_name,keyType ="SYMBOL",ont="MF" ,qvalueCutoff = 0.05 )
#收集前10行条目：
if (  nrow(as.data.frame(erich.go.MF)) > 10 ){
  ego_result_MF <- na.omit(as.data.frame(erich.go.MF))[1:10, ]
  MF_number<- 10
} else {
  ego_result_MF <- na.omit(as.data.frame(erich.go.MF)[1:nrow(as.data.frame(erich.go.MF)), ])
  MF_number<- nrow(as.data.frame(erich.go.MF))
  }


go_enrich_df <- data.frame(ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
                           Description=c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description),
                           GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
                           type=factor(c(rep("biological process", BP_number), rep("cellular component", CC_number),
                                         rep("molecular function", MF_number)), levels=c("molecular function", "cellular component", "biological process")))

## numbers as data on x axis
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
## shorten the names of GO terms
shorten_names <- function(x, n_word=4, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
  {
    if (nchar(x) > 40) x <- substr(x, 1, 40)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  } 
  else
  {
    return(x)
  }
}

labels=(sapply(
  levels(go_enrich_df$Description)[as.numeric(go_enrich_df$Description)],
  shorten_names))
names(labels) = rev(1:nrow(go_enrich_df))

## colors for bar // green, blue, orange
CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5")
library(ggplot2)
p <- ggplot(data=go_enrich_df, aes(x=number, y=GeneNumber, fill=type)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() +   #coord_flip为对图像进行翻转
  scale_fill_manual(values = CPCOLS) + theme_bw() + 
  scale_x_discrete(labels=labels) +
  xlab("GO term") + 
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "The Most Enriched GO Terms")

pdf("GO_enrichment.pdf")
p
dev.off()


##########################################################################################
#第五步,KEGG分析：

erich_kegg <- enrichKEGG(gene = gene_entrez,organism  = species_kegg_name ,qvalueCutoff = 1)  #做KEGG分析
pdf("KEGG-barplot.pdf",width = 12)
barplot(erich_kegg,color="qvalue")  ##x轴为基因counts数，y为信号通路，颜色为p值
dev.off()	
pdf("KEGG-dotplot.pdf",width = 12)
dotplot(erich_kegg, showCategory=12,color="qvalue")     	 ##x为geneRatio数
dev.off()	
write.csv(as.data.frame(erich_kegg),"KEGG.csv",row.names=F )  


#生成pathview图，针对有差异倍数的：
{
genelist4<- as.data.frame(genelist2)
genelist4$gene<- rownames(genelist2)
kegg_list<- cbind(erich_kegg$ID,erich_kegg$geneID)
ke<-erich_kegg$ID
for (i in ke) {
  print(i)
  kegg_single_gene<-strsplit(kegg_list[match(i,kegg_list[,1]),2],"/")
  genelist3<- genelist4[match(kegg_single_gene[[1]],genelist4[,2]),]
  write.csv(genelist3,paste(i,'.csv',sep=''),row.names = F)
  x.inv <- try(pathview(gene.data = genelist2,pathway.id = i,species = species_kegg_name,limit=list(gene=max(abs(genefoldchange)),cpd=1),kegg.native=TRUE), silent=TRUE)
  if ('try-error' %in% class(x.inv)) next
  pathview(gene.data = genelist2,pathway.id = i,species = species_kegg_name,limit=list(gene=max(abs(genelist3$genefoldchange)),cpd=1),kegg.native=TRUE)  
}
}

#生成pathview图，针对没有差异倍数的输入：
{
ke<-erich_kegg$ID
for (i in ke) {
  print(i)
  x.inv <- try(pathview(gene.data = gene_entrez,pathway.id = i,species = "hsa",kegg.native=TRUE), silent=TRUE)
  if ('try-error' %in% class(x.inv)) next
  pathview(gene.data = gene_entrez,pathway.id = i,species = "hsa",kegg.native=TRUE)  ###将部分基因ID转换成gene name;
  
}
}