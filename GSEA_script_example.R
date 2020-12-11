#################################################
#                                               #
#         Gene Set Enrichment Analysis          #
#                                               #
#################################################

########### FILE DESCRIPTION

# NAME: GSEA_script.R
# DATE: 5 december 2020
# AUTHOR: Pernille
#
# DESCRIPTION: Perform GSEA for Manuscript Aregs


library(clusterProfiler)
library(pathview)
library(org.Mm.eg.db)

##---------------------------------------------##
##-----------------Functions-------------------##
##---------------------------------------------##

#' f_gseGO 
#' @description function to run gseGO from clusterprofiler. Perform GSEA on GO terms
#' @param DE_data a data.frame of the differential expression analysis results
#' @param gene_id_col name of the column containing the gene names, if == to "rownames" will use the rownames()
#' @param col_name_stat name of the column containing the values to order the gene by. Recomended to be statistics (t.value or )
#' @param kT keyType argument to give to gseGO, indicates the type of gene names used. Use keytypes(org.Mm.eg.db) to know the accecpted kT
f_gseGO <- function(DE_data, gene_id_col = "rownames", col_name_stat = "stat",  kT = "ENSEMBL",  ontology ="BP"){
  LG <- DE_data[, col_name_stat]
  if(gene_id_col == "rownames"){
    names(LG) <- rownames(DE_data)
  }else{
    names(LG) <- as.character(DE_data[,gene_id_col])
  }
  LG <- sort(LG, decreasing = T)
  #print(head(LG))
  setseed(05122020)
  GSEA_res <- gseGO(geneList = LG, ont = ontology, minGSSize=10, OrgDb = org.Mm.eg.db, keyType = kT,  seed = T)
  
  return(GSEA_res)
}

#' f_gseKEGG 
#' @description function to run gseKEGG from clusterprofiler. Perform GSEA on KEGG pathways
#' @param DE_data a data.frame of the differential expression analysis results
#' @param gene_id_col name of the column containing the gene names, if == to "rownames" will use the rownames()
#' @param col_name_stat name of the column containing the values to order the gene by. Recomended to be statistics (t.value or )
#' @param need_to_convert gseKEGG only accept ENTREZID or UNIPROT (as well as "kegg" and "ncib-proteinid"). 
#' If need_to_convert == T will convert the gene name from genetype to kT
#' @param genetype indicates the type of gene names used in our data that will be converted to kT if need_to_convert == T
#' @param kT can be "ENTREZID", "UNIPROT" (if need_to_convert == F can also be "kegg" or "ncib-proteinid")
f_gseKEGG <- function(DE_data,  gene_id_col = "rownames", col_name_stat = "stat", need_to_convert = T, genetype = "ENSEMBL", kT = "ENTREZID"){
  if(gene_id_col == "rownames"){ g <- rownames(DE_data)}else{g <- as.character(DE_data[,gene_id_col])}

  if(need_to_convert){
    LG <- data.frame(stats = DE_data[, col_name_stat],
                     gene = g,
                     stringsAsFactors = F)
    gene_convert<-bitr(g,
                       fromType = genetype, 
                       toType = kT,
                       OrgDb = "org.Mm.eg.db")
    colnames(gene_convert) <- c("gene", kT)
    LG <- merge(LG, gene_convert)
    LG_kegg <- LG$stats; names(LG_kegg) <- LG[, kT]
  }else{
    LG_kegg <- DE_data[, col_name_stat]; names(LG_kegg) <- g
  }
  LG_kegg <- sort(LG_kegg, decreasing = T)
  if(kT == "ENTREZID"){kT = "ncbi-geneid"}
  if(kT == "UNIPROT"){kT = "uniprot" }
  setseed(05122020)
  GSEA_res <- gseKEGG(geneList = LG_kegg, minGSSize=10, organism = "mmu", keyType = kT,  seed = T)
  
  return(GSEA_res)
}

#' f_GSEA 
#' @description function to run GSEA on specific pathways of interest
#' @param DE_data a data.frame of the differential expression analysis results
#' @param gene_id_col name of the column containing the gene names, if == to "rownames" will use the rownames()
#' @param col_name_stat name of the column containing the values to order the gene by. Recomended to be statistics (t.value or )
#' @param kT keyType argument to give to gseGO, indicates the type of gene names used. Use keytypes(org.Mm.eg.db) to know the accecpted kT
#' @param Term2Gene Terms to test in the Term2Gene format
f_GSEA <- function(DE_data, gene_id_col = "rownames",col_name_stat = "stat",  kT = "ENSEMBL", Term2Gene){
  LG <- DE_data[, col_name_stat]
  if(gene_id_col == "rownames"){
    names(LG) <- rownames(DE_data)
  }else{
    names(LG) <- as.character(DE_data[,gene_id_col])
  }
  LG <- sort(LG, decreasing = T)
  #print(head(LG))
  setseed(05122020)
  GSEA_res <- GSEA(geneList = LG, minGSSize = 1, TERM2GENE = Term2Gene, pvalueCutoff = 1)
  return(GSEA_res)
}

##---------------------------------------------##
##---------------Main Example------------------##
##---------------------------------------------##

  ## -- RNA-seq
head(myDEresults) #In my case DeSeq result, but basically need to be a data.frame or at least accessible as a data.frame
# log2 fold change (MLE): simpleCategories F3pos vs F3neg 
# Wald test p-value: simpleCategories F3pos vs F3neg 
# DataFrame with 6 rows and 7 columns
#                     baseMean log2FoldChange     lfcSE      stat      pvalue        padj      geneID
#                     <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric> <character>
# ENSMUSG00000050830   7.78713        6.69927  0.840312   7.97236 1.55668e-15 3.60756e-14        Vwc2
# ENSMUSG00000023043  14.98143        6.34412  1.104946   5.74157 9.38018e-09 1.05402e-07       Krt18
# ENSMUSG00000034486   5.71999        6.27795  0.893967   7.02257 2.17828e-12 3.72156e-11        Gbx2
# ENSMUSG00000030732   9.86125        6.22052  0.833568   7.46253 8.48794e-14 1.68039e-12      Chrdl2
# ENSMUSG00000049382   5.51357        6.21556  1.391327   4.46736 7.91908e-06 5.63064e-05        Krt8
# ENSMUSG00000051029  10.47732        6.10635  0.618245   9.87691 5.24250e-23 2.05627e-21   Serpinb1b

  ## --Mass-spec
head(myMSresults)
#   Accession Gene.names                             Protein.Desc                                                Protein.Group.Accessions Sorted_AregvsF3Neg_logFC
# 1    E9Q616      Ahnak         AHNAK nucleoprotein (desmoyokin)                                                       E9Q616;A0A494BBD5             1.1171715466
# 2    Q9QXS1       Plec                                  Plectin Q9QXS1;E9Q3W4;A0A3B2W7J8;A0A0R4J218;A0A0R4J221;A0A0R4J223;E9Q9J6;E9PW24             0.5287499078
# 3    A3KGU5     Sptan1 Spectrin alpha chain, non-erythrocytic 1                                                           A3KGU5;E9Q447             0.3484873183
# 4    P16546     Sptan1 Spectrin alpha chain, non-erythrocytic 1                                                                  P16546             0.6384291779
# 5    Q62261     Sptbn1  Spectrin beta chain, non-erythrocytic 1                                  Q62261;A0A0A0MQG2;P15508;Q3UGX2;E9Q397             0.7061508767
# 6    Q8BTM8       Flna                                Filamin-A                               Q8BTM8;B7FAU9;B7FAV1;J3JS91;F6Z2C0;F7AVL7             0.0009865701

# - Prepare mass spec results
uniprots <- data.table::fread("/Users/pernillerainer/Projects/Manuscript_Areg2/MassSpec/uniprot-filtered-organism__Mus+musculus+(Mouse)+[10090]_+AND+revie--.tab", 
                              header = T, sep = "\t" ) 
f <- function(numbers){
  numbers <- strsplit(numbers, split = ";")[[1]]
  isSwissProt <- which(numbers %in% uniprots$Entry)
  
  if(length(isSwissProt)>=1){
    output <- numbers[isSwissProt[1]]
  }else{
    output <- numbers[1]
    print(output)
  }
  return(output)
}
myMSresults$MyAccession <- unlist(lapply(1:nrow(myMSresults), function(x) f(myMSresults[x, "Protein.Group.Accessions"])))
table(myMSresults$MyAccession %in% uniprots$Entry)
# FALSE  TRUE 
#   284  5259 

#### --------------------------
####         GSEA GO

## -- Run GSEA
  # RNA-seq
GSEA_GO <- f_gseGO(myDEresults)
  # Mass spec
GSEA_GO_MS <- f_gseGO(myMSresults, col_name_stat = "Sorted_AregvsF3Neg_logFC", gene_id_col = "MyAccession", kT = "UNIPROT")

## -- Explore data
GSEA_GO@result[grep("lipid",GSEA_GO@result$Description ), 1:7]


#### --------------------------
####        GSEA KEGG
# /!\ KEGG contains a lot of disease associated terms

## -- Run GSEA 
  # RNA-seq
GSEA_KEGG <- f_gseKEGG(myDEresults)
  # Mass specc
GSEA_KEGG_MS <- f_gseKEGG(myMSresults, col_name_stat = "Sorted_AregvsF3Neg_logFC", gene_id_col = "MyAccession", kT = "UNIPROT", need_to_convert = F)

## -- Explore data
GSEA_KEGG@result[ grep("lipid", GSEA_KEGG@result$Description),]

#### --------------------------
####  GSEA of specific terms

# Need to create Term2Gene
  #1. Can download .gmt file (from msigdb for example)
myT2G <- read.gmt("my_path_to_gmt_filt.gmt")

  #2. Can create data.frame 
myT2G <- rbind( data.frame(term = "White fat cell diff - GO_0050872",
                           gene = GSEA_GO@geneSets$`GO:0050872`), #GOOD TO KNOW: GSEA output contain lists of genes from many gene sets!
                data.frame(term = "My custom gene list",
                           gene = my_gene_list))
myT2G <- myT2G[!is.na(myT2G$gene),]

  #(3. To cconvert to another KeyType)
conv <- bitr(myT2G$gene, fromType = "ENSEMBL", toType = "UNIPROT", OrgDb = org.Mm.eg.db, drop = F)
colnames(conv) <- c("gene", "Uniprot")
myT2Uniprot <- merge(myT2G, conv)[, c("term", "Uniprot")]
colnames(myT2Uniprot) <- c("term", "gene"); myT2Uniprot <- myT2Uniprot[!is.na(myT2Uniprot$gene),]

## -- Run gsea
  # RNA-seq
GSEA_myTermsOfInterest <- f_GSEA(DE_data = myDEresults, Term2Gene = myT2G) 
  #Mass spec
GSEA_myTermsOfInterest_MS <- f_GSEA(DE_data = myMSresults, Term2Gene = myT2Uniprot, col_name_stat = "Sorted_AregvsF3Neg_logFC", gene_id_col = "MyAccession") 


#### --------------------------
####       PLOT RESULTS

#1. GSEA Plot
gseaplot(x = GSEA_GO, geneSetID ="GO:0045444", title = paste("NES:", round(GSEA_GO@result["GO:0045444", "NES"],3),
                                                             "padj:", round(GSEA_GO@result["GO:0045444", "p.adjust"], 5)))

#2. Emapplot 
emapplot(GSEA_GO, showCategory = 20)
  #Pathways of interests
selected_pathways <- GSEA_GO@result$Description[GSEA_GO@result$NES > 0][1:50]
emapplot(GSEA_GO, showCategory = selected_pathways, color = "NES")

#3. Dotplot
dotplot(GSEA_GO, showCategory = 20)
  #Dot plot of selected pathways and colored by NES
y <- as.data.frame(GSEA_GO@result[, c("Description", "NES", "p.adjust","setSize")])
y <- subset(y, NES > 0)
y <- y[order(y$NES, decreasing = F),]; y$Description <- factor(y$Description, levels = y$Description)
p <- ggplot(y, 
            aes(x = NES, y = Description)) + 
  geom_point(aes(size = setSize, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.05), low="red") +
  ylab(NULL) +
  ggtitle("GSEA NES > 0")

#4. Ridge plot
ridgeplot(GSEA_GO)

#5. pathview
  #Visualize kegg pathway
myDEresults$logFC_0<-ifelse(myDEresults$padj>0.05 | is.na(myDEresults$padj), 0, myDEresults$log2FoldChange)
genePW <- myDEresults$logFC_0
names(genePW) <- data.annot[rownames(myDEresults), "gene_short_name"] #convert to symbol

pathview(gene.data = genePW,
         pathway.id = "mmu00830",
         species="mmu",
         gene.idtype = "SYMBOL") #data(gene.idtype.list)

#5. heatmap
  # Find genes and lead genes of pathway
All_genes_RetinolM <- GSEA_GO@geneSets[[GSEA_GO@result$ID[grep("retinol metabolic pathway",GSEA_GO@result$Description)]]]
lead_genes_RetinolM <- strsplit(GSEA_GO@result$core_enrichment[grep("retinol metabolic pathway",GSEA_GO@result$Description)], "/")[[1]]
lead_genes_RetinolM <- data.annot[lead_genes_RetinolM, "gene_short_name"]

norm_data <- readRDS("myNormalizedCPMdata.Rds")
norm_data <- norm_data[rownames(norm_data) %in% All_genes_RetinolM,]
rownames(d) <- data.annot[rownames(d),"gene_short_name"]

my_annot_c <- data.frame(row.names = colnames(d),
                         cellType = sapply(strsplit(colnames(d),"_"), `[`, 1)) 
my_annot_r <- data.frame(row.names = rownames(d),
                         lead = rep("No", nrow(d)))
my_annot_r[lead_genes_RetinolM, "lead"] <- "lead"
my_colour = list(
  lead = c("No" = "gray", "lead" = "firebrick"),
  cellType = c("F3pos" = "#08519C", "F3neg" = "#FAE418")
)
pheatmap::pheatmap(norm_data, scale = "row", 
                   clustering_method = "ward.D2", 
                   annotation_colors = my_colour,
                   #color = colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(100),
                   annotation = my_annot_c,
                   annotation_row = my_annot_r,
                   border_color = "NA",
                   cutree_col = 2, cellwidth = 10)




