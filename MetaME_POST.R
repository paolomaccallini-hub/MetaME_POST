# file name: MetaME_POST
#
#-------------------------------------------------------------------------------
# Packages
#-------------------------------------------------------------------------------
#
library(readxl)  
library(httr) 
library(R.utils) 
library(data.table)
library(TissueEnrich)
library(magick)
library(stringr)
library(org.Hs.eg.db)
library(openxlsx)
library(archive)
#
#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------
#
# This function reads MAGMA GSA results
#
read_magma_gsa<-function(path) {
  df <- fread(path, skip = "VARIABLE")
  df <- df[TYPE == "SET"]
  df[, P := as.numeric(P)]
  df[, BETA := as.numeric(BETA)]
  df[, P_bon := pmin(P * nrow(df), 1)]
  df[, P_bh := p.adjust(P, method = "BH")] # Benjamini Hochberg
  df<-df[order(df$P)] 
  return(df)
}
#
# This function reads MAGMA tissue results
#
read_magma_tissue<-function(path) {
  df <- fread(path, skip = "VARIABLE")
  df[, P := as.numeric(P)]
  df[, BETA := as.numeric(BETA)]
  df[, P_bon := pmin(P * nrow(df), 1)]
  df[, P_bh := p.adjust(P, method = "BH")] # Benjamini Hochberg
  df<-df[order(df$P)]
  return(df)
}
#
# This function reads MAGMA cell results
#
read_magma_cell<-function(path) {
  df <- fread(path)
  df[, P := as.numeric(P)]
  df[, BETA := as.numeric(BETA)]
  df[, P_bon := pmin(P * nrow(df), 1)]
  df[, P_bh := p.adjust(P, method = "BH")] # Benjamini Hochberg
  df<-df[order(df$P)]
  return(df)
}
#
# This function reads and edit "MSigDB_20231Hs_MAGMA.txt"
#
read_gmt_ensembl<-function(path) {
  lines <- readLines(path)
  rbindlist(lapply(lines, function(line) {
    fields <- strsplit(line, "\t")[[1]]
    genes  <- fields[-(1:2)]        # skip name and URL
    genes  <- genes[genes != ""]
    data.table(gs_name=fields[1], ensembl_id=genes)
  }))
}
#
# This function reads and edit "gtex_v8_ts_DEG.txt"
#
read_tissue_ensembl<-function(path) {
  df<-fread(path)
  df<-df[grepl("\\.up",df$V1)]
  colnames(df)<-c("FULL_NAME","NGENES","ensembl_id")
  for (i in 1:nrow(df)) df$FULL_NAME[i]<-gsub("\\.up","",df$FULL_NAME[i])
  results_list<-list()
  for (i in 1:nrow(df)) {
    genes<-stringr::str_split_1(df$ensembl_id[i],":")
    results_list[[i]]<-data.table(gs_name=rep(df$FULL_NAME[i],length(genes)),
                            ensembl_id=genes) 
  }
  #
  # Collapse list into a single data table
  #
  results<-rbindlist(results_list)
  return(results)
}
#
# This function performs ORA on "MSigDB_20231Hs_MAGMA.txt" 
# and on edited "gtex_v8_ts_DEG.txt"
#
ora_fuma<-function(my_genes_ensembl,background,msigdb) {
  #
  # Remove genes not in background
  #
  index<-which(my_genes_ensembl%in%background)
  gene_list<-my_genes_ensembl[index]
  #
  # Remove genes not in background
  #
  index<-which(msigdb$ensembl_id%in%background)
  msigdb<-msigdb[index,]
  #
  # Retrieve number of elements in background and in query list
  #
  N1<-length(gene_list)
  NB<-length(background)
  #
  # Retrieve unique gene sets
  #
  gene_sets<-unique(msigdb$gs_name)
  #
  # Where we store the results?
  #
  results_list<-list()
  #
  # Perform ORA
  #
  for (gs in gene_sets) {
    #
    # Current gene set and its size
    #
    gene_set<-msigdb[gs_name==gs,ensembl_id]
    N2<-length(gene_set)
    #
    # Overlap?
    #
    Overlap<-gene_list[gene_list%in%gene_set]
    NR<-length(Overlap)
    #
    # Calculate overlap pval P(X>NR-1)=P(X>=NR)
    #
    if (NR==0) {
      P_val<-1
    } else {
      #
      # We use the notation you find in R documentation.
      # The urn is the set of background genes. The white balls are the genes of
      # the gene set under scrutiny (gs), the black balls are all the remaining
      # genes of the background, the white balls extracted are the overlapping
      # genes between gene_list and gs. Since phyper computes P(X>q) but we want
      # P(X>=q), therefore we consider q-1 instead of q in the function.
      # All that said we have:
      # 
      q<-NR # number of white balls extracted: overlap between gs and gene_list
      m<-N2 # number of white balls in the urn: genes in gs
      n<-NB-N2 # number of black balls in the urn: genes of background not included in gs
      k<-N1 # number of balls extracted: size of gene_list
      P_val<-as.numeric(phyper(q-1,m,n,k,lower.tail=F)) # P(X>NR-1)=P(X>=NR)
    }
    #
    # Save results
    #
    results_list[[gs]]<-data.table(
      FULL_NAME   = gs,
      n_genes_set = N2,
      n_overlap   = NR,
      P           = P_val,
      genes       = paste(Overlap,collapse=";")
    )
  }
  #
  # Collapse list into a single data table
  #
  results<-rbindlist(results_list)
  #
  results[,P_bh:=p.adjust(P,method="BH")]
  results[,P_bon:=p.adjust(P,method="bonferroni")]
  setorder(results, P)
  return(results)
}
#
# This function reads and edit cell type datasets.
# It derives SEGs for each cell type using top decile expression proprotion
# as described in: https://www.nature.com/articles/s41467-024-55611-1
# in Eq.1 and 2
# 
#
make_celltype_genesets<-function(archive_path,file_name,top_quantile=0.9,log2fc_cutoff=1) {
  #
  # Read dataset
  #
  con<-archive_read(archive_path, file = file_name)
  dt<-fread(text = paste(readLines(con), collapse = "\n"))
  close(con)
  #
  # Write on screen what you are doing 
  #
  print(paste("Working on",file_name))
  #
  # Retrieve cell types
  #
  gene_col<-names(dt)[1]
  cell_types<-setdiff(names(dt),c(gene_col, "Average"))
  #
  # Convert to linear scale from log scale
  # 
  mat_lin<-as.matrix(dt[, ..cell_types])
  mat_lin<-2^mat_lin - 1 
  mat_lin[mat_lin<0]<-0 # Handle potential precision noise
  #
  # Calculate expression proportions
  # 
  row_sums<-rowSums(mat_lin)
  prop_mat<-mat_lin/ifelse(row_sums==0,1,row_sums) # avoid division by zero
  #
  mylist<-list() # results here, one element per cell type
  counter<-1
  for (j in seq_along(cell_types)) {
    #
    scores<-prop_mat[, j]
    threshold<-quantile(scores,top_quantile,na.rm=TRUE)
    lfc_scores<-dt[[cell_types[j]]]-dt$Average
    idx<-which(scores>threshold&lfc_scores>log2fc_cutoff) 
    #
    if (length(idx)>0) {
      mylist[[counter]]<-data.table(
        gs_name = cell_types[j],
        ensembl_id = dt[[gene_col]][idx],
        specificity_score = scores[idx],
        log2_fc = lfc_scores[idx]
      )
      counter<-counter+1
    }
  }
  #
  return(rbindlist(mylist))
}
#
# This function add annotations to DrpViz results
#
Annot_DV<-function(dt) {
  for (i in 1:nrow(dt)) {
    tissue<-gsub("DropViz_","",dt$Dataset[i])
    tissue<-gsub("_level2","",tissue)
    full_name<-gsub("_",".",dt$Cell_type[i])
    full_name<-gsub("\\.[^.]+\\.[^.]+$","",full_name)
    index<-which(DropViz_anno$tissue==tissue&DropViz_anno$full_name==full_name)
    if (length(index)>0) {
      dt$class_marker[i]<-DropViz_anno$class_marker[index]
      dt$type_marker[i]<-DropViz_anno$type_marker[index]
      dt$common_name[i]<-DropViz_anno$common_name[index]
    } else {
      dt$class_marker[i]<-NA
      dt$type_marker[i]<-NA
      dt$common_name[i]<-NA
    }
  }
  return(dt)
}
#
#-------------------------------------------------------------------------------
# Create folder if absent
#-------------------------------------------------------------------------------
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Output")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
folder_path<-file.path(current_dir,"Data")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
folder_path<-file.path(current_dir,"Replication")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
folder_path<-file.path(current_dir,"Data/Cell_Type")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
#-----------------------------------------------------------------------------
# Download cell annotations from DropViz
#-----------------------------------------------------------------------------
#
current_dir<-getwd()
url<-paste0("https://storage.googleapis.com/dropviz-downloads/static/annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS")
destfile<-"DropViz_anno.RDS"
file_path<-file.path(current_dir,"Data",destfile)
if(!file.exists(file_path)) {
  print("Downloading data from Drop Viz Atlas")
  RETRY(
    verb = "GET",
    url = url,
    write_disk(file_path, overwrite = TRUE),
    times = 5,           # up to 5 attempts
    pause_min = 5,       # wait 5s between attempts
    terminate_on = c(404, 403) # don't retry on these errors
  )
}  
DropViz_anno<-readRDS("Data/DropViz_anno.RDS")
#
#-----------------------------------------------------------------------------
# Edit scRNA-seq archive (download from FUMA webstite, place it in Data)
#-----------------------------------------------------------------------------
#
# Read files in archive
#
arc<-archive("Data/preprocessed_scrnaseq.tar.gz")
file_list<-arc$path
#
# Filter to L2 resolution, Siletti + Seeker + DropViz only
#
file_list_filtered<-file_list[
  grepl("level2",file_list,ignore.case=T) &
    grepl("Siletti|Seeker|DropViz|Saunders",file_list,ignore.case=T)
]
#
# For each data set build a file to use for ORA
#
for (i in 1:length(file_list_filtered)){
  archive_path<-"Data/preprocessed_scrnaseq.tar.gz"
  file_name<-file_list_filtered[i]
  top_quantile<-0.90 
  log2fc_cutoff<-1
  results<-make_celltype_genesets(archive_path,file_name,top_quantile,log2fc_cutoff)
  file_name<-gsub("celltype/","",file_name)
  write.table(results,paste0("Data/Cell_type/",file_name),sep=",",row.names=F)
}
#
#-----------------------------------------------------------------------------
# Download ME/CFS module form Zhang S. 2025 
# We use supplementary table 2.
#-----------------------------------------------------------------------------
#
current_dir<-getwd()
url<-paste0("https://www.medrxiv.org/content/medrxiv/early/2025/05/11/2025.04.15.25325899/DC9/embed/media-9.xlsx")
destfile<-"media-9.xlsx"
file_path<-file.path(current_dir,"Data",destfile)
if(!file.exists(file_path)) {
  print("Downloading data from Zhang S. et al. 2025 (https://pmc.ncbi.nlm.nih.gov/articles/PMC12047926/)...")
  RETRY(
    verb = "GET",
    url = url,
    write_disk(file_path, overwrite = TRUE),
    times = 5,           # up to 5 attempts
    pause_min = 5,       # wait 5s between attempts
    terminate_on = c(404, 403) # don't retry on these errors
  )
}   
#
#-----------------------------------------------------------------------------
# Read ME/CFS module form Zhang S. 2025 (https://pmc.ncbi.nlm.nih.gov/articles/PMC12047926/)
# We use supplementary table 2.
#-----------------------------------------------------------------------------
#
Module<-read_xlsx("Data/media-9.xlsx") # all 17759 genes of Zhang's study
Zhang_module<-subset.data.frame(Module,q_value<0.02) # 115 genes associated with ME/CFS  
#
#-----------------------------------------------------------------------------
# Retrieve ensembl IDs for Zhang_Module and for Module
#-----------------------------------------------------------------------------
#
mapped<-AnnotationDbi::select(
  org.Hs.eg.db,
  keys=Zhang_module$Gene,
  columns=c("SYMBOL","ENSEMBL"),
  keytype="SYMBOL"
)
my_genes_ensembl<-unique(mapped$ENSEMBL[!is.na(mapped$ENSEMBL)])
#
mapped<-AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = Module$Gene,
  columns = c("SYMBOL", "ENSEMBL"),
  keytype = "SYMBOL"
)
my_universe_ensembl<-unique(mapped$ENSEMBL[!is.na(mapped$ENSEMBL)])
#
#-------------------------------------------------------------------------------
# Discovery: DME_1_MVP
# Replication: Zhang
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
# Gene-set enrichment analysis
#-------------------------------------------------------------------------------
#
replication_path<-"Zhang"
#
# Perform ORA on Zhang module
#
msigdb<-read_gmt_ensembl("MSigDB_20231Hs_MAGMA.txt")
background<-my_universe_ensembl
Replication<-ora_fuma(my_genes_ensembl,background,msigdb)
write.table(Replication,paste0("Replication/",replication_path,"_to_Replication_GSEA.csv"),
            sep=",",row.names=F)
#
#-------------------------------------------------------------------------------
# Tissue enrichment analysis
#-------------------------------------------------------------------------------
#
replication_path<-"Zhang"
#
# Perform ORA on Zhang module
#
GTExv8<-read_tissue_ensembl("gtex_v8_ts_DEG.txt")
background<-my_universe_ensembl
Replication<-ora_fuma(my_genes_ensembl,background,GTExv8)
write.table(Replication,paste0("Replication/",replication_path,"_to_Replication_Tissue.csv"),
            sep=",",row.names=F)
#
#-------------------------------------------------------------------------------
# DropViz cell type analysis
#-------------------------------------------------------------------------------
#
discovery_path<-"DME_1_MVP"
replication_path<-"Zhang"
#
# Read Discovery 
#
Discovery<-read_magma_cell("FUMA/DME_1_MVP/DropViz_L2/magma_celltype_step1.txt")
Discovery<-unique(Discovery)
#
# Perform ORA on Zhang module
#
Replication_list<-list()
background<-my_universe_ensembl
unique_ds<-unique(Discovery$Dataset)
for (i in 1:length(unique_ds)) {
  file_name<-unique_ds[i]
  dataset<-fread(paste0("Data/Cell_type/",file_name,".txt"))
  Replication_list[[i]]<-ora_fuma(my_genes_ensembl,background,dataset)
  Replication_list[[i]]$Dataset<-rep(file_name,nrow(Replication_list[[i]]))
}
Replication<-rbindlist(Replication_list)
for (j in 1:ncol(Replication)) {
  if (colnames(Replication)[j]=="FULL_NAME") {
    colnames(Replication)[j]<-"Cell_type"
  }
}
Replication<-unique(Replication)
#
# Calculate P_bon and P_bh
#
Replication[, P_bon := pmin(P * nrow(Replication), 1)]
Replication[, P_bh := p.adjust(P, method = "BH")] 
write.table(Replication,paste0("Replication/",replication_path,"_to_Replication_DropViz_L2.csv"),
            sep=",",row.names=F)
#
#-------------------------------------------------------------------------------
# Siletti-Seeker cell type analysis with Bonferroni correction
# Note: FUMA v1.8.3 produces a duplicate entry for
# 16_Siletti_CerebralCortex.IFG.A44-A45_Human_2022_level2
# in the dataset list. The duplicate is present in the output
# but Bonferroni correction is applied correctly by FUMA.
# The duplicate has been removed from the dataset list in this
# script to avoid confusion in downstream processing.
#-------------------------------------------------------------------------------
#
discovery_path<-"DME_1_MVP"
replication_path<-"Zhang"
#
# Read Discovery
#
Discovery<-read_magma_cell("FUMA/DME_1_MVP/Siletti_Seeker_L2/magma_celltype_step1.txt")
Discovery<-unique(Discovery) # I found duplicates
#
# Perform ORA on Zhang module
#
Replication_list<-list()
background<-my_universe_ensembl
unique_ds<-unique(Discovery$Dataset)
for (i in 1:length(unique_ds)) {
  file_name<-unique_ds[i]
  dataset<-fread(paste0("Data/Cell_type/",file_name,".txt"))
  Replication_list[[i]]<-ora_fuma(my_genes_ensembl,background,dataset)
  Replication_list[[i]]$Dataset<-rep(file_name,nrow(Replication_list[[i]]))
}
Replication<-rbindlist(Replication_list)
for (j in 1:ncol(Replication)) {
  if (colnames(Replication)[j]=="FULL_NAME") {
    colnames(Replication)[j]<-"Cell_type"
  }
}
Replication<-unique(Replication)
#
# Calculate P_bon and P_bh
#
Replication[, P_bon := pmin(P * nrow(Replication), 1)]
Replication[, P_bh := p.adjust(P, method = "BH")] 
write.table(Replication,paste0("Replication/",replication_path,"_to_Replication_Siletti_Seeker_L2.csv"),
            sep=",",row.names=F)
#
#-------------------------------------------------------------------------------
# Output and Supplementary material
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
# README
#-------------------------------------------------------------------------------
#
readme <- data.frame(
  Field = c(
    #---------------------------------------------------------------------------
    "Table S1",
    #---------------------------------------------------------------------------
    "topLeadSNP", "GRCh38", "GRCh37", "start", "end",
    "META_Z", "META_MAF", "META_Log10P",
    "DME_1_Z", "DME_1_MAF", "DME_1_Log10P",
    "MVP_Z", "MVP_MAF", "MVP_Log10P",
    "IndSigSNPs",
    "",
    #---------------------------------------------------------------------------
    "Table S2",
    #---------------------------------------------------------------------------
    "topLeadSNP", "GRCh38", "GRCh37", "start", "end",
    "DME_1_Z", "DME_1_MAF", "DME_1_Log10P",
    "META_Z", "META_MAF", "META_Log10P",
    "MVP_Z", "MVP_MAF", "MVP_Log10P",
    "",
    #---------------------------------------------------------------------------
    "Table S3",
    #---------------------------------------------------------------------------
    "FULL_NAME", "NGENES", "BETA", "SE", "P", "P_bon", "P_bh",
    "BETA_DME_1", "SE_DME_1", "P_DME_1", "P_bon_DME_1", "P_bh_DME_1",
    "BETA_MVP", "SE_MVP", "P_MVP", "P_bon_MVP", "P_bh_MVP",
    "P_Zhang", "P_bon_Zhang", "P_bh_Zhang",
    "n_gene_set_Zhang","n_overlap_Zhang","genes_Zhang",
    "",
    #---------------------------------------------------------------------------
    "Table S4",
    #---------------------------------------------------------------------------
    "FULL_NAME", "NGENES", "BETA", "SE", "P", "P_bon", "P_bh",
    "BETA_DME_1", "SE_DME_1", "P_DME_1", "P_bon_DME_1", "P_bh_DME_1",
    "BETA_MVP", "SE_MVP", "P_MVP", "P_bon_MVP", "P_bh_MVP",
    "P_Zhang", "P_bon_Zhang", "P_bh_Zhang",
    "n_gene_set_Zhang","n_overlap_Zhang","genes_Zhang",
    "",
    #---------------------------------------------------------------------------
    "Table S5",
    #---------------------------------------------------------------------------
    "Dataset", "Cell_type", "BETA", "SE", "P", "P_bon", "P_bh",
    "BETA_DME_1", "SE_DME_1", "P_DME_1", "P_bon_DME_1", "P_bh_DME_1",
    "BETA_MVP", "SE_MVP", "P_MVP", "P_bon_MVP", "P_bh_MVP",
    "P_Zhang", "P_bon_Zhang", "P_bh_Zhang",
    "n_gene_set_Zhang","n_overlap_Zhang","genes_Zhang",
    "",
    #---------------------------------------------------------------------------
    "Table S6",
    #---------------------------------------------------------------------------
    "Dataset", "Cell_type", "BETA", "SE", "P", "P_bon", "P_bh",
    "BETA_DME_1", "SE_DME_1", "P_DME_1", "P_bon_DME_1", "P_bh_DME_1",
    "BETA_MVP", "SE_MVP", "P_MVP", "P_bon_MVP", "P_bh_MVP",
    "P_Zhang", "P_bon_Zhang", "P_bh_Zhang",
    "n_gene_set_Zhang","n_overlap_Zhang","genes_Zhang",
    ""
  ),
  Description = c(
    #---------------------------------------------------------------------------
    "Genomic risk loci for the DME-1+MVP meta-analysis. Lead SNPs and locus boundaries.",
    #---------------------------------------------------------------------------
    "rsID of the top lead SNP in the genomic risk locus",
    "Genomic coordinates of the lead SNP (GRCh38): CHR:BP:A1:A2",
    "Genomic coordinates of the lead SNP (GRCh37): CHR:BP:A1:A2",
    "Start position of the genomic risk locus (GRCh37)",
    "End position of the genomic risk locus (GRCh37)",
    "Z-score of the lead SNP in the DME-1+MVP meta-analysis",
    "Effect allele frequency in the DME-1+MVP meta-analysis (weighted average)",
    "-log10(P) of the lead SNP in the DME-1+MVP meta-analysis",
    "Z-score of the lead SNP in the DecodeME (DME-1) GWAS",
    "Effect allele frequency in DME-1",
    "-log10(P) of the lead SNP in DME-1",
    "Z-score of the lead SNP in the Million Veteran Program (MVP) GWAS",
    "Effect allele frequency in MVP",
    "-log10(P) of the lead SNP in MVP",
    "All independent significant SNPs in the locus (r2 < 0.6)",
    "",
    #---------------------------------------------------------------------------
    "Genomic risk loci for DecodeME (DME-1) alone. Lead SNPs and locus boundaries.",
    #---------------------------------------------------------------------------
    "rsID of the top lead SNP in the genomic risk locus",
    "Genomic coordinates of the lead SNP (GRCh38): CHR:BP:A1:A2",
    "Genomic coordinates of the lead SNP (GRCh37): CHR:BP:A1:A2",
    "Start position of the genomic risk locus (GRCh37)",
    "End position of the genomic risk locus (GRCh37)",
    "Z-score of the lead SNP in DME-1",
    "Effect allele frequency in DME-1",
    "-log10(P) of the lead SNP in DME-1",
    "Z-score of the lead SNP in the DME-1+MVP meta-analysis",
    "Effect allele frequency in the DME-1+MVP meta-analysis (weighted average)",
    "-log10(P) of the lead SNP in the DME-1+MVP meta-analysis",
    "Z-score of the lead SNP in MVP",
    "Effect allele frequency in MVP",
    "-log10(P) of the lead SNP in MVP",
    "",
    #---------------------------------------------------------------------------
    "MAGMA competitive gene-set enrichment analysis (MSigDB 2023.1.Hs).",
    #---------------------------------------------------------------------------
    "Full gene set name from MSigDB 2023.1.Hs",
    "Number of genes in the gene set",
    "MAGMA regression coefficient (DME-1+MVP meta-analysis)",
    "Standard error of the regression coefficient (DME-1+MVP meta-analysis)",
    "Raw P-value (DME-1+MVP meta-analysis)",
    "Bonferroni-corrected P-value (DME-1+MVP meta-analysis)",
    "Benjamini-Hochberg-adjusted P-value (DME-1+MVP meta-analysis)",
    "MAGMA regression coefficient (DME-1 alone)",
    "Standard error (DME-1 alone)",
    "Raw P-value (DME-1 alone)",
    "Bonferroni-corrected P-value (DME-1 alone)",
    "Benjamini-Hochberg-adjusted P-value (DME-1 alone)",
    "MAGMA regression coefficient (MVP alone)",
    "Standard error (MVP alone)",
    "Raw P-value (MVP alone)",
    "Bonferroni-corrected P-value (MVP alone)",
    "Benjamini-Hochberg-adjusted P-value (MVP alone)",
    "Raw P-value from ORA of 115 Zhang et al. (2025) ME/CFS candidate genes against the gene set; background: 17,759 STRING network genes",
    "Bonferroni-corrected P-value for Zhang ORA",
    "Benjamini-Hochberg-adjusted P-value for Zhang ORA",
    "size of the gene set in Zhang ORA",
    "size of the overlap in Zhang ORA",
    "IDs of overlapping genes",
    "",
    #---------------------------------------------------------------------------
    "MAGMA gene-tissue expression analysis (GTEx v8; 54 tissues; 17,280 genes). All 54 tissues shown.",
    #---------------------------------------------------------------------------
    "Full tissue name from GTEx v8",
    "Number of genes tested (17,280 for all tissues)",
    "MAGMA regression coefficient (DME-1+MVP meta-analysis)",
    "Standard error of the regression coefficient (DME-1+MVP meta-analysis)",
    "Raw P-value (DME-1+MVP meta-analysis)",
    "Bonferroni-corrected P-value (DME-1+MVP meta-analysis; k=54)",
    "Benjamini-Hochberg-adjusted P-value (DME-1+MVP meta-analysis)",
    "MAGMA regression coefficient (DME-1 alone)",
    "Standard error (DME-1 alone)",
    "Raw P-value (DME-1 alone)",
    "Bonferroni-corrected P-value (DME-1 alone; k=54)",
    "Benjamini-Hochberg-adjusted P-value (DME-1 alone)",
    "MAGMA regression coefficient (MVP alone)",
    "Standard error (MVP alone)",
    "Raw P-value (MVP alone)",
    "Bonferroni-corrected P-value (MVP alone; k=54)",
    "Benjamini-Hochberg-adjusted P-value (MVP alone)",
    "Raw P-value from ORA of 115 Zhang et al. (2025) ME/CFS candidate genes against upregulated genes per tissue (gtex_v8_ts_DEG.txt); background: 17,759 STRING network genes",
    "Bonferroni-corrected P-value for Zhang ORA (k=54)",
    "Benjamini-Hochberg-adjusted P-value for Zhang ORA",
    "size of the gene set in Zhang ORA",
    "size of the overlap in Zhang ORA",
    "IDs of overlapping genes",
    "",
    #---------------------------------------------------------------------------
    "FUMA cell-type analysis, DropViz mouse brain atlas (9 regions; level-2 resolution; 565 cell types total). All tested cell types shown.",
    #---------------------------------------------------------------------------
    "DropViz dataset identifier (region)",
    "Cell type label at level-2 resolution (DropViz nomenclature: CellClass.ClassMarker.TypeMarker)",
    "MAGMA regression coefficient (DME-1+MVP meta-analysis)",
    "Standard error of the regression coefficient (DME-1+MVP meta-analysis)",
    "Raw P-value from step 1 of the FUMA cell-type pipeline (DME-1+MVP meta-analysis)",
    "Bonferroni-corrected P-value (DME-1+MVP meta-analysis, all cell types across all DropViz datasets)",
    "Benjamini-Hochberg-adjusted P-value (DME-1+MVP meta-analysis)",
    "MAGMA regression coefficient (DME-1 alone)",
    "Standard error (DME-1 alone)",
    "Raw P-value (DME-1 alone)",
    "Bonferroni-corrected P-value (DME-1 alone)",
    "Benjamini-Hochberg-adjusted P-value (DME-1 alone)",
    "MAGMA regression coefficient (MVP alone)",
    "Standard error (MVP alone)",
    "Raw P-value (MVP alone)",
    "Bonferroni-corrected P-value (MVP alone)",
    "Benjamini-Hochberg-adjusted P-value (MVP alone)",
    "Raw P-value from ORA of 115 Zhang et al. (2025) ME/CFS candidate genes against cell-type foreground gene set; background: 17,759 STRING network genes",
    "Bonferroni-corrected P-value for Zhang ORA (k=565)",
    "Benjamini-Hochberg-adjusted P-value for Zhang ORA",
    "size of the gene set in Zhang ORA",
    "size of the overlap in Zhang ORA",
    "IDs of overlapping genes",
    "",
    #---------------------------------------------------------------------------
    "FUMA cell-type analysis, Siletti et al. (2023) + Seeker et al. (2023) human brain atlases (107 datasets; level-2 resolution; 2,099 cell types total). All tested cell types shown.",
    #---------------------------------------------------------------------------
    "Dataset identifier (Siletti or Seeker dissection label)",
    "Cell type label at level-2 resolution",
    "MAGMA regression coefficient (DME-1+MVP meta-analysis)",
    "Standard error of the regression coefficient (DME-1+MVP meta-analysis)",
    "Raw P-value from step 1 of the FUMA cell-type pipeline (DME-1+MVP meta-analysis)",
    "Bonferroni-corrected P-value (DME-1+MVP meta-analysis, all cell types across all Siletti/Seeker datasets)",
    "Benjamini-Hochberg-adjusted P-value (DME-1+MVP meta-analysis)",
    "MAGMA regression coefficient (DME-1 alone)",
    "Standard error (DME-1 alone)",
    "Raw P-value (DME-1 alone)",
    "Bonferroni-corrected P-value (DME-1 alone)",
    "Benjamini-Hochberg-adjusted P-value (DME-1 alone)",
    "MAGMA regression coefficient (MVP alone)",
    "Standard error (MVP alone)",
    "Raw P-value (MVP alone)",
    "Bonferroni-corrected P-value (MVP alone)",
    "Benjamini-Hochberg-adjusted P-value (MVP alone)",
    "Raw P-value from ORA of 115 Zhang et al. (2025) ME/CFS candidate genes against cell-type foreground gene set; background: 17,759 STRING network genes",
    "Bonferroni-corrected P-value for Zhang ORA",
    "Benjamini-Hochberg-adjusted P-value for Zhang ORA",
    "size of the gene set in Zhang ORA",
    "size of the overlap in Zhang ORA",
    "IDs of overlapping genes",
    ""
  ),
  stringsAsFactors = FALSE
)
#
# Write to Excel
#
wb_name<-"S1_File"
wb<-createWorkbook()
sheet_name<-"README" # name the sheet
#
addWorksheet(wb,sheet_name)
writeData(wb,sheet_name,readme, colNames = FALSE)
table_rows <- which(grepl("^Table S", readme$Field))
for (r in table_rows) {
  addStyle(wb,sheet_name,
           createStyle(textDecoration = "bold"),
           rows = r, cols = 1:2, stack = TRUE)
}
setColWidths(wb,sheet_name,cols=1,widths=35)
setColWidths(wb,sheet_name,cols=2,widths=80)
saveWorkbook(wb,paste0(wb_name,".xlsx"),overwrite=T)
#
#-------------------------------------------------------------------------------
# Supplementary material: Risk loci for DME_1_MVP
#-------------------------------------------------------------------------------
#
populations<-c("DME_1","MVP")
myriskloci<-fread("FUMA/DME_1_MVP/MAGMA/GenomicRiskLoci.txt")
head(myriskloci)
#
# Read meta-analysis
#
file_name<-paste0(current_dir,"/Output/GWAS_METAL_DME_1_MVP_GRCh38.tsv.gz")
mydata<-fread(file_name,sep="\t")
#
# Find lead SNPs in meta analysis
#
mybest<-list()
index<-which(mydata$SNP%in%myriskloci$rsID)
mybest[[1]]<-mydata[index,]
names(mybest)[1]<-"METAL"
#
# Recover the best SNPs from the input GWAS
#
for (i in 1:length(populations)) {
  file_name_gz<-paste0(current_dir,"/Munged/",populations[i],"_GRCh38.tsv.gz")
  mydata<-fread(file_name_gz)
  df<-mybest[[1]][,c("SNP","A1","A2")]
  mybest[[1+i]]<-merge(mydata,df,by=c("SNP","A1","A2"),all.x=F)
  mybest[[1+i]]<-mybest[[1+i]][mybest[[1]]$SNP,on="SNP"] # ask for the same order in SNPs
  names(mybest)[1+i]<-populations[i]
  remove(mydata)
}
#
# Build a custom table for risk loci
#
myRL<-data.frame(topLeadSNP=myriskloci$rsID)
for (i in 1:nrow(myRL)) {
  myRL$GRCh38[i]<-paste0(c(mybest[[1]]$CHR[i],mybest[[1]]$BP[i],mybest[[1]]$A1[i],mybest[[1]]$A2[i]),collapse=":")
}
myRL$GRCh37<-myriskloci$uniqID
myRL$start<-myriskloci$start
myRL$end<-myriskloci$end
myRL$META_Z<-mybest[[1]]$Z
myRL$META_MAF<-mybest[[1]]$FRQ
myRL$META_Log10P<-round(-log10(mybest[[1]]$P),2)
myRL$DME_1_Z<-mybest[[2]]$Z
myRL$DME_1_MAF<-mybest[[2]]$FRQ
myRL$DME_1_Log10P<-round(-log10(mybest[[2]]$P),2)
myRL$MVP_Z<-mybest[[3]]$Z
myRL$MVP_MAF<-mybest[[3]]$FRQ
myRL$MVP_Log10P<-round(-log10(mybest[[3]]$P),2)
myRL$IndSigSNPs<-myriskloci$IndSigSNPs
#
# Write to Excel
#
wb_name<-"S1_File"
wb<-loadWorkbook(paste0(wb_name,".xlsx"))
myresult<-myRL
sheet_name<-"S1 Risk Loci Meta" # name the sheet
#
addWorksheet(wb,sheet_name)
writeData(wb,sheet_name,myresult)
addStyle(wb,sheet_name,createStyle(textDecoration ="bold",border="Bottom"),
         rows=1,cols=1:ncol(myresult))
setColWidths(wb,sheet_name,cols=1:ncol(myresult),widths="auto")
freezePane(wb,sheet_name,firstRow=T)
saveWorkbook(wb,paste0(wb_name,".xlsx"),overwrite=T)
#
#-------------------------------------------------------------------------------
# Supplementary table: Risk loci form DME_1
#-------------------------------------------------------------------------------
#
myriskloci<-fread("FUMA/DME_1/MAGMA/GenomicRiskLoci.txt")
head(myriskloci)
#
# Read DME_1
#
file_name<-paste0(current_dir,"/Munged/DME_1_GRCh38.tsv.gz")
mydata<-fread(file_name,sep="\t")
#
# Find lead SNPs in DME_1
#
mybest<-list()
index<-which(mydata$SNP%in%myriskloci$rsID)
mybest[[1]]<-mydata[index,]
names(mybest)[1]<-"METAL"
#
# Recover the best SNPs from the meta analysis and MVP
#
i<-1
file_name<-paste0(current_dir,"/Output/GWAS_METAL_DME_1_MVP_GRCh38.tsv.gz")
mydata<-fread(file_name)
df<-mybest[[1]][,c("SNP","A1","A2")]
mybest[[1+i]]<-merge(mydata,df,by=c("SNP","A1","A2"),all.x=F)
mybest[[1+i]]<-mybest[[1+i]][mybest[[1]]$SNP,on="SNP"] # ask for the same order in SNPs
names(mybest)[1+i]<-"DME_1_MVP"
remove(mydata)
i<-2
file_name<-paste0(current_dir,"/Munged/MVP_GRCh38.tsv.gz")
mydata<-fread(file_name)
df<-mybest[[1]][,c("SNP","A1","A2")]
mybest[[1+i]]<-merge(mydata,df,by=c("SNP","A1","A2"),all.x=F)
mybest[[1+i]]<-mybest[[1+i]][mybest[[1]]$SNP,on="SNP"] # ask for the same order in SNPs
names(mybest)[1+i]<-"MVP"
remove(mydata)
#
# Build a custom table for risk loci
#
myRL<-data.frame(topLeadSNP=myriskloci$rsID)
for (i in 1:nrow(myRL)) {
  myRL$GRCh38[i]<-paste0(c(mybest[[1]]$CHR[i],mybest[[1]]$BP[i],mybest[[1]]$A1[i],mybest[[1]]$A2[i]),collapse=":")
}
myRL$GRCh37<-myriskloci$uniqID
myRL$start<-myriskloci$start
myRL$end<-myriskloci$end
myRL$DME_1_Z<-mybest[[1]]$Z
myRL$DME_1_MAF<-mybest[[1]]$FRQ
myRL$DME_1_Log10P<-round(-log10(mybest[[1]]$P),2)
myRL$META_Z<-mybest[[2]]$Z
myRL$META_MAF<-mybest[[2]]$FRQ
myRL$META_Log10P<-round(-log10(mybest[[2]]$P),2)
myRL$MVP_Z<-mybest[[3]]$Z
myRL$MVP_MAF<-mybest[[3]]$FRQ
myRL$MVP_Log10P<-round(-log10(mybest[[3]]$P),2)
#
# Write to Excel
#
wb<-loadWorkbook(paste0(wb_name,".xlsx"))
sheet_name<-"S2 Risk Loci DME_1" # name the sheet
myresult<-myRL
#
addWorksheet(wb,sheet_name)
writeData(wb,sheet_name,myresult)
addStyle(wb,sheet_name,
         createStyle(textDecoration ="bold",border="Bottom"),
         rows=1,cols=1:ncol(myresult))
setColWidths(wb,sheet_name,cols=1:ncol(myresult),widths="auto")
freezePane(wb,sheet_name,firstRow=T)
saveWorkbook(wb,paste0(wb_name,".xlsx"),overwrite=T)
#
#-------------------------------------------------------------------------------
# Supplementary table: Tissue analysis
#-------------------------------------------------------------------------------
#
DME_1_MVP<-read_magma_tissue("FUMA/DME_1_MVP/MAGMA/magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out")
DME_1<-read_magma_tissue("FUMA/DME_1/MAGMA/magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out")
MVP<-read_magma_tissue("FUMA/MVP/MAGMA/magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out")
Zhang<-fread("Replication/Zhang_to_Replication_Tissue.csv")
#
# Are there duplicated lines? All data tables?
#
setDT(DME_1_MVP); DME_1_MVP <- unique(DME_1_MVP)
setDT(DME_1);     DME_1     <- unique(DME_1)
setDT(MVP);       MVP       <- unique(MVP)
setDT(Zhang);     Zhang     <- unique(Zhang)
#
DME_1_MVP <- merge(DME_1_MVP, 
                   DME_1[, .(FULL_NAME, BETA_DME_1 = BETA, SE_DME_1 = SE, 
                             P_DME_1 = P, P_bon_DME_1 = P_bon, P_bh_DME_1 = P_bh)], 
                   by = "FULL_NAME", all.x = TRUE)

DME_1_MVP <- merge(DME_1_MVP, 
                   MVP[, .(FULL_NAME, BETA_MVP = BETA, SE_MVP = SE, 
                           P_MVP = P, P_bon_MVP = P_bon, P_bh_MVP = P_bh)], 
                   by = "FULL_NAME", all.x = TRUE)

DME_1_MVP <- merge(DME_1_MVP, 
                   Zhang[, .(FULL_NAME, P_Zhang = P, P_bon_Zhang = P_bon, P_bh_Zhang = P_bh, 
                             n_gene_set_Zhang = n_genes_set, n_overlap_Zhang = n_overlap, 
                             genes_Zhang = genes)], 
                   by = "FULL_NAME", all.x = TRUE)
#
DME_1_MVP<-subset.data.frame(DME_1_MVP,select=c("FULL_NAME","NGENES","BETA","SE",
                                                "P","P_bon","P_bh",
                                                "BETA_DME_1","SE_DME_1",
                                                "P_DME_1","P_bon_DME_1","P_bh_DME_1",
                                                "BETA_MVP","SE_MVP",
                                                "P_MVP","P_bon_MVP","P_bh_MVP",
                                                "P_Zhang","P_bon_Zhang","P_bh_Zhang",
                                                "n_gene_set_Zhang","n_overlap_Zhang","genes_Zhang"))
DME_1_MVP<-DME_1_MVP[order(DME_1_MVP$P),]
#
# Write a table with only significant results (BH)
#
dt<-subset.data.frame(DME_1_MVP,P_bh<=0.05)
dt[, P_bon_rep_Zhang := pmin(P_Zhang * nrow(dt), 1)] 
write.table(dt,"Table_Tissue_BH.csv",sep=",",row.names=F,quote=F)
#
# Write a table with only significant results (Bon)
#
dt<-subset.data.frame(DME_1_MVP,P_bon<=0.05)
dt[, P_bon_rep_Zhang := pmin(P_Zhang * nrow(dt), 1)] 
write.table(dt,"Table_Tissue_Bon.csv",sep=",",row.names=F,quote=F)
#
# Load existing workbook
#
wb<-loadWorkbook(paste0(wb_name,".xlsx"))
sheet_name<-"S4 Tissue Analysis" # name the sheet
myresult<-DME_1_MVP
#
addWorksheet(wb,sheet_name)
writeData(wb,sheet_name,myresult)
addFilter(wb,sheet_name,row=1,cols=1:ncol(myresult))
addStyle(wb,sheet_name,
         createStyle(textDecoration ="bold",border="Bottom"),
         rows=1,cols=1:ncol(myresult))
setColWidths(wb,sheet_name,cols=1:ncol(myresult),widths="auto")
freezePane(wb,sheet_name,firstRow=T)
saveWorkbook(wb,paste0(wb_name,".xlsx"),overwrite=T)
#
#-------------------------------------------------------------------------------
# Supplementary table: Gene Set analyses
#-------------------------------------------------------------------------------
#
DME_1_MVP<-read_magma_gsa("FUMA/DME_1_MVP/MAGMA/magma.gsa.out")
DME_1<-read_magma_gsa("FUMA/DME_1/MAGMA/magma.gsa.out")
MVP<-read_magma_gsa("FUMA/MVP/MAGMA/magma.gsa.out")
Zhang<-fread("Replication/Zhang_to_Replication_GSEA.csv")
#
# Are there duplicated lines? All data tables?
#
setDT(DME_1_MVP); DME_1_MVP <- unique(DME_1_MVP)
setDT(DME_1);     DME_1     <- unique(DME_1)
setDT(MVP);       MVP       <- unique(MVP)
setDT(Zhang);     Zhang     <- unique(Zhang)
#
# Sequential Merges (Left Joins)
# 
res <- merge(DME_1_MVP, 
             DME_1[, .(FULL_NAME, BETA_DME_1 = BETA, SE_DME_1 = SE, P_DME_1 = P, 
                       P_bon_DME_1 = P_bon, P_bh_DME_1 = P_bh)], 
             by = "FULL_NAME", all.x = TRUE)

res <- merge(res, 
             MVP[, .(FULL_NAME, BETA_MVP = BETA, SE_MVP = SE, P_MVP = P, 
                     P_bon_MVP = P_bon, P_bh_MVP = P_bh)], 
             by = "FULL_NAME", all.x = TRUE)

res <- merge(res, 
             Zhang[, .(FULL_NAME, P_Zhang = P, P_bon_Zhang = P_bon, P_bh_Zhang = P_bh, 
                       n_gene_set_Zhang = n_genes_set, n_overlap_Zhang = n_overlap, 
                       genes_Zhang = genes)], 
             by = "FULL_NAME", all.x = TRUE)
#
DME_1_MVP<-res
DME_1_MVP<-subset.data.frame(DME_1_MVP,select=c("FULL_NAME","NGENES","BETA","SE",
                                                "P","P_bon","P_bh",
                                                "BETA_DME_1","SE_DME_1",
                                                "P_DME_1","P_bon_DME_1","P_bh_DME_1",
                                                "BETA_MVP","SE_MVP",
                                                "P_MVP","P_bon_MVP","P_bh_MVP",
                                                "P_Zhang","P_bon_Zhang","P_bh_Zhang",
                                                "n_gene_set_Zhang","n_overlap_Zhang","genes_Zhang"))
#
DME_1_MVP<-DME_1_MVP[order(DME_1_MVP$P),]
#
# Write a table with only significant results (BH)
#
dt<-subset.data.frame(DME_1_MVP,P_bh<=0.05)
dt[, P_bon_rep_Zhang := pmin(P_Zhang * nrow(dt), 1)] 
write.table(dt,"Table_Gene_Set_BH.csv",sep=",",row.names=F,quote=F)
#
# Write a table with only significant results (Bon)
#
dt<-subset.data.frame(DME_1_MVP,P_bon<=0.05)
dt[, P_bon_rep_Zhang := pmin(P_Zhang * nrow(dt), 1)] 
write.table(dt,"Table_Gene_Set_Bon.csv",sep=",",row.names=F,quote=F)
#
# Write to Excel
#
wb<-loadWorkbook(paste0(wb_name,".xlsx"))
sheet_name<-"S3 Gene Set Analysis" # name the sheet
myresult<-DME_1_MVP
#
addWorksheet(wb,sheet_name)
writeData(wb,sheet_name,myresult)
addFilter(wb,sheet_name,row=1,cols=1:ncol(myresult))
addStyle(wb,sheet_name,
         createStyle(textDecoration ="bold",border="Bottom"),
         rows=1,cols=1:ncol(myresult))
setColWidths(wb,sheet_name,cols=1:ncol(myresult),widths="auto")
freezePane(wb,sheet_name,firstRow=T)
saveWorkbook(wb,paste0(wb_name,".xlsx"),overwrite=T)
#
#-------------------------------------------------------------------------------
# Supplementary table: DropViz L2
#-------------------------------------------------------------------------------
#
DME_1_MVP<-fread("FUMA/DME_1_MVP/DropViz_L2/magma_celltype_step1.txt")
DME_1<-fread("FUMA/DME_1/DropViz_L2/magma_celltype_step1.txt")
MVP<-fread("FUMA/MVP/DropViz_L2/magma_celltype_step1.txt")
Zhang<-fread("Replication/Zhang_to_Replication_DropViz_L2.csv")
#
# Are there duplicated lines? All data tables?
#
setDT(DME_1_MVP); DME_1_MVP <- unique(DME_1_MVP)
setDT(DME_1);     DME_1     <- unique(DME_1)
setDT(MVP);       MVP       <- unique(MVP)
setDT(Zhang);     Zhang     <- unique(Zhang)
#
# Calculate P_bh and P_bon 
#
DME_1_MVP[, P_bon := pmin(P * nrow(DME_1_MVP), 1)]
DME_1_MVP[, P_bh := p.adjust(P, method = "BH")]
DME_1[, P_bon := pmin(P * nrow(DME_1), 1)]
DME_1[, P_bh := p.adjust(P, method = "BH")]
MVP[, P_bon := pmin(P * nrow(MVP), 1)]
MVP[, P_bh := p.adjust(P, method = "BH")]
Zhang[, P_bon := pmin(P * nrow(Zhang), 1)]
Zhang[, P_bh := p.adjust(P, method = "BH")]
#
DME_1_MVP <- merge(DME_1_MVP, 
                   DME_1[, .(Cell_type, Dataset, 
                             BETA_DME_1 = BETA, SE_DME_1 = SE, P_DME_1 = P, 
                             P_bon_DME_1 = P_bon, P_bh_DME_1 = P_bh)], 
                   by = c("Cell_type", "Dataset"), all.x = TRUE)

DME_1_MVP <- merge(DME_1_MVP, 
                   MVP[, .(Cell_type, Dataset, 
                           BETA_MVP = BETA, SE_MVP = SE, P_MVP = P, 
                           P_bon_MVP = P_bon, P_bh_MVP = P_bh)], 
                   by = c("Cell_type", "Dataset"), all.x = TRUE)

DME_1_MVP <- merge(DME_1_MVP, 
                   Zhang[, .(Cell_type, Dataset, 
                             P_Zhang = P, P_bon_Zhang = P_bon, P_bh_Zhang = P_bh, 
                             n_gene_set_Zhang = n_genes_set, n_overlap_Zhang = n_overlap, 
                             genes_Zhang = genes)], 
                   by = c("Cell_type", "Dataset"), all.x = TRUE)
#
DME_1_MVP<-subset.data.frame(DME_1_MVP,select=c("Dataset","Cell_type","BETA","SE",
                                                "P","P_bon","P_bh",
                                                "BETA_DME_1","SE_DME_1",
                                                "P_DME_1","P_bon_DME_1","P_bh_DME_1",
                                                "BETA_MVP","SE_MVP",
                                                "P_MVP","P_bon_MVP","P_bh_MVP",
                                                "P_Zhang","P_bon_Zhang","P_bh_Zhang",
                                                "n_gene_set_Zhang","n_overlap_Zhang","genes_Zhang"))
DME_1_MVP<-DME_1_MVP[order(DME_1_MVP$P),]
#
# Write a table with only significant results (BH)
#
step3<-fread("FUMA/DME_1_MVP/DropViz_L2_BH/step1_2_summary.txt")
step3<-subset.data.frame(step3,step3==1)
index<-which(DME_1_MVP$Cell_type%in%step3$Cell_type&DME_1_MVP$Dataset%in%step3$Dataset)
dt<-DME_1_MVP[index,]
dt[, P_bon_rep_Zhang := pmin(P_Zhang * nrow(dt), 1)] 
dt<-Annot_DV(dt) # Add DropViz annotations
write.table(dt,"Table_DropViz_BH.csv",sep=",",row.names=F,quote=F)
#
# Write a table with only significant results (Bonferroni + Step2)
#
step3<-fread("FUMA/DME_1_MVP/DropViz_L2/step1_2_summary.txt")
step3<-subset.data.frame(step3,step3==1)
index<-which(DME_1_MVP$Cell_type%in%step3$Cell_type&DME_1_MVP$Dataset%in%step3$Dataset)
dt<-DME_1_MVP[index,]
dt[, P_bon_rep_Zhang := pmin(P_Zhang * nrow(dt), 1)] 
dt<-Annot_DV(dt) # Add DropViz annotations
write.table(dt,"Table_DropViz_Bon.csv",sep=",",row.names=F,quote=F)
#
# Load existing workbook
#
wb<-loadWorkbook(paste0(wb_name,".xlsx"))
sheet_name<-"S5 DropViz L2 Analysis" # name the sheet
myresult<-DME_1_MVP
#
addWorksheet(wb,sheet_name)
writeData(wb,sheet_name,myresult)
addFilter(wb,sheet_name,row=1,cols=1:ncol(myresult))
addStyle(wb,sheet_name,
         createStyle(textDecoration ="bold",border="Bottom"),
         rows=1,cols=1:ncol(myresult))
setColWidths(wb,sheet_name,cols=1:ncol(myresult),widths="auto")
freezePane(wb,sheet_name,firstRow=T)
saveWorkbook(wb,paste0(wb_name,".xlsx"),overwrite=T)
#
#-------------------------------------------------------------------------------
# Supplementary table: Siletti_Seeker L2
#-------------------------------------------------------------------------------
#
DME_1_MVP<-fread("FUMA/DME_1_MVP/Siletti_Seeker_L2/magma_celltype_step1.txt")
DME_1<-fread("FUMA/DME_1/Siletti_Seeker_L2/magma_celltype_step1.txt")
MVP<-fread("FUMA/MVP/Siletti_Seeker_L2/magma_celltype_step1.txt")
Zhang<-fread("Replication/Zhang_to_Replication_Siletti_Seeker_L2.csv")
#
# Are there duplicated lines? All data tables?
#
setDT(DME_1_MVP); DME_1_MVP <- unique(DME_1_MVP)
setDT(DME_1);     DME_1     <- unique(DME_1)
setDT(MVP);       MVP       <- unique(MVP)
setDT(Zhang);     Zhang     <- unique(Zhang)
#
# Calculate P_bh and P_bon 
#
DME_1_MVP[, P_bon := pmin(P * nrow(DME_1_MVP), 1)]
DME_1_MVP[, P_bh := p.adjust(P, method = "BH")]
DME_1[, P_bon := pmin(P * nrow(DME_1), 1)]
DME_1[, P_bh := p.adjust(P, method = "BH")]
MVP[, P_bon := pmin(P * nrow(MVP), 1)]
MVP[, P_bh := p.adjust(P, method = "BH")]
Zhang[, P_bon := pmin(P * nrow(Zhang), 1)]
Zhang[, P_bh := p.adjust(P, method = "BH")]
#
DME_1_MVP <- merge(DME_1_MVP, 
                   DME_1[, .(Cell_type, Dataset, 
                             BETA_DME_1 = BETA, SE_DME_1 = SE, P_DME_1 = P, 
                             P_bon_DME_1 = P_bon, P_bh_DME_1 = P_bh)], 
                   by = c("Cell_type", "Dataset"), all.x = TRUE)

DME_1_MVP <- merge(DME_1_MVP, 
                   MVP[, .(Cell_type, Dataset, 
                           BETA_MVP = BETA, SE_MVP = SE, P_MVP = P, 
                           P_bon_MVP = P_bon, P_bh_MVP = P_bh)], 
                   by = c("Cell_type", "Dataset"), all.x = TRUE)

DME_1_MVP <- merge(DME_1_MVP, 
                   Zhang[, .(Cell_type, Dataset, 
                             P_Zhang = P, P_bon_Zhang = P_bon, P_bh_Zhang = P_bh, 
                             n_gene_set_Zhang = n_genes_set, n_overlap_Zhang = n_overlap, 
                             genes_Zhang = genes)], 
                   by = c("Cell_type", "Dataset"), all.x = TRUE)
#
DME_1_MVP<-subset.data.frame(DME_1_MVP,select=c("Dataset","Cell_type","BETA","SE",
                                                "P","P_bon","P_bh",
                                                "BETA_DME_1","SE_DME_1",
                                                "P_DME_1","P_bon_DME_1","P_bh_DME_1",
                                                "BETA_MVP","SE_MVP",
                                                "P_MVP","P_bon_MVP","P_bh_MVP",
                                                "P_Zhang","P_bon_Zhang","P_bh_Zhang",
                                                "n_gene_set_Zhang","n_overlap_Zhang","genes_Zhang"))
DME_1_MVP<-DME_1_MVP[order(DME_1_MVP$P),]
#
# Write a table with only significant results
#
step3<-fread("FUMA/DME_1_MVP/Siletti_Seeker_L2/step1_2_summary.txt")
step3<-subset.data.frame(step3,step3==1)
index<-which(DME_1_MVP$Cell_type%in%step3$Cell_type&DME_1_MVP$Dataset%in%step3$Dataset)
dt<-DME_1_MVP[index,]
dt<-subset.data.frame(dt,P_bh<=0.05)
dt[,P_bon_rep_Zhang:=pmin(P_Zhang*nrow(dt),1)]  
write.table(dt,"Table_Siletti_Seeker_BH.csv",sep=",",row.names=F,quote=F)
#
dt<-subset.data.frame(DME_1_MVP,P_bon<=0.05)
dt[, P_bon_rep_Zhang := pmin(P_Zhang * nrow(dt), 1)]  
write.table(dt,"Table_Siletti_Seeker_Bon.csv",sep=",",row.names=F,quote=F)
#
# Load existing workbook
#
wb<-loadWorkbook(paste0(wb_name,".xlsx"))
sheet_name<-"S6 Siletti Seeker L2 Analysis" # name the sheet
myresult<-DME_1_MVP
#
addWorksheet(wb,sheet_name)
writeData(wb,sheet_name,myresult)
addFilter(wb,sheet_name,row=1,cols=1:ncol(myresult))
addStyle(wb,sheet_name,
         createStyle(textDecoration ="bold",border="Bottom"),
         rows=1,cols=1:ncol(myresult))
setColWidths(wb,sheet_name,cols=1:ncol(myresult),widths="auto")
freezePane(wb,sheet_name,firstRow=T)
saveWorkbook(wb,paste0(wb_name,".xlsx"),overwrite=T)

