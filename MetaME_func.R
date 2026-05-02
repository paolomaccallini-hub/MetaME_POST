# file name: MetaME_func
#
library(httr) 
library(R.utils) 
library(data.table)
library(MungeSumstats)
library(yaml) 
library(ggplot2) 
library(scattermore)
library(patchwork)
#
#-------------------------------------------------------------------------------
# Add a data folder, if absent
#-------------------------------------------------------------------------------
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Data")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
#-------------------------------------------------------------------------------
# This function works on DecodeME summary statistics
#-------------------------------------------------------------------------------
#
DME_func<-function(sumstats) {
  #
  #-------------------------------------------------------------------------------
  # Build data base for DecodeME
  #-------------------------------------------------------------------------------
  #
  current_dir<-getwd()
  folder_path<-file.path(current_dir,"Data/DecodeME")  
  if(!dir.exists(folder_path)) {
    dir.create(folder_path) 
  }
  #
  # Download summary statistics, if not present
  #
  url<-"https://files.de-1.osf.io/v1/resources/rgqs3/providers/osfstorage/688881b3c8c15886a1e7cc93/?zip="
  file_path<-file.path("Data/DecodeME/DecodeME_summary.zip")
  if(!file.exists(file_path)) {
    print("Downloading the summary statistics...")
    RETRY(
      verb = "GET",
      url = url,
      write_disk(file_path, overwrite = TRUE),
      times = 5,           # up to 5 attempts
      pause_min = 5,       # wait 5s between attempts
      terminate_on = c(404,403) # don't retry on these errors
    )
  }   
  #
  #---------------------------------------------------------------------------
  # Munge
  #---------------------------------------------------------------------------
  #
  # Path to munged file
  #
  munge_path<-paste0("Munged/DME_",sumstats,"_GRCh38.tsv.gz")
  #
  if (!file.exists(munge_path)) {
    #
    # Read variants that passed quality filter
    #
    input_phenotype<-"gwas_qced.var.gz"
    zip_path<-file.path("Data/DecodeME","DecodeME_summary.zip")
    tmpdir<-tempdir()
    unzip(zip_path,files=input_phenotype,exdir=tmpdir)
    myQCEDvariants<-fread(file.path(tmpdir,input_phenotype))
    #
    # Read summary statistics
    #
    input_phenotype<-paste0("gwas_",sumstats,".regenie.gz")
    zip_path<-file.path("Data/DecodeME","DecodeME_summary.zip")
    tmpdir<-tempdir()
    unzip(zip_path,files=input_phenotype,exdir=tmpdir)
    mydata<-fread(file.path(tmpdir,input_phenotype))
    head(mydata)
    colnames(mydata)
    #
    # Remove columns we wont use
    #
    mydata<-mydata[,-"EXTRA"]
    mydata<-mydata[,-"A1FREQ_CASES"]
    mydata<-mydata[,-"A1FREQ_CONTROLS"]
    mydata<-mydata[,-"TEST"]
    mydata<-mydata[,-"CHISQ"]
    #
    # Keep only variants that passed INFO quality filter
    #
    mydata<-mydata[mydata$ID %in% myQCEDvariants$ID, ]
    remove(myQCEDvariants)
    #
    # Filter by MAF (use cut-off of uncommon variants)
    #
    mydata<-mydata[mydata$A1FREQ<=1-MAFco_uc,]
    mydata<-mydata[mydata$A1FREQ>=MAFco_uc,]
    #
    # Specify effect allele: in format_sumstats A1 is the non-effect allele, A2 is the effect allele,
    # FRQ is the frequency of the effect allele 
    # (https://www.bioconductor.org/packages/devel/bioc/vignettes/MungeSumstats/inst/doc/MungeSumstats.html).
    # In DecodeME, ALLELE0 is	Non effect alleles (reference), ALLELE1	is Effect alleles (alternate),
    # A1FREQ is	Effect allele frequency (from DecodeME readme)
    #
    for (j in 1:ncol(mydata)) {
      if (colnames(mydata)[j]=="ALLELE0") colnames(mydata)[j]<-"other_allele"
      if (colnames(mydata)[j]=="ALLELE1") colnames(mydata)[j]<-"effect_allele"
      if (colnames(mydata)[j]=="A1FREQ") colnames(mydata)[j]<-"effect_allele_frequency"
    }
    #
    # Munge 
    #
    format_sumstats(mydata,ref_genome="GRCh38",
                    compute_z="BETA",
                    save_path=munge_path, # where to save the sumstat
                    bi_allelic_filter=F, # keep non-biallelic SNPs
                    flip_frq_as_biallelic=T, # keep non-biallelic SNPs even when they need flipping 
                    log_mungesumstats_msgs=T, # keep a log of all the messages and warnings
                    log_folder="Munged") # were to save the log
    #
    # Check for errors
    # 
    mymunged<-fread(munge_path)
    test<-0
    while(test<=30) {
      n1<-sample(c(1:nrow(mymunged)),1)
      n2<-which(mydata$CHROM==mymunged$CHR[n1]&mydata$GENPOS==mymunged$BP[n1])
      if (length(n2)==1) {
        if (((mydata[n2,"effect_allele",with=F]==mymunged[n1,"A2",with=F])&
             abs(mydata[n2,"BETA",with=F]-mymunged[n1,"BETA",with=F])<1e-6)|
            ((mydata[n2,"effect_allele",with=F]==mymunged[n1,"A1",with=F])&
             abs(mydata[n2,"BETA",with=F]+mymunged[n1,"BETA",with=F])<1e-6)) {
          test<-test+1
          next
        } else {
          stop("Munging generated an error!")
        }  
      }
    }
    remove(mymunged,mydata)
  } 
  if (T) {
    #
    # Lift over to GRCh37 
    #
    munge_path<-paste0("Munged/DME_",sumstats,"_GRCh37.tsv.gz")
    #
    if (!file.exists(munge_path)) {
      #
      # Read GRCh38
      #
      munge_path<-paste0("Munged/DME_",sumstats,"_GRCh38.tsv.gz")
      mydata<-fread(munge_path)
      #
      # Assign GRCh37 path
      #
      munge_path<-paste0("Munged/DME_",sumstats,"_GRCh37.tsv.gz")
      #
      format_sumstats(mydata,
                      ref_genome="GRCh38",
                      convert_ref_genome="GRCh37",
                      save_path=munge_path, # where to save the sumstat
                      bi_allelic_filter=F, # keep non-biallelic SNPs
                      flip_frq_as_biallelic=T, # keep non-biallelic SNPs even when they need flipping 
                      log_mungesumstats_msgs=T, # keep a log of all the messages and warnings
                      log_folder="Munged") # were to save the log
      #
      # Check for errors
      # 
      mydata<-fread(paste0("Munged/DME_",sumstats,"_GRCh38.tsv.gz"))
      mymunged<-fread(paste0("Munged/DME_",sumstats,"_GRCh37.tsv.gz"))
      test<-0
      while(test<=30) {
        n1<-sample(c(1:nrow(mymunged)),1)
        n2<-which(mydata$SNP==mymunged$SNP[n1])
        if (length(n2)==1) {
          if (((mydata[n2,"A2",with=F]==mymunged[n1,"A2",with=F])&
               abs(mydata[n2,"BETA",with=F]-mymunged[n1,"BETA",with=F])<1e-6)|
              ((mydata[n2,"A2",with=F]==mymunged[n1,"A1",with=F])&
               abs(mydata[n2,"BETA",with=F]+mymunged[n1,"BETA",with=F])<1e-6)) {
            test<-test+1
            next
          } else {
            stop("Munging generated an error!")
          }  
        }
      }
      remove(mymunged,mydata)  
      gc() # free unused memory
    } 
  }
}
#
#-------------------------------------------------------------------------------
# This function works on UK Biobank summary statistics
#-------------------------------------------------------------------------------
#
UKBNL_func<-function(sumstats) {
  #
  #-------------------------------------------------------------------------------
  # Build database for UK Biobank (Neale Lab)
  #-------------------------------------------------------------------------------
  #
  current_dir<-getwd()
  folder_path<-file.path(current_dir,"Data/NealeLab")  
  if(!dir.exists(folder_path)) {
    dir.create(folder_path) 
  }
  #
  # Download phenotype data, version 2, if not present
  #
  for (sex in c("both_sexes","female","male")) {
    url<-paste0("https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/phenotypes.",sex,".v2.tsv.bgz")
    destfile<-paste0("phenotypes.",sex,".v2.tsv.bgz")
    file_path<-file.path(current_dir,"Data/NealeLab/",destfile)
    if(!file.exists(file_path)) {
      print("Downloading the list of phenotypes...")
      RETRY(
        verb = "GET",
        url = url,
        write_disk(file_path, overwrite = TRUE),
        times = 5,           # up to 5 attempts
        pause_min = 5,       # wait 5s between attempts
        terminate_on = c(404,403) # don't retry on these errors
      )
    }   
  }
  #
  # Download complete list of imputed variants, if not present
  #
  url<-"https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz"
  destfile<-"variants.tsv.bgz"
  file_path<-file.path(current_dir,"Data/NealeLab",destfile)
  if(!file.exists(file_path)) {
    print("Downloading and editing the full set of imputed variants from Neale Lab (UK Biobank)")
    RETRY(
      verb = "GET",
      url = url,
      write_disk(file_path, overwrite = TRUE),
      times = 5,           # up to 5 attempts
      pause_min = 5,       # wait 5s between attempts
      terminate_on = c(404, 403) # don't retry on these errors
    ) 
    #
    # Edit the file keeping only the columns: myvariants and p_hwe (this takes a while)
    #
    all_variants<-fread(file_path,header="auto",sep="\t")
    all_variants<-all_variants[,c("variant","rsid","p_hwe","info")]
    write.table(all_variants,file=file_path,sep="\t",col.names=T,row.names=F)
  }  
  #
  # Download summary statistics for the selected phenotypes
  #
  for (phenotype in c("20002_1482")) {
    for (sex in c("both_sexes","male","female")) {
      url<-paste0("https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/",phenotype,".gwas.imputed_v3.",sex,".tsv.bgz")
      destfile<-paste0(phenotype,".gwas.imputed_v3.",sex,".tsv.bgz")
      file_path<-file.path(current_dir,"Data/NealeLab/",destfile)
      if(!file.exists(file_path)) {
        print("Downloading summary statistics")
        RETRY(
          verb = "GET",
          url = url,
          write_disk(file_path, overwrite = TRUE),
          times = 5,           # up to 5 attempts
          pause_min = 5,       # wait 5s between attempts
          terminate_on = c(404, 403) # don't retry on these errors
        )
      }  
    }
  } 
  gc() # free unused memory
  #
  #---------------------------------------------------------------------------
  # Munge to GRCh38
  #---------------------------------------------------------------------------
  #
  # Path to munged file
  #
  munge_path<-paste0("Munged/UKBNL_",sumstats,"_GRCh38.tsv.gz")
  #
  if (!file.exists(munge_path)) {
    #
    # Read summary statistics
    #
    file_name<-paste0("Data/NealeLab/20002_1482.gwas.imputed_v3.",sumstats,".tsv.bgz") 
    mydata<-fread(file_name)
    head(mydata)
    colnames(mydata)
    #
    # Read variant annotations
    #
    all_variants<-fread("Data/NealeLab/variants.tsv.bgz",header=TRUE,sep="\t")
    #
    # Remove columns we wont use
    #
    mydata<-mydata[,-"AC"]
    mydata<-mydata[,-"ytx"]
    mydata<-mydata[,-"expected_case_minor_AC"]
    #
    # Add annotation
    #
    mydata<-merge(mydata,all_variants,by="variant",all.x=T)
    remove(all_variants)
    #
    # Filter by MAF
    #
    mydata<-mydata[mydata$minor_AF>=MAFco_uc,]
    #
    # Filter by info
    #
    mydata<-mydata[mydata$info>=INFOco,]
    mydata<-mydata[,-"info"]
    #
    # Remove low confidence variants
    #
    mydata<-mydata[low_confidence_variant==F,]
    mydata<-mydata[,-"low_confidence_variant"]
    #
    # Remove variants outside HWE
    #
    mydata<-mydata[p_hwe>pco_HWE]
    mydata<-mydata[,-"p_hwe"]
    #
    # Split columns with coordinates and alleles
    # Note that in this sumstat the alternative allele is the effect allele:
    # https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?gid=227859291#gid=227859291
    #
    mydata[,c("CHR","BP","other_allele","effect_allele"):=tstrsplit(variant,":",fixed=T)]
    mydata[,BP:=as.numeric(BP)]
    mydata<-mydata[,-"variant"]
    #
    # Add N
    #
    phenotypes<-fread(paste0("Data/NealeLab/phenotypes.",sumstats,".v2.tsv.bgz"))
    phenotypes<-phenotypes[phenotype=="20002_1482"]
    N_cases<-phenotypes$n_cases
    N_controls<-phenotypes$n_controls
    N<-phenotypes$n_non_missing
    mydata$N<-rep(N,nrow(mydata))
    mydata$N_cases<-rep(N_cases,nrow(mydata))
    mydata$N_controls<-rep(N_controls,nrow(mydata))
    #
    # Munge 
    #
    format_sumstats(mydata,ref_genome="GRCh37",
                    convert_ref_genome="GRCh38",
                    compute_z="BETA",
                    bi_allelic_filter=F,
                    flip_frq_as_biallelic=T,
                    save_path=munge_path,
                    log_mungesumstats_msgs=T, # keep a log of all the messages and warnings
                    log_folder="Munged") # were to save the log
    #
    # Check for errors
    # 
    mymunged<-fread(munge_path)
    test<-0
    while(test<=30) {
      n1<-sample(c(1:nrow(mymunged)),1)
      n2<-which(mydata$rsid==mymunged$SNP[n1])
      if (length(n2)==1) {
        if (((mydata[n2,"effect_allele",with=F]==mymunged[n1,"A2",with=F])&
             abs(mydata[n2,"beta",with=F]-mymunged[n1,"BETA",with=F])<1e-6)|
            ((mydata[n2,"effect_allele",with=F]==mymunged[n1,"A1",with=F])&
             abs(mydata[n2,"beta",with=F]+mymunged[n1,"BETA",with=F])<1e-6)) {
          test<-test+1
          next
        } else {
          stop("Munging generated an error!")
        }  
      }
    }
    remove(mydata,mymunged)
    gc() # free unused memory
  } 
  #
  #---------------------------------------------------------------------------
  # Munge to GRCh37
  #---------------------------------------------------------------------------
  #
  # Path to munged file
  #
  munge_path<-paste0("Munged/UKBNL_",sumstats,"_GRCh37.tsv.gz")
  #
  if (!file.exists(munge_path)) {
    #
    # Read summary statistics
    #
    file_name<-paste0("Data/NealeLab/20002_1482.gwas.imputed_v3.",sumstats,".tsv.bgz") 
    mydata<-fread(file_name)
    head(mydata)
    colnames(mydata)
    #
    # Read variant annotations
    #
    all_variants<-fread("Data/NealeLab/variants.tsv.bgz",header=TRUE,sep="\t")
    #
    # Remove columns we wont use
    #
    mydata<-mydata[,-"AC"]
    mydata<-mydata[,-"ytx"]
    mydata<-mydata[,-"expected_case_minor_AC"]
    #
    # Add annotation
    #
    mydata<-merge(mydata,all_variants,by="variant",all.x=T)
    remove(all_variants)
    #
    # Filter by MAF
    #
    mydata<-mydata[mydata$minor_AF>=MAFco_uc,]
    #
    # Filter by info
    #
    mydata<-mydata[mydata$info>=INFOco,]
    mydata<-mydata[,-"info"]
    #
    # Remove low confidence variants
    #
    mydata<-mydata[low_confidence_variant==F,]
    mydata<-mydata[,-"low_confidence_variant"]
    #
    # Remove variants outside HWE
    #
    mydata<-mydata[p_hwe>pco_HWE]
    mydata<-mydata[,-"p_hwe"]
    #
    # Split columns with coordinates and alleles
    # Note that in this sumstat the alternative allele is the effect allele:
    # https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?gid=227859291#gid=227859291
    #
    mydata[,c("CHR","BP","other_allele","effect_allele"):=tstrsplit(variant,":",fixed=T)]
    mydata[,BP:=as.numeric(BP)]
    mydata<-mydata[,-"variant"]
    #
    # Add N
    #
    phenotypes<-fread(paste0("Data/NealeLab/phenotypes.",sumstats,".v2.tsv.bgz"))
    phenotypes<-phenotypes[phenotype=="20002_1482"]
    N_cases<-phenotypes$n_cases
    N_controls<-phenotypes$n_controls
    N<-phenotypes$n_non_missing
    mydata$N<-rep(N,nrow(mydata))
    mydata$N_cases<-rep(N_cases,nrow(mydata))
    mydata$N_controls<-rep(N_controls,nrow(mydata))
    #
    # Munge 
    #
    format_sumstats(mydata,ref_genome="GRCh37",
                    convert_ref_genome="GRCh37",
                    compute_z="BETA",
                    bi_allelic_filter=F,
                    flip_frq_as_biallelic=T,
                    save_path=munge_path,
                    log_mungesumstats_msgs=T, # keep a log of all the messages and warnings
                    log_folder="Munged") # were to save the log
    #
    # Check for errors
    # 
    mymunged<-fread(munge_path)
    test<-0
    while(test<=30) {
      n1<-sample(c(1:nrow(mymunged)),1)
      n2<-which(mydata$rsid==mymunged$SNP[n1])
      if (length(n2)==1) {
        if (((mydata[n2,"effect_allele",with=F]==mymunged[n1,"A2",with=F])&
             abs(mydata[n2,"beta",with=F]-mymunged[n1,"BETA",with=F])<1e-6)|
            ((mydata[n2,"effect_allele",with=F]==mymunged[n1,"A1",with=F])&
             abs(mydata[n2,"beta",with=F]+mymunged[n1,"BETA",with=F])<1e-6)) {
          test<-test+1
          next
        } else {
          stop("Munging generated an error!")
        }  
      }
    }
    remove(mydata,mymunged)
    gc() # free unused memory
  } 
}
#
#-------------------------------------------------------------------------------
# This function retrieves total number of variants and number of genome-wide 
# significant variants from MungeSumstats output
#-------------------------------------------------------------------------------
#
parse_munge_log <- function(log_file) {
  txt <- readLines(log_file, warn=FALSE)
  txt <- gsub("\033\\[[^m]*m|\033G[0-9]*;|\033g", "", txt)
  txt <- gsub("\r", "", txt)
  #
  # Find lines from the second (post-munging) summary report
  # The second report follows the "Writing in tabular format" line
  #
  write_line <- grep("Writing in tabular format", txt)
  if (length(write_line)==0) return(list(rows=NA_integer_, gws=NA_integer_))
  txt_after <- txt[write_line:length(txt)]
  #
  rows_line <- grep("[0-9,]+ rows \\(", txt_after, value=TRUE)[1]
  gws_line  <- grep("[0-9,]+ genome-wide", txt_after, value=TRUE)[1]
  #
  rows_after <- as.integer(gsub(",","", regmatches(rows_line, regexpr("[0-9,]+", rows_line))))
  gws_after  <- as.integer(gsub(",","", regmatches(gws_line,  regexpr("[0-9,]+", gws_line))))
  #
  return(list(rows=rows_after, gws=gws_after))
}
#
#-------------------------------------------------------------------------------
# This function calculate lambda and plot the QQ plot. It is inspired
# by a function by Matti Pirinen - Chapter 6, page 7 of his GWAS lectures
#-------------------------------------------------------------------------------
#
lambda_func<-function(mydata,name) {
  #
  # Calculate lambda 
  #
  expect.stats<-qchisq(ppoints(nrow(mydata)),df=1,lower.tail=FALSE)
  obs.stats<-qchisq(mydata$P,df=1,lower.tail=FALSE)
  lambda<-median(obs.stats)/median(expect.stats)
  #
  # Generate QQ plot for meta analysis: ChiSquared values
  #
  tiff(filename=paste0("Output/ChiSquared_QQ_",name,".tiff"),
       width=7,height=7,units="in",res=300,compression="lzw")
  #
  par(mar = c(5, 5, 4, 2))
  scattermoreplot(sort(expect.stats),sort(obs.stats),
                  xlab = expression(Expected ~ chi[1]^2),
                  ylab = expression(Observed ~ chi[1]^2),
                  cex=5,col="black",xaxt="n",yaxt="n",
                  cex.lab=1.3)
  axis(1, at=seq(0, 30, by=2))
  axis(2, at=seq(0, 50, by=5))
  abline(0,1,lwd=2,col="red",lty=2)
  legend("topleft",legend=bquote(lambda[GC]==.(signif(lambda, 3))),bty="n",cex=1.2)
  dev.off()  
  #
  # Generate QQ plot for meta analysis: -Log10P
  #
  tiff(filename=paste0("Output/Log10P_QQ_",name,".tiff"),
       width=7,height=7,units="in",res=300,compression="lzw")
  #
  par(mar = c(5, 5, 4, 2))
  scattermoreplot(sort(-log10(ppoints(nrow(mydata)))),sort(-log10(mydata$P)),
                  xlab = expression(Expected ~ -Log[10](P)),
                  ylab = expression(Observed ~ -Log[10](P)),
                  cex=5,col="black",xaxt="n",yaxt="n",
                  cex.lab=1.3)
  axis(1, at = seq(0,8,by=1))    
  axis(2, at = seq(0,12,by=1))
  abline(0,1,lwd=2,col="red",lty=2)
  legend("topleft",legend=bquote(lambda[GC]==.(signif(lambda, 3))),bty="n",cex=1.2)
  dev.off()  
  #
  # Return result
  #
  return(lambda)
}
#
#-------------------------------------------------------------------------------
# As lambda_func() but specific for meta analyses: it compares the meta gwas
# with the component sumstats
#-------------------------------------------------------------------------------
#
meta_lambda_func<-function(mydata,name) {
  #
  # retrieve composing sustats
  #
  names<-c()
  h<-1
  for (p in populations) {
    if (grepl(p,name)) {
      names[h]<-p
      h<-h+1
    }
  }
  #
  # Vectors and lists
  #
  expect.stats<-list()
  obs.stats<-list()
  lambda<-c()
  log10P<-list()
  #
  # Calculate lambda for meta analysis
  #
  expect.stats[[1]]<-qchisq(ppoints(nrow(mydata)),df=1,lower.tail=FALSE)
  obs.stats[[1]]<-qchisq(mydata$P,df=1,lower.tail=FALSE)
  lambda[1]<-median(obs.stats[[1]])/median(expect.stats[[1]])
  log10P[[1]]<--log10(mydata$P)
  #
  # Calculate lambda for composing populations
  #
  for (k in 1:length(names)) {
    #
    # Read sumstat
    #
    file_name_gz<-paste0(current_dir,"/Munged/",names[k],"_GRCh37.tsv.gz")
    mydata<-fread(file_name_gz)
    expect.stats[[k+1]]<-qchisq(ppoints(nrow(mydata)),df=1,lower.tail=FALSE)
    obs.stats[[k+1]]<-qchisq(mydata$P,df=1,lower.tail=FALSE)
    lambda[k+1]<-median(obs.stats[[k+1]])/median(expect.stats[[k+1]])
    log10P[[k+1]]<--log10(mydata$P)
  }
  #
  # Generate QQ plot for meta analysis: ChiSquared values
  #
  tiff(filename=paste0("Output/ChiSquared_QQ_",name,".tiff"),
       width=7,height=7,units="in",res=300,compression="lzw")
  #
  par(mar=c(5,5,4,2))
  x_max<-max(sapply(expect.stats, max))
  y_max<-max(sapply(obs.stats, max))
  colors<-c("red","blue","forestgreen","darkorange","purple")
  legends<-list()
  lableg<-c("META",names)
  plot(NULL,xlim=c(0,x_max),ylim=c(0,y_max),xlab=expression(Expected~chi[1]^2),
       ylab=expression(Observed~chi[1]^2),xaxt="n",yaxt="n",cex.lab=1.3)
  axis(1,at=seq(0,x_max,by=2))
  axis(2,at=seq(0,y_max,by=5))
  #
  for (k in 1:length(obs.stats)) {
    par(new=T)
    scattermoreplot(sort(expect.stats[[k]]),sort(obs.stats[[k]]),cex=8,col=colors[k],
                    xlab="",ylab="",xlim=c(0, x_max),ylim=c(0, y_max),xaxt="n",yaxt="n")
    legends[[k]] <- bquote(lambda[.(lableg[k])] == .(signif(lambda[k],4)))
  }
  abline(0,1,lwd=2,col="black",lty=2)
  legend("topleft",legend=legends,col=colors,lwd=5,bty="n",cex=1.3)
  dev.off()
  #
  # Generate QQ plot for meta analysis: -Log10(P)
  #
  tiff(filename=paste0("Output/Log10P_QQ_",name,".tiff"),
       width=7,height=7,units="in",res=300,compression="lzw")
  #
  par(mar=c(5,5,4,2))
  x_max<-8
  y_max<-max(sapply(log10P,max))
  colors<-c("red","blue","forestgreen","darkorange","purple")
  legends<-list()
  lableg<-c("META",names)
  plot(NULL,xlim=c(0,x_max),ylim=c(0,y_max),
       xlab = expression(Expected ~ -Log[10](P)),
       ylab = expression(Observed ~ -Log[10](P)),xaxt="n",yaxt="n",cex.lab=1.3)
  axis(1,at=seq(0,x_max,by=1))
  axis(2,at=seq(0,y_max,by=1))
  #
  for (k in 1:length(obs.stats)) {
    par(new=T)
    N<-length(log10P[[k]])
    scattermoreplot(sort(-log10(ppoints(N))),sort(log10P[[k]]),
                    cex=8,col=colors[k],ylab="",xlab="",
                    xlim=c(0,x_max),ylim=c(0,y_max),xaxt="n",yaxt="n")
    legends[[k]] <- bquote(lambda[.(lableg[k])] == .(signif(lambda[k],4)))
  }
  abline(0,1,lwd=2,col="black",lty=2)
  legend("topleft",legend=legends,col=colors,lwd=5,bty="n",cex=1.3)
  dev.off()
  #
  # Return result
  #
  return(lambda[1])
}

