#Taking Rick's Code and putting it in an Rscript

#load libaries
library(Matrix)
library(irlba)
library(threejs)
library(tidyverse)
library(lfa)
library(furrr)
library(glm2)
library(cowplot)
library(drake)
options(future.globals.maxSize= 10000000000)
future::plan(multicore)

#For running in terminal
args = commandArgs(trailingOnly=T)
Chr <- args[1]
#Chr <- '22'
#set directories and paths
GenoPath='/Users/luke/genomes/genomes/hg19/phase3/'
VCFname = '.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
Dir = '/Users/luke/Documents/QualityPaper/rick/'
setwd(Dir)

#Get ID names
if(!file.exists("ids.txt")){
ids <- readLines(pipe(paste("gzcat ",GenoPath,"ALL.chr",Chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz  |
                      sed -n /^#CHROM/p | 
                      tr '\t' '\n' | 
                      tail -n +10",sep='')))

writeLines(ids, "ids.txt")}

if(file.exists("ids.txt")){
  ids <- unlist(read.table("ids.txt"))
  head(ids)
}

#Get Quality Scores
if(!file.exists("bas_dt")){
read_basfiles <- function(bas_filepath) {
  
  bas_file <- read_tsv(bas_filepath,
                  col_types = cols_only('sample' = 'c',        
                                       'average_quality_of_mapped_bases' = col_guess())) 
  return( bas_file)
}

bas_files <- list.files(path = "/Users/luke/genomes/bas_file/1000GenomesBAS/",
                        pattern="*bam.bas", full.names = TRUE)

bas_dt <- 
  map_dfr(bas_files,read_basfiles ) %>%
  rename(avg_qual="average_quality_of_mapped_bases") %>%
  group_by(sample) %>% summarize(avg_qual = mean(avg_qual))

#keep samples in 1kGP VCF files
bas_dt <- bas_dt[which((bas_dt$sample %in% ids)==T),]

#get populations of each ID
ped <- 
  read_tsv(url("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped"))%>% 
  rename(sample="Individual ID")

bas_dt <- inner_join( ped  %>% 
                       select(sample, Population),
                       bas_dt,
                     by="sample")
save(bas_dt, file = '~/Documents/QualityPaper/rick/bas_dt')
}
if(file.exists("bas_dt")){
load('~/Documents/QualityPaper/rick/bas_dt')
}
                   
#Load Genotype Data
#Load VCF file into a sparse matrix
to_datalist <- function(i, data_genotype, data_iaf){
  
    data <- tibble(genotype=data_genotype[i,],
               iaf=data_iaf[i,])
               
     return(data)       
}
get_data <- function(){
  
  p <- pipe(paste("gzcat ",GenoPath,"ALL.chr",Chr, VCFname," |
          sed /^#/d  |
          cut  -f '10-' |
          ./a.out | 
          cut -f '1-2'",sep=''))
  x <- read.table(p,
                colClasses=c("integer","integer"), 
                fill=TRUE,
                row.names=NULL)

  
  data <- sparseMatrix(i=x[,2], j=x[,1], x=1L)
  
  
  return(data)
  
}
                     
get_af_data <- function(data, d=4){
  
  # compute the first d factors including the intercept
  LF <- lfa(data, d)
  # compute the individual-specific allele frequencies
  data_iaf <- af(data, LF)
  loci <- c(1:dim(data)[1])
  names(loci) <- loci
  
  data_list <- future_map(loci, 
                 to_datalist,
                 data_genotype=data,
                 data_iaf=data_iaf) 
  
  data_list <- 
  enframe(data_list,
          name = "locus", 
          value = "dataset")
  
  return(data_list)
  
}

get_predicted_mean <- function(model_fit) {
 predicted <- 
   augment(model_fit, type.predict = "response") %>%
   select("genotype","offset.iaf.",".fitted") %>%
   rename(iaf="offset.iaf.", mu=".fitted" )
  return(predicted)
}

get_stats <- function(model_fit){
  stats_dt <- 
    tidy(model_fit)  %>%
    filter(term=="avg_qual_mean") %>%
    select(estimate, p.value) %>%
    mutate(deviance= anova(model_fit)["Deviance"][[1]][2])
  return(stats_dt)
}

fit_model <- function(dataset){
  
   fit <- glm(genotype ~ avg_qual + offset(iaf),
                           family ="binomial",
                           data=dataset)
   
   stats <- get_stats(fit)
   predicted <- get_predicted_mean(fit)
   
   return(list(stats, predicted))
}

get_plot <- function(data){
  
 p <- ggplot(bind_cols(data, 
                   bas_dt %>% 
                     select(Population,
                            avg_qual_mean))) +
    geom_point(aes(avg_qual_mean, 
                   genotype),
               size=.5,
               alpha= 0.4) +
    geom_line(aes(avg_qual_mean,
                  iaf, 
                  colour="af"),
              linetype="dashed") +
    geom_line(aes(avg_qual_mean,
                  mu,
                  colour="mu")) +
      theme_bw() +
    facet_wrap(~Population)
  
  return(p)
}

run_workflow <- function(loci){

  data <- get_data()
  data <- t(as.matrix(data[,loci]))
  data_list <- get_af_data(data) 
  
  rm(data)
  gc()
  
  data_list <- 
    data_list %>%
  #slice(1:10) %>%
  transmute(
    locus=locus,
    fit1 = future_map(dataset, fit_model))  %>%
  transmute(
    locus=locus,
    stats = future_map(fit1, c(1)),
    predicted =future_map(fit1, c(2))) %>%
    unnest(stats) %>%
    arrange(-deviance) 
  
  return(data_list)
  
}

data_list <- run_workflow (loci=c(1:50000))

save(data_list, file = paste('~/genomes/genomes/hg19/DataListsDataLists/Chr_',Chr,'_DataList.Rdata', sep=''))
