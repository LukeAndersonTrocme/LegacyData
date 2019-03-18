## ------------------------------------------------------------------------
library(tidyverse)
library(lfa)
library(furrr)
library(glm2)
library(tictoc)
library(data.table)
options(future.globals.maxSize= 1000000000000)
plan(multicore)

## ------------------------------------------------------------------------

args <- commandArgs(trailingOnly=T)

path <- args[1]
setwd(path)

plinkPath <- args[2] 

chromosome <- args[3]

data_genotype <- read.bed(paste0(plinkPath,'chr', chromosome))

message("plink data loaded")

## ------------------------------------------------------------------------
data_bim <- fread(paste0(plinkPath,'chr', chromosome,'.bim'), 
                  col.names = c('Chr','rsID','X',
                                'Pos','Ref','Alt'))

rownames(data_genotype) <- data_bim$rsID

## ----lfa-----------------------------------------------------------------
tic()
#LF <- lfa(data_genotype, 5)
toc()
message("lfa run complete")

## ----save----------------------------------------------------------------
tic()
#saveRDS(LF, "./data/preprocessed/lf5.rds")
toc()

## ----read lf-------------------------------------------------------------
LF <- readRDS( "./data/preprocessed/lf5.rds")

## ------------------------------------------------------------------------
#this runs the logistic regression
lreg <- function(y, X){
  X <- rbind(X, X)
  y1 <- as.numeric((y==1) | (y==2))
  y2 <- as.numeric(y==2)
  y <- c(y1,y2)
  b <- glm2(cbind(y, 2-y) ~ -1 + X, family="binomial")$coef
  
  return(b)
  
}
#modified from gcat
fit_glm <- function(..., LF, trait){
  snp <- c(...)
  ind <- !is.na(snp) #logical index vector
  snp_no_na <- snp[ind]
  LF_no_na <- LF[ind,]
  p0 <- rep(NA, length(snp))
  p1 <- rep(NA, length(snp))
  trait_no_na <- trait[ind]
  b0 <- lreg(snp_no_na, LF_no_na) #coefficients from logreg 
  b1 <- lreg(snp_no_na, cbind(LF_no_na, trait_no_na))
  est0 <- .Call("mv", LF_no_na, b0)
  est1 <- .Call("mv", cbind(LF_no_na, trait_no_na), b1)
  p0[ind] <- exp(est0)/(1+exp(est0))
  p1[ind] <- exp(est1)/(1+exp(est1))
  
  dataset <-
    tibble(snp=snp,
           p0 = p0,
           p1 = p1)
  
  stat <-
    dataset %>%
    summarize(
      dev0 = -2 * sum(snp * log(p0) +
                        (2 - snp) * log(1 - p0),
                      na.rm = TRUE),
      dev1 = -2 * sum(snp * log(p1) +
                        (2 - snp) * log(1 - p1),
                      na.rm = TRUE),
      devdiff = dev0 - dev1,
      dev_logp = -log10(pchisq(devdiff,
                               1,
                               lower.tail = FALSE))
    ) %>% 
    mutate(num_nas = sum(!ind),
	   num_zeros= sum(snp_no_na ==0),
           b=b1[length(b1)])
  
  return(list(dataset=dataset, stat=stat))
         
}

## ------------------------------------------------------------------------
gc()

## ------------------------------------------------------------------------
fam_file <-fread(paste0(plinkPath,'chr', chromosome,'.fam'), 
                  col.names=c('sample','a','b','c','d','e'))
                 
bas_filepath <- "./data/metadata/average_quality_dt.txt"
bas_dt <- read_tsv(bas_filepath)
bas_dt <- merge(fam_file, 
                bas_dt,
                by = 'sample') 
avg_qual_mean <- bas_dt$avg_qual_mean

## ------------------------------------------------------------------------
head(avg_qual_mean)

## ------------------------------------------------------------------------
message("fitting model")
tic("fitting model")
labels <- rownames(data_genotype)
data_fit <- 
	data_genotype %>%
	as_tibble() %>%
	transmute(
	fit = future_pmap(., fit_glm,
	LF=LF, 
	trait=avg_qual_mean),
		rsID=labels)  %>%
	unnest(fit) %>%
	mutate(l= rep(c("data","stat"), n()/2))%>%
	spread(l, fit)%>%
	unnest(stat)
   toc()


## ------------------------------------------------------------------------
message("Saving Data")
saveRDS(data_fit, 
	file.path("./data/preprocessed/subsets",
	paste0("fits_chrom_", chromosome, ".rds")))

write.table(data_fit[-2], 
              file.path("./data/preprocessed/subsets", 
                        paste0("fits_chrom_", chromosome, ".txt")),
              quote = FALSE,
              row.names = FALSE)

## ------------------------------------------------------------------------
gc()

## ------------------------------------------------------------------------
sessionInfo()

