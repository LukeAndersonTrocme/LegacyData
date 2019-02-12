#Make QQ
library(data.table)
library(ggplot2)

Repeat_Indels<-fread('/Users/luke/Documents/QualityPaper/Misc/Repeat_Indels.csv')
NORepeat_Indels <-fread('/Users/luke/Documents/QualityPaper/Misc/NORepeat_Indels.csv')
Repeat_SNPs <-fread('/Users/luke/Documents/QualityPaper/Misc/Repeat_SNPs.csv')
NORepeat_SNPs <-fread('/Users/luke/Documents/QualityPaper/Misc/NORepeat_SNPs.csv')

makeQQ <- function(pvalues, name){
set.seed(1234)
pvalue <- runif(length(pvalues), min=0, max=1)
inflation <- function(ps) {
  chisq <- qchisq(1 - pvalue, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}
inflation(pvalue)

gg_qqplot(pvalues, name) +
  theme_bw(base_size = 24) +
  annotate(
    geom = "text",
    x = -Inf,
    y = Inf,
    hjust = -0.15,
    vjust = 1 + 0.15 * 3,
    label = sprintf("λ = %.2f", inflation(ps)),
    size = 8
  ) + ggtitle(name) +
  theme(
    axis.ticks = element_line(size = 0.5),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
 }
 ##modified from : https://slowkow.com/notes/ggplot2-qqplot/
 gg_qqplot <- function(ps,name) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)))
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}

makeQQ(Repeat_Indels$p, 'Repeat Indels')
ggsave('~/Documents/QualityPaper/Figures/Repeat_Indels_QQ.jpg',height=8,width=12)
makeQQ(NORepeat_Indels$p, 'Non Repeat Indels')
ggsave('~/Documents/QualityPaper/Figures/NORepeat_Indels_QQ.jpg',height=8,width=12)
makeQQ(Repeat_SNPs$p, 'Repeat SNPs')
ggsave('~/Documents/QualityPaper/Figures/Repeat_SNPs_QQ.jpg',height=8,width=12)
makeQQ(NORepeat_SNPs $p, 'Non Repeat SNPs')
ggsave('~/Documents/QualityPaper/Figures/NORepeat_SNPs_QQ.jpg',height=8,width=12)
