'''Running edgeR on Hough samples
author: Carsten Knutsen
Date: March 18 2023
conda environment: edgeR
'''
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("edgeR",force=TRUE)
#BiocManager::install("biomaRt")




library(edgeR)
# read in cleaned counts
x <- read.delim("/home/carsten/alvira_bioinformatics/hough_sequencing/data/output_data/clean_counts.csv", row.names = 'symbol', sep = ',')

# run lung
group <- factor(c('AP',
                  'AP',
                  'AP',
                  'AS',
                  'AS', 
                  'AS',
                  'JP', 
                  'JP',
                  'JP',
                  'JS',
                  'JS',
                  'JS')
                )
y <- DGEList(counts=x,group=group)
keep <- filterByExpr(y, design = group)
y <- y[keep, ,keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method = 'TMM')
tmm <- cpm(y)
write.table(tmm, sep = ",", file="/home/carsten/alvira_bioinformatics/hough_sequencing/data/output_data/tmm.csv",col.names=NA)
# Run QLF test for every level
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
y <- estimateDisp(y,design)
fit <- glmQLFit(y, design)
contrasts <- makeContrasts( 
  AS_v_AP = AP-AS,
  JS_v_JP = JP-JS,
  JS_v_AS = AS-JS,
  JP_v_AP = AP-JP,
  juvenile_v_adult_treatment = (AP-AS)-(JP-JS),
  saline_v_pseudomonas_age = (AP-JP)-(AS-JS),
  levels = design )
for (i in 1:ncol(contrasts)){
  current.glmQLFTest <- glmQLFTest(fit, contrast = contrasts[,i])
  out <- topTags(current.glmQLFTest, n=Inf)
  col = colnames(contrasts)[i]
  output <- sprintf("/home/carsten/alvira_bioinformatics/hough_sequencing/data/output_data/deg_tests/%s_glmQLFTest.csv", col)
  write.table(out, sep = ",", file=output,col.names=NA)
}


