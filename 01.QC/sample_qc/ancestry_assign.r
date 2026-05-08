#!/usr/bin/env Rscript

pcadir <- "/stanley/huang_lab/home/myu/IBD/preprocess/QC/ancestry_assign"
pcplot0 <- TRUE      #TRUE or FALSE
npc <- 10
pred_prob <- 0.8

###
setwd(pcadir)
library(randomForest)
library(tidyverse)
library(ggsci)
library(RColorBrewer)


print("Read in data")
# Read in data
pca_w_ref <- read.table("all_samples.pca.txt", header=T, sep='\t', stringsAsFactors=F)
pca_w_ref$pop <- pca_w_ref$FID


pca_w_ref$superpop <- "gnomad"
# if ( refname == "1kg" ) {
pca_w_ref$superpop <- ifelse(pca_w_ref$pop=="CHB" | pca_w_ref$pop=="JPT" | pca_w_ref$pop=="CHS" | pca_w_ref$pop=="CDX" | pca_w_ref$pop=="KHV", "EAS", pca_w_ref$superpop)
pca_w_ref$superpop <- ifelse(pca_w_ref$pop=="CEU" | pca_w_ref$pop=="TSI" | pca_w_ref$pop=="FIN" | pca_w_ref$pop=="GBR" | pca_w_ref$pop=="IBS", "EUR", pca_w_ref$superpop)
pca_w_ref$superpop <- ifelse(pca_w_ref$pop=="YRI" | pca_w_ref$pop=="LWK" | pca_w_ref$pop=="GWD" | pca_w_ref$pop=="MSL" | pca_w_ref$pop=="ESN" | pca_w_ref$pop=="ASW" | pca_w_ref$pop=="ACB", "AFR.AMR", pca_w_ref$superpop)
pca_w_ref$superpop <- ifelse(pca_w_ref$pop=="MXL" | pca_w_ref$pop=="PUR" | pca_w_ref$pop=="CLM" | pca_w_ref$pop=="PEL", "AFR.AMR", pca_w_ref$superpop)
pca_w_ref$superpop <- ifelse(pca_w_ref$pop=="GIH" | pca_w_ref$pop=="PJL" | pca_w_ref$pop=="BEB" | pca_w_ref$pop=="STU" | pca_w_ref$pop=="ITU", "SAS", pca_w_ref$superpop)

pca_w_ref$superpop <- factor(pca_w_ref$superpop,
                              levels=c("gnomad","AFR.AMR","EAS","EUR","SAS"))

# }


print("Make PC plots: gnomad + ref samples, colored by super population")
# PC plots: gnomad + ref samples, colored by super population
if ( pcplot0 == "TRUE" ){
    for(i in 1:2){
        p = ggplot(pca_w_ref, aes_string(x=paste0('PC',i), y=paste0('PC',i+1))) +
        geom_point(aes(color=superpop), alpha=0.57, size=0.97) +
        scale_color_manual(values = c("grey25",pal_d3("category20")(20)[c(1:5)])) + 
        theme_bw() +
        guides(color = guide_legend(override.aes = list(size=2))) +
        labs(x=paste0('PC',i), y=paste0('PC',i+1), title="gnomad + 1KG samples",
        color="Super population")
    # print(p)
        ggsave(paste0("gnomad_ref.PC",i,"_PC",i+1,"_by_superpop.png"), p, width=7.5, height=6)
    }
}

print(paste0("Use RF to predict ancestry based on top ",npc," PCs"))
# RF function to predict ancestry using PCs:
pop_forest <- function(training_data, data, ntree=100, seed=42, pcs=1:npc) {
  set.seed(seed)
  form <- formula(paste('as.factor(known_pop) ~', paste0('PC', pcs, collapse = ' + ')))
  forest <- randomForest(form,
                        data = training_data,
                        importance = T,
                        ntree = ntree)
  print(forest)
  fit_data <- data.frame(predict(forest, data, type='prob'), sample = data$sample)
  fit_data %>%
    gather(predicted_pop, probability, -sample) %>%
    group_by(sample) %>%
    slice(which.max(probability))
}

# Prep data
# Explanation:
# `trdat` and `tedat` are training and testing data. training data has to have a 
# column `known_pop` and `PC1` to `PC10` or so. `tedat` is expected to have a 
# column `sample` which is just the sample ID, and also PC columns.


trdat <- pca_w_ref %>%
  filter(superpop != 'gnomad') %>%
  select(superpop, PC1:PC10) %>%
  rename(known_pop = superpop)
trdat$known_pop <- as.character(trdat$known_pop)

tedat <- pca_w_ref %>%
  filter(superpop=='gnomad') %>%
  select(IID, PC1:PC10)
names(tedat)[1] <- "sample"
tedat$sample <- as.character(tedat$sample)



print("Prediction results:")
# Make prediction
pop_pred <- as.data.frame(pop_forest(training_data = trdat, data = tedat))
# when number of PCs (npc = 6)
# Type of random forest: classification
#                      Number of trees: 100
# No. of variables tried at each split: 2

#         OOB estimate of  error rate: 0.16%
# Confusion matrix:
#     AFR AMR EAS EUR SAS class.error
# AFR 659   2   0   0   0 0.003025719
# AMR   2 345   0   0   0 0.005763689
# EAS   0   0 504   0   0 0.000000000
# EUR   0   0   0 503   0 0.000000000
# SAS   0   0   0   0 489 0.000000000


print("Overall ancestry assignment:")
table(pop_pred$predicted_pop)
 #  AFR   AMR   EAS   EUR   SAS
 # 1861  5396   582 28266   319

# summary(pop_pred$probability)
# dim(pop_pred %>% filter(probability<0.5))
# dim(pop_pred %>% filter(probability<0.9))

print(paste0("Subset to prediction prob. > ", pred_prob, ":"))
pop_pred.sub <- pop_pred %>%
  filter(probability > pred_prob)
table(pop_pred.sub$predicted_pop)
# when pred_prob = 0.9
 #  AFR   AMR   EAS   EUR   SAS
 # 1607  1840   504 26677   297

# when pred_prob = 0.8
 #  AFR   AMR   EAS   EUR   SAS
 # 1657  2056   511 27090   299



# Merge back with pca_w_ref to udpate PC plots (colored by pred. super population)
names(pop_pred)[1] <- "IID"
pca_w_ref.gnomad <- subset(pca_w_ref, superpop=="gnomad")
pca_w_ref.pred <- merge(pca_w_ref.gnomad, pop_pred, by="IID")
pca_w_ref.pred$predicted_pop <- factor(pca_w_ref.pred$predicted_pop, levels=c("AFR.AMR","EAS","EUR","SAS"))


print("Make PC plots: gnomad samples only, colored by predicted super population/ancestry")
# PC plots: gnomad samples only, colored by predicted super population/ancestry
for(i in 1:2){
  p = ggplot(pca_w_ref.pred, aes_string(x=paste0('PC',i), y=paste0('PC',i+1))) +
    geom_point(aes(color=predicted_pop), alpha=0.57, size=0.97) +
    scale_color_manual(values = pal_d3("category20")(20)[1:5]) +
    theme_bw() +
    guides(color = guide_legend(override.aes = list(size=2))) +
    labs(x=paste0('PC',i), y=paste0('PC',i+1), title="gnomad samples (in 1KG space)",
         color="Predicted population")
  # print(p)
  ggsave(paste0("gnomad_ref.PC",i,"_PC",i+1,"_by_predPop.png"),
         p, width=7.75, height=6)
}


print(paste0("Make PC plots: gnomad samples only, colored by predicted superpop w. a pred. prob > ",pred_prob))
# PC plots: gnomad samples only, colored by predicted super population w. a pred. prob > pred_prob
pca_w_ref.pred.sub <- subset(pca_w_ref.pred, probability > pred_prob)
for(i in 1:2){
  p = ggplot(pca_w_ref.pred.sub, aes_string(x=paste0('PC',i), y=paste0('PC',i+1))) +
    geom_point(aes(color=predicted_pop), alpha=0.57, size=0.97) +
    scale_color_manual(values = pal_d3("category20")(20)[1:5]) +
    theme_bw() +
    guides(color = guide_legend(override.aes = list(size=2))) +
    labs(x=paste0('PC',i), y=paste0('PC',i+1), title="gnomad samples (in 1KG space)",
         color="Predicted population")
  # print(p)
  ggsave(paste0("gnomad_ref.PC",i,"_PC",i+1,"_by_predPop",pred_prob,".png"),
         p, width=7.75, height=6)
}


print("Export predicted ancestral groups for study samples:")
# Save predicted population labels and EUR family and individual IDs
write.table(pca_w_ref.pred, "gnomad_ref.PC.predPop.tsv", quote=F,
  col.names=T, row.names=F, sep="\t")

# eur_inds <- pca_w_ref.pred.sub[pca_w_ref.pred.sub$predicted_pop=="EUR", c("FID","IID")]
# write.table(eur_inds, paste0("gnomad_ref.PC.predPop",pred_prob,".EUR.indlist"), quote=F,
#   col.names=F, row.names=F, sep="\t")

