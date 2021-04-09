library(dplyr)
library(stringr)
library(ggplot2)
library(mgcv)
library(cowplot)
library(MuMIn)
library(corrplot)
library(lme4)
library(scatterpie)
setwd("~/Documents/Berkeley/Virulence_traits")
rm(list=ls())

## load data
dat <- read.csv(file="Extension/code_and_data/data/stringent_data.csv", header = T, stringsAsFactors = F)
dat_sub <- dat[ ,c("hOrder","SppName_ICTV_MSL2018b", "CFR_avg",
                   "phylo_dist",
                   "VirusFamilyPublicationCount",
                   "vFamily",
                   "ReservoirPublicationCount",
                   "IsVectorBorne",
                   "spill_type")]

sero <- read.csv(file="Extension/code_and_data/data/loose_data.csv", header = T, stringsAsFactors = F)
unique(sero$SppName_ICTV_MSL2018b) #24 viruses, each with one transmission chain
sero_sub <- sero[ ,c("hOrder","SppName_ICTV_MSL2018b", "CFR_avg",
                    "phylo_dist",
                    "VirusFamilyPublicationCount",
                    "vFamily",
                    "ReservoirPublicationCount",
                    "IsVectorBorne",
                    "spill_type")]

dat <- rbind(sero_sub, dat_sub)

## convert all categorical predictors to factors
dat$SppName_ICTV_MSL2018b <- as.factor(dat$SppName_ICTV_MSL2018b)
dat$hOrder <- as.factor(dat$hOrder)
dat$spill_type <- as.factor(dat$spill_type)
dat$phylo_dist <- as.factor(dat$phylo_dist) #also modeled as categorical random effect in PNAS (only 5 unique values)
dat$vFamily <- as.factor(dat$vFamily)
dat$IsVectorBorne <- as.factor(dat$IsVectorBorne)
str(dat)
dat[complete.cases(dat)==FALSE,  ]

#final set of model variables
dat.mort <- dat %>% distinct(SppName_ICTV_MSL2018b, CFR_avg, hOrder, spill_type, .keep_all = TRUE)

#how many viruses have multiple entries?
dup <- dat.mort[duplicated(dat.mort$SppName_ICTV_MSL2018b)==TRUE, ] #only 9 viruses
dup_dat <- dat.mort %>%
  filter(SppName_ICTV_MSL2018b %in% dup$SppName_ICTV_MSL2018b)


#global model w/ automated term selection: global CFR estimates
gam_mort <- gam(CFR_avg ~
                  s(vFamily, bs = 're'  ) +
                  s(hOrder, bs = 're'  ) +
                  s(VirusFamilyPublicationCount, k=7, bs="tp") +
                  s(spill_type, bs = 're' ) +
                  s(IsVectorBorne, bs = 're'),
                data = dat.mort,
                select=TRUE,
                family = gaussian)
summary(gam_mort) #no variables dropped
gam.check(gam_mort)

#selected model
#only one variable drops
gam_mort <- gam(CFR_avg ~
                  s(vFamily, bs = 're'  ) +
                  s(hOrder, bs = 're'  ) +
                 # s(VirusFamilyPublicationCount, k=7, bs="tp") +
                  s(spill_type, bs = 're' ) +
                  s(IsVectorBorne, bs = 're'),
                data = dat.mort,
                select=TRUE,
                family = gaussian)
summary(gam_mort) #no variables dropped
gam.check(gam_mort)

effects_dat <- plot.gam(gam_mort, residuals=TRUE, select = 2)[[2]]
effects <- as.factor(levels(effects_dat$raw))
effects <- cbind.data.frame(effects, effects_dat$fit)
print(effects)

effects_dat <- plot.gam(gam_mort, residuals=TRUE, select = 2)[[3]]
effects <- as.factor(levels(effects_dat$raw))
effects <- cbind.data.frame(effects, effects_dat$fit)
print(effects)

effects_dat <- plot.gam(gam_mort, residuals=TRUE, select = 2)[[4]]
effects <- as.factor(levels(effects_dat$raw))
effects <- cbind.data.frame(effects, effects_dat$fit)
print(effects)




