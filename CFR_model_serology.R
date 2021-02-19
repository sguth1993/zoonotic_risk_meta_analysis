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

## Model selection
Terms <- list(
  hOrder = c(NA, "s(hOrder, bs='re')"),
  phylo_dist = c(NA, "s(phylo_dist, bs='re')"),
  VirusFamilyPublicationCount  = c(NA, "s(VirusFamilyPublicationCount, k=7, bs='tp')"),
  vFamily  = c(NA, "s(vFamily, bs='re')"),
  ReservoirPublicationCount  = c(NA, "s(ReservoirPublicationCount, k=7, bs='tp')"),
  IsVectorBorne  = c(NA, "s(IsVectorBorne, bs='re')"),
  spill_type  = c(NA, "s(spill_type, bs = 're')")
)

## All possible combinations of these terms:
CompetingFullModels <- expand.grid(Terms)

## Build formulas
CompetingFullModels <- apply(CompetingFullModels, 1, function(row) paste(na.omit(row), collapse = ' + ')) %>% 
  data.frame(Formula = ., stringsAsFactors = F) %>% 
  mutate(Formula = ifelse(nchar(Formula) == 0, '1', Formula),
         Formula = paste('CFR_avg ~', Formula))

## Model fit
CompetingFits <- CompetingFullModels %>% 
  group_by(Formula) %>% 
  do(ModelFit = try(gamm(as.formula(.$Formula), 
                         family = 'gaussian',
                         data = dat.mort,
                         method = 'REML')))

removeRows <- lapply(CompetingFits$ModelFit, function(x) 'try-error' %in% class(x)) %>% 
  unlist()
FailedFormulas <- CompetingFits[removeRows, ]
stopifnot(nrow(FailedFormulas) == 0)

CompetingFits <- CompetingFits[!removeRows, ]

## Add AIC:
RankedModels <- CompetingFits %>% 
  mutate(AIC = AIC(ModelFit),
         DevianceExplained = summary(ModelFit$gam)$d) %>% 
  ungroup() %>% 
  arrange(AIC) %>% 
  mutate(DeltaAIC = AIC - AIC[1])

## run model with lowest AIC to make plots
lowest_AIC <- RankedModels$Formula[which.min(RankedModels$AIC)]

gamm_mort_lowest <- gamm(as.formula(lowest_AIC),
                         data = dat.mort,
                         family = gaussian)
summary(gamm_mort_lowest$gam)
gam.check(gamm_mort_lowest$gam)

#is there a model with more predictors and a low deltaAIC?
models <- RankedModels[str_detect(RankedModels$Formula, "IsVectorBorne")==TRUE & str_detect(RankedModels$Formula, "spill_type")==TRUE, ]
best_fit <- models[which.min(models$AIC), ]
best_fit$DeltaAIC #deltaAIC is only 0.83
best_fit$Formula

#### SELECTED MODEL ####
gamm_mort <- gamm(as.formula(best_fit$Formula),
                  data = dat.mort,
                  family = gaussian)
summary(gamm_mort$gam)
gam.check(gamm_mort$gam)


#### HOST ORDER EFFECTS ####
order_dat <- plot(gamm_mort$gam, residuals = TRUE, pages =1)[[1]]
orders <- as.factor(levels(order_dat$raw))
order_effects <- cbind.data.frame(orders, order_dat$fit)
names(order_effects) <- c("hOrder", "effect")
print(order_effects) 

#### VIRUS FAMILY EFFECTS ####
vfam_dat <- plot(gamm_mort$gam, residuals = TRUE, pages =1)[[3]]
vfams <- as.factor(levels(vfam_dat$raw))

vfam_effects <- cbind.data.frame(vfams, vfam_dat$fit)
names(vfam_effects) <- c("vFamily", "effect")
print(vfam_effects) 

#### VECTOR ####
vector_dat <- plot(gamm_mort$gam, residuals = TRUE, pages =1)[[5]]
vectors <- as.factor(levels(vector_dat$raw))

vector_effects <- cbind.data.frame(vectors, vector_dat$fit)
names(vector_effects) <- c("vector", "effect")
print(vector_effects) 

#### SPILL TYPE ####
spill_dat <- plot(gamm_mort$gam, residuals = TRUE, pages =1)[[6]]
spills <- as.factor(levels(spill_dat$raw))

spill_effects <- cbind.data.frame(spills, spill_dat$fit)
names(spill_effects) <- c("spill", "effect")
print(spill_effects) 

