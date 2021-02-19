library(dplyr)
library(stringr)
library(ggplot2)
library(mgcv)
library(cowplot)
library(MuMIn)
library(corrplot)
library(lme4)
library("ape")
library(VGAM)
library(glmmTMB)
library(DHARMa)
setwd("~/Documents/Berkeley/Virulence_traits")
rm(list=ls())

## load data
dat <- read.csv(file="Extension/code_and_data/data/stringent_data.csv", header = T, stringsAsFactors = F)

## convert categorical predictors to factors
dat$SppName_ICTV_MSL2018b <- as.factor(dat$SppName_ICTV_MSL2018b)
dat$hOrder <- as.factor(dat$hOrder)
dat$spill_type <- as.factor(dat$spill_type)
dat$phylo_dist <- as.factor(dat$phylo_dist) #also modeled as categorical random effect in PNAS (only 5 unique values)
dat$vFamily <- as.factor(dat$vFamily)
dat$human.trans <- as.integer(dat$human.trans)
dat$IsVectorBorne <- as.factor(dat$IsVectorBorne)

# note that unlike the CFR models, we're keeping viruses with only one recorded case
dat_sub <- dat
## round cases and deaths (for poisson model) that we had to proxy by dividing by number of years
dat_sub$deaths_per_outbreak_rounded <- round(dat_sub$deaths_per_outbreak)
dat_sub$cases_per_outbreak_rounded <- round(dat_sub$cases_per_outbreak)

## bin outbreak year by decade
# and calculate predictor variable: years since outbreak
dat_sub$outbreak_year <- as.numeric(dat_sub$outbreak_year)
dat_sub$outbreak_decade <-floor(dat_sub$outbreak_year/10)*10
dat_sub$years_ago <- 2020 - dat_sub$outbreak_year

## global model
global_fit <- glmmTMB(cases_per_outbreak_rounded ~ 
                  #log(GDP) +
                  #log(pop_size) +
                  #phylo_dist + 
                  #years_ago + 
                  hOrder + 
                  vFamily + 
                  human.trans + 
                  IsVectorBorne +
                  #log(ReservoirPublicationCount) +
                  log(VirusFamilyPublicationCount) +
                  (1 | spill_type),
                data = dat_sub,
                family = "truncated_nbinom2")
summary(global_fit)


  

## Model selection
Terms <- list(
  GDP = c(NA, "log(GDP)"),
  pop_size = c(NA, "log(pop_size)"),
  phylo_dist = c(NA, "phylo_dist"),
  years_ago = c(NA, "years_ago"),
  hOrder = c(NA, "hOrder"),
  vFamily = c(NA, "vFamily"),
  human.trans = c(NA, "human.trans"),
  IsVectorBorne = c(NA, "IsVectorBorne"),
  ReservoirPublicationCount = c(NA, "log(ReservoirPublicationCount)"),
  VirusFamilyPublicationCount= c(NA, "log(VirusFamilyPublicationCount)"),
  spill_type = c(NA, "(1 | spill_type)"))

## All possible combinations of these terms:
CompetingFullModels <- expand.grid(Terms)

## Build formulas
CompetingFullModels <- apply(CompetingFullModels, 1, function(row) paste(na.omit(row), collapse = ' + ')) %>% 
  data.frame(Formula = ., stringsAsFactors = F) %>% 
  mutate(Formula = ifelse(nchar(Formula) == 0, '1', Formula),
         Formula = paste('cases_per_outbreak_rounded ~', Formula))

## Model fit
CompetingFits <- CompetingFullModels %>% 
  group_by(Formula) %>% 
  do(ModelFit = try(glmmTMB(as.formula(.$Formula), 
                            data = dat_sub,
                         family = 'truncated_nbinom2')))

removeRows <- lapply(CompetingFits$ModelFit, function(x) 'try-error' %in% class(x)) %>% 
  unlist()
FailedFormulas <- CompetingFits[removeRows, ]
stopifnot(nrow(FailedFormulas) == 0)

CompetingFits <- CompetingFits[!removeRows, ]

## Add AIC:
RankedModels <- CompetingFits %>% 
  mutate(AIC = AIC(ModelFit)) %>% 
        
  ungroup() %>% 
  arrange(AIC) %>% 
  mutate(DeltaAIC = AIC - AIC[1])

## run model with lowest AIC to make plots
lowest_AIC <- RankedModels$Formula[which.min(RankedModels$AIC)]
best_fit <- RankedModels$Formula[3]
RankedModels$DeltaAIC[3]
#"cases_per_outbreak_rounded ~ log(GDP) + log(pop_size) + years_ago + hOrder + vFamily"

fit1 <- glmmTMB(as.formula(best_fit),
        data = dat_sub,
        family = "truncated_nbinom2")
summary(fit1)

effects <- allEffects(fit1)
effects
## Diagnostics:

# Model fit
ggplot(NULL, aes(x = log1p(fit1$frame$cases_per_outbreak_rounded), y = log1p(predict(fit1, type = "response")))) +
  geom_jitter(shape = 1, width = 0.05, height = 0.05) +
  geom_abline(linetype = 2) #+ 
  labs(x = "log(1 + observed)", y = "log(1 + predicted)")

# Simulated residuals
simulation_output <- simulateResiduals(fittedModel = fit1)

plot(simulation_output)
testUniformity(simulation_output)
testZeroInflation(simulation_output)
testDispersion(simulation_output)

plotResiduals(simulation_output, log(dat_sub$GDP))
plotResiduals(simulation_output, log(dat_sub$ReservoirPublicationCount))
plotResiduals(simulation_output, log(dat_sub$VirusFamilyPublicationCount))
plotResiduals(simulation_output, dat_sub$IsVectorBorne)
plotResiduals(simulation_output, dat_sub$outbreak_start)
plotResiduals(simulation_output, dat_sub$outbreak_country)
plotResiduals(simulation_output, dat_sub$spill_type)
plotResiduals(simulation_output, dat_sub$hOrder)
plotResiduals(simulation_output, dat_sub$vFamily)
plotResiduals(simulation_output, dat_sub$SppName_ICTV_MSL2018b)


