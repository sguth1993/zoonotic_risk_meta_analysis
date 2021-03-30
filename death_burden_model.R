library(dplyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(mgcv)
library(cowplot)
library(MuMIn)
library(corrplot)
library(lme4)
library("ape")
library(VGAM)
library(scatterpie)
library(glmmTMB)
setwd("~/Documents/Berkeley/Virulence_traits")
rm(list=ls())

## load data
dat <- read.csv(file="Extension/code_and_data/data/stringent_data_new_3-27.csv", header = T, stringsAsFactors = F)

## convert predictors to the correct class
dat$SppName_ICTV_MSL2018b <- as.factor(dat$SppName_ICTV_MSL2018b)
dat$hOrder <- as.factor(dat$hOrder)
dat$phylo_dist <- as.factor(dat$phylo_dist)
dat$spill_type <- as.factor(dat$spill_type)
dat$vFamily <- as.factor(dat$vFamily)
dat$IsVectorBorne <- as.factor(dat$IsVectorBorne)
dat$death_burden_since_1950 <- as.numeric(dat$death_burden_since_1950)
dat$vaccination_effects <- as.factor(dat$vaccination_effects)
dat$Tr.primary <- as.factor(dat$Tr.primary)

#calculate years since death count start
dat$death_period <- dat$death_count_end - dat$death_count_start

# note that unlike the CFR models, we're keeping viruses with only one recorded case
dat_sub <- dat

#remove Rotavirus A
dat_sub$SppName_ICTV_MSL2018b[is.na(dat_sub$death_burden_since_1950)==TRUE]
dat_sub <- dat_sub[is.na(dat_sub$death_burden_since_1950)==FALSE, ]

#final set of model variables
dat.death <- dat_sub %>% distinct(SppName_ICTV_MSL2018b, death_burden_since_1950, 
                                  hOrder, 
                                  #spill_type, 
                                  .keep_all = TRUE)
unique(dat.death$SppName_ICTV_MSL2018b)

#check sample size for reservoir orders and virus families
order_sum <- dat.death%>%
  group_by(hOrder) %>%
  summarise(count = n())
#should drop Eulipotyphla, Pilosa, and Perissodactyla

vir_sum <- dat.death%>%
  group_by(vFamily) %>%
  summarise(count = n())
#should drop Peribunyaviridae,  Hepeviridae 

#  Create dummy variables for orders and virus families to use as random effects
dummy_orders = as.data.frame(model.matrix(~hOrder, 
                                          data=dat.death, 
                                          contrasts.arg=list(hOrder=contrasts(dat.death$hOrder, contrasts=F))))[,-1]
colnames(dummy_orders)
dummy_orders <- dummy_orders[ ,-c(6:8)]

dummy_virfams = as.data.frame(model.matrix(~vFamily, 
                                           data=dat.death, 
                                           contrasts.arg=list(vFamily=contrasts(dat.death$vFamily, contrasts=F))))[,-1]
colnames(dummy_virfams)
dummy_virfams <- dummy_virfams[ ,-c(7,11)]

dummy_spill = as.data.frame(model.matrix(~spill_type, 
                                         data=dat.death, 
                                         contrasts.arg=list(spill_type=contrasts(dat.death$spill_type, contrasts=F))))[,-1]


dummy_vector = as.data.frame(model.matrix(~IsVectorBorne, 
                                          data=dat.death, 
                                          contrasts.arg=list(IsVectorBorne=contrasts(dat.death$IsVectorBorne, contrasts=F))))[,-1]
dummy_vaccine = as.data.frame(model.matrix(~vaccination_effects, 
                                           data=dat.death, 
                                           contrasts.arg=list(vaccination_effects=contrasts(dat.death$vaccination_effects, contrasts=F))))[,-1]
dummy_trans = as.data.frame(model.matrix(~Tr.primary,
                                                 data=dat.death,
                                                 contrasts.arg=list(Tr.primary=contrasts(dat.death$Tr.primary, contrasts=F))))[,-1]
dummys = cbind(dummy_orders, dummy_virfams, dummy_vector, dummy_vaccine, 
               dummy_trans,
               dummy_spill)
dummys <- lapply(dummys, as.factor)

data_set = cbind(dat.death, dummys)

#print for pasting into terms and models
paste0("s(", names(dummys), ", bs = 're')")

#global model with automated term selection
gam_death <- gam(death_burden_since_1950 ~  
                   s(vFamilyArenaviridae , bs = 're') +
                   s(vFamilyBornaviridae , bs = 're') +
                   s(vFamilyBunyaviridae , bs = 're') +
                   s(vFamilyCoronaviridae , bs = 're') +
                   s(vFamilyFiloviridae , bs = 're') +
                   s(vFamilyOrthomyxoviridae , bs = 're') +
                   s(vFamilyFlaviviridae , bs = 're') +
                   s(vFamilyHerpesviridae , bs = 're') +
                   s(vFamilyParamyxoviridae , bs = 're') +
                   s(vFamilyPicornaviridae , bs = 're') +
                   s(vFamilyPoxviridae , bs = 're') +
                   s(vFamilyReoviridae , bs = 're') +
                   s(vFamilyRetroviridae , bs = 're') +
                   s(vFamilyRhabdoviridae , bs = 're') +
                   s(vFamilyTogaviridae , bs = 're') +
                   
                   s(hOrderCARNIVORA , bs = 're') +
                   s(hOrderCETARTIODACTYLA , bs = 're') +
                   s(hOrderCHIROPTERA , bs = 're') +
                   s(hOrderDIPROTODONTIA , bs = 're') +
                   s(hOrderPRIMATES , bs = 're') +
                   s(hOrderRODENTIA , bs = 're') +
                   
                   s(Tr.primarybiting, bs = 're') +                         
                   s(Tr.primarybodily_fluids, bs = 're')   +              
                   s(Tr.primarydirect_contact, bs = 're')   +              
                   s(Tr.primaryfecal_oral, bs = 're')        +            
                   s(Tr.primaryinhalation_aerosolized_excreta, bs = 're') +
                   s(Tr.primaryrespiratory, bs = 're')                   +
                   s(Tr.primaryvector, bs = 're') +
     
                   s(spill_typebridged, bs = 're') +
                   
                   s(vaccination_effectsY, bs = 're') +
                  
                   s(VirusSppPublicationCount, k =7, bs ='tp'),
                   #s(VirusFamilyPublicationCount, k =7, bs ='tp'),
                 offset = scale(death_period),
                 data = data_set,
                 select = TRUE,
                 family = nb())
summary(gam_death)
AIC(gam_death)
gam.check(gam_death)
plot(gam_death)

#selected model
selected_gam_death <- gam(death_burden_since_1950 ~  
                            #s(vFamilyArenaviridae , bs = 're') +
                            #s(vFamilyBornaviridae , bs = 're') +
                            s(vFamilyBunyaviridae , bs = 're') +
                            s(vFamilyCoronaviridae , bs = 're') +
                            s(vFamilyFiloviridae , bs = 're') +
                            #s(vFamilyOrthomyxoviridae , bs = 're') +
                            s(vFamilyFlaviviridae , bs = 're') +
                            #s(vFamilyHerpesviridae , bs = 're') +
                            s(vFamilyParamyxoviridae , bs = 're') +
                            s(vFamilyPicornaviridae , bs = 're') +
                            s(vFamilyPoxviridae , bs = 're') +
                            s(vFamilyReoviridae , bs = 're') +
                            s(vFamilyRetroviridae , bs = 're') +
                            #s(vFamilyRhabdoviridae , bs = 're') +
                            #s(vFamilyTogaviridae , bs = 're') +
                            
                            #s(hOrderCARNIVORA , bs = 're') +
                            #s(hOrderCETARTIODACTYLA , bs = 're') +
                            #s(hOrderCHIROPTERA , bs = 're') +
                            s(hOrderDIPROTODONTIA , bs = 're') +
                            s(hOrderPRIMATES , bs = 're') +
                            #s(hOrderRODENTIA , bs = 're') +
                            
                            #s(Tr.primarybiting, bs = 're') +                         
                            #s(Tr.primarybodily_fluids, bs = 're')   +              
                            s(Tr.primarydirect_contact, bs = 're')   +              
                            s(Tr.primaryfecal_oral, bs = 're')        +            
                            # s(Tr.primaryinhalation_aerosolized_excreta, bs = 're') +
                            #s(Tr.primaryrespiratory, bs = 're')                   +
                            s(Tr.primaryvector, bs = 're') +
                            
                            #s(spill_typebridged , bs = 're'), # +
                          
                          s(vaccination_effectsY , bs = 're') + 
                          #s(VirusFamilyPublicationCount, k =7, bs ='tp'),
                          s(VirusSppPublicationCount, k =7, bs ='tp'),
                            
                          offset = scale(death_period),
                          data = data_set,
                          select = TRUE,
                          family = nb())
summary(selected_gam_death)
coef(selected_gam_death)
gam.check(selected_gam_death)

#### FIGURES ####
plotData <- get_partial_effects_binary(fit = gam_death, 
                                       var = c('vFamilyBunyaviridae',
                                               'vFamilyCoronaviridae',
                                               'vFamilyFiloviridae',
                                               'vFamilyFlaviviridae',
                                               'vFamilyParamyxoviridae',
                                               'vFamilyPicornaviridae',
                                               'vFamilyPoxviridae',
                                               'vFamilyReoviridae',
                                               'vFamilyRetroviridae',
                                               'hOrderDIPROTODONTIA',
                                               'hOrderPRIMATES',
                                        
                                               'Tr.primarydirect_contact',
                                               'Tr.primaryfecal_oral',
                                               'Tr.primaryvector',
                                               
                                               'vaccination_effectsY'),
                                       seWithMean = TRUE, fixedEffect = FALSE, removeNegatives = TRUE) 


ggplot(plotData$effects, aes(x = variable, y = y, colour = IsSignificant, fill = IsSignificant)) +
  geom_hline(yintercept = 0, linetype = 2, colour = 'grey20') +
  geom_boxplot(aes(middle = y, lower = ylower, upper = yupper, ymin = ylower, ymax = yupper), 
               stat = 'identity', alpha = 0.5, colour = NA) +
  geom_boxplot(aes(middle = y, lower = y, upper = y, ymin = y, ymax = y), 
               stat = 'identity', alpha = 0.5, size = 0.5) +
  
  geom_jitter(aes(y = Residual), alpha = 0.6, size = 0.8, width = 0.35, height = 0,
              data = plotData$partialResiduals) +
  scale_colour_manual(values = c(No = 'grey60', Yes = '#225ea8'), guide = F) +
  scale_fill_manual(values = c(No = 'grey60', Yes = '#225ea8'), guide = F) +
  scale_y_continuous(labels = function(x) sprintf('%1.1f', x) ) +
  labs(x = NULL, y = 'Effect on death burden') +
  #scale_x_discrete(labels = gam_mort_labs) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())


#pdf file
ggsave(file = "Extension/code_and_data/figures/Figure_2.pdf",
       units="mm",
       width=80,
       height=80,
       scale=2,
       dpi=300)



#more diagnostics
# Model fit
ggplot(NULL, aes(x = log1p(selected_gam_death$y), y = log1p(selected_gam_death$fitted.values))) +
  geom_jitter(shape = 1, width = 0.05, height = 0.05) +
  geom_abline(linetype = 2) + 
  labs(x = "log(1 + observed)", y = "log(1 + predicted)")


