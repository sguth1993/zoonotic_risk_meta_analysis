library(dplyr)
library(stringr)
library(ggplot2)
library(mgcv)
library(cowplot)
library(MuMIn)
library(corrplot)
library(lme4)
library(scatterpie)
library(reshape2)
setwd("~/Documents/Berkeley/Virulence_traits")
rm(list=ls())

## load data
dat <- read.csv(file="Extension/code_and_data/data/stringent_data.csv", header = T, stringsAsFactors = F)

## convert all categorical predictors to factors
dat$SppName_ICTV_MSL2018b <- as.factor(dat$SppName_ICTV_MSL2018b)
dat$hOrder <- as.factor(dat$hOrder)
dat$spill_type <- as.factor(dat$spill_type)
dat$vFamily <- as.factor(dat$vFamily)
dat$phylo_dist <- as.factor(dat$phylo_dist)
dat$IsVectorBorne <- as.factor(dat$IsVectorBorne)

# remove viruses with only one human case recorded in history
dat_rare <- dat[dat$tot_cases==1, ]
dat_suf <- dat[dat$tot_cases>1, ]

#final set of model variables
dat.mort <- dat_suf %>% distinct(SppName_ICTV_MSL2018b, CFR_avg, hOrder, spill_type, .keep_all = TRUE)
length(unique(dat.mort$SppName_ICTV_MSL2018b))

#how many viruses have multiple entries?
dup <- dat.mort[duplicated(dat.mort$SppName_ICTV_MSL2018b)==TRUE, ] #only 9 viruses
dup_dat <- dat.mort %>%
  filter(SppName_ICTV_MSL2018b %in% dup$SppName_ICTV_MSL2018b)

#check sample size for reservoir orders and virus families
order_sum <- dat.mort%>%
  group_by(hOrder) %>%
  summarise(count = n())
#should drop Eulipotyphla and Perissodactyla

vir_sum <- dat.mort%>%
  group_by(vFamily) %>%
  summarise(count = n())
#should drop Peribunyaviridae,  Hepeviridae 

#  Create dummy variables for orders and virus families to use as random effects
dummy_orders = as.data.frame(model.matrix(~hOrder, 
                                    data=dat.mort, 
                                    contrasts.arg=list(hOrder=contrasts(dat.mort$hOrder, contrasts=F))))[,-1]
colnames(dummy_orders)
dummy_orders <- dummy_orders[ ,-c(6:7)]

dummy_virfams = as.data.frame(model.matrix(~vFamily, 
                                          data=dat.mort, 
                                          contrasts.arg=list(vFamily=contrasts(dat.mort$vFamily, contrasts=F))))[,-1]
colnames(dummy_virfams)
dummy_virfams <- dummy_virfams[ ,-c(7,11)]

dummy_spill = as.data.frame(model.matrix(~spill_type, 
                                         data=dat.mort, 
                                         contrasts.arg=list(spill_type=contrasts(dat.mort$spill_type, contrasts=F))))[,-1]
dummy_vector = as.data.frame(model.matrix(~IsVectorBorne, 
                                         data=dat.mort, 
                                         contrasts.arg=list(IsVectorBorne=contrasts(dat.mort$IsVectorBorne, contrasts=F))))[,-1]
dummys = cbind(dummy_orders, dummy_virfams, dummy_spill, dummy_vector)
data_set = cbind(dat.mort, dummys)

#print for pasting into terms and models
paste0("s(", names(dummys), ", bs = 're')")

#global model: use automated term selection by double penalty smoothing to identify variables to remove
#let's first see if GDPPC is significant 
country_gam_mort <- gam(CFR_country_specific ~ 
                         s(weighted_GDP, k=7, bs='tp') +
                         s(phylo_dist, bs = 're') +
                         s(ReservoirPublicationCount, k=7, bs="tp") +
                         s(VirusFamilyPublicationCount, k=7, bs="tp") +
                         
                         #spillover type random effects
                         s(spill_typebridged, bs = 're') +
                         #s(spill_typedirect, bs = 're') +       
                         
                         #vector-borne random effects
                         #s(IsVectorBorne0, bs = 're') +
                         s(IsVectorBorne1, bs = 're') +
                         
                         #host order random effects
                         s(hOrderAVES, bs = 're') +
                         s(hOrderCARNIVORA, bs = 're') +
                         s(hOrderCETARTIODACTYLA, bs = 're') +
                         s(hOrderCHIROPTERA, bs = 're') +
                         s(hOrderDIPROTODONTIA, bs = 're') +
                        # s(hOrderEULIPOTYPHLA, bs = 're') +
                        # s(hOrderPERISSODACTYLA, bs = 're') +
                         s(hOrderPRIMATES, bs = 're') +
                         s(hOrderRODENTIA, bs = 're') +
                         
                         #virus family random effects
                         s(vFamilyArenaviridae, bs = 're') +
                         s(vFamilyBornaviridae, bs = 're') +
                         s(vFamilyBunyaviridae, bs = 're') +
                         s(vFamilyCoronaviridae, bs = 're') +
                         s(vFamilyFiloviridae, bs = 're') +
                         s(vFamilyFlaviviridae, bs = 're') +
                        # s(vFamilyHepeviridae, bs = 're') +
                         s(vFamilyHerpesviridae, bs = 're') +
                         s(vFamilyOrthomyxoviridae, bs = 're') +
                         s(vFamilyParamyxoviridae, bs = 're') +
                       #  s(vFamilyPeribunyaviridae, bs = 're') +
                         s(vFamilyPicornaviridae, bs = 're') +
                         s(vFamilyPoxviridae, bs = 're') +
                         s(vFamilyReoviridae, bs = 're') +
                         s(vFamilyRetroviridae, bs = 're') +
                         s(vFamilyRhabdoviridae, bs = 're') +
                         s(vFamilyTogaviridae, bs = 're'),
                       data = data_set,
                       select=TRUE,
                       family = gaussian)

summary(country_gam_mort)
gam.check(country_gam_mort)
AIC(country_gam_mort) 

#selected country-specific CFR gam
country_gam_mort <- gam(CFR_country_specific ~ 
                         # s(weighted_GDP, k=7, bs='tp') +
                        #  s(phylo_dist, bs = 're') +
                         # s(ReservoirPublicationCount, k=7, bs="tp") +
                        #  s(VirusFamilyPublicationCount, k=7, bs="tp") +
                          
                          #spillover type random effects
                          s(spill_typebridged, bs = 're') +
                         # s(spill_typedirect, bs = 're') +       
                          
                          #vector-borne random effects
                        #  s(IsVectorBorne0, bs = 're') +
                          s(IsVectorBorne1, bs = 're') +
                          
                          #host order random effects
                          s(hOrderAVES, bs = 're') +
                          s(hOrderCARNIVORA, bs = 're') +
                          s(hOrderCETARTIODACTYLA, bs = 're') +
                          s(hOrderCHIROPTERA, bs = 're') +
                        #  s(hOrderDIPROTODONTIA, bs = 're') +
                        #  s(hOrderEULIPOTYPHLA, bs = 're') +
                        #  s(hOrderPERISSODACTYLA, bs = 're') +
                        #  s(hOrderPRIMATES, bs = 're') +
                        #  s(hOrderRODENTIA, bs = 're') +
                          
                          #virus family random effects
                          s(vFamilyArenaviridae, bs = 're') +
                          s(vFamilyBornaviridae, bs = 're') +
                        #  s(vFamilyBunyaviridae, bs = 're') +
                          s(vFamilyCoronaviridae, bs = 're') +
                         # s(vFamilyFiloviridae, bs = 're') +
                        #  s(vFamilyFlaviviridae, bs = 're') +
                        #  s(vFamilyHepeviridae, bs = 're') +
                        #  s(vFamilyHerpesviridae, bs = 're') +
                        #  s(vFamilyOrthomyxoviridae, bs = 're') +
                          s(vFamilyParamyxoviridae, bs = 're') +
                        #  s(vFamilyPeribunyaviridae, bs = 're') +
                        #  s(vFamilyPicornaviridae, bs = 're') +
                          s(vFamilyPoxviridae, bs = 're') +
                      #    s(vFamilyReoviridae, bs = 're') +
                          s(vFamilyRetroviridae, bs = 're') +
                          s(vFamilyRhabdoviridae, bs = 're'), # +
                          #s(vFamilyTogaviridae, bs = 're'),
                        data = data_set,
                        select=TRUE,
                        family = gaussian)

summary(country_gam_mort)
gam.check(country_gam_mort)
AIC(country_gam_mort)  
#nope, GDPPC drops--so let's switch to the global CFR estimates

#global model w/ automated term selection: global CFR estimates
global_gam_mort <- gam(CFR_avg ~ 
                         s(phylo_dist, bs = 're') +
                           s(ReservoirPublicationCount, k=7, bs="tp") +
                           s(VirusFamilyPublicationCount, k=7, bs="tp") +
                         
                         #spillover type random effects
                         s(spill_typebridged, bs = 're') +
                         #s(spill_typedirect, bs = 're') +       
                         
                         #vector-borne random effects
                         #s(IsVectorBorne0, bs = 're') +
                         s(IsVectorBorne1, bs = 're') +
                          
                          #host order random effects
                          s(hOrderAVES, bs = 're') +
                          s(hOrderCARNIVORA, bs = 're') +
                          s(hOrderCETARTIODACTYLA, bs = 're') +
                          s(hOrderCHIROPTERA, bs = 're') +
                            s(hOrderDIPROTODONTIA, bs = 're') +
                            # s(hOrderEULIPOTYPHLA, bs = 're') +
                            # s(hOrderPERISSODACTYLA, bs = 're') +
                             s(hOrderPRIMATES, bs = 're') +
                             s(hOrderRODENTIA, bs = 're') +
                          
                          #virus family random effects
                          s(vFamilyArenaviridae, bs = 're') +
                          s(vFamilyBornaviridae, bs = 're') +
                           s(vFamilyBunyaviridae, bs = 're') +
                          s(vFamilyCoronaviridae, bs = 're') +
                            s(vFamilyFiloviridae, bs = 're') +
                             s(vFamilyFlaviviridae, bs = 're') +
                           #  s(vFamilyHepeviridae, bs = 're') +
                          s(vFamilyHerpesviridae, bs = 're') +
                             s(vFamilyOrthomyxoviridae, bs = 're') +
                          s(vFamilyParamyxoviridae, bs = 're') +
                           #  s(vFamilyPeribunyaviridae, bs = 're') +
                             s(vFamilyPicornaviridae, bs = 're') +
                          s(vFamilyPoxviridae, bs = 're') +
                             s(vFamilyReoviridae, bs = 're') +
                          s(vFamilyRetroviridae, bs = 're') +
                          s(vFamilyRhabdoviridae, bs = 're') +
                            s(vFamilyTogaviridae, bs = 're'),
                        data = data_set,
                        select=TRUE,
                        family = gaussian)

summary(global_gam_mort)
gam.check(global_gam_mort)
AIC(global_gam_mort) 

#selected model
selected_gam_mort <- gam(CFR_avg ~ 
                         #  s(phylo_dist, bs = 're') +
                         #  s(ReservoirPublicationCount, k=7, bs="tp") +
                         #  s(VirusFamilyPublicationCount, k=7, bs="tp") +
                           
                           #spillover type random effects
                           s(spill_typebridged, bs = 're') +
                        #   s(spill_typedirect, bs = 're') +       
                           
                           #vector-borne random effects
                        #   s(IsVectorBorne0, bs = 're') +
                           s(IsVectorBorne1, bs = 're') +
                           
                          #host order random effects
                           s(hOrderAVES, bs = 're') +
                           s(hOrderCARNIVORA, bs = 're') +
                           s(hOrderCETARTIODACTYLA, bs = 're') +
                           s(hOrderCHIROPTERA, bs = 're') +
                         #  s(hOrderDIPROTODONTIA, bs = 're') +
                        #   s(hOrderEULIPOTYPHLA, bs = 're') +
                        #   s(hOrderPERISSODACTYLA, bs = 're') +
                        #   s(hOrderPRIMATES, bs = 're') +
                        #   s(hOrderRODENTIA, bs = 're') +

                        #virus family random effects
                         s(vFamilyArenaviridae, bs = 're') +
                         s(vFamilyBornaviridae, bs = 're') +
                        # s(vFamilyBunyaviridae, bs = 're') +
                         s(vFamilyCoronaviridae, bs = 're') +
                       #  s(vFamilyFiloviridae, bs = 're') +
                      #   s(vFamilyFlaviviridae, bs = 're') +
                      #   s(vFamilyHepeviridae, bs = 're') +
                         s(vFamilyHerpesviridae, bs = 're') +
                      #   s(vFamilyOrthomyxoviridae, bs = 're') +
                         s(vFamilyParamyxoviridae, bs = 're') +
                      #   s(vFamilyPeribunyaviridae, bs = 're') +
                      #   s(vFamilyPicornaviridae, bs = 're') +
                         s(vFamilyPoxviridae, bs = 're') +
                      #   s(vFamilyReoviridae, bs = 're') +
                         s(vFamilyRetroviridae, bs = 're') +
                         s(vFamilyRhabdoviridae, bs = 're'),# +
                     #    s(vFamilyTogaviridae, bs = 're'),
                         data = data_set,
                        select=TRUE,
                         family = gaussian)

summary(selected_gam_mort)
gam.check(selected_gam_mort)
AIC(selected_gam_mort) 

#### FIGURES ####
plotData <- get_partial_effects_binary(fit = selected_gam_mort, 
                                       var = c('spill_typebridged',
                                               'IsVectorBorne1',
                                               'hOrderAVES',
                                               'hOrderCARNIVORA',
                                               'hOrderCETARTIODACTYLA',
                                               'hOrderCHIROPTERA',
                                               'vFamilyArenaviridae',
                                               'vFamilyBornaviridae',
                                               'vFamilyCoronaviridae',
                                               'vFamilyHerpesviridae',
                                               'vFamilyParamyxoviridae',
                                               'vFamilyRetroviridae',
                                               'vFamilyRhabdoviridae'),
                                       seWithMean = TRUE, fixedEffect = FALSE, removeNegatives = TRUE) 


gam_mort_labs <- c(
                   "Aves",
                    "Carnivora",
                     "Cetartiodactyla",
                      "Chiroptera",
                   
                   "vector-borne",
                   'bridged spillover',
                       "Arenaviridae",
                       'Bornaviridae',
                      'Coronaviridae',
                     'Herpesviridae',
                    'Paramyxoviridae',
                     'Retroviridae',
               'Rhabdoviridae')

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
  labs(x = NULL, y = 'Effect on human CFR') +
  scale_x_discrete(labels = gam_mort_labs) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())


#pdf file
ggsave(file = "Extension/code_and_data/figures/Figure_1.pdf",
       units="mm",
       width=80,
       height=80,
       scale=2,
       dpi=300)




