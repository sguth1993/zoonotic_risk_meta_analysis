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
dat$vFamily <- as.factor(dat$vFamily)
dat$phylo_dist <- as.factor(dat$phylo_dist)
dat$IsVectorBorne <- as.factor(dat$IsVectorBorne)

#convert spill_type to binary
dat$spill_type[dat$spill_type=="bridged"] <- 1
dat$spill_type[dat$spill_type=="direct"] <- 0
dat$spill_type <- as.factor(dat$spill_type)

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

#global model: use automated term selection by double penalty smoothing to identify variables to remove
#let's first see if GDPPC is significant 
country_gam_mort <- gam(CFR_country_specific ~ 
                         s(weighted_GDP, k=7, bs='tp') +
                          s(vFamily, bs = 're'  ) +
                          s(hOrder, bs = 're'  ) +
                          s(VirusSppPublicationCount, k=7, bs="tp") +
                          s(spill_type, bs = 're' ) +
                          s(IsVectorBorne, bs = 're'),
                       data = dat.mort,
                       select=TRUE,
                       family = gaussian)

summary(country_gam_mort) #weighted_GDP is the only variable that drops completely
gam.check(country_gam_mort)
AIC(country_gam_mort) 

#selected country-specific CFR gam
country_gam_mort <- gam(CFR_country_specific ~ 
                          #s(weighted_GDP, k=7, bs='tp') +
                          s(vFamily, bs = 're'  ) +
                          s(hOrder, bs = 're'  ) +
                          s(VirusSppPublicationCount, k=7, bs="tp") +
                          s(spill_type, bs = 're' ) +
                          s(IsVectorBorne, bs = 're'),
                        data = dat.mort,
                        select=TRUE,
                        family = gaussian)

summary(country_gam_mort)
gam.check(country_gam_mort)
AIC(country_gam_mort) 

#nope, GDPPC drops--so let's switch to the global CFR estimates

#global model w/ automated term selection: global CFR estimates
gam_mort <- gam(CFR_avg ~
                         s(vFamily, bs = 're'  ) +
                         s(hOrder, bs = 're'  ) +
                         s(VirusSppPublicationCount, k=7, bs="tp") +
                         s(spill_type, bs = 're' ) +
                         s(IsVectorBorne, bs = 're'),
                       data = dat.mort,
                       select=TRUE,
                       family = gaussian)
summary(gam_mort) #no variables dropped
gam.check(gam_mort)

#### FIGURES ####

## host order
plotData <- get_partial_effects(gam_mort,
                                var = 'hOrder')
orders <- c("Aves", 
            "Carnivora",
            "Cetartiodactyla",
            "Chiroptera",
            "Diprotodontia",
            "Eulipotyphla",
            "Perissodactyla",
            "Primates",
            "Rodentia")

hOrder <- ggplot(plotData$effects, aes(x = hOrder, y = y, colour = IsSignificant, fill = IsSignificant)) +
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
  labs(x = NULL, y = 'Effect on CFR in humans') +
  scale_x_discrete(labels = orders) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 14),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())

#pdf file
ggsave(file = "Extension/code_and_data/figures/CFR_hOrder.pdf",
       units="mm",
       width=80,
       height=80,
       scale=2,
       dpi=300)


## virus family
plotData <- get_partial_effects(global_gam_mort,
                                var = 'vFamily')

vFamily <- ggplot(plotData$effects, aes(x = vFamily, y = y, colour = IsSignificant, fill = IsSignificant)) +
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
  labs(x = NULL, y = 'Effect on CFR in humans') +
  #scale_x_discrete(labels = gam_mort_labs) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 14),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())

#pdf file
ggsave(file = "Extension/code_and_data/figures/CFR_vFamily.pdf",
       units="mm",
       width=80,
       height=80,
       scale=2,
       dpi=300)

#spill type and vector-borne
plotData <- get_partial_effects_binary(gam_mort, 
                                       c('IsVectorBorne',
                                         'spill_type'),
                                       fixedEffect = FALSE,
                                       removeNegatives = TRUE)


binary <- ggplot(plotData$effects, aes(x = variable, y = y, colour = IsSignificant, fill = IsSignificant)) +
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
  labs(x = NULL, y = 'Effect on CFR in humans') +
  scale_x_discrete(labels = c("vector-\nborne", "bridged\nspillover")) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())

print(binary)
ggsave(file = "Extension/code_and_data/figures/CFR_binary.pdf",
       units="mm",
       width=80,
       height=80,
       scale=2,
       dpi=300)

plotData_pub_count <- get_partial_effects_continuous(gam_mort, 
                                                     var = 'VirusSppPublicationCount')

## virus species publication count
pub_count <- ggplot(plotData_pub_count$effects, aes(x = log(VirusSppPublicationCount), y = y, colour = IsSignificant, fill = IsSignificant)) +
  geom_hline(yintercept = 0, linetype = 2, colour = 'grey20') +
  geom_jitter(aes(y = Residual), alpha = 0.6, size = 0.8,
              width = 3, height = 0,
              data = plotData_pub_count$partialResiduals) +
  geom_ribbon(aes(ymin = ylower, ymax = yupper, colour = NA), alpha = 0.5) +
 # geom_rug(aes(y = Residual), colour = 'grey80', data = plotData_pub_count$partialResiduals) +
  geom_line(size = 1) +
  scale_colour_manual(values = c(No = 'grey60', Yes = '#225ea8'), guide = F) +
  scale_fill_manual(values = c(No = 'grey60', Yes = '#225ea8'), guide = F) +
  scale_y_continuous(labels = function(x) sprintf('%1.1f', x) ) +
  labs(x = "log(virus species publication count)", y = 'Effect on CFR in humans') +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())
print(pub_count)

ggsave(file = "Extension/code_and_data/figures/CFR_pub_count.pdf",
       units="mm",
       width=80,
       height=80,
       scale=2,
       dpi=300)

Fig1a <- cowplot::plot_grid(hOrder, vFamily, 
                           ncol=2, nrow=1, align = "h",
                           labels = c( "A", "B"),
                           label_size = 20)
print(Fig1a)

Fig1b <- cowplot::plot_grid(binary, pub_count, 
                           ncol=1, nrow=2, align = "v",
                           labels = c( "C", "D"),
                           label_size = 20)
print(Fig1b)

Fig1 <- cowplot::plot_grid(Fig1a, Fig1b, 
                            ncol=2, nrow=1, align = "h",
                          rel_widths = c(3,1)) 
print(Fig1)

ggsave(file = "Extension/code_and_data/figures/Figure_1.pdf",
       units="mm",
       width=200,
       height=80,
       scale=2,
       dpi=300)


