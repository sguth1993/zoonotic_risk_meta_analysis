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
dat <- read.csv(file="Extension/code_and_data/data/stringent_data.csv", header = T, stringsAsFactors = F)

## convert predictors to the correct class
dat$SppName_ICTV_MSL2018b <- as.factor(dat$SppName_ICTV_MSL2018b)
dat$hOrder <- as.factor(dat$hOrder)
dat$phylo_dist <- as.factor(dat$phylo_dist)
dat$vFamily <- as.factor(dat$vFamily)
dat$IsVectorBorne <- as.factor(dat$IsVectorBorne)
dat$death_burden_since_1950 <- as.numeric(dat$death_burden_since_1950)
dat$vaccination_effects <- as.factor(dat$vaccination_effects)
dat$Tr.primary <- as.factor(dat$Tr.primary)

#convert spill_type to binary
dat$spill_type[dat$spill_type=="bridged"] <- 1
dat$spill_type[dat$spill_type=="direct"] <- 0
dat$spill_type <- as.factor(dat$spill_type)

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
                                  spill_type, 
                                  .keep_all = TRUE)
unique(dat.death$SppName_ICTV_MSL2018b)

#global model with automated term selection
gam_death <- gam(death_burden_since_1950 ~  
                   s(vFamily, bs = 're' ) + 
                   s(hOrder, bs = 're' ) + 
                   s(Tr.primary, bs = 're' ) + 
                   s(spill_type, bs = 're' ) +
                   s(VirusSppPublicationCount, k =7, bs ='tp'),
                 offset = scale(death_period),
                 data = dat.death,
                 select = TRUE,
                 family = nb())
summary(gam_death)
gam.check(gam_death)

#selected model
#host order and spillover type dropped
gam_death <- gam(death_burden_since_1950 ~  
                   s(vFamily, bs = 're' ) + 
                  # s(hOrder, bs = 're' ) + 
                   s(Tr.primary, bs = 're' ) + 
                  # s(spill_type, bs = 're' ) +
                   s(VirusSppPublicationCount, k =7, bs ='tp'),
                 offset = scale(death_period),
                 data = dat.death,
                 select = TRUE,
                 family = nb())
summary(gam_death)
gam.check(gam_death)

#### FIGURES ####
plotData <- get_partial_effects(gam_death,
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
  labs(x = NULL, y = 'Strength of effect on death burden') +
  #scale_x_discrete(labels = gam_mort_labs) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 14),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())
print(vFamily)

#pdf file
ggsave(file = "Extension/code_and_data/figures/death_vFamily.pdf",
       units="mm",
       width=80,
       height=80,
       scale=2,
       dpi=300)

## transmission route
plotData_Tr.primary <- get_partial_effects(gam_death,
                                var = 'Tr.primary')

routes <- c("biting",
            "bodily fluids",
            "direct contact",
            "fecal-oral",
            "aerosolized animal excreta",
            "respiratory droplets",
            "unknown",
            "vector-borne")

route <- ggplot(plotData_Tr.primary$effects, aes(x = Tr.primary, y = y, colour = IsSignificant, fill = IsSignificant)) +
  geom_hline(yintercept = 0, linetype = 2, colour = 'grey20') +
  geom_boxplot(aes(middle = y, lower = ylower, upper = yupper, ymin = ylower, ymax = yupper), 
               stat = 'identity', alpha = 0.5, colour = NA) +
  geom_boxplot(aes(middle = y, lower = y, upper = y, ymin = y, ymax = y), 
               stat = 'identity', alpha = 0.5, size = 0.5) +
  
  geom_jitter(aes(y = Residual), alpha = 0.6, size = 0.8, width = 0.35, height = 0,
              data = plotData_Tr.primary$partialResiduals) +
  scale_colour_manual(values = c(No = 'grey60', Yes = '#225ea8'), guide = F) +
  scale_fill_manual(values = c(No = 'grey60', Yes = '#225ea8'), guide = F) +
  scale_y_continuous(labels = function(x) sprintf('%1.1f', x) ) +
  labs(x = NULL, y = 'Strength of effect on death burden') +
  scale_x_discrete(labels = routes) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 14),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())
print(route)

#pdf file
ggsave(file = "Extension/code_and_data/figures/death_route.pdf",
       units="mm",
       width=80,
       height=80,
       scale=2,
       dpi=300)

## virus species publication count
plotData_pub_count <- get_partial_effects_continuous(gam_death, 
                               var = 'VirusSppPublicationCount')

pub_count <- ggplot(plotData_pub_count$effects, aes(x = log(VirusSppPublicationCount), y = y, colour = IsSignificant, fill = IsSignificant)) +
  geom_hline(yintercept = 0, linetype = 2, colour = 'grey20') +
  geom_jitter(aes(y = Residual), alpha = 0.6, size = 0.8,
              width = 3, height = 0,
              data = plotData_pub_count$partialResiduals) +
  geom_ribbon(aes(ymin = ylower, ymax = yupper, colour = NA), alpha = 0.5) +
  #geom_rug(aes(y = Residual), colour = 'grey80', data = plotData_pub_count$partialResiduals) +
  geom_line(size = 1) +
  scale_colour_manual(values = c(No = 'grey60', Yes = '#225ea8'), guide = F) +
  scale_fill_manual(values = c(No = 'grey60', Yes = '#225ea8'), guide = F) +
  scale_y_continuous(labels = function(x) sprintf('%1.1f', x) ) +
  labs(x = NULL, y = 'Strength of effect on death burden') +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())
print(pub_count)

#pdf file
ggsave(file = "Extension/code_and_data/figures/death_pub_count.pdf",
       units="mm",
       width=80,
       height=80,
       scale=2,
       dpi=300)

Fig3 <- cowplot::plot_grid(vFamily, route, pub_count, 
                           ncol=3, nrow=1, align = "h",
                           rel_widths = c(2,1,1)) 
print(Fig3)

ggsave(file = "Extension/code_and_data/figures/Figure_3.pdf",
       units="mm",
       width=170,
       height=80,
       scale=2,
       dpi=300)



## death burden by order, raw data
orders <- c("Aves", 
            "Carnivora",
            "Cetartiodactyla",
            "Chiroptera",
            "Diprotodontia",
            "Eulipotyphla",
            "Perissodactyla",
            "Pilosa",
            "Primates",
            "Rodentia")

dat.death$deaths_per_year <- dat.death$death_burden_since_1950/dat.death$death_period

ggplot(dat.death, aes(x = hOrder, y = deaths_per_year)) +
  geom_boxplot() +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))  +
  scale_x_discrete(labels = orders) +
  labs(x= NULL, y = "deaths per year") +
  #scale_y_continuous(trans='log10')  +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 14),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())

#pdf file
ggsave(file = "Extension/code_and_data/figures/deaths_per_year_hOrder.pdf",
       units="mm",
       width=80,
       height=80,
       scale=2,
       dpi=300)
