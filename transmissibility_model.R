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
dat$human.trans <- as.numeric(dat$human.trans)
dat$IsVectorBorne <- as.factor(dat$IsVectorBorne)
dat$Tr.primary <- as.factor(dat$Tr.primary)

#convert spill_type to binary
dat$spill_type[dat$spill_type=="bridged"] <- 1
dat$spill_type[dat$spill_type=="direct"] <- 0
dat$spill_type <- as.factor(dat$spill_type)

# and viruses with only one human case recorded in history
dat_rare <- dat[dat$tot_cases==1, ]
dat_suf <- dat[dat$tot_cases>1, ]

#final set of model variables
dat.trans <- dat_suf %>% distinct(SppName_ICTV_MSL2018b, human.trans, hOrder, spill_type, .keep_all = TRUE)
length(unique(dat.trans$SppName_ICTV_MSL2018b))

#how many viruses have multiple entries?
dup <- dat.trans[duplicated(dat.trans$SppName_ICTV_MSL2018b)==TRUE, ] #only 9 viruses
dup_dat <- dat.trans %>%
  filter(SppName_ICTV_MSL2018b %in% dup$SppName_ICTV_MSL2018b)

#global model: use automated term selection by double penalty smoothing to identify variables to remove
gam_trans <- gam(human.trans ~
                         s(vFamily, bs = 're'  ) +
                         s(hOrder, bs = 're'  ) +
                         s(VirusSppPublicationCount, k=7, bs="tp") +
                         s(Tr.primary, bs = 're') +
                         s(spill_type, bs = 're' ),
                       data = dat.trans,
                       select=TRUE,
                 family = ocat(theta = c(1,2,3,4)))
summary(gam_trans)
gam.check(gam_trans)

#selected model
#only virus family drops--effects likely overwhelmed by transmission route and virus publication count
gam_trans <- gam(human.trans ~
                  # s(vFamily, bs = 're'  ) + 
                   s(hOrder, bs = 're'  ) +
                   s(VirusSppPublicationCount, k=7, bs="tp") +
                   s(Tr.primary, bs = 're') +
                   s(spill_type, bs = 're' ),
                 data = dat.trans,
                 select=TRUE,
                 family = ocat(theta = c(1,2,3,4)))
summary(gam_trans)
gam.check(gam_trans)

#### FIGURES ####

## host order
plotData <- get_partial_effects(gam_trans,
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
  labs(x = NULL, y = 'Effect on between-human\n transmissibility') +
  scale_x_discrete(labels = orders) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 14),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())
print(hOrder)

#pdf file
ggsave(file = "Extension/code_and_data/figures/Tran_hOrder.pdf",
       units="mm",
       width=80,
       height=80,
       scale=2,
       dpi=300)

## transmission route
plotData <- get_partial_effects(gam_trans,
                                var = 'Tr.primary')

routes <- c("biting",
            "bodily fluids",
            "direct contact",
            "fecal-oral",
            "aerosolized animal excreta",
            "respiratory droplets",
            "unknown",
            "vector-borne")


route <- ggplot(plotData$effects, aes(x = Tr.primary, y = y, colour = IsSignificant, fill = IsSignificant)) +
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
  labs(x = NULL, y = 'Effect on between-human\n transmissibility') +
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
ggsave(file = "Extension/code_and_data/figures/Tran_route.pdf",
       units="mm",
       width=80,
       height=80,
       scale=2,
       dpi=300)

## spillover type
plotData <- get_partial_effects_binary(gam_trans, 
                                       var = 'spill_type',
                                       seWithMean = TRUE, 
                                       fixedEffect = FALSE, removeNegatives = TRUE)


spill <- ggplot(plotData$effects, aes(x = spill_type, y = y, colour = IsSignificant, fill = IsSignificant)) +
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
  labs(x = NULL, y = 'Effect on between-human\n transmissibility') +
  scale_x_discrete(labels = "bridged spillover") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())
print(spill)

#pdf file
ggsave(file = "Extension/code_and_data/figures/Tran_spill.pdf",
       units="mm",
       width=80,
       height=80,
       scale=2,
       dpi=300)

## virus species publication count
plotData_pub_count <- get_partial_effects_continuous(gam_trans, 
                                                     var = 'VirusSppPublicationCount')

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
  labs(x = "log(virus species publication count)", y = 'Effect on between-human\n transmissibility') +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 14),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())
print(pub_count)

#pdf file
ggsave(file = "Extension/code_and_data/figures/Tran_pub_count.pdf",
       units="mm",
       width=80,
       height=80,
       scale=2,
       dpi=300)

Fig2a <- cowplot::plot_grid(hOrder, route, 
                            ncol=2, nrow=1, align = "h",
                            labels = c( "A", "B"),
                            label_size = 20)
print(Fig2a)

Fig2b <- cowplot::plot_grid(spill, pub_count, 
                            ncol=1, nrow=2, align = "v",
                            labels = c( "C", "D"),
                            label_size = 20,
                            hjust = c(-0.5,-0.5))
print(Fig2b)

Fig2 <- cowplot::plot_grid(Fig2a, Fig2b, 
                           ncol=2, nrow=1, align = "h",
                           rel_widths = c(3,1)) 
print(Fig2)

ggsave(file = "Extension/code_and_data/figures/Figure_2.pdf",
       units="mm",
       width=200,
       height=80,
       scale=2,
       dpi=300)
