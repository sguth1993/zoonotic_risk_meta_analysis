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
library(scales)
library(ggbeeswarm)
library(ggtree)
setwd("~/Documents/Berkeley/Virulence_traits")
#rm(list=ls())

## load data
dat <- read.csv(file="Extension/code_and_data/data/stringent_data.csv", header = T, stringsAsFactors = F)

## convert predictors to the correct class
dat$SppName_ICTV_MSL2018b <- as.factor(dat$SppName_ICTV_MSL2018b)
dat$hOrder <- as.factor(dat$hOrder)
dat$phylo_dist <- as.numeric(dat$phylo_dist)
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
dat.death <- dat.death[is.na(dat.death$SppName_ICTV_MSL2018b)==FALSE, ]

dat.death$VirusSppPublicationCount_log <- log(dat.death$VirusSppPublicationCount)
dat.death$ReservoirNspecies_log <- log(dat.death$ReservoirNspecies)

#automated term selection by double penalty smoothing
gam_death_global <- gam(death_burden_since_1950 ~  
                          s(vFamily, bs = 're'  ) +
                          s(hOrder, bs = 're'  ) +
                          s(phylo_dist, k=5, bs='ts') +
                          s(ReservoirNspecies, k=7, bs='ts') +
                          #s(VirusSppPublicationCount, k=7, bs='ts') +
                          s(spill_type, bs = 're' ) +
                          s(IsVectorBorne, bs = 're'),
                        offset = scale(death_period),
                        data = dat.death,
                        select = TRUE,
                        family = nb())
summary(gam_death_global)
#gam.check(gam_death_global)
AIC(gam_death_global)

#record effect sizes
effects_dat <- plot.gam(gam_mort, residuals=TRUE, select = 1)[[1]]
effects <- as.factor(levels(effects_dat$raw))
effects_vfam <- cbind.data.frame(effects, effects_dat$fit)
effects_vfam <- effects_vfam[order(effects_vfam$`effects_dat$fit`, decreasing = TRUE), ]
print(effects_vfam)

effects_dat <- plot.gam(gam_mort, residuals=TRUE, select = 5)[[5]]
effects <- as.factor(levels(effects_dat$raw))
effects_spill <- cbind.data.frame(effects, effects_dat$fit)
effects_spill$effects <- recode_factor(effects_spill$effects, "0" = "Direct spillover", "1" = "Bridged spillover")
effects_spill <- effects_spill[effects_spill$effects=="Bridged spillover", ]
print(effects_spill)

effect_sizes <- rbind(effects_vfam, effects_spill)
write.csv(effect_sizes,"Extension/code_and_data/data/effect_sizes/death_effects_no_pub.csv")


#### AIC-maximization model selection ####
Terms <- list(
  vFamily = c(NA, "s(vFamily, bs = 're' )"),
  IsVectorBorne = c(NA, "s(IsVectorBorne, bs = 're' )"),
  phylo_dist = c(NA, "s(phylo_dist, k=6, bs='tp')"),
  ReservoirNspecies = c(NA, "s(ReservoirNspecies, k=7, bs='tp')"),
  hOrder = c(NA,  "s(hOrder, bs = 're' )"), 
  spill_type = c(NA, "s(spill_type, bs = 're' )")) 

## All possible combinations of these terms:
CompetingFullModels <- expand.grid(Terms)

## Build formulas
CompetingFullModels <- apply(CompetingFullModels, 1, function(row) paste(na.omit(row), collapse = ' + ')) %>% 
  data.frame(Formula = ., stringsAsFactors = F) %>% 
  mutate(Formula = ifelse(nchar(Formula) == 0, '1', Formula),
         Formula = paste('death_burden_since_1950 ~', Formula))

## Model fit
CompetingFits <- CompetingFullModels %>% 
  group_by(Formula) %>% 
  do(ModelFit = try(gam(as.formula(.$Formula),
                        offset = scale(death_period),
                        family = nb(),
                        data = dat.death, 
                        select = FALSE,    # Not using internal model selection - mixing selection strategies makes presentation akward
                        method = 'REML')))

removeRows <- lapply(CompetingFits$ModelFit, function(x) 'try-error' %in% class(x)) %>% 
  unlist()
FailedFormulas <- CompetingFits[removeRows, ]
stopifnot(nrow(FailedFormulas) == 0)

CompetingFits <- CompetingFits[!removeRows, ]

## Add AIC:
RankedModels <- CompetingFits %>% 
  mutate(AIC = AIC(ModelFit),
         DevianceExplained = summary(ModelFit)$dev.expl) %>% 
  ungroup() %>% 
  arrange(AIC) %>% 
  mutate(DeltaAIC = AIC - AIC[1])

RankedModels$Formula[which.min(RankedModels$AIC)]

gam_death <- gam(as.formula(RankedModels$Formula[1]),
                 offset = scale(death_period),
                 data = dat.death,
                 select = FALSE,
                 family = nb())
summary(gam_death)
#gam.check(gam_death)
AIC(gam_death)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## ---- Plot models rankings: -----------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Find models to display:
#  - Top 15
library(tidyr)
TOP_N_MODELS <- 15

TERM_NAMES <- tribble(~ Term,                   ~PlotLabel,             ~TermCategory,
                      'vFamily',          'virus family',                'Virus effects',
                      'IsVectorBorne',    'vector-borne transmission',   'Virus effects',
                      'phylo_dist',       'phylogenetic distance',       'Reservoir effects',
                      'ReservoirNspecies', 'reservoir species richness', '',
                      'hOrder',           'reservoir host order',        'Reservoir effects',
                      'spill_type',       'bridged spillover',              'Virus effects')





displayModels <- c(1:8,10:15) 

# Get deviance explained, etc for these models:
rankedSummary <- summarise_ranked_models_death(RankedModels[displayModels, ], cores = 8)

# Add ranks:
ranks <- RankedModels %>% 
  select(Formula) %>% 
  mutate(Rank = row_number())

rankedSummary <- rankedSummary %>% 
  select(-Rank) %>% 
  left_join(ranks, by = 'Formula')


rankedSummary_sub <- rankedSummary[ ,c("Rank", "AIC", 
                                       #"PlotLabel" ,
                                       "PercentDevianceExplained",  "DevianceExplained",
                                       "sig", 'termClean', 'DeltaAIC')]
ninth <- data.frame(Rank = 9,
                    AIC = "",
                    #PlotLabel = "virus family",
                    PercentDevianceExplained=016.95879,
                    DevianceExplained=51.9,
                    sig = "Y",
                    termClean = "vFamily",
                    DeltaAIC = 2.258)

rankedSummary_plot <- rbind(rankedSummary_sub, ninth)

# Expand summary to include missing terms in each model:
# - full join ensures rows are created for missing terms
termDf <- lapply(1:15, function(x) data.frame(Rank = x, TERM_NAMES)) %>% 
  bind_rows()

rankedSummary_plot <- rankedSummary_plot %>% 
  full_join(termDf, by = c('termClean' = 'Term',
                           'Rank' = 'Rank'))



# Calculate model labels:
AICs <- rankedSummary_plot %>% 
  distinct(Rank, DeltaAIC, DevianceExplained) %>% 
  na.omit() %>% 
  mutate(DeltaAIC = sprintf("%.2f", DeltaAIC),
         DeltaAIC = str_pad(DeltaAIC, width = 5),
         DevianceExplained = sprintf("%.1f", DevianceExplained*100),
         Label = paste0(' ', DeltaAIC, '       ', DevianceExplained, '%'))

AICs$DevianceExplained[15] <- "54.6"
AICs$Label <- paste0(' ', AICs$DeltaAIC, '       ', AICs$DevianceExplained, '%')


AICvals <- AICs$Label
names(AICvals) <- AICs$Rank

# Plot:
rankedSummary_plot$sig[is.na(rankedSummary_plot$sig)==TRUE] <- "not_app"


death_models <- ggplot(rankedSummary_plot, aes(x = PlotLabel, y = Rank, fill = PercentDevianceExplained)) +
  geom_tile(aes(colour = sig, width=0.7, height=0.7), size = 1) +
  scale_color_manual(values = c(N = 'grey60', Y = 'black', not_app = 'grey60'), guide = F) +
  scale_y_continuous(breaks = displayModels,
                     sec.axis = sec_axis(trans = ~.,
                                         breaks = 1:15,
                                         labels = AICvals),
                     trans = 'reverse', expand = expansion(add = c(0.1, 0.1))) +
  scale_x_discrete(expand = c(0.055, 0.055)) +
  scale_fill_viridis_c(option = 'viridis', direction = -1, na.value = 'white',
                       breaks = c(1, seq(5, 50, by = 5))) +
  labs(x = NULL, y = NULL, fill = 'Relative deviance\nexplained (%)') + 
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
        axis.ticks.y = element_blank(),
        legend.position = 'top', 
        legend.title = element_text(vjust = 0.8, hjust = 1),
        legend.key.height = unit(0.8, 'lines'),
        legend.key.width = unit(1.5, 'lines'),
        legend.margin = margin(0, 0, 0, 0))

print(death_models)


#### FIGURES ####

## virus family
plotData <- get_partial_effects(gam_death,
                                var = 'vFamily')

vFamily <- ggplot(plotData$effects, aes(x = vFamily, y = y, colour = IsSignificant, fill = IsSignificant)) +
  coord_cartesian(clip="off") +
  geom_hline(yintercept = 0, linetype = 2, colour = 'grey20') +
  geom_boxplot(aes(middle = y, lower = ylower, upper = yupper, ymin = ylower, ymax = yupper), 
               stat = 'identity', alpha = 0.5, colour = NA) +
  geom_boxplot(aes(middle = y, lower = y, upper = y, ymin = y, ymax = y), 
               stat = 'identity', alpha = 0.5, size = 0.5) +
  
  geom_jitter(aes(y = Residual), alpha = 0.6, size = 0.8, width = 0.35, height = 0,
              data = plotData$partialResiduals) +
  scale_colour_manual(values = c(No = 'grey60', Yes = 'mediumorchid4'), guide = F) +
  scale_fill_manual(values = c(No = 'grey60', Yes = 'mediumorchid4'), guide = F) +
  scale_y_continuous(labels = function(x) sprintf('%1.1f', x) ) +
  labs(x = NULL, y = 'Effect on death burden') +
  #scale_x_discrete(labels = gam_mort_labs) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 14),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())
print(vFamily)

# phylogenetic distance
plotData_phylo_dist <- get_partial_effects_continuous(gam_death,
                                                      var = 'phylo_dist',
                                                      resolution = 1,
                                                      seWithMean = TRUE)

phylo_dist <- ggplot(plotData_phylo_dist$effects, aes(x = phylo_dist, y = y, colour = IsSignificant, fill = IsSignificant)) +
  coord_cartesian(clip="off") +
  geom_hline(yintercept = 0, linetype = 2, colour = 'grey20') +
  geom_jitter(aes(y = Residual), alpha = 0.6, size = 0.8,
              width = 3, height = 0,
              data = plotData_phylo_dist$partialResiduals) +
  geom_ribbon(aes(ymin = ylower, ymax = yupper, colour = NA), alpha = 0.5) +
  # geom_rug(aes(y = Residual), colour = 'grey80', data = plotData_pub_count$partialResiduals) +
  geom_line(size = 1) +
  scale_colour_manual(values = c(No = 'grey60', Yes = 'mediumorchid4'), guide = F) +
  scale_fill_manual(values = c(No = 'grey60', Yes = 'mediumorchid4'), guide = F) +
  scale_y_continuous(labels = function(x) sprintf('%1.1f', x) ) +
  labs(x = "phylogenetic distance from primates", y = 'Effect on death burden') +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())
print(phylo_dist)

## reservoir species richness
plotData_sp_rich <- get_partial_effects_continuous(gam_death,
                                                   var = 'ReservoirNspecies')


sp_rich <- ggplot(plotData_sp_rich$effects, aes(x = ReservoirNspecies, y = y, colour = IsSignificant, fill = IsSignificant)) +
  coord_cartesian(clip="off") +
  geom_hline(yintercept = 0, linetype = 2, colour = 'grey20') +
  geom_jitter(aes(y = Residual), alpha = 0.6, size = 0.8,
              width = 3, height = 0,
              data = plotData_sp_rich$partialResiduals) +
  geom_ribbon(aes(ymin = ylower, ymax = yupper, colour = NA), alpha = 0.5) +
  # geom_rug(aes(y = Residual), colour = 'grey80', data = plotData_pub_count$partialResiduals) +
  geom_line(size = 1) +
  scale_colour_manual(values = c(No = 'grey60', Yes = 'mediumorchid4'), guide = F) +
  scale_fill_manual(values = c(No = 'grey60', Yes = 'mediumorchid4'), guide = F) +
  scale_y_continuous(labels = function(x) sprintf('%1.1f', x) ) +
  labs(x = "reservoir species richness", y = 'Effect on death burden') +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())
print(sp_rich)

#vector-borne
plotData <- get_partial_effects_binary_single(gam_death,
                                              c('IsVectorBorne'),
                                              fixedEffect = FALSE,
                                              removeNegatives = TRUE)


vector <- ggplot(plotData$effects, aes(x = variable, y = y, colour = IsSignificant, fill = IsSignificant)) +
  coord_cartesian(clip="off") +
  geom_hline(yintercept = 0, linetype = 2, colour = 'grey20') +
  geom_boxplot(aes(middle = y, lower = ylower, upper = yupper, ymin = ylower, ymax = yupper),
               stat = 'identity', alpha = 0.5, colour = NA) +
  geom_boxplot(aes(middle = y, lower = y, upper = y, ymin = y, ymax = y),
               stat = 'identity', alpha = 0.5, size = 0.5) +

  geom_jitter(aes(y = Residual), alpha = 0.6, size = 0.8, width = 0.35, height = 0,
              data = plotData$partialResiduals) +
  scale_colour_manual(values = c(No = 'grey60', Yes = 'mediumorchid4'), guide = F) +
  scale_fill_manual(values = c(No = 'grey60', Yes = 'mediumorchid4'), guide = F) +
  scale_y_continuous(labels = function(x) sprintf('%1.1f', x) ) +
  labs(x = NULL, y = 'Effect on death burden') +
  scale_x_discrete(labels = c("vector-borne transmission")) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())

print(vector)

#combine into one figure
Fig_S9a <- cowplot::plot_grid(death_models, vFamily,
                             ncol=2, nrow=1, align = "h",
                             labels = c("A", "B"),
                             label_size = 20)

Fig_S9b <- cowplot::plot_grid(phylo_dist, sp_rich, vector,
                              ncol=1, nrow=3, align = "h",
                              labels = c( "C", "D", "E"),
                              label_size = 20)
Fig_S9 <- cowplot::plot_grid(Fig_S9a, Fig_S9b,
                              ncol=2, nrow=1, align = "h",
                              rel_widths = c(3,1),
                              label_size = 20)
print(Fig_S9)

ggdraw(Fig_S9) + 
  draw_plot_label(
    c("AIC", "Deviance\nexplained"),
    c(0.311,0.323),
    c(0.938,0.972),
    size = 8
  )

ggsave(file = "Extension/code_and_data/figures/Figure_3.pdf",
       units="mm",
       width=200,
       height=110,
       scale=2,
       dpi=300)




