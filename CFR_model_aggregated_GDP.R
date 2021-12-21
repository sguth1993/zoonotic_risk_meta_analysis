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
library("ggpubr")
setwd("~/Documents/Berkeley/Virulence_traits")
rm(list=ls())

## load data
dat <- read.csv(file="Extension/code_and_data/data/stringent_data.csv", header = T, stringsAsFactors = F)

## convert all categorical predictors to factors
dat$SppName_ICTV_MSL2018b <- as.factor(dat$SppName_ICTV_MSL2018b)
dat$hOrder <- as.factor(dat$hOrder)
dat$Tr.primary <- as.factor(dat$Tr.primary)
dat$vFamily <- as.factor(dat$vFamily)
dat$phylo_dist <- as.numeric(dat$phylo_dist)
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
dat.mort <- dat.mort[is.na(dat.mort$SppName_ICTV_MSL2018b)==FALSE, ]

#how many viruses have multiple entries?
dup <- dat.mort[duplicated(dat.mort$SppName_ICTV_MSL2018b)==TRUE, ] #only 9 viruses
dup_dat <- dat.mort %>%
  filter(SppName_ICTV_MSL2018b %in% dup$SppName_ICTV_MSL2018b)

#convert CFRs to proportions to model as beta distribution
dat.mort$CFR_country_specific <- dat.mort$CFR_country_specific/100

#transformation for extreme 0 and 1 values
dat.mort$CFR_country_specific <- (dat.mort$CFR_country_specific*(87-1) + 0.5)/87

#automated term selection by double penalty smoothing
gam_mort <- gam(CFR_country_specific ~
                  s(vFamily, bs = 're'  ) +
                  s(hOrder, bs = 're'  ) +
                  s(phylo_dist, k=5, bs='ts') +
                  s(ReservoirNspecies, k=7, bs='ts') +
                  s(VirusSppPublicationCount, k=7, bs='ts') +
                  s(spill_type, bs = 're' ) +
                  s(IsVectorBorne, bs = 're') +
                  s(weighted_GDP, k=7, bs='ts'),
                data = dat.mort,
                select=TRUE,
                family = betar(link = "logit"))
summary(gam_mort) 
#gam.check(gam_mort)

#record effect sizes
effects_dat <- plot.gam(gam_mort, residuals=TRUE, select = 1)[[1]]
effects <- as.factor(levels(effects_dat$raw))
effects_vfam <- cbind.data.frame(effects, effects_dat$fit)
effects_vfam <- effects_vfam[order(effects_vfam$`effects_dat$fit`, decreasing = TRUE), ]
print(effects_vfam)

effects_dat <- plot.gam(gam_mort, residuals=TRUE, select = 2)[[2]]
effects <- as.factor(levels(effects_dat$raw))
effects_hOrd <- cbind.data.frame(effects, effects_dat$fit)
effects_hOrd <- effects_hOrd[order(effects_hOrd$`effects_dat$fit`, decreasing = TRUE), ]
print(effects_hOrd)

effects_dat <- plot.gam(gam_mort, residuals=TRUE, select = 6)[[6]]
effects <- as.factor(levels(effects_dat$raw))
effects_spill <- cbind.data.frame(effects, effects_dat$fit)
effects_spill$effects <- recode_factor(effects_spill$effects, "0" = "Direct spillover", "1" = "Bridged spillover")
effects_spill <- effects_spill[effects_spill$effects=="Bridged spillover", ]
print(effects_spill)

effects_dat <- plot.gam(gam_mort, residuals=TRUE, select = 7)[[7]]
effects <- as.factor(levels(effects_dat$raw))
effects_vec <- cbind.data.frame(effects, effects_dat$fit)
effects_vec$effects <- recode_factor(effects_vec$effects, "0" = "Direct transmission", "1" = "Vector-borne transmission")
effects_vec <- effects_vec[effects_vec$effects=="Vector-borne transmission", ]
print(effects_vec)

effect_sizes <- rbind(effects_vfam, effects_hOrd, effects_spill, effects_vec)
write.csv(effect_sizes,"Extension/code_and_data/data/effect_sizes/CFR_effects_country_aggregated.csv")


#let's first see if GDPPC is significant 
#### AIC-maximization model selection ####
Terms <- list(
  weighted_GDP = c(NA, "s(weighted_GDP, k=7, bs='ts')"),
  vFamily = c(NA, "s(vFamily, bs = 're' )"),
  IsVectorBorne = c(NA, "s(IsVectorBorne, bs = 're' )"),
  phylo_dist = c(NA, "s(phylo_dist, k=5, bs='ts')"),
  ReservoirNspecies = c(NA, "s(ReservoirNspecies, k=7, bs='ts')"),
  hOrder = c(NA,  "s(hOrder, bs = 're' )"), 
  spill_type = c(NA, "s(spill_type, bs = 're' )"), 
  VirusSppPublicationCount = c(NA, "s(VirusSppPublicationCount, k=7, bs='ts')"))

## All possible combinations of these terms:
CompetingFullModels <- expand.grid(Terms)

## Build formulas
CompetingFullModels <- apply(CompetingFullModels, 1, function(row) paste(na.omit(row), collapse = ' + ')) %>% 
  data.frame(Formula = ., stringsAsFactors = F) %>% 
  mutate(Formula = ifelse(nchar(Formula) == 0, '1', Formula),
         Formula = paste('CFR_country_specific ~', Formula))

## Model fit
CompetingFits <- CompetingFullModels %>% 
  group_by(Formula) %>% 
  do(ModelFit = try(gam(as.formula(.$Formula),
                        data = dat.mort,
                        select=FALSE,
                        family = betar(link = "logit"))))

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


country_gam_mort <- gam(as.formula(RankedModels$Formula[1]),
                        data = dat.mort,
                        select=FALSE,
                        family = betar(link = "logit"))

summary(country_gam_mort) 
gam.check(country_gam_mort)
AIC(country_gam_mort) 



#GDP is not significant
#plot for the supplement

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## ---- Plot models rankings: -----------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Find models to display:
#  - Top 15
library(tidyr)
TOP_N_MODELS <- 15

TERM_NAMES <- tribble(~ Term,                   ~PlotLabel,             ~TermCategory,
                      'weighted_GDP', 'aggregated GDP',                    '',
                      'vFamily',          'virus family',                'Virus effects',
                      'IsVectorBorne',    'vector-borne transmission',   'Virus effects',
                      'phylo_dist',       'phylogenetic distance',       'Reservoir effects',
                      'ReservoirNspecies', 'reservoir species richness', '',
                      'hOrder',           'reservoir host order',        'Reservoir effects',
                      'spill_type',       'bridged spillover',              'Virus effects',
                      'VirusSppPublicationCount', 'virus species publication count', 'Virus effects')




displayModels <- c(1:15) 

# Get deviance explained, etc for these models:
rankedSummary <- summarise_ranked_models_mort(RankedModels[displayModels, ], cores = 8)

# Add ranks:
ranks <- RankedModels %>% 
  select(Formula) %>% 
  mutate(Rank = row_number())

rankedSummary <- rankedSummary %>% 
  select(-Rank) %>% 
  left_join(ranks, by = 'Formula')


# Expand summary to include missing terms in each model:
# - full join ensures rows are created for missing terms
termDf <- lapply(displayModels, function(x) data.frame(Rank = x, TERM_NAMES)) %>% 
  bind_rows()

rankedSummary <- rankedSummary %>% 
  full_join(termDf, by = c('termClean' = 'Term',
                           'Rank' = 'Rank'))

# Calculate model labels:
AICs <- rankedSummary %>% 
  distinct(Rank, DeltaAIC, DevianceExplained) %>% 
  na.omit() %>% 
  mutate(DeltaAIC = sprintf("%.2f", DeltaAIC),
         DeltaAIC = str_pad(DeltaAIC, width = 5),
         DevianceExplained = sprintf("%.1f", DevianceExplained*100),
         Label = paste0(' ', DeltaAIC, '       ', DevianceExplained, '%'))

AICvals <- AICs$Label
names(AICvals) <- AICs$Rank

# Plot:
rankedSummary$sig[is.na(rankedSummary$sig)==TRUE] <- "not_app"

CFR_country_models <- ggplot(rankedSummary, aes(x = PlotLabel, y = Rank, fill = PercentDevianceExplained)) +
  geom_tile(aes(colour = sig, width=0.7, height=0.7), size = 1) +
  scale_color_manual(values = c(N = 'grey60', Y = 'black', not_app = 'grey60'), guide = F) +
  scale_y_continuous(breaks = displayModels,
                     sec.axis = sec_axis(trans = ~.,
                                         breaks = displayModels,
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
print(CFR_country_models)

#### FIGURES ####
## host order
plotData <- get_partial_effects(fit = country_gam_mort,
                                var = 'hOrder',
                                seWithMean = TRUE)

plotData$effects$hOrder <- factor(plotData$effects$hOrder, levels = c("PRIMATES",
                                                                      "RODENTIA",
                                                                      "CETARTIODACTYLA",
                                                                      "PERISSODACTYLA",
                                                                      "CARNIVORA",
                                                                      "EULIPOTYPHLA",
                                                                      "CHIROPTERA",
                                                                      #"PILOSA",
                                                                      "DIPROTODONTIA",
                                                                      "AVES" ))
orders <- c("Primates       ",
            "Rodentia       ",
            "Cetartiodactyla       ",
            "Perissodactyla       ",
            "Carnivora       ",
            "Eulipotyphla       ",
            "Chiroptera       ",
            "Diprotodontia       ",
            "Aves       ")

d <- ggimage::phylopic_uid(c("Gorilla_gorilla",
                             "Mus_musculus_domesticus",
                             "Sus_scrofa",
                             "Equus_asinus",
                             "Ailurus_fulgens",
                             "Suncus_murinus",
                             "Acerodon_celebensis",
                             "Phascolarctos_cinereus",
                             "Passer_domesticus"))




d$order <- c("Aves", 
             "Carnivora",
             "Cetartiodactyla",
             "Chiroptera",
             #"Pilosa",
             "Diprotodontia",
             "Eulipotyphla",
             "Perissodactyla",
             "Primates",
             "Rodentia")

d$order <- ordered(d$order, levels = c("Primates",
                                       "Rodentia",
                                       "Cetartiodactyla",
                                       "Perissodactyla",
                                       "Carnivora",
                                       "Eulipotyphla",
                                       "Chiroptera",
                                       #"Pilosa",
                                       "Diprotodontia",
                                       "Aves"))
d$x <- 1:9
d$color <- 'black'

reservoirs <- unique(dat.mort$hOrder)                  
reservoirs <- ordered(reservoirs, levels = c("PRIMATES",
                                             "RODENTIA",
                                             "CETARTIODACTYLA",
                                             "PERISSODACTYLA",
                                             "CARNIVORA",
                                             "EULIPOTYPHLA",
                                             "CHIROPTERA",
                                             #"PILOSA",
                                             "DIPROTODONTIA",
                                             "AVES"))

hOrder_country <- ggplot(plotData$effects, aes(x = hOrder, y = y, colour = IsSignificant, fill = IsSignificant)) +
  coord_cartesian(ylim = c(-2.5,2.57), 
                  clip="off") +
  geom_phylopic(image=d$uid, alpha = 1, colour = 'gray10', x = d$x, y = -3) +
  geom_hline(yintercept = 0, linetype = 2, colour = 'grey20') +
  geom_boxplot(aes(middle = y, lower = ylower, upper = yupper, ymin = ylower, ymax = yupper), 
               stat = 'identity', alpha = 0.5, colour = NA) +
  geom_boxplot(aes(middle = y, lower = y, upper = y, ymin = y, ymax = y), 
               stat = 'identity', alpha = 0.5, size = 0.5) +
  geom_jitter(aes(y = Residual), alpha = 0.6, size = 0.8, width = 0.35, height = 0,
              data = plotData$partialResiduals) +
  geom_bracket(aes(xmin = "PRIMATES", xmax = "PRIMATES", label = "0"), data = plotData$effects, y.position = 2.9, colour = "black") + 
  geom_bracket(aes(xmin = "RODENTIA", xmax = "RODENTIA", label = "180"), data = plotData$effects, y.position = 2.9, colour = "black") + 
  geom_bracket(aes(xmin = "CETARTIODACTYLA", xmax = "CHIROPTERA", label = "193"), data = plotData$effects, y.position = 2.9, colour = "black") + 
  geom_bracket(aes(xmin = "DIPROTODONTIA", xmax = "DIPROTODONTIA", label = "317"), data = plotData$effects, y.position = 2.9, colour = "black") + 
  geom_bracket(aes(xmin = "AVES", xmax = "AVES", label = "624"), data = plotData$effects, y.position = 2.9, colour = "black") + 
  scale_colour_manual(values = c(No = 'grey60', Yes = 'mediumorchid4'), guide = F) +
  scale_fill_manual(values = c(No = 'grey60', Yes = 'mediumorchid4'), guide = F) +
  scale_y_continuous(labels = function(x) sprintf('%1.1f', x) ) +
  labs(x = NULL, y = 'Effect on CFR in humans') +
  scale_x_discrete(labels = orders) +
  theme_bw() +
  theme(plot.margin = unit(c(1,0,0,0), "cm"),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 14),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())

print(hOrder_country)

## virus family
plotData <- get_partial_effects(country_gam_mort,
                                var = 'vFamily')

vFamily_country <- ggplot(plotData$effects, aes(x = vFamily, y = y, colour = IsSignificant, fill = IsSignificant)) +
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
  labs(x = NULL, y = 'Effect on CFR in humans') +
  #scale_x_discrete(labels = gam_mort_labs) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 14),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())

print(vFamily_country)

#spill type and vector-borne
plotData <- get_partial_effects_binary(country_gam_mort, 
                                       c('IsVectorBorne',
                                         'spill_type'),
                                       fixedEffect = FALSE,
                                       removeNegatives = TRUE)


binary_country <- ggplot(plotData$effects, aes(x = variable, y = y, colour = IsSignificant, fill = IsSignificant)) +
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
  labs(x = NULL, y = 'Effect on CFR in humans') +
  scale_x_discrete(labels = c("vector-borne\ntransmission", "bridged\nspillover")) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())

print(binary_country)

Fig1 <- cowplot::plot_grid(CFR_country_models, hOrder_country, vFamily_country, binary_country,
                            ncol=4, nrow=1, align = "h",
                            labels = c("A", "B", "C", "D"),
                            label_size = 20)
print(Fig1)

Fig1 <- ggdraw(Fig1) + 
  draw_plot_label(
    c("AIC", "Deviance\nexplained", "Phylogenetic distance from Primates"),
    c(0.196,0.21, 0.28),
    c(0.94,0.98, 0.999),
    size = 8
  )
print(Fig1)

ggsave(file = "Extension/code_and_data/figures/Figure_S8.pdf",
       units="mm",
       width=240,
       height=90,
       scale=2,
       dpi=300)

