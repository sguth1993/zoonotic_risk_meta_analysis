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
dat$phylo_dist <- as.numeric(dat$phylo_dist)
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
dat.trans <- dat.trans[is.na(dat.trans$SppName_ICTV_MSL2018b)==FALSE, ]

#how many viruses have multiple entries?
dup <- dat.trans[duplicated(dat.trans$SppName_ICTV_MSL2018b)==TRUE, ] #only 9 viruses
dup_dat <- dat.trans %>%
  filter(SppName_ICTV_MSL2018b %in% dup$SppName_ICTV_MSL2018b)

#automated term selection by double penalty smoothing
gam_trans_global <- gam(human.trans ~
                          s(vFamily, bs = 're'  ) +
                          s(hOrder, bs = 're'  ) +
                          s(phylo_dist, k=5, bs='ts') +
                          s(ReservoirNspecies, k=7, bs='ts') +
                          s(VirusSppPublicationCount, k=7, bs='ts') +
                          s(spill_type, bs = 're' ) +
                          s(IsVectorBorne, bs = 're'),
                       data = dat.trans,
                       select=TRUE,
                 family = ocat(theta = c(1,2,3,4)))
summary(gam_trans_global)
#gam.check(gam_trans_global)

#record effect sizes
effects_dat <- plot.gam(gam_mort, residuals=TRUE, select = 2)[[2]]
effects <- as.factor(levels(effects_dat$raw))
effects_hOrd <- cbind.data.frame(effects, effects_dat$fit)
effects_hOrd <- effects_hOrd[order(effects_hOrd$`effects_dat$fit`, decreasing = TRUE), ]
print(effects_hOrd)

effects_dat <- plot.gam(gam_mort, residuals=TRUE, select = 7)[[7]]
effects <- as.factor(levels(effects_dat$raw))
effects_vec <- cbind.data.frame(effects, effects_dat$fit)
effects_vec$effects <- recode_factor(effects_vec$effects, "0" = "Direct transmission", "1" = "Vector-borne transmission")
effects_vec <- effects_vec[effects_vec$effects=="Vector-borne transmission", ]
print(effects_vec)

effect_sizes <- rbind(effects_hOrd, effects_vec)
write.csv(effect_sizes,"Extension/code_and_data/data/effect_sizes/trans_effects.csv")

#### AIC-maximization model selection ####
Terms <- list(
  vFamily = c(NA, "s(vFamily, bs = 're' )"),
  IsVectorBorne = c(NA, "s(IsVectorBorne, bs = 're' )"),
  phylo_dist = c(NA, "s(phylo_dist, k=5, bs='tp')"),
  ReservoirNspecies = c(NA, "s(ReservoirNspecies, k=7, bs='tp')"),
  hOrder = c(NA,  "s(hOrder, bs = 're' )"), 
  spill_type = c(NA, "s(spill_type, bs = 're' )"), 
  VirusSppPublicationCount = c(NA, "s(VirusSppPublicationCount, k=7, bs='tp')"))

## All possible combinations of these terms:
CompetingFullModels <- expand.grid(Terms)

## Build formulas
CompetingFullModels <- apply(CompetingFullModels, 1, function(row) paste(na.omit(row), collapse = ' + ')) %>% 
  data.frame(Formula = ., stringsAsFactors = F) %>% 
  mutate(Formula = ifelse(nchar(Formula) == 0, '1', Formula),
         Formula = paste('human.trans ~', Formula))

## Model fit
CompetingFits <- CompetingFullModels %>% 
  group_by(Formula) %>% 
  do(ModelFit = try(gam(as.formula(.$Formula),
                        data = dat.trans,
                        select=FALSE,
                        family = ocat(theta = c(1,2,3,4)))))

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

gam_tran <- gam(as.formula(RankedModels$Formula[1]),
                 data = dat.trans,
                 select = FALSE,
                 family = ocat(theta = c(1,2,3,4)))
summary(gam_tran)
#gam.check(gam_tran)
AIC(gam_tran)

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
                      'spill_type',       'bridged spillover',              'Virus effects',
                      'VirusSppPublicationCount', 'virus species publication count', 'Virus effects')




displayModels <- c(1:15) 

# Get deviance explained, etc for these models:
rankedSummary <- summarise_ranked_models_tran(RankedModels[displayModels, ], cores = 8)

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

tran_models <- ggplot(rankedSummary, aes(x = PlotLabel, y = Rank, fill = PercentDevianceExplained)) +
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


print(tran_models)

#selected model
gam_tran <- gam(as.formula(RankedModels$Formula[1]),
                data = dat.trans,
                select = FALSE,
                family = ocat(theta = c(1,2,3,4)))
summary(gam_tran)
gam.check(gam_tran)
AIC(gam_tran)



## virus family
plotData <- get_partial_effects(gam_tran,
                                var = 'vFamily')

fam <- ggplot(plotData$effects, aes(x = vFamily, y = y, colour = IsSignificant, fill = IsSignificant)) +
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
  labs(x = NULL, y = 'Effect on transmissibility\n within human populations') +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 14),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())
print(fam)

# phylogenetic distance
plotData_phylo_dist <- get_partial_effects_continuous(gam_tran, 
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
  labs(x = "phylogenetic distance from Primates", y = 'Effect on transmissibility\n within human populations') +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())
print(phylo_dist)


## vector-borne
plotData <- get_partial_effects_binary(gam_tran, 
                                       var = 'IsVectorBorne',
                                       seWithMean = TRUE, 
                                       fixedEffect = FALSE, removeNegatives = TRUE)


vector <- ggplot(plotData$effects, aes(x = IsVectorBorne, y = y, colour = IsSignificant, fill = IsSignificant)) +
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
  labs(x = NULL, y = 'Effect on transmissibility\n within human populations') +
  scale_x_discrete(labels = "vector-borne transmission") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())
print(vector)

## virus species publication count
plotData_pub_count <- get_partial_effects_continuous(gam_tran, 
                                                     var = 'VirusSppPublicationCount')

pub_count <- ggplot(plotData_pub_count$effects, aes(x = VirusSppPublicationCount, y = y, colour = IsSignificant, fill = IsSignificant)) +
  coord_cartesian(clip="off") +
  geom_hline(yintercept = 0, linetype = 2, colour = 'grey20') +
  geom_point(aes(x = VirusSppPublicationCount, y = Residual), alpha = 0.6, size = 0.8,
             data = plotData_pub_count$partialResiduals) +
  geom_ribbon(aes(ymin = ylower, ymax = yupper, colour = NA), alpha = 0.5) +
  geom_line(size = 1) +
  scale_colour_manual(values = c(No = 'grey60', Yes = 'mediumorchid4'), guide = F) +
  scale_fill_manual(values = c(No = 'grey60', Yes = 'mediumorchid4'), guide = F) +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 100, 500,1000,5000, 10000, 40000)) +
  labs(x = "virus species publication count", y = 'Effect on transmissibility\n within human populations') +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 11),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())
print(pub_count)

Fig2a <- cowplot::plot_grid(tran_models, fam, phylo_dist, 
                            ncol=3, nrow=1, align = "h",
                            labels = c( "A", "B", "C"),
                            label_size = 20)
print(Fig2a)

Fig2b <- cowplot::plot_grid(vector, pub_count, 
                            ncol=1, nrow=2, align = "v",
                            labels = c( "D", "E"),
                            label_size = 20,
                            hjust = c(-0.5,-0.5))
print(Fig2b)

Fig2 <- cowplot::plot_grid(Fig2a, Fig2b, 
                           ncol=2, nrow=1, align = "h",
                           rel_widths = c(3,1)) 
print(Fig2)

ggdraw(Fig2) + 
  draw_plot_label(
    c("AIC", "Deviance\nexplained"),
    c(0.196,0.21),
    c(0.94,0.98),
    size = 8
  )

ggsave(file = "Extension/code_and_data/figures/Figure_2.pdf",
       units="mm",
       width=240,
       height=90,
       scale=2,
       dpi=300)


#for supplement: model with hOrder instead of phylo_dist
gam_trans_hOrder <- gam(human.trans ~
                         s(vFamily, bs = 're'  ) +
                         s(IsVectorBorne, bs = 're'  ) +
                         s(hOrder, bs = 're'  ) +
                         s(VirusSppPublicationCount, k=7, bs="tp"), 
                       data = dat.trans,
                       select=FALSE,
                       family = ocat(theta = c(1,2,3,4)))
summary(gam_trans_hOrder)
#gam.check(gam_trans_hOrder)
AIC(gam_trans_hOrder) #93.83

#### FIGURES ####

## host order
plotData <- get_partial_effects(gam_trans_hOrder,
                                var = 'hOrder')
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

hOrder <- ggplot(plotData$effects, aes(x = hOrder, y = y, colour = IsSignificant, fill = IsSignificant)) +
  coord_cartesian(ylim = c(-6,6), 
                  clip="off") +
  geom_phylopic(image=d$uid, alpha = 1, colour = 'gray10', x = d$x, y = -7) +
  geom_hline(yintercept = 0, linetype = 2, colour = 'grey20') +
  geom_boxplot(aes(middle = y, lower = ylower, upper = yupper, ymin = ylower, ymax = yupper), 
               stat = 'identity', alpha = 0.5, colour = NA) +
  geom_boxplot(aes(middle = y, lower = y, upper = y, ymin = y, ymax = y), 
               stat = 'identity', alpha = 0.5, size = 0.5) +
  geom_jitter(aes(y = Residual), alpha = 0.6, size = 0.8, width = 0.35, height = 0,
              data = plotData$partialResiduals) +
  geom_bracket(aes(xmin = "PRIMATES", xmax = "PRIMATES", label = "0"), data = plotData$effects, y.position = 7, colour = "black") + 
  geom_bracket(aes(xmin = "RODENTIA", xmax = "RODENTIA", label = "180"), data = plotData$effects, y.position = 7, colour = "black") + 
  geom_bracket(aes(xmin = "CETARTIODACTYLA", xmax = "CHIROPTERA", label = "193"), data = plotData$effects, y.position = 7, colour = "black") + 
  geom_bracket(aes(xmin = "DIPROTODONTIA", xmax = "DIPROTODONTIA", label = "317"), data = plotData$effects, y.position = 7, colour = "black") + 
  geom_bracket(aes(xmin = "AVES", xmax = "AVES", label = "624"), data = plotData$effects, y.position = 7, colour = "black") + 
  scale_colour_manual(values = c(No = 'grey60', Yes = 'mediumorchid4'), guide = F) +
  scale_fill_manual(values = c(No = 'grey60', Yes = 'mediumorchid4'), guide = F) +
  scale_y_continuous(labels = function(x) sprintf('%1.1f', x) ) +
  labs(x = NULL, y = 'Effect on transmissibility\n within human populations') +
  scale_x_discrete(labels = orders) +
  theme_bw() +
  theme(plot.margin = unit(c(1,0,0,0), "cm"),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 14),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())

print(hOrder)

## virus family
plotData <- get_partial_effects(gam_trans_hOrder,
                                var = 'vFamily')

fam <- ggplot(plotData$effects, aes(x = vFamily, y = y, colour = IsSignificant, fill = IsSignificant)) +
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
  labs(x = NULL, y = 'Effect on transmissibility\n within human populations') +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 14),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())
print(fam)



## vector-borne
plotData <- get_partial_effects_binary(gam_trans_hOrder, 
                                       var = 'IsVectorBorne',
                                       seWithMean = TRUE, 
                                       fixedEffect = FALSE, removeNegatives = TRUE)


vector <- ggplot(plotData$effects, aes(x = IsVectorBorne, y = y, colour = IsSignificant, fill = IsSignificant)) +
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
  labs(x = NULL, y = 'Effect on transmissibility\n within human populations') +
  scale_x_discrete(labels = "vector-borne transmission") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())
print(vector)

## virus species publication count
plotData_pub_count <- get_partial_effects_continuous(gam_trans_hOrder, 
                                                     var = 'VirusSppPublicationCount')

pub_count <- ggplot(plotData_pub_count$effects, aes(x = VirusSppPublicationCount, y = y, colour = IsSignificant, fill = IsSignificant)) +
  coord_cartesian(clip="off") +
  geom_hline(yintercept = 0, linetype = 2, colour = 'grey20') +
  geom_point(aes(x = VirusSppPublicationCount, y = Residual), alpha = 0.6, size = 0.8,
             data = plotData_pub_count$partialResiduals) +
  geom_ribbon(aes(ymin = ylower, ymax = yupper, colour = NA), alpha = 0.5) +
  geom_line(size = 1) +
  scale_colour_manual(values = c(No = 'grey60', Yes = 'mediumorchid4'), guide = F) +
  scale_fill_manual(values = c(No = 'grey60', Yes = 'mediumorchid4'), guide = F) +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 100, 500,1000,5000, 10000, 40000)) +
  labs(x = "virus species publication count", y = 'Effect on transmissibility\n within human populations') +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 11),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())
print(pub_count)

Fig2a <- cowplot::plot_grid(hOrder, fam, 
                            ncol=2, nrow=1, align = "h",
                            labels = c( "A", "B"),
                            label_size = 20)
print(Fig2a)

Fig2b <- cowplot::plot_grid(vector, pub_count, 
                            ncol=1, nrow=2, align = "v",
                            labels = c( "C", "D"),
                            label_size = 20,
                            hjust = c(-0.5,-0.5))
print(Fig2b)

Fig2 <- cowplot::plot_grid(Fig2a, Fig2b, 
                           ncol=2, nrow=1, align = "h",
                           rel_widths = c(3,1)) 
print(Fig2)

ggsave(file = "Extension/code_and_data/figures/Figure_S10.pdf",
       units="mm",
       width=200,
       height=80,
       scale=2,
       dpi=300)


#### raw transmissibility data ####
unique(length(dat.trans$SppName_ICTV_MSL2018b[dat.trans$human.trans==1]))/
  unique(length(dat.trans$SppName_ICTV_MSL2018b))

unique(length(dat.trans$SppName_ICTV_MSL2018b[dat.trans$human.trans==4]))/
  unique(length(dat.trans$SppName_ICTV_MSL2018b))

unique(length(dat.trans$SppName_ICTV_MSL2018b[dat.trans$human.trans==4 & dat.trans$hOrder=="PRIMATES"]))/
  unique(length(dat.trans$SppName_ICTV_MSL2018b[dat.trans$human.trans==4]))

#grouped by host order
colors= c('AVES'="navyblue",
          'CARNIVORA'="forestgreen",
          'CETARTIODACTYLA'=  "yellowgreen", 
          'CHIROPTERA'="purple", 
          'DIPROTODONTIA'= "violet", 
          'EULIPOTYPHLA' ="maroon",
          'PERISSODACTYLA'= "darkorange",
          'PRIMATES' ="firebrick1",
          'RODENTIA'=  "gray40") 


orders <- c("Aves", 
            "Carnivora",
            "Cetartiodactyla",
            "Chiroptera",
            "Diprotodontia",
            "Eulipotyphla",
            "Perissodactyla",
            "Primates",
            "Rodentia")


ggplot(dat.trans, aes(x = human.trans)) +
  geom_histogram(aes(fill = hOrder, colour = 'white'), binwidth = 1) +
  scale_color_manual(values = "white", guide = FALSE) +
  scale_fill_manual(values =colors, labels = orders) +
  labs(x="transmissibility within human populations",
       y="number of virus species") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = c(0.78,0.78),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(.15,"cm"),
        legend.text = element_text(size = 15),
        legend.background = element_rect(fill="gray97", color="black")) 

#pdf file
ggsave(file = "Extension/code_and_data/figures/Figure_S9.pdf",
       units="mm",
       width=70,
       height=70,
       scale=2,
       dpi=300)


#### CFR vs. transmissibility tradeoff ####
library(ggbeeswarm)
ggplot(dat.trans, aes(x = CFR_avg, y = as.factor(human.trans))) +
  geom_point(aes(color = hOrder)) +
  geom_boxplot() +
  #geom_beeswarm(aes(color = hOrder)) +
  scale_color_manual(values =colors, labels = orders) +
  #labs(x= "transmissibility within human populations", y = "CFR in humans") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        # legend.position = c(0.78,0.78),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(.15,"cm"),
        legend.text = element_text(size = 15),
        legend.background = element_rect(fill="gray97", color="black"))

chisq.test(x = dat.trans$CFR_avg, y = dat.trans$human.trans,  simulate.p.value = TRUE)

ggplot(dat.trans, aes(x = CFR_avg, y = as.factor(human.trans))) +
  geom_boxplot() +
  geom_beeswarm(aes(color = hOrder), groupOnX = FALSE) +
  scale_color_manual(values =colors, labels = orders) +
  labs(x= "CFR in humans", y = "transmissibility within human populations") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        # legend.position = c(0.78,0.78),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(.15,"cm"),
        legend.text = element_text(size = 15),
        legend.background = element_rect(fill="gray97", color="black"))


#pdf file
ggsave(file = "Extension/code_and_data/figures/Figure_S11.pdf",
       units="mm",
       width=90,
       height=70,
       scale=2,
       dpi=300)


