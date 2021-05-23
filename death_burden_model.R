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
library(ggrepel)
setwd("~/Documents/Berkeley/Virulence_traits")
rm(list=ls())

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

#global model
gam_death_global <- gam(death_burden_since_1950 ~  
                   s(vFamily, bs = 're' ) + 
                   s(IsVectorBorne, bs = 're' ) + 
                   s(phylo_dist, k=6, bs="tp") +
                    s(ReservoirNspecies_log, k=7, bs="tp") +
                   s(hOrder, bs = 're' ) + 
                   s(spill_type, bs = 're' ) +
                   s(VirusSppPublicationCount, k=7, bs ='tp'),
                 offset = scale(death_period),
                 data = dat.death,
                 select = FALSE,
                 family = nb())
summary(gam_death_global)
#gam.check(gam_death_global)
AIC(gam_death_global)

#### model selection ####
Terms <- list(
  vFamily = c(NA, "s(vFamily, bs = 're' )"),
  IsVectorBorne = c(NA, "s(IsVectorBorne, bs = 're' )"),
  phylo_dist = c(NA, "s(phylo_dist, k=6, bs='tp')"),
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
                      'spill_type',       'bridged spillover',              'Virus effects',
                      'VirusSppPublicationCount', 'virus species publication count', 'Virus effects')




displayModels <- c(1:15) 

# Get deviance explained, etc for these models:
rankedSummary <- summarise_ranked_models_death(RankedModels[displayModels, ], cores = 8)

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

death_models <- ggplot(rankedSummary, aes(x = PlotLabel, y = Rank, fill = PercentDevianceExplained)) +
  geom_tile(aes(colour = sig, width=0.7, height=0.7), size = 1) +
  scale_color_manual(values = c(N = 'grey60', Y = 'black', not_app = 'grey60'), guide = F) +
  scale_y_continuous(breaks = displayModels,
  									 sec.axis = sec_axis(trans = ~.,
  									 										breaks = displayModels,
  									 										labels = AICvals),
  									 trans = 'reverse', expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_viridis_c(option = 'viridis', direction = -1, na.value = 'white',
                       breaks = c(1, seq(5, 50, by = 5))) +
  labs(x = NULL, y = NULL, fill = 'Deviance explained (%)') + 
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

#selected model
gam_death <- gam(as.formula(RankedModels$Formula[1]),
                 offset = scale(death_period),
                 data = dat.death,
                 select = FALSE,
                 family = nb())
summary(gam_death)
#gam.check(gam_death)
AIC(gam_death)

#### FIGURES ####
## host order
plotData <- get_partial_effects(fit = gam_death,
                                var = 'hOrder',
                                seWithMean = TRUE)
orders <- c("Aves      ",
            "Carnivora      ",
            "Cetartiodactyla      ",
            "Chiroptera      ",
            "Diprotodontia      ",
            "Eulipotyphla      ",
            "Perissodactyla      ",
            "Pilosa      ",
            "Primates      ",
            "Rodentia      ")

d <- phylopic_uid(c("Passer_domesticus",
                             "Ailurus_fulgens",
                             "Sus_scrofa",
                             "Acerodon_celebensis",
                             "Phascolarctos_cinereus",
                             "Suncus_murinus",
                             "Equus_asinus",
                             "ground_sloths",
                             "Gorilla_gorilla",
                             "Mus_musculus_domesticus"))

d$order <- c("Aves",
             "Carnivora",
             "Cetartiodactyla",
             "Chiroptera",
             "Diprotodontia",
             "Eulipotyphla",
             "Perissodactyla",
             "Pilosa",
             "Primates",
             "Rodentia")
d$order <- as.factor(d$order)
d$x <- 1:10
d$color <- 'black'

reservoirs <- as.factor(c("AVES",
                          "CARNIVORA",
                          "CETARTIODACTYLA",
                          "CHIROPTERA",
                          "DIPROTODONTIA",
                          "EULIPOTYPHLA" ,
                          "PERISSODACTYLA",
                          "PILOSA",
                          "PRIMATES",
                          "RODENTIA"))


hOrder <- ggplot(plotData$effects, aes(x = hOrder, y = y, colour = IsSignificant, fill = IsSignificant)) +
  coord_cartesian(ylim = c(-2.5,2.5),
                  clip="off") +
  geom_phylopic(image=d$uid, alpha = 1, colour = 'gray10', x = d$x, y = -3) +
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
  scale_x_discrete(labels = orders) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 14),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())

print(hOrder)


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


## virus species publication count
plotData_pub_count <- get_partial_effects_continuous(gam_death, 
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
  labs(x = "virus species publication count", y = 'Effect on death burden') +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 11),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())
print(pub_count)

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

#spill type
plotData <- get_partial_effects_binary_single(gam_death,
                                              c('spill_type'),
                                              fixedEffect = FALSE,
                                              removeNegatives = TRUE)


spill_type <- ggplot(plotData$effects, aes(x = variable, y = y, colour = IsSignificant, fill = IsSignificant)) +
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
  scale_x_discrete(labels = c("bridged spillover")) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())

print(spill_type)

#combine into one figure
Fig3a <- cowplot::plot_grid(death_models,vFamily, hOrder,
                            ncol=3, nrow=1, align = "h",
                            labels = c( "A", "B", "C"),
                            label_size = 20)
print(Fig3a)

Fig3b <- cowplot::plot_grid(pub_count, spill_type, phylo_dist, sp_rich,
                            ncol=4, nrow=1, align = "h",
                            labels = c( "D", "E", "F", "G"),
                            label_size = 20)
print(Fig3b)

Fig3 <- cowplot::plot_grid(Fig3a, Fig3b,
                           ncol=1, nrow=2, align = "h",
                           rel_heights = c(3,2))

Fig3 <- ggdraw(Fig3) + 
  draw_plot_label(
    c("AIC", "Deviance\nexplained"),
    c(0.272,0.291),
    c(0.96,0.99),
    size = 8
  )
print(Fig3)

ggsave(file = "Extension/code_and_data/figures/Figure_S9.pdf",
       units="mm",
       width=210,
       height=110,
       scale=2,
       dpi=300)

#### plots of raw death data ####

## raw data
length(unique(dat.death$SppName_ICTV_MSL2018b[dat.death$death_burden_since_1950==0]))/
  length(unique(dat.death$SppName_ICTV_MSL2018b))

length(unique(dat.death$SppName_ICTV_MSL2018b[dat.death$death_burden_since_1950<50]))/
  length(unique(dat.death$SppName_ICTV_MSL2018b))

sum(unique(dat.death$death_burden_since_1950[dat.death$death_burden_since_1950>100000]))/sum(unique(dat.death$death_burden_since_1950))

length(unique(dat.death$SppName_ICTV_MSL2018b[dat.death$death_burden_since_1950>100000]))/length(unique(dat.death$SppName_ICTV_MSL2018b))

sum(unique(dat.death$death_burden_since_1950[dat.death$Tr.primary=="respiratory"]))/sum(unique(dat.death$death_burden_since_1950))

dat.death$death_burden_since_1950[dat.death$SppName_ICTV_MSL2018b=="Hantaan orthohantavirus"]/sum(dat.death$death_burden_since_1950[dat.death$vFamily=="Bunyaviridae" & 
                                    dat.death$hOrder=="RODENTIA" &
                                    dat.death$SppName_ICTV_MSL2018b!="Hantaan orthohantavirus"])

#raw data histogram for supplement
colors= c('AVES'="navyblue",
          'CARNIVORA'="forestgreen",
          'CETARTIODACTYLA'=  "yellowgreen", 
          'CHIROPTERA'="purple", 
          'DIPROTODONTIA'= "violet", 
          'EULIPOTYPHLA' ="maroon",
          'PERISSODACTYLA'= "darkorange",
          'PILOSA'="lightpink",
          'PRIMATES' ="firebrick1",
          'RODENTIA'=  "gray40") 


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


ggplot(dat.death, aes(x = death_burden_since_1950)) +
  geom_histogram(aes(fill = hOrder, colour = 'white'), boundary = -0.5) +
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), 
                     labels = scales::label_comma(accuracy = 1),
                     breaks = c(0,10,50,100,1000,100000,1000000,2000000)) +
  #stat_bin(aes(position=position_stack(vjust=0.5))) +
  scale_color_manual(values = "white", guide = FALSE) +
  scale_fill_manual(values =colors, labels = orders) +
  labs(x="post-1950 cumulative death count",
       y="number of virus species") +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1,size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = c(0.72,0.72),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(.15,"cm"),
        legend.text = element_text(size = 15),
        legend.background = element_rect(fill="gray97", color="black")) 

#pdf file
ggsave(file = "Extension/code_and_data/figures/Figure_S8.pdf",
       units="mm",
       width=70,
       height=70,
       scale=2,
       dpi=300)

## death burden per year by order
dat.death$deaths_per_year <- dat.death$death_burden_since_1950/dat.death$death_period

orders <- c("Aves      ", 
            "Carnivora      ",
            "Cetartiodactyla      ",
            "Chiroptera      ",
            "Diprotodontia      ",
            "Eulipotyphla      ",
            "Perissodactyla      ",
            "Pilosa    ",
            "Primates      ",
            "Rodentia      ")

d <- ggimage::phylopic_uid(c("Passer_domesticus",
                             "Ailurus_fulgens",
                             "Sus_scrofa",
                             "Acerodon_celebensis",
                             "Phascolarctos_cinereus",
                             "Suncus_murinus",
                             "Equus_asinus",
                             "ground_sloths",
                             "Gorilla_gorilla",
                             "Mus_musculus_domesticus"))

d$order <- c("Aves", 
             "Carnivora",
             "Cetartiodactyla",
             "Chiroptera",
             "Diprotodontia",
             "Eulipotyphla",
             "Perissodactyla",
             "Pilosa",
             "Primates",
             "Rodentia")
d$order <- as.factor(d$order)
d$x <- 1:10
d$color <- 'black'

reservoirs <- as.factor(c("AVES",
                          "CARNIVORA",
                          "CETARTIODACTYLA",
                          "CHIROPTERA",
                          "DIPROTODONTIA",
                          "EULIPOTYPHLA" , 
                          "PERISSODACTYLA",
                          "PILOSA",
                          "PRIMATES",
                          "RODENTIA"))


plot.death <- dat.death
plot.death <- dat.death %>% distinct(SppName_ICTV_MSL2018b, spill_type,
                                     .keep_all = TRUE)

plot.death$human.trans <- factor(plot.death$human.trans, levels = c('4','3','2','1'))
plot.death$deaths_per_year <- plot.death$death_burden_since_1950/plot.death$death_period

colors = c("midnightblue",
           "darkorchid3",
           #"#ff4d4d",
           "magenta",
           "cornflowerblue")

plot.death$names <- plot.death$SppName_ICTV_MSL2018b
plot.death$names <- as.character(plot.death$names)
plot.death <- plot.death[order(-plot.death$deaths_per_year), ]
plot.death$names[c(10:88)] <- NA

plot.death$names <- gsub(" virus", "", plot.death$names)
plot.death$names[plot.death$names=="Severe acute respiratory syndrome-related coronavirus-2"] <- "SARS-CoV-2"
plot.death$names[plot.death$names=="Rabies lyssavirus"] <- "Rabies"
plot.death$names[plot.death$names=="Ebolavirus"] <- "Ebola"
plot.death$names[plot.death$names=="Influenza A"] <- "Influenza A (bridged spillover)"

colors = c("midnightblue",
           "darkorchid3",
           #"#ff4d4d",
           "magenta",
           "cornflowerblue")

death_order <- ggplot(plot.death, aes(x = hOrder, y = deaths_per_year)) +
  coord_cartesian(ylim = c(0,2000000), 
                  clip="off") +
  geom_phylopic(data = d, image=d$uid, alpha = 1, colour = 'gray10', x = d$x, y = -0.55) +
  geom_boxplot() +
  geom_beeswarm(aes(color = human.trans)) +
  scale_color_manual(values = colors, labels = c("4 (high)", "3", "2", "1 (low)")) + 
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     labels = scales::label_comma(accuracy = 1),
                     breaks = c(0,10,100,1000,100000,1000000)) +
  scale_x_discrete(labels = orders) +
  labs(x= NULL, 
       y = "death burden per year",
       color = "transmissibility\nbetween humans") +
  theme_bw() +
  guides(ncol = 4,
         nrow = 1,
    col = guide_legend(title.position = "top", override.aes = list(size = 3))) +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1,size = 12),
        axis.text.y = element_text(size = 12),
        panel.grid = element_blank(),
        legend.title.align=0.5,
        legend.key.size = unit(.3,"cm"),
        legend.text = element_text(size = 12),
        legend.background = element_rect(fill="gray97", color="black"),
        #legend.position = "top",
        legend.position = c(0.8, 0.8)) 
print(death_order)

# death burden per year by virus family
death_fam <- ggplot(plot.death, aes(x = vFamily, y = deaths_per_year)) +
  coord_cartesian(clip="off") +
  geom_boxplot() +
  geom_beeswarm(aes(color = human.trans)) +
  scale_color_manual(values = colors) + 
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), 
                     labels = scales::label_comma(accuracy = 1),
                     breaks = c(0,10,100,1000,100000,1000000)) +
  labs(x= NULL, 
       y = "death burden per year",
       color = "virus family") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 14),
        axis.text.y = element_text(size = 10),
        legend.position = "none",
        panel.grid = element_blank())
print(death_fam)

#death burden per year by transmission route
routes <- c("biting",
            "bodily fluids",
            "direct contact",
            "fecal-oral",
            "inhalation of\naerosolized excreta",
            "respiratory",
            "unknown",
            "vector")
death_route <- ggplot(plot.death, aes(x = Tr.primary, y = deaths_per_year)) +
  coord_cartesian(clip="off") +
  geom_boxplot() +
  geom_beeswarm(aes(color = human.trans)) +
  scale_color_manual(values = colors) + 
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), 
                     labels = scales::label_comma(accuracy = 1),
                     breaks = c(0,10,100,1000,100000,1000000)) +
  scale_x_discrete(labels = routes) +
  labs(x= NULL, 
       y = "death burden per year",
       color = "transmissibility\nbetween humans") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 14),
        axis.text.y = element_text(size = 10),
        legend.position = "none",
        panel.grid = element_blank()) 

print(death_route)

CFR_death <- ggplot(plot.death, aes(x = CFR_avg, y = deaths_per_year, color = human.trans)) + 
  coord_cartesian(clip="off") +
  geom_point(size = 2) +
  scale_color_manual(values = colors) + 
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), 
                     labels = scales::label_comma(accuracy = 1),
                     breaks = c(0,10,100,1000,10000,100000,1000000,2000000)) +
  labs(x="CFR in humans",
       y="death burden per year",
       color = "transmissibility\nbetween humans") +
  geom_text_repel(aes(label=names),
                  show.legend = FALSE,
                  nudge_x = 0.1) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none",
        panel.grid = element_blank())
print(CFR_death)

## combine in cowplot
Fig4b <- cowplot::plot_grid(CFR_death,
                           ncol=1, nrow=1, align = "h")
                          
print(Fig4b)

Fig4a <- cowplot::plot_grid(death_route, Fig4b,
                           ncol=2, nrow=1, align = "h",
                           labels = c("C", "D"),
                           rel_heights = c(1,2),
                           #hjust = -0.1,
                           label_size = 20)
print(Fig4a)

Fig4 <- cowplot::plot_grid(death_order, death_fam, death_route, Fig4b,
                            ncol=2, nrow=2, align = "h",
                           axis = "b",
                            labels = c( "A","B", "C", "D"),
                            label_size = 20)
print(Fig4)

ggsave(file = "Extension/code_and_data/figures/Figure_4.pdf",
       units="mm",
       width=160,
       height=160,
       scale=2,
       dpi=300)

