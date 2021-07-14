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
library(zoib)
library("ggpubr")


setwd("~/Documents/Berkeley/Virulence_traits")
rm(list=ls())

## GDPPC wasn't significant in the country-specific GAM
#so let's use global estimates

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

#how many viruses have multiple entries?
dup <- dat.mort[duplicated(dat.mort$SppName_ICTV_MSL2018b)==TRUE, ] #only 9 viruses
dup_dat <- dat.mort %>%
  filter(SppName_ICTV_MSL2018b %in% dup$SppName_ICTV_MSL2018b)

#convert CFRs to proportions to model as beta distribution
dat.mort$CFR <- dat.mort$CFR_avg/100

#transformation for extreme 0 and 1 values
dat.mort$CFR <- (dat.mort$CFR*(86-1) + 0.5)/86

#global model w/ global CFR estimates
gam_mort <- gam(CFR ~
                  s(vFamily, bs = 're'  ) +
                  s(phylo_dist, k=5, bs="tp") +
                  s(ReservoirNspecies, k=7, bs="tp") +
                  s(hOrder, bs = 're'  ) +
                  s(VirusSppPublicationCount, k=7, bs="tp") +
                  s(spill_type, bs = 're' )  +
                  s(IsVectorBorne, bs = 're'),
                data = dat.mort,
                select=FALSE,
                family = betar(link = "logit"))
summary(gam_mort)
#gam.check(gam_mort)
AIC(gam_mort)

#### model selection ####
Terms <- list(
  vFamily = c(NA, "s(vFamily, bs = 're' )"),
  phylo_dist = c(NA, "s(phylo_dist, k=5, bs='tp')"),
  ReservoirNspecies = c(NA, "s(ReservoirNspecies, k=7, bs='tp')"),
  hOrder = c(NA,  "s(hOrder, bs = 're' )"), 
  spill_type = c(NA, "s(spill_type, bs = 're' )"), 
  IsVectorBorne = c(NA, "s(IsVectorBorne, bs = 're' )"),
  VirusSppPublicationCount = c(NA, "s(VirusSppPublicationCount, k=7, bs='tp')"))

## All possible combinations of these terms:
CompetingFullModels <- expand.grid(Terms)


## Build formulas
CompetingFullModels <- apply(CompetingFullModels, 1, function(row) paste(na.omit(row), collapse = ' + ')) %>% 
  data.frame(Formula = ., stringsAsFactors = F) %>% 
  mutate(Formula = ifelse(nchar(Formula) == 0, '1', Formula),
         Formula = paste('CFR ~', Formula))


## Model fit
CompetingFits <- CompetingFullModels %>% 
  group_by(Formula) %>% 
  do(ModelFit = try(gam(as.formula(.$Formula),
                        family = betar(link = "logit"),
                        data = dat.mort, 
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

gam_mort <- gam(as.formula(RankedModels$Formula[1]),
                  data = dat.mort,
                family = betar(link = "logit"),
                method = "REML")
summary(gam_mort)
gam.check(gam_mort)
AIC(gam_mort)

effects_dat <- plot.gam(gam_mort, residuals=TRUE, select = 2)[[2]]
effects <- as.factor(levels(effects_dat$raw))
effects <- cbind.data.frame(effects, effects_dat$fit)
print(effects)

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Find models to display:
#  - Top 15
#  - Top 5 models following the best one containing separate virus and reservoir effects 
#    (i.e. not containing family-specific reservoir effect)
#  - Top 5 models not accounting for publication bias
library(tidyr)
TOP_N_MODELS <- 15

TERM_NAMES <- tribble(~ Term,                   ~PlotLabel,             ~TermCategory,
                      'vFamily',          'virus family',                '',
                      'phylo_dist',       'phylogenetic distance',       'Reservoir effects',
                      'ReservoirNspecies', 'reservoir species richness', '',
                      'hOrder',           'reservoir host order',        'Reservoir effects',
                      'IsVectorBorne',       'vector-borne transmission', '',
                      'spill_type',       'bridged spillover',                '',
                      'VirusSppPublicationCount', 'virus species publication count', '')



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

mort_models <- ggplot(rankedSummary, aes(x = PlotLabel, y = Rank, fill = PercentDevianceExplained)) +
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

print(mort_models)

#selected model
gam_mort <- gam(as.formula(RankedModels$Formula[1]),
                data = dat.mort,
                select=FALSE,
                family = betar(link = "logit"))
summary(gam_mort)
#plot(gam_mort)
gam.check(gam_mort)
AIC(gam_mort)

#### FIGURES ####

## host order
plotData <- get_partial_effects(fit = gam_mort,
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

hOrder <- ggplot(plotData$effects, aes(x = hOrder, y = y, colour = IsSignificant, fill = IsSignificant)) +
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

print(hOrder)

## virus family
plotData <- get_partial_effects(gam_mort,
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
  labs(x = NULL, y = 'Effect on CFR in humans') +
  #scale_x_discrete(labels = gam_mort_labs) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 14),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())

print(vFamily)

#vector-borne and bridged spillover
plotData <- get_partial_effects_binary(gam_mort, 
                                       c('IsVectorBorne', 'spill_type'),
                                       fixedEffect = FALSE,
                                       removeNegatives = TRUE)


binary <- ggplot(plotData$effects, aes(x = variable, y = y, colour = IsSignificant, fill = IsSignificant)) +
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
  scale_x_discrete(labels = c("vector-borne transmission", "bridged spillover")) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())

print(binary)


Fig1 <- cowplot::plot_grid(mort_models, hOrder, vFamily, binary,
                           ncol=4, nrow=1, align = "h", 
                           labels = c( "A", "B", "C", "D"),
                           label_size = 20)
print(Fig1)

Fig1 <- ggdraw(Fig1) + 
  draw_plot_label(
    c("AIC", "Deviance\nexplained", "Phylogenetic distance from Primates"),
    c(0.196,0.21, 0.288),
    c(0.94,0.98, 0.999),
    size = 8
  )
print(Fig1)

ggsave(file = "Extension/code_and_data/figures/Figure_1.pdf",
       units="mm",
       width=240,
       height=90,
       scale=2,
       dpi=300)


#### predictions ####
#build empty dataset to house predictions
#include all predictor variables in top selected model
gam.dat <- dat.mort %>%
  group_by(hOrder) %>%
  summarise(Nobs=length(hOrder))
  
gam.dat$hOrder <- as.factor(gam.dat$hOrder)

gam.dat$spill_type <- 0
gam.dat$spill_type <- as.factor(gam.dat$spill_type)
gam.dat$IsVectorBorne <- 0
gam.dat$IsVectorBorne <- as.factor(gam.dat$IsVectorBorne)
#fill in a place holder for virus family but then exclude it in the predictions
gam.dat$vFamily <- "Poxviridae"
gam.dat$CFR = NA

gam.dat$CFR <- predict.gam(gam_mort, newdata=gam.dat, exclude = "s(vFamily)", type = "response")

gam.dat$CFR_lci <- gam.dat$CFR - (1.96*predict(gam_mort, newdata=gam.dat, type="response", se.fit=TRUE)$se)
gam.dat$CFR_uci <-gam.dat$CFR + (1.96*predict(gam_mort, newdata=gam.dat, type="response", se.fit=TRUE)$se)
gam.dat$CFR_lci[gam.dat$CFR_lci<0] <- 0
gam.dat$CFR[gam.dat$CFR<0] <- 0
gam.dat$hOrder <- str_to_title(gam.dat$hOrder)

gam.dat <- arrange(gam.dat, desc(CFR))

#gam.dat$hOrder <- factor(gam.dat$hOrder, levels=unique(gam.dat$hOrder))

gam.dat$CFR <- 100*gam.dat$CFR
gam.dat$CFR_lci <- 100*gam.dat$CFR_lci
gam.dat$CFR_uci <- 100*gam.dat$CFR_uci

colors= c('PRIMATES' ="firebrick1", 'PERISSODACTYLA'= "darkorange", 'CETARTIODACTYLA'=  "yellowgreen", 
          'CARNIVORA'="forestgreen",  'SCANDENTIA'="darkseagreen1" , 'LAGOMORPHA'="cornflowerblue", 
          'RODENTIA'=  "gray40",  'CHIROPTERA'="purple", 'DIPROTODONTIA'= "violet", 
          'PROBOSCIDEA'="magenta", 'EULIPOTYPHLA' ="maroon", 'PILOSA'="lightpink",
          'AVES'="navyblue")
gam.dat$hOrder <- as.factor(gam.dat$hOrder)
gam.dat$hOrder <- toupper(gam.dat$hOrder)
gam.dat$hOrder <- factor(gam.dat$hOrder, levels = c("PRIMATES",
                                                    "RODENTIA",
                                                    "CETARTIODACTYLA",
                                                    "PERISSODACTYLA",
                                                    "CARNIVORA",
                                                    "EULIPOTYPHLA",
                                                    "CHIROPTERA",
                                                    #"PILOSA",
                                                    "DIPROTODONTIA",
                                                    "AVES" ))


#and plot it...
p1 <- ggplot(data=gam.dat) +
  coord_cartesian(ylim = c(0,100), 
                  clip="off") +
  geom_phylopic(image=d$uid, alpha = 1, colour = 'gray10', x = d$x, y = -10) +
  geom_errorbar(aes(x=hOrder, ymin=CFR_lci, ymax=CFR_uci, color=hOrder),  width=0, linetype=3) +
  geom_point(aes(hOrder, CFR, color=hOrder, size=Nobs)) + 
  theme_bw() +
  geom_bracket(aes(xmin = "PRIMATES", xmax = "PRIMATES", label = "0"), data = gam.dat, y.position = 107, colour = "black") + 
  geom_bracket(aes(xmin = "RODENTIA", xmax = "RODENTIA", label = "180"), data = gam.dat, y.position = 107, colour = "black") + 
  geom_bracket(aes(xmin = "CETARTIODACTYLA", xmax = "CHIROPTERA", label = "193"), data = gam.dat, y.position = 107, colour = "black") + 
  geom_bracket(aes(xmin = "DIPROTODONTIA", xmax = "DIPROTODONTIA", label = "317"), data = gam.dat, y.position = 107, colour = "black") + 
  geom_bracket(aes(xmin = "AVES", xmax = "AVES", label = "624"), data = gam.dat, y.position = 107, colour = "black") + 
  scale_color_manual(values=colors, guide = "none") +
  labs(x = NULL, y = 'Predicted CFR in humans') +
  scale_x_discrete(labels = orders) +
  theme(plot.margin = unit(c(1,0,0,0), "cm"),
    panel.grid = element_blank(), 
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 14),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 7),
        legend.title = element_blank())

print(p1)

p1 <- ggdraw(p1) + 
  draw_plot_label(
    c("Phylogenetic distance from Primates"),
    c(0.188),
    c(0.999),
    size = 8
  )
print(p1)

ggsave(file = "Extension/code_and_data/figures/Figure_S2.pdf",
       units="mm",
       width=80,
       height=80,
       scale=2,
       dpi=300)

#### raw CFR data ####
length(unique(dat.mort$SppName_ICTV_MSL2018b[dat.mort$CFR_avg==0]))/
  length(unique(dat.mort$SppName_ICTV_MSL2018b))

length(unique(dat.mort$SppName_ICTV_MSL2018b[dat.mort$CFR_avg<10]))/
  length(unique(dat.mort$SppName_ICTV_MSL2018b))

length(unique(dat.mort$SppName_ICTV_MSL2018b[dat.mort$CFR_avg>50 & dat.mort$hOrder=="CHIROPTERA"]))/
  length(unique(dat.mort$SppName_ICTV_MSL2018b[dat.mort$CFR_avg>50]))


colors= c('PRIMATES' ="firebrick1", 'PERISSODACTYLA'= "darkorange", 'CETARTIODACTYLA'=  "yellowgreen", 
          'CARNIVORA'="forestgreen",  'SCANDENTIA'="darkseagreen1" , 'LAGOMORPHA'="cornflowerblue", 
          'RODENTIA'=  "gray40",  'CHIROPTERA'="purple", 'DIPROTODONTIA'= "violet", 
          'PROBOSCIDEA'="magenta", 'EULIPOTYPHLA' ="maroon", 'PILOSA'="lightpink",
          'AVES'="navyblue")

orders <- c("Aves", 
            "Carnivora",
            "Cetartiodactyla",
            "Chiroptera",
            "Diprotodontia",
            "Eulipotyphla",
            "Perissodactyla",
            "Primates",
            "Rodentia")


ggplot(dat.mort, aes(x = CFR_avg)) +
  geom_histogram(aes(fill = hOrder, colour = 'white'), binwidth =10, boundary = -0.5) +
  #stat_bin(aes(position=position_stack(vjust=0.5))) +
  scale_color_manual(values = "white", guide = FALSE) +
  scale_fill_manual(values =colors, labels = orders) +
  labs(x="human case fatality rate",
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
ggsave(file = "Extension/code_and_data/figures/Figure_S1.pdf",
       units="mm",
       width=70,
       height=70,
       scale=2,
       dpi=300)

