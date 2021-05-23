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

#global model w/ global CFR estimates
gam_mort <- gam(CFR_avg ~
                  s(vFamily, bs = 're'  ) +
                  s(phylo_dist, k=5, bs="tp") +
                  s(ReservoirNspecies, k=7, bs="tp") +
                  s(hOrder, bs = 're'  ) +
                  s(Tr.primary, bs = 're'  ) +
                  s(VirusSppPublicationCount, k=7, bs="tp") +
                  s(spill_type, bs = 're' )  +
                  s(IsVectorBorne, bs = 're'),
                data = dat.mort,
                select=FALSE,
                family = gaussian)
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
         Formula = paste('CFR_avg ~', Formula))


## Model fit
CompetingFits <- CompetingFullModels %>% 
  group_by(Formula) %>% 
  do(ModelFit = try(gam(as.formula(.$Formula),
                        family = gaussian(),
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
                  family = gaussian)
summary(gam_mort)

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

print(mort_models)

#selected model
gam_mort <- gam(as.formula(RankedModels$Formula[1]),
                data = dat.mort,
                select=FALSE,
                family = gaussian)
summary(gam_mort)
#plot(gam_mort)
#gam.check(gam_mort)
AIC(gam_mort)

#### FIGURES ####

## host order
plotData <- get_partial_effects(fit = gam_mort,
                                var = 'hOrder',
                                seWithMean = TRUE)
orders <- c("Aves       ", 
            "Carnivora       ",
            "Cetartiodactyla       ",
            "Chiroptera       ",
            "Diprotodontia       ",
            "Eulipotyphla       ",
            "Perissodactyla       ",
            "Primates       ",
            "Rodentia       ")

d <- ggimage::phylopic_uid(c("Passer_domesticus",
                           "Ailurus_fulgens",
                           "Sus_scrofa",
                           "Acerodon_celebensis",
                           "Phascolarctos_cinereus",
                           "Suncus_murinus",
                           "Equus_asinus",
                          "Gorilla_gorilla",
                           "Mus_musculus_domesticus"))

d$order <- c("Aves", 
"Carnivora",
"Cetartiodactyla",
"Chiroptera",
"Diprotodontia",
"Eulipotyphla",
"Perissodactyla",
"Primates",
"Rodentia")
d$order <- as.factor(d$order)
d$x <- 1:9
d$color <- 'black'
                           
reservoirs <- as.factor(c("AVES",
                          "CARNIVORA",
                          "CETARTIODACTYLA",
                          "CHIROPTERA",
                          "DIPROTODONTIA",
                          "EULIPOTYPHLA" , 
                          "PERISSODACTYLA",
                          #"PILOSA",
                          "PRIMATES",
                          "RODENTIA"))

hOrder <- ggplot(plotData$effects, aes(x = hOrder, y = y, colour = IsSignificant, fill = IsSignificant)) +
  coord_cartesian(ylim = c(-50,50), 
                  clip="off") +
  geom_phylopic(image=d$uid, alpha = 1, colour = 'gray10', x = d$x, y = -60) +
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

#vector-borne
plotData <- get_partial_effects_binary_single(gam_mort, 
                                       c('IsVectorBorne'),
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
  scale_x_discrete(labels = c("vector-borne transmission")) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())

print(binary)

## virus species publication count
plotData_pub_count <- get_partial_effects_continuous(gam_mort, 
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
  labs(x = "virus species publication count", y = 'Effect on CFR in humans') +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 11),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())
print(pub_count)


Fig1a <- cowplot::plot_grid(mort_models, hOrder, vFamily, 
                            ncol=3, nrow=1, align = "h", 
                            labels = c( "A", "B", "C"),
                            label_size = 20)
print(Fig1a)


Fig1b <- cowplot::plot_grid(binary, pub_count, 
                            ncol=1, nrow=2, align = "v",
                            labels = c( "D", "E"),
                            label_size = 20)
print(Fig1b)

Fig1 <- cowplot::plot_grid(Fig1a, Fig1b, 
                           ncol=2, nrow=1, align = "h",
                           rel_widths = c(3,1)) 
print(Fig1)


Fig1 <- ggdraw(Fig1) + 
  draw_plot_label(
    c("AIC", "Deviance\nexplained"),
    c(0.196,0.21),
    c(0.94,0.98),
    size = 8
  )
print(Fig1)

ggsave(file = "Extension/code_and_data/figures/Figure_1.pdf",
       units="mm",
       width=240,
       height=90,
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

#### alternative model forms ####
## with reservoir species richness
gam_mort_sp <- gam(as.formula(RankedModels$Formula[9]),
                data = dat.mort,
                select=FALSE,
                family = gaussian)
summary(gam_mort_sp)
#plot(gam_mort_sp)
#gam.check(gam_mort_sp)
AIC(gam_mort_sp)

## virus family
plotData <- get_partial_effects(gam_mort_sp,
                                var = 'vFamily')

vFamily <- ggplot(plotData$effects, aes(x = vFamily, y = y, colour = IsSignificant, fill = IsSignificant)) +
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

## reservoir species richness
plotData_sp_rich <- get_partial_effects_continuous(gam_mort_sp, 
                                                     var = 'ReservoirNspecies')


sp_rich <- ggplot(plotData_sp_rich$effects, aes(x = ReservoirNspecies, y = y, colour = IsSignificant, fill = IsSignificant)) +
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
  labs(x = "reservoir species richness", y = 'Effect on CFR in humans') +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())
print(sp_rich)

#vector-borne
plotData <- get_partial_effects_binary_single(gam_mort_sp, 
                                              c('IsVectorBorne'),
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
  scale_colour_manual(values = c(No = 'grey60', Yes = 'mediumorchid4'), guide = F) +
  scale_fill_manual(values = c(No = 'grey60', Yes = 'mediumorchid4'), guide = F) +
  scale_y_continuous(labels = function(x) sprintf('%1.1f', x) ) +
  labs(x = NULL, y = 'Effect on CFR in humans') +
  scale_x_discrete(labels = c("vector-borne\ntransmission")) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())

print(binary)

## virus species publication count
plotData_pub_count <- get_partial_effects_continuous(gam_mort_sp, 
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
  labs(x = "virus species publication count", y = 'Effect on CFR in humans') +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 11),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())
print(pub_count)

Fig1a1 <- cowplot::plot_grid(vFamily,
                            ncol=1, nrow=1, align = "h",
                            labels = c( "A"),
                            label_size = 20)
print(Fig1a1)

Fig1a2 <- cowplot::plot_grid(sp_rich,
                            ncol=1, nrow=1, align = "h",
                            labels = c( "B"),
                            label_size = 20)
print(Fig1a2)

Fig1a <- cowplot::plot_grid(Fig1a1, Fig1a2,
                             ncol=2, nrow=1, align = "h",
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

ggsave(file = "Extension/code_and_data/figures/Figure_S2.pdf",
       units="mm",
       width=200,
       height=85,
       scale=2,
       dpi=300)


