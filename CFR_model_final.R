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
dat$phylo_dist <- as.factor(dat$phylo_dist) #also modeled as categorical random effect in PNAS (only 5 unique values)
dat$vFamily <- as.factor(dat$vFamily)
dat$IsVectorBorne <- as.factor(dat$IsVectorBorne)
dat$human.trans <- as.integer(dat$human.trans)

# and viruses with only one human case recorded in history
dat_rare <- dat[dat$tot_cases==1, ]
dat_suf <- dat[dat$tot_cases>1, ]

#### GDP was not significant, so let's instead run the models with global CFR statistics ####

#final set of model variables
dat.mort <- dat_suf %>% distinct(SppName_ICTV_MSL2018b, CFR_avg, hOrder, spill_type, .keep_all = TRUE)
length(unique(dat.mort$SppName_ICTV_MSL2018b))

#how many viruses have multiple entries?
dup <- dat.mort[duplicated(dat.mort$SppName_ICTV_MSL2018b)==TRUE, ] #only 9 viruses
dup_dat <- dat.mort %>%
  filter(SppName_ICTV_MSL2018b %in% dup$SppName_ICTV_MSL2018b)

#global model
global_gamm_mort <- gamm(CFR_avg ~ 
                           s(hOrder, bs="re") +
                           s(vFamily, bs="re") +
                           s(phylo_dist, bs="re") +
                           s(ReservoirPublicationCount, k=7, bs="tp") +
                           s(VirusFamilyPublicationCount, k=7, bs="tp") +
                           s(IsVectorBorne,bs="re") +
                           s(spill_type, bs = "re"),
                         data = dat.mort,
                         family = gaussian)
summary(global_gamm_mort$gam)
gam.check(global_gamm_mort$gam)
plot(global_gamm_mort$gam)

## Model selection
Terms <- list(
  hOrder = c(NA, "s(hOrder, bs='re')"),
  phylo_dist = c(NA, "s(phylo_dist, bs='re')"),
  VirusFamilyPublicationCount  = c(NA, "s(VirusFamilyPublicationCount, k=7, bs='tp')"),
  vFamily  = c(NA, "s(vFamily, bs='re')"),
  ReservoirPublicationCount  = c(NA, "s(ReservoirPublicationCount, k=7, bs='tp')"),
  IsVectorBorne  = c(NA, "s(IsVectorBorne, bs='re')"),
  spill_type  = c(NA, "s(spill_type, bs = 're')")
)

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
  do(ModelFit = try(gamm(as.formula(.$Formula), 
                         family = 'gaussian',
                         data = dat.mort,
                         method = 'REML')))

removeRows <- lapply(CompetingFits$ModelFit, function(x) 'try-error' %in% class(x)) %>% 
  unlist()
FailedFormulas <- CompetingFits[removeRows, ]
stopifnot(nrow(FailedFormulas) == 0)

CompetingFits <- CompetingFits[!removeRows, ]

## Add AIC:
RankedModels <- CompetingFits %>% 
  mutate(AIC = AIC(ModelFit),
         DevianceExplained = summary(ModelFit$gam)$d) %>% 
  ungroup() %>% 
  arrange(AIC) %>% 
  mutate(DeltaAIC = AIC - AIC[1])

## run model with lowest AIC to make plots
lowest_AIC <- RankedModels$Formula[which.min(RankedModels$AIC)]

gamm_mort_lowest <- gamm(as.formula(lowest_AIC),
                         data = dat.mort,
                         family = gaussian)
summary(gamm_mort_lowest$gam)
plot(gamm_mort_lowest$gam)
gam.check(gamm_mort_lowest$gam)

#is there a model within 2 AIC where phylo_dist is significant?
phylo <- RankedModels[str_detect(RankedModels$Formula, "phylo_dist")==TRUE, ]
phylo_fit <- phylo[which.min(phylo$AIC), ]
phylo_fit$DeltaAIC #deltaAIC is 2

gamm_phylo <- gamm(as.formula(phylo_fit$Formula),
                  data = dat.mort,
                  family = gaussian)
summary(gamm_phylo$gam)
#nope

#is there a best fit model with both host order and vector-borne?
models <- RankedModels[str_detect(RankedModels$Formula, "hOrder")==TRUE & str_detect(RankedModels$Formula, "IsVectorBorne")==TRUE, ]
models <- RankedModels[str_detect(RankedModels$Formula, "hOrder")==TRUE, ]

best_fit <- models[which.min(models$AIC), ]
best_fit$DeltaAIC #deltaAIC is only 1.2
best_fit$Formula

#### BEST FIT MODEL ####
gamm_mort <- gamm(as.formula(best_fit$Formula),
                  data = dat.mort,
                  family = gaussian)
summary(gamm_mort$gam)
gam.check(gamm_mort$gam)

#### HOST ORDER PLOT ####
phylo_dist <- data.frame(dat$phylo_dist, dat$hOrder)
colnames(phylo_dist) <- c("phylo_dist", "hOrder")
phylo_dist <- phylo_dist %>% distinct(hOrder, phylo_dist)

order_dat <- plot(gamm_mort$gam, residuals = TRUE, pages =1)[[1]]
orders <- as.factor(levels(order_dat$raw))

order_effects <- cbind.data.frame(orders, order_dat$fit)
names(order_effects) <- c("hOrder", "effect")
print(order_effects) 
order_effects <- merge(order_effects, phylo_dist, by = "hOrder")

order_resid <- cbind.data.frame(order_dat$raw, order_dat$p.resid)
names(order_resid) <- c("hOrder", "residuals")
order_resid <- merge(order_resid, phylo_dist, by = "hOrder")

colors= c('PRIMATES' ="firebrick1", 'PERISSODACTYLA'= "darkorange", 'CETARTIODACTYLA'=  "yellowgreen", 
          'CARNIVORA'="forestgreen",  'SCANDENTIA'="darkseagreen1" , 'LAGOMORPHA'="cornflowerblue", 
          'RODENTIA'=  "gray40",  'CHIROPTERA'="purple", 'DIPROTODONTIA'= "violet", 
          'PROBOSCIDEA'="magenta", 'EULIPOTYPHLA' ="maroon", 'DIDELPHIMORPHIA'="lightpink",
          'AVES'="navyblue")

## try visualizing in cytochrome b distance order
order_effects$hOrder <- factor(order_effects$hOrder, levels=c("PRIMATES",
                                                              "RODENTIA","CETARTIODACTYLA",
                                                              "PERISSODACTYLA","CARNIVORA",
                                                              "CHIROPTERA", "EULIPOTYPHLA",
                                                              #"SCANDENTIA","LAGOMORPHA",
                                                              "DIPROTODONTIA",
                                                              #"PROBOSCIDEA",
                                                              #"DIDELPHIMORPHIA", 
                                                              "AVES"))

mort_hOrder <- ggplot() + 
  geom_point(data=order_effects, aes(x=hOrder, y = effect, col = hOrder), size=10) +
  geom_point(data=order_resid, aes(x=hOrder, y = residuals), col="black", alpha=.3, size=3) +
  geom_hline(aes(yintercept=0), linetype=2, alpha=.5) +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        panel.grid = element_blank(),
        plot.margin = unit(c(.1,.1,0,.1), "cm"),
        #legend.position = "none",
        #legend.position = c(.8,.75),
        legend.key.size = unit(.15,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.background = element_rect(fill="gray97", color="black")) +
  scale_color_manual(values =colors) +
  guides(color=guide_legend(ncol=2)) +
  labs(x="reservoir order",
       y="effect on human case fatality rate")

print(mort_hOrder)

#pdf file
ggsave(file = "Extension/code_and_data/figures/CFR_hOrder.pdf",
       units="mm",
       width=140,
       height=80,
       scale=2,
       dpi=300)

#### VIRUS FAMILY PLOT ####
vfam_dat <- plot(gamm_mort$gam, residuals = TRUE, pages =1)[[3]]
vfams <- as.factor(levels(vfam_dat$raw))

vfam_effects <- cbind.data.frame(vfams, vfam_dat$fit)
names(vfam_effects) <- c("vFamily", "effect")
print(vfam_effects) 

#summarize host orders per virus for pie charts
#to do: 1) more spacing between x tics
#2) manually assign virus family labels on the x axis
#3) make sure hOrder colors match hOrder plot
#4) cowplot to merge hOrder and vFam plots
vir_order <- dat.mort[ ,c("vFamily", "hOrder")]
vir_sum <- dcast(vir_order, vFamily ~hOrder, fun.aggregate = length)
vfam_effects_pie <- merge(vfam_effects, vir_sum, by = "vFamily")

#now factor viral family by order phylogenetic distance
melted <- melt(vfam_effects_pie[ , c(1,3:11)])
top <- melted %>%
  group_by(vFamily) %>%
  top_n(1, value)
top$ID <- 1:nrow(top)
top <- top[top$ID!=15, ]
top <- top[top$ID!=16, ]
top <- top[top$ID!=19, ]
top$variable <- as.character(top$variable)

vfam_effects_pie$max_order <- NA
for(i in 1:nrow(top)){
  vfam_effects_pie$max_order[vfam_effects_pie$vFamily==top$vFamily[i]] <- top$variable[i]
}

host_ranking <- data.frame(hOrder = c("PRIMATES",
                                      "RODENTIA","CETARTIODACTYLA",
                                      "PERISSODACTYLA","CARNIVORA",
                                      "CHIROPTERA", "EULIPOTYPHLA",
                                      "DIPROTODONTIA",
                                      "AVES"),
                           rank = 1:9)
host_ranking$hOrder <- as.character(host_ranking$hOrder)

vfam_effects_pie$max_order_rank <- NA
for(i in 1:nrow(host_ranking)){
  vfam_effects_pie$max_order_rank[vfam_effects_pie$max_order==host_ranking$hOrder[i]] <- host_ranking$rank[i]
}
vfam_effects_pie <- arrange(vfam_effects_pie, max_order_rank)
vfam_effects_pie$axis <- seq(from = 10, to = 180, by = 10) 


vfam_resid <- cbind.data.frame(vfam_dat$raw, vfam_dat$p.resid)
names(vfam_resid) <- c("vFamily", "residuals")
vfam_resid$axis <- NA
for(i in 1:nrow(vfam_effects_pie)){
  vfam_resid$axis[vfam_resid$vFamily==vfam_effects_pie$vFamily[i]] <- vfam_effects_pie$axis[i]
}


colors= c('PRIMATES' ="firebrick1", 'PERISSODACTYLA'= "darkorange", 'CETARTIODACTYLA'=  "yellowgreen", 
          'CARNIVORA'="forestgreen",  'SCANDENTIA'="darkseagreen1" , 'LAGOMORPHA'="cornflowerblue", 
          'RODENTIA'=  "gray40",  'CHIROPTERA'="purple", 'DIPROTODONTIA'= "violet", 
          'PROBOSCIDEA'="magenta", 'EULIPOTYPHLA' ="maroon", 'DIDELPHIMORPHIA'="lightpink",
          'AVES'="navyblue")

mort_vfam <- ggplot() + 
  geom_scatterpie(data = vfam_effects_pie,
                  aes(axis, effect, r=4),
                  color = NA, #gets rid of black outlines in pie chart
                  cols = c("AVES", 
                           "CARNIVORA", 
                           "CETARTIODACTYLA",
                           "CHIROPTERA",
                           "DIPROTODONTIA",
                           "EULIPOTYPHLA",
                           "PERISSODACTYLA",
                           "PRIMATES",
                           "RODENTIA")) +
  coord_fixed() +
  scale_x_continuous(breaks = vfam_effects_pie$axis, labels = vfam_effects_pie$vFamily) +
  geom_point(data=vfam_resid, aes(x=axis, y = residuals), col="black", alpha=.3, size=3) +
  geom_hline(aes(yintercept=0), linetype=2, alpha=.5) +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, size = 15),
        axis.text.y = element_text(size = 12),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        plot.margin = unit(c(.1,.1,0,.1), "cm"),
        legend.position = "none",
        legend.key.size = unit(.15,"cm"),
        legend.title = element_blank(),
        legend.background = element_rect(fill="gray97", color="black")) +
  scale_fill_manual(values = colors) +
  guides(color=guide_legend(ncol=2)) +
  labs(x="virus family",
       y="effect on human case fatality rate")

print(mort_vfam)

#pdf file
ggsave(file = "Extension/code_and_data/figures/CFR_vfam.pdf",
       units="mm",
       width=140,
       height=80,
       scale=2,
       dpi=300)

#### VECTOR ####
vector_dat <- plot(gamm_mort$gam, residuals = TRUE, pages =1)[[5]]
vectors <- as.factor(levels(vector_dat$raw))

vector_effects <- cbind.data.frame(vectors, vector_dat$fit)
names(vector_effects) <- c("vector", "effect")
print(vector_effects) 

#### SPILL TYPE ####
spill_dat <- plot(gamm_mort$gam, residuals = TRUE, pages =1)[[6]]
spills <- as.factor(levels(spill_dat$raw))

spill_effects <- cbind.data.frame(spills, spill_dat$fit)
names(spill_effects) <- c("spill", "effect")
print(spill_effects) 




