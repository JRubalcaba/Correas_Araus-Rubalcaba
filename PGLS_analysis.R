require(ape)
require(MCMCglmm)

#### Load FMR database ####

dir <- "..."
data <- read.csv(paste0(dir,"/FMR_database.csv"), sep=";")

#### Phylogenetic analyses ####

phylo <- read.tree(paste0(dir,"/FritzTree_mammals_consensus.tre")) # Consensus tree for mammals (Fritz et al. 2009)

data_mammals_pgls <- data.frame(Species=data$Species[which(data$Class=="Mammalia")], 
                                 FMR_Watt=data$FMR_Watt[which(data$Class=="Mammalia")], 
                                 FMR_M=data$FMR_M[which(data$Class=="Mammalia")], 
                                 FMR_M_ind_simul=FMR_M_ind_simul, 
                                 FMR_ind_simul=FMR_ind_simul,
                                 Mass=data$Mass[which(data$Class=="Mammalia")],
                                 LocationCode=data$LocationCode[which(data$Class=="Mammalia")], 
                                 TALOCmean=data$TALOCmean[which(data$Class=="Mammalia")],
                                 TMEAN=data$TMEAN[which(data$Class=="Mammalia")], 
                                 Tanomalies=data$Tanomalies[which(data$Class=="Mammalia")])
data_mammals_pgls <- data_mammals_pgls[complete.cases(data_mammals_pgls),] # Remove NAs

data_mammals_pgls$animal <- "not in tree"
for(i in 1:nrow(data_mammals_pgls)){
  species <- data_mammals_pgls$Species[i]
  if(rlang::is_empty(phylo$tip.label[grepl(species,phylo$tip.label)])){ # if sp is not in tree leave "not in tree"
  } else {
    data_mammals_pgls$animal[i]<-phylo$tip.label[grepl(species,phylo$tip.label)]} # else put the sp name from the tree
}

length(unique(data_mammals_pgls$animal[which(!data_mammals_pgls$animal=="not in tree")])) # species in tree
length(unique(data_mammals_pgls$Species[which(data_mammals_pgls$animal=="not in tree")])) # not in tree

nameslist <- phylo$tip.label
treenameslist <- as.data.frame(table(data_mammals_pgls$animal))
Speciestoretain <- intersect(treenameslist$Var1, nameslist)
pruned.tree <- drop.tip(phylo,phylo$tip.label[-match(Speciestoretain,phylo$tip.label)])
plot(pruned.tree, cex=0.5)
tree <- pruned.tree
tree <- compute.brlen(tree)
tree$node.label<-NULL
is.ultrametric(tree) # checks
is.rooted(tree)
any(duplicated(tree$node.label))

# Reordenamos los data para que coincidan especies en la matriz de data y en el árbol
data_mammals_pgls_sub <- data_mammals_pgls[(data_mammals_pgls$animal %in% tree$tip.label),]
data_mammals_pgls_sub <- data_mammals_pgls[which(!data_mammals_pgls$animal=="not in tree"),]
length(unique(data_mammals_pgls_sub$animal)) #  40 sp
nrow(data_mammals_pgls_sub)

########## Mean T vs T anomaly -----
# Phylogenetic mixed-effect models 

mod_mamm_T_Tanom_inter <- MCMCglmm(FMR_M ~ TMEAN * Tanomalies + I(TMEAN^2), 
                                          random=~Species, data=data_mammals_pgls_sub, 
                                          pedigree=tree, nitt=10000)
mod_mamm_T_Tanom <- MCMCglmm(FMR_M ~ TMEAN + I(TMEAN^2) +Tanomalies, 
                                   random=~Species, data=data_mammals_pgls_sub, 
                                   pedigree=tree, nitt=10000)
mod_mamm_T_Tanom_Tanom2 <- MCMCglmm(FMR_M ~ TMEAN + I(TMEAN^2) + Tanomalies + I(Tanomalies^2), 
                                          random=~Species, data=data_mammals_pgls_sub, 
                                          pedigree=tree, nitt=10000)
mod_mamm_T_Tanom_Tanom2_inter <- MCMCglmm(FMR_M ~ TMEAN * Tanomalies + I(Tanomalies^2) + I(TMEAN^2), 
                                          random=~Species, data=data_mammals_pgls_sub, 
                                          pedigree=tree, nitt=10000)
mod_mamm_Tanom <- MCMCglmm(FMR_M ~ Tanomalies, 
                       random=~Species, data=data_mammals_pgls_sub, 
                       pedigree=tree, nitt=10000)
mod_mamm_Tanom2 <- MCMCglmm(FMR_M ~ Tanomalies+I(Tanomalies^2), 
                           random=~Species, data=data_mammals_pgls_sub, 
                           pedigree=tree, nitt=10000)
mod_mamm_T <- MCMCglmm(FMR_M ~ TMEAN, 
                            random=~Species, data=data_mammals_pgls_sub, 
                            pedigree=tree, nitt=10000)
mod_mamm_T2 <- MCMCglmm(FMR_M ~ TMEAN + I(TMEAN^2), 
                       random=~Species, data=data_mammals_pgls_sub, 
                       pedigree=tree, nitt=10000)

# DIC comparison
mod_mamm_T$DIC
mod_mamm_T2$DIC
mod_mamm_Tanom$DIC
mod_mamm_Tanom2$DIC
mod_mamm_T_Tanom_inter$DIC
mod_mamm_T_Tanom$DIC
mod_mamm_T_Tanom_Tanom2_inter$DIC # Best model (lowest DIC)
mod_mamm_T_Tanom_Tanom2$DIC

###### Results -------

require(ggplot2)
require(RColorBrewer)

# Check model MCMC
plot(mod_mamm_T_Tanom_Tanom2_inter)
summary(mod_mamm_T_Tanom_Tanom2_inter) 

# Parameters of the best model
best_model <- function(Tmean, Tanom){
  0.6360113-0.0671558*Tmean+0.1096303*Tanom+0.0181288*Tanom^2+0.0011452*Tmean^2-0.0100275*Tanom*Tmean
}

summary(data_mammals_pgls_sub$TMEAN) # Observed ranges of mean T
Tmean1 = 10 # Mean Temperature for data below median T
Tmean3 = 23 # Mean T for data above median T

# FMR vs mean T
ggplot(data_mammals_pgls_sub, aes(y=FMR_M, x=TMEAN)) + theme_classic() +
  ylab("Mass-specific Field Metabolic Rate") + xlab("Mean Temperature (ºC)") +
  geom_point(size=2) + 
  scale_color_gradient(low="black", high="grey80") +
  geom_function(fun=function(x) best_model(Tanom=0,Tmean=x), lwd=1, col="black") 
ggplot() + theme_classic() +
  ylab("Mass-specific Field Metabolic Rate") + xlab("Temperature anomaly (ºC)") +
  geom_point(data_mammals_pgls_sub[which(data_mammals_pgls_sub$TMEAN<19),], mapping=aes(y=FMR_M, x=Tanomalies), size=2, col="black") + 
  geom_point(data_mammals_pgls_sub[which(data_mammals_pgls_sub$TMEAN>19),], mapping=aes(y=FMR_M, x=Tanomalies), size=2, col="grey80") + 
  scale_color_gradient(low="black", high="grey80") +
  geom_function(fun=function(x) best_model(Tanom=x,Tmean=Tmean1), lwd=1, col="black") + 
  geom_function(fun=function(x) best_model(Tanom=x,Tmean=Tmean3), lwd=1, col="grey80")

#### Obseved vs Predicted FMR -----

data_mammals_pgls_warm <- data_mammals_pgls_sub[data_mammals_pgls_sub$TMEAN>19,]
data_mammals_pgls_cold <- data_mammals_pgls_sub[data_mammals_pgls_sub$TMEAN<19,]

data_mammals_pgls_warm$animal <- "not in tree"
for(i in 1:nrow(data_mammals_pgls_warm)){
  species <- data_mammals_pgls_warm$Species[i]
  if(rlang::is_empty(phylo$tip.label[grepl(species,phylo$tip.label)])){ # if sp is not in tree leave "not in tree"
  } else {
    data_mammals_pgls_warm$animal[i]<-phylo$tip.label[grepl(species,phylo$tip.label)]} # else put the sp name from the tree
}

length(unique(data_mammals_pgls_warm$animal[which(!data_mammals_pgls_warm$animal=="not in tree")])) # species in tree
length(unique(data_mammals_pgls_warm$Species[which(data_mammals_pgls_warm$animal=="not in tree")])) # not in tree

data_mammals_pgls_cold$animal <- "not in tree"
for(i in 1:nrow(data_mammals_pgls_cold)){
  species <- data_mammals_pgls_cold$Species[i]
  if(rlang::is_empty(phylo$tip.label[grepl(species,phylo$tip.label)])){ # if sp is not in tree leave "not in tree"
  } else {
    data_mammals_pgls_cold$animal[i]<-phylo$tip.label[grepl(species,phylo$tip.label)]} # else put the sp name from the tree
}

length(unique(data_mammals_pgls_cold$animal[which(!data_mammals_pgls_cold$animal=="not in tree")])) # species in tree
length(unique(data_mammals_pgls_cold$Species[which(data_mammals_pgls_cold$animal=="not in tree")])) # not in tree

nameslist <- phylo$tip.label
treenameslist <- as.data.frame(table(data_mammals_pgls$animal))
Speciestoretain <- intersect(treenameslist$Var1, nameslist)
pruned.tree <- drop.tip(phylo,phylo$tip.label[-match(Speciestoretain,phylo$tip.label)])
plot(pruned.tree, cex=0.5)
tree_warm <- pruned.tree
tree_warm <- compute.brlen(tree_warm)
tree_warm$node.label<-NULL
is.ultrametric(tree_warm) # checks
is.rooted(tree_warm)
any(duplicated(tree_warm$node.label))

nameslist <- phylo$tip.label
treenameslist <- as.data.frame(table(data_mammals_pgls$animal))
Speciestoretain <- intersect(treenameslist$Var1, nameslist)
pruned.tree <- drop.tip(phylo,phylo$tip.label[-match(Speciestoretain,phylo$tip.label)])
plot(pruned.tree, cex=0.5)
tree_cold <- pruned.tree
tree_cold <- compute.brlen(tree_cold)
tree_cold$node.label<-NULL
is.ultrametric(tree_cold) # checks
is.rooted(tree_cold)
any(duplicated(tree_cold$node.label))

data_mammals_pgls_warm_sub <- data_mammals_pgls_warm[(data_mammals_pgls_warm$animal %in% tree_warm$tip.label),]
data_mammals_pgls_warm_sub <- data_mammals_pgls_warm_sub[which(!data_mammals_pgls_warm_sub$animal=="not in tree"),]
length(unique(data_mammals_pgls_warm_sub$animal)) #  19 sp

data_mammals_pgls_cold_sub <- data_mammals_pgls_cold[(data_mammals_pgls_cold$animal %in% tree_cold$tip.label),]
data_mammals_pgls_cold_sub <- data_mammals_pgls_cold_sub[which(!data_mammals_pgls_cold_sub$animal=="not in tree"),]
length(unique(data_mammals_pgls_cold_sub$animal)) #  24 sp

###### Results -------

OBSvsPRED <- MCMCglmm(log(FMR_Watt) ~ log(FMR_ind_simul) * TMEAN, 
                      random=~Species, data=data_mammals_pgls_sub, 
                      pedigree=tree, nitt=10000)
summary(OBSvsPRED)  

OBSvsPRED_warm <- MCMCglmm(log(FMR_Watt) ~ log(FMR_ind_simul), 
                      random=~Species, data=data_mammals_pgls_warm_sub, 
                      pedigree=tree, nitt=10000)
summary(OBSvsPRED_warm)  

OBSvsPRED_cold <- MCMCglmm(log(FMR_Watt) ~ log(FMR_ind_simul), 
                           random=~Species, data=data_mammals_pgls_cold_sub, 
                           pedigree=tree, nitt=10000)
summary(OBSvsPRED_cold) 

## Mass-specific FMR vs Temperature for observed and predicted

FMR_TempOBS <- MCMCglmm(FMR_M ~ TALOCmean, 
                        random=~Species, data=data_mammals_pgls_sub, 
                        pedigree=tree, nitt=10000)
FMR_TempPRED <- MCMCglmm(FMR_M ~ TALOCmean, 
                        random=~Species, data=data_mammals_pgls_sub, 
                        pedigree=tree, nitt=10000)
summary(FMR_TempOBS) 
summary(FMR_TempPRED) 

FMR_TempOBS_warm <- MCMCglmm(FMR_M ~ TALOCmean, 
                      random=~Species, data=data_mammals_pgls_warm_sub, 
                      pedigree=tree, nitt=10000)
FMR_TempPRED_warm <- MCMCglmm(FMR_M_ind_simul ~ TALOCmean, 
                     random=~Species, data=data_mammals_pgls_warm_sub, 
                     pedigree=tree, nitt=10000)
summary(FMR_TempOBS_warm)
summary(FMR_TempPRED_warm)  

FMR_TempOBS_cold <- MCMCglmm(FMR_M ~ TALOCmean, 
                             random=~Species, data=data_mammals_pgls_cold_sub, 
                             pedigree=tree, nitt=10000)
FMR_TempPRED_cold <- MCMCglmm(FMR_M_ind_simul ~ TALOCmean, 
                              random=~Species, data=data_mammals_pgls_cold_sub, 
                              pedigree=tree, nitt=10000)
summary(FMR_TempOBS_cold)
summary(FMR_TempPRED_cold) 

## Plots

p1 <- ggplot() + 
  theme_classic() + theme(axis.text = element_text(colour="black")) +
  ylab("log Field Metabolic Rate (W)") + xlab("log Predicted Metabolic Rate (W)") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12, colour="black")) +
  geom_point(data_mammals_pgls_cold_sub, mapping=aes(y=log(FMR_Watt), x=log(FMR_ind_simul))) + 
  geom_point(data_mammals_pgls_warm_sub, mapping=aes(y=log(FMR_Watt), x=log(FMR_ind_simul)), col="grey80") + 
  
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  geom_abline(intercept = 0.9950, slope = 0.7995, colour="grey80", size=1) +
  geom_abline(intercept = 0.6137, slope = 0.9648, colour="black", size=1)

p2 <- ggplot() + theme_classic() +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12, colour="black")) +
  ylab("Modelled mass-specific Field Metabolic Rate") + xlab("Temperature (ºC)") +
  theme(axis.text = element_text(colour="black")) +
  geom_point(data_mammals_pgls_cold_sub, mapping=aes(y=FMR_M, x=TALOCmean), size=2,  pch=1, col="black") +
  geom_abline(intercept=0.53982, slope=-0.03145, col="black", size=1, linetype="dashed")+

  geom_point(data_mammals_pgls_cold_sub, mapping=aes(y=FMR_M_ind_simul, x=TALOCmean), size=2, col="black") + 
  geom_abline(intercept=1.12911, slope=-0.06963, col="black", size=1, linetype="solid")

p3 <- ggplot() + theme_classic() +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12, colour="black")) +
  ylab("Modelled mass-specific Field Metabolic Rate") + xlab("Temperature (ºC)") +
  geom_point(data_mammals_pgls_warm_sub, mapping=aes(y=FMR_M, x=TALOCmean), size=2, pch=1, col="black") +
  geom_abline(intercept=-0.047444, slope=-0.007519, col="black", size=1, linetype="dashed")+
  
  geom_point(data_mammals_pgls_warm_sub, mapping=aes(y=FMR_M_ind_simul, x=TALOCmean), size=2,  col="black") + 
  geom_abline(intercept=5.4283, slope=-0.2298, col="black", size=1, linetype="solid")

require(ggpubr)
ggarrange(p1,p2,p3,ncol=3) #1000 x 300
