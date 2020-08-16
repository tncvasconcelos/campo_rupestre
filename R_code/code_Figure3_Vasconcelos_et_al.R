#################################################################
############## Code to reproduce plots in Figure 2 ##############
#                     (among other things)
#################################################################

# This is part of:
# Vasconcelos TNC, Alcantara S, Andrino CO, Forest F, Reginato M, Simon MF, Pirani JR 
# Fast diversification through a mosaic of evolutionary histories 
# characterizes the endemic flora of ancient Neotropical mountains. Proceedings B.


# rm(list=ls())

library(ape)
library(hisse)
library(parallel)
library(tidyverse)
library(phytools)
library(PMCMR)
library(dplyr)

######## RATES #########

setwd("ESM_data/GeoHiSSE/")
# set working directory to the GeoHiSSE folder in the ESM data

list.files(pattern="models") -> models
list.files(pattern="recon") -> recon

sub("_GeoHiSSE_models.RData", "", models) -> group

sapply(models, readRDS, simplify=F) -> list.geohisse
sapply(recon, readRDS, simplify=F) -> recon.models

#### This will extract speciation, extinction and net-diversification rates for each group
# based on the 3 best models among the 36 GeoHiSSE models;
# it will also plot boxplots for each group
# and save tables that will be used to plot Figure 2

for(i in 1:length(group)){
  ## extract models and recons
  l0<-list() 
  l1<-list() 
  list.geohisse[[i]]->l0
  recon.models[[i]]->l1
  assign(paste(group[i], "models", sep="_"), l0) 
  assign(paste(group[i], "recon", sep="_"), l1)
  
  if(group[i]=="Chamaecrista"){
    # excluding model 7 # - this model had to be excluded in Chamaecrista 
    # due to extreme values of AIC
    
  ## getting 3 best models - Akaike weights ##
    aic<-GetAICWeights(l0[-7], criterion="AIC")
    round(tail(sort(aic),3),3)->x
    which(aic==x[1])->m1
    which(aic==x[2])->m2
    which(aic==x[3])->m3
    
  }
  else {
  ## getting 3 best models - Akaike weights ##

  aic<-GetAICWeights(l0, criterion="AIC")
  tail(sort(aic),3)->x
  which(aic==x[1])->m1
  which(aic==x[2])->m2
  which(aic==x[3])->m3
}
  ## reconstruct areas and rates ##
  
  names<-c("speciation","extinction","net.div")
  
  for(u in 1:length(names)){

    xx <- plot.geohisse.states(x = l1[c(m1,m2,m3)], 
                               rate.param = names[u], type = "phylogram",  
                               show.tip.label = F, fsize=0.2 ,
                               legend = T,  state.colors=c("blue","red","purple")) 
    
    xx[1]-> rate_tree
    xx[2]-> area_tree
  
    ## boxplots for each group ##
    
    rate_tree$rate.tree->obj
    tree=rate_tree$rate.tree$tree
    
    ii<-c(as.numeric(names(obj$tree$maps[[1]])[1]),
          sapply(obj$tree$maps,function(x) as.numeric(names(x)[length(x)])))
    
    a<-setNames(ii/(length(obj$cols)-1)*diff(obj$lims)+obj$lims[1],
                c(obj$tree$edge[1,1],obj$tree$edge[,2]))
    
    l0$model1->mod1
    dat=mod1$data
    a<-a[as.character(1:Ntip(mod1$phy))]
    as.data.frame(a)->a
    rownames(a)<-mod1$phy$tip.label
    rownames(dat)<-dat[,1]
    
    merge(a,dat,by="row.names",all=TRUE)->states_rates
    for (z in 1:length(states_rates$area)){
      if (states_rates[z, "area"] == 0 ) {
        states_rates[z, "area"] = "Both (W)" 
      }
      if (states_rates[z, "area"] == 1 ) {
        states_rates[z, "area"] = "Non_cr" 
      }
      if (states_rates[z, "area"] == 2 ) {
        states_rates[z, "area"] = "Cr" 
      }
    }
    
    pdf(file=paste(names[u], group[i], "boxplot.pdf", sep="_"), height = 8, width=4)
    boxplot(states_rates$a~states_rates$area, notch=F, pch=20, col=c("darkorchid","tomato2","dodgerblue2"))
    dev.off()
    
    write.csv(states_rates, file=paste(names[u],"rates", group[i],".csv", sep="_"))
  }
}  


######################################################################
######################################################################

# Figure 2 plots and conover tests #


pdf(file="Figure2_base.pdf")
par(mfrow=c(3,3))

#### ALL GROUPS ####
list.files(pattern="speciation_rates") -> speciation
list.files(pattern="extinction_rates") -> extinction
list.files(pattern="net.div_rates") -> net.div

sapply(speciation, read.csv, simplify=F) -> speciation_list
sapply(extinction, read.csv, simplify=F) -> extinction_list
sapply(net.div, read.csv, simplify=F) -> net.div_list

as.data.frame(data.table::rbindlist(speciation_list, fill=T))->speciation_all
speciation_all %>% filter(a != 0, a != "NA") -> speciation_all
as.data.frame(data.table::rbindlist(extinction_list, fill=T))->extinction_all
extinction_all %>% filter(a != 0, a != "NA") -> extinction_all
as.data.frame(data.table::rbindlist(net.div_list, fill=T))->net.div_all
net.div_all %>% filter(a != 0, a != "NA") -> net.div_all

boxplot(log(speciation_all$a)~speciation_all$area, notch=T, ylim=c(-5, 3))
posthoc.kruskal.conover.test(speciation_all$a,speciation_all$area)
boxplot(log(extinction_all$a)~extinction_all$area, notch=T, ylim=c(-10, 3))
posthoc.kruskal.conover.test(extinction_all$a,extinction_all$area)
boxplot(log(net.div_all$a)~net.div_all$area, notch=T, ylim=c(-5, 4))
posthoc.kruskal.conover.test(net.div_all$a,net.div_all$area)


#### EUDICOTS ####

eudicots<-c("Calliandra","Chamaecrista","Lychnophorinae","Metastelmatinae","Myrcia","Mimosa","Marcetieae", "Melastomateae","Diplusodon")

eud_spe<-c()
eud_ext<-c()
eud_net<-c()
for (h in 1:length(eudicots)){
  grep(eudicots[h], speciation)->x0
  eud_spe<-c(speciation[x0], eud_spe)
  eud_ext<-c(extinction[x0], eud_ext)
  eud_net<-c(net.div[x0], eud_net)
}

sapply(eud_spe, read.csv, simplify=F) -> speciation_list
sapply(eud_ext, read.csv, simplify=F) -> extinction_list
sapply(eud_net, read.csv, simplify=F) -> net.div_list

as.data.frame(data.table::rbindlist(speciation_list, fill=T))->speciation_all
speciation_all %>% filter(a != 0, a != "NA") -> speciation_all
as.data.frame(data.table::rbindlist(extinction_list, fill=T))->extinction_all
extinction_all %>% filter(a != 0, a != "NA") -> extinction_all
as.data.frame(data.table::rbindlist(net.div_list, fill=T))->net.div_all
net.div_all %>% filter(a != 0, a != "NA") -> net.div_all

boxplot(log(speciation_all$a)~speciation_all$area, notch=T, ylim=c(-5, 3))
posthoc.kruskal.conover.test(speciation_all$a,speciation_all$area)
boxplot(log(extinction_all$a)~extinction_all$area, notch=T, ylim=c(-10, 3))
posthoc.kruskal.conover.test(extinction_all$a,extinction_all$area)
boxplot(log(net.div_all$a)~net.div_all$area, notch=T, ylim=c(-5, 4))
posthoc.kruskal.conover.test(net.div_all$a,net.div_all$area)


#### MONOCOTS ####
monocots<-c("Cattleya","Habenaria","Dyckia","Paepalanthus","Velloziaceae","Trimezieae")

mon_spe<-c()
mon_ext<-c()
mon_net<-c()
for (h in 1:length(monocots)){
  grep(monocots[h], speciation)->x0
  mon_spe<-c(speciation[x0], mon_spe)
  mon_ext<-c(extinction[x0], mon_ext)
  mon_net<-c(net.div[x0], mon_net)
}

sapply(mon_spe, read.csv, simplify=F) -> speciation_list
sapply(mon_ext, read.csv, simplify=F) -> extinction_list
sapply(mon_net, read.csv, simplify=F) -> net.div_list

as.data.frame(data.table::rbindlist(speciation_list, fill=T))->speciation_all
speciation_all %>% filter(a != 0, a != "NA") -> speciation_all
as.data.frame(data.table::rbindlist(extinction_list, fill=T))->extinction_all
extinction_all %>% filter(a != 0, a != "NA") -> extinction_all
as.data.frame(data.table::rbindlist(net.div_list, fill=T))->net.div_all
net.div_all %>% filter(a != 0, a != "NA") -> net.div_all

boxplot(log(speciation_all$a)~speciation_all$area, notch=T, ylim=c(-5, 4))
posthoc.kruskal.conover.test(speciation_all$a,speciation_all$area)
boxplot(log(extinction_all$a)~extinction_all$area, notch=T, ylim=c(-10, 3))
posthoc.kruskal.conover.test(extinction_all$a,extinction_all$area)
boxplot(log(net.div_all$a)~net.div_all$area, notch=T, ylim=c(-5, 4))
posthoc.kruskal.conover.test(net.div_all$a,net.div_all$area)

###############################################################

dev.off()
