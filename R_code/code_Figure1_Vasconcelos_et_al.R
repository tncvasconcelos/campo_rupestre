#################################################################
############## Code to reproduce trees in Figure 1 ##############
#                     (among other things)
#################################################################

# This is part of:
# Vasconcelos TNC, Alcantara S, Andrino CO, Forest F, Reginato M, Simon MF, Pirani JR 
# Fast diversification through a mosaic of evolutionary histories 
# characterizes the endemic flora of ancient Neotropical mountains. Proceedings B.

# rm(list=ls()) # uncomment to clean your environment 

library(ape)
library(hisse)
library(parallel)
library(tidyverse)
library(phytools)
library(dplyr)
library(geiger)
library(magrittr)
library(BioGeoBEARS)
library(optimx)
library(phangorn)

###############################################################
###################### GEOSSE CODE ############################
###############################################################

setwd("/ESM_data/GeoHiSSE/")
# set working directory to the GeoHiSSE folder in the ESM data

list.files(pattern="models") -> models
list.files(pattern="recon") -> recon

sub("_GeoHiSSE_models.RData", "", models) -> group

sapply(models, readRDS, simplify=F) -> list.geohisse
sapply(recon, readRDS, simplify=F) -> recon.models

trees_geohisse <- list() # store trees 
for(i in 1:length(group)){
  ## extract models and recons files

  l0<-list() 
  l1<-list() 
  list.geohisse[[i]]->l0
  recon.models[[i]]->l1
  assign(paste(group[i], "models", sep="_"), l0) 
  assign(paste(group[i], "recon", sep="_"), l1)
  trees_geohisse[[i]] <- l0$model1$phy
  
  if(group[i]=="Chamaecrista"){
    # excluding model 7 # - this model had to be excluded in Chamaecrista 
    # due to extreme values of AIC
  
  # aic plot
  pdf(file=paste("akaikeweight",group[i],".pdf", sep="_"), width=10, height=5)
  lab=1:36
  plot(GetAICWeights(l0[-7], criterion="AIC"), pch=20, ylab="akaike weight", xlab=lab[-7])
  
  dev.off()
  
  # contmap
  aic<-GetAICWeights(l0[-7], criterion="AIC")
  tail(sort(aic),3)->x
  which(aic==x[1])->m1
  which(aic==x[2])->m2
  which(aic==x[3])->m3

  xx <- plot.geohisse.states(x = l1[c(m1,m2,m3)], 
                             rate.param = "speciation", type = "phylogram",  
                             show.tip.label = F, fsize=0.2 ,
                             legend = T,  state.colors=c("blue","red","purple")) 
  xx[2]-> area_tree
  
  
  ## range reconstruction - using original GeoSSE full model (model number 2) 
  l1[[2]]->rec1
  rec1<-rec1$node.mat[,2:4]
  colnames(rec1)<-c("V1","V2","V3")
  marginal_final<-rec1
  write.csv(marginal_final, file=paste(group[i], "marginal_prob_geohisse.csv", sep="_"))
  
  ##
  
}
  else {
    # aic plot
    pdf(file=paste("akaikeweight",group[i],".pdf", sep="_"), width=10, height=5)
    plot(GetAICWeights(l0, criterion="AIC"), pch=20, ylab="akaike weight")
    dev.off()
    
    # contmap
    aic<-GetAICWeights(l0, criterion="AIC")
    tail(sort(aic),3)->x
    which(aic==x[1])->m1
    which(aic==x[2])->m2
    which(aic==x[3])->m3

    
    xx <- plot.geohisse.states(x = l1[c(m1,m2,m3)], 
                               rate.param = "speciation", type = "phylogram",  
                               show.tip.label = F, fsize=0.2 ,
                               legend = T,  state.colors=c("blue","red","purple")) 
    xx[2]-> area_tree
    
  
    ## range reconstruction - using original GeoSSE full model (model number 2) 
    l1[[2]]->rec1
    rec1<-rec1$node.mat[,2:4]
    colnames(rec1)<-c("V1","V2","V3")
    marginal_final<-rec1
    write.csv(marginal_final, file=paste(group[i], "marginal_prob_geohisse.csv", sep="_"))
    marginal_final<-as.matrix(marginal_final)
  }
  
  
## PLOTS   
  colors_list_for_states<-c("blue", "red", "purple")
  
  pdf(file=paste("rec_trees", group[i],"tree.pdf", sep="_"))
  plot(area_tree$state.tree, fsize=0.2, type="phylogram", lwd=0.75, outline=F)
  title(main = "Range reconstruction (contmap) all models")
  
  plot(l0$model1$phy, show.tip.label=F)  
  nodelabels(pie=marginal_final[1:length(marginal_final[,1]),],cex=0.5, 
             piecol=colors_list_for_states)
  title(main = "Range reconstruction (marginal probabilities) GeoSSE only")
  
  dev.off()
  
  
  #### Extracting age of in situ speciation events ####
  
  Ntip(l0$model1$phy) -> n.tips
  Nnode(l0$model1$phy) -> n.nodes
  c((n.tips++1):(n.tips++n.nodes)) -> nodes
  
  max(nodeHeights(l0$model1$phy)) -> tree.length
  sapply(nodes, nodeheight, tree=l0$model1$phy) -> n
  (n-tree.length)*-1 -> ages
  
  #
  
  cbind(marginal_final, ages)->ranges1
  insitu1 <- as.data.frame(ranges1) %>% filter(V2 > V1 + V3) # thresholds for insitu
  insitu1[,c("V1","V3")]<-NULL
  rep(group[i], length(insitu1$V2))->n1
  cbind(n1, insitu1)->insitu1
  
  write.csv(insitu1, file=paste(group[i], "nodes_insitu_geosse.csv", sep="_"))
  
  #

}

list.files(pattern="prob_geohisse.csv") -> prob
sapply(prob, read.csv, simplify=F) -> probs


###############################################################
################### BIOGEOBEARS CODE ##########################
###############################################################


setwd("ESM_data/BioGeoBEARS")
# set working directory to the BioGeoBEARS folder in the ESM data

list.files(pattern="dist.csv") -> files.mat
sub("dist.csv", "clean.tre", files.mat) -> files.mcc
sapply(files.mcc, read.tree, simplify=F) -> trees.mcc
sub("_clean.tre", "", files.mcc) -> labs

keep = c("cr", "nonCr")
areanames = c("E", "N")

sapply(files.mat, read.csv, row.names=1, simplify=F) -> mats
lapply(mats, subset, select=keep) -> mats

### MCC #####
for (i in 1:length(labs)) {
  labs[i] -> l0
  mats[[i]] -> m0
  colnames(m0) <- areanames
  trees.mcc[[i]] -> mcc0
  rownames(m0)[which(m0$E == "OUT")] -> drop
  if (length(drop) > 0) {
    drop.tip(mcc0, drop) -> mcc0
    m0[-which(m0$E == "OUT"),] -> m0
    mcc0 -> trees.mcc[[i]]
  }
  m0[mcc0$tip.label,] -> m0
  m0 -> mats[[i]]
  define_tipranges_object(m0) -> m0
  save_tipranges_to_LagrangePHYLIP(m0, lgdata_fn = paste(l0, "_phylip.txt", sep=""))
}

#

### setup geral 
max_range_size = 2
n.cores = 4
c("Camporupestre","nonCamporupestre") -> actual_names

#############################################
### Compare DEC and DEC+j
#############################################

vector("list", length=length(labs)) -> best.models -> model.comps -> ml.models
names(ml.models) <- names(best.models) <- names(model.comps) <- labs

results<-list()
for (i in 1:length(labs)) {
  labs[i] -> l0
  trees.mcc[[i]] -> t0
  
  if (is.ultrametric(t0) == F) {
    force.ultrametric(t0) -> t0
    t0$edge.length[t0$edge.length == 0] <- 0.01
    force.ultrametric(t0) -> t0
  }
  if (min(t0$edge.length)==0) {
    force.ultrametric(t0) -> t0
    t0$edge.length[t0$edge.length == 0] <- 0.01
    force.ultrametric(t0) -> t0
  }
  
  tree.file = paste(l0, "_biogeo.tre", sep="")
  matrix.file = paste(l0, "_phylip.txt", sep="")
  write.tree(t0, tree.file)
  
  ### DEC
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = tree.file
  BioGeoBEARS_run_object$geogfn = matrix.file
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$speedup = F          
  BioGeoBEARS_run_object$use_optimx = TRUE     
  BioGeoBEARS_run_object$num_cores_to_use = n.cores
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE    
  resDEC = bears_optim_run(BioGeoBEARS_run_object)
  ### DEC+J
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = tree.file
  BioGeoBEARS_run_object$geogfn = matrix.file
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$speedup = F          
  BioGeoBEARS_run_object$use_optimx = TRUE     
  BioGeoBEARS_run_object$num_cores_to_use = n.cores
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE    
  dstart = resDEC$outputs@params_table["d","est"]
  estart = resDEC$outputs@params_table["e","est"]
  jstart = 0.0001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
  resDECj = bears_optim_run(BioGeoBEARS_run_object)
  ### Stats
  restable = NULL
  teststable = NULL
  LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
  LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECj)
  numparams1 = 3
  numparams2 = 2
  stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
  res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  tmp_tests = conditional_format_table(stats)
  restable = rbind(restable, res2, res1)
  teststable = rbind(teststable, tmp_tests)
  row.names(restable) = c("DEC", "DEC+J")
  teststable$alt = "DEC+J"
  teststable$null = "DEC"
  restable$AIC <- c(teststable$AIC2, teststable$AIC1)
  restable -> model.comps[[i]]
  write.table(as.matrix(restable), file=paste(l0, "_biogeobears_pars_AIC.txt", sep=""))
  rownames(restable)[which.min(restable$AIC)] -> best.model
  if (best.model == "DEC") {
    resDEC -> biogeo
  } else {
    resDECj -> biogeo
  }
  
  
  biogeo -> ml.models[[i]]
  best.model -> best.models[[i]]
  

  results[[i]]<-biogeo$ML_marginal_prob_each_state_at_branch_top_AT_node
  
  ### plot 
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  pdf(file=paste(l0, "_biogeobears.pdf", sep=""), width=5, height=10)

  colors_list_for_states<-c("green", "red", "blue", "darkorchid")
  plot_BioGeoBEARS_results(biogeo, paste(l0, "-", best.model), addl_params=list("j"), plotwhat="pie", colors_list_for_states=colors_list_for_states, label.offset=0.001, tipcex=0.001, statecex=1.1, splitcex=0.4, titlecex=0.8, plotsplits=F,  cornercoords_loc=scriptdir, include_null_range=T, tr=t0, show.tip.label = F, tipranges = NULL)
  
  dev.off()
  try(dev.off())
} 

saveRDS(results, file="BioGeoBears_marginal_probs.RData")

###############################################################
####################### PIE CHARTS ############################
###############################################################



setwd("/")

#### BioGeoBEARS ####

pdf(file="Figure1_base.pdf", width=30, height=200)

# layout 
bla<-c()
for(z in 1:length(trees.mcc)){
  length(trees.mcc[[z]]$tip.label)/40->b
  bla<-c(bla, b)
}

layout(matrix(1:30,ncol=2), heights=bla, TRUE) 

# plot

for(i in 1:length(trees.mcc)){
  
  plot(trees.mcc[[i]], show.tip.label=F, x.lim=c(1,80))  
  x<-results[[i]]
  npieMax<-trees.mcc[[i]]$Nnode+length(trees.mcc[[i]]$tip.label)
  npieMin<-length(trees.mcc[[i]]$tip.label)+1
  colors_list_for_states<-c("green", "red", "blue", "darkorchid")
  nodelabels(pie=x[npieMin:npieMax,],cex=0.4, piecol=colors_list_for_states)
  axisPhylo()
  
}

for(i in 1:length(trees_geohisse)){
  
  plot(trees_geohisse[[i]], show.tip.label=F, x.lim=c(1,80))  
  x<-probs[[i]]
  x[,1]<-NULL
  x<-as.matrix(x)
  colors_list_for_states<-c("blue", "red", "darkorchid")
  nodelabels(pie=x[1:length(x[,1]),],cex=0.4, piecol=colors_list_for_states)
  axisPhylo()
  
}

dev.off()




