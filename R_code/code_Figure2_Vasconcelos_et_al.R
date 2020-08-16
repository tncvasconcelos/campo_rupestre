#################################################################
############## Code to reproduce plots in Figure 3 ##############
#################################################################

# This is part of:
# Vasconcelos TNC, Alcantara S, Andrino CO, Forest F, Reginato M, Simon MF, Pirani JR 
# Fast diversification through a mosaic of evolutionary histories 
# characterizes the endemic flora of ancient Neotropical mountains. Proceedings B.


###############################################################
######################### BOXPLOTS ############################
###############################################################


pdf(file="Figure3a_base.pdf",width=10, height=10)

par(mfrow=c(1,2))
setwd("/ESM_data/GeoHiSSE/")

list.files(pattern="nodes_insitu_geosse.csv") -> geosse
sapply(geosse, read.csv, simplify=F) -> total_geosse
data.table::rbindlist(total_geosse) -> total_geosse
boxplot(total_geosse$ages~ total_geosse$n1,lwd=0.2, las=2,
        at=rank(tapply(total_geosse$ages, total_geosse$n1, max)), 
        col="tomato2", pch=20, ylim=c(0,30)) # 

wd <- "/Users/thaisvasconcelos/Desktop/revision_cr/files to resubmit"
setwd(paste(wd, "/ESM_data/BioGeoBEARS", sep=""))

# the following table contains node ages with highest likelihood of campo rupestre endemic 
# from all BioGeoBEARS analyses, including those on the whole seed-plant phylogeny

total<-read.csv("boxplot_biogeobears.csv") 

boxplot(total$ages~ total$group,lwd=0.2, las=2,
        at=sort(rank(tapply(total$ages, total$group, max))), col="tomato2", pch=20, ylim=c(0,30)) 

dev.off()



######################################################################
######################################################################

##### Rate through time from BAMM ####

library(ape)
library(BAMMtools)
library(RColorBrewer)

setwd("ESM_data/BAMM_files")

# Calliandra 
tree_calliandra<-read.tree("BAMM_calliandra_bamm_tree.txt")
ed_calliandra <- getEventData(tree_calliandra, "BAMM_calliandra_event_data.txt", burnin = 0.2)
endemic_calliandra <- getCladeRates(ed_calliandra, node = 168) # node 168 for the MRCA of Calliandra endemic radiation (radiation number 5 in fig.1)
mean(endemic_calliandra$lambda) # 0.76431 # mean speciation rate 

# Chamaecrista 
tree_chamaecrista <- read.tree("BAMM_chamaecrista_bamm_tree.txt")
ed_chamaecrista <- getEventData(tree_chamaecrista, "BAMM_chamaecrista_event_data.txt", burnin = 0.2)
endemic_chamaecrista1 <- getCladeRates(ed_chamaecrista, node = 172) # node 172 for the MRCA of Chamaecrista endemic radiation (radiation number 6 in fig. 1)
mean(endemic_chamaecrista1$lambda) # 0.90564
endemic_chamaecrista2 <- getCladeRates(ed_chamaecrista, node = 114) # node 114 for the MRCA of Chamaecrista endemic radiation (radiation number 7 in fig. 1)
mean(endemic_chamaecrista2$lambda) # 1.24245

# Trimezieae 
tree_trimezieae <- read.tree("BAMM_trimezieae_bamm_tree.txt")
ed_trimezieae  <- getEventData(tree_trimezieae , "BAMM_trimezieae_event_data.txt", burnin = 0.2)
endemic_trimezieae  <- getCladeRates(ed_trimezieae , node = 61) # node 61 for the MRCA of Pseudotrimezia
mean(endemic_trimezieae$lambda) # 0.5391305

# Mimosa 
tree_mimosa <- read.tree("BAMM_mimosa_bamm_tree.txt")
ed_mimosa  <- getEventData(tree_mimosa , "BAMM_mimosa_event_data.txt", burnin = 0.2)
endemic_mimosa  <- getCladeRates(ed_mimosa , node = 389) # node 389 for the MRCA of Mimosa clade O
mean(endemic_mimosa$lambda) # 1.401293

# Diplusodon 
tree_diplusodon <- read.tree("BAMM_diplusodon_bamm_tree.txt")
ed_diplusodon  <- getEventData(tree_diplusodon , "BAMM_diplusodon_event_data.txt", burnin = 0.2)
endemic_diplusodon  <- getCladeRates(ed_diplusodon) # the whole tree is a radiation
mean(endemic_diplusodon$lambda) # 0.4844791

# Lychnophorinae 
tree_lychnophorinae <- read.tree("BAMM_lychnophorinae_bamm_tree.txt")
ed_lychnophorinae  <- getEventData(tree_lychnophorinae, "BAMM_lychnophorinae_event_data.txt", burnin = 0.2)
endemic_lychnophorinae  <- getCladeRates(ed_lychnophorinae, node=80) # node 80 for the MRCA of Lychnophorinae
mean(endemic_lychnophorinae$lambda) # 1.084128

# Marcetia 
tree_marcetia <- read.tree("BAMM_marcetia_bamm_tree.txt")
ed_marcetia  <- getEventData(tree_marcetia, "BAMM_marcetia_event_data.txt", burnin = 0.2)
endemic_marcetia  <- getCladeRates(ed_marcetia, node=93) # node 93 for the MRCA of Marcetia
mean(endemic_marcetia$lambda) # 0.5166609

# Paepalanthus s.l.
tree_paepalanthus <- read.tree("BAMM_paepalanthus_bamm_tree.txt")
ed_paepalanthus  <- getEventData(tree_paepalanthus, "BAMM_paepalanthus_event_data.txt", burnin = 0.2)
endemic_paepalanthus1  <- getCladeRates(ed_paepalanthus, node=185) # node 185 for the MRCA of Paepalanthus s.s.
mean(endemic_paepalanthus1$lambda) # 1.046901
endemic_paepalanthus2  <- getCladeRates(ed_paepalanthus, node=313) # node 313 for the MRCA of Actinocephalus
mean(endemic_paepalanthus2$lambda) # 1.066719

# Cattleya
tree_cattleya <- read.tree("BAMM_cattleya_bamm_tree.txt")
ed_cattleya  <- getEventData(tree_cattleya, "BAMM_cattleya_event_data.txt", burnin = 0.2)
endemic_cattleya  <- getCladeRates(ed_cattleya, node=103) # node 103 for the MRCA of Cattleya sect. Parviflora
mean(endemic_cattleya$lambda) # 1.499022

# Velloziaceae 
tree_velloziaceae <- read.tree("BAMM_velloziaceae_bamm_tree.txt")
ed_velloziaceae  <- getEventData(tree_velloziaceae, "BAMM_velloziaceae_event_data.txt", burnin = 0.2)
endemic_velloziaceae1  <- getCladeRates(ed_velloziaceae, node=157) # Vellozia
mean(endemic_velloziaceae1$lambda) # 0.1786386
endemic_velloziaceae2  <- getCladeRates(ed_velloziaceae, node=254) # Barbacenia
mean(endemic_velloziaceae2$lambda) # 0.299503

# Rate through time from the stem node
# (this can take some time) 

pdf(file="Figure3b_base")
plotRateThroughTime(ed_mimosa, node = 388, avgCol="orange", intervalCol = "orange", opacity=0.002, start.time = 36, ylim=c(0, 2))
plotRateThroughTime(ed_calliandra, node = 167, avgCol="navy", intervalCol = "navy", opacity=0.002,  start.time = 36, ylim=c(0, 2), add=T)
plotRateThroughTime(ed_lychnophorinae, node = 79, avgCol="limegreen", intervalCol = "limegreen",opacity=0.002, start.time = 36, ylim=c(0, 2), add=T)
plotRateThroughTime(ed_chamaecrista, node = 170, avgCol="chocolate1", intervalCol = "chocolate1", opacity=0.002, start.time = 36, ylim=c(0, 2), add=T)
plotRateThroughTime(ed_chamaecrista, node = 113,avgCol="red", intervalCol = "red", opacity=0.002, start.time = 36, ylim=c(0, 2), add=T)
plotRateThroughTime(ed_trimezieae, node = 60, avgCol="darkred", intervalCol = "darkred", opacity=0.002, start.time = 36, ylim=c(0, 2),add=T)
plotRateThroughTime(ed_velloziaceae, node = 156, avgCol="darkorchid3", intervalCol = "darkorchid3", opacity=0.002, start.time = 36, ylim=c(0, 2),add=T)
plotRateThroughTime(ed_velloziaceae, node = 253, avgCol="darkolivegreen3", intervalCol = "darkolivegreen3", opacity=0.002, start.time = 36, ylim=c(0, 2), add=T)
plotRateThroughTime(ed_marcetia, node = 92, avgCol="violetred1", intervalCol = "violetred1", opacity=0.002, start.time = 36, ylim=c(0, 2), add=T)
plotRateThroughTime(ed_paepalanthus, node = 184, avgCol="darkgreen", intervalCol = "darkgreen",opacity=0.002, start.time = 36, ylim=c(0, 2), add=T)
plotRateThroughTime(ed_paepalanthus, node = 312, avgCol="steelblue", intervalCol = "steelblue", opacity=0.002, start.time = 36, ylim=c(0, 2), add=T)
plotRateThroughTime(ed_cattleya, node = 102, avgCol="bisque3", intervalCol = "bisque3", opacity=0.002, start.time = 36, ylim=c(0, 2), add=T)
plotRateThroughTime(ed_diplusodon, avgCol="yellow1", intervalCol = "yellow1", opacity=0.002, start.time = 36, ylim=c(0, 2), add=T)

dev.off()

