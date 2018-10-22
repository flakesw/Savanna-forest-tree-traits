## Phylogenetic analysis


#import species list and remove some duplicates
spList_orig <- read.csv("./raw data/lista_spp_plantas_families_sf_011618.csv")[-c(113, 130, 138, 176), ]
spList_orig[!(spList_orig$FG %in% c("S", "F", "G")), "FG"] <- "G"
spList_orig$FG <- droplevels(spList_orig$FG)

#make species list in proper format
spList <- data.frame(species = paste(spList_orig$Genus, spList_orig$Species, sep = " "),
                     genus = spList_orig$Genus,
                     family = spList_orig$Family)

library("phytools")                       # load the "phytools" package.
library("ape")
source("./S.PhyloMaker-master/R_codes for S.PhyloMaker")
example<-read.csv("./S.PhyloMaker-master/example.splist",header=T, sep = "\t")       # read in the example species list.    
phylo<-read.tree("./S.PhyloMaker-master/PhytoPhylo")      # read in the megaphylogeny.    
nodes<-read.csv("./S.PhyloMaker-master/nodes",header=T, sep = "\t")     # read in the nodes information of the megaphylogeny.    
result<-S.PhyloMaker(spList=spList, tree=phylo, nodes=nodes)      # run the function S.PhyloMaker.    
str(result)       # the structure of the ouput of S.PhyloMaker.    
par(mfrow=c(1,3),mar=c(0,0,1,0))       # show the phylogenies of the three scenarios.
plot(result$Scenario.1,cex=1.1,main="Scenario One")
plot(result$Scenario.2,cex=1.1,main="Scenario Two")
plot(result$Scenario.3,cex=1.1,main="Scenario Three")

tree_out <- result$Scenario.3

plot(tree_out, type = "fan")


#-------------------------------------------------------------------
# playing around with ancestral states -- does not work great 
x <- spList_orig$FG
names(x) <- paste(spList_orig$Genus, spList_orig$Species, sep = "_")

fit<-fastAnc(tree_out, x ,vars=TRUE,CI=TRUE)
fit

make.simmap(tree_out, x = x, model = "empirical")


cols <- ifelse(spList_orig$FG == "F", "blue", 
               ifelse(spList_orig$FG == "G", "green", "red"))

plotTree(tree_out,type="fan",fsize=.6, ftype="i", tip.color = cols, no.margin = TRUE)
## setEnv=TRUE for this type is experimental. please be patient with bugs
cols<-setNames(palette()[1:length(unique(x))],sort(unique(x)))
tiplabels(pie=to.matrix(x,sort(unique(x))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(tree)),fsize=0.8)

fitER<-rerootingMethod(tree_out, x, model="ER", type="discrete")

plotTree(tree_out,type="fan",fsize=.6, ftype="i", tip.color = cols, no.margin = TRUE)
nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, 
           piecol = c("blue", "red", "yellow"), cex = 0.6)
tiplabels(pie = to.matrix(x, sort(unique(x))), piecol = c("blue", "red", "yellow"), 
          cex = 0.3)

#simmap tutorial from phytools
anole.tree<-read.tree("http://www.phytools.org/eqg2015/data/anole.tre")
## plot tree
plotTree(anole.tree,type="fan",ftype="i")

svl<-read.csv("http://www.phytools.org/eqg2015/data/svl.csv",
              row.names=1)
## change this into a vector
svl<-as.matrix(svl)[,1]
svl

fit<-fastAnc(anole.tree,svl,vars=TRUE,CI=TRUE)
fit

obj<-contMap(anole.tree,svl,plot=FALSE)
plot(obj,type="fan",legend=0.7*max(nodeHeights(anole.tree)),
     fsize=c(0.7,0.9))


#example for discrete characters
data(anoletree)
## this is just to pull out the tip states from the
## data object - normally we would read this from file
x<-getStates(anoletree,"tips")
tree<-anoletree
rm(anoletree)
tree

plotTree(tree,type="fan",fsize=0.8,ftype="i")
## setEnv=TRUE for this type is experimental. please be patient with bugs
cols<-setNames(palette()[1:length(unique(x))],sort(unique(x)))
tiplabels(pie=to.matrix(x,sort(unique(x))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(tree)),fsize=0.8)

## estimate ancestral states under a ER model
fitER<-ace(x,tree,model="ER",type="discrete")
fitER
