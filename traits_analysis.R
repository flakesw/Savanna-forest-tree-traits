## Bivariate and ordination analysis

# load libraries
library("ggplot2")
library("plyr")
library("vegan")
library("multcomp")
library("multcompView")
library("effects")
library("maptools")
library("ape")
library("cluster")
library("phytools")
library("Hotelling")

# setwd("C:\\Users\\Sam\\Google Drive\\Projects\\Savanna traits")
set.seed(45750765)


#----------------------------------------------------------------------------
# Function to plot elliptical hulls for phylogenetic PCA,
# modified from vegan::ordiellipse. I removed a bunch of functionality,
# so now it only draws elliptical hulls. I really forget why I had to butcher it.
#----------------------------------------------------------------------------
phyellipse <- function (ord, groups, display = "sites",
                        col = NULL, 
                        show.groups, choices = c(1,2), ...) 
{

  pts <- ord$S
  pts <- as.matrix(pts)
  pts <- pts[, choices, drop = FALSE]
  take <- groups %in% show.groups
  pts <- pts[take, , drop = FALSE]

  tmp <- ellipsoidhull(pts)
  lines(predict(tmp), col = col)
  
}


#----------------------------------------------------------------------------
# Import and wrangle data
#----------------------------------------------------------------------------
#read data
clean_species <- read.csv("./clean data/clean_species2020-06-09.csv")
spList_orig <- read.csv("./raw data/sp_info_with_splink.csv")
spList_orig <- spList_orig[spList_orig$Code %in% clean_species$Code, ]
write.csv(spList_orig, "./Clean data/species_list_for_appendix.csv")


plot_data <- read.csv("./Clean data/plot_data2020-06-09.csv")
plot_data <- subset(plot_data, !is.na(CSA_30))
names(plot_data)[21] <- "FG"
bark_model <- readRDS("./clean data/bark_model.RDS")

clean_species$FG <- factor(clean_species$FG,levels(clean_species$FG)[c(3,2,1)])
# clean_species[clean_species$Code %in% c("GOPO", "TAOC", "RAGU", "MYLI"),]$FG <- "S" #to test influence of these species
clean_species$Leaf_size <- clean_species$Leaf_size/100 
clean_species$Total_leaf_size <- clean_species$Total_leaf_size/100 

clean_species_reduced <- clean_species[complete.cases(clean_species), ]
# clean_species_reduced <- clean_species_reduced[clean_species_reduced$Habit %in% c("Tree", "tree", "Treelet", "treelet"), ]

clean_species_trans <- clean_species
clean_species_reduced_orig <- clean_species_reduced
clean_species_reduced_trans <- clean_species_reduced

#vars to log-trans
log_trans <- c(4, 5, 6, 9, 11)

#log-trans predictor variables
clean_species_trans[, log_trans] <- log(clean_species[, log_trans])
clean_species_reduced_trans[, log_trans] <- log(clean_species_reduced_trans[, log_trans])


# Bring in phylogeny
tree <- read.tree("./clean data/phylogeny_v")

#relabel tree with code instead of names, to use later
tree_code <- tree
tree_code$tip.label  <- spList_orig[match(tree$tip.label, gsub(" ", "_", spList_orig$New.name)), "Code"] 


#plot the tree
tiff(filename = "./plots/phylogeny_pruned.tiff",
    width = 7,
    height = 7,
    units = "in",
    pointsize = 7.8,
    compression = "lzw",
    res = 600)
  
  #deal with the weird name problems
  fg <- spList_orig[match(tree$tip.label, gsub(" ", "_", spList_orig$New.name)), "classification66"]

  # fg <- fg[!is.na(fg)]
  
  color_groups <- data.frame(group = as.character(c("F", "G", "S")),
                             col = as.character(c("black", "#1b9e77", "#E69F00")))
  color_points <- as.character(color_groups[match(fg, color_groups$group), "col"])
  
  par(mar = c(0, 0, 0, 0))
  
  plot(tree, type = "fan",
       show.tip.label = TRUE,
       tip.color = color_points,
       show.node.label = TRUE)
  legend(x = par()$usr[1], y = par()$usr[4], 
         legend = c("Forest", "Generalist", "Savanna"), 
         col = c("black", "#1b9e77", "#E69F00"), 
         pt.bg = c("black", "#1b9e77", "#E69F00"), 
         pch = 21,
         cex = 1.3,
         bty = "n")
dev.off()

#------------------------------------------------------------------------------
# What proportion of tree species sampled? For results section
#------------------------------------------------------------------------------
#calculate using both clean_species and clean_species_reduced
total_ba <- sum(plot_data$CSA_30_expand, na.rm = TRUE)
sampled_ba <- sum(plot_data[plot_data$Code %in% clean_species$Code, ]$CSA_30_expand, na.rm = TRUE)
sampled_ba/total_ba

plot_ba <- aggregate(plot_data[, c("P", "CSA_30_expand")], by = list(plot_data$P), FUN = function(x){sum(x, na.rm=TRUE)})
sampled_plot_ba <- aggregate(plot_data[plot_data$Code %in% clean_species$Code, c("P", "CSA_30_expand")], 
                             by = list(plot_data[plot_data$Code %in% clean_species$Code, ]$P), FUN = function(x){sum(x, na.rm=TRUE)})
sampled_plot_ba$CSA_30_expand / plot_ba$CSA_30_expand

#------------------------------------------------------------------------------
# Bivariate analyses
#------------------------------------------------------------------------------
traits_names <- names(clean_species[-c(1:3, 5, 15, 16)])
traits_names <- traits_names[c(1,2,7,3,4,5,8,9,10,6)]
# Calculate phylogenetic signal, ANOVA, Tukey's post-hoc, phylogenetic ANOVA
log_trans_traits <- c(1, 2, 6, 3) # which traits to test transformed

#intialize a table to catch all the results
pvals_table <- data.frame(matrix(ncol = 11, nrow = 9))
#t for tukey, p for phylogenetic comparison
names(pvals_table) <- c("Trait", "tG-S", "tF-S", "tF-G", "pG-S", "pF-S", "pF-G", "ANOVA", "phylANOVA", "lambda", "n")

for (i in 1:length(traits_names)){

  #subset the tree to only match species with data for the trait
  tree_pruned <- drop.tip(tree_code, setdiff(tree_code$tip.label, clean_species_trans[!is.na(clean_species_trans[, traits_names[i]]), ]$Code))
  tree_pruned$tip.label <- droplevels(tree_pruned$tip.label)
  
  #reorder species traits to match phylogeny
  clean_species_ordered <- clean_species_trans[match(tree_pruned$tip.label, clean_species_trans$Code), ] #reorder and subset to match phylogeny
  
  pvals_table[i, 11] <- nrow(clean_species_ordered)
  
  #trait name
  pvals_table[i, 1] <- traits_names[i]
  #fit an ANOVA
  model <- lm(clean_species_ordered[, traits_names[i]] ~ as.factor(clean_species_ordered$FG))
  aov <- aov(model)
  #ANOVA p-value from F-test
  pvals_table[i, 8] <- summary(aov)[[1]][[5]][[1]] #this is dumb; is this really the right way to extract the p-value?!
  
  #calculate post-hoc phylogenetically-corrected differences and add letters to figure
  tuke <- TukeyHSD(aov)
  pvals <- tuke$`as.factor(clean_species_ordered$FG)`[, 4]
  pvals_table[i, c(2:4)] <- pvals
  
  #phylogenetic ANOVA
  phyl_aov_results <- phylANOVA(tree = tree_pruned, x = clean_species_ordered$FG, 
                                y = clean_species_ordered[, traits_names[i]], 
                                nsim=1000, posthoc=TRUE, p.adj="holm")
  
  pvals_table[i, 9] <- phyl_aov_results$Pf
  
  pvals <- c(phyl_aov_results$Pt[2,1], phyl_aov_results$Pt[3,1], phyl_aov_results$Pt[3,2]) 
  
  pvals_table[i, c(5:7)] <- pvals
  
  
  #calculate the phylogenetic signal in the trait
  sig <- phylosig(tree = tree_pruned, x = clean_species_ordered[, traits_names[i]], 
                  method="lambda", test=TRUE, nsim=1000)
  pvals_table[i, 10] <- sig$lambda
}

write.csv(pvals_table, "./Model output/FG_traits_pvals.csv")

#---------------------------
# Figure 2: trait differences between functional groups

traits_names_clean <- c(expression(paste("Leaf size (cm"^"2", ")")),
                        "Leaf thickness (mm)",
                        expression(paste("Specific leaf area (cm"^"2", " g"^"-1", ")")),
                        "Max height (m)",
                        "Height at 5 cm dia. (m)", 
                        "       Crown ratio\nat 5cm dia. (unitless)",
                        "  Bark thickness\nat 5cm dia. (mm)",
                        "  Bark thickness\nat 8mm dia. (mm)",
                        expression(paste("Wood density (g cm"^"-3",")")),
                        "         Light code\nat 5 cm dia. (unitless)")

tiff(filename="./plots/traits_by_FG.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 7, 
     height=4, 
     pointsize=9, 
     res=600)

par(mfrow = c(2, 5),
    mar = c(2,6,2,0.5),
    oma = c(1, 0, 0, 0))

for (i in 1:length(traits_names)){
  # model <- lm(clean_species_reduced[, traits_names[i]] ~ as.factor(clean_species_reduced$FG))
  noise <- rnorm(nrow(clean_species_trans), 0, .02)

  means <- aggregate(clean_species_trans[, traits_names[i]], by = list(as.factor(clean_species_trans$FG)),
                       FUN = mean, na.rm = TRUE)
    
  
  if(!(i %in% log_trans_traits)){
    ylim <- c((0.8 * min(clean_species_trans[, traits_names[i]], na.rm = TRUE, na.rm = TRUE)), 
              (1.2 * max(clean_species_trans[, traits_names[i]], na.rm = TRUE)))
    yrange <- abs(ylim[2]-ylim[1])
    
    plot(clean_species_trans[, traits_names[i]] ~ I(as.numeric(clean_species_trans$FG) + noise),
         xlab = "",
         ylab = "",
         xaxt = 'n',
         yaxt = 'n',
         xlim = c(0.8, 3.2),
         ylim = ylim,
         cex.lab = 1.3
       )
    
    axis(side = 2, las = 2)
    
    
  }else{ # for log-trans traits

    means[, 2] <- (means[, 2])
    
    ylim <- c(log(0.8 * min((clean_species[, traits_names[i]]), na.rm = TRUE)), 
              log(1.3 * max((clean_species[, traits_names[i]]), na.rm = TRUE)))
    
    if(i == 1){
      ylim <- c(log(0.8 * min((clean_species[, traits_names[i]]), na.rm = TRUE)), 
              log(1.8 * max((clean_species[, traits_names[i]]), na.rm = TRUE)))
    }
    
    yrange <- abs(ylim[2]-ylim[1])
    
    plot((clean_species_trans[, traits_names[i]]) ~ I(as.numeric(clean_species_trans$FG) + noise),
         xlab = "",
         ylab = "",
         xaxt = 'n',
         xlim = c(0.8, 3.2),
         ylim = ylim,
         cex.lab = 1.3,
         yaxt = "n"
    )
    
    ax_labs <- axisTicks(log10(range(exp(clean_species_trans[, traits_names[i]]), na.rm = TRUE)), log = TRUE, n = 5)

    axis(side = 2, at = log(ax_labs), labels = ax_labs, las = 2)
    
  }
  
  points(means[, 2] ~ c(1,2,3), pch = 18, cex = 2)
  lines(means[, 2] ~ c(1,2,3), lty = 1, lwd = 2)
  
  mtext(text = levels(clean_species_trans$FG), side = 1, at = c(1,2,3), line = 1)
  mtext(text = paste0("(", letters[i], ")") , side = 3, at = 0.4, line = .3)
  mtext(text = traits_names_clean[i], side = 2, line = 2.6)
  
  pvals <- as.numeric(pvals_table[i, c(5:7)])
  names(pvals) <- c("S-G", "S-F", "G-F")
  pvals <- ifelse(pvals < (0.05), TRUE, FALSE)
  tuke_letters <- multcompLetters(pvals, reversed=FALSE)$Letters
  
  text(x = c(1,2,3), y = ylim[2] - (ylim[2] - ylim[1])*.03, labels = tuke_letters, cex = 1.2)

}

dev.off()



#------------------------------------------------------------------------------
# Same figure as a boxplot
#------------------------------------------------------------------------------

tiff(filename="./plots/Figure_2_traits_by_FG_boxplot.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 7, 
     height=4, 
     pointsize=9, 
     res=600)

par(mfrow = c(2, 5),
    mar = c(2,5.2,2,0.5),
    oma = c(1, .3, 0, 0))

for (i in 1:length(traits_names)){
  # model <- lm(clean_species_reduced[, traits_names[i]] ~ as.factor(clean_species_reduced$FG))
  noise <- rnorm(nrow(clean_species_trans), 0, .02)
  
  means <- aggregate(clean_species_trans[, traits_names[i]], by = list(as.factor(clean_species_trans$FG)),
                     FUN = mean, na.rm = TRUE)
  
  
  if(!(i %in% log_trans_traits)){
    ylim <- c((0.8 * min(clean_species_trans[, traits_names[i]], na.rm = TRUE, na.rm = TRUE)), 
              (1.2 * max(clean_species_trans[, traits_names[i]], na.rm = TRUE)))
    yrange <- abs(ylim[2]-ylim[1])
    
    boxplot(clean_species_trans[, traits_names[i]] ~ clean_species_trans$FG,
         xlab = "",
         ylab = "",
         xaxt = 'n',
         yaxt = 'n',
         xlim = c(0.5, 3.5),
         ylim = ylim,
         cex.lab = 1.3
    )
    
    axis(2, las = 2)
    
  }else{ # for log-trans traits
    
    # means[, 2] <- (means[, 2])
    
    ylim <- c(log(0.8 * min((clean_species[, traits_names[i]]), na.rm = TRUE)), 
              log(1.3 * max((clean_species[, traits_names[i]]), na.rm = TRUE)))
    
    if(i == 1){
      ylim <- c(log(0.8 * min((clean_species[, traits_names[i]]), na.rm = TRUE)), 
                log(1.8 * max((clean_species[, traits_names[i]]), na.rm = TRUE)))
    }
    
    yrange <- abs(ylim[2]-ylim[1])
    
    boxplot(clean_species_trans[, traits_names[i]] ~ clean_species_trans$FG,
         xlab = "",
         ylab = "",
         xaxt = 'n',
         xlim = c(0.5, 3.5),
         ylim = ylim,
         cex.lab = 1.3,
         yaxt = "n"
    )
    
    ax_labs <- axisTicks(log10(range(exp(clean_species_trans[, traits_names[i]]), na.rm = TRUE)), log = TRUE, n = 5)
    
    axis(side = 2, at = log(ax_labs), labels = ax_labs, las = 2)
  }
  # 
  # points(means[, 2] ~ c(1,2,3), pch = 18, cex = 2)
  # lines(means[, 2] ~ c(1,2,3), lty = 1, lwd = 2)
  
  mtext(text = levels(clean_species_trans$FG), side = 1, at = c(1,2,3), line = 1)
  mtext(text = paste0("(", letters[i], ")") , side = 3, at = 0.4, line = .3)
  mtext(text = traits_names_clean[i], side = 2, line = 2.3)
  
  pvals <- as.numeric(pvals_table[i, c(5:7)])
  names(pvals) <- c("S-G", "S-F", "G-F")
  pvals <- ifelse(pvals < (0.05), TRUE, FALSE)
  tuke_letters <- multcompLetters(pvals, reversed=FALSE)$Letters
 
  text(x = c(1,2,3), y = ylim[2] - (ylim[2] - ylim[1])*.03, labels = tuke_letters, cex = 1.2)
  # }
  
}

dev.off()


#------------------------------------------------------------------------------
#PCA on traits
#------------------------------------------------------------------------------

pca_data <- scale(clean_species_reduced_trans[, c("Leaf_size", "Leaf_thickness", "Max_height",
                                            "Height_at_5cm", "Crown_ratio", "SLA",
                                            "Bark_at_5cm", "Bark_at_8mm", "Wood_density", "Light_at_5cm")])

rownames(pca_data) <- clean_species_reduced_trans$Code

pca_groups <- clean_species_reduced_trans[, c("FG")]

pca_sp_names <- droplevels(spList_orig[match(clean_species_reduced_trans$Code, spList_orig$Code), "NewCode"])

#run the PCA
pca_results <- rda(pca_data, scale = TRUE)
trait_names <- gsub("_", " ", rownames(pca_results$CA$v))
sp_scores <- vegan::scores(pca_results, choices = c(1:10), display = c("sites"))
trait_scores <- vegan::scores(pca_results, choices = c(1:10), display = c("species"))

summary(pca_results)

biplot(pca_results, choices = c(1,2))

set.seed(45750766)

tiff(filename="./plots/Figure_3_traits_PCA_axes_1_2.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 7, 
     height=7, 
     pointsize=9, 
     res=600)
color_groups <- data.frame(group = as.character(c("F", "G", "S")),
                           col = as.character(c("black", "#1b9e77", "#E69F00")))
color_points <- as.character(color_groups[match(pca_groups, color_groups$group), "col"])

plot(NA,
     xlab = "PC1", ylab = "PC2",
     xlim = c( (min(c(sp_scores[, 1], trait_scores[, 1]) - .5)), 
               (max(c(sp_scores[, 1], trait_scores[, 1])) + 0.4)),
     ylim = c( (min(c(sp_scores[, 2], trait_scores[, 2]) - .5)), 
               (max(c(sp_scores[, 2], trait_scores[, 2])) + 0.4)))

abline(h = 0)
abline(v = 0)

x<-pointLabel(x = sp_scores[, 1], y = sp_scores[, 2], labels = pca_sp_names, col = color_points, cex = 0.7)

for(i in 1:3){
  ordiellipse(pca_results, pca_groups, show.groups = color_groups[i, 1], 
              col = as.character(color_groups[i, 2]), choices = c(1,2),
              kind = "sd", conf = 0.90)
}

legend(x = par()$usr[1], y = par()$usr[4], 
       legend = c("Forest", "Generalist", "Savanna"), 
       col = c("black", "#1b9e77", "#E69F00"), 
       pt.bg = c("black", "#1b9e77", "#E69F00"), 
       pch = 21,
       cex = 1.5,
       bty = "n")

segments(x0 = 0, y0 = 0, x1 = trait_scores[, 1], y1 = trait_scores[, 2])
yoffsets <- c(0.12, -.05, 0.03, .05, -.06, .05, -.08, .05, -.15, 0.02)
xoffsets <- c(-0.15, -.33, 0, 0.2, -.05, 0, -.25, -.3, .27, -.3)

text(x = trait_scores[, 1]+xoffsets, y = trait_scores[, 2]+yoffsets, labels = trait_names, cex = 1.3)

dev.off()

#################
# PCs 1 + 3


tiff(filename="./plots/Figure_S2_traits_PCA_axes_1_3.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 7, 
     height=7, 
     pointsize=9, 
     res=600)

color_groups <- data.frame(group = as.character(c("F", "G", "S")),
                           col = as.character(c("black", "#1b9e77", "#E69F00")))
color_points <- as.character(color_groups[match(pca_groups, color_groups$group), "col"])

plot(NA,
     xlab = "PC1", ylab = "PC3",
     xlim = c( (min(c(sp_scores[, 1], trait_scores[, 1]) - .5)), 
               (max(c(sp_scores[, 1], trait_scores[, 1])) + 0.4)),
     ylim = c( (min(c(sp_scores[, 3], trait_scores[, 3]) - .5)), 
               (max(c(sp_scores[, 3], trait_scores[, 3])) + 0.4)),
     yaxt = "n")

abline(h = 0)
abline(v = 0)

axis(2, las = 2)

pointLabel(x = sp_scores[, 1], y = sp_scores[, 3], labels = pca_sp_names, col = color_points, cex = 0.7)



for(i in 1:3){
  ordiellipse(pca_results, pca_groups, show.groups = color_groups[i, 1], 
              col = as.character(color_groups[i, 2]), choices = c(1,3),
              kind = "sd", conf = 0.90)
}

legend(x = par()$usr[1], y = par()$usr[4], 
       legend = c("Forest", "Generalist", "Savanna"), 
       col = c("black", "#1b9e77", "#E69F00"), 
       pt.bg = c("black", "#1b9e77", "#E69F00"), 
       pch = 21,
       cex = 1.5,
       bty = "n")

segments(x0 = 0, y0 = 0, x1 = trait_scores[, 1], y1 = trait_scores[, 3])
yoffsets <- c(0, -.05, -.05, 0, .05, 0, -.04, -.04, -.1, 0)
xoffsets <- c(0, 0, 0.23, 0, -.05, 0, -.25, -.25, 0, -0.2)

text(x = trait_scores[, 1]+xoffsets, y = trait_scores[, 3]+yoffsets, labels = trait_names, cex = 1.3)

dev.off()


# #---------------------------------------------------------------------------------------------
# Tables for supplement
#--------------------------------------------------------------------------------------------

species_scores_table <- as.data.frame(cbind(as.character(pca_sp_names), 
                                            as.character(clean_species_reduced$FG),
                                            sp_scores))
names(species_scores_table)[c(1,2)] <- c("Species Code", "Functional Type")
write.csv(species_scores_table, "./Model output/species_score_table.csv")

#--------------------------------------------------------------------------
# Checking differences between functional groups
#---------------------------------------------------------------------------
#Hotelling's T2 for pairwise differences

hotelling_results <- data.frame("H2" = numeric(3))

groups <- matrix(data = c("S", "G", "F", "G", "F", "S"), nrow = 3, ncol = 2, byrow = TRUE)

pca_results_out <- cbind(as.data.frame(sp_scores), pca_groups)
names(pca_results_out)[11] <- "FG"

for(i in 1:3){
  h2_pca <- hotelling.test(x = pca_results_out[pca_results_out$FG == groups[i,1], c(1:2)], 
                 y = pca_results_out[pca_results_out$FG == groups[i,2], c(1:2)],
                 perm = FALSE)
  hotelling_results[i, 1] <- h2_pca$pval
}

#PERMANOVA analysis
anosim_data <- as.matrix(clean_species_reduced_trans[, c("Leaf_size", "Leaf_thickness", "Max_height",
                                         "Height_at_5cm", "Crown_ratio", "SLA",
                                         "Bark_at_5cm", "Bark_at_8mm", "Wood_density", "Light_at_5cm")])
functional_group_anosim <- adonis(anosim_data ~ clean_species_reduced_trans$FG, method = "euclidean")
summary(functional_group_anosim)
print(functional_group_anosim)


#*******************************************************************************************
# Estimating trait distributions for plots
#*******************************************************************************************
#create new columns to fill
plot_data$Leaf_size <- NA
plot_data$Leaf_thickness <- NA
plot_data$Height_at_5cm <- NA
plot_data$Max_height <- NA
plot_data$Crown_ratio <- NA
plot_data$SLA <- NA
plot_data$Bark_at_5cm <- NA
plot_data$Bark_at_8mm <- NA
plot_data$est_bark_thickness <- NA
plot_data$Wood_density <- NA
plot_data$Light_at_5cm <- NA


was_measured <- plot_data$Code %in% clean_species$Code
# estimate traits for entries in the plot data

for(i in 1:nrow(plot_data)){#this is surely the worst way to do this
  
  if(was_measured[i]){
    plot_data$Leaf_size[i] <- clean_species[clean_species$Code == as.character(plot_data$Code[i]), "Total_leaf_size"]
    plot_data$Leaf_thickness[i] <- clean_species[clean_species$Code == as.character(plot_data$Code[i]), "Leaf_thickness"]
    plot_data$Height_at_5cm[i] <- clean_species[clean_species$Code == as.character(plot_data$Code[i]), "Height_at_5cm"]
    plot_data$Max_height[i] <- clean_species[clean_species$Code == as.character(plot_data$Code[i]), "Max_height"]
    plot_data$Crown_ratio[i] <- clean_species[clean_species$Code == as.character(plot_data$Code[i]), "Crown_ratio"]
    plot_data$SLA[i] <- clean_species[clean_species$Code == as.character(plot_data$Code[i]), "SLA"]
    plot_data$Bark_at_5cm[i] <- clean_species[clean_species$Code == as.character(plot_data$Code[i]), "Bark_at_5cm"]
    plot_data$Bark_at_8mm[i] <- clean_species[clean_species$Code == as.character(plot_data$Code[i]), "Bark_at_8mm"]
    plot_data$est_bark_thickness[i] <- predict(bark_model, newdata = list(max_d30 = plot_data$Max_D30[i], Code = plot_data$Code[i]))
    plot_data$Wood_density[i] <- clean_species[clean_species$Code == as.character(plot_data$Code[i]), "Wood_density"]
    plot_data$Light_at_5cm[i] <- clean_species[clean_species$Code == as.character(plot_data$Code[i]), "Light_at_5cm"]
  } 
  
}
# names(plot_data)
# trait_summary <- aggregate(plot_data[, c(9, 21:32)], by = list(plot_data$Classification66, plot_data$Life.Form), FUN = mean, na.rm = TRUE)

#--------------------------------------------------------------------------
# Aggregate means for each plot, weighted by CSA of trees
#--------------------------------------------------------------------------
#all trees, not just shrubs (had previously excluded shrubs)

plot_agg_traits <- data.frame(Plot = as.integer(seq(1:30)),
                              BA = numeric(30),
                              Leaf_size = numeric(30),
                              Leaf_thickness = numeric(30),
                              Max_height = numeric(30),
                              Height_at_5cm = numeric(30),
                              Crown_ratio = numeric(30),
                              SLA = numeric(30),
                              Bark_at_5cm = numeric(30),
                              Bark_at_8mm = numeric(30),
                              est_bark_thickness = numeric(30),
                              Wood_density = numeric(30),
                              Light_at_5cm = numeric(30),
                              Light_code = numeric(30),
                              Height_mean = numeric(30),
                              Density = numeric(30),
                              QMD = numeric(30),
                              SavannaBA = numeric(30),
                              SavannaDens = numeric(30),
                              SavannaQMD = numeric(30),
                              GeneralistBA = numeric(30),
                              GeneralistDens = numeric(30),
                              GeneralistQMD = numeric(30),
                              ForestBA = numeric(30),
                              ForestDens = numeric(30),
                              ForestQMD = numeric(30)
                              )


for ( i in 1:30){
  trees_select <- plot_data[plot_data$P == i & !is.na(plot_data$CSA_30_expand), ]

  plot_agg_traits$Leaf_size[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Leaf_size), ]$Leaf_size, w = trees_select[!is.na(trees_select$Leaf_size), ]$CSA_30_expand, na.rm = TRUE)
  plot_agg_traits$Leaf_thickness[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Leaf_thickness), ]$Leaf_thickness, w = trees_select[!is.na(trees_select$Leaf_thickness), ]$CSA_30_expand, na.rm = TRUE)
  plot_agg_traits$Height_at_5cm[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Height_at_5cm), ]$Height_at_5cm, w = trees_select[!is.na(trees_select$Height_at_5cm), ]$CSA_30_expand, na.rm = TRUE)
  plot_agg_traits$Max_height[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Max_height), ]$Max_height, w = trees_select[!is.na(trees_select$Max_height), ]$CSA_30_expand, na.rm = TRUE)
  plot_agg_traits$Crown_ratio[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Crown_ratio), ]$Crown_ratio, w = trees_select[!is.na(trees_select$Crown_ratio), ]$CSA_30_expand, na.rm = TRUE)
  plot_agg_traits$SLA[i] <- weighted.mean(x = trees_select[!is.na(trees_select$SLA), ]$SLA, w = trees_select[!is.na(trees_select$SLA), ]$CSA_30_expand, na.rm = TRUE)
  plot_agg_traits$Bark_at_5cm[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Bark_at_5cm), ]$Bark_at_5cm, w = trees_select[!is.na(trees_select$Bark_at_5cm), ]$CSA_30_expand, na.rm = TRUE)
  plot_agg_traits$Bark_at_8mm[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Bark_at_8mm), ]$Bark_at_8mm, w = trees_select[!is.na(trees_select$Bark_at_8mm), ]$CSA_30_expand, na.rm = TRUE)
  plot_agg_traits$est_bark_thickness[i] <- weighted.mean(x = trees_select[!is.na(trees_select$est_bark_thickness), ]$est_bark_thickness, w = trees_select[!is.na(trees_select$est_bark_thickness), ]$CSA_30_expand, na.rm = TRUE)
  plot_agg_traits$Wood_density[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Wood_density), ]$Wood_density, w = trees_select[!is.na(trees_select$Wood_density), ]$CSA_30_expand, na.rm = TRUE)
  plot_agg_traits$Light_at_5cm[i] = weighted.mean(x = trees_select[!is.na(trees_select$Light_at_5cm), ]$Light_at_5cm, w = trees_select[!is.na(trees_select$Light_at_5cm), ]$CSA_30_expand, na.rm = TRUE)
  plot_agg_traits$Light_code[i] = weighted.mean(x = trees_select[!is.na(trees_select$Light.code), ]$Light.code, w = trees_select[!is.na(trees_select$Light.code), ]$CSA_30_expand, na.rm = TRUE)
  plot_agg_traits$Height_mean[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Ht), ]$Ht, w = trees_select[!is.na(trees_select$Ht), ]$CSA_30_expand, na.rm = TRUE)

  #some other useful structure variables
  plot_agg_traits$BA[i] <- sum(trees_select$CSA_BH_expand, na.rm = TRUE) / 10000 * 10
  plot_agg_traits$Density[i] <- sum(trees_select$Density_expand, na.rm = TRUE) * 10
  plot_agg_traits$QMD[i] <- sqrt((plot_agg_traits$BA[i] / plot_agg_traits$Density[i]) / 0.00007854)
  plot_agg_traits$SavannaBA[i] <- sum(trees_select[trees_select$classification66 == "S", ]$CSA_BH_expand, na.rm = TRUE) / 10000 * 10
  plot_agg_traits$SavannaDens[i] <- sum(trees_select[trees_select$classification66 == "S", "Density_expand"], na.rm= TRUE)  * 10 
  plot_agg_traits$SavannaQMD[i] <- sqrt((plot_agg_traits$SavannaBA[i] / plot_agg_traits$SavannaDens[i]) / 0.00007854)
  plot_agg_traits$GeneralistBA[i] <- sum(trees_select[trees_select$classification66 == "G", ]$CSA_BH_expand, na.rm = TRUE) / 10000 * 10
  plot_agg_traits$GeneralistDens[i] <- sum(trees_select[trees_select$classification66 == "G", "Density_expand"], na.rm= TRUE) * 10 
  plot_agg_traits$GeneralistQMD[i] <- sqrt((plot_agg_traits$GeneralistBA[i] / plot_agg_traits$GeneralistDens[i]) / 0.00007854)
  plot_agg_traits$ForestBA[i] <- sum(trees_select[trees_select$classification66 == "F", ]$CSA_BH_expand, na.rm = TRUE) / 10000 * 10
  plot_agg_traits$ForestDens[i] <- sum(trees_select[trees_select$classification66 == "F", "Density_expand"], na.rm= TRUE) * 10
  plot_agg_traits$ForestQMD[i] <- sqrt((plot_agg_traits$ForestBA[i] / plot_agg_traits$ForestDens[i]) / 0.00007854)
  
}


#--------------------------------------------------------------------------------------
# CWM analysis

traits_to_plot <- names(plot_agg_traits)[c(3:10,12,13)]

traits_to_plot <- traits_to_plot[c(1, 2, 6, 3, 4, 5, 7, 8, 9, 10)]

traits_names_clean <- c(expression(paste("Leaf size (cm"^"2", ")")),
                        "Leaf thickness (mm)",
                        expression(atop("Specific leaf area", "(cm"^"2"*" g"^"-1"*")")),
                        "Max height (m)",
                        "Height at 5 cm dia. (m)", 
                        "       Crown ratio\nat 5cm dia. (unitless)",
                        "  Bark thickness\nat 5cm dia. (mm)",
                        "  Bark thickness\nat 8mm dia. (mm)",
                        expression(paste("Wood density (g cm"^"-3",")")),
                        "         Light code\nat 5 cm dia. (unitless)")

#Generate Figure 4

tiff(filename="./plots/Figure_4_community_weighted_traits.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 5, 
     height=8, 
     pointsize=12, 
     res=600)

par(mfrow = c(5,2))
par(oma = c(1,1,3,0), mar = c(4,5,1,1), family = "sans")

for(i in 1:length(traits_to_plot)){

  trait <- traits_to_plot[i]
  
  plot(plot_agg_traits[, eval(substitute(trait))] ~ plot_agg_traits$BA,
       xlab = "",
       ylab = "",
       ylim = c(min(plot_agg_traits[, eval(substitute(trait))]), max(plot_agg_traits[, eval(substitute(trait))]) * 1.05),
       pch = 21, 
       bg = "grey",
       yaxt = "n")
  
  axis(2, las = 2)

  temp <- summary(lm( plot_agg_traits[, eval(substitute(trait))] ~ plot_agg_traits$BA))
  
  if(!(i %in% c(3,4,5,9))){
  text(x = par()$usr[2] - .4*(par()$usr[2] - par()$usr[1]), 
       y = par()$usr[4] - .08*(par()$usr[4] - par()$usr[3]),
       pos = 4, 
       labels = ifelse(coef(temp)[2, 4] > 0.001,
                       paste("p =", round(coef(temp)[2, 4], 3)),
                       "p < 0.001"))
  text(x = par()$usr[2] - .4*(par()$usr[2] - par()$usr[1]), 
       y = par()$usr[4] - .2*(par()$usr[4] - par()$usr[3]),  
       pos = 4, 
       labels = paste("r =", (round(sqrt(temp$r.squared), 2) * sign(coef(temp)[2,1])) ))
  } else{
    text(x = par()$usr[2] - .4*(par()$usr[2] - par()$usr[1]), 
         y = par()$usr[3] + .08*(par()$usr[4] - par()$usr[3]),
         pos = 4, 
         labels = paste("r =", (round(sqrt(temp$r.squared), 2) * sign(coef(temp)[2,1])) ))
    text(x = par()$usr[2] - .4*(par()$usr[2] - par()$usr[1]), 
         y = par()$usr[3] + .2*(par()$usr[4] - par()$usr[3]),  
         pos = 4, 
         labels = ifelse(coef(temp)[2, 4] > 0.001,
                         paste("p =", round(coef(temp)[2, 4], 3)),
                         "p < 0.001"))
  }
  
  mtext(side = 2, text = traits_names_clean[i], line = 2.6, cex = 0.8)
  mtext(side = 3, text = paste0("(", letters[i], ")"), line = 0.7, at = 2, cex = 0.8)
  if(i %in% c(9,10)){
    mtext(side = 1, text = expression(paste("Basal area (m"^"2", " ha"^"-1", ")")), line = 2.8)
  }
  
  abline(coef(temp)[, 1], lty = ifelse(temp$coefficients[2, 4] < (0.05), 1, 2))
  
}

dev.off()

#---------------------------------------------------------------------------------------
# SNC analysis
#---------------------------------------------------------------------------------------
traits <- clean_species_trans[, -c(1, 3, 5, 15, 16)] 
rownames(traits) <- traits$Code
traits <- traits[, -1]

traits_to_plot <- names(traits)

traits_to_plot <- traits_to_plot[c(1, 2, 7, 3, 4, 5, 8, 9, 10, 6)] #reorder

traits_names_clean <- c(expression(paste("Leaf size (cm"^"2", ")")),
                        "Leaf thickness (mm)",
                        expression(paste("Specific leaf area (cm"^"2"," g"^"-1",")")),
                        "Max height (m)",
                        "Height at 5 cm dia. (m)", 
                        "Crown ratio at 5cm dia. (unitless)",
                        "Bark thickness at 5cm dia. (mm)",
                        "Bark thickness at 8mm dia. (mm)",
                        expression(paste("Wood density (g cm"^"-3",")")),
                        "Light code at 5 cm dia. (unitless)")


#initialize species x plot abundance array
abund <- array(0, dim = c(30, length(unique(plot_data$Code))),
               dimnames = list(sites = seq(1:30),
                               species = unique(plot_data$Code)))

#sum BA per species per plot
for(i in 1:30){
  for(j in 1:length(unique(plot_data$Code))){
    abund[i, j] <- sum(plot_data[plot_data$Code == colnames(abund)[j] &
                                   plot_data$P == i, "CSA_30_expand"], na.rm = TRUE)
  }
}


#drop species with no occurrences (or, more accurately, no rows with measured diameter at 30)
n.sites <- apply(abund, 2, function(x){sum(x > 0)})
abund <- abund[, which(n.sites > 0)]

#match species that have trait data and abundances in the plot
abund <- abund[, attr(abund, "dimnames")$species %in% rownames(traits)]
abund <- abund[, order(attr(abund, "dimnames")$species)] #alphabetize
# apply(abund, 2, FUN = sum)

#reduce species to those that have abundance data
traits <- traits[rownames(traits) %in% attr(abund, "dimnames")$species, ]
traits <- traits[order(rownames(traits)), ] #alphabetize

#relative abundance matrix
L <- as.matrix(abund)
for(i in 1:ncol(L)){
  L[, i] <- L[, i] / colSums(L)[i]
}

# calculate weighted SNC, by multiplying relative abundance in the plot
# by site BA
SNC <- as.numeric(as.character(plot_agg_traits$BA %*% L)) 

color_groups <- data.frame(group = as.character(c("F", "G", "S")),
                           col = as.character(c("black", "#1b9e77", "#E69F00")))
fg <- clean_species[match(rownames(traits), clean_species$Code), "FG"]
color_points <- as.character(color_groups[match(fg, color_groups$group), "col"])

# Generate supplemental figure xx
tiff(filename="./plots/species_niche_centroids.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 4.5, 
     height=7.5, 
     pointsize=12, 
     res=600)

par(mfrow = c(5,2))
par(oma = c(1,2.5,1,0), mar = c(4,2,1,3), family = "sans")

for(i in 1:length(traits_to_plot)){
  
  trait <- traits_to_plot[i]
  
  not_na <- !is.na(traits[, eval(substitute(trait))])
  
  plot(SNC[not_na] ~ traits[not_na, eval(substitute(trait))],
       xlab = "",
       ylab = "",
       ylim = c(0, 30),
       pch = 21, 
       bg = color_points[not_na])
  
  temp <- summary(lm( SNC ~ traits[, eval(substitute(trait))]))
  if(i == 1){
    text(x = par()$usr[1], 
         y = par()$usr[4] - .08*(par()$usr[4] - par()$usr[3]),
         pos = 4, 
         labels = ifelse(coef(temp)[2, 4] > 0.001,
                         paste("p =", round(coef(temp)[2, 4], 3)),
                         "p < 0.001"))
    text(x = par()$usr[1], 
         y = par()$usr[4] - .2*(par()$usr[4] - par()$usr[3]),  
         pos = 4, 
         labels = paste("r =", (round(sqrt(temp$r.squared), 2) * sign(coef(temp)[2,1])) ))
  }
  if(i %in% c(2,7,8)){ #plot p-values in different locations for different plots. This is a dumb way to do this sorry future me
    text(x = par()$usr[2] - .4*(par()$usr[2] - par()$usr[1]), 
         y = par()$usr[4] - .08*(par()$usr[4] - par()$usr[3]),
         pos = 4, 
         labels = ifelse(coef(temp)[2, 4] > 0.001,
                         paste("p =", round(coef(temp)[2, 4], 3)),
                         "p < 0.001"))
    text(x = par()$usr[2] - .4*(par()$usr[2] - par()$usr[1]), 
         y = par()$usr[4] - .2*(par()$usr[4] - par()$usr[3]),  
         pos = 4, 
         labels = paste("r =", (round(sqrt(temp$r.squared), 2) * sign(coef(temp)[2,1])) ))
  }
  if(i %in% c(3,4,5,9)){
    text(x = par()$usr[2] - .4*(par()$usr[2] - par()$usr[1]), 
         y = par()$usr[3] + .08*(par()$usr[4] - par()$usr[3]),
         pos = 4, 
         labels = paste("r =", (round(sqrt(temp$r.squared), 2) * sign(coef(temp)[2,1])) ))
    text(x = par()$usr[2] - .4*(par()$usr[2] - par()$usr[1]), 
         y = par()$usr[3] + .2*(par()$usr[4] - par()$usr[3]),  
         pos = 4, 
         labels = ifelse(coef(temp)[2, 4] > 0.001,
                         paste("p =", round(coef(temp)[2, 4], 3)),
                         "p < 0.001"))
  }
  if(i %in% c(6,10)){
    text(x = par()$usr[1], 
         y = par()$usr[3] + .08*(par()$usr[4] - par()$usr[3]),
         pos = 4, 
         labels = paste("r =", (round(sqrt(temp$r.squared), 2) * sign(coef(temp)[2,1])) ))
    text(x = par()$usr[1], 
         y = par()$usr[3] + .2*(par()$usr[4] - par()$usr[3]),  
         pos = 4, 
         labels = ifelse(coef(temp)[2, 4] > 0.001,
                         paste("p =", round(coef(temp)[2, 4], 3)),
                         "p < 0.001"))
  }
  
  mtext(side = 1, text = traits_names_clean[i], line = 2.3, cex = 0.7)
  
  mtext(side = 3, text = paste0("(", letters[i], ")"), 
        at = (par()$usr[1] + .1*(par()$usr[2] - par()$usr[1])), 
        line = 0.6, cex = 0.8)
  
  # if(i %in% c(1,3,5,7,9)){
    mtext(side = 2, outer = TRUE,
          text = expression(paste("SNC for BA (m"^"2", " ha"^"-1", ")")), cex = 1.15, line = 0)
  # }
  #make line dashed or solid depending upon the p-value
  abline(coef(temp)[, 1], lty = ifelse(temp$coefficients[2, 4] < (0.05), 1, 2)) 
  
}

dev.off()


# Plot PCA scores against SNC


scores_test <- sp_scores[, c(1, 2)]
scores_test <- scores_test[rownames(scores_test) %in% rownames(traits), ]
scores_test <- scores_test[order(rownames(scores_test)), ]

SNC <- as.numeric(as.character(plot_agg_traits$BA %*% L))
names(SNC) <- colnames(abund)
SNC <- SNC[names(SNC) %in% rownames(scores_test)]


tiff(filename="./plots/pca_against_SNC.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 3, 
     height=6, 
     pointsize=12, 
     res=600)

par(mfrow = c(2,1))
par(mar = c(5,5,1,1))
par(oma = c(``,0,0,0))

plot(SNC ~ scores_test[, 1],
     xlab = "PC1")
abline(lm(SNC ~ scores_test[, 1]))
summary(lm(SNC ~ scores_test[, 1]))

plot(SNC ~ scores_test[, 2],
     xlab = "PC2")
abline(lm(SNC ~ scores_test[, 2]), lty = 2)
summary(lm(SNC ~ scores_test[, 2]))

dev.off()


# average SNC for each FT

tiff(filename="./plots/SNC by FG.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 3, 
     height=3, 
     pointsize=12, 
     res=600)

SNC <- as.numeric(as.character(plot_agg_traits$BA %*% L))
names(SNC) <- colnames(abund)
SNC <- SNC[names(SNC) %in% rownames(scores_test)]


fg_test <- spList_orig[match(names(SNC), spList_orig$Code), "classification66"]
fg_test <- factor(fg_test, levels=rev(levels(fg_test)))

par(oma = c(0,0,0,0))
par(mar = c(4.1,5.1,2.1,1.1))

boxplot(SNC ~ fg_test, xlab = "Functional type", ylab = expression(paste("SNC for BA (m"^"2", " ha"^"-1", ")")))
aggregate(SNC, by = list(fg_test), FUN = mean)
aggregate(SNC, by = list(fg_test), FUN = function(x){sd(x)/sqrt(length(x))})

dev.off()



#----------------------------------------------------------------------------------
# Figure 1: Functional type ~ BA
#----------------------------------------------------------------------------------
tiff(filename="./plots/Figure_1_change_in_FTs.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 4, 
     height=4, 
     pointsize=12, 
     res=600)

par(mfrow = c(2,1))
par(oma = c(2,2.4,0,1))
par(mar = c(1,2,1,0))
cols <- c("black", "#1b9e77", "#E69F00")
pchs <- c(15,16,17)

plot(NA, ylim = c(0, 30), xlim = c(0, 30), xaxt = "n", yaxt = "n", ylab = "", xlab = "")
legend(x = -1, y = 31, legend = c("Forest", "Generalist", "Savanna"), 
       col = cols[c(1,2,3)], pch = pchs[c(1,2,3)], bty = "n")
mtext(side = 2, text = "Functional type", line = 3.4, adj = 0.5)
mtext(side = 2, text = expression(paste("basal area (m"^"2", " ha"^"-1", ")")), line = 2,
      adj = 0.5)
mtext(side = 3, text = "(a)", at = 0)
axis(1, labels = FALSE)
axis(2, las = 2)

i <- 1
for (type in c("ForestBA", "GeneralistBA", "SavannaBA")){
  points(plot_agg_traits[[type]] ~ BA, data = plot_agg_traits,
         pch = pchs[i], bg = cols[i], col = cols[i])
  newdata <- seq(0, 28, length.out = 100)
  mod <- lm(plot_agg_traits[[type]] ~ poly(BA,2), data = plot_agg_traits)
  # print(coef(summary(mod))[2, 4])
  preds <- predict(mod, newdata = list(BA = newdata))
  lines(preds ~ newdata, col = cols[i], lwd = 2)
  i <- i+1
}


plot(NA, ylim = c(0, 1), xlim = c(0, 30), xaxt = "n", yaxt = "n", ylab = "", xlab = "")
mtext(side = 2, text = "Proportion of BA", line = 2.5)
mtext(side = 1, text = expression(paste("Plot basal area (m"^"2", " ha"^"-1", ")")), line = 1, outer = TRUE,
      adj = 0.5)
mtext(side = 3, text = "(b)", at = 0)
axis(1, labels = TRUE, padj = -.5)
axis(2, las = 2)

i <- 1
for (type in c("ForestBA", "GeneralistBA", "SavannaBA")){
  points(I(plot_agg_traits[[type]]/BA) ~ BA, data = plot_agg_traits,
         pch = pchs[i], bg = cols[i], col = cols[i])
  newdata <- seq(0, 28, length.out = 100)
  mod <- lm(I(plot_agg_traits[[type]]/plot_agg_traits$BA) ~ poly(BA,2), data = plot_agg_traits)
  # print(coef(summary(mod))[2, 4])
  preds <- predict(mod, newdata = list(BA = newdata))
  if(type != "GeneralistBA"){lines(preds ~ newdata, col = cols[i], lwd = 2)}
  i <- i+1
}

dev.off()


tiff(filename="./plots/change_in_FTs_for_talk.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 4, 
     height=4, 
     pointsize=12, 
     res=600)

par(oma = c(2,1.7,0,1))
par(mar = c(1,2,1,0))
cols <- c("black", "#1b9e77", "#E69F00")
pchs <- c(15,16,17)

plot(NA, ylim = c(0, 8000), xlim = c(0, 30), xaxt = "n", ylab = "", xlab = "")
mtext(side = 2, text = expression(paste("Density (", "ha"^"-1", ")")), line = 2)
mtext(side = 1, text = expression(paste("Basal area (m"^"2", " ha"^"-1", ")")), line = 1, outer = TRUE)
axis(1, labels = TRUE, padj = -.5)

i <- 1
for (type in c("ForestDens", "GeneralistDens", "SavannaDens")){
  points(plot_agg_traits[[type]] ~ BA, data = plot_agg_traits, 
         pch = pchs[i], bg = cols[i], col = cols[i])
  newdata <- seq(0, 28, length.out = 100)
  mod <- lm(plot_agg_traits[[type]] ~ poly(BA,2), data = plot_agg_traits)
  print(coef(summary(mod))[2, 4])
  preds <- predict(mod, newdata = list(BA = newdata))
  lines(preds ~ newdata, col = cols[i], lwd = 2)
  i <- i+1
}

dev.off()

tiff(filename="./plots/change_in_FTs_for_talk.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 4, 
     height=4, 
     pointsize=12, 
     res=600)

par(oma = c(2,1.7,0,1))
par(mar = c(1,2,1,0))
cols <- c("#1b9e77", "#E69F00")
pchs <- c(16,17)

plot(NA, ylim = c(0, 8000), xlim = c(0, 30), xaxt = "n", ylab = "", xlab = "")
mtext(side = 2, text = expression(paste("Density (", "ha"^"-1", ")")), line = 2)
mtext(side = 1, text = expression(paste("Basal area (m"^"2", " ha"^"-1", ")")), line = 1, outer = TRUE)
axis(1, labels = TRUE, padj = -.5)


points(plot_agg_traits[["SavannaDens"]] ~ BA, data = plot_agg_traits, 
       pch = pchs[2], bg = cols[2], col = cols[2])
newdata <- seq(0, 28, length.out = 100)
mod <- lm(plot_agg_traits[["SavannaDens"]] ~ poly(BA,2), data = plot_agg_traits)
print(coef(summary(mod))[2, 4])
preds <- predict(mod, newdata = list(BA = newdata))
lines(preds ~ newdata, col = cols[2], lwd = 2)

points(I(plot_agg_traits[["ForestDens"]] + plot_agg_traits[["GeneralistDens"]])~ BA, data = plot_agg_traits, 
       pch = pchs[1], bg = cols[1], col = cols[1])

newdata <- seq(0, 28, length.out = 100)
mod <- lm(I(plot_agg_traits[["ForestDens"]] + plot_agg_traits[["GeneralistDens"]]) ~ poly(BA,2), data = plot_agg_traits)
print(coef(summary(mod))[2, 4])
preds <- predict(mod, newdata = list(BA = newdata))
lines(preds ~ newdata, col = cols[1], lwd = 2)


dev.off()



#------------------------------------------------------------------------------
# Change in FTs over stand BA
#
# TODO: make tables and figures for supplement
#------------------------------------------------------------------------------
summary(lm(plot_agg_traits$SavannaBA ~ plot_agg_traits$BA))
summary(lm(plot_agg_traits$SavannaBA ~ poly(plot_agg_traits$BA,2)))
plot(plot_agg_traits$SavannaBA ~ plot_agg_traits$BA)

summary(lm(plot_agg_traits$SavannaDens ~ plot_agg_traits$Density))
plot(plot_agg_traits$SavannaDens ~ plot_agg_traits$Density)

summary(lm(plot_agg_traits$GeneralistBA ~ plot_agg_traits$BA))
plot(plot_agg_traits$GeneralistBA ~ plot_agg_traits$BA)
summary(lm(plot_agg_traits$GeneralistDens ~ plot_agg_traits$Density))
plot(plot_agg_traits$GeneralistDens ~ plot_agg_traits$Density)

summary(lm(plot_agg_traits$ForestBA ~ plot_agg_traits$BA))
plot(plot_agg_traits$ForestBA ~ plot_agg_traits$BA)
summary(lm(plot_agg_traits$ForestDens ~ plot_agg_traits$Density))
plot(plot_agg_traits$ForestDens ~ plot_agg_traits$Density)


summary(lm(I(plot_agg_traits$SavannaBA / plot_agg_traits$BA)  ~ poly(plot_agg_traits$BA,2)))
summary(lm(I(plot_agg_traits$GeneralistBA / plot_agg_traits$BA)  ~ poly(plot_agg_traits$BA,2)))
summary(lm(I(plot_agg_traits$ForestBA / plot_agg_traits$BA)  ~ poly(plot_agg_traits$BA,2)))

plot(I(plot_agg_traits$SavannaBA / plot_agg_traits$BA)  ~ plot_agg_traits$BA)
plot(I(plot_agg_traits$GeneralistBA / plot_agg_traits$BA)  ~ plot_agg_traits$BA)
plot(I(plot_agg_traits$ForestBA / plot_agg_traits$BA)  ~ plot_agg_traits$BA)

summary(lm(I(plot_agg_traits$SavannaDens / plot_agg_traits$Density)  ~ plot_agg_traits$Density))
summary(lm(I(plot_agg_traits$GeneralistDens / plot_agg_traits$Density)  ~ plot_agg_traits$Density))
summary(lm(I(plot_agg_traits$ForestDens / plot_agg_traits$Density)  ~ plot_agg_traits$Density))

plot(I(plot_agg_traits$SavannaDens / plot_agg_traits$Density)  ~ plot_agg_traits$Density)
plot(I(plot_agg_traits$GeneralistDens / plot_agg_traits$Density)  ~ plot_agg_traits$Density)
plot(I(plot_agg_traits$ForestDens / plot_agg_traits$Density)  ~ plot_agg_traits$Density)
