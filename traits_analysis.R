## Bivariate and ordination analysis

# TODO

# load libraries
library("ggplot2")
library("plyr")
library("vegan")
library("multcomp")
library("multcompView")
library("effects")
library("plotrix")
library("ape")
library("cluster")
library("FD")
library("phytools")
library("Hotelling")

setwd("C:\\Users\\Sam\\Google Drive\\Projects\\Savanna traits")
set.seed(45750762)

#----------------------------------------------------------------------------
# Function to plot elliptical hulls for phylogenetic PCA,
# modified from vegan::ordiellipse. I removed a bunch of functionality,
# so now it only draws elliptical hulls
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
clean_species <- read.csv("./clean data/clean_species2018-11-07.csv")
spList_orig <- read.csv("./raw data/lista_spp_plantas_families_sf_101618.csv")[-c(113, 130, 138, 176), ]
spList_orig$Sp_name <- paste(spList_orig$Genus, " ", spList_orig$Species)
plot_data <- read.csv("./Clean data/plot_data2018-11-07.csv")
plot_data <- subset(plot_data, !is.na(CSA_30))
bark_model <- readRDS("./clean data/bark_model.RDS")
bark_model_fg <- readRDS("./clean data/bark_model_fg.RDS")

clean_species$FG <- factor(clean_species$FG,levels(clean_species$FG)[c(3,2,1)])
clean_species$Leaf_size <- clean_species$Leaf_size/100 
clean_species$Total_leaf_size <- clean_species$Total_leaf_size/100 

clean_species_reduced <- clean_species[complete.cases(clean_species), ]
clean_species_reduced <- clean_species_reduced[clean_species_reduced$Habit %in% c("Tree", "tree", "Treelet", "treelet"), ]

clean_species_reduced_orig <- clean_species_reduced
clean_species_reduced_trans <- clean_species_reduced

#vars to log-trans
log_trans <- c(4, 5, 6, 9, 11, 12, 13)

#log-trans predictor variables
clean_species_reduced_trans[, log_trans] <- log(clean_species_reduced_trans[, log_trans])

clean_species_reduced <- clean_species_reduced_trans


# Bring in phylogeny
tree <- readRDS("phylogeny_code.RDS")

tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, clean_species$Code))

#plot the tree
tiff(filename = "./plots/phylogeny_pruned.tiff",
    width = 7,
    height = 7,
    units = "in",
    pointsize = 8,
    compression = "lzw",
    res = 600)
  
  tree_pruned_plot <- tree_pruned
  tree_pruned_plot$tip.label <- as.character(spList_orig[, "Sp_name"][match(tree_pruned$tip.label, spList_orig$Code)])
  
  sp_codes <- tree_pruned$tip.label
  fg <- spList_orig[match(sp_codes, spList_orig$Code), "FG"]
  fg <- fg[!is.na(fg)]
  
  color_groups <- data.frame(group = as.character(c("F", "G", "S")),
                             col = as.character(c("black", "green", "dark orange")))
  color_points <- as.character(color_groups[match(fg, color_groups$group), "col"])
  
  plot(tree_pruned_plot, type = "fan",
       show.tip.label = TRUE,
       tip.color = color_points)
  legend(x = par()$usr[1], y = par()$usr[4], 
         legend = c("Forest", "Generalist", "Savanna"), 
         col = c("black", "green", "dark orange"), 
         pt.bg = c("black", "green", "dark orange"), 
         pch = 21,
         cex = 1.3,
         bty = "n")
dev.off()

#------------------------------------------------------------------------------
# What proportion of tree species sampled? For results section
#------------------------------------------------------------------------------
#calculate using both clean_species and clean_species_reduced
total_ba <- sum(plot_data$Density_expand, na.rm = TRUE)
sampled_ba <- sum(plot_data[plot_data$Code %in% clean_species$Code, ]$Density_expand, na.rm = TRUE)
sampled_ba/total_ba

plot_ba <- aggregate(plot_data[, c("P", "CSA_BH_expand")], by = list(plot_data$P), FUN = function(x){sum(x, na.rm=TRUE)})
sampled_plot_ba <- aggregate(plot_data[plot_data$Code %in% clean_species$Code, c("P", "CSA_BH_expand")], 
                             by = list(plot_data[plot_data$Code %in% clean_species$Code, ]$P), FUN = function(x){sum(x, na.rm=TRUE)})
sampled_plot_ba$CSA_BH_expand / plot_ba$CSA_BH_expand

#------------------------------------------------------------------------------
# Bivariate analyses
#------------------------------------------------------------------------------
traits_names <- names(clean_species[-c(1:4, 15, 16)])

# Calculate phylogenetic signal, ANOVA, Tukey's post-hoc, phylogenetic ANOVA
log_trans <- c(1, 2, 5, 7, 8, 9) # which traits to test transformed

#fix the tree
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, clean_species_reduced$Code)) # for phylogenetic ANOVA
#reorder species traits to match phylogeny
clean_species_ordered_orig <- clean_species_reduced_orig[match(tree_pruned$tip.label, clean_species_reduced_orig$Code), ] #reorder to match phylogeny

#intialize a table to catch all the results
pvals_table <- data.frame(matrix(ncol = 10, nrow = 9))
#t for tukey, p for phylogenetic comparison
names(pvals_table) <- c("Trait", "tG-S", "tF-S", "tF-G", "pG-S", "pF-S", "pF-G", "ANOVA", "phylANOVA", "lambda")

for (i in 1:length(traits_names)){
  #trait name
  pvals_table[i, 1] <- traits_names[i]
  #fit an ANOVA
  model <- lm(clean_species_reduced[, traits_names[i]] ~ as.factor(clean_species_reduced$FG))
  aov <- aov(model)
  #ANOVA p-value from F-test
  pvals_table[i, 8] <- summary(aov)[[1]][[5]][[1]] #this is dumb; is this really the right way to extract the p-value?!
  
  #calculate post-hoc phylogenetically-corrected differences and add letters to figure
  tuke <- TukeyHSD(aov)
  pvals <- tuke$`as.factor(clean_species_reduced$FG)`[, 4]
  pvals_table[i, c(2:4)] <- pvals
  
  #phylogenetic ANOVA
  phyl_aov_results <- phylANOVA(tree = tree_pruned, x = clean_species_ordered_orig$FG, 
                                y = clean_species_ordered_orig[, traits_names[i]], 
                                nsim=1000, posthoc=TRUE, p.adj="holm")
  pvals_table[i, 9] <- phyl_aov_results$Pf
  
  pvals <- c(phyl_aov_results$Pt[2,1], phyl_aov_results$Pt[3,1], phyl_aov_results$Pt[3,2]) 
  
  pvals_table[i, c(5:7)] <- pvals
  
  
  #calculate the phylogenetic signal in the trait
  sig <- phylosig(tree = tree_pruned, x = clean_species_ordered_orig[, traits_names[i]], 
                  method="lambda", test=TRUE, nsim=1000)
  pvals_table[i, 10] <- sig$lambda
}

write.csv(pvals_table, "./Model output/FG_traits_pvals.csv")

#---------------------------
# Figure 2: trait differences between functional groups

traits_names_clean <- c(expression(paste("Leaf size (cm"^"2", ")")),
                        "Leaf thickness (mm)",
                        "Max height (m)",
                        "Height at 5 cm dia. (m)", 
                        "       Crown ratio\nat 5cm dia. (unitless)",
                        "         Light code\nat 5 cm dia. (unitless)",
                        expression(paste("Specific leaf area (cm"^"2", " g"^"-1", ")")),
                        "  Bark thickness\nat 5cm dia. (mm)",
                        "  Bark thickness\nat 8mm dia. (mm)",
                        expression(paste("Wood density (g cm"^"-3",")")))

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
  noise <- rnorm(nrow(clean_species_reduced), 0, .02)
  
  means <- aggregate(clean_species_reduced_orig[, traits_names[i]], by = list(as.factor(clean_species_reduced$FG)),
                     FUN = mean)
  
  
  if(!(i %in% log_trans)){
    plot(clean_species_reduced[, traits_names[i]] ~ I(as.numeric(clean_species_reduced$FG) + noise),
         xlab = "",
         ylab = "",
         xaxt = 'n',
         xlim = c(0.8, 3.2),
         ylim = c((0.8 * min(clean_species_reduced[, traits_names[i]], na.rm = TRUE, na.rm = TRUE)), 
                  (1.2 * max(clean_species_reduced[, traits_names[i]], na.rm = TRUE))),
         cex.lab = 1.3
       )
  }else{
    plot(exp(clean_species_reduced[, traits_names[i]]) ~ I(as.numeric(clean_species_reduced$FG) + noise),
         xlab = "",
         ylab = "",
         xaxt = 'n',
         xlim = c(0.8, 3.2),
         ylim = c((0.8 * min(exp(clean_species_reduced[, traits_names[i]]), na.rm = TRUE)), 
                  (1.2 * max(exp(clean_species_reduced[, traits_names[i]]), na.rm = TRUE))),
         cex.lab = 1.3
    )
  }
  points(means[, 2] ~ c(1,2,3), pch = 18, cex = 2)
  lines(means[, 2] ~ c(1,2,3), lty = 1, lwd = 2)
  
  mtext(text = levels(clean_species_reduced$FG), side = 1, at = c(1,2,3), line = 1)
  mtext(text = paste0("(", letters[i], ")") , side = 3, at = 0.4, line = .3)
  mtext(text = traits_names_clean[i], side = 2, line = 1.9)
  
  pvals <- as.numeric(pvals_table[i, c(5:7)])
  names(pvals) <- c("G-S", "F-S", "F-G")
  pvals <- ifelse(pvals < (0.05), TRUE, FALSE)
  tuke_letters <- multcompLetters(pvals, reversed=TRUE)$Letters[c(3,1,2)]
  
  #fixes letters sometimes being in the wrong order -- not sure why this happens?
  if(tuke_letters[1] != "a"){
    tuke_letters <- multcompLetters(pvals, reversed=FALSE)$Letters[c(3,1,2)]
  }
  
  #plot the tuke letters for phylogenetic ANOVA
  ylim <- par("usr")[c(3,4)]
  yrange <- ylim[2]-ylim[1]
  text(x = c(1,2,3), y = ylim[2] - yrange/10, labels = tuke_letters, cex = 1.2)
  
}

dev.off()

#------------------------------------------------------------------------------
#PCA on traits
#------------------------------------------------------------------------------

pca_data <- scale(clean_species_reduced[, c("Total_leaf_size", "Leaf_thickness", "Max_height",
                                            "Height_at_5cm", "Crown_ratio", "SLA",
                                            "Bark_at_5cm", "Bark_at_8mm", "Wood_density", "Light_at_5cm")])


rownames(pca_data) <- clean_species_reduced$Code

pca_groups <- clean_species_reduced[, c("FG")]

pca_sp_names <- clean_species_reduced$Code

# pca_data <- pca_data[complete.cases(pca_data), ]
# pca_results <- prcomp(pca_data, scale. = T)

#run the PCA
pca_results <- rda(pca_data, scale = TRUE)
trait_names <- gsub("_", " ", rownames(pca_results$CA$v))
sp_scores <- scores(pca_results,  choices = c(1:10), display = c("sites"))
trait_scores <- scores(pca_results, choices = c(1:10), display = c("species"))

summary(pca_results)

biplot(pca_results, choices = c(1,3))
# ordihull(ord = pca_results, groups = factor(clean_species_reduced$FG))

tiff(filename="./plots/traits_PCA_axes_1_2.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 7, 
     height=7, 
     pointsize=9, 
     res=600)
color_groups <- data.frame(group = as.character(c("F", "G", "S")),
                           col = as.character(c("black", "green", "dark orange")))
color_points <- as.character(color_groups[match(pca_groups, color_groups$group), "col"])

plot(NA,
     xlab = "PC1", ylab = "PC2",
     xlim = c( (min(c(sp_scores[, 1], trait_scores[, 1]) - .5)), 
               (max(c(sp_scores[, 1], trait_scores[, 1])) + 0.4)),
     ylim = c( (min(c(sp_scores[, 2], trait_scores[, 2]) - .5)), 
               (max(c(sp_scores[, 2], trait_scores[, 2])) + 0.4)))
text(y = sp_scores[, 2], x = sp_scores[, 1], labels = pca_sp_names, col = color_points, cex = .7)
# points(y = sp_scores[, 2], x = sp_scores[, 1], col = color_points, bg = color_points, pch = 21, cex = .7)

abline(h = 0)
abline(v = 0)

for(i in 1:3){
  ordiellipse(pca_results, pca_groups, show.groups = color_groups[i, 1], 
              col = as.character(color_groups[i, 2]), choices = c(1,2),
              kind = "sd", conf = 0.90)
}

legend(x = par()$usr[1], y = par()$usr[4], 
       legend = c("Forest", "Generalist", "Savanna"), 
       col = c("black", "green", "dark orange"), 
       pt.bg = c("black", "green", "dark orange"), 
       pch = 21,
       cex = 1.5,
       bty = "n")

segments(x0 = 0, y0 = 0, x1 = trait_scores[, 1], y1 = trait_scores[, 2])
# yoffsets <- c(-0.05, .05, -.07, -.08, .05, -.05, -.05, -.04, 0.05, -.1)
# xoffsets <- c(0, -.05, 0, 0.2, -.05, 0, -.2, -.2, 0, 0.2)
yoffsets <- 0
xoffesets <- 0
text(x = trait_scores[, 1], y = trait_scores[, 2], labels = trait_names, cex = 1.3)

dev.off()

#################
# PCs 1 + 3

tiff(filename="./plots/traits_PCA_axes_1_3.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 7, 
     height=7, 
     pointsize=9, 
     res=600)
color_groups <- data.frame(group = as.character(c("F", "G", "S")),
                           col = as.character(c("black", "green", "dark orange")))
color_points <- as.character(color_groups[match(pca_groups, color_groups$group), "col"])

plot(NA,
     xlab = "PC1", ylab = "PC3",
     xlim = c( (min(c(sp_scores[, 1], trait_scores[, 1]) - .5)), 
               (max(c(sp_scores[, 1], trait_scores[, 1])) + 0.4)),
     ylim = c( (min(c(sp_scores[, 3], trait_scores[, 3]) - .5)), 
               (max(c(sp_scores[, 3], trait_scores[, 3])) + 0.4)))
text(y = sp_scores[, 3], x = sp_scores[, 1], labels = pca_sp_names, col = color_points, cex = .7)
# points(y = sp_scores[, 2], x = sp_scores[, 1], col = color_points, bg = color_points, pch = 21, cex = .7)

abline(h = 0)
abline(v = 0)

for(i in 1:3){
  ordiellipse(pca_results, pca_groups, show.groups = color_groups[i, 1], 
              col = as.character(color_groups[i, 2]), choices = c(1,3),
              kind = "sd", conf = 0.90)
}

legend(x = par()$usr[1], y = par()$usr[4], 
       legend = c("Forest", "Generalist", "Savanna"), 
       col = c("black", "green", "dark orange"), 
       pt.bg = c("black", "green", "dark orange"), 
       pch = 21,
       cex = 1.5,
       bty = "n")

segments(x0 = 0, y0 = 0, x1 = trait_scores[, 1], y1 = trait_scores[, 3])
# yoffsets <- c(-0.05, .05, -.07, -.08, .05, -.05, -.05, -.04, 0.05, -.1)
# xoffsets <- c(0, -.05, 0, 0.2, -.05, 0, -.2, -.2, 0, 0.2)
yoffsets <- 0
xoffesets <- 0
text(x = trait_scores[, 1], y = trait_scores[, 3], labels = trait_names, cex = 1.3)

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

# phy_results_out <- cbind(as.data.frame(phy_sp_scores), clean_species_ordered$FG)
# names(phy_results_out)[11] <- "FG"

for(i in 1:3){
  h2_pca <- hotelling.test(x = pca_results_out[pca_results_out$FG == groups[i,1], c(1:2)], 
                 y = pca_results_out[pca_results_out$FG == groups[i,2], c(1:2)],
                 perm = FALSE)
  hotelling_results[i, 1] <- h2_pca$pval
}

#ANOSIM analysis
anosim_data <- clean_species_reduced[, c("Leaf_size", "Leaf_thickness", "Max_height",
                                         "Height_at_5cm", "Crown_ratio", "SLA",
                                         "Bark_at_5cm", "Bark_at_8mm", "Wood_density", "Light_at_5cm")]

functional_group_anosim <- anosim(anosim_data, grouping = clean_species_reduced$FG, distance = "euclidean")
# functional_group_anosim <- adonis(anosim_data ~ clean_species_reduced$FG, method = "euclidean")
summary(functional_group_anosim)
print(functional_group_anosim)
plot(functional_group_anosim)

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

trait_summary <- aggregate(plot_data[, c(6, 19:30)], by = list(plot_data$FG, plot_data$Life.Form), FUN = mean, na.rm = TRUE)

#--------------------------------------------------------------------------
# Aggregate means for each plot, weighted by CSA of trees
#--------------------------------------------------------------------------
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
  trees_select <- plot_data[plot_data$P == i & !is.na(plot_data$Density_expand) & plot_data$Life.Form %in% c("tree", "treelet"), ]

  plot_agg_traits$Leaf_size[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Leaf_size), ]$Leaf_size, w = trees_select[!is.na(trees_select$Leaf_size), ]$Density_expand, na.rm = TRUE)
  plot_agg_traits$Leaf_thickness[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Leaf_thickness), ]$Leaf_thickness, w = trees_select[!is.na(trees_select$Leaf_thickness), ]$Density_expand, na.rm = TRUE)
  plot_agg_traits$Height_at_5cm[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Height_at_5cm), ]$Height_at_5cm, w = trees_select[!is.na(trees_select$Height_at_5cm), ]$Density_expand, na.rm = TRUE)
  plot_agg_traits$Max_height[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Max_height), ]$Max_height, w = trees_select[!is.na(trees_select$Max_height), ]$Density_expand, na.rm = TRUE)
  plot_agg_traits$Crown_ratio[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Crown_ratio), ]$Crown_ratio, w = trees_select[!is.na(trees_select$Crown_ratio), ]$Density_expand, na.rm = TRUE)
  plot_agg_traits$SLA[i] <- weighted.mean(x = trees_select[!is.na(trees_select$SLA), ]$SLA, w = trees_select[!is.na(trees_select$SLA), ]$Density_expand, na.rm = TRUE)
  plot_agg_traits$Bark_at_5cm[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Bark_at_5cm), ]$Bark_at_5cm, w = trees_select[!is.na(trees_select$Bark_at_5cm), ]$Density_expand, na.rm = TRUE)
  plot_agg_traits$Bark_at_8mm[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Bark_at_8mm), ]$Bark_at_8mm, w = trees_select[!is.na(trees_select$Bark_at_8mm), ]$Density_expand, na.rm = TRUE)
  plot_agg_traits$est_bark_thickness[i] <- weighted.mean(x = trees_select[!is.na(trees_select$est_bark_thickness), ]$est_bark_thickness, w = trees_select[!is.na(trees_select$est_bark_thickness), ]$Density_expand, na.rm = TRUE)
  plot_agg_traits$Wood_density[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Wood_density), ]$Wood_density, w = trees_select[!is.na(trees_select$Wood_density), ]$Density_expand, na.rm = TRUE)
  plot_agg_traits$Light_at_5cm[i] = weighted.mean(x = trees_select[!is.na(trees_select$Light_at_5cm), ]$Light_at_5cm, w = trees_select[!is.na(trees_select$Light_at_5cm), ]$Density_expand, na.rm = TRUE)
  plot_agg_traits$Light_code[i] = weighted.mean(x = trees_select[!is.na(trees_select$Light.code), ]$Light.code, w = trees_select[!is.na(trees_select$Light.code), ]$Density_expand, na.rm = TRUE)
  plot_agg_traits$Height_mean[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Ht), ]$Ht, w = trees_select[!is.na(trees_select$Ht), ]$Density_expand, na.rm = TRUE)

  #some other useful structure variables
  plot_agg_traits$BA[i] <- sum(trees_select$CSA_BH_expand, na.rm = TRUE) / 10000 * 10
  plot_agg_traits$Density[i] <- sum(trees_select$Density_expand, na.rm = TRUE) * 10
  plot_agg_traits$QMD[i] <- sqrt((plot_agg_traits$BA[i] / plot_agg_traits$Density[i]) / 0.00007854)
  plot_agg_traits$SavannaBA[i] <- sum(trees_select[trees_select$FG == "S", ]$CSA_BH_expand, na.rm = TRUE) / 10000 * 10
  plot_agg_traits$SavannaDens[i] <- sum(trees_select[trees_select$FG == "S", "Density_expand"], na.rm= TRUE)  * 10 
  plot_agg_traits$SavannaQMD[i] <- sqrt((plot_agg_traits$SavannaBA[i] / plot_agg_traits$SavannaDens[i]) / 0.00007854)
  plot_agg_traits$GeneralistBA[i] <- sum(trees_select[trees_select$FG == "G", ]$CSA_BH_expand, na.rm = TRUE) / 10000 * 10
  plot_agg_traits$GeneralistDens[i] <- sum(trees_select[trees_select$FG == "G", "Density_expand"], na.rm= TRUE) * 10 
  plot_agg_traits$GeneralistQMD[i] <- sqrt((plot_agg_traits$GeneralistBA[i] / plot_agg_traits$GeneralistDens[i]) / 0.00007854)
  plot_agg_traits$ForestBA[i] <- sum(trees_select[trees_select$FG == "F", ]$CSA_BH_expand, na.rm = TRUE) / 10000 * 10
  plot_agg_traits$ForestDens[i] <- sum(trees_select[trees_select$FG == "F", "Density_expand"], na.rm= TRUE) * 10
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

#Generate Figure xx

tiff(filename="./plots/community_weighted_traits.tiff", 
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
       bg = "grey")

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
  
  mtext(side = 2, text = traits_names_clean[i], line = 2, cex = 0.8)
  mtext(side = 3, text = paste0("(", letters[i], ")"), line = 0.7, at = 2, cex = 0.8)
  if(i %in% c(9,10)){
    mtext(side = 1, text = expression(paste("Basal area (m"^"2", " ha"^"-1", ")")), line = 2.8)
  }
  
  abline(coef(temp)[, 1], lty = ifelse(temp$coefficients[2, 4] < (0.05/10), 1, 2))
  
}

dev.off()

#---------------------------------------------------------------------------------------
# SNC analysis
#---------------------------------------------------------------------------------------
traits <- clean_species_reduced_orig[, -c(1, 3, 4, 15, 16)] 
rownames(traits) <- traits$Code
traits <- traits[, -1]

traits_to_plot <- names(traits)

traits_to_plot <- traits_to_plot[c(1, 2, 7, 3, 4, 5, 8, 9, 10, 6)] #reorder

traits_names_clean <- c(expression(paste("Leaf size (cm"^"2", ")")),
                        "Leaf thickness (mm)",
                        expression(paste("Specific leaf area (cm"^"2"," g"^"-1",")")),
                        "Max height (m)",
                        "Height at 5 cm dia. (m)", 
                        "       Crown ratio at 5cm dia. (unitless)",
                        "  Bark thickness at 5cm dia. (mm)",
                        "  Bark thickness at 8mm dia. (mm)",
                        expression(paste("Wood density (g cm"^"-3",")")),
                        "         Light code at 5 cm dia. (unitless)")


#initialize species x plot abundance array
abund <- array(0, dim = c(30, length(unique(plot_data$Code))),
               dimnames = list(sites = seq(1:30),
                               species = unique(plot_data$Code)))

#weight abundance by density
for(i in 1:30){
  for(j in 1:length(unique(plot_data$Code))){
    abund[i, j] <- sum(plot_data[plot_data$Code == colnames(abund)[j] &
                                   plot_data$P == i, "Density_expand"], na.rm = TRUE)
  }
}

#match species that have trait data and abundances in the plot
abund <- abund[, attr(abund, "dimnames")$species %in% rownames(traits)]
abund <- abund[, order(attr(abund, "dimnames")$species)]
apply(abund, 2, FUN = sum)

traits <- traits[rownames(traits) %in% attr(abund, "dimnames")$species, ]
traits <- traits[order(rownames(traits)), ]

L <- as.matrix(abund)
cs <- colSums(L)

for(i in 1:ncol(L)){
  L[, i] <- L[, i] / colSums(L)[i]
}

SNC <- as.numeric(as.character(plot_agg_traits$BA %*% L)) #calculate weighted SNC

color_groups <- data.frame(group = as.character(c("F", "G", "S")),
                           col = as.character(c("black", "green", "dark orange")))
color_points <- as.character(color_groups[match(pca_groups, color_groups$group), "col"]) #species are the same as in the pca

# Generate supplemental figure xx
tiff(filename="./plots/species_niche_centroids.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 5, 
     height=8, 
     pointsize=12, 
     res=600)

par(mfrow = c(5,2))
par(oma = c(1,1,1,0), mar = c(4,4,1,3.5), family = "sans")

for(i in 1:length(traits_to_plot)){
  
  trait <- traits_to_plot[i]
  
  plot(SNC ~ traits[, eval(substitute(trait))],
       xlab = "",
       ylab = "",
       ylim = c(0, 30),
       pch = 21, 
       bg = color_points)
  
  temp <- summary(lm( SNC ~ traits[, eval(substitute(trait))]))
  
  if(!(i %in% c(3,4,5,9))){ #plot p-values in different locations for different plots
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
  
  mtext(side = 1, text = traits_names_clean[i], line = 2.3, cex = 0.8)
  
  mtext(side = 3, text = paste0("(", letters[i], ")"), 
        at = (par()$usr[1] + .1*(par()$usr[2] - par()$usr[1])), 
        line = 0.6, cex = 0.8)
  
  # if(i %in% c(1,3,5,7,9)){
    mtext(side = 2, outer = TRUE,
          text = expression(paste("Species niche centroid for basal area (m"^"2", " ha"^"-1", ")")), cex = 1.15, line = -2)
  # }
  #make line dashed or solid depending upon the p-value, with bonferroni correction
  abline(coef(temp)[, 1], lty = ifelse(temp$coefficients[2, 4] < (0.05/10), 1, 2)) 
  
}

dev.off()


#----------------------------------------------------------------------------------------
# Functional diversity
#----------------------------------------------------------------------------------------

traits <- clean_species[, -c(1, 3, 4, 15, 16)] 
rownames(traits) <- traits$Code
traits <- traits[, -1]

#initialize species x plot abundance array
abund <- array(0, dim = c(30, length(unique(plot_data$Code))),
                          dimnames = list(sites = seq(1:30),
                                          species = unique(plot_data$Code)))

#weight abundance by basal area
for(i in 1:30){
  for(j in 1:length(unique(plot_data$Code))){
    abund[i, j] <- sum(plot_data[plot_data$Code == colnames(abund)[j] &
                                         plot_data$P == i, "Density_expand"], na.rm = TRUE)
  }
}


abund <- abund[, attr(abund, "dimnames")$species %in% rownames(traits)]
abund <- abund[, order(attr(abund, "dimnames")$species)]

traits <- traits[rownames(traits) %in% attr(abund, "dimnames")$species, ]
traits <- traits[order(rownames(traits)), ]

traits_std <- scale(traits[, sapply(traits, is.numeric)])
  
traits <- traits_std


FD_plots <- dbFD(traits, abund, w.abun = TRUE, corr = "cailliez")

plot(FD_plots$RaoQ ~ plot_agg_traits$BA)
summary(lm(FD_plots$RaoQ ~ plot_agg_traits$BA))
plot(FD_plots$FDis ~ plot_agg_traits$BA)
summary(lm(FD_plots$FDis ~ plot_agg_traits$BA))
plot(FD_plots$FRic ~ plot_agg_traits$BA)
summary(lm(FD_plots$FRic ~ plot_agg_traits$BA))
plot(FD_plots$FDiv ~ plot_agg_traits$BA)
summary(lm(FD_plots$FDiv ~ plot_agg_traits$BA))
plot(FD_plots$FEve ~ plot_agg_traits$BA)
summary(lm(FD_plots$FEve ~ plot_agg_traits$BA))

plot(FD_plots$FRic ~ plot_agg_traits$Light_code)
summary(lm(FD_plots$FRic ~ plot_agg_traits$Light_code))
plot(FD_plots$FDiv ~ plot_agg_traits$Light_code)
summary(lm(FD_plots$FDiv ~ plot_agg_traits$Light_code))
plot(FD_plots$FEve ~ plot_agg_traits$Light_code)
summary(lm(FD_plots$FEve ~ plot_agg_traits$Light_code))
plot(FD_plots$RaoQ ~ plot_agg_traits$Light_code)
summary(lm(FD_plots$RaoQ ~ plot_agg_traits$Light_code))
plot(FD_plots$FDis ~ plot_agg_traits$Light_code)
summary(lm(FD_plots$FDis ~ plot_agg_traits$Light_code))

## Multipanel figure
## sorry that this code is especially WET

tiff(filename="./plots/functional_diversity.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 4, 
     height=5, 
     pointsize=12, 
     res=600)

par(mfrow = c(3, 2))
par(oma = c(2,0,0,0), mar = c(1.8,4,2.5,1), family = "sans")

FRic_lm <- summary((lm(FD_plots$FRic ~ plot_agg_traits$BA)))
  plot(FD_plots$FRic ~ plot_agg_traits$BA, 
       xlab = "",
       ylab = "",
       pch = 21,
       bg = "grey")
  abline(coef(FRic_lm)[, 1])
  text(x = 15, y = par("usr")[3] + ((par("usr")[4] - par("usr")[3]) * .3), 
       pos = 4, 
       labels = paste("r =", (round(sqrt(FRic_lm$r.squared), 2) * sign(coef(FRic_lm)[2,1])) ))
  text(x = 15, y = par("usr")[3] + ((par("usr")[4] - par("usr")[3]) * .2),
       pos = 4, 
       labels = ifelse(coef(FRic_lm)[2, 4] > 0.001,
                       paste("p =", round(coef(FRic_lm)[2, 4], 3)),
                       "p < 0.001"))
  mtext(side = 3, text = "(a)", 
        at = (par()$usr[1] + .1*(par()$usr[2] - par()$usr[1])), 
        line = 0.6, cex = 0.8)
  mtext(side = 2, text = "FRic", line = 2.2)

FDiv_lm <- summary(lm(FD_plots$FDiv ~ plot_agg_traits$BA))
  plot(FD_plots$FDiv ~ plot_agg_traits$BA,
       xlab = "",
       ylab = "",
       pch = 21,
       bg = "grey")
  abline(coef(FDiv_lm)[, 1])
  text(x = 15, y = par("usr")[4] - ((par("usr")[4] - par("usr")[3]) * .2), 
       pos = 4, 
       labels = paste("r =", (round(sqrt(FDiv_lm$r.squared), 2) * sign(coef(FDiv_lm)[2,1])) ))
  text(x = 15, y = par("usr")[4] - ((par("usr")[4] - par("usr")[3]) * .3),
       pos = 4, 
       labels = ifelse(coef(FDiv_lm)[2, 4] > 0.001,
                       paste("p =", round(coef(FDiv_lm)[2, 4], 3)),
                       "p < 0.001"))
  mtext(side = 3, text = "(b)", 
        at = (par()$usr[1] + .1*(par()$usr[2] - par()$usr[1])), 
        line = 0.6, cex = 0.8)
  mtext(side = 2, text = "FDiv", line = 2.2)

FEve_lm <- summary(lm(FD_plots$FEve ~ plot_agg_traits$BA))
  plot(FD_plots$FEve ~ plot_agg_traits$BA,
       xlab = "",
       ylab = "",
       pch = 21,
       bg = "grey")
  abline(coef(FEve_lm)[, 1], lty = 2)
  text(x = 15, y = par("usr")[4] - ((par("usr")[4] - par("usr")[3]) * .2), 
       pos = 4, 
       labels = paste("r =", (round(sqrt(FEve_lm$r.squared), 2) * sign(coef(FEve_lm)[2,1])) ))
  text(x = 15, y = par("usr")[4] - ((par("usr")[4] - par("usr")[3]) * .3),
       pos = 4, 
       labels = ifelse(coef(FEve_lm)[2, 4] > 0.001,
                       paste("p =", round(coef(FEve_lm)[2, 4], 3)),
                       "p < 0.001"))
  mtext(side = 3, text = "(c)", 
        at = (par()$usr[1] + .1*(par()$usr[2] - par()$usr[1])), 
        line = 0.6, cex = 0.8)
  mtext(side = 2, text = "FEve", line = 2.2)

FDis_lm <- summary(lm(FD_plots$FDis ~ plot_agg_traits$BA))
  plot(FD_plots$FDis ~ plot_agg_traits$BA,
       xlab = "",
       ylab = "",
       pch = 21,
       bg = "grey")
  abline(coef(FDis_lm)[, 1])
  text(x = 0, y = par("usr")[3] + ((par("usr")[4] - par("usr")[3]) * .22), 
       pos = 4, 
       labels = paste("r =", (round(sqrt(FRic_lm$r.squared), 2) * sign(coef(FRic_lm)[2,1])) ))
  text(x = 0, y = par("usr")[3] + ((par("usr")[4] - par("usr")[3]) * .12),
       pos = 4, 
       labels = ifelse(coef(FRic_lm)[2, 4] > 0.001,
                       paste("p =", round(coef(FRic_lm)[2, 4], 3)),
                       "p < 0.001"))
  mtext(side = 3, text = "(d)", 
        at = (par()$usr[1] + .1*(par()$usr[2] - par()$usr[1])), 
        line = 0.6, cex = 0.8)
  mtext(side = 2, text = "FDis", line = 2.2)

div <- diversity(abund, index = "invsimpson")
  Simpson_lm <- summary(lm(div ~ plot_agg_traits$BA))
  plot(div ~ plot_agg_traits$BA,
       xlab = "",
       ylab = "",
       pch = 21,
       bg = "grey")
  abline(coef(Simpson_lm)[, 1])
  text(x = 0, y = par("usr")[4] - ((par("usr")[4] - par("usr")[3]) * .2), 
       pos = 4, 
       labels = paste("r =", (round(sqrt(Simpson_lm$r.squared), 2) * sign(coef(Simpson_lm)[2,1])) ))
  text(x = 0, y = par("usr")[4] - ((par("usr")[4] - par("usr")[3]) * .3),
       pos = 4, 
       labels = ifelse(coef(Simpson_lm)[2, 4] > 0.001,
                       paste("p =", round(coef(Simpson_lm)[2, 4], 3)),
                       "p < 0.001"))
  mtext(side = 3, text = "(e)", 
        at = (par()$usr[1] + .1*(par()$usr[2] - par()$usr[1])), 
        line = 0.6, cex = 0.8)
  mtext(side = 2, text = "Inverse Simpson's", line = 2.2)

mtext(side = 1, text = expression(paste("Basal area (m"^"2", " ha"^"-1", ")")), line = 2.8)



dev.off()


#*******************************************************************************************
# Classify species
# analysis for appendix xx
#*******************************************************************************************
habitats <- as.vector(rep("S", times = 30))
habitats[seq(2,24, by = 2)] <- "F"

comm_mat_counts <- array(0, dim = c(30, length(unique(plot_data$Code))),
                         dimnames = list(sites = seq(1:30),
                                         species = unique(plot_data$Code)))

for(i in 1:30){
  for(j in 1:length(unique(plot_data$Code))){
    comm_mat_counts[i, j] <- nrow(plot_data[plot_data$Code == colnames(comm_mat_counts)[j] &
                                              plot_data$P == i,])
  }
}

clam <- clamtest(comm_mat_counts, habitats, coverage.limit = 1, alpha = 0.05, specialization = 2/3)
clam[order(clam$Species), ]

clam_merge <- clam
names(clam_merge)[1] <- "Code"
clam_merge <- join(clam_merge, clean_species[, c(2, 16)], by = c("Code"))
table(clam_merge$FG, clam_merge$Classes)[, c(2, 1, 3, 4)]

#------------------------------------------------------------------------------
# Change in FTs over stand BA
#
# TODO: make tables and figures for supplement
#------------------------------------------------------------------------------
summary(lm(plot_agg_traits$SavannaBA ~ plot_agg_traits$BA))
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


summary(lm(I(plot_agg_traits$SavannaBA / plot_agg_traits$BA)  ~ plot_agg_traits$BA))
summary(lm(I(plot_agg_traits$GeneralistBA / plot_agg_traits$BA)  ~ plot_agg_traits$BA))
summary(lm(I(plot_agg_traits$ForestBA / plot_agg_traits$BA)  ~ plot_agg_traits$BA))

summary(lm(I(plot_agg_traits$SavannaDens / plot_agg_traits$Density)  ~ plot_agg_traits$Density))
summary(lm(I(plot_agg_traits$GeneralistDens / plot_agg_traits$Density)  ~ plot_agg_traits$Density))
summary(lm(I(plot_agg_traits$ForestDens / plot_agg_traits$Density)  ~ plot_agg_traits$Density))


#----------------------------------------------------------------------------------
# Figure 1: Functional type ~ BA
#----------------------------------------------------------------------------------
tiff(filename="./plots/change_in_FTs.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 4, 
     height=4, 
     pointsize=12, 
     res=600)

par(mfrow = c(2,1))
par(oma = c(2,1.7,0,1))
par(mar = c(1,2,1,0))
cols <- c("Green", "Orange", "Black")
pchs <- c(15,16,17)

plot(NA, ylim = c(0, 35), xlim = c(0, 30), xaxt = "n", ylab = "", xlab = "")
legend(x = -1, y = 36, legend = c("Savanna", "Generalist", "Forest"), 
       col = cols[c(2,1,3)], pch = pchs[c(2,1,3)], bty = "n")
mtext(side = 2, text = expression(paste("Basal area (m"^"2", " ha"^"-1", ")")), line = 2)
mtext(side = 3, text = "(a)", at = 0)
axis(1, labels = FALSE)

for (type in c("GeneralistBA", "SavannaBA", "ForestBA")){
  points(plot_agg_traits[[type]] ~ BA, data = plot_agg_traits, 
         pch = pchs[i], bg = cols[i], col = cols[i])
  newdata <- seq(0, 28, length.out = 100)
  mod <- lm(plot_agg_traits[[type]] ~ poly(BA,2), data = plot_agg_traits)
  print(coef(summary(mod))[2, 4])
  preds <- predict(mod, newdata = list(BA = newdata))
  lines(preds ~ newdata, col = cols[i], lwd = 2)
  i <- i+1
}


plot(NA, ylim = c(0, 10000), xlim = c(0, 30), xaxt = "n", ylab = "", xlab = "")
mtext(side = 2, text = expression(paste("Density (", "ha"^"-1", ")")), line = 2)
mtext(side = 1, text = expression(paste("Basal area (m"^"2", " ha"^"-1", ")")), line = 1, outer = TRUE)
mtext(side = 3, text = "(b)", at = 0)
axis(1, labels = TRUE, padj = -.5)

for (type in c("GeneralistDens", "SavannaDens", "ForestDens")){
  points(plot_agg_traits[[type]] ~ BA, data = plot_agg_traits, 
         pch = pchs[i], bg = cols[i], col = cols[i])
  newdata <- seq(0, 28, length.out = 100)
  mod <- lm(plot_agg_traits[[type]] ~ poly(BA,2), data = plot_agg_traits)
  print(coef(summary(mod))[2, 4])
  preds <- predict(mod, newdata = list(BA = newdata))
  lines(preds ~ newdata, col = cols[i], lwd = 2)
  i <- i+1
}


# i <- 1
# plot(NA, ylim = c(0, 25), xlim = c(0, 36), ylab = "", xlab = "")
# mtext(side = 1, text = expression(paste("Basal area (m"^"2", " ha"^"-1", ")")), outer = TRUE)
# mtext(side = 2, text = "Quadratic mean \ndiameter (cm)", line = 2)
# for (type in c("GeneralistQMD", "SavannaQMD", "ForestQMD")){
#   points(plot_agg_traits[[type]] ~ BA, data = plot_agg_traits, 
#          pch = 21, bg = cols[i], col = cols[i])
#   newdata <- seq(0, 36, length.out = 100)
#   mod <- lm(plot_agg_traits[[type]] ~ poly(BA,2), data = plot_agg_traits)
#   preds <- predict(mod, newdata = list(BA = newdata))
#   lines(preds ~ newdata, col = cols[i], lwd = 2)
#   i <- i+1
# }
  
dev.off()
#############################################################################################
# PCA with weighted circle sizes
#############################################################################################
abund <- array(0, dim = c(30, length(unique(plot_data$Code))),
               dimnames = list(sites = seq(1:30),
                               species = unique(plot_data$Code)))

#weight abundance by basal area
for(i in 1:30){
  for(j in 1:length(unique(plot_data$Code))){
    abund[i, j] <- sum(plot_data[plot_data$Code == colnames(abund)[j] &
                                   plot_data$P == i, "Density_expand"], na.rm = TRUE)
  }
}


abund <- abund[, colnames(abund) %in% pca_sp_names]

abund_tot <- abund[c(29, 13, 12), ]

abund_tot <- abund_tot[, order(colnames(abund_tot))]

hull <- chull(y = sp_scores[, 2], x = sp_scores[, 1])
FRic_max <- polyarea(x = sp_scores[hull, 1], y = sp_scores[hull, 2])


tiff(filename="./plots/traits_PCA_axes_1_2_circles.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 7, 
     height=5, 
     pointsize=9, 
     res=600)

par(mfrow = c(1,3))
par(oma = c(3,3,0,0))
color_groups <- data.frame(group = as.character(c("F", "G", "S")),
                           col = as.character(c("black", "green", "dark orange")))

for(i in 1:3){
color_points <- as.character(color_groups[match(pca_groups, color_groups$group), "col"])[abund_tot[i, ] > 0]
sp_scores_red <- sp_scores[abund_tot[i, ] > 0, ]

abund_red <- abund_tot[i, abund_tot[i, ] > 0]

plot(NA,
     xlab = "", ylab = "", xaxt = "n", yaxt = "n",
     xlim = c( (min(c(sp_scores[, 1], trait_scores[, 1]))),
               (max(c(sp_scores[, 1], trait_scores[, 1])))),
     ylim = c( (min(c(sp_scores[, 2], trait_scores[, 2]))),
               (max(c(sp_scores[, 2], trait_scores[, 2])))))

mtext(outer = FALSE, cex = 1.25, text = ifelse(i == 1, expression(paste("Savanna, BA = 0.56 m"^"2", " ha"^"-1")),
                                            ifelse(i == 2, expression(paste("Woodland, BA = 11.9 m"^"2", " ha"^"-1")),
                                                   expression(paste("Forest, BA = 30.4 m"^"2", " ha"^"-1")))))

axis(1, cex = 1.2)
axis(2, cex = 1.2)

# points(y = sp_scores[, 2], x = sp_scores[, 1], col = color_points, bg = color_points, pch = 21, cex = .7)

symbols(y = sp_scores_red[, 2], x = sp_scores_red[, 1], 
        circles=sqrt(abund_red/pi), inches=0.15,
        bg=color_points,
        add = TRUE)
# text(y = sp_scores_red[, 2], x = sp_scores_red[, 1], labels = rownames(sp_scores_red), col = "black", cex = .7)
abline(h = 0)
abline(v = 0)

#Calculate and plot FRic
hull <- chull(y = sp_scores_red[, 2], x = sp_scores_red[, 1])
polygon(sp_scores_red[hull, c(1,2)])
text(x = -0.7, y = -1.1, 
     labels = paste("FRic = ",
                    round(polyarea(x = sp_scores_red[hull, 1], y = sp_scores_red[hull, 2])/FRic_max, 3)),
     cex = 1.4)


#calculate weighted centroid of traits
center <- c(weighted.mean(x = sp_scores_red[, 1], w = abund_red),
            weighted.mean(x = sp_scores_red[, 2], w = abund_red))
points(x = center[1], y = center[2],  cex = 2, pch = 21, bg = "blue")

#Draw weighted segments
segments(x0 = center[1], x1 = sp_scores_red[, 1], 
         y0 = center[2], y1 = sp_scores_red[, 2],
         lwd = (3 * abund_red/ max(abund_red)))

#Calculate and plot FDis
dists <- apply(sp_scores_red[, c(1,2)], 1, FUN = function(x){dist(rbind(x, center))})

FDis <- sum(abund_red * dists) / sum(abund_red)

text(x = -0.7, y = -1.2, labels = paste("FDis = ", round(FDis, digits = 3)),
     cex = 1.4)

if(i == 2){
  
  legend(x = par()$usr[1], y = par()$usr[4], 
         legend = c("Forest", "Generalist", "Savanna", "Centroid"), 
         col = c("black", "green", "dark orange", "dark blue"), 
         pt.bg = c("black", "green", "dark orange", "dark blue"), 
         pch = 21,
         cex = 1.5,
         bty = "n")
}

}


mtext(side = 1, outer = TRUE, text = "PC1", cex = 1.2)
mtext(side = 2, outer = TRUE, text = "PC2", cex = 1.2)
dev.off()
