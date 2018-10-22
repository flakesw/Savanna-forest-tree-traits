## Bivariate and ordination analysis

# TODO

# load libraries
library("ggplot2")
library("plyr")
library("vegan")
library("multcomp")
library("multcompView")
library("pca3d")
library("effects")
library("plotrix")
library("ape")
library("cluster")

setwd("C:\\Users\\Sam\\Google Drive\\Projects\\Savanna traits")
set.seed(45750762)
#----------------------------------------------------------------------------
# Function to plot elliptical hulls for phylogenetic PCA,
# modified from vegan::ordiellipse
#----------------------------------------------------------------------------
phyellipse <- function (ord, groups, display = "sites", 
                        kind = c("sd", "se", "ehull"), 
                        conf, draw = c("lines", "polygon", "none"), 
                        w = weights(ord, display), 
                        col = NULL, alpha = 127, 
                        show.groups, label = FALSE, 
                        border = NULL, lty = NULL, lwd = NULL, ...) 
{
  weights.default <- function(object, ...) NULL
  kind <- match.arg(kind)
  draw <- match.arg(draw)
  pts <- ord$S
  # pts <- scores(pca_results, display = "sites")
  pts <- as.matrix(pts)
  if (ncol(pts) > 2) 
    pts <- pts[, 1:2, drop = FALSE]
  if (ncol(pts) < 2) 
    stop("ordiellipse needs two dimensions")
  w <- eval(w)
  if (length(w) == 1) 
    w <- rep(1, nrow(pts))
  if (is.null(w)) 
    w <- rep(1, nrow(pts))
  if (alpha < 1) 
    alpha <- round(alpha * 255)
  if (draw == "polygon" && !is.null(col)) 
    col <- rgb(t(col2rgb(col)), alpha = alpha, maxColorValue = 255)
  if (!missing(show.groups)) {
    take <- groups %in% show.groups
    pts <- pts[take, , drop = FALSE]
    groups <- groups[take]
    w <- w[take]
  }
  out <- seq(along = groups)
  inds <- names(table(groups))
  col.new <- border.new <- lty.new <- lwd.new <- NULL
  for (arg in c("col", "border", "lty", "lwd")) {
    tmp <- mget(arg, ifnotfound = list(NULL))[[1]]
    if (is.null(tmp)) 
      tmp <- ifelse(suppressWarnings(is.null(par(arg))), 
                    par("fg"), par(arg))
    if (length(inds) != length(tmp)) {
      tmp <- rep_len(tmp, length(inds))
    }
    assign(paste(arg, ".new", sep = ""), tmp)
  }
  if (is.null(col) && draw == "polygon") 
    col.new <- rep_len("transparent", length(inds))
  else if (is.null(col) && draw == "lines") 
    col.new <- rep_len(par("fg"), length(inds))
  res <- list()
  if (label) {
    cntrs <- matrix(NA, nrow = length(inds), ncol = 2)
    rownames(cntrs) <- inds
  }
  kk <- complete.cases(pts) & !is.na(groups)
  for (is in inds) {
    gr <- out[groups == is & kk]
    if (length(gr)) {
      X <- pts[gr, , drop = FALSE]
      W <- w[gr]
      if (kind == "ehull") {
        tmp <- ellipsoidhull(X)
        mat <- list(cov = tmp$cov, center = tmp$loc, 
                    n.obs = nrow(X))
      }
      else mat <- cov.wt(X, W)
      if (mat$n.obs == 1) 
        mat$cov[] <- 0
      if (kind == "se") 
        mat$cov <- mat$cov * sum(mat$wt^2)
      if (kind == "ehull") 
        t <- sqrt(tmp$d2)
      else {
        if (missing(conf)) 
          t <- 1
        else t <- sqrt(qchisq(conf, 2))
      }
      if (mat$n.obs > 1) 
        xy <- vegan:::veganCovEllipse(mat$cov, mat$center, t)
      else xy <- X
      if (draw == "lines") 
        vegan:::ordiArgAbsorber(xy, FUN = lines, col = if (is.null(col)) 
          par("fg")
          else col.new[match(is, inds)], lty = lty.new[match(is, 
                                                             inds)], lwd = lwd.new[match(is, inds)])
      else if (draw == "polygon") 
        vegan:::ordiArgAbsorber(xy[, 1], xy[, 2], col = col.new[match(is, 
                                                                      inds)], border = border.new[match(is, inds)], 
                                lty = lty.new[match(is, inds)], lwd = lwd.new[match(is, 
                                                                                    inds)], FUN = polygon)
      if (label && draw != "none") {
        cntrs[is, ] <- mat$center
      }
      mat$scale <- t
      res[[is]] <- mat
    }
  }
  if (label && draw != "none") {
    if (draw == "lines") 
      ordiArgAbsorber(cntrs[, 1], cntrs[, 2], labels = rownames(cntrs), 
                      col = col.new, FUN = text, ...)
    else ordiArgAbsorber(cntrs, col = NULL, FUN = ordilabel, 
                         ...)
  }
  class(res) <- "ordiellipse"
  invisible(res)
}


#----------------------------------------------------------------------------
# Import and wrangle data
#----------------------------------------------------------------------------

#read data
clean_species <- read.csv("./clean data/clean_species2018-10-08.csv")
plot_data <- read.csv("./Clean data/plot_data2018-10-08.csv")
plot_data <- subset(plot_data, !is.na(CSA_30))
bark_model <- readRDS("./clean data/bark_model.RDS")

clean_species$FG <- factor(clean_species$FG,levels(clean_species$FG)[c(3,2,1)])
clean_species$Leaf_size <- clean_species$Leaf_size/100 
clean_species$Total_leaf_size <- clean_species$Total_leaf_size/100 

clean_species_reduced <- clean_species[complete.cases(clean_species), ]
clean_species_reduced <- clean_species_reduced[clean_species_reduced$Habit %in% c("Tree", "tree", "Treelet", "treelet"), ]

#### Transformations
#exploration
hist(clean_species_reduced$Leaf_size)
hist(log(clean_species_reduced$Leaf_size)) #yes
hist(clean_species_reduced$Total_leaf_size)
hist(log(clean_species_reduced$Total_leaf_size)) #yes
hist(clean_species_reduced$Leaf_thickness)
hist(log(clean_species_reduced$Leaf_thickness)) #yes
hist(clean_species_reduced$Max_height)#no
hist(clean_species_reduced$Height_at_5cm)#no
hist(clean_species_reduced$Crown_ratio)#no
hist(clean_species_reduced$Light_at_5cm)#no
hist(clean_species_reduced$SLA)
hist(log(clean_species_reduced$SLA))#yes
hist(clean_species_reduced$Bark_at_5cm)
hist(log(clean_species_reduced$Bark_at_5cm)) #yes
hist(clean_species_reduced$Bark_at_8mm)
hist(log(clean_species_reduced$Bark_at_8mm)) #yes
hist(clean_species_reduced$Wood_density)#no


clean_species_reduced_trans <- clean_species_reduced

#vars to log-trans
log_trans <- c(4, 5, 6, 9, 11, 12, 13)

#log-trans predictor variables
clean_species_reduced_trans[, log_trans] <- log(clean_species_reduced_trans[, log_trans])

clean_species_reduced <- clean_species_reduced_trans


# Bring in phylogeny
tree <- readRDS("phylogeny_code.RDS")

tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, clean_species$Code))

# plot(tree_pruned, type = "fan")


#------------------------------------------------------------------------------
# Bivariate analyses
#------------------------------------------------------------------------------
# Figure 1: trait differences between functional groups
log_trans <- c(1, 2, 5, 7, 8)

traits_names <- names(clean_species[-c(1:4, 13, 15, 16)])
traits_names_clean <- c(expression(paste("Leaf size (cm"^"2", ")")),
                        "Leaf thickness (mm)",
                        "Max height (m)",
                        "Height at 5 cm dia. (m)", 
                        "       Crown ratio\nat 5cm dia. (unitless)",
                        "         Light code\nat 5 cm dia. (unitless)",
                        expression(paste("Specific leaf area (cm"^"2", " g"^"-1", ")")),
                        "  Bark thickness\nat 5cm dia. (mm)",
                        expression(paste("Wood density (g cm"^"-3",")")))

tiff(filename="./plots/traits_by_FG.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 5, 
     height=6, 
     pointsize=9, 
     res=600)

par(mfrow = c(3, 3),
    mar = c(2,7,2,2),
    oma = c(1, 1, 0, 0))

for (i in 1:length(traits_names)){
  model <- lm(clean_species_reduced[, traits_names[i]] ~ as.factor(clean_species_reduced$FG))
  noise <- rnorm(nrow(clean_species_reduced), 0, .02)
  
  
  if(!(i %in% log_trans)){
    plot(clean_species_reduced[, traits_names[i]] ~ I(as.numeric(clean_species_reduced$FG) + noise),
         xlab = "",
         ylab = traits_names_clean[i],
         xaxt = 'n',
         xlim = c(0.8, 3.2),
         ylim = c((0.8 * min(clean_species_reduced[, traits_names[i]], na.rm = TRUE, na.rm = TRUE)), 
                  (1.2 * max(clean_species_reduced[, traits_names[i]], na.rm = TRUE))),
         cex.lab = 1.3
       )
  }else{
    plot(exp(clean_species_reduced[, traits_names[i]]) ~ I(as.numeric(clean_species_reduced$FG) + noise),
         xlab = "",
         ylab = traits_names_clean[i],
         xaxt = 'n',
         xlim = c(0.8, 3.2),
         ylim = c((0.8 * min(exp(clean_species_reduced[, traits_names[i]]), na.rm = TRUE)), 
                  (1.2 * max(exp(clean_species_reduced[, traits_names[i]]), na.rm = TRUE))),
         cex.lab = 1.3
    )
  }
  
  mtext(text = levels(clean_species_reduced$FG), side = 1, at = c(1,2,3), line = 1)
  mtext(text = paste0("(", letters[i], ")") , side = 3, at = 0.4, line = .3)
  
  
  #calculte post-hoc differences and add letters to figure
  tuke <- TukeyHSD(aov(model))
  pvals <- tuke$`as.factor(clean_species_reduced$FG)`[, 4]
  
  tuke_letters <- multcompLetters(pvals, reversed=TRUE)$Letters[c(3,1,2)]
  
  ylim <- par("usr")[c(3,4)]
  yrange <- ylim[2]-ylim[1]
  
  text(x = c(1,2,3), y = ylim[2] - yrange/10, labels = tuke_letters, cex = 1.2)
}

dev.off()

#------------------------------------------------------------------------------
# Phylogenetic ANOVA
#------------------------------------------------------------------------------

tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, clean_species_reduced$Code))

length(tree_pruned$tip.label)

clean_species_ordered <- clean_species_reduced[match(tree_pruned$tip.label, clean_species_reduced$Code), ]

phylANOVA(tree = tree_pruned, x = clean_species_ordered$FG, y = clean_species_ordered$Total_leaf_size, nsim=1000, posthoc=TRUE, p.adj="holm")
phylosig(tree = tree_pruned, x = clean_species_ordered$Total_leaf_size, method="lambda", test=FALSE, nsim=1000, se=NULL, start=NULL,
         control=list())
phylANOVA(tree = tree_pruned, x = clean_species_ordered$FG, y = clean_species_ordered$Leaf_thickness, nsim=1000, posthoc=TRUE, p.adj="holm")
phylosig(tree = tree_pruned, x = clean_species_ordered$Leaf_thickness, method="lambda", test=FALSE, nsim=1000, se=NULL, start=NULL,
         control=list())
phylANOVA(tree = tree_pruned, x = clean_species_ordered$FG, y = clean_species_ordered$Max_height, nsim=1000, posthoc=TRUE, p.adj="holm")
phylosig(tree = tree_pruned, x = clean_species_ordered$Max_height, method="lambda", test=FALSE, nsim=1000, se=NULL, start=NULL,
         control=list())
phylANOVA(tree = tree_pruned, x = clean_species_ordered$FG, y = clean_species_ordered$Height_at_5cm, nsim=1000, posthoc=TRUE, p.adj="holm")
phylosig(tree = tree_pruned, x = clean_species_ordered$Height_at_5cm, method="lambda", test=FALSE, nsim=1000, se=NULL, start=NULL,
         control=list())
phylANOVA(tree = tree_pruned, x = clean_species_ordered$FG, y = clean_species_ordered$Crown_ratio, nsim=1000, posthoc=TRUE, p.adj="holm")
phylosig(tree = tree_pruned, x = clean_species_ordered$Crown_ratio, method="lambda", test=FALSE, nsim=1000, se=NULL, start=NULL,
         control=list())
phylANOVA(tree = tree_pruned, x = clean_species_ordered$FG, y = clean_species_ordered$Light_at_5cm, nsim=1000, posthoc=TRUE, p.adj="holm")
phylosig(tree = tree_pruned, x = clean_species_ordered$Light_at_5cm, method="lambda", test=FALSE, nsim=1000, se=NULL, start=NULL,
         control=list())
phylANOVA(tree = tree_pruned, x = clean_species_ordered$FG, y = clean_species_ordered$SLA, nsim=1000, posthoc=TRUE, p.adj="holm")
phylosig(tree = tree_pruned, x = clean_species_ordered$SLA, method="lambda", test=FALSE, nsim=1000, se=NULL, start=NULL,
         control=list())
phylANOVA(tree = tree_pruned, x = clean_species_ordered$FG, y = clean_species_ordered$Bark_at_5cm, nsim=1000, posthoc=TRUE, p.adj="holm")
phylosig(tree = tree_pruned, x = clean_species_ordered$Bark_at_5cm, method="lambda", test=FALSE, nsim=1000, se=NULL, start=NULL,
         control=list())
phylANOVA(tree = tree_pruned, x = clean_species_ordered$FG, y = clean_species_ordered$Bark_at_8mm, nsim=1000, posthoc=TRUE, p.adj="holm")
phylosig(tree = tree_pruned, x = clean_species_ordered$Bark_at_8mm, method="lambda", test=FALSE, nsim=1000, se=NULL, start=NULL,
         control=list())
phylANOVA(tree = tree_pruned, x = clean_species_ordered$FG, y = clean_species_ordered$Wood_density, nsim=1000, posthoc=TRUE, p.adj="holm")
phylosig(tree = tree_pruned, x = clean_species_ordered$Wood_density, method="lambda", test=FALSE, nsim=1000, se=NULL, start=NULL,
         control=list())


pvals_test <- c(test$Pt[2,1], test$Pt[3,1], test$Pt[3,2]) 
names(pvals_test) <- c("G-S", "F-S", "F-G")
multcompLetters(pvals_test, reversed=TRUE)$Letters[c(3,1,2)]


phylosig(tree = tree_pruned, x = clean_species_ordered$Wood_density, method="lambda", test=FALSE, nsim=1000, se=NULL, start=NULL,
         control=list())
#------------------------------------------------------------------------------
# What proportion of tree species sampled? For results section
#------------------------------------------------------------------------------
#calculate using both clean_species and clean_species_reduced
total_ba <- sum(plot_data$CSA_BH_expand, na.rm = TRUE)/10000
sampled_ba <- sum(plot_data[plot_data$Code %in% clean_species$Code, ]$CSA_BH_expand, na.rm = TRUE)/10000
sampled_ba/total_ba

plot_ba <- aggregate(plot_data[, c("P", "CSA_BH_expand")], by = list(plot_data$P), FUN = function(x){sum(x, na.rm=TRUE)})
sampled_plot_ba <- aggregate(plot_data[plot_data$Code %in% clean_species$Code, c("P", "CSA_BH_expand")], 
                             by = list(plot_data[plot_data$Code %in% clean_species$Code, ]$P), FUN = function(x){sum(x, na.rm=TRUE)})
sampled_plot_ba$CSA_BH_expand / plot_ba$CSA_BH_expand

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

# summary(pca_results)

biplot(pca_results, choices = c(1,2))
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

legend(x = -1.7, y = 1.65, 
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


#-----------------------------------------------------------------------------
# Phylogenetic PCA
#-----------------------------------------------------------------------------
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, clean_species_reduced$Code))

clean_species_ordered <- clean_species_reduced[match(tree_pruned$tip.label, clean_species_reduced$Code), ]

phy_pca_data <- scale(clean_species_ordered[, c("Total_leaf_size", "Leaf_thickness", "Max_height",
                                            "Height_at_5cm", "Crown_ratio", "SLA",
                                            "Bark_at_5cm", "Bark_at_8mm", "Wood_density", "Light_at_5cm")])
rownames(phy_pca_data) <- clean_species_ordered$Code
phy_pca_groups <- clean_species_ordered[, c("FG")]

phy_pca_sp_names <- clean_species_ordered$Code



#run the PCA
phy_pca_results <- phyl.pca(tree = tree_pruned, Y = phy_pca_data[match(tree_pruned$tip.label, rownames(phy_pca_data)), ], method="BM", mode="cov")
trait_names <- gsub("_", " ", rownames(phy_pca_results$Evec))
sp_scores <- phy_pca_results$S
trait_scores <- phy_pca_results$L * 3

summary(phy_pca_results)
biplot(phy_pca_results)




tiff(filename="./plots/traits_phylo_PCA_axes_1_2.tiff",
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
color_points <- as.character(color_groups[match(phy_pca_groups, color_groups$group), "col"])

plot(NA,
     xlab = "PC1", ylab = "PC2",
     xlim = c( (min(c(sp_scores[, 1], trait_scores[, 1]) - .5)),
               (max(c(sp_scores[, 1], trait_scores[, 1])) + 0.4)),
     ylim = c( (min(c(sp_scores[, 2], trait_scores[, 2]) - .5)),
               (max(c(sp_scores[, 2], trait_scores[, 2])) + 0.4)))
text(y = sp_scores[, 2], x = sp_scores[, 1], labels = phy_pca_sp_names, col = color_points, cex = .7)
# points(y = sp_scores[, 2], x = sp_scores[, 1], col = color_points, bg = color_points, pch = 21, cex = .7)

abline(h = 0)
abline(v = 0)

for(i in 1:3){
  phyellipse(ord = phy_pca_results, groups = factor(clean_species_ordered$FG),
             show.groups = color_groups[i, 1],
             col = as.character(color_groups[i, 2]),
             kind = "ehull")
}

legend(x = -1.7, y = 1.65,
       legend = c("Forest", "Generalist", "Savanna"),
       col = c("black", "green", "dark orange"),
       pt.bg = c("black", "green", "dark orange"),
       pch = 21,
       cex = 1.5,
       bty = "n")

segments(x0 = 0, y0 = 0, x1 = trait_scores[, 1], y1 = trait_scores[, 2])
yoffsets <- c(-0.05, .05, -.07, -.08, .05, -.05, -.05, -.04, 0.05, -.1)
xoffsets <- c(0, -.05, 0, 0.2, -.05, 0, -.2, -.2, 0, 0.2)
text(x = trait_scores[, 1], y = trait_scores[, 2], labels = trait_names, cex = 1.3)

dev.off()

anosim(phy_pca_results$S, grouping = clean_species_reduced$FG, distance = "euclidean")



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
library("Hotelling")
pca_results_out <- cbind(as.data.frame(sp_scores), clean_species_ordered$FG)
names(pca_results_out)[11] <- "FG"
test <- hotelling.test(x = sp_scores[clean_species_ordered$FG == "F", c(1,2)], 
               y = sp_scores[clean_species_ordered$FG == "G", c(1,2)],
               perm = TRUE,
               B = 10000, progBar = FALSE)
summary(test)
test$stats
test$pval
plot(test)


anosim_data <- scale(clean_species_reduced[, c("Leaf_size", "Leaf_thickness", "Max_height",
                                               "Height_at_5cm", "Crown_ratio", "SLA",
                                               "Bark_at_5cm", "Bark_at_8mm", "Wood_density", "Light_at_5cm")])

functional_group_anosim <- anosim(anosim_data, grouping = clean_species_reduced$FG, distance = "euclidean")
functional_group_anosim <- anosim(scores(pca_results,  choices = c(1:10), display = c("sites")), grouping = clean_species_reduced$FG, distance = "euclidean")

summary(functional_group_anosim)
plot(functional_group_anosim)


adonis_data <- scale(clean_species_reduced[, c("Leaf_size", "Leaf_thickness", "Max_height",
                                               "Height_at_5cm", "Crown_ratio", "SLA",
                                               "Bark_at_5cm", "Bark_at_8mm", "Wood_density", "Light_at_5cm")])

functional_group_adonis <- adonis(adonis_data ~ clean_species_reduced$FG, method = "euclidean")
print(functional_group_adonis)


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
    plot_data$est_bark_thickness[i] <- predict(bark_model, newdata = list(D30_eff = plot_data$D30_eff[i], Code = plot_data$Code[i]))
    plot_data$Wood_density[i] <- clean_species[clean_species$Code == as.character(plot_data$Code[i]), "Wood_density"]
    plot_data$Light_at_5cm[i] <- clean_species[clean_species$Code == as.character(plot_data$Code[i]), "Light_at_5cm"]
  } 
  
}

trait_summary <- aggregate(plot_data[, c(6, 19:30)], by = list(plot_data$FG, plot_data$Life.Form), FUN = mean, na.rm = TRUE)


for (i in 1:nrow(plot_data)){
  #fill in non-measured trees with the mean value for their funcional type and habit
  if(!was_measured[i]){
  plot_data$Leaf_size[i] <- trait_summary[trait_summary$Group.1 == as.character(plot_data$FG[i]) &
                                            trait_summary$Group.2 == as.character(plot_data$Life.Form[i]), "Leaf_size"]
  plot_data$Leaf_thickness[i] <- trait_summary[trait_summary$Group.1 == as.character(plot_data$FG[i]) &
                                                 trait_summary$Group.2 == as.character(plot_data$Life.Form[i]), "Leaf_thickness"]
  plot_data$Height_at_5cm[i] <- trait_summary[trait_summary$Group.1 == as.character(plot_data$FG[i]) &
                                                trait_summary$Group.2 == as.character(plot_data$Life.Form[i]), "Height_at_5cm"]
  plot_data$Max_height[i] <- trait_summary[trait_summary$Group.1 == as.character(plot_data$FG[i]) &
                                             trait_summary$Group.2 == as.character(plot_data$Life.Form[i]), "Max_height"]
  plot_data$Crown_ratio[i] <- trait_summary[trait_summary$Group.1 == as.character(plot_data$FG[i]) &
                                                     trait_summary$Group.2 == as.character(plot_data$Life.Form[i]), "Crown_ratio"]
  plot_data$SLA[i] <- trait_summary[trait_summary$Group.1 == as.character(plot_data$FG[i]) &
                                      trait_summary$Group.2 == as.character(plot_data$Life.Form[i]), "SLA"]
  plot_data$Bark_at_5cm[i] <- trait_summary[trait_summary$Group.1 == as.character(plot_data$FG[i]) &
                                              trait_summary$Group.2 == as.character(plot_data$Life.Form[i]), "Bark_at_5cm"]
  plot_data$Bark_at_8mm[i] <- trait_summary[trait_summary$Group.1 == as.character(plot_data$FG[i]) &
                                              trait_summary$Group.2 == as.character(plot_data$Life.Form[i]), "Bark_at_8mm"]

  plot_data$est_bark_thickness[i] <- trait_summary

  plot_data$Wood_density[i] <- trait_summary[trait_summary$Group.1 == as.character(plot_data$FG[i]) &
                                              trait_summary$Group.2 == as.character(plot_data$Life.Form[i]), "Wood_density"]
  plot_data$Light_at_5cm[i] <- trait_summary[trait_summary$Group.1 == as.character(plot_data$FG[i]) &
                                              trait_summary$Group.2 == as.character(plot_data$Life.Form[i]), "Light_at_5cm"]
  }
}

#--------------------------------------------------------------------------

# Aggregate means for each plot, weighted by CSA of trees

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
                              PercentSavannaBA = numeric(30),
                              PercentSavannaDens = numeric(30)
                              )


for ( i in 1:30){
  trees_select <- plot_data[plot_data$P == i & !is.na(plot_data$CSA_BH_expand) & plot_data$Life.Form %in% c("tree", "treelet"), ]
  
  plot_agg_traits$Leaf_size[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Leaf_size), ]$Leaf_size, w = trees_select[!is.na(trees_select$Leaf_size), ]$CSA_BH_expand, na.rm = TRUE)
  plot_agg_traits$Leaf_thickness[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Leaf_thickness), ]$Leaf_thickness, w = trees_select[!is.na(trees_select$Leaf_thickness), ]$CSA_BH_expand, na.rm = TRUE)
  plot_agg_traits$Height_at_5cm[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Height_at_5cm), ]$Height_at_5cm, w = trees_select[!is.na(trees_select$Height_at_5cm), ]$CSA_BH_expand, na.rm = TRUE)
  plot_agg_traits$Max_height[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Max_height), ]$Max_height, w = trees_select[!is.na(trees_select$Max_height), ]$CSA_BH_expand, na.rm = TRUE)
  plot_agg_traits$Crown_ratio[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Crown_ratio), ]$Crown_ratio, w = trees_select[!is.na(trees_select$Crown_ratio), ]$CSA_BH_expand, na.rm = TRUE)
  plot_agg_traits$SLA[i] <- weighted.mean(x = trees_select[!is.na(trees_select$SLA), ]$SLA, w = trees_select[!is.na(trees_select$SLA), ]$CSA_BH_expand, na.rm = TRUE)
  plot_agg_traits$Bark_at_5cm[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Bark_at_5cm), ]$Bark_at_5cm, w = trees_select[!is.na(trees_select$Bark_at_5cm), ]$CSA_BH_expand, na.rm = TRUE)
  plot_agg_traits$Bark_at_8mm[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Bark_at_8mm), ]$Bark_at_8mm, w = trees_select[!is.na(trees_select$Bark_at_8mm), ]$CSA_BH_expand, na.rm = TRUE)
  plot_agg_traits$est_bark_thickness[i] <- weighted.mean(x = trees_select[!is.na(trees_select$est_bark_thickness), ]$est_bark_thickness, w = trees_select[!is.na(trees_select$est_bark_thickness), ]$CSA_BH_expand, na.rm = TRUE)
  plot_agg_traits$Wood_density[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Wood_density), ]$Wood_density, w = trees_select[!is.na(trees_select$Wood_density), ]$CSA_BH_expand, na.rm = TRUE)
  plot_agg_traits$Light_at_5cm[i] = weighted.mean(x = trees_select[!is.na(trees_select$Light_at_5cm), ]$Light_at_5cm, w = trees_select[!is.na(trees_select$Light_at_5cm), ]$CSA_BH_expand, na.rm = TRUE)
  plot_agg_traits$Light_code[i] = weighted.mean(x = trees_select[!is.na(trees_select$Light.code), ]$Light.code, w = trees_select[!is.na(trees_select$Light.code), ]$CSA_BH_expand, na.rm = TRUE)
  plot_agg_traits$Height_mean[i] <- weighted.mean(x = trees_select[!is.na(trees_select$Ht), ]$Ht, w = trees_select[!is.na(trees_select$Ht), ]$CSA_BH_expand, na.rm = TRUE)

  #some other useful variables
  plot_agg_traits$BA[i] <- sum(trees_select$CSA_BH_expand, na.rm = TRUE) / 10000 * 10
  plot_agg_traits$Density[i] <- sum(trees_select$Density_expand, na.rm = TRUE)
  plot_agg_traits$QMD[i] <- sqrt((plot_agg_traits$BA[i] / plot_agg_traits$Density[i]) / 0.00007854)
  ba_savanna <- sum(trees_select[trees_select$FG == "S", ]$CSA_BH_expand, na.rm = TRUE) / 10000 * 10
  plot_agg_traits$PercentSavannaBA[i] <- ba_savanna / plot_agg_traits$BA[i]
  plot_agg_traits$PercentSavannaDens[i] <- sum(trees_select[trees_select$FG == "S", "Density_expand"], na.rm= TRUE) / plot_agg_traits$Density[i]
}


tiff(filename="./plots/community_weighted_traits.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 7, 
     height=7, 
     pointsize=12, 
     res=600)



traits_to_plot <- names(plot_agg_traits)[3:16]

par(mfrow = c(4, 4))
par(oma = c(2,2,0,0), mar = c(4,4,1,1), family = "sans")

for(i in 1:length(traits_to_plot)){
  
  # there has got to be a better way to substitute the trait name than this,
  # maybe one that keeps the data = plot_agg_traits structure. But I can't
  # find it
  trait <- traits_to_plot[i]
  
  plot(plot_agg_traits[, eval(substitute(trait))] ~ plot_agg_traits$BA,
       xlab = "BA",
       ylab = substitute(trait))
  
  summary(lm(plot_agg_traits[, eval(substitute(trait))] ~ plot_agg_traits$BA))
  
  
  temp <- summary(lm( plot_agg_traits[, eval(substitute(trait))] ~ plot_agg_traits$BA))
  
  text(x = par()$usr[2] - .1*(par()$usr[2] - par()$usr[1]), 
       y = par()$usr[4] - .05*(par()$usr[4] - par()$usr[3]), 
       labels = paste("r =", round(sqrt(temp$r.squared), 2)))
  text(x = par()$usr[2] - .1*(par()$usr[2] - par()$usr[1]), 
       y = par()$usr[4] - .13*(par()$usr[4] - par()$usr[3]),  
       labels = paste("p =", round(coef(temp)[2, 4], 6)))
  
  abline(coef(temp)[, 1], lty = ifelse(temp$coefficients[2, 4] < 0.05, 1, 2))
  
  
}

dev.off()


plot(Leaf_thickness ~ BA, data = plot_agg_traits)
plot(Height_at_5cm ~ BA, data = plot_agg_traits)
plot(Max_height ~ BA, data = plot_agg_traits)
plot(Crown_ratio ~ BA, data = plot_agg_traits)
plot(SLA ~ BA, data = plot_agg_traits)
plot(Bark_at_5cm ~ BA, data = plot_agg_traits)
plot(est_bark_thickness ~ BA, data = plot_agg_traits)
plot(Bark_at_8mm ~ BA, data = plot_agg_traits)
plot(Wood_density ~ BA, data = plot_agg_traits)
plot(Height_mean ~ BA, data = plot_agg_traits)
plot(Light_at_5cm ~ BA, data = plot_agg_traits)
plot(Light_code ~ BA, data = plot_agg_traits)
plot(Density ~ BA, data = plot_agg_traits)
plot(QMD ~ BA, data = plot_agg_traits)
plot(PercentSavannaBA ~ BA, data = plot_agg_traits)
plot(PercentSavannaDens ~ BA, data = plot_agg_traits)

# tiff(filename="./plots/figure_for_proposal.tiff", 
#      type = "cairo",
#      antialias = "gray",
#      compression = "lzw",
#      units="in", 
#      width = 4, 
#      height=4, 
#      pointsize=11, 
#      res=600)
# 
# par(mfrow = c(2,1))
# par(mar = c(4,6,.5,.5),
#     oma = c(0,0,0,0))
# plot(SLA ~ BA, data = plot_agg_traits, 
#      xlab = expression(paste("Basal area (m"^"2", " ha"^"-1", ")")), 
#      ylab = "",
#      cex = 1.5,
#      pch = 21, 
#      bg = "grey")
# mtext("Specific leaf area", side = 2, line = 3.2)
# mtext(expression(paste("(cm"^"2", " g"^"-1", ")")), side = 2, line = 2)
# abline(coef(lm(SLA ~ BA, data = plot_agg_traits)), lwd = 2.5)
# plot(Bark_at_5cm ~ BA, data = plot_agg_traits, 
#      ylab = "",
#      xlab = expression(paste("Basal area (m"^"2", " ha"^"-1", ")")),
#      cex = 1.5,
#      pch = 21, 
#      bg = "grey")
# mtext("Bark thickness\nat 5cm dia. (mm)", side = 2, line = 2.5)
# abline(coef(lm(Bark_at_5cm ~ BA, data = plot_agg_traits)), lwd = 2.5)
# 
# dev.off()
# 
# 
# tiff(filename="./plots/figure_for_proposal2.tiff", 
#      type = "cairo",
#      antialias = "gray",
#      compression = "lzw",
#      units="in", 
#      width = 4, 
#      height=4, 
#      pointsize=11, 
#      res=600)
# 
# par(mfrow = c(2,1))
# par(mar = c(4,6,.5,.5),
#     oma = c(0,0,0,0))
# plot(Height_at_5cm ~ BA, data = plot_agg_traits, 
#      xlab = expression(paste("Basal area (m"^"2", " ha"^"-1", ")")), 
#      ylab = "",
#      cex = 1.5,
#      pch = 21, 
#      bg = "grey")
# mtext("Height at 5 cm (m)", side = 2, line = 3.2)
# abline(coef(lm(Height_at_5cm ~ BA, data = plot_agg_traits)), lwd = 2.5)
# 
# plot(Light_at_5cm ~ BA, data = plot_agg_traits, 
#      ylab = "",
#      xlab = expression(paste("Basal area (m"^"2", " ha"^"-1", ")")),
#      cex = 1.5,
#      pch = 21, 
#      bg = "grey")
# mtext("Light code at 5 cm", side = 2, line = 2.5)
# abline(coef(lm(Light_at_5cm ~ BA, data = plot_agg_traits)), lwd = 2.5)
# 
# dev.off()

# #---------------------------------------------------------------------------------
# # Community analysis
# #---------------------------------------------------------------------------------
# comm_mat_ba <- array(0, dim = c(30, length(unique(plot_data$Code))),
#                      dimnames = list(sites = seq(1:30),
#                                      species = unique(plot_data$Code)))
# 
# i <- j <- 1
# # abundance scaled to % of plot -- not sure why I did this
# for(i in 1:30){
#   for(j in 1:length(unique(plot_data$Code))){
#     BA <- plot_agg_traits$BA[i]
#     comm_mat_ba[i, j] <- sum(plot_data[plot_data$Code == colnames(comm_mat_ba)[j] &
#                                          plot_data$P == i, "CSA_30"]) /10000 / BA
#   }
# }
# 
# 
# 
# # metaMDS(comm, distance = "bray", k = 2, try = 20, trymax = 20,
# #         engine = c("monoMDS", "isoMDS"), autotransform =TRUE,
# #         noshare = (engine == "isoMDS"), wascores = TRUE, expand = TRUE,
# #         trace = 1, plot = FALSE, previous.best, ...)
# 
# nms <- metaMDS(comm_mat_ba)
# 
# nms
# stressplot(nms)
# 
# plot(nms)
# 
# ordiplot(nms,type="n")
# orditorp(nms,display="species",col="red",air=0.01)
# orditorp(nms,display="sites",cex=1.25,air=0.01)
# 
# 
# 
# # Ordination of traits
# 
# community_trait_matrix <- array(0, dim = c(30, 8),
#                                 dimnames = list(sites = seq(1:30),
#                                                 species = seq(1:8)))
# 
# 
# 
# for(i in 1:30){
#   for(j in 1:8){
#     community_trait_matrix[i, j] <- plot_agg_traits[i, j+2]
#   }
# }
# 
# nms_traits <- metaMDS(community_trait_matrix)
# 
# 
# nms_traits
# stressplot(nms_traits)
# 
# plot(nms_traits)
# 
# ordiplot(nms_traits,type="n")
# orditorp(nms_traits,display="species",col="red",air=0.01)
# orditorp(nms_traits,display="sites",cex=1.25,air=0.01)
# 
# 
# library("vegan")

#----------------------------------------------------------------------------------------
# Functional diversity
#----------------------------------------------------------------------------------------
library("FD")

traits <- clean_species[, -c(1, 3)] #this looks like the same as PCA and ANOSIM dataset above,
  #clean up to make one spp x trait dataframe to use every time
rownames(traits) <- traits$Code
traits <- traits[, -1]

abund <- array(0, dim = c(30, length(unique(plot_data$Code))),
                          dimnames = list(sites = seq(1:30),
                                          species = unique(plot_data$Code)))
for(i in 1:30){
  for(j in 1:length(unique(plot_data$Code))){
    abund[i, j] <- sum(plot_data[plot_data$Code == colnames(abund)[j] &
                                         plot_data$P == i, "Density_expand"])
  }
}


abund <- abund[, attr(abund, "dimnames")$species %in% rownames(traits)]
abund <- abund[, order(attr(abund, "dimnames")$species)]

traits <- traits[rownames(traits) %in% attr(abund, "dimnames")$species, ]
traits <- traits[order(rownames(traits)), ]

traits_std <- scale(traits[, sapply(traits, is.numeric)])
  
traits <- traits_std
#avoid "trivially correlated" traits (Mason 2005)
#traits <- traits[, c("Leaf_size", "Height_at_5cm", "Crown_ratio_at_5cm", "SLA", "Bark_at_5cm", "Wood_density")]

#using FD package
# FD calculates CWD traits differently -- figure out why?
test <- functcomp(traits, abund)
plot(test$Leaf_size ~ plot_agg_traits$Leaf_size)
plot(test$Bark_at_5cm ~ plot_agg_traits$Bark_at_5cm)
plot(test$Bark_at_5cm ~ plot_agg_traits$BA)
plot(plot_agg_traits$Bark_at_5cm ~ plot_agg_traits$BA)

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
plot(FD_plots$RaoQ ~ plot_agg_traits$BA)


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
  

plot(FD_plots$RaoQ ~ div)

##MUltipanel figure

tiff(filename="./plots/functional_diversity.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 7, 
     height=7, 
     pointsize=12, 
     res=600)

par(mfrow = c(3, 2))
par(oma = c(2,2,0,0), mar = c(4,4,1,1), family = "sans", bty = 'y')

FRic_lm <- (lm(FD_plots$FRic ~ plot_agg_traits$BA))
summary(FRic_lm)
plot(FD_plots$FRic ~ plot_agg_traits$BA, 
     xlab = "BA",
     ylab = "FRic")
abline(coef(FRic_lm))
text(x = 23, y = par("usr")[4] - ((par("usr")[4] - par("usr")[3]) * .2), 
     labels = paste0("r = ", round(sqrt(summary(FRic_lm)$r.squared), 2)))
text(x = 3, y = par("usr")[4] - ((par("usr")[4] - par("usr")[3]) * .1), labels = "(a)", cex = 1.7)


FDiv_lm <- lm(FD_plots$FDiv ~ plot_agg_traits$BA)
summary(FDiv_lm)
plot(FD_plots$FDiv ~ plot_agg_traits$BA,
     xlab = "BA",
     ylab = "FDiv")
abline(coef(FDiv_lm))
text(x = 23, y = par("usr")[4] - ((par("usr")[4] - par("usr")[3]) * .2), 
     labels = paste0("r = ", round(sqrt(summary(FDiv_lm)$r.squared), 2)))
text(x = 3, y = par("usr")[4] - ((par("usr")[4] - par("usr")[3]) * .1), labels = "(b)", cex = 1.7)


FEve_lm <- lm(FD_plots$FEve ~ plot_agg_traits$BA)
summary(FEve_lm)
plot(FD_plots$FEve ~ plot_agg_traits$BA,
     xlab = "BA",
     ylab = "FEve")
abline(coef(FEve_lm))
text(x = 23, y = par("usr")[4] - ((par("usr")[4] - par("usr")[3]) * .2), 
     labels = paste0("r = ", round(sqrt(summary(FEve_lm)$r.squared), 2)))
text(x = 3, y = par("usr")[4] - ((par("usr")[4] - par("usr")[3]) * .1), labels = "(c)", cex = 1.7)

FDis_lm <- lm(FD_plots$FDis ~ plot_agg_traits$BA)
summary(FDis_lm)
plot(FD_plots$FDis ~ plot_agg_traits$BA,
     xlab = "BA",
     ylab = "FDis")
abline(coef(FDis_lm))
text(x = 23, y = par("usr")[4] - ((par("usr")[4] - par("usr")[3]) * .2), 
     labels = paste0("r = ", round(sqrt(summary(FDis_lm)$r.squared), 2)))
text(x = 3, y = par("usr")[4] - ((par("usr")[4] - par("usr")[3]) * .1), labels = "(d)", cex = 1.7)


div <- diversity(abund, index = "invsimpson")
Simpson_lm <- lm(div ~ plot_agg_traits$BA)
summary(Simpson_lm)
plot(div ~ plot_agg_traits$BA,
     xlab = "BA",
     ylab = "Inverse Simpson's Index")
abline(coef(Simpson_lm))
text(x = 23, y = par("usr")[4] - ((par("usr")[4] - par("usr")[3]) * .2), 
     labels = paste0("r = ", round(sqrt(summary(Simpson_lm)$r.squared), 2)))
text(x = 3, y = par("usr")[4] - ((par("usr")[4] - par("usr")[3]) * .1), labels = "(e)", cex = 1.7)



dev.off()

# 
# #-----------------------------------------------------------------------------------------
# # Digging deeper in the light and fire data
# #-----------------------------------------------------------------------------------------
# #savanna plots
# hist(plot_data$Light.code[plot_data$P %in% c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25:30)])
# #forest plots
# hist(plot_data$Light.code[plot_data$P %in% c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)])
# 
# hist(plot_data[!(plot_data$P %in% c(1,2, 9, 10, 17, 18)), ]$Char.height)
# 
# light_plot_data  <- plot_data[!is.na(plot_data$Light.code) & !is.na(plot_data$Char.height), ]
# newdata <- data.frame(Light.code = seq(1, 5, length.out = 100))
# 
# 
# 
# #all burned plots
# all_light_data <- light_plot_data[!(light_plot_data$P %in% c(1, 2, 9, 10, 17, 18)), ]
# names(all_light_data)[2] <- "Plot"
# 
# plot(Char.height ~ Light.code, data = all_light_data)
# 
# light_char_lm <- lm(Char.height ~ poly(Light.code, 3), data = all_light_data)
# 
# summary(light_char_lm)
# 
# newdata$pred1 <- predict(light_char_lm, newdata = newdata)
# 
# lines(newdata$pred1 ~ newdata$Light.code)
# 
# #just cerradao burned plots
# cerradao_light_data <- light_plot_data[light_plot_data$P %in% c(4, 6, 8, 12, 16, 20, 22, 24), ]
# 
# plot(Char.height ~ Light.code, data = cerradao_light_data)
# 
# light_char_lm <- lm(Char.height ~ poly(Light.code, 3), data = cerradao_light_data)
# 
# summary(light_char_lm)
# 
# newdata$pred2 <- predict(light_char_lm, newdata = newdata)
# 
# lines(newdata$pred2 ~ newdata$Light.code)
# 
# #just cerrado burned plots
# cerrado_light_data <- light_plot_data[light_plot_data$P %in% c(3, 5, 7, 11, 13, 15, 19, 21, 23, 25:30), ]
# 
# plot(Char.height ~ Light.code, data = cerrado_light_data)
# 
# light_char_lm <- lm(Char.height ~ poly(Light.code, 3), data = cerrado_light_data)
# 
# summary(light_char_lm)
# 
# newdata$pred3 <- predict(light_char_lm, newdata = newdata)
# 
# lines(newdata$pred3 ~ newdata$Light.code)
# 
# 
# 
# 
# 
# 
# ### lm for light with interaction with BA
# all_light_data <- join(all_light_data, plot_agg_traits[, c("Plot", "BA")], by = "Plot")
# 
# light_char_lm <- lm(log(Char.height+.1) ~ Light.code + BA, data = all_light_data)
# 
# summary(light_char_lm)
# 
# plot(allEffects(light_char_lm, partial.residuals = TRUE))
# 
# 
# library(earth)
# ba_char_lm <- lm(Char.height ~ BA, data = all_light_data)
# summary(ba_char_lm)
# plot(Char.height ~ BA, data = all_light_data)
# abline(coef(ba_char_lm))
# 
# ba_char_mars <- earth(Char.height ~ BA, data = all_light_data)
# plot(ba_char_mars)
# summary(ba_char_mars)
# 
# plot(Char.height ~ BA, data = all_light_data)
# # abline(coef(ba_char_lm))
# preds <- predict(ba_char_mars, newdata = data.frame("BA" = seq(0, 28, length.out = 200)))
# lines(preds ~ seq(0, 28, length.out = 200), lwd = 2)
# 
# ## Aggregate data
# 
# light_agg <- aggregate(all_light_data, by = list(all_light_data$Plot), FUN = mean)
# light_se <- aggregate(all_light_data, by = list(all_light_data$Plot), FUN = sd)
# plot(light_agg$Char.height ~ light_agg$BA)
# arrows(x0 = light_agg$BA, y0 = light_agg$Char.height,
#        x1 = light_agg$BA, y1 = I(light_agg$Char.height + light_se$Char.height) )
# 
# #Combined figure
# library(RColorBrewer)
# col_map <- densCols(x = all_light_data$BA, y = all_light_data$Char.height, nbin = 128,
#                     colramp = colorRampPalette(brewer.pal(9, "YlGnBu")[-c(1:2)], bias = 2))
# 
# tiff(filename="./plots/char_height_BA.tiff", 
#      type = "cairo",
#      antialias = "gray",
#      compression = "lzw",
#      units="in", 
#      width = 5, 
#      height=3, 
#      pointsize=9, 
#      res=600)
# 
# plot(Char.height/100 ~ BA, data = all_light_data,
#      ylim = c(0, 4),
#      xlab = expression(paste("Tree Basal Area (m"^"2", "ha"^"-1", ")")),
#      ylab = "Char height (m)",
#      pch = 21,
#      cex = 0.7,
#      bg = col_map,
#      col = col_map)
# # abline(coef(ba_char_lm))
# lines(preds/100 ~ seq(0, 27.5, length.out = 200), lwd = 2)
# points(light_agg$Char.height/100 ~ light_agg$BA, pch = 150,
#        bg = "black",
#        cex = 1.5)
# 
# dev.off()
# 
# #*******************************************************************************************
# # Classify species
# #*******************************************************************************************
# habitats <- as.vector(rep("S", times = 30))
# habitats[seq(2,24, by = 2)] <- "F"
# 
# comm_mat_counts <- array(0, dim = c(30, length(unique(plot_data$Code))),
#                          dimnames = list(sites = seq(1:30),
#                                          species = unique(plot_data$Code)))
# 
# for(i in 1:30){
#   for(j in 1:length(unique(plot_data$Code))){
#     comm_mat_counts[i, j] <- nrow(plot_data[plot_data$Code == colnames(comm_mat_counts)[j] &
#                                               plot_data$P == i,])
#   }
# }
# 
# clam <- clamtest(comm_mat_counts, habitats, coverage.limit = 1, alpha = 0.05, specialization = 2/3)
# clam[order(clam$Species), ]
# 
# clam_merge <- clam
# names(clam_merge)[1] <- "Code"
# clam_merge <- join(clam_merge, clean_species[, c(2, 16)], by = c("Code"))
# table(clam_merge$FG, clam_merge$Classes)[, c(2, 1, 3, 4)]


# #--------------------------------------------------------------------------------------------
# # Simulate fire removal
# #--------------------------------------------------------------------------------------------
# 
# trees_after_fire <- plot_data[plot_data$, ] #this causes some problems 
# #probably by filtering out some species without bark estimates
# 
# 
# trees_after_agg <- data.frame(Plot = as.integer(seq(1:30)),
#                               BA = numeric(30),
#                               FG = numeric(30),
#                               PercentSavannaBA = numeric(30),
#                               PercentSavannaDens = numeric(30)
# )
# 
# for(i in 1:30){
#   trees_select <- trees_after_fire[trees_after_fire$P == i, ]
#   trees_after_agg$BA[i] <- sum(trees_select$CSA_BH_expand, na.rm = TRUE) / 10000
#   trees_after_agg$Density[i] <- sum(trees_select$Density_expand, na.rm = TRUE)
#   
#   
#   ba_savanna <- sum(trees_select[trees_select$FG == "S", ]$CSA_BH_expand, na.rm = TRUE) / 10000
#   trees_after_agg$PercentSavannaBA[i] <- ba_savanna / trees_after_agg$BA[i]
#   trees_after_agg$PercentSavannaDens[i] <- sum(trees_select[trees_select$FG == "S", "Density_expand"], na.rm= TRUE) / trees_after_agg$Density[i]
# }
# 
# plot(PercentSavannaDens ~ BA, data = trees_after_agg)
# 
# plot(PercentSavannaDens ~ BA, data = plot_agg_traits)
# 
# plot(NA,
#      xlim = c(0, 3),
#      ylim = c(0, 1),
#      xlab = "Basal Area (m2 per .1 ha)",
#      ylab = "Percent Savanna Species")
# 
# points(PercentSavannaBA ~ BA, data = plot_agg_traits, col = "red")
# abline(coef(lm(PercentSavannaBA ~ BA, data = plot_agg_traits)), col = "red")
# points(PercentSavannaBA ~ BA, data = trees_after_agg, col = "blue")
# abline(coef(lm(PercentSavannaBA ~ BA, data = trees_after_agg)), col = "blue")
# 
# legend(x = 2.3, y = .8, 
#        legend = c("Before Selection (Bark)", "After Selection (Bark)"), 
#        col = c("red", "blue"), 
#        pch = 1)
# arrows(x0 = plot_agg_traits$BA, y0 = plot_agg_traits$PercentSavannaBA,
#        x1 = trees_after_agg$BA, y1 = trees_after_agg$PercentSavannaBA)
# 
# plot(I(trees_after_agg$PercentSavannaDens - plot_agg_traits$PercentSavannaDens) ~ plot_agg_traits$BA)
# 
# plot(I((trees_after_agg$BA - plot_agg_traits$BA)/plot_agg_traits$BA)  ~ plot_agg_traits$BA)
