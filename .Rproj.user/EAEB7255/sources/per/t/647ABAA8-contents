
#----------------------------------------------------------------------------------------
# Functional diversity
#----------------------------------------------------------------------------------------

traits <- clean_species_reduced_trans[, -c(1, 3, 4, 15, 16)] 
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
                                   plot_data$P == i, "CSA_30_expand"], na.rm = TRUE)
  }
}

abund_notraits <- abund

setdiff(rownames(traits), attr(abund, "dimnames")$species)
setdiff(attr(abund, "dimnames")$species, rownames(traits))

abund <- abund[, attr(abund, "dimnames")$species %in% rownames(traits)]
abund <- abund[, order(attr(abund, "dimnames")$species)]

traits <- traits[rownames(traits) %in% attr(abund, "dimnames")$species, ]
traits <- traits[order(rownames(traits)), ]

FD_plots <- dbFD(traits, abund, w.abun = TRUE, corr = "sqrt", stand.x = TRUE, stand.FRic = TRUE)

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
abline(coef(FRic_lm)[, 1], lty = ifelse(coef(FRic_lm)[2, 4] < 0.05, 1, 2))
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
abline(coef(FDiv_lm)[, 1], lty = ifelse(coef(FDiv_lm)[2, 4] < 0.05, 1, 2))
text(x = 15, y = par("usr")[3] + ((par("usr")[4] - par("usr")[3]) * .3), 
     pos = 4, 
     labels = paste("r =", (round(sqrt(FDiv_lm$r.squared), 2) * sign(coef(FDiv_lm)[2,1])) ))
text(x = 15, y = par("usr")[3] + ((par("usr")[4] - par("usr")[3]) * .2),
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
abline(coef(FEve_lm)[, 1], lty = ifelse(coef(FEve_lm)[2, 4] < 0.05, 1, 2))
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
abline(coef(FDis_lm)[, 1], lty = ifelse(coef(FDis_lm)[2, 4] < 0.05, 1, 2))
text(x = 15, y = par("usr")[3] + ((par("usr")[4] - par("usr")[3]) * .3), 
     pos = 4, 
     labels = paste("r =", (round(sqrt(FDis_lm$r.squared), 2) * sign(coef(FDis_lm)[2,1])) ))
text(x = 15, y = par("usr")[3] + ((par("usr")[4] - par("usr")[3]) * .2),
     pos = 4, 
     labels = ifelse(coef(FDis_lm)[2, 4] > 0.001,
                     paste("p =", round(coef(FDis_lm)[2, 4], 3)),
                     "p < 0.001"))
mtext(side = 3, text = "(d)", 
      at = (par()$usr[1] + .1*(par()$usr[2] - par()$usr[1])), 
      line = 0.6, cex = 0.8)
mtext(side = 2, text = "FDis", line = 2.2)

div <- diversity(abund_notraits, index = "invsimpson")
Simpson_lm <- summary(lm(div ~ plot_agg_traits$BA))
plot(div ~ plot_agg_traits$BA,
     xlab = "",
     ylab = "",
     pch = 21,
     bg = "grey")
abline(coef(Simpson_lm)[, 1], lty = ifelse(coef(Simpson_lm)[2, 4] < 0.05, 1, 2))
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

J <- diversity(abund_notraits, index = "shannon")/log(specnumber(abund))
pielou_lm <- summary(lm(J ~ plot_agg_traits$BA))
plot(J ~ plot_agg_traits$BA,
     xlab = "",
     ylab = "",
     pch = 21,
     bg = "grey")
abline(coef(pielou_lm)[, 1], lty = ifelse(coef(pielou_lm)[2, 4] < 0.05, 1, 2))
text(x = 0, y = par("usr")[4] - ((par("usr")[4] - par("usr")[3]) * .2), 
     pos = 4, 
     labels = paste("r =", (round(sqrt(pielou_lm$r.squared), 2) * sign(coef(pielou_lm)[2,1])) ))
text(x = 0, y = par("usr")[4] - ((par("usr")[4] - par("usr")[3]) * .3),
     pos = 4, 
     labels = ifelse(coef(pielou_lm)[2, 4] > 0.001,
                     paste("p =", round(coef(pielou_lm)[2, 4], 3)),
                     "p < 0.001"))
mtext(side = 3, text = "(e)", 
      at = (par()$usr[1] + .1*(par()$usr[2] - par()$usr[1])), 
      line = 0.6, cex = 0.8)
mtext(side = 2, text = "Pielou's evenness", line = 2.2)

mtext(side = 1, text = expression(paste("Basal area (m"^"2", " ha"^"-1", ")")), 
      line = 1.1, outer = TRUE)



dev.off()



#############################################################################################
# PCA with weighted circle sizes
#############################################################################################
abund <- array(0, dim = c(30, length(unique(plot_data$Code))),
               dimnames = list(sites = seq(1:30),
                               species = unique(plot_data$Code)))

#weight abundance by basal area. We already did this once but doing it again for good measure
for(i in 1:30){
  for(j in 1:length(unique(plot_data$Code))){
    abund[i, j] <- sum(plot_data[plot_data$Code == colnames(abund)[j] &
                                   plot_data$P == i, "CSA_30_expand"], na.rm = TRUE)
  }
}

abund <- abund[, attr(abund, "dimnames")$species %in% rownames(traits)]
abund <- abund[, order(attr(abund, "dimnames")$species)]

traits <- traits[rownames(traits) %in% attr(abund, "dimnames")$species, ]
traits <- traits[order(rownames(traits)), ]

#pull out pca data that also has abundance data
pca_data <- pca_data[rownames(pca_data) %in% colnames(abund), ]
#match functional type for each species
pca_groups <- clean_species_reduced_trans[match(rownames(pca_data), clean_species_reduced_trans$Code), c("FG")]

hull <- chull(y = pca_data[, 2], x = pca_data[, 1])
FRic_max <- polyarea(x = pca_data[hull, 1], y = pca_data[hull, 2])

abund_tot <- abund[c(26, 13, 4), ]

abund_tot <- abund_tot[, order(colnames(abund_tot))]
pca_data <- pca_data[order(rownames(pca_data)), ]

#match scores with reduced dataset
sp_scores <- sp_scores[rownames(sp_scores) %in% rownames(pca_data), ]
sp_scores <- sp_scores[order(rownames(sp_scores)), ]


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
par(oma = c(2,2,0,0))
par(mar = c(3, 2, 4, 1))
color_groups <- data.frame(group = as.character(c("F", "G", "S")),
                           col = as.character(c("black", "#1b9e77", "#E69F00")))

lets <- c("(a)", "(b)", "(c)")

for(i in 1:3){
  color_points <- as.character(color_groups[match(pca_groups, color_groups$group), "col"])[abund_tot[i, ] > 0]
  sp_scores_red <- sp_scores[which(abund_tot[i, ] > 0), ]
  
  abund_red <- abund_tot[i, which(abund_tot[i, ] > 0)]
  
  plot(NA,
       xlab = "", ylab = "", xaxt = "n", yaxt = "n",
       xlim = c(min(sp_scores[, 1]),
                max(sp_scores[, 1])),
       ylim = c(min(sp_scores[, 2]),
                max(sp_scores[, 2])))
  
  mtext(outer = FALSE, cex = 1.25, text = ifelse(i == 1, expression(paste("Savanna, BA = 6.47 m"^"2", " ha"^"-1")),
                                                 ifelse(i == 2, expression(paste("Woodland, BA = 9.05 m"^"2", " ha"^"-1")),
                                                        expression(paste("Forest, BA = 24.1 m"^"2", " ha"^"-1")))))
  
  axis(1, cex = 1.2)
  axis(2, cex = 1.2)
  
  # points(y = sp_scores[, 2], x = sp_scores[, 1], col = color_points, bg = color_points, pch = 21, cex = .7)
  
  symbols(y = sp_scores_red[, 2], x = sp_scores_red[, 1], 
          circles=sqrt(abund_red/pi), inches=0.15,
          bg = color_points,
          fg = color_points,
          add = TRUE)
  # text(y = sp_scores_red[, 2], x = sp_scores_red[, 1], labels = rownames(sp_scores_red), col = "black", cex = .7)
  abline(h = 0)
  abline(v = 0)
  
  #Calculate and plot FRic
  hull <- chull(y = sp_scores_red[, 2], x = sp_scores_red[, 1])
  hull2 <- convhulln(as.matrix(sp_scores_red[, c(1,2)]), "FA") #same result as chull
  polygon(sp_scores_red[hull, c(1,2)])
  text(x = -0.65, y = -1.2, 
       labels = paste("FRic = ",
                      round(polyarea(x = sp_scores_red[hull, 1], y = sp_scores_red[hull, 2])/FRic_max, 3)),
       cex = 1.4)
  
  #calculate FDiv
  #calculate g (center of gravity of hull vertices)
  x <- sp_scores_red[hull, c(1,2)]
  g <- apply(x, 2, FUN = mean)
  
  #calculate distance from g
  dG <- apply(sp_scores_red[, c(1,2)], 1, FUN = function(x){dist(rbind(x, g), method = "euclidean")})
  dG_bar <- mean(dG)
  
  # distances from mean distance from g
  dist_from_g_bar <- dG - dG_bar
  weights <- abund_red / sum(abund_red)
  
  delta_d <- weights %*% dist_from_g_bar
  
  abs_delta_d <- weights %*% abs(dist_from_g_bar)
  
  FDiv <- (delta_d + dG_bar)/(abs_delta_d + dG_bar)
  
  #plot centroid
  points(x = g[1], y = g[2],  cex = 2, pch = 21, bg = "blue")
  #circle depicting center of gravity
  ts <- seq(0, 2*pi, length.out = 100)
  polygon(x = cos(ts)*dG_bar + g[1], y = sin(ts)*dG_bar + g[2], lty = 2) 
  
  #lines depicting distance from center of gravity
  x_dist <- sp_scores_red[, 1] - g[1]
  y_dist <- sp_scores_red[, 2] - g[2]
  r <- sqrt((sp_scores_red[, 1] - g[1])^2 + (sp_scores_red[, 2] - g[2])^2)
  
  for(j in 1:nrow(sp_scores_red)){
    segments(x0 = sp_scores_red[j, 1], x1 = (x_dist[j] / r[j])*dG_bar + g[1], 
             y0 = sp_scores_red[j, 2], y1 = (y_dist[j] / r[j])*dG_bar + g[2],
             lwd = (weights[j])^(1/4) * 4)
  }
  
  text(x = -0.65, y = -1.1, labels = paste("FDiv = ", round(FDiv, digits = 3)),
       cex = 1.4)
  
  text(x = -1, y = 1.2, cex = 2, labels = lets[i])
  
  if(i == 2){
    
    legend(x = -0.05, y = -.7, 
           legend = c("Forest", "Generalist", "Savanna", "Centroid"), 
           col = c("black", "#1b9e77", "#E69F00", "dark blue"), 
           pt.bg = c("black", "#1b9e77", "#E69F00", "dark blue"), 
           pch = 21,
           cex = 1.5,
           bty = "n")
  }
  
}


mtext(side = 1, outer = TRUE, text = "PC1", cex = 1.2)
mtext(side = 2, outer = TRUE, text = "PC2", cex = 1.2)
dev.off()
