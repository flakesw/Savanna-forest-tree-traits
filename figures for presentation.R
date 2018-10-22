# Figures for talk

# Light code at a reference diameter
code_select <- "VOTU"
data_select <- clean_data[clean_data$Code == code_select, ]
light_data <- join(light_data, plot_data_no_equals, by = c("Number"), type = "left")
light_model <- lm(as.numeric(Light.code) ~ log(Max_D30), data = light_data)

plot(as.numeric(Light.code) ~ Max_D30, data = light_data,
     xlab = "Diameter (cm)",
     ylab = "Light code",
     main = "Vochysia tucanorum",
     cex.lab = 1.2)

lines(seq(0.05, 55, length.out = 500), predict(light_model, newdata = list("Max_D30" = seq(0.05, 55, length.out = 500))), lwd = 2)
abline(v = 5, lty = 2, lwd = 3)



# shifts with succession
par(oma = c(.5, .5, .5, 0),
    mar = c(5.1, 4.6, 4.1, 2.1))

plot(SLA ~ BA, data = plot_agg_traits,
     ylab = expression(paste("Specific leaf area (cm"^"2", " g"^"-1", ")")),
     xlab = expression(paste("Basal area (m"^"2", " ha"^"-1", ")")),
     main = "Specific leaf area",
     cex.lab = 1.3)
abline(coef(lm(SLA ~ BA, data = plot_agg_traits)))
summary(lm(SLA ~ BA, data = plot_agg_traits))


plot(Light_at_5cm ~ BA, data = plot_agg_traits,
     ylab = "Light code at 5cm",
     xlab = expression(paste("Basal area (m"^"2", " ha"^"-1", ")")),
     main = "Light code",
     cex.lab = 1.3)
abline(coef(lm(Light_at_5cm ~ BA, data = plot_agg_traits)))
summary(lm(Light_at_5cm ~ BA, data = plot_agg_traits))

plot(est_bark_thickness ~ BA, data = plot_agg_traits,
     ylab = "Bark thickness (mm)",
     xlab = expression(paste("Basal area (m"^"2", " ha"^"-1", ")")),
     main = "Estimated bark thickness",
     cex.lab = 1.3)
summary(lm(est_bark_thickness ~ poly(BA, 2), data = plot_agg_traits))
b_mod <- lm(est_bark_thickness ~ poly(BA, 2), data = plot_agg_traits)
eff <- predict(b_mod, list("BA" = seq(0, 30)))
lines(seq(0,30), eff)


plot(Height_at_5cm ~ BA, data = plot_agg_traits,
     ylab = "Height (m) at 5cm diam",
     xlab = expression(paste("Basal area (m"^"2", " ha"^"-1", ")")),
     main = "Specific leaf area",
     cex.lab = 1.3)
abline(coef(lm(Height_at_5cm ~ BA, data = plot_agg_traits)))
summary(lm(Height_at_5cm ~ BA, data = plot_agg_traits))




#PCA

tiff(filename="./plots/PCA_for_presentation.tiff", 
     type = "cairo",
     antialias = "gray",
     compression = "lzw",
     units="in", 
     width = 7, 
     height=7, 
     pointsize=9, 
     res=600)
color_groups <- data.frame(group = as.character(c("F", "G", "S")),
                           col = c("#1b9e77", "#d95f02", "#7570b3"),
                           pch = c(17, 18, 19))
color_points <- as.character(color_groups[match(pca_groups, color_groups$group), "col"])
pch_points <- as.character(color_groups[match(pca_groups, color_groups$group), "pch"])

plot(NA,
     xlab = "PC1", ylab = "PC2",
     xlim = c( (min(c(sp_scores[, 1], trait_scores[, 1]) - .5)), 
               (max(c(sp_scores[, 1], trait_scores[, 1])) + 0.4)),
     ylim = c( (min(c(sp_scores[, 2], trait_scores[, 2]) - .5)), 
               (max(c(sp_scores[, 2], trait_scores[, 2])) + 0.4)))


points(y = sp_scores[, 2], x = sp_scores[, 1], col = color_points, pch = as.numeric(pch_points), cex = 1)
abline(h = 0)
abline(v = 0)

for(i in 1:3){
  ordiellipse(pca_results, pca_groups, show.groups = color_groups[i, 1], 
              col = as.character(color_groups[i, 2]), choices = c(1,2),
              kind = "sd", conf = 0.90)
}

legend(x = -1.7, y = 2.5, legend = c("Forest", "Generalist", "Savanna"), 
       col = color_groups$col[c(2,3,1)], pch = color_groups$pch,
       cex = 1.5)

segments(x0 = 0, y0 = 0, x1 = trait_scores[c(1,2,3,5,6,7,8,10), 1], y1 = trait_scores[c(1,2,3,5,6,7,8,10), 2])
text(x = trait_scores[c(1,2,3,5,6,7,8,10), 1], y = trait_scores[c(1,2,3,5,6,7,8,10), 2], labels = trait_names[c(1,2,3,5,6,7,8,10)], cex = 1.3)

dev.off()
