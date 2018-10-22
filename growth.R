library("plyr")
library("effects")
plot_data <- read.csv("./Clean data/plot_data2018-01-24.csv")
growth_data <- read.csv("./raw data/tree_growth_041718_plots_13_14_25_26.csv", stringsAsFactors = FALSE)
names(growth_data)[4] <- "Number"
growth_data <- growth_data[!(is.na(as.numeric(as.character(growth_data$cresc)))), ]
growth_data <- growth_data[!(growth_data$H.2 > 0), ]

growth_data$cresc <- as.numeric(as.character(growth_data$cresc))
growth_data$For.DAP2 <- as.numeric(as.character(growth_data$For.DAP2))
growth_data$Form.DAP1 <- as.numeric(as.character(growth_data$Form.DAP1))
growth_data <- growth_data[!(abs(growth_data$cresc) > 3), ]

growth_data$bai <- 3.1416 * growth_data$cresc^2 * sign(growth_data$cresc)
growth_data$old_ba <- 3.1416 * growth_data$For.DAP2^2
growth_data$rel_growth <- growth_data$bai/growth_data$old_ba

growth_data <- growth_data[growth_data$rel_growth > -0.1, ]


all_plot_data <- join(plot_data, growth_data, by = "Number", type = "inner")

boxplot(cresc ~ Light.code, data = all_plot_data)
boxplot(bai ~ Light.code, data = all_plot_data)
boxplot(rel_growth ~ Light.code, data = all_plot_data) #still need to remove some trees that got topkilled/broken 

growth_lm <- lm(cresc ~ Light.code, data = all_plot_data[all_plot_data$Code == "MYFA", ])
summary(growth_lm)
plot(allEffects(growth_lm))



plot(cresc ~ Light.code, data = all_plot_data[all_plot_data$Code == "GOPO", ],
     pch = 21, col = "red", bg = "red")
abline(coef(lm(cresc ~ Light.code, data = all_plot_data[all_plot_data$Code == "GOPO", ])),
       col = "red", lwd = 2)
points(cresc ~ Light.code, data = all_plot_data[all_plot_data$Code == "MYFA", ],
       pch = 21, col = "blue", bg = "blue")
abline(coef(lm(cresc ~ Light.code, data = all_plot_data[all_plot_data$Code == "MYFA", ])),
       col = "blue", lwd = 2)



#self-thinning
plot(log(plot_agg_traits$QMD) ~ log(plot_agg_traits$Density))

self_thin_rq <- rq(log(QMD) ~ log(Density), tau = 0.95, data = plot_agg_traits)
summary(self_thin_rq)
lines(predict(self_thin_rq) ~ plot_agg_traits$Density)

     