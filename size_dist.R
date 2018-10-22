## exploring tree size distributions
library("plyr")

trees_in_plots <- read.csv("./clean data/plot_data2018-01-24.csv")
plot_ba <- aggregate(trees_in_plots[, c("P", "CSA_BH_expand")], by = list(trees_in_plots$P), FUN = function(x){sum(x, na.rm=TRUE)})
plot_ba$CSA_BH_expand <- plot_ba$CSA_BH_expand/1000
trees_in_plots <- join(trees_in_plots, plot_ba, by = "P")

occo <- trees_in_plots[trees_in_plots$Code == "OCCO", ]
hist(2*sqrt(occo$CSA_BH/3.14))
for(i in 1:30){
  occo_plot <- occo[occo$P == i, ]
  if(nrow(occo_plot) > 0){
    hist(2*sqrt(occo_plot$CSA_BH/3.14), main = paste0("Plot", i))
  }
}


gopo <- trees_in_plots[trees_in_plots$Code == "GOPO", ]
hist(2*sqrt(gopo$CSA_BH/3.14))
for(i in 1:30){
  gopo_plot <- gopo[gopo$P == i, ]
  if(nrow(gopo_plot) > 0){
    hist(2*sqrt(gopo_plot$CSA_BH/3.14), main = paste0("Plot", i))
  }
}


votu <- trees_in_plots[trees_in_plots$Code == "VOTU", ]

for(i in 1:30){
  votu_plot <- votu[votu$P == i, ]
  if(nrow(votu_plot) > 0){
    hist(2*sqrt(votu_plot$CSA_BH/3.14), main = paste0("Plot", i))
  }
}

peob <- trees_in_plots[trees_in_plots$Code == "PEOB", ]

for(i in 1:30){
  peob_plot <- peob[peob$P == i, ]
  if(nrow(peob_plot) > 0){
    hist(2*sqrt(peob_plot$CSA_BH/3.14), main = paste0("Plot", i))
  }
}

stro <- trees_in_plots[trees_in_plots$Code == "STRO", ]

for(i in 1:30){
  stro_plot <- stro[stro$P == i, ]
  if(nrow(stro_plot) > 0){
    hist(2*sqrt(stro_plot$CSA_BH/3.14), main = paste0("Plot", i))
  }
}

maac <- trees_in_plots[trees_in_plots$Code == "MAAC", ]
hist(2*sqrt(maac$CSA_BH/3.14))

for(i in 1:30){
  maac_plot <- maac[maac$P == i, ]
  if(nrow(maac_plot) > 0){
    hist(2*sqrt(maac_plot$CSA_BH/3.14), main = paste0("Plot", i))
  }
}

