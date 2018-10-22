library("plyr")

trees_in_plots <- read.csv("./clean data/plot_data2018-01-24.csv")

species_counts <- read.csv(".\\raw data\\lista_spp_plantas_edited_sf082317.csv", stringsAsFactors=FALSE)
  
tab <- data.frame(table(trees_in_plots$Code))
names(tab)[1] <- "Code"

species_counts <- join(species_counts, tab, by = c("Code"))
species_counts$BA <- NA


for(i in 1:nrow(species_counts)){
  species_counts$BA[i] <- sum(trees_in_plots[trees_in_plots$Code == species_counts$Code[i], ]$CSA_30_expand)
}

write.csv(species_counts, "species_counts.csv")
