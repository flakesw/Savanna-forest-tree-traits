# reduce species list
spList_orig <- read.csv("./raw data/lista_spp_plantas_families_sf_2020_04_08.csv", stringsAsFactors = FALSE)
# data <- read.csv("./raw data/field_data.csv", stringsAsFactors = FALSE)

# spList_orig <- spList_orig[spList_orig$Code %in% data$Code, ]

species_for_splink <- c(spList_orig$New.name, paste(spList_orig$Genus, spList_orig$Species, sep = " "))

species_for_splink <- species_for_splink[!(duplicated(species_for_splink))]

write.csv(species_for_splink, file = "./raw data/species_list_splink.csv")
