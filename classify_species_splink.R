# Classify species from Specieslink data
# Sam Flake, sflake@gmail.com
#
# inputs: raw SpeciesLink records (./raw data/SaoPaulo_spLink_data)
#         list of species and some species information (./raw data/lista_spp_plantas_families_sf_2020_04_08)
# outputs: species classified by their occurrence (./raw data/species_link_classification)
# This script processes data downloaded manually from SpeciesLink and classifies each species
# according to its habitat preferences, as indicated in the "notes" and "localities" field 
# in herbarium records.

library(plyr)
library(flora)



splink_results <- list.files(path = "./raw data/SaoPaulo_spLink_data", pattern = "\\.txt$", full.names = TRUE)
splink_data <- ldply(splink_results, .fun = function(x){
                          read.table(x, quote = "", header = TRUE, 
                          sep = "\t", encoding = "UTF-8", dec = ".", fill = TRUE)})
sp_info <- read.csv(".\\raw data\\lista_spp_plantas_families_sf_2020_04_08.csv", stringsAsFactors = FALSE)

# check to make sure names are up to date and spelled right. 
# I revised species list to match REFLORA on 2020-06-05. 
# TODO: automate name correction, flora package does some of this
# flora_names <- get.taxa(sp_info$New.name)
# flora_names[!(flora_names$search.str %in% sp_info$New.name), ]

# reformat data to lump subspecies and clean up errors
# remove special characters
splink_data$species <- gsub("[^[:alnum:] ]", "", splink_data$species) 
# take first word from each species row
splink_data$species_epi <- gsub("([A-Za-z]+).*", "\\1", splink_data$species) 
# recombine genus + species
splink_data$Species_name <- paste(splink_data$genus, splink_data$species_epi, sep = " ") 

sp_count <- table(splink_data$Species_name)

# remove species with low counts, to remove misspellings
splink_data <- splink_data[splink_data$Species_name %in% names(sp_count[(sp_count > 3)]), ]

# do some manual proofing :(
# TODO automate

splink_data[splink_data$Species_name == "Cybistax antisyphillitica", "Species_name"] <- "Cybistax antisyphilitica"
splink_data[splink_data$Species_name == "Copaifera langsdorfii", "Species_name"] <- "Copaifera langsdorffii"
splink_data[splink_data$Species_name == "Duguetia furfuraceae", "Species_name"] <- "Duguetia furfuracea"
splink_data[splink_data$Species_name == "Kielmeyera coriaceae", "Species_name"] <- "Kielmeyera coriacea"
splink_data[splink_data$Species_name == "Myrsine coriaceae", "Species_name"] <- "Myrsine coriacea"
splink_data[splink_data$Species_name == "Chromolaena maximiliani", "Species_name"] <- "Chromolaena maximilianii"
splink_data[splink_data$Species_name == "Pinus elliotti", "Species_name"] <- "Pinus elliottii"

# replace species names that aren't in the new name list with their new names. This is to correct for entries
# still using the old name
splink_data[!(splink_data$Species_name %in% sp_info$New.name), "Species_name"] <- 
  sp_info[match(splink_data[!(splink_data$Species_name %in% sp_info$New.name), "Species_name"], sp_info$Old.name), "New.name"]

min(table(splink_data$Species_name))
max(table(splink_data$Species_name))


# words to search for for each habitat type
savanna_words <- c("campo", "cerrado", "savanna", "savannah", "grassland", "campis", "campina")
s_unlisted <- paste(unlist(savanna_words), collapse = "|")
forest_words <- c("mata", "floresta", "forest", "cerradao", "cerradão", "bosque")
f_unlisted <- paste(unlist(forest_words), collapse = "|")


# some place names and other phrases including "campo" to remove
campo_remove <- c("Campo de Congomhas", "Campo Congomhas", "Campo de Congonhas", "Campo Congonhas",
                  "Sao Jose dos Campos", "Sao José dos Campos", "São José dos Campos", 
                  "Campos do Jordao", "Campos do Jordão", "Campos de Jordao", "Campos de Jordão", 
                  "Campos das Sete Lagôas", "Campos das Sete Lagoas", "Campos da Sete Lagôas", "Campos da Sete Lagoas",
                  "Campos de Sete Lagôas", "Campos de Sete Lagoas", "Campo de Sete Lagôas", "Campo de Sete Lagoas",
                  "Amrérico de Campos", 
                  "Campos ad Yapanema", 
                  "São Bernardo do Campo", "Sao Bernardo do Campo",
                  "Luzia do Campo",
                  "Fazenda Campo Grande",
                  "Campos do Butantã",
                  "Bairro Campo Novo",
                  "Campo Alegre",
                  "campo de futbol",
                  "campo e nautica", "Campo e Naútica", "Campo e Náutica",
                  "Campo Alegre",
                  "Novo Campo Elísios",
                  "José de Campos Novaes",
                  "Mandioqueiro-do-campo", "Mandioquinha do campo",
                  "Campo Limpo Paulista",
                  "Curso Taxonomia de Campo de Espécies Vegetais",
                  "Campos, M.C.R. 2008",
                  "cambará-do-campo",
                  "Açuquinha do campo",
                  "araticum-do-campo",
                  "peroba-do-campo",
                  "alecrim do campo", "alecrim-do-campo", "alegrim do campo",
                  "Campos Novos Paulista",
                  "A vegetação dos remanescentes de cerrado no Estado de São Paulo",
                  "Viabilidade da Conservação dos Remanescentes de Cerrado",
                  "campinas",
                  "Amrérico de Campos",
                  "cedro do campo", "cedro de campo", "pimenteira do cerrado",
                  "pimenta do campo",
                  "angico do cerrado", "angico branco do cerrado",
                  "Análise florística e fitossociológica do estrato herbáceo-subarbustivo do cerrado na reserva biológica de Mogi Guaçu e em Itirapina, SP",
                  "Abundance and distribution of native and alien grasses in 'cerrado' (brazilian savanna) biological reserve",
                  "Composição florísrtica de uma área de cerrado em Moji Mirim (SP)",
                  "Composição florística e fitossociológica da vegetação de cerrado no município de Luís Antonio",
                  "Florística e fitossociologia de um cerrado marginal brasileiro",
                  "Florística do Cerrado de Emas (Pirassununga, SP)",
                  "Composição florísrtica de uma área de cerrado em Moji Mirim (SP)",
                  "Reserva Biológica Campininha",
                  "Trilha do Campo",
                  "Estação do Campo",
                  "campo de Aviação",
                  "elementos de cerrado"
                  )

forest_remove <- c("horto florestal",
                   "instituto florestal", "Inst. Florestal",
                   "floresta estadual",
                   "Museu Florestal",
                   "reserva florestal",
                   "usina da mata", "Usina DaMata",
                   "Mata da Cumbica",
                   "Fazenda da Mata",
                   "Jardim do Bosque",
                   "mata do viveiro",
                   "matado viveiro",
                   "matadouro",
                   "mata da mariana",
                   "corta a mata",
                   "Bosque São José",
                   "Sub-bosque", "Subosque", 
                   "Bosque dos Jequitibás",
                   "rodeada antigamente por mata",
                   "microtopográficas e edáficas da Floresta Ombrófila Densa do Núcleo Picinguaba/PESM, Ubatuba",
                   "Estrutura da vegetação arbórea e regeneração natural em remanescentes de mata ciliar do Rio Mogi Guaçu - SP",
                   "Florística e fitossociologia do estrato arbóreo de um remanescente de mata ciliar do rio Jacaré-Pepira, Brotas, SP",
                   "Fitossociologia do componente arbóreo de um trecho de mata em São Paulo, SP"
                   )

# combine into one string
find.string <- tolower(paste(paste(unlist(campo_remove), collapse = "|"), paste(unlist(forest_remove), collapse = "|"), sep = "|")) 

# remove confounding words from notes and locality and combine into one string
splink_data$clean_notes <- tolower(paste(gsub(find.string, "", splink_data$notes, ignore.case = TRUE), gsub(find.string, "", splink_data$locality, ignore.case = TRUE), sep = " "))

splink_data$savanna <- NA
splink_data$forest <- NA

# look for savanna and forest words in the notes/localities
splink_data[grep(s_unlisted, splink_data$clean_notes, ignore.case = TRUE), "savanna"] <- 1
splink_data[grep(f_unlisted, splink_data$clean_notes, ignore.case = TRUE), "forest"] <- 1

# remove duplicates and species that do not match
splink_data <- splink_data[!duplicated(paste(splink_data$collectioncode, splink_data$catalognumber)), ]
splink_data <- splink_data[!is.na(splink_data$Species_name), ]

# summarize data
splink_aggregate <- data.frame(species = unique(splink_data$Species_name)[order(unique(splink_data$Species_name))],
                               savanna_count = numeric(length(unique(splink_data$Species_name))),
                               forest_count = numeric(length(unique(splink_data$Species_name))))

splink_aggregate$savanna_count <- aggregate(splink_data$savanna, 
                                            by = list(splink_data$Species_name), 
                                            FUN = function(x){sum(x, na.rm = TRUE)})[, 2]

splink_aggregate$forest_count <- aggregate(splink_data$forest, 
                                            by = list(splink_data$Species_name), 
                                            FUN = function(x){sum(x, na.rm = TRUE)})[, 2]

# classify based upon thresholds of habitat preference
threshold <- 0.8

splink_aggregate$total_count <- splink_aggregate$savanna_count + splink_aggregate$forest_count
splink_aggregate$classification80 <- ifelse(splink_aggregate$savanna_count/splink_aggregate$total_count > threshold,
                                          "S", ifelse(splink_aggregate$forest_count/splink_aggregate$total_count > threshold,
                                                      "F", "G"))

threshold <- 0.7
splink_aggregate$classification70 <- ifelse(splink_aggregate$savanna_count/splink_aggregate$total_count > threshold,
                                            "S", ifelse(splink_aggregate$forest_count/splink_aggregate$total_count > threshold,
                                                        "F", "G"))

threshold <- 0.66
splink_aggregate$classification66 <- ifelse(splink_aggregate$savanna_count/splink_aggregate$total_count > threshold,
                                            "S", ifelse(splink_aggregate$forest_count/splink_aggregate$total_count > threshold,
                                                        "F", "G"))

#remove data with low sample size
for(i in 1:nrow(splink_aggregate)){
  if(splink_aggregate$total_count[i] <= 5){
    splink_aggregate$classification80[i] <- NA
    splink_aggregate$classification70[i] <- NA
    splink_aggregate$classification66[i] <- NA
  }
}

write.csv(splink_aggregate, file = "./raw data/species_link_classification.csv")
