library(plyr)
library(rvest)
library(dplyr)
library(curl)
library(V8)

##Classify species from Specieslink data
splink_results <- list.files(path = "./raw data/SaoPaulo_spLink_data", pattern = "\\.txt$", full.names = TRUE)
splink_data <- ldply(splink_results, .fun = function(x){
                          read.table(x, quote = "", header = TRUE, 
                          sep = "\t", encoding = "UTF-8", dec = ".", fill = TRUE)})


#reformat data to lump subspecies and clean up errors
#remove special characters
splink_data$species <- gsub("[^[:alnum:] ]", "", splink_data$species) 
#take first word from each species row
splink_data$species_epi <- gsub("([A-Za-z]+).*", "\\1", splink_data$species) 
#recombine genus + species
splink_data$Species_name <- paste(splink_data$genus, splink_data$species_epi, sep = " ") 

sp_count <- table(splink_data$Species_name)

#remove species with low counts, which are all due to misspellings or unclear ID
splink_data <- splink_data[splink_data$Species_name %in% names(sp_count[(sp_count > 3)]), ]

#do some manual proofing :(
splink_data <- splink_data[-grep("cf", splink_data$Species_name), ]

splink_data[splink_data$Species_name == "Cybistax antisyphillitica", "Species_name"] <- "Cybistax antisyphilitica"
splink_data[splink_data$Species_name == "Copaifera langsdorfii", "Species_name"] <- "Copaifera langsdorffii"
splink_data[splink_data$Species_name == "Duguetia furfuraceae", "Species_name"] <- "Duguetia furfuracea"
splink_data[splink_data$Species_name == "Kielmeyera coriaceae", "Species_name"] <- "Kielmeyera coriacea"
splink_data[splink_data$Species_name == "Myrsine coriaceae", "Species_name"] <- "Myrsine coriacea"
splink_data[splink_data$Species_name == "Chromolaena maximiliani", "Species_name"] <- "Chromolaena maximilianii"
splink_data[splink_data$Species_name == "Pinus elliotti", "Species_name"] <- "Pinus elliottii"


#words to search for for each habitat type
savanna_words <- c("campo", "cerrado", "savanna", "savannah", "grassland", "campis", "campina")
s_unlisted <- paste(unlist(savanna_words), collapse = "|")
forest_words <- c("mata", "floresta", "forest", "cerradao", "cerradão")
f_unlisted <- paste(unlist(forest_words), collapse = "|")

#------------------------------------------------------------------------------
# #follow links and scrape html from links to other herbaria
# 
# 
# #todo: extract URLs from specieslink data
# #transform URLs into proper form for smithsonian 
# #finish XML extraction from smithonian
# #merge data back into splink dataframe
# 
# #if link is from smithsonian
# 
# 
# add <- "http://n2t.net/ark:/65665/3e6b15c49-86e0-41bd-b947-3025e8590b55"
# add2 <- "https://collections.nmnh.si.edu/search/botany/?ark=ark:/65665/3e6b15c4986e041bdb9473025e8590b55"
# suffix <- gsub("-", "", url_parse(add)$path)
# add3 <- url(paste0("https://collections.nmnh.si.edu/search/botany/?ark=", suffix))
# link <- read_html(add3)
# # 
# # link2 <- read_html(url(add2))
# # link1 <- read_html(url(add))
# # 
# # node1 <- link1 %>% html_nodes('body')
# # html_text(node1)
# 
# # -----------------------------------------------------------------------------
# #htmlunit method -- doesn't work
# # library(htmlunit)
# # library(rvest)
# # library(purrr)
# # library(tibble)
# 
# # library(rvest)
# 
# # js_pg <- hu_read_html("https://collections.nmnh.si.edu/search/botany/?ark=ark:/65665/3e6b15c4986e041bdb9473025e8590b55",
# #                       emulate = "ie", js_delay = 10, ignore_ssl_errors = FALSE)
# # 
# # desc_nodes <- xml_find_all(js_pg, "//tr") %>% html_children()
# # grep("Locality", html_text(js_pg))
# 
# #------------------------------------------------------------------------------
# # # try with splashr
# # install.packages("splashr")
# # library("splashr")
# # library("stevedore")
# # 
# # install_splash()
# # splash_container <- start_splash()
# 
# # library("splashr")
# # docker <- stevedore::docker_client()
# # install_splash()
# # splash_container <- start_splash()
# 
# #------------------------------------------------------------------------------
# #with decapitated
# # devtools::install_github("hrbrmstr/decapitated")
# # library("decapitated")
# # decapitated::download_chromium("C:/Users/Sam/Documents/R/chromium")
# # C:/Users/Sam/Documents/R/chromium/chrome-win32/chrome.exe
# 
# # chrome_read_html("https://collections.nmnh.si.edu/search/botany/?ark=ark:/65665/3e6b15c4986e041bdb9473025e8590b55")
# # chrome_read_html("http://n2t.net/ark:/65665/3e6b15c49-86e0-41bd-b947-3025e8590b55")
# # doesn't work -- locks up when trying to load website
# 
# 
# -------------------------------------------------------------------------------
# #try with RSelenium
# # devtools::install_github("ropensci/RSelenium")
# # install.packages("RSelenium")
# # library("RSelenium")
# # 
# RSelenium::checkForServer()
# rD <- RSelenium::rsDriver(browser = "firefox", port = 4444L)
# 
# remDr <- rD[["client"]]
# remDr$navigate("http://www.google.com/ncr")
# remDr$navigate("http://www.bbc.com")
# remDr$close()
# 
# system("taskkill /im java.exe /f > nul 2>&1", intern = FALSE, ignore.stdout=TRUE, ignore.stderr=TRUE)
# 
# # 
# # 
# # remDr <- remoteDriver(remoteServerAddr = "localhost"
# #                       , port = 4444
# #                       , browserName = "firefox"
# #                       )
# # remDr$open()
# 
# 
# # https://thatdatatho.com/2019/01/22/tutorial-web-scraping-rselenium/
# 
# driver <- RSelenium::rsDriver(browser = "chrome",
#                               chromever =
#                                 system2(command = "wmic",
#                                         args = 'datafile where name="C:\\\\Program Files (x86)\\\\Google\\\\Chrome\\\\Application\\\\chrome.exe" get Version /value',
#                                         stdout = TRUE,
#                                         stderr = TRUE) %>%
#                                 stringr::str_extract(pattern = "(?<=Version=)\\d+\\.\\d+\\.\\d+\\.") %>%
#                                 magrittr::extract(!is.na(.)) %>%
#                                 stringr::str_replace_all(pattern = "\\.",
#                                                          replacement = "\\\\.") %>%
#                                 paste0("^",  .) %>%
#                                 stringr::str_subset(string =
#                                                       binman::list_versions(appname = "chromedriver") %>%
#                                                       dplyr::last()) %>%
#                                 as.numeric_version() %>%
#                                 max() %>%
#                                 as.character())
# 
# remote_driver <- driver[["client"]] 
# remote_driver$navigate("https://www.latlong.net/convert-address-to-lat-long.html")

# 
# #------------------------------------------------------------------------------
# #if link comes from new york bg, easy with rvest
# add <- "http://sweetgum.nybg.org/science/vh/specimen_details.php?irn=148244"
# page <- read_html(add)
# 
# desc_nodes <- xml_find_all(page, "//ul[contains(@class, 'occurrence-terms')]") %>% html_children()
# 
# loc_node <- xml_node(page, xpath = "/html/body/div[2]/section[5]/div/div/div/div/div[1]/div/div/div/ul/li[3]/p")
# hab_node <- xml_node(page, xpath = "/html/body/div[2]/section[5]/div/div/div/div/div[1]/div/div/div/ul/li[4]/p")
# paste(html_text(loc_node), html_text(hab_node))

#todo: check on GBIF

#some place names and other phrases including "campo" to remove
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
                  "cedro do campo", 
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
                  "campo de Aviação"
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
                   "rodeada antigamente por mata",
                   "microtopográficas e edáficas da Floresta Ombrófila Densa do Núcleo Picinguaba/PESM, Ubatuba",
                   "Estrutura da vegetação arbórea e regeneração natural em remanescentes de mata ciliar do Rio Mogi Guaçu - SP",
                   "Florística e fitossociologia do estrato arbóreo de um remanescente de mata ciliar do rio Jacaré-Pepira, Brotas, SP",
                   "Fitossociologia do componente arbóreo de um trecho de mata em São Paulo, SP"
                   )

find.string <- tolower(paste(paste(unlist(campo_remove), collapse = "|"), paste(unlist(forest_remove), collapse = "|"), sep = "|")) 

splink_data$clean_notes <- tolower(paste(gsub(find.string, "", splink_data$notes, ignore.case = TRUE), gsub(find.string, "", splink_data$locality, ignore.case = TRUE), sep = " "))

splink_data$savanna <- NA
splink_data$forest <- NA

#look for savanna and forest words in the notes/localities
splink_data[grep(s_unlisted, splink_data$clean_notes, ignore.case = TRUE), "savanna"] <- 1
splink_data[grep(f_unlisted, splink_data$clean_notes, ignore.case = TRUE), "forest"] <- 1

splink_aggregate <- data.frame(species = unique(splink_data$Species_name)[order(unique(splink_data$Species_name))],
                               savanna_count = numeric(length(unique(splink_data$Species_name))),
                               forest_count = numeric(length(unique(splink_data$Species_name))))

splink_aggregate$savanna_count <- aggregate(splink_data$savanna, 
                                            by = list(splink_data$Species_name), 
                                            FUN = function(x){sum(x, na.rm = TRUE)})[, 2]

splink_aggregate$forest_count <- aggregate(splink_data$forest, 
                                            by = list(splink_data$Species_name), 
                                            FUN = function(x){sum(x, na.rm = TRUE)})[, 2]

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
