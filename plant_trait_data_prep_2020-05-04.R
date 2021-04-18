# Process raw data and produce analysis-ready data
# Sam Flake, sflake@gmail.com
#
# This script processes raw data on species traits and stand inventory and outputs 
# species-level mean trait data and cleaned stand inventory data. The inventory data includes
# some kludged together data from other data sources, especially light codes collected from other
# projects at the same study area. It also creates a phylogeny for the study species. 
#
# inputs: raw data on individual traits (./raw data/leaf.csv, wood.csv, field.csv)
#         species list and some ancillary information (./raw data/lista_spp_plantas_families_sf_2020_04_08.csv)
#         species classifications created by classify_species_splink.R (./raw data/species_link_classification.csv)
#         stand inventory data (./raw data/30_parcelas_arvores_2015_complet_out_2015_sf_edits.csv)
#         ancillary information including light codes (./raw data/CharLight_final_102218.csv and FrostDamageAll.csv)
#         maximum heights of species gleaned from the literature (./raw data/species_max_heights.csv)
#
# outputs: species info with classifications (./raw data/sp_info_with_splink.csv)
#          clean individual trait data (.\\clean data\\clean_data.csv)
#          clean species trait data  (./clean data/clean_species.csv)
#          clean species trait data with sample sizes (./clean data/clean_species_with_sample_size.csv)
#          clean plot inventory data (./clean data/plot_data)
#          species phylogeny (./clean data/phylogeny_v)


# load libraries
library("plyr")
# TODO update to dplyr
library("effects")
library("lme4")
library("phytools")
library("ape")
source("./S.PhyloMaker-master/R_codes for S.PhyloMaker") #phylogeny tools
library("devtools")

# import data
leaf <- read.csv(".\\raw data\\leaf_traits.csv", stringsAsFactors = FALSE)
wood <- read.csv(".\\raw data\\twig_traits.csv", stringsAsFactors = FALSE)
field <- read.csv(".\\raw data\\field_data.csv", stringsAsFactors = FALSE)
sp_info <- read.csv(".\\raw data\\lista_spp_plantas_families_sf_2020_04_08.csv", stringsAsFactors = FALSE)
splink_class <- read.csv("./raw data/species_link_classification.csv", stringsAsFactors = FALSE)
plot_data <- read.csv("./raw data/30_parcelas_arvores_2015_complet_out_2015_sf_edits.csv", stringsAsFactors = FALSE)
charlight <- read.csv("./raw data/CharLight_final_102218.csv")
frost <- read.csv("./raw data/FrostDamageAll.csv")
max_heights <- read.csv("./raw data/species_max_heights.csv")

#------------------------------------------------------------------------------
# Clean up trait data
#------------------------------------------------------------------------------
# just pull in the pertinent data from each separate datasheet

leaf_reduced <- leaf[!(leaf$Individual %in% c("Omit", "omit")), c("Individual", "Leaflet.or.leaf.size..mm2.", "Total.leaf.size..mm2.",
                                                                  "Leaf.thickness.1..mm.", "Leaf.thickness.2..mm.", "Leaf.thickness.3..mm.", "SLA..cm2.g.1.", "Method")]
leaf_reduced$Leaf.thickness.1..mm. <- as.numeric(leaf_reduced$Leaf.thickness.1..mm.)
leaf_reduced$Leaf.thickness.2..mm. <- as.numeric(leaf_reduced$Leaf.thickness.2..mm.)
leaf_reduced$Leaf.thickness.3..mm. <- as.numeric(leaf_reduced$Leaf.thickness.3..mm.)
leaf_reduced$SLA..cm2.g.1. <- as.numeric(leaf_reduced$SLA..cm2.g.1.)

wood_reduced <- wood[!(wood$Individual %in% c("Omit", "omit")), c("Individual", "Outer.diameter", "Twig.bark.thickness", "Relative.bark.thickness", "Density..g.mL.")]

field_reduced <- field[!(field$Individual %in% c("Omit", "omit")), c("Individual", "X", "Code", "D30", "DBH", "Ht", "Lower.leaf", "Bark1", "Bark2", "Bark3", "Light", "Date")]
names(field_reduced)[2] <- "Bonus"

#get max DBH and max D30
#split diams at comma, then take the maximum
field_reduced$max_d30 <- vapply(strsplit(field_reduced$D30, split = ","), 
                                1, FUN = function(x){return(as.numeric(max(unlist(x))))})
field_reduced$max_DBH <- vapply(strsplit(field_reduced$DBH, split = ","), 
                                1, FUN = function(x){return(as.numeric(max(unlist(x))))})
field_reduced$D30_eff <- vapply(strsplit(field_reduced$D30, split = ","), 
                           1, FUN = function(x){return(sqrt(sum(((as.numeric(unlist(x)))^2), na.rm = TRUE)))})

sp_info$Life.Form <- as.factor(tolower(sp_info$Life.Form))

#add specieslink classification
names(splink_class)[2] <- "New.name"
sp_info <- join(sp_info, splink_class[, c(2, 6, 7, 8)], by = c("New.name"))
write.csv(sp_info, "./raw data/sp_info_with_splink.csv")

# merge trait data
# from past sam: is there a more elegant way to do this than joinjoinjoinjoin?
# present sam: yes, it could use pipes
clean_data <- join(join(join(join(field_reduced, leaf_reduced, by = c("Individual"), type = "full"),
                   wood_reduced, by = c("Individual"), type = "full"),
                   sp_info, by = c("Code"), type = "left"), 
                   max_heights[, c(2,3)], by = c("Code"), type = "left")
  
# Adding/calculating some variables
clean_data$Mean.Leaf.thickness <- NA
clean_data$Mean.Bark.thickness <- NA
clean_data$Crown.ratio <- NA
clean_data$Taper <- NA
clean_data$Inner.diam <- NA
clean_data$Stem.rel.bark.thick <- NA

for(i in 1:nrow(clean_data)){ #why is this all in a for loop? It seems like easily vectorized stuff
  clean_data$Mean.Leaf.thickness[i] <- mean(unlist(clean_data[i, c("Leaf.thickness.1..mm.", "Leaf.thickness.2..mm.", "Leaf.thickness.3..mm.")]), na.rm=TRUE)
  clean_data$Mean.Bark.thickness[i] <- mean(unlist(clean_data[i, c("Bark1", "Bark2", "Bark3")]), na.rm=TRUE)
  clean_data$Crown.ratio[i] <- (clean_data$Ht[i] - clean_data$Lower.leaf[i]) / clean_data$Ht[i]
  clean_data$Taper[i] <- clean_data$Ht[i] / clean_data$max_d30[i]
  clean_data$Inner.diam[i] <- clean_data$max_d30[i] - clean_data$Mean.Bark.thickness[i]/5
  clean_data$Stem.rel.bark.thick[i] <- (clean_data$Mean.Bark.thickness[i]/5)/ clean_data$max_d30[i]  
}

clean_data$Ht <- clean_data$Ht/100

#Write data
write.csv(clean_data, ".\\clean data\\clean_data.csv")

#------------------------------------------------------------------------------
# Clean up plot data
#------------------------------------------------------------------------------

plot_data <- plot_data[, !(colnames(plot_data) %in% c("X", "X.1", "X.2", "X.3", "X.4"))]

#fix some special characters in diameter data
plot_data$D30 <- gsub("=", "", plot_data$D30)
plot_data$D30 <- gsub("\\(", "", plot_data$D30)
plot_data$D30 <- gsub("\\)", "", plot_data$D30)
plot_data$DBH <- gsub("=", "", plot_data$DBH)
plot_data$DBH <- gsub("\\(", "", plot_data$DBH)
plot_data$DBH <- gsub("\\)", "", plot_data$DBH)


#clean up weird data entry
plot_data_no_equals <- plot_data[plot_data$Species != "=", ]

for(i in 1:nrow(plot_data_no_equals)){
  tree_number <- plot_data_no_equals[i, "Number"]
  
  if(length(plot_data[plot_data$Number == tree_number, "Number"])>1){ #for those tree with  more than one entry
    
    plot_data_no_equals[i, "DBH"] <- paste(plot_data[plot_data$Number == tree_number, "DBH"], collapse = "+") #add other entry to DBH
    if(length(plot_data[plot_data$Number == tree_number, "D30"])>0){ # do the same for D30
      which_d30 <- sapply(plot_data[plot_data$Number == tree_number, "D30"], FUN = function(x){x!=""})
      
      plot_data_no_equals[i, "D30"] <- paste(plot_data[plot_data$Number == tree_number, "D30"][which_d30], collapse = "+")
            
    }
  }
}

#convert species name to a 4-letter code
plot_data_no_equals$Code <- as.character(sapply(plot_data_no_equals$Species, FUN = function(x){toupper(paste(substr(unlist(strsplit(x , split = " ")), 1, 2), collapse = ""))}))

#Clean up codes and some other issues
plot_data_no_equals[plot_data_no_equals$Species == "Campomanesia guaviroba", "Code"] <- "CAGU2"
plot_data_no_equals[plot_data_no_equals$Species == "Chromolaena maximilianii", "Code"] <- "CHMA2"
plot_data_no_equals[plot_data_no_equals$Species == "Ocotea puberula", "Code"] <- "OCPU2"
plot_data_no_equals[plot_data_no_equals$Species == "Ocotea velutina", "Code"] <- "OCVE2"
plot_data_no_equals[plot_data_no_equals$Species == "Diospyros hispida", "Code"] <- "DIBU"
plot_data_no_equals[plot_data_no_equals$Species == "Lippia sidoides", "Code"] <- "LISA"

#fix some data entry errors
plot_data_no_equals[plot_data_no_equals$Number == 4237, "Ht"] <- 5.5
plot_data_no_equals[plot_data_no_equals$Number == 4393, "Ht"] <- 1.9
plot_data_no_equals[plot_data_no_equals$Number == 4319, "Ht"] <- 1.8
plot_data_no_equals[plot_data_no_equals$Number == 1996, "Ht"] <- 4.5

plot_data_no_equals[plot_data_no_equals$Number == 4413, "DBH"] <- 19.3
plot_data_no_equals[plot_data_no_equals$Number == 4204, "DBH"] <- 13.3
plot_data_no_equals[plot_data_no_equals$Number == 4462, "DBH"] <- 14.4
plot_data_no_equals[plot_data_no_equals$Number == 5761, "DBH"] <- 7.8
plot_data_no_equals[plot_data_no_equals$Number == 2275, "DBH"] <- 7.4
plot_data_no_equals[plot_data_no_equals$Number == 3302, "DBH"] <- 7.2
plot_data_no_equals[plot_data_no_equals$Number == 7166, "DBH"] <- 5.8
plot_data_no_equals[plot_data_no_equals$Number == 4204, "DBH"] <- 13.3
plot_data_no_equals[plot_data_no_equals$Number == 4413, "DBH"] <- 19.3
plot_data_no_equals[plot_data_no_equals$Number == 4462, "DBH"] <- 11.4
plot_data_no_equals[plot_data_no_equals$Number == 5289, "DBH"] <- paste("13+11.1")
plot_data_no_equals[plot_data_no_equals$Number == 4477, "DBH"] <- 1.8
plot_data_no_equals[plot_data_no_equals$Number == 1809, "DBH"] <- 1.3
plot_data_no_equals[plot_data_no_equals$Number == 4672, "DBH"] <- 1.1 
plot_data_no_equals[plot_data_no_equals$Number == 5156, "DBH"] <- paste("2+1.5")
plot_data_no_equals[plot_data_no_equals$Number == 5467, "DBH"] <- paste("19+2.2")
plot_data_no_equals[plot_data_no_equals$Number == 2398, "DBH"] <- paste("2.4+1.3+1.2")
plot_data_no_equals[plot_data_no_equals$Number == 1601, "DBH"] <- paste("1.9+1.4+1.4")
plot_data_no_equals[plot_data_no_equals$Number == 7003, "DBH"] <- paste("1.2+1.3+1.1")
plot_data_no_equals[plot_data_no_equals$Number == 2250, "DBH"] <- paste("1.1+1.1+1.5")
plot_data_no_equals[plot_data_no_equals$Number == 7010, "DBH"] <- 11.3
plot_data_no_equals[plot_data_no_equals$Number == 5378, "DBH"] <- paste("6.9+6.3")
plot_data_no_equals[plot_data_no_equals$Number == 5378, "DBH"] <- paste("2.2+2.0+2.3+1.9+3+2.2")

plot_data_no_equals[plot_data_no_equals$Number == 3577, "D30"] <- 8.6
plot_data_no_equals[plot_data_no_equals$Number == 1162, "D30"] <- 8.5
plot_data_no_equals[plot_data_no_equals$Number == 5777, "D30"] <- 7.3
plot_data_no_equals[plot_data_no_equals$Number == 4488, "D30"] <- 7.1
plot_data_no_equals[plot_data_no_equals$Number == 3469, "D30"] <- 4.7
plot_data_no_equals[plot_data_no_equals$Number == 5269, "D30"] <- 4.5
plot_data_no_equals[plot_data_no_equals$Number == 2879, "D30"] <- 4.3
plot_data_no_equals[plot_data_no_equals$Number == 4999, "D30"] <- NA
plot_data_no_equals[plot_data_no_equals$Number == 3453, "D30"] <- 12
plot_data_no_equals[plot_data_no_equals$Number == 979, "D30"] <- 24
plot_data_no_equals[plot_data_no_equals$Number == 2071, "D30"] <- NA
plot_data_no_equals[plot_data_no_equals$Number == 5309, "D30"] <- 2.3
plot_data_no_equals[plot_data_no_equals$Number == 3496, "D30"] <- 1.6

#calculate max diameters and cross-sectional areas for all trees; replace NaNs and Infs with NA
#don't worry about the errors

plot_data_no_equals$Max_DBH <- as.numeric(sapply(plot_data_no_equals$DBH, FUN = function(x){max(as.numeric(unlist(strsplit(x, split = "+", 
                                                                                                                           fixed = TRUE, perl = FALSE, useBytes = FALSE))))}))
plot_data_no_equals$Max_DBH <- ifelse(!is.finite(plot_data_no_equals$Max_DBH), NA, plot_data_no_equals$Max_DBH)

plot_data_no_equals$Max_D30 <- as.numeric(sapply(plot_data_no_equals$D30, FUN = function(x){max(as.numeric(unlist(strsplit(x, split = "+", 
                                                                                                                           fixed = TRUE, perl = FALSE, useBytes = FALSE))))}))
plot_data_no_equals$Max_D30 <- ifelse(!is.finite(plot_data_no_equals$Max_D30), NA, plot_data_no_equals$Max_D30)

plot_data_no_equals$CSA_BH <- as.numeric(sapply(plot_data_no_equals$DBH, FUN = function(x){sum((((as.numeric(unlist(strsplit(x, split = "+", 
                                                                                                                             fixed = TRUE, perl = FALSE, useBytes = FALSE))))/2)^2)*3.14159)}))
plot_data_no_equals$CSA_BH <- ifelse(!is.finite(plot_data_no_equals$CSA_BH), NA, plot_data_no_equals$CSA_BH)

plot_data_no_equals$CSA_30 <- as.numeric(sapply(plot_data_no_equals$D30, FUN = function(x){sum((((as.numeric(unlist(strsplit(x, split = "+", 
                                                                                                                             fixed = TRUE, perl = FALSE, useBytes = FALSE))))/2)^2)*3.14159)}))
plot_data_no_equals$CSA_30 <- ifelse(!is.finite(plot_data_no_equals$CSA_30), NA, plot_data_no_equals$CSA_30)

# if trees are smaller than 5 cm, give them an expansion factor of 4 to account for subplot sampling
plot_data_no_equals$Density_expand <- ifelse(!is.na(plot_data_no_equals$Max_D30), 
                                             ifelse(plot_data_no_equals$Max_D30 > 5, 1, 4), 
											 ifelse(plot_data_no_equals$Max_DBH > 5, 1, 4)) 
                                             
plot_data_no_equals$CSA_30_expand <- plot_data_no_equals$CSA_30 * plot_data_no_equals$Density_expand

plot_data_no_equals$CSA_BH_expand <- plot_data_no_equals$CSA_BH * plot_data_no_equals$Density_expand

plot_data_no_equals$DBH_eff <- sqrt(plot_data_no_equals$CSA_BH/pi)

plot_data_no_equals$D30_eff <- sqrt(plot_data_no_equals$CSA_30/pi)

plot_data_no_equals <- join(plot_data_no_equals, sp_info[, c(4,5,6,7,8,15)], by = c("Code"), type = "left")

#------------------------------------------------------------------------------
#Doing some bark data manipulation
#------------------------------------------------------------------------------
## Get relative bark thickness for stems

bark_table <- clean_data[, c("Code", "max_d30", "Mean.Bark.thickness")]

## Get relative bark thickness for twigs
bark_table_2 <- clean_data[, c("Code", "Outer.diameter", "Twig.bark.thickness")]
bark_table_2$Outer.diameter <- bark_table_2$Outer.diameter / 10
names(bark_table_2) <- names(bark_table)[1:3]

#combine all bark measurements into one data frame
bark_table <- rbind(bark_table, bark_table_2)
bark_table$Code <- as.factor(bark_table$Code)

# fit models of bark ~ diameter and for mean bark thickness per species
bark_model_global  <- lm(Mean.Bark.thickness ~  Code*log(max_d30), data = bark_table)
bark_model_fg <- lm(Mean.Bark.thickness ~  Code, data = bark_table)
saveRDS(bark_model_global, "./clean data/bark_model.RDS")
saveRDS(bark_model_fg, "./clean data/bark_model_fg.RDS")
#-------------------------------------------------------------------------
# Spruce up frost and light data
#-------------------------------------------------------------------------
frost_temp <- frost[, c("placa", "light")]
names(frost_temp)[1] <- "N."
frost_temp$"N." <- as.character(frost_temp$"N.")
frost_temp <- frost_temp[complete.cases(frost_temp), ]
frost_temp <- frost_temp[ frost_temp$"N." != "", ]

#fix one row
frost_temp[frost_temp$N. == 5566, "light"] <- 2

charlight <- join(charlight, frost_temp, by = c("N."), type = "full")

#move data from frost light to charlight for those where charlight is missing light
charlight[is.na(charlight$Light.code), "Light.code"] <- charlight[is.na(charlight$Light.code), "light"]

frost$placa <- as.character(frost$placa)
charlight$"N." <- as.character(charlight$"N.")
charlight_temp <- charlight[charlight$"N." != "" & !(is.na(charlight$"N.")), ]
names(charlight_temp)[5] <- "placa"
charlight_temp <- charlight_temp[, c("placa", "Light.code")]

charlight$Code <- as.character(sapply(as.character(charlight$Especie), FUN = function(x){toupper(paste(substr(unlist(strsplit(x , split = " ")), 1, 2), collapse = ""))}))

charlight[charlight$Code == "DIHI", "Code"] <- "DIBU"
names(charlight)[5] <- "Number"

#add light codes and char heights from char and frost data to census data
plot_data_no_equals <- join(plot_data_no_equals, charlight[, c("Number", "Light.code", "Char.height")], by = c("Number"), type = "left" )

write.csv(plot_data_no_equals, paste0("./Clean data/plot_data", Sys.Date(), ".csv"))

#------------------------------------------------------------------------
# Create species-level data
#------------------------------------------------------------------------
#initialize species-level dataframe
clean_species <- data.frame(Code = unique(clean_data$Code),
                            N = numeric(length(unique(clean_data$Code))),
                            N_leaf = numeric(length(unique(clean_data$Code))),
                            N_wood = numeric(length(unique(clean_data$Code))),
                            N_bark = numeric(length(unique(clean_data$Code))),
                            N_height = numeric(length(unique(clean_data$Code))),
                            Leaf_size = numeric(length(unique(clean_data$Code))),
                            Total_leaf_size = numeric(length(unique(clean_data$Code))),
                            Leaf_thickness = numeric(length(unique(clean_data$Code))),
                            Max_height = numeric(length(unique(clean_data$Code))),
                            Height_at_5cm = numeric(length(unique(clean_data$Code))),
                            Crown_ratio = numeric(length(unique(clean_data$Code))),
                            Light_at_5cm = numeric(length(unique(clean_data$Code))),
                            SLA = numeric(length(unique(clean_data$Code))),
                            Bark_at_5cm = numeric(length(unique(clean_data$Code))),
                            Bark_at_8mm = numeric(length(unique(clean_data$Code))),
                            Wood_density = numeric(length(unique(clean_data$Code))),
                            FG = character(length(unique(clean_data$Code))),
                            Habit = character(length(unique(clean_data$Code))),
                            stringsAsFactors = FALSE
)

# fill the dataframe for each trait; some stuff is just transferred directly and others
# require some calculation, e.g. bark data
for(i in 1:length(unique(clean_data$Code))){
  code_select <- unique(clean_data$Code)[i]
  data_select <- clean_data[clean_data$Code == code_select, ]
  
  clean_species$N[i] <- nrow(data_select)
  clean_species$N_wood[i] <- nrow(data_select[!is.na(data_select$Density..g.mL.), ])
  clean_species$N_leaf[i] <- nrow(data_select[!is.na(data_select$SLA..cm2.g.1.), ])
  clean_species$N_bark[i] <- nrow(data_select[!is.na(data_select$Mean.Bark.thickness), ])
  clean_species$N_height[i] <- nrow(plot_data_no_equals[plot_data_no_equals$Code == code_select & !is.na(plot_data_no_equals$H), ]) + 
                                nrow(data_select[!is.na(data_select$Ht), ])
  
  if(sum(!is.na(data_select$Leaflet.or.leaf.size..mm2.))>=3){
    clean_species$Leaf_size[i] <- mean(data_select$Leaflet.or.leaf.size..mm2., na.rm = TRUE)}
  
  if(sum(!is.na(data_select$Total.leaf.size..mm2.))>=3){
  clean_species$Total_leaf_size[i] <- mean(data_select$Total.leaf.size..mm2., na.rm = TRUE)}
  
  if(sum(!is.na(data_select$Mean.Leaf.thickness))>=3){
  clean_species$Leaf_thickness[i] <- mean(data_select$Mean.Leaf.thickness, na.rm = TRUE)}
  
  if(sum(!is.na(data_select$SLA..cm2.g.1.))>=3){
  clean_species$SLA[i] <-mean(data_select$SLA..cm2.g.1., na.rm = TRUE)}
  
  if(sum(!is.na(data_select$Density..g.mL.))>=3){
  clean_species$Wood_density[i] <- mean(data_select$Density..g.mL., na.rm = TRUE)}
  
  #create regressions for predicting bark thicknesses
  if(sum(!is.na(data_select$Twig.bark.thickness))>=3){
      
      clean_species$Bark_at_8mm[i] <- mean(data_select$Relative.bark.thickness, na.rm = TRUE) * 8
      }else{
    clean_species$Bark_at_8mm[i] <- NA}
    
  if(sum(!is.na(data_select$Mean.Bark.thickness))>=3 & 
     max(data_select$max_d30, na.rm = TRUE) >= 5){
    clean_species$Bark_at_5cm[i] <- predict(bark_model_global, newdata = list(max_d30 = 5, Code = code_select))
    }
  
  
  #Height at 5cm
  # use data from both plots and field measurements, whichever is larger
  ht_subset <- subset(plot_data_no_equals, Code == code_select)[, c("Ht", "D30_eff")]
  ht_subset <- rbind(ht_subset, data_select[, c("Ht", "D30_eff")])
  ht_subset <- ht_subset[ht_subset$D30_eff > 0 & !(is.na(ht_subset$D30_eff)) & !(is.na(ht_subset$Ht)), ]

  if(nrow(ht_subset) >= 5 &
     min(ht_subset$D30_eff, na.rm = TRUE) <= 5 &
     max(ht_subset$D30_eff, na.rm = TRUE) >=5){
    
    height_model <- lm(Ht ~ log(D30_eff), data = ht_subset)
    clean_species$Height_at_5cm[i] <- predict(height_model, newdata = list(D30_eff = 5))
  } else{clean_species$Height_at_5cm[i] <- NA}
  
  
  #Crown ratio at 5cm
  if(nrow(data_select) >= 3){
    clean_species$Crown_ratio[i] <- mean(data_select$Crown.ratio, na.rm = TRUE)
  } else{clean_species$Crown_ratio[i] <- NA}
  
  
  #FG from specieslink data
  clean_species$FG[i] <- as.character(data_select$classification66[1])
  
  #Light at 5cm
  #don't worry about the warnings -- it returns NAs at the end anyway like it should
  light_subset <- subset(plot_data_no_equals, Code == code_select)[, c("Light.code", "D30_eff")]
  light_select <- data_select[, c("Light", "D30_eff")]; names(light_select)[1] <- c("Light.code")
  light_subset <- rbind(light_subset, light_select)
  light_subset <- light_subset[light_subset$D30_eff > 0 & !(is.na(light_subset$D30_eff)) & !(is.na(light_subset$Light.code)), ]
  
  #only return a number if there's at least 5 individuals per species
  if(nrow(light_subset) >= 5 &
     min(light_subset$D30_eff, na.rm = TRUE) < 6 &
     max(light_subset$D30_eff, na.rm = TRUE) >4){
    light_model <- lm(Light.code ~ log(D30_eff), data = light_subset)
    clean_species$Light_at_5cm[i] <- predict(light_model, newdata = list(D30_eff = 5))
  } else{ clean_species$Light_at_5cm[i] <- NA}
  
  
  #Habit
  clean_species$Habit[i] <- as.character(data_select$Life.Form[1])
  
  clean_species$Max_height[i] <- as.numeric(data_select$Max.height[1])
  
}

clean_species_sample_size <- clean_species
clean_species_sample_size <- join(clean_species_sample_size, sp_info, by = c("Code"))
# includes data on sample size for paper
write.csv(clean_species_sample_size, "./clean data/clean_species_with_sample_size.csv")

#get rid of some extraneous data
clean_species <- clean_species[, -c(3:6)]

#sample size per trait? 
N_spp <- apply(clean_species, c(2), FUN = function(x){sum(!is.na(x))})
N_ind <- apply(clean_data, c(2), FUN = function(x){sum(!is.na(x))})


#replace zeroes with NAs
clean_species[, -c(1,2,14,15)] <- apply(clean_species[, -c(1,2,14,15)], c(1,2), FUN = function(x){ifelse(x == 0, NA, x)})
clean_species[, -c(1,2,14,15)] <- apply(clean_species[, -c(1,2,14,15)], c(1,2), FUN = function(x){ifelse(x<0, NA, x)})

write.csv(clean_species, paste0("./clean data/clean_species", Sys.Date(), ".csv"))

#------------------------------------------------------------------------------
# Bring in phylogeny
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Tree using V.Phylomaker
# code borrowed from jinyizju/V.PhyloMaker

# devtools::install_github("jinyizju/V.PhyloMaker")
library("V.PhyloMaker")

#import species list and remove some duplicates
spList_orig <- read.csv("./raw data/lista_spp_plantas_families_sf_2020_04_08.csv", stringsAsFactors = FALSE)

#remove all the other species
spList_orig <- spList_orig[spList_orig$Code %in% clean_species$Code, ]


# split up the binomial
spList_orig$new_genus <- unlist(lapply(spList_orig$New.name, FUN = function(x){strsplit(x, split = " ")[[1]][1]}))
spList_orig$new_species <- unlist(lapply(spList_orig$New.name, FUN = function(x){strsplit(x, split = " ")[[1]][2]}))

#make species list in proper format
spList <- data.frame(species = spList_orig$New.name,
                     genus = spList_orig$new_genus,
                     family = spList_orig$Family,
                     species.relative = spList_orig$species.relative,
                     genus.relative = spList_orig$genus.relative,
                     stringsAsFactors = FALSE)


#Update one genus that has old name in the V.Phylomaker
nodes.info.1[nodes.info.1$genus == "Gochnatia", ]$genus <- "Moquiniastrum"

#how many genera and species are in the  megaphylogeny?
# sum(spList$genus %in% nodes.info.1$genus)
# sum(gsub(" ", "_", spList$species) %in% GBOTB.extended$tip.label)

# generate the phylogeny 
rel <- bind.relative(sp.list=spList, tree=GBOTB.extended, nodes=nodes.info.1)
tree_rel <- phylo.maker(sp.list=rel$species.list, tree=rel$phylo,
                      nodes=rel$nodes.info, scenarios="S3")

write.tree(tree_rel$scenario.3, "./clean data/phylogeny_v")
