
## To do
## Enter light data
## Get species mean light levels
## For each species: get max heights from plot data or from lit; 

#load libraries
library("plyr")

leaf <- read.csv(".\\raw data\\leaf_traits.csv", stringsAsFactors = FALSE)
wood <- read.csv(".\\raw data\\twig_traits.csv", stringsAsFactors = FALSE)
field <- read.csv(".\\raw data\\field_data.csv", stringsAsFactors = FALSE)
sp_info <- read.csv(".\\raw data\\lista_spp_plantas_edited_sf082317.csv", stringsAsFactors=FALSE)
plot_data <- read.csv("./raw data/30_parcelas_arvores_2015_complet_out_2015_sf_edits.csv", stringsAsFactors = FALSE)
charlight <- read.csv("./raw data/CharLight_for_data_entry_sf_092417.csv")
frost <- read.csv("./raw data/FrostDamageAll.csv")

leaf_reduced <- leaf[!(leaf$Individual %in% c("Omit", "omit")), c("Individual", "Leaflet.or.leaf.size..mm2.", "Total.leaf.size..mm2.",
                                                                  "Leaf.thickness.1..mm.", "Leaf.thickness.2..mm.", "Leaf.thickness.3..mm.", "SLA..cm2.g.1.")]
leaf_reduced$Leaf.thickness.1..mm. <- as.numeric(leaf_reduced$Leaf.thickness.1..mm.)
leaf_reduced$Leaf.thickness.2..mm. <- as.numeric(leaf_reduced$Leaf.thickness.2..mm.)
leaf_reduced$Leaf.thickness.3..mm. <- as.numeric(leaf_reduced$Leaf.thickness.3..mm.)
leaf_reduced[leaf_reduced$SLA..cm2.g.1. == "#VALUE", "SLA..cm2.g.1."] <- NA
leaf_reduced[leaf_reduced$SLA..cm2.g.1. == 0, "SLA..cm2.g.1."] <- NA
leaf_reduced$SLA..cm2.g.1. <- as.numeric(leaf_reduced$SLA..cm2.g.1.)


wood_reduced <- wood[!(wood$Individual %in% c("Omit", "omit")), c("Individual", "Outer.diameter", "Twig.bark.thickness", "Relative.bark.thickness", "Density..g.mL.")]

field_reduced <- field[!(field$Individual %in% c("Omit", "omit")), c("Individual", "X", "Code", "D30", "DBH", "Ht", "Lower.leaf", "Bark1", "Bark2", "Bark3", "Light")]
names(field_reduced)[2] <- "Bonus"

#get max DBH and max D30
field_reduced$max_d30 <- vapply(strsplit(field_reduced$D30, split = ","), 
                                1, FUN = function(x){return(as.numeric(max(unlist(x))))})
field_reduced$max_DBH <- vapply(strsplit(field_reduced$DBH, split = ","), 
                                1, FUN = function(x){return(as.numeric(max(unlist(x))))})
# field_reduced$BA <- vapply(strsplit(field_reduced$D30, split = ","), 
#                            1, FUN = function(x){return(as.numeric( sum(3.14159 * (as.numeric(unlist(x))^2 ))))})

clean_data <- join(join(join(leaf_reduced, wood_reduced, by = c("Individual"), type = "full"),
                        field_reduced, by = c("Individual"), type = "full"),
                   sp_info, by = c("Code"), type = "left")



#### Adding/calculating some variables
clean_data$Mean.Leaf.thickness <- NA
clean_data$Mean.Bark.thickness <- NA
clean_data$Crown.ratio <- NA
clean_data$Taper <- NA
clean_data$Inner.diam <- NA
clean_data$Stem.rel.bark.thick <- NA

for(i in 1:nrow(clean_data)){
  clean_data$Mean.Leaf.thickness[i] <- mean(unlist(clean_data[i, c("Leaf.thickness.1..mm.", "Leaf.thickness.2..mm.", "Leaf.thickness.3..mm.")]), na.rm=TRUE)
  clean_data$Mean.Bark.thickness[i] <- mean(unlist(clean_data[i, c("Bark1", "Bark2", "Bark3")]), na.rm=TRUE)
  clean_data$Crown.ratio[i] <- (clean_data$Ht[i] - clean_data$Lower.leaf[i]) / clean_data$Ht[i]
  clean_data$Taper[i] <- clean_data$Ht[i] / clean_data$max_d30[i]
  clean_data$Inner.diam[i] <- clean_data$max_d30[i] - clean_data$Mean.Bark.thickness[i]/5
  clean_data$Stem.rel.bark.thick[i] <- (clean_data$Mean.Bark.thickness[i]/5)/ clean_data$max_d30[i]  
}



#Write data
write.csv(clean_data, ".\\clean data\\clean_data.csv")



#-------------------------------------------------------------------------------------
# Clean up plot data

plot_data <- plot_data[, !(colnames(plot_data) %in% c("X", "X.1", "X.2", "X.3", "X.4"))]


#fix some special characters in diameter data
plot_data$D30 <- gsub("=", "", plot_data$D30)
plot_data$D30 <- gsub("\\(", "", plot_data$D30)
plot_data$D30 <- gsub("\\)", "", plot_data$D30)
plot_data$DBH <- gsub("=", "", plot_data$DBH)
plot_data$DBH <- gsub("\\(", "", plot_data$DBH)
plot_data$DBH <- gsub("\\)", "", plot_data$DBH)

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

#calculate max diameters and cross-sectional areas for all trees; replace NaNs and Infs with NA
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


#convert species name to a 4-letter code
plot_data_no_equals$Code <- as.character(sapply(plot_data_no_equals$Species, FUN = function(x){toupper(paste(substr(unlist(strsplit(x , split = " ")), 1, 2), collapse = ""))}))

write.csv(plot_data_no_equals, paste0("./Clean data/plot_data", Sys.Date(), ".csv"))
#---------------------------------------------------------------
#Doing some bark data manipulation

## Get relative bark thickness for stems

bark_table <- clean_data[, c("Code", "max_d30", "Mean.Bark.thickness")]

## Get relative bark thickness for twigs

bark_table_2 <- clean_data[, c("Code", "Outer.diameter", "Twig.bark.thickness")]
bark_table_2$Outer.diameter <- bark_table_2$Outer.diameter / 10
names(bark_table_2) <- names(bark_table)
bark_table <- rbind(bark_table, bark_table_2) #combine all bark measurements into one data frame

#-------------------------------------------------------------------------
# Spruce up frost and light data
#-------------------------------------------------------------------------

frost_temp <- frost[, c("placa", "light")]
names(frost_temp)[1] <- "N."
frost_temp$"N." <- as.character(frost_temp$"N.")
frost_temp <- frost_temp[complete.cases(frost_temp), ]
frost_temp <- frost_temp[ frost_temp$"N." != "", ]

charlight <- join(charlight, frost_temp, by = c("N."), type = "full")

charlight[is.na(charlight$Light.code), "Light.code"] <- charlight[is.na(charlight$Light.code), "light"]

frost$placa <- as.character(frost$placa)
charlight$"N." <- as.character(charlight$"N.")
charlight_temp <- charlight[charlight$"N." != "" & !(is.na(charlight$"N.")), ]
names(charlight_temp)[5] <- "placa"
charlight_temp <- charlight_temp[, c("placa", "Light.code")]

frost_new <- join(frost, charlight_temp, by = c("placa"), type = "left")

write.csv(frost_new, paste0("./clean data/FrostDamageMerged", Sys.Date(), ".csv"))


charlight$Code <- as.character(sapply(as.character(charlight$Espécie), FUN = function(x){toupper(paste(substr(unlist(strsplit(x , split = " ")), 1, 2), collapse = ""))}))

charlight[charlight$Code == "DIHI", "Code"] <- "DIBU"
names(charlight)[5] <- "Number"

#------------------------------------------------------------------------
# Create species-level data
#------------------------------------------------------------------------

#initialize species-level dataframe
clean_species <- data.frame(Code = unique(clean_data$Code),
                            N = numeric(length(unique(clean_data$Code))),
                            Leaf_size = numeric(length(unique(clean_data$Code))),
                            Leaf_thickness = numeric(length(unique(clean_data$Code))),
                            Max_height = numeric(length(unique(clean_data$Code))),
                            Height_at_5cm = numeric(length(unique(clean_data$Code))),
                            Crown_ratio_at_5cm = numeric(length(unique(clean_data$Code))),
                            Light_at_5cm = numeric(length(unique(clean_data$Code))),
                            SLA = numeric(length(unique(clean_data$Code))),
                            Bark_at_5cm = numeric(length(unique(clean_data$Code))),
                            Bark_at_8mm = numeric(length(unique(clean_data$Code))),
                            #Leaf_N = numeric(length(unique(clean_data$Code))),
                            Wood_density = numeric(length(unique(clean_data$Code))),
                            FG = character(length(unique(clean_data$Code))),
                            Habit = character(length(unique(clean_data$Code))),
                            stringsAsFactors = FALSE
)

i <- 3
for(i in 1:length(unique(clean_data$Code))){
  code_select <- unique(clean_data$Code)[i]
  data_select <- clean_data[clean_data$Code == code_select, ]
  
  clean_species$N[i] <- nrow(data_select)#length(which(complete.cases(data_select)))
  clean_species$Leaf_size[i] <- mean(data_select$Leaflet.or.leaf.size..mm2., na.rm = TRUE)
  clean_species$Leaf_thickness[i] <- mean(data_select$Mean.Leaf.thickness, na.rm = TRUE)
  clean_species$SLA[i] <-mean(data_select$SLA..cm2.g.1., na.rm = TRUE)
  clean_species$Wood_density[i] <- mean(data_select$Density..g.mL., na.rm = TRUE)
  
  #create regressions for predicting bark thicknesses
  twig_model  <- lm(Twig.bark.thickness ~ log(Outer.diameter), data = data_select)
  clean_species$Bark_at_8mm[i] <- predict(twig_model, newdata = list(Outer.diameter = 8))
  
  bark_model  <- lm(Mean.Bark.thickness ~ log(max_d30), data = bark_table[bark_table$Code == code_select, ])
  clean_species$Bark_at_5cm[i] <- predict(bark_model, newdata = list(max_d30 = 5))
  
  #Height at 5cm
  # for species with data, use data from the plots. Else use data from field
  if(nrow(subset(plot_data_no_equals, Code == code_select)) > 5){
    height_model <- lm(Ht ~ as.numeric(Max_D30), data = subset(plot_data_no_equals, Code == code_select & !(is.na(Max_D30))))
    clean_species$Height_at_5cm[i] <- predict(height_model, newdata = list(Max_D30 = 5))
  } else{ #build regression from field data
    height_model <- lm(Ht ~ as.numeric(max_d30), data = data_select)
    clean_species$Height_at_5cm[i] <- predict(height_model, newdata = list(max_d30 = 5)) / 100 #divide by 100 to convert cm to m
  }
  
  #Crown ratio at 5cm
  cr_model <- lm(Crown.ratio ~ max_d30, data = data_select)
  clean_species$Crown_ratio_at_5cm[i] <- predict(cr_model, newdata = list(max_d30 = 5))
  
  #FG
  clean_species$FG[i] <- as.character(data_select$FG[1])
  
  #Light at 5cm
  light_data <- charlight[charlight$Code == code_select, ]
  
  if(nrow(subset(plot_data_no_equals, Code == code_select)) > 5){
    light_data <- join(light_data, plot_data_no_equals, by = c("Number"), type = "left")
    light_model <- lm(as.numeric(Light.code) ~ Max_D30, data = light_data)
    
    plot(as.numeric(Light.code) ~ Max_D30, data = light_data)
    
    clean_species$Light_at_5cm[i] <- predict(light_model, newdata = list(Max_D30 = 5))
    } #else{ #build regression from field data
      #light_data <- join(light_data, data_select, by = c("Number"), type = "left")
      #light_model <- lm(as.numeric(Light.code) ~ max_d30, data = light_data)
      
      #clean_species$Light_at_5cm[i] <- predict(light_model_model, newdata = list(max_d30 = 5))
  #}
  
  
  #Habit
  clean_species$Habit[i] <- as.character(data_select$Life.Form[1])
  
  # clean_species$Max_height[i] <-
}

write.csv(clean_species, paste0("./clean data/clean_species", Sys.Date(), ".csv"))






# 
# 
# 
# 
# 
# 
# 
# 
# 
# #-----------------------------------------------------------------------------------------------------
# ### Not needed
# 
# 
# #################
# ## Compare codes
# #################
# 
# names(leaf)
# leaf$Individual
# names(wood)
# wood$Individual
# names(field)
# field$Individual
# 
# check_names <- leaf[!(leaf$Individual %in% wood$Individual), ]
# 
# check_names2 <- wood[!(wood$Individual %in% leaf$Individual), ]
# 
# check_names3 <- leaf[!(leaf$Individual %in% field$Individual), ]
# 
# check_names4 <- field[!(field$Individual %in% leaf$Individual), ]
# 
# check_names5 <- wood[!(wood$Individual %in% field$Individual), ]
# 
# check_names6 <- field[!(field$Individual %in% wood$Individual), ]
# 
# #################################
# ## Find the max diameter for each species and currently sampled diameters
# ################################
# tree_data <- read.csv("30_parcelas_arvores_2015_complet_out_2015.csv", stringsAsFactors = FALSE)
# names(tree_data)[c(5,6)] <- c("Tree.number", "Genus.species")
# species_table <- read.csv("species_table.csv")
# 
# tree_data$max_D30 <- NA
# 
# for(i in (1:nrow(tree_data))){
#   
#   tree_data$max_D30[i] <- max(as.numeric(unlist(strsplit(tree_data[i, "D..30cm."], split = "+", 
#                                                          fixed = TRUE, perl = FALSE, useBytes = FALSE))))
#   
# }
# 
# tree_data_merge <- merge(tree_data, species_table, sort = FALSE)
# 
# need_sample <- data.frame(Code = unique(field$Code),
#                           max_D30 = numeric(length(unique(field$Code))),
#                           max_sampled = numeric(length(unique(field$Code))))
# 
# for(i in 1:nrow(need_sample)){
#   need_sample$max_D30[i] <- max(tree_data_merge[as.character(tree_data_merge$Code) == 
#                                                   as.character(need_sample$Code[i]), ]$max_D30,
#                                 na.rm=TRUE)
#   need_sample$max_sampled[i] <- max(max(as.numeric(unlist(strsplit(field[as.character(field$Code) == 
#                                             as.character(need_sample$Code[i]), ]$D30, 
#                                             split = ",", fixed=TRUE)))))
#   
# }
# 
# 
# need_sample_with_abund <- merge(need_sample, species_table)
# write.csv(need_sample_with_abund, "need_bark_sample.csv")
