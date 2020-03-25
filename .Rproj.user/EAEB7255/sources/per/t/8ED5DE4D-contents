## To do

#load libraries
library("plyr")
library("effects")
library("lme4")
library("phytools")
library("ape")
source("./S.PhyloMaker-master/R_codes for S.PhyloMaker") #phylogeny tools

leaf <- read.csv(".\\raw data\\leaf_traits.csv", stringsAsFactors = FALSE)
wood <- read.csv(".\\raw data\\twig_traits.csv", stringsAsFactors = FALSE)
field <- read.csv(".\\raw data\\field_data.csv", stringsAsFactors = FALSE)
sp_info <- read.csv(".\\raw data\\lista_spp_plantas_families_sf_09252019.csv", stringsAsFactors=FALSE)
plot_data <- read.csv("./raw data/30_parcelas_arvores_2015_complet_out_2015_sf_edits.csv", stringsAsFactors = FALSE)
charlight <- read.csv("./raw data/CharLight_final_102218.csv")
frost <- read.csv("./raw data/FrostDamageAll.csv")
max_heights <- read.csv("./raw data/species_max_heights.csv")

#------------------------------------------------------------------------------
# Clean up trait data
#------------------------------------------------------------------------------
leaf_reduced <- leaf[!(leaf$Individual %in% c("Omit", "omit")), c("Individual", "Leaflet.or.leaf.size..mm2.", "Total.leaf.size..mm2.",
                                                                  "Leaf.thickness.1..mm.", "Leaf.thickness.2..mm.", "Leaf.thickness.3..mm.", "SLA..cm2.g.1.")]
leaf_reduced$Leaf.thickness.1..mm. <- as.numeric(leaf_reduced$Leaf.thickness.1..mm.)
leaf_reduced$Leaf.thickness.2..mm. <- as.numeric(leaf_reduced$Leaf.thickness.2..mm.)
leaf_reduced$Leaf.thickness.3..mm. <- as.numeric(leaf_reduced$Leaf.thickness.3..mm.)
# leaf_reduced[leaf_reduced$SLA..cm2.g.1. == "#VALUE", "SLA..cm2.g.1."] <- NA
# leaf_reduced[leaf_reduced$SLA..cm2.g.1. == 0, "SLA..cm2.g.1."] <- NA
leaf_reduced$SLA..cm2.g.1. <- as.numeric(leaf_reduced$SLA..cm2.g.1.)


wood_reduced <- wood[!(wood$Individual %in% c("Omit", "omit")), c("Individual", "Outer.diameter", "Twig.bark.thickness", "Relative.bark.thickness", "Density..g.mL.")]

field_reduced <- field[!(field$Individual %in% c("Omit", "omit")), c("Individual", "X", "Code", "D30", "DBH", "Ht", "Lower.leaf", "Bark1", "Bark2", "Bark3", "Light")]
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

#is there a more elegant way to do this than joinjoinjoinjoin?
clean_data <- join(join(join(join(field_reduced, leaf_reduced, by = c("Individual"), type = "full"),
                   wood_reduced, by = c("Individual"), type = "left"),
                   sp_info, by = c("Code"), type = "left"), 
                   max_heights[, c(2,3)], by = c("Code"), type = "left")
  
#### Adding/calculating some variables
clean_data$Mean.Leaf.thickness <- NA
clean_data$Mean.Bark.thickness <- NA
clean_data$Crown.ratio <- NA
clean_data$Taper <- NA
clean_data$Inner.diam <- NA
clean_data$Stem.rel.bark.thick <- NA

for(i in 1:nrow(clean_data)){ #why is this all in a for loop?
  clean_data$Mean.Leaf.thickness[i] <- mean(unlist(clean_data[i, c("Leaf.thickness.1..mm.", "Leaf.thickness.2..mm.", "Leaf.thickness.3..mm.")]), na.rm=TRUE)
  clean_data$Mean.Bark.thickness[i] <- mean(unlist(clean_data[i, c("Bark1", "Bark2", "Bark3")]), na.rm=TRUE)
  clean_data$Crown.ratio[i] <- (clean_data$Ht[i] - clean_data$Lower.leaf[i]) / clean_data$Ht[i]
  clean_data$Taper[i] <- clean_data$Ht[i] / clean_data$max_d30[i]
  clean_data$Inner.diam[i] <- clean_data$max_d30[i] - clean_data$Mean.Bark.thickness[i]/5
  clean_data$Stem.rel.bark.thick[i] <- (clean_data$Mean.Bark.thickness[i]/5)/ clean_data$max_d30[i]  
}

clean_data$Ht <- clean_data$Ht/100

setdiff(sp_info$Code, clean_data$Code)

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

plot_data_no_equals <- plot_data[plot_data$Species != "=", ]

#clean up weird data entry
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
plot_data_no_equals[plot_data_no_equals$Species == "Lippia sidoides", "Code"] <- "LISA" #use same values for LISA and LISI


#fix some mis-entered data
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

plot_data_no_equals$Density_expand <- ifelse(!is.na(plot_data_no_equals$Max_D30), 
                                             ifelse(plot_data_no_equals$Max_D30 > 5, 1, 4), ifelse(plot_data_no_equals$Max_DBH > 5, 1, 4)) 
                                             
plot_data_no_equals$CSA_30_expand <- plot_data_no_equals$CSA_30 * plot_data_no_equals$Density_expand

plot_data_no_equals$CSA_BH_expand <- plot_data_no_equals$CSA_BH * plot_data_no_equals$Density_expand

plot_data_no_equals$DBH_eff <- sqrt(plot_data_no_equals$CSA_BH/pi)

plot_data_no_equals$D30_eff <- sqrt(plot_data_no_equals$CSA_30/pi)

plot_data_no_equals <- join(plot_data_no_equals, sp_info[, c(5,6,7,8)], by = c("Code"), type = "left")

write.csv(plot_data_no_equals, paste0("./Clean data/plot_data", Sys.Date(), ".csv"))

#------------------------------------------------------------------------------
#Doing some bark data manipulation
#------------------------------------------------------------------------------

## Get relative bark thickness for stems

bark_table <- clean_data[, c("Code", "max_d30", "Mean.Bark.thickness")]

## Get relative bark thickness for twigs

bark_table_2 <- clean_data[, c("Code", "Outer.diameter", "Twig.bark.thickness")]
bark_table_2$Outer.diameter <- bark_table_2$Outer.diameter / 10
names(bark_table_2) <- names(bark_table)[1:3]
bark_table <- rbind(bark_table, bark_table_2) #combine all bark measurements into one data frame
bark_table$Code <- as.factor(bark_table$Code)

bark_model_global  <- lmer(Mean.Bark.thickness ~  Code*log(max_d30) + (1|Code), data = bark_table)

saveRDS(bark_model_global, "./clean data/bark_model.RDS")

bark_table_fg <- join(bark_table, clean_data[, c("Code", "FG", "Life.Form")], by = c("Code"))

bark_model_fg <- lm(Mean.Bark.thickness ~  FG*Life.Form*log(max_d30), data = bark_table_fg)
saveRDS(bark_model_global, "./clean data/bark_model_fg.RDS")
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

# frost_new <- join(frost, charlight_temp, by = c("placa"), type = "left")

# write.csv(frost_new, paste0("./clean data/FrostDamageMerged", Sys.Date(), ".csv"))


charlight$Code <- as.character(sapply(as.character(charlight$Especie), FUN = function(x){toupper(paste(substr(unlist(strsplit(x , split = " ")), 1, 2), collapse = ""))}))

charlight[charlight$Code == "DIHI", "Code"] <- "DIBU"
names(charlight)[5] <- "Number"

plot_data_no_equals <- join(plot_data_no_equals, charlight[, c("Number", "Light.code", "Char.height")], by = c("Number"), type = "left" )

write.csv(plot_data_no_equals, paste0("./Clean data/plot_data", Sys.Date(), ".csv"))

#------------------------------------------------------------------------
# Create species-level data
#------------------------------------------------------------------------

#initialize species-level dataframe
clean_species <- data.frame(Code = unique(clean_data$Code),
                            N = numeric(length(unique(clean_data$Code))),
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
                            #Leaf_N = numeric(length(unique(clean_data$Code))),
                            Wood_density = numeric(length(unique(clean_data$Code))),
                            FG = character(length(unique(clean_data$Code))),
                            Habit = character(length(unique(clean_data$Code))),
                            stringsAsFactors = FALSE
)

for(i in 1:length(unique(clean_data$Code))){
  code_select <- unique(clean_data$Code)[i]
  data_select <- clean_data[clean_data$Code == code_select, ]
  
  clean_species$N[i] <- nrow(data_select)#length(which(complete.cases(data_select)))
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
  # use data from both plots and field measurements
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
  
  
  #FG
  clean_species$FG[i] <- as.character(data_select$FG[1])
  
  #Light at 5cm
  #don't worry about the warnings -- it returns NAs at the end anyway like it should
  light_subset <- subset(plot_data_no_equals, Code == code_select)[, c("Light.code", "D30_eff")]
  light_select <- data_select[, c("Light", "D30_eff")]; names(light_select)[1] <- c("Light.code")
  light_subset <- rbind(light_subset, light_select)
  light_subset <- light_subset[light_subset$D30_eff > 0 & !(is.na(light_subset$D30_eff)) & !(is.na(light_subset$Light.code)), ]
  
  # print(paste(i, ": ", nrow(light_subset)))
  
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

#replace zeroes with NAs
clean_species[, -c(1,2,14,15)] <- apply(clean_species[, -c(1,2,14,15)], c(1,2), FUN = function(x){ifelse(x == 0, NA, x)})

clean_species[, -c(1,2,14,15)] <- apply(clean_species[, -c(1,2,14,15)], c(1,2), FUN = function(x){ifelse(x<0, NA, x)})

#sample size per trait? 
N_spp <- apply(clean_species, c(2), FUN = function(x){sum(!is.na(x))})
N_ind <- apply(clean_data, c(2), FUN = function(x){sum(!is.na(x))})
  
write.csv(clean_species, paste0("./clean data/clean_species", Sys.Date(), ".csv"))

#------------------------------------------------------------------------------
# Bring in phylogeny
#------------------------------------------------------------------------------

## Phylogenetic analysis


#import species list and remove some duplicates
spList_orig <- read.csv("./raw data/lista_spp_plantas_families_sf_09252019.csv")
spList_orig[!(spList_orig$FG %in% c("S", "F", "G")), "FG"] <- "G"
spList_orig$FG <- droplevels(spList_orig$FG)

#make species list in proper format
spList <- data.frame(species = paste(spList_orig$Genus, spList_orig$Species, sep = " "),
                     genus = spList_orig$Genus,
                     family = spList_orig$Family,
                     stringsAsFactors = FALSE)



phylo<-read.tree("./S.PhyloMaker-master/PhytoPhylo")      # read in the megaphylogeny.    
nodes<-read.csv("./S.PhyloMaker-master/nodes",header=T, sep = "\t")     # read in the nodes information of the megaphylogeny.    
# phylo <- drop.tip(phylo, 8791) #drop Ocotea pulchella to fix the myrtaceae
# phylo <- drop.tip(phylo, 13609) #drop cupania vernalis to fix cupania/matayba
# phylo <- drop.tip(phylo, c())#drop myrcia tomentosa and multiflora



# tree_test <- drop.tip(phylo, phylo$tip.label[!(phylo$tip.label %in% sub(" ", "_", spList$species))])

result<-S.PhyloMaker(spList=spList, tree=phylo, nodes=nodes)      # run the function S.PhyloMaker.    

# plot(result$Scenario.3,cex=1.1,main="Scenario Three")

tree <- result$Scenario.3

#check on myrcias, pouteria ramiflora, 
unique(clean_data$Code)[!(unique(clean_data$Code) %in% spList_orig$Code)]

splist_test <- sub(" ", "_", spList$species)
# splist_test[!(splist_test %in% phylo$tip.label)]
phylo$tip.label[grep("Leptolobium", phylo$tip.label)]

# spList_orig[spList_orig$Code %in% clean_species$Code, "Authority"]

# str(tree)
# 
# saveRDS(tree, file = "phylogeny.RDS")

#relabel tips with code instead of species name
tree$tip.label <- as.character(spList_orig[, "Code"][match(tree$tip.label, sub(" ", "_", spList$species))])

saveRDS(tree, file = "phylogeny_code.RDS")
