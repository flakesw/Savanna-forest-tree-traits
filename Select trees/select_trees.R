setwd("C:\\Users\\Sam\\Documents\\Research\\Plant traits\\Select trees")
tree_data <- read.csv("30_parcelas_arvores_2015_complet_out_2015.csv", stringsAsFactors = FALSE)

species_table <- table(tree_data$Espécie)[rev(order(table(tree_data$Espécie)))]

tree_data$Habitat <- ifelse(tree_data$P %in% c("1", "3", "5", "7", "9", "11", "13", "15", "17",
                                                  "19", "21", "23"), "Savanna", 
                            ifelse(tree_data$P %in% c("2", "4", "6", "8", "10", "12", "14",
                                                      "16", "18", "20", "22", "24"), "Forest", 
                                   ifelse(tree_data$P %in% c("25", "26", "27", "28", "29", "30"), "Campo", NA)))

open_trees <- table(tree_data[tree_data$Habitat %in% c("Savanna", "Campo"), ]$Espécie)
open_trees <- open_trees[order(open_trees)]

forest_trees <- table(tree_data[tree_data$Habitat == "Forest", ]$Espécie)
forest_trees <- forest_trees[order(forest_trees)]


tree_data$Unique_subplot <- paste0(tree_data$P, tree_data$Sub)
tree_data$BA <- NA
for( i in 1:nrow(tree_data)){
  trees_select <- tree_data[tree_data$Unique_subplot == tree_data$Unique_subplot[i], ]
  tree_data[i, "BA"] <- sum(as.numeric(trees_select$sum.DAP), na.rm = TRUE)
}

abundant_species <- names(subset(species_table, species_table>10))

abundant_open <- names(subset(open_trees, open_trees > 20))


opar <- par()
par(oma = c(7, 1,1,1))
boxplot(BA ~ Espécie, data = tree_data[tree_data$Espécie %in% abundant_open & 
                                         tree_data$Habitat %in% c("Savanna", "Campo"), ],
        las = 2)


write.csv(species_table, "species_table.csv")

###################

sla <- read.csv("field_data_traits_summer_2017_temp.csv")

table(sla$Code)

length(table(sla$Code))

######################
