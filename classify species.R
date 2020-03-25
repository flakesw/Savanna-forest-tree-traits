#*******************************************************************************************
# Classify species
# analysis for appendix xx
# todo: add Fridley analysis?
#*******************************************************************************************
habitats <- as.vector(rep("S", times = 30))
habitats[seq(2,24, by = 2)] <- "F"

comm_mat_counts <- array(0, dim = c(30, length(unique(plot_data$Code))),
                         dimnames = list(sites = seq(1:30),
                                         species = unique(plot_data$Code)))

for(i in 1:30){
  for(j in 1:length(unique(plot_data$Code))){
    comm_mat_counts[i, j] <- nrow(plot_data[plot_data$Code == colnames(comm_mat_counts)[j] &
                                              plot_data$P == i,])
  }
}

clam <- clamtest(comm_mat_counts, habitats, coverage.limit = 1, alpha = 0.05, specialization = 2/3)
clam[order(clam$Species), ]

clam_merge <- clam
names(clam_merge)[1] <- "Code"
clam_merge <- join(clam_merge, clean_species[, c(2, 16)], by = c("Code"))
table(clam_merge$FG, clam_merge$Classes)[, c(2, 1, 3, 4)]
