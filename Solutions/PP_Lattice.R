
# Draw and save three lattice graphs by feeding interaction type
#		for predator, prey and predator/prey mass, before calcualting the
#		mean and median for each.

library(lattice)

MyDF <- read.csv("../Data/EcolArchives-E089-51-D1.csv")

# plot three lattice graphs and save as separate .pdf files:
pdf("../Results/Pred_Lattice.pdf", 10, 10)
densityplot(~log(Predator.mass) | Type.of.feeding.interaction, data=MyDF)
dev.off()

pdf("../Results/Prey_Lattice.pdf", 10, 10)
densityplot(~log(Prey.mass) | Type.of.feeding.interaction, data=MyDF)
dev.off()

pdf("../Results/Size_Ratio_Lattice.pdf", 10, 10)
densityplot(~log(Prey.mass/Predator.mass) | Type.of.feeding.interaction, data=MyDF)
dev.off()

# Using the "Smart" R way with tapply:
Prey_mean <- tapply(log(MyDF$Prey.mass), MyDF$Type.of.feeding.interaction, mean)
Prey_median <- tapply(log(MyDF$Prey.mass), MyDF$Type.of.feeding.interaction, median)
Predator_mean <- tapply(log(MyDF$Prey.mass), MyDF$Type.of.feeding.interaction, mean)
Predator_median <- tapply(log(MyDF$Prey.mass), MyDF$Type.of.feeding.interaction, median)
Size_Ratio_mean <- tapply(log(MyDF$Predator.mass / MyDF$Prey.mass), MyDF$Type.of.feeding.interaction, mean)
Size_Ratio_median <- tapply(log(MyDF$Predator.mass / MyDF$Prey.mass), MyDF$Type.of.feeding.interaction, median)
a <- data.frame(levels(MyDF$Type.of.feeding.interaction), Prey_mean, Prey_median, Predator_mean, Predator_median, Size_Ratio_mean, Size_Ratio_median)

# Save the results to csv
write.csv(a, "../Results/PP_Results.csv", row.names=FALSE) # save data frame to a .csv file
