
library(ggplot2)

gametologs <- read.csv("/Users/idarolti/Desktop/Folders/Physalia/practicals/day3/gametologs_divergence/6.plot/gametologs_dS_position.txt", row.names=1,header=F,sep="\t")

ggplot(gametologs, aes(x=V3, y=V2)) + 
	geom_point(size=2, alpha=0.7) +
	labs(title="Gametologs divergence",
			x="Chromosomal position (Mb)",
			y="Pairwise divergence dSxy") +
	theme_minimal()
