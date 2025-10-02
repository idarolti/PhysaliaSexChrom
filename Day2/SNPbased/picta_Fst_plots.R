library(tidyverse)

Fst <- read_tsv("picta.fst.10kb.windowed.weir.fst") %>% replace_na(list(WEIGHTED_FST=0)) %>% rename(scaf = CHROM, start = BIN_START, end = BIN_END) %>% mutate(start=as.numeric(start)) %>% filter(WEIGHTED_FST >=0) %>% replace_na(list(WEIGHTED_FST=0))
scaffold_lengths <- read_tsv("Poecilia_picta.fna.fai", col_names = c("scaf","length")) %>% filter(grepl("CM",scaf))

cutoff1 <- round(nrow(Fst)*0.05)
Fst_sort_cut1 <- head(Fst %>% arrange(-WEIGHTED_FST), n=cutoff1) 
min(Fst_sort_cut1$`WEIGHTED_FST`)

# the lowest value is the 5% cutoff

cutoff2 <- round(nrow(Fst)*0.01)
Fst_sort_cut2 <- head(Fst %>% arrange(-WEIGHTED_FST), n=cutoff2) 
min(Fst_sort_cut2$`WEIGHTED_FST`)

# the lowest value is the 1% cutoff

table(Fst_sort_cut2$scaf)

quantile(Fst$`WEIGHTED_FST`, c(0.95, 0.99), na.rm = T)

mean(Fst$WEIGHTED_FST)

Fst_join <- left_join(scaffold_lengths,Fst) %>% filter(WEIGHTED_FST>=0)

#for x axis, we want cumulative bases for each position in the genome for a continuous axis
nCHR <- length(unique(Fst_join$scaf))
Fst_join$BPcum <- 0
s <- 0
nbp <- c()
for (i in unique(Fst_join$scaf)){
  nbp[i] <- max(Fst_join[Fst_join$scaf == i,]$start)
  Fst_join[Fst_join$scaf == i,"BPcum"] <- Fst_join[Fst_join$scaf == i,"start"] + s
  s <- s + nbp[i]
}

axis.set <- Fst_join %>% 
  group_by(scaf) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)

Fst_join %>% ggplot(aes(x=BPcum,y=WEIGHTED_FST,color=as.factor(scaf))) +
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis.set$scaf, breaks = axis.set$center) + 
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x="CHROMOSOME",color="") + 
  guides(color = guide_legend(ncol = 1, byrow = F)) + 
  geom_hline(yintercept = 0.336323, linetype="dotted", color = "black", size=0.75) + 
  geom_hline(yintercept = 0.470588, linetype="dotted", color = "black", size=1.5)  + 
  labs(title="Fst for P. picta with 5% and 1% cutoff")

