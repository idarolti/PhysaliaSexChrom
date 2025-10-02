library(tidyverse)
library(ggpubr)

# given the command run for SNP density, we have an easy way to isolate males and females. want to load in the individual SNP densities calculated for each individual, get the mean male and mean female SNP Density in that window, run permutation tests to see which windows are significantly different.

# load in the male data:
myFiles <- list.files(pattern="^MALE*")

# build a backbone
backbone <- read_delim(file=myFiles[1],delim = "\t",col_names = T) 
firstsplit <- strsplit(myFiles[1], "_")
secondsplit <- strsplit(firstsplit[[1]][3], "\\.")
column_ID <- paste(firstsplit[[1]][2],secondsplit[[1]][1],sep="_")
backbone.upgrade <- backbone %>% unite(LOCATION,c("CHROM","BIN_START"),sep=":") %>% dplyr::rename(!!column_ID := "VARIANTS/KB") %>% select(-SNP_COUNT)

#loop the other files

for(i in 2:length(myFiles)){
  file <- read_delim(myFiles[i],delim = "\t",col_names = T) 
  firstsplit <- strsplit(myFiles[i], "_")
  secondsplit <- strsplit(firstsplit[[1]][3], "\\.")
  column_ID <- paste(firstsplit[[1]][2],secondsplit[[1]][1],sep="_")
  file.upgrade <- file %>% unite(LOCATION,c("CHROM","BIN_START"),sep=":") %>% dplyr::rename(!!column_ID := "VARIANTS/KB") %>% select(-SNP_COUNT)
  backbone.upgrade <- full_join(backbone.upgrade,file.upgrade,by="LOCATION")
}

SNPdensity.males <- backbone.upgrade %>% separate(LOCATION, c("scaf","base"),sep=":")

# load the female data:

myFiles <- list.files(pattern="^FEMALE*")

# build a backbone
backbone <- read_delim(file=myFiles[1],delim = "\t",col_names = T) 
firstsplit <- strsplit(myFiles[1], "_")
secondsplit <- strsplit(firstsplit[[1]][3], "\\.")
column_ID <- paste(firstsplit[[1]][2],secondsplit[[1]][1],sep="_")
backbone.upgrade <- backbone %>% unite(LOCATION,c("CHROM","BIN_START"),sep=":") %>% dplyr::rename(!!column_ID := "VARIANTS/KB") %>% select(-SNP_COUNT)

#loop the other files

for(i in 2:length(myFiles)){
  file <- read_delim(myFiles[i],delim = "\t",col_names = T) 
  firstsplit <- strsplit(myFiles[i], "_")
  secondsplit <- strsplit(firstsplit[[1]][3], "\\.")
  column_ID <- paste(firstsplit[[1]][2],secondsplit[[1]][1],sep="_")
  file.upgrade <- file %>% unite(LOCATION,c("CHROM","BIN_START"),sep=":") %>% dplyr::rename(!!column_ID := "VARIANTS/KB") %>% select(-SNP_COUNT)
  backbone.upgrade <- full_join(backbone.upgrade,file.upgrade,by="LOCATION")
}


SNPdensity.females <- backbone.upgrade %>% separate(LOCATION, c("scaf","base"),sep=":")

#join them
SNPdensity <- full_join(SNPdensity.males,SNPdensity.females)

#make new tibble, change type for base, and replace NAs
SNPdensity.rows <- SNPdensity
SNPdensity.rows$base <- as.numeric(as.character(SNPdensity.rows$base))
SNPdensity.rows <- SNPdensity.rows %>% replace_na(list(Males = 0, Females = 0)) %>% replace(is.na(.), 0)

# use a subsetter to get the male and female average densities 
male_n <- ncol(SNPdensity.rows[ , grepl( "_male" , names( SNPdensity.rows ) ) ])
males_true <- bind_cols(SNPdensity.rows %>% select(1:2),SNPdensity.rows[ , grepl( "_male" , names( SNPdensity.rows ) ) ] %>% mutate(mean_Males = rowMeans(.))) %>% select(scaf,base,mean_Males)

female_n <- ncol(SNPdensity.rows[ , grepl( "_female" , names( SNPdensity.rows ) ) ])
females_true <- bind_cols(SNPdensity.rows %>% select(1:2),SNPdensity.rows[ , grepl( "_female" , names( SNPdensity.rows ) ) ] %>% mutate(mean_Females = rowMeans(.))) %>% select(scaf,base,mean_Females)

true_SNPdensity <- full_join(males_true,females_true) %>% mutate(mean_MvF_dif = mean_Males - mean_Females)
table(true_SNPdensity$mean_MvF_dif > 0)
table(true_SNPdensity$mean_MvF_dif < 0)
table(true_SNPdensity$mean_MvF_dif == 0)

# so, now we would like to permute!
# this means shuffling - we can shuffle the data columns and then use the code above.

# need an initial run of the permutation:

perm <- sample(3:length(SNPdensity.rows))

SNPdensity.rows.perm <- SNPdensity.rows %>% select(1:2,perm)

males_perm <- bind_cols(SNPdensity.rows.perm %>% select(1:2),SNPdensity.rows.perm %>% select(3:(male_n+2)) %>% mutate(mean_Males = rowMeans(.))) %>% select(scaf,base,mean_Males)

females_perm <- bind_cols(SNPdensity.rows.perm %>% select(1:2),SNPdensity.rows.perm %>% select((male_n+3):(male_n+2+female_n)) %>% mutate(mean_Females = rowMeans(.))) %>% select(scaf,base,mean_Females)

true_SNPdensity_position_perm <- full_join(males_perm,females_perm) %>% mutate(mean_MvF_dif = mean_Males - mean_Females)

perm_backbone <- true_SNPdensity_position_perm %>% select(scaf,base,mean_MvF_dif) %>% rename(p1=mean_MvF_dif)

# then we want to run the rest of the loops (since seed wasn't set, your run of this will likely produce a slightly different final file than my own on the GitHub repo, but it should not change the overall results)

#this loop takes maybe 5 minutes to run for 1000 permutations on my iMac. I've written something faster in python for the 100k permutations discussed in the MS, but this output is used for the main SexFindR plot (and the 100k permutations vs these 1000 permutations didn't provide us any new insight for fugu anyway).

for(i in 2:1000){
  perm <- sample(3:length(SNPdensity.rows))
  SNPdensity.rows.perm <- SNPdensity.rows %>% select(1:2,perm)
  males_perm <- bind_cols(SNPdensity.rows.perm %>% select(1:2),SNPdensity.rows.perm %>% select(3:(male_n+2)) %>% mutate(mean_Males = rowMeans(.))) %>% select(scaf,base,mean_Males)
  females_perm <- bind_cols(SNPdensity.rows.perm %>% select(1:2),SNPdensity.rows.perm %>% select((male_n+3):(male_n+2+female_n)) %>% mutate(mean_Females = rowMeans(.))) %>% select(scaf,base,mean_Females)
  column_ID <- paste0("p",i)
  perm_upgrade <- full_join(males_perm,females_perm) %>% mutate(mean_MvF_dif = mean_Males - mean_Females) %>% select(scaf,base,mean_MvF_dif) %>% dplyr::rename(!!column_ID := "mean_MvF_dif")
  perm_backbone <- full_join(perm_backbone,perm_upgrade)
}

perm_with_true <- full_join(true_SNPdensity %>% select(scaf,base,mean_MvF_dif),perm_backbone)

perm_with_true_p <- perm_with_true %>% 
  mutate(Pvalue = case_when(mean_MvF_dif > 0 ~ (rowSums(perm_with_true[,grep("p", names(perm_with_true))] > perm_with_true$mean_MvF_dif)+1)/(length(perm_with_true[,grep("p", names(perm_with_true))])+1),
                            mean_MvF_dif < 0 ~ (rowSums(perm_with_true[,grep("p", names(perm_with_true))] < perm_with_true$mean_MvF_dif)+1)/(length(perm_with_true[,grep("p", names(perm_with_true))])+1),
                            mean_MvF_dif == 0 ~ 1)) # here we are using a conditional mutate where if the MvF difference is > 0 we look for outliers above, if the MvF difference is < 0  we look for outliers below, and if the MvF difference is ==0, we say p value is 1.

write_tsv(perm_with_true_p,"SNPdensity_perm_with_true_p_picta.txt")

#for step 3, we only need a few columns out of this file and a subset of the genome, so we'll parse it here to save space on GitHub
write_tsv(perm_with_true_p %>% select(scaf,base,mean_MvF_dif,Pvalue), "SNPdensity_SexFindR_picta.txt")
# adding window size for consistency between analyses

# want a look at proportion of scaffold - find the outlier
scaffold_length <- read_tsv("Poecilia_picta.fna.fai",col_names = F) %>% rename(scaf=X1,length=X2)

perm_with_true_p001_snp <- perm_with_true_p %>% filter(Pvalue <= 0.001) 
count_snp_p001_scaffolds <- perm_with_true_p001_snp %>% select(scaf) %>% count(scaf)
proportion_count_snp_p001_scaffolds <- left_join(count_snp_p001_scaffolds,scaffold_length) %>% mutate(proportion = (n*10000)/length)
proportion_count_snp_p001_scaffolds %>% filter(length > 10000000) %>% ggplot(aes(x=scaf,y=proportion)) + geom_col() + theme(axis.text.x = element_text(angle = 90))


##### SCATTER PLOT

# SNP density
# SNP density with permutations, remove any that don't have p <= 0.05.
getwd()
SNP <- read_tsv("SNPdensity_SexFindR_picta.txt") %>% 
  # filter out unplaced and mitochondrial data
  filter(Pvalue <= 0.001)

# SNP density

# find points to highlight
SNPchr <- left_join(SNP, scaffold_length)

(
  site <-
    SNPchr %>%
    filter(scaf == "CM065364.1") %>%
    ggplot(aes(x = base,
               y = mean_MvF_dif)) +
    geom_point(size = 1) +
    labs(x = "", color = "") +
    labs(y = "10kb SNP density \n(male mean - female mean)", color = "") +
    labs(
      title = bquote("Density of sex-associated SNPs for P. picta"),
      subtitle = "Per window",
      x = "CM065364.1",
      color = ""
    ) +
    scale_color_manual(values = cols,
                       guide = "none") +
    scale_y_continuous(
      limits = c(
        min(SNPchr$mean_MvF_dif, na.rm = TRUE) * 1.25,
        max(SNPchr$mean_MvF_dif, na.rm = TRUE) * 1.25
      ),
      expand = c(0, 0)
    ) +
    theme_bw() +
    theme(legend.position = "none") +
    theme(
      axis.text.x = element_text(size = 10, vjust = .5),
      axis.title.x = element_text(size = 15),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(
        size = 15,
        angle = 90,
        hjust = .5,
        vjust = .5,
        face = "plain"
      ),
      title = element_text(size = 15)
    )
)
