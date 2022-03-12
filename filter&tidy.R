library(dplyr)
library(tidyr)
options(stringsAsFactors = F)
setwd("~/temp/")
datafiles <- dir(pattern = "result_021822PJX")
samples <- gsub('(.*)\\.csv', '\\1', datafiles)
data_cum <- data.frame()
data_db <- data.frame()
for(i in 1:length(datafiles)){
  data <- read.csv(datafiles[i], row.names = 1)
  
  # 18O > 2 is false positive signal
  lab18ind <- grepl("18O", data$Annotation)
  labnum <- as.numeric(gsub("18O_([0-9]+)", "\\1", data$Annotation[lab18ind]))
  lab18filt <- labnum <= 2
  filtind <- data$Annotation == "16O"
  filtind[filtind == F] <- lab18filt
  subdata <- data[filtind, ]
  # signal contain 18O but no 16O is FP
  lab18ind2 <- grepl("18O", subdata$Annotation)
  sig16 <- subdata$Lipid[!lab18ind2]
  sig18 <- subdata$Lipid[lab18ind2]
  filtind2 <- sig18 %in% sig16
  subdata18 <- subdata[which(lab18ind2)[filtind2], ]
  subdata16 <- subdata[!lab18ind2, ]
  subdata2 <- rbind(subdata16, subdata18)
  
  # Fix error signal: PC(36:4)+H
  fixind <- subdata2$Lipid %in% "PC(36:4)+H"
  subdata2$Lipid[fixind] <- paste0(subdata2$lipidClass[fixind], 
                                   subdata2$fattyAcid[fixind], "+H")
  # Tidy the data
  data_tidyi <- subdata2 %>%
    select(Lipid, Intensity, Annotation) %>%
    #group_by(Lipid) %>%
    mutate(ord = factor(Annotation, levels = unique(c('16O', Annotation)))) %>%
    arrange(Lipid, ord) %>%
    group_by(Lipid) %>%
    mutate(Intensity_norm = Intensity / Intensity[1], 
           sample = samples[i]) %>%
    select(-ord)
  data_dbi <- subset(subdata2, 
                     select = c("Lipid", "fattyAcid", "lipidClass", "lipidForm", "Theoretical_mz", "Adduct"))
  if(i == 1){
    data_cum <- data_tidyi
    data_db <- data_dbi
  }else{
    data_cum <- rbind(data_cum, data_tidyi)
    data_db <- rbind(data_db, data_dbi)
  }
}

data_tidy <- data_cum %>%
  group_by(Lipid, Annotation) %>%
  mutate(allInten = paste(Intensity, Intensity_norm, sep = '/')) %>%
  select(-c(Intensity, Intensity_norm)) %>%
  spread(sample, allInten)
for(i in 1:length(samples)){
  data_tidy <- data_tidy %>%
    separate(samples[i], sep = '/', convert = T,
             into = c(samples[i], paste0(samples[i], "_norm")))
}
dbind <- match(data_tidy$Lipid, data_db$Lipid)
data_tidy <- cbind(as.data.frame(data_tidy), 
                   subset(data_db, select = -Lipid)[dbind, ])

## Do some filter
# filter unlabeled signals greater than labeled signals
data_tidy[is.na(data_tidy)] <- 0
sigind <- 
  apply(data_tidy, 1, function(x) mean(as.numeric(c(x[3], x[5]))) <= mean(as.numeric(c(x[7], x[9]))))
data_tidy <- data_tidy[sigind, ]
# filter only have 16O signals / 18O signals
data_tidy <- data_tidy %>%
  group_by(Lipid) %>%
  filter(n() > 1)
write.csv(data_tidy, "data_tidy.csv", row.names = F)


