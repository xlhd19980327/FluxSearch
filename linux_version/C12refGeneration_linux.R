args=commandArgs(T)
wd <- args[1]
out <- args[2]

## For Windows path
#wd <- gsub("\\", "/", wd, fixed = T)
#out <- gsub("\\", "/", out, fixed = T)
#print(xmlin)

setwd(wd)
ref12CFile <- dir("./", 
                  pattern = ".*?\\.txt")
ref12C <- data.frame(stringsAsFactors = F)
for (i in ref12CFile){
  file <- read.table(paste0("./", i), header = T, sep = '\t', skip = 5)
  file <- subset(file, 
                 select = c("LipidIon", "Class", "FattyAcid", "Ion", 
                            "Formula", "ObsMz", "Rt", "Pol.", "z", "ObsMass", 
                            "CalcMz", "Grade"))
  ref12C <- rbind(ref12C, file)
}
write.csv(ref12C, paste0(out, "/ref_12C.csv"))
