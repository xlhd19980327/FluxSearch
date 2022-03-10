#source("/home/lifs/myWork/mycode/13Chandle/lipidFluxSearch_functions_0.5.1.R")
source("/home/lifs/myWork/mycode/13Chandle/lipidFluxSearch_lfsRevised.R")

args=commandArgs(T)
xmlfiles <- args[1]
C12files <- args[2]
resultout <- args[3]

#pos <- dir(paste0(xmlfiles, "/positive/"), pattern = ".*?\\.csv")
### In this case, neg and pos reaults should have the same name, or this may fail
### Note: score = 0.5
neg <- dir(paste0(xmlfiles, "/negative/"), pattern = ".*?\\.csv")
for(i in neg){
  print(paste0(i, " starts: ", date()))
  lipids <- Flux_result(paste0(xmlfiles, "/negative/", i), 
                        paste0(xmlfiles, "/positive/", i), 
                        paste0(C12files, "/ref_12C.csv"), score=0.5)
  write.csv(lipids, paste0(resultout, "/result_", i))
  print(paste0(i, " ends: ", date()))
}
