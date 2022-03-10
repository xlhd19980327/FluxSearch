require(fuzzyjoin)
require(data.table)
require(dplyr)


resultWrap <- function(inputFile, refInput, ion_mode, nC = 40){
  inputPeaklist <- fread(inputFile)
  reference_list <- ref_13C(refInput, ion_mode, nC = nC)
  reference <- reference_list[[1]]
  reference_mz <- reference_list[[2]]
  
  # enviPick result info
  mz1 <- as.numeric(inputPeaklist$"m/z")
  rt1 <- (as.numeric(inputPeaklist$"RT")/60)
  Intensity <- inputPeaklist$"max_int"
  mzList <- data.frame(mz1, rt1, Intensity)
  
  result <- data.frame()
  for(i in 0:nC){
    reference_calc <- reference
    reference_calc$mz2 <- reference_mz[[i+1]] 
    reference_calc$Annotation <- names(reference_mz)[i+1]
    result <- rbind(result, 
                    msMatch(mzList, reference_calc))
  }
  return(result)
}

msMatch <- function(mzFile, refFile) {
  resultFile <- fuzzy_join(
    mzFile,
    refFile,
    multi_by = c("mz1" = "mz2", "rt1" = "rt2"),
    multi_match_fun = mmf,
    mode = "left"
  )
  
  if(length(resultFile$Score) != 0) {
    res_order <- resultFile[order(-resultFile$Score,
                                  resultFile$delmz_ppm), ]
    #resultReducedFile <- na.omit(res_order)
    uniResult <- res_order[!duplicated(res_order[, c("fattyAcid", "lipidForm", "Theoretical_mz")]), ]
    return(uniResult)
  } else {
    Score <- NA
    delmz_ppm <- NA
    uniResult_na <- data.frame(resultFile, Score, delmz_ppm)
    return(uniResult_na)
  }
  
}


ref_13C <- function(refInput, ion_mode, nC){
  # O16: 15.99491461957, O18: 17.99915961287
  Odist <- 2.0042449933
  ref_list <- refLib_ls(refInput, ion_mode)
  ref_13Cmzlist <- lapply(0:nC, function(x) ref_list$mz2 + (x*(Odist/ref_list$z)))
  #length(ref_13Clist) <- nC + 1
  names(ref_13Cmzlist) <- c("16O", 
                            paste0("18O_", 
                                   sapply(1:nC, function(x) ifelse(x < 10, paste0("0", x), as.character(x))))
  )
  return(list(ref_list, ref_13Cmzlist))
}

refLib_ls <- function(refInput, ion_mode){
  refPeaklist <- fread(refInput)
  # extract the negative peaks for reference negative files
  refPeaklist_neg <- subset(refPeaklist, refPeaklist$Pol == switch(ion_mode, pos = "P", neg = "N"))
  # re-organize the refPeaklist_neg
  mz2 <- refPeaklist_neg$ObsMz
  rt2 <- refPeaklist_neg$Rt
  Lipid <- refPeaklist_neg$LipidIon
  fattyAcid <- refPeaklist_neg$FattyAcid
  lipidClass <- refPeaklist_neg$Class
  lipidForm <- refPeaklist_neg$Formula
  Theoretical_mz <- refPeaklist_neg$CalcMz
  Adduct <- refPeaklist_neg$Ion
  Charge <- refPeaklist_neg$Pol
  z <- refPeaklist_neg$z
  Grade_ls <- refPeaklist_neg$Grade
  
  ref <- data.frame(mz2, rt2, Lipid, fattyAcid, lipidClass, lipidForm,
                    Theoretical_mz, Adduct, Charge, z, Grade_ls)
  ref$mz2 <- as.numeric(ref$mz2)
  ref$rt2 <- as.numeric(ref$rt2)
  ref$z <- as.numeric(ref$z)
  return(ref)
}

mmf <- function(x, y) {
  # The differnce between data vs. reference
  delta_mz <- abs(x[,1] - y[,1])
  mz_dist <- 1000000*abs(x[,1] - y[,1])/abs(x[,1])
  rt_dist <- abs(x[,2] - y[,2])
  out <- data_frame(merge = rt_dist <= 2 & mz_dist <= 10,
                    Score = 1 - sqrt((mz_dist*0.1)^2 + (rt_dist)^2),
                    delmz_ppm = mz_dist)
  return(out)
}


#===================================================================================================
# Wrap all positive and negative files together

Flux_result <- function(input_negative, input_positive, referInput, score = 0.8) {
  flux_neg <- resultWrap(input_negative, referInput, ion_mode = "neg")
  flux_pos <- resultWrap(input_positive, referInput, ion_mode = "pos")
  output_all <- rbind(flux_neg, flux_pos)
  output_all <- Sgrade(output_all)
  # select feature with good grades
  output_all <- subset(output_all, Score >= score)
  
  return(output_all)
}


# Grading function
Sgrade <- function(inputResult) {
  inputResult$Grades <- cut(inputResult$Score,
                            breaks = c(-2, 0, 0.5, 0.7, 0.8, 0.9, 1),
                            labels = c("F", "D", "C", "B", "A", "A+"))
  return(inputResult)
}
