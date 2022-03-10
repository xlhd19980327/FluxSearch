library(enviPick)
drtgap <- 300
dmzdens <- 20
minpeak <- 3
drtsmall2 <- 20  ##20 
drtfill <- 10
drtdens2 <- 300  ##300 
minint <- 3  ##use 3 before
SN <- 3  
SB <- 1
recurs <- 3
ended <- 1
weight <- 1
maxint <- 6 ##use 6 before
ion_mode <- c("positive", "negative") #switch ion mode: c("positive", "negative")
MSlevel <- 1
dmzgap <-  dmzdens * 2 + 1

args=commandArgs(T)
xmlin <- args[1]
xmlout <- args[2]

## For Windows path
#xmlin <- gsub("\\", "/", xmlin, fixed = T)
#xmlout <- gsub("\\", "/", xmlout, fixed = T)
#print(xmlin)

## Modify the NEG/POS when switching ion mode
mzXML <- dir(xmlin, 
             pattern = ".*?\\.mzXML")
#csvFile <- gsub("(.*?)\\.mzXML", "\\1", mzXML)
for(k in ion_mode){
	for(i in mzXML){
		MSlist<-readMSdata(paste0(xmlin, "/", i),
							MSlevel=c(MSlevel), ion_mode = k)
		MSlist<-mzagglom(MSlist,dmzgap=dmzgap,ppm=TRUE,drtgap=drtgap,minpeak=minpeak,maxint=10^maxint)
		MSlist<-mzclust(MSlist,dmzdens=dmzdens,ppm=TRUE,drtdens=60,minpeak=minpeak,maxint=10^maxint)
		MSlist<-mzpick(MSlist, minpeak = minpeak, drtsmall = drtsmall2, drtfill = drtfill, 
						drttotal = drtdens2, recurs = recurs, weight = weight, SB = SB, SN=SN, minint = 10^minint, maxint = 10^maxint, ended = ended)
		#writePeaklist(MSlist,"C:/Users/lfs/Desktop/myLearning/lipGroup/riverGroup/practice/100919QltmData/NEG",
		#              "NEG-13_1.txt")
		j <- gsub("(.*?)\\.mzXML", "\\1", i)
		write.csv(MSlist[[8]], 
					paste0(xmlout, "/", k, "/", 
						j, ".csv"))
		print(paste0(k, " is ok!\n"))
	}
}




