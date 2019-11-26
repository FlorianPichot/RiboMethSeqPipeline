#IDEM (if required for UCounts name of final file converted to Count_)


#SCRIPT TO READ 5'counts and 3'counts and create UCount_RNA file

extDataDir <- Path/to/Samples
setwd(extDataDir)

DataDirs <- list.dirs(extDataDir, full.names = TRUE, recursive = FALSE)
DataDirs

#READING local DATA into separate dataframes

for (z in DataDirs) {
list5<-list.files(z, pattern ="^UCount5prime", full.names = TRUE, recursive = TRUE)  #TAKE only UCount 
RNA_list <- substr (list5,unlist(regexpr( '5prime_', list5)) + 7, unlist(regexpr('.csv', list5)) - 1)
type <- as.list (RNA_list)

for (y in type) {
count5file <- list.files(z, pattern = paste0("^UCount5prime_",y,".csv"), full.names = TRUE, recursive = FALSE)
count3file <- list.files(z, pattern = paste0("^UCount3prime_",y,".csv"), full.names = TRUE, recursive = FALSE)

sample5r<-read.csv(file = count5file, sep=" ", header = FALSE)
sample3r<-read.csv(file = count3file, sep=" ", header = FALSE)

sample5r$V4 <- NULL
sample5r$V1 <- NULL
sample3r$V4 <- NULL
sample3r$V1 <- NULL

sample5r<-na.omit(sample5r)						#remove evantual NA values
sample3r<-na.omit(sample3r) 						#remove evantual NA values

names(sample5r)[1] <- "position"
names(sample5r)[2] <- "counts5"
names(sample3r)[1] <- "position"
names(sample3r)[2] <- "counts3"

sample <-merge(sample5r, sample3r, "position", all = TRUE)
sample[is.na(sample)] <- 0
sample$counts<-sample$counts5+sample$counts3
sample$position<-sample$position+1

sample$counts5<-NULL
sample$counts3<-NULL

sample$x<-""
sample$y<-""
sample <- sample[c(3,1,2,4)]    
write.table (sample, file = paste0(z, "/Count_", y, ".csv"), quote =FALSE, sep = " ", row.names = FALSE, col.names=FALSE)
}
}
