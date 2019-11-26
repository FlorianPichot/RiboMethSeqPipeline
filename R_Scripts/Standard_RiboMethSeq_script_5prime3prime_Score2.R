rm(list=ls(all=TRUE))

####### Input Area #######

type <- list("5.8S", "18S", "28S")		# rRNA list to study
extDataDir <-"path/to/Counts/Directory"		# Path to the samples
ResultsDir <-"path/to/Results/Directory"	# Path to the results directory (a directory is created if needed)

#PLACE Hsapiens_rRNA.csv (table with known methylated positions) in the working folder
setwd(ResultsDir)

####### PROCESSING #######

list.files(extDataDir)

DataDirs <- list.dirs(extDataDir, full.names = FALSE, recursive = FALSE)
DataDirs

ResultsDir <- paste0(ResultsDir,"/Calcul_5prime")
dir.create (ResultsDir, showWarnings=FALSE)

# Import the reference table
ref_Hs <- read.csv2("Hsapiens_rRNA.csv", sep=" ", blank.lines.skip= TRUE, header = FALSE, col.names=(c("type","position","identity")))
ref_Hs$identity <- paste0 (ref_Hs$identity, "_Hs") 
ref_Hs <- na.omit(ref_Hs)


#READING local DATA into separate dataframes


for (y in type) {
for (z in DataDirs) {
SampleDir <- paste0(ResultsDir,"/", z)
dir.create (SampleDir, showWarnings = FALSE)


results <-NA
miR_local_dir <- paste(extDataDir, "/",z, sep = "")
#print (miR_local_dir)
#print (z)
miR_file <- paste (miR_local_dir,"/", "Count_",y,"_hs_rRNA.csv", sep ="")
#print (miR_file)

#reading file from csv with names of colums and rows
sample=read.csv(miR_file, sep=" ", header = FALSE)
sample$V1 <- NULL
sample$V4 <- NULL
names(sample)[1] <- "position"
names(sample)[2] <- "counts"
assign(paste0(z),sample)
print (head(z))

sample$position <- sample$position+1    # Position shift applied to bed files where first base = pos '0'
rnames <- sample$position
row.names(sample) <- rnames		

n <- max(sample$position,na.rm = TRUE)	

table <- data.frame (position = 1:n)	# Completing positions where no 5 prime counts has been found
results <- merge (sample, table, by="position", all=TRUE)


results$counts [is.na(results$counts)] <- 1	# Change all NA counts values by 1 (otherwise, the script crashes)

results$ratio <- NA

for (i in 2:n) {results$ratio[i] = (results$counts[i]/results$counts[i-1])}		

results$around <- NA
for (i in 3:n) {
  j <- i-2
  k <- i-1 
  l <- i+1 
  m <- i+2
  results$around[i] = (1 - results$ratio[i]/(mean(results$ratio[j:k])+ mean (results$ratio[l:m])))}		# Quantification de la variabilité de notre résultat par rapport à ses 12 valeurs voisines


results$deviation <- NA
for (i in 4:n) {
  j <- i-3 
  k <- i-1 
  l <- i+1 
  m <- i+3
  results$deviation[i] = (1 - results$around[i]/(sd(results$around[j:k])+ sd(results$around[l:m])))}	# Mesure de la variance de nos résultats pour chaque position

results$ratio2 <- NA	
for (i in 2:n) {results$ratio2[i] = 1/results$ratio[i]}		# Inverse du premier ratio


results$around2 <- NA
for (i in 3:n) {
  j <- i-2
  k <- i-1 
  l <- i+1 
  m <- i+2
  results$around2[i] = (1 - results$ratio2[i]/(mean(results$ratio2[j:k])+ mean (results$ratio2[l:m])))}
  
results$cumul <- NA
for (i in 3:n) {results$cumul[i] = (results$around[i]+results$around2[i+1])/2}

results$mean <- NA
for (i in 3:n) {results$mean[i] = mean(results$around[i],results$cumul[i])}	

results$scoreA <- NA
for (i in 3:n) {
  j = i-2
  k = i-1
  l = i+1
  m = i+2
  results$scoreA[i] = 1-(2*results$counts[i]+1)/(0.5*abs(mean(results$counts[j:k])-sd(results$counts[j:k]))+results$counts[i]+0.5*abs(mean(results$counts[l:m])-sd(results$counts[l:m])))
}	

results$scoreB <- NA
for (i in 3:n) {
  e = i-2
  f = i-1
  j = i+1
  k = i+2
  results$scoreB[i] = abs(results$counts[i]-0.5*((results$counts[e]*0.9+results$counts[f]*1)/(0.9+1)+(results$counts[j]*1+results$counts[k]*0.9)/(0.9+1)))/(results$counts[i]+1)}	

results$scoreC <- NA
for (i in 3:n) {
  e = i-2
  f = i-1
  j = i+1
  k = i+2
  results$scoreC[i] = 1-results$counts[i]/(0.5*((results$counts[e]*0.9+results$counts[f]*1)/(0.9+1)+(results$counts[j]*1+results$counts[k]*0.9)/(0.9+1)))}
  
write.table (results, file=paste(SampleDir,"/Treatment_",z,"_", y,".csv",sep=""), sep = " ", row.names = FALSE)

results <- na.omit (results)

# Filtering positions with low 5 prime coverage to minimize false positives hits

mean <- mean (results$counts)
seuil <- 0.05 * mean

n <- length (results$counts)
i <- 1 
j <- 1
while (i <= n){
	if (results$counts[i] < seuil){
		j <- i
		while ((results$counts[j] < seuil) & (j < n)){
			j <- j+1
		}
		if (j-i > 4) {
			results$counts[i:j] <- NA
			i <- j+1
			print (paste0("DELETE_", z, "_", y))
		}
		else {
			i <- i+1
	}	}
	else {
	i <- i+1
	}
} 

results <- na.omit(results)

results_mean <- data.frame (position = results$position, mean = results$mean)
results_mean <- results_mean [order (results_mean$mean, decreasing = TRUE),]
results_scoreA <- data.frame (position = results$position, scoreA = results$scoreA)
results_scoreA <- results_scoreA [order (results_scoreA$scoreA, decreasing = TRUE),]
results_scoreB <- data.frame (position = results$position, scoreB = results$scoreB)
results_scoreB <- results_scoreB [order (results_scoreB$scoreB, decreasing = TRUE),]
results_scoreC <- data.frame (position = results$position, scoreC = results$scoreC)
results_scoreC <- results_scoreC [order (results_scoreC$scoreC, decreasing = TRUE),]

ref_Hs_type <- ref_Hs[grep (paste0('^',y), ref_Hs$type), ]
score_mean <- merge (results_mean, ref_Hs_type, by="position", all = TRUE)
score_mean <- score_mean[order(score_mean$mean, decreasing = TRUE),]
scoreA <- merge (results_scoreA, ref_Hs_type, by="position", all = TRUE)
scoreA <- scoreA[order(scoreA$scoreA, decreasing = TRUE),]
scoreC <- merge (results_scoreC, ref_Hs_type, by="position")

names (results_mean) [2] <- paste0 ("mean.", z)
names (results_scoreA) [2] <- paste0 ("ScoreA.", z)

write.table (score_mean, file = paste0(SampleDir, "/Score_mean_", y, ".csv"), sep = " ", row.names=FALSE)
write.table (scoreA, file = paste0(SampleDir, "/Score_A_", y, ".csv"), sep = " ", row.names=FALSE)
write.table (scoreC, file = paste0(SampleDir, "/Score_C_", y, ".csv"), sep = " ", row.names=FALSE)

n = length (score_mean$position)
score_mean$mean_rank <- 1:n
score_mean$type <- NULL
score_mean$identity <- NULL
scoreA$scoreA_rank <- 1:n
scoreA$type <- NULL
scoreA$identity <- NULL


graphe <- merge (score_mean, scoreA, by = "position")
graphe <- merge (graphe, ref_Hs_type, by= "position", all = TRUE)
Nm <- graphe[grep("[AGUCP]m", graphe$identity),]
Psi <- graphe[grep ("Psi", graphe$identity),]
MetBase <- graphe[grep ("^m",graphe$identity),]
N <- graphe[is.na(graphe$identity),]

setwd(SampleDir)

library(ggplot2)

windows()
plotmean1 <- ggplot (graphe, aes(mean)) + xlim(0,1) + geom_histogram(alpha = 0.3, fill = "green", bins = 250) + theme_bw()
ggsave (filename = paste0 ("Histogramme_global_score_mean_2_", y, "_", z, ".pdf"), plot = plotmean1)

windows()
plotmean2 <- ggplot (Nm, aes(mean)) + xlim(0,1) + geom_histogram(alpha = 0.75, fill = "deepskyblue", bins = 250) + geom_histogram(data = N, alpha = 0.5, fill = "grey", bins = 250) + geom_histogram(data = Psi, alpha = 0.5, fill = "red", bins = 250) + theme_bw()
ggsave (filename = paste0 ("Histogramme_colour_score_mean_2_", y, "_", z, ".pdf"), plot = plotmean2)

windows()
plotA1 <- ggplot (graphe, aes(scoreA)) + xlim(0,1) + geom_histogram(alpha = 0.3, fill = "green", bins = 250) + theme_bw()
ggsave (filename = paste0 ("Histogramme_global_score_A_2_", y, "_", z, ".pdf"), plot = plotA1)

windows()
plotA2 <- ggplot (Nm, aes(scoreA)) + xlim(0,1) + geom_histogram(alpha = 0.75, fill = "deepskyblue", bins = 250) + geom_histogram(data = N, alpha = 0.5, fill = "grey", bins = 250) + geom_histogram(data = Psi, alpha = 0.5, fill = "red", bins = 250) + theme_bw()
ggsave (filename = paste0 ("Histogramme_colour_score_A_2_", y, "_", z, ".pdf"), plot = plotA2)


graphics.off()


library(ggplot2)
library(gridExtra)

ggplot (graphe, aes(mean)) + 
xlim(0,1) + 
geom_density(data = Nm, alpha = 0.75, fill = "deepskyblue", colour="deepskyblue") + 
geom_density(data = N, alpha = 0.5, fill = "grey", colour="grey") + 
geom_density(data = Psi, alpha = 0.5, fill = "red", colour="red") +
theme_bw()
ggsave (filename = paste0 ("Density_score_mean_2_", y, "_", z, ".pdf"), plot = last_plot())


ggplot (Nm, aes(scoreA)) + 
xlim(0,1) + 
geom_density(alpha = 0.75, fill = "deepskyblue", colour = "deepskyblue") + 
geom_density(data = N, alpha = 0.5, fill = "grey", colour = "grey") + 
geom_density(data = Psi, alpha = 0.5, fill = "red", colour = "red") +
theme_bw()
ggsave (filename = paste0 ("Density_score_score_A_2_", y, "_", z, ".pdf"), plot = last_plot())

theme0 <- function(...) theme( legend.position = "none",
                               panel.background = element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.margin = unit(0,"null"),
                               axis.ticks = element_blank(),
                               axis.text.x = element_blank(),
                               axis.text.y = element_blank(),
                               axis.title.x = element_blank(),
                               axis.title.y = element_blank(),
                               axis.ticks.length = unit(0,"null"),
                               axis.ticks.margin = unit(0,"null"),
                               panel.border=element_rect(color=NA),...)


p1 <- ggplot (Nm, aes(x = mean_rank, y = scoreA_rank)) + 
geom_point(alpha = 0.75, colour = "deepskyblue") + 
geom_point(data = N, alpha = 0.5, colour = "grey", size = 0.2) + 
geom_point(data = Psi, alpha = 0.5, colour = "red") +
scale_x_continuous(expand=c(0.02,0)) +
scale_y_continuous(expand=c(0.02,0)) +
theme_bw()

p2 <- ggplot (Nm, aes(x = mean_rank)) + 
geom_density(alpha = 0.75, fill = "deepskyblue", colour = "deepskyblue") + 
geom_density(data = N, alpha = 0.5, fill = "grey", colour="grey") + 
geom_density(data = Psi, alpha = 0.5, fill = "red", colour="red") +
theme_bw() +
theme0(plot.margin = unit(c(0,-0.3,-0.4,2.8),"lines")) 

p3 <- ggplot (Nm, aes(x = scoreA_rank)) + 
geom_density(alpha = 0.75, fill = "deepskyblue", colour = "deepskyblue") + 
geom_density(data = N, alpha = 0.5, fill = "grey", colour = "grey") + 
geom_density(data = Psi, alpha = 0.5, fill = "red", colour = "red") +
coord_flip() +
theme_bw() +
theme0(plot.margin = unit(c(-0.3,0,1.9,0),"lines"))

grid.arrange(arrangeGrob(p2,ncol=2,widths=c(3,1)),
             arrangeGrob(p1,p3,ncol=2,widths=c(3,1)),
             heights=c(1,3))

ggsave (filename = paste0 ("Density+plot_rank_2_", y, "_", z, ".pdf"), plot = grid.arrange(arrangeGrob(p2,ncol=2,widths=c(3,1)),
             arrangeGrob(p1,p3,ncol=2,widths=c(3,1)),
             heights=c(1,3)))



p1 <- ggplot (Nm, aes(x = scoreA, y = mean)) + 
xlim(0,1) +
ylim(0,1) +
geom_point(alpha = 0.75, colour = "deepskyblue") + 
geom_point(data = N, alpha = 0.5, colour = "grey", size = 0.2) + 
geom_point(data = Psi, alpha = 0.5, colour = "red") +
theme_bw()

p2 <- ggplot (Nm, aes(scoreA)) + 
xlim(0,1) + 
geom_density(alpha = 0.75, fill = "deepskyblue", colour = "deepskyblue") + 
geom_density(data = N, alpha = 0.5, fill = "grey", colour = "grey") + 
geom_density(data = Psi, alpha = 0.5, fill = "red", colour = "red") +
theme_bw() +
theme0(plot.margin = unit(c(0,0.4,-0.4,3.3),"lines")) 

p3 <- ggplot (graphe, aes(mean)) + 
xlim(0,1) + 
coord_flip() +
geom_density(data = Nm, alpha = 0.75, fill = "deepskyblue", colour="deepskyblue") + 
geom_density(data = N, alpha = 0.5, fill = "grey", colour="grey") + 
geom_density(data = Psi, alpha = 0.5, fill = "red", colour="red") +
theme_bw() +
theme0(plot.margin = unit(c(0.4,0,2.6,0),"lines"))

grid.arrange(arrangeGrob(p2,ncol=2,widths=c(3,1)),
             arrangeGrob(p1,p3,ncol=2,widths=c(3,1)),
             heights=c(1,3))

ggsave (filename = paste0 ("Density+plot_2_score_", y, "_", z, ".pdf"), plot = grid.arrange(arrangeGrob(p2,ncol=2,widths=c(3,1)),
             arrangeGrob(p1,p3,ncol=2,widths=c(3,1)),
             heights=c(1,3)))

graphics.off()

}


}
