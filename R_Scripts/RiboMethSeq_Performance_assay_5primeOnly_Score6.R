####### Input Area #######

type <- list("5.8S", "18S", "28S")		      # rRNA list to study
extDataDir <-"Path/to/Samples/Directory"		# Path to the directory with the 'counts' files
ResultsDir <-"Path/to/Results/Directory"		# Path to the results directory (A directory is created if needed)

#PLACE Hsapiens_rRNA.csv in the working folder
setwd(ResultsDir)

## ROCCurve Functions ##
#FUNCTIONS TO CALCULATE THE nb OF TRUE POSITIVES
fctTP <- function (identity) {
  i <- ifelse (grepl("[AGCUP]m", identity), 1, 0)
  return (cumsum(i))
}
#FUNCTIONS TO CALCULATE THE nb OF FALSE POSITIVES
fctFP <- function (identity) {
  i <- ifelse (grepl("[AGCUP]m", identity), 0, 1)
  return (cumsum(i))
}
#FUNCTIONS TO CALCULATE THE nb OF FALSE NEGATIVES
fctFN <- function (identity, TP) {
  i <- ifelse (grepl("[AGCUP]m", identity), 1, 0)
  return (sum(i) - TP)
}
#FUNCTIONS TO CALCULATE THE nb OF TRUE NEGATIVES
fctTN <- function (identity, FP) {
  i <- ifelse (grepl("[AGCUP]m", identity), 0, 1)
  return (sum(i) - FP)
}
#FUNCTIONS TO CALCULATE THE Matthews Correlation Coefficient (MCC)
fctMCC <- function (TP, TN, FP, FN) {
  return (((TP*TN)-(FP*FN))/(sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))))
}
#FUNCTIONS TO CALCULATE THE True Positive Rate (TPR)
fctTPR <- function (TP, FN) {
  return (TP/(TP+FN))
}
#FUNCTIONS TO CALCULATE THE True Negative Rate (TNR)
fctTNR <- function (TN, FP) {
  return (TN/(TN+FP))
}
#FUNCTIONS TO CALCULATE THE Positive Predictive Value (PPV)
fctPPV <- function (TP, FP) {
  return (TP/(TP+FP))
}
#FUNCTIONS TO CALCULATE THE Negative Predictive Value (NPV)
fctNPV <- function (TN, FN) {
  return (TN/(TN+FN))
}
#FUNCTIONS TO CALCULATE THE False Negative Rate (FNR)
fctFNR <- function (FN, TP) {
  return (FN/(FN+TP))
}
#FUNCTIONS TO CALCULATE THE False Positive Rate (FPR)
fctFPR <- function (TN, FP) {
  return (FP/(FP+TN))
}
#FUNCTIONS TO CALCULATE THE False Discovery Rate (FDR)
fctFDR <- function (TP, FP) {
  return (FP/(FP+TP))
}
#FUNCTIONS TO CALCULATE THE False Omission Rate (FOR)
fctFOR <- function (TN, FN) {
  return (FN/(FN+TN))
}
#FUNCTIONS TO CALCULATE THE Accuracy (ACC)
fctACC <- function (TP, TN, FP, FN) {
  return ((TP+TN)/(TP+TN+FP+FN))
}


####### PROCESSING #######
list.files(extDataDir)

DataDirs <- list.dirs(extDataDir, full.names = FALSE, recursive = FALSE)

ResultsDir <- paste0(ResultsDir,"/2DensityPlot(ScoreMAX6)_MCC_", basename(extDataDir))
dir.create (ResultsDir, showWarnings=FALSE)

ref_Hs <- read.csv2("Hsapiens_rRNA.csv", sep=" ", blank.lines.skip= TRUE, header = FALSE, col.names=(c("type","position","identity")))
ref_Hs$identity <- paste0 (ref_Hs$identity, "_Hs") 
ref_Hs <- na.omit(ref_Hs)


#READING local DATA into separate dataframes

memory <- data.frame ()
total <- data.frame ()

for (z in DataDirs) {
  
  tot_score_max <- data.frame()
  tot_scoreA <- data.frame()
  
  for (y in type) {
    SampleDir <- paste0(ResultsDir,"/", z)
    dir.create (SampleDir, showWarnings = FALSE)
    
    
    results <-NA
    miR_local_dir <- paste(extDataDir, "/",z, sep = "")
    #print (miR_local_dir)
    #print (z)
    miR_file <- paste (miR_local_dir,"/", "UCount5prime_",y,"_hs_rRNA.csv", sep ="")
    #print (miR_file)
    
    #reading file from csv with names of colums and rows
    sample=read.csv(miR_file, sep=" ", header = FALSE)
    sample$V1 <- NULL
    sample$V4 <- NULL
    names(sample)[1] <- "position"
    names(sample)[2] <- "counts"
    assign(paste0(z),sample)
    print (head(z))
    
    rnames <- sample$position
    row.names(sample) <- rnames
    
    n <- max(sample$position,na.rm = TRUE)	
    
    table <- data.frame (position = 1:n)	
    results <- merge (sample, table, by="position", all=TRUE)	
    
    results$position <- results$position+1
    results$counts [is.na(results$counts)] <- 1	  # Fill NA values (otherwise, the script crashes)
    
   
    results$ratio <- NA	
    
    for (i in 2:n) {results$ratio[i] = (results$counts[i]/results$counts[i-1])}
    
    results$around <- NA
    for (i in 7:n) {
      j <- i-6
      k <- i-1 
      l <- i+1 
      m <- i+6
      results$around[i] = (1 - results$ratio[i]/(mean(results$ratio[j:k])+ mean (results$ratio[l:m])))}
    
    
    results$deviation <- NA
    for (i in 7:n) {
      j <- i-3 
      k <- i-1 
      l <- i+1 
      m <- i+3
      results$deviation[i] = (1 - results$around[i]/(sd(results$around[j:k])+ sd(results$around[l:m])))}
    
    results$ratio2 <- NA	
    for (i in 2:n) {results$ratio2[i] = 1/results$ratio[i]}	
    
    
    results$around2 <- NA
    for (i in 7:n) {
      j <- i-6
      k <- i-1 
      l <- i+1 
      m <- i+6
      results$around2[i] = (1 - results$ratio2[i]/(mean(results$ratio2[j:k])+ mean (results$ratio2[l:m])))}	
    
    results$cumul <- NA
    for (i in 7:n) {results$cumul[i] = (results$around[i]+results$around2[i+1])/2}
    
    results$max <- NA
    for (i in 7:n) {results$max[i] = max(results$around[i],results$cumul[i])}
    
    results$scoreA <- NA
    for (i in 7:n) {j = i-6
    k = i-1
    l = i+1
    m = i+6
    results$scoreA[i] = 1-(2*results$counts[i]+1)/(0.5*abs(mean(results$counts[j:k])-sd(results$counts[j:k]))+results$counts[i]+0.5*abs(mean(results$counts[l:m])-sd(results$counts[l:m])))
    }
    
    results$scoreB <- NA
    for (i in 7:n) {a = i-6
    b = i-5
    c = i-4
    d = i-3
    e = i-2
    f = i-1
    j = i+1
    k = i+2
    l = i+3
    m = i+4
    o = i+5
    p = i+6
    results$scoreB[i] = abs(results$counts[i]-0.5*((results$counts[a]*0.5+results$counts[b]*0.6+results$counts[c]*0.7+results$counts[d]*0.8+results$counts[e]*0.9+results$counts[f]*1)/(0.5+0.6+0.7+0.8+0.9+1)+(results$counts[j]*1+results$counts[k]*0.9+results$counts[l]*0.8+results$counts[m]*0.7+results$counts[o]*0.6+results$counts[p]*0.5)/(0.5+0.6+0.7+0.8+0.9+1)))/(results$counts[i]+1)}
    
    results$scoreC <- NA
    for (i in 7:n) {a = i-6
    b = i-5
    c = i-4
    d = i-3
    e = i-2
    f = i-1
    j = i+1
    k = i+2
    l = i+3
    m = i+4
    o = i+5
    p = i-6
    results$scoreC[i] = 1-results$counts[i]/(0.5*((results$counts[a]*0.5+results$counts[b]*0.6+results$counts[c]*0.7+results$counts[d]*0.8+results$counts[e]*0.9+results$counts[f]*1)/(0.5+0.6+0.7+0.8+0.9+1)+(results$counts[j]*1+results$counts[k]*0.9+results$counts[l]*0.8+results$counts[m]*0.7+results$counts[o]*0.6+results$counts[p]*0.5)/(0.5+0.6+0.7+0.8+0.9+1)))}
    
    write.table (results, file=paste(SampleDir,"/Treatment_",z,"_", y,".csv",sep=""), sep = " ", row.names = FALSE)
    
    
    
    results <- na.omit (results)
    
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
    
    results_max <- data.frame (position = results$position, max = results$max)
    results_max <- results_max [order (results_max$max, decreasing = TRUE),]
    results_scoreA <- data.frame (position = results$position, scoreA = results$scoreA)
    results_scoreA <- results_scoreA [order (results_scoreA$scoreA, decreasing = TRUE),]
    
    
    ref_Hs_type <- ref_Hs[grep (paste0('^',y), ref_Hs$type), ]
    score_max <- merge (results_max, ref_Hs_type, by="position", all = TRUE)
    score_max$identity [is.na(score_max$identity)] <- "No_Met"
    score_max <- score_max[order(score_max$max, decreasing = TRUE),]
    scoreA <- merge (results_scoreA, ref_Hs_type, by="position", all = TRUE)
    scoreA$identity [is.na(scoreA$identity)] <- "No_Met"
    scoreA <- scoreA[order(scoreA$scoreA, decreasing = TRUE),]

    write.table (score_max, file = paste0(SampleDir, "/Score_max_", y, ".csv"), sep = " ", row.names=FALSE)
    write.table (scoreA, file = paste0(SampleDir, "/Score_A_", y, ".csv"), sep = " ", row.names=FALSE)
    
    fusion <- merge(scoreA, score_max, by = c("position","type","identity"))
    fusion$type <- y
    
    memory <- rbind (memory,fusion)
    tot_score_max <- rbind(tot_score_max, fusion)
    tot_scoreA <- rbind(tot_scoreA, fusion)
    
  }
  ###### ROC Curves calculation (with functions) ######
  
  ###### ROC curve ######
  # max #
  tot_score_max <- tot_score_max[order(tot_score_max$max, decreasing = TRUE),]
  
  
  tot_score_max$TP <- fctTP(identity = tot_score_max$identity)
  tot_score_max$FP <- fctFP(identity = tot_score_max$identity)
  tot_score_max$FN <- fctFN(identity = tot_score_max$identity, TP = tot_score_max$TP)
  tot_score_max$TN <- fctTN(identity = tot_score_max$identity, FP = tot_score_max$FP)
  tot_score_max$MCC <- fctMCC(TP = tot_score_max$TP, TN = tot_score_max$TN, FP = tot_score_max$FP, FN = tot_score_max$FN)
  tot_score_max$TPR <- fctTPR(TP = tot_score_max$TP, FN = tot_score_max$FN)
  tot_score_max$TNR <- fctTNR(TN = tot_score_max$TN, FP = tot_score_max$FP)
  tot_score_max$PPV <- fctPPV(TP = tot_score_max$TP, FP = tot_score_max$FP)
  tot_score_max$NPV <- fctNPV(TN = tot_score_max$TN, FN = tot_score_max$FN)
  tot_score_max$FNR <- fctFNR(TP = tot_score_max$TP, FN = tot_score_max$FN)
  tot_score_max$FPR <- fctFPR(TN = tot_score_max$TN, FP = tot_score_max$FP)
  tot_score_max$FDR <- fctFDR(TP = tot_score_max$TP, FP = tot_score_max$FP)
  tot_score_max$FOR <- fctFOR(TN = tot_score_max$TN, FN = tot_score_max$FN)
  tot_score_max$ACC <- fctACC(TP = tot_score_max$TP, TN = tot_score_max$TN, FP = tot_score_max$FP, FN = tot_score_max$FN)
  
  write.table (tot_score_max, file=paste0(SampleDir, "/ROC_max.csv"), sep = " ", row.names = FALSE)  
  
  # ROCurve score max #
  dev.new()
  plot (tot_score_max$FPR, tot_score_max$TPR, type = "l", col = "blue", xlim = c(0,0.05), main = "ROCurve for Score max", sub=paste0("MCC = ", max(unlist(tot_score_max$MCC))), xlab=" ", ylab=" ")
  lines(tot_score_max$FPR, tot_score_max$MCC, type = "l", col = "orange")
  savePlot (filename = paste0 (SampleDir, "/ROCurve_Score_max"), "pdf")
  graphics.off()
  
  # scoreA #
  tot_scoreA <- tot_scoreA[order(tot_scoreA$scoreA, decreasing = TRUE),]
  
  tot_scoreA$TP <- fctTP(identity = tot_scoreA$identity)
  tot_scoreA$FP <- fctFP(identity = tot_scoreA$identity)
  tot_scoreA$FN <- fctFN(identity = tot_scoreA$identity, TP = tot_scoreA$TP)
  tot_scoreA$TN <- fctTN(identity = tot_scoreA$identity, FP = tot_scoreA$FP)
  tot_scoreA$MCC <- fctMCC(TP = tot_scoreA$TP, TN = tot_scoreA$TN, FP = tot_scoreA$FP, FN = tot_scoreA$FN)
  tot_scoreA$TPR <- fctTPR(TP = tot_scoreA$TP, FN = tot_scoreA$FN)
  tot_scoreA$TNR <- fctTNR(TN = tot_scoreA$TN, FP = tot_scoreA$FP)
  tot_scoreA$PPV <- fctPPV(TP = tot_scoreA$TP, FP = tot_scoreA$FP)
  tot_scoreA$NPV <- fctNPV(TN = tot_scoreA$TN, FN = tot_scoreA$FN)
  tot_scoreA$FNR <- fctFNR(TP = tot_scoreA$TP, FN = tot_scoreA$FN)
  tot_scoreA$FPR <- fctFPR(TN = tot_scoreA$TN, FP = tot_scoreA$FP)
  tot_scoreA$FDR <- fctFDR(TP = tot_scoreA$TP, FP = tot_scoreA$FP)
  tot_scoreA$FOR <- fctFOR(TN = tot_scoreA$TN, FN = tot_scoreA$FN)
  tot_scoreA$ACC <- fctACC(TP = tot_scoreA$TP, TN = tot_scoreA$TN, FP = tot_scoreA$FP, FN = tot_scoreA$FN)
  
  write.table (scoreA, file=paste0(SampleDir, "/ROC_scoreA.csv"), sep = " ", row.names = FALSE)
  
  # ROCurve score A #
  dev.new()
  plot (tot_score_A$FPR, tot_score_A$TPR, type = "l", col = "blue", xlim = c(0,0.05), main = "ROCurve for Score A", sub=paste0("MCC = ", max(unlist(tot_score_A$MCC))), xlab=" ", ylab=" ")
  lines(tot_score_A$FPR, tot_score_A$MCC, type = "l", col = "orange")
  savePlot (filename = paste0 (SampleDir, "/ROCurve_Score_A"), "pdf")
  graphics.off()
  
  # Compilation of MCC results #
  
  MCC_scoreA <- tot_scoreA[grep(max(tot_scoreA$MCC, na.rm = TRUE),tot_scoreA$MCC),]
  MCC_score_max <- tot_score_max[grep(max(tot_score_max$MCC, na.rm = TRUE),tot_score_max$MCC),]
  
  MCC_total <-cbind(z, MCC_score_max$MCC, MCC_score_max$TPR, MCC_score_max$PPV, MCC_score_max$FDR, MCC_scoreA$MCC, MCC_scoreA$TPR, MCC_scoreA$PPV, MCC_scoreA$FDR)
  colnames (MCC_total) <- c("sample", "MCC_max", "TPR_max", "PPV_max", "FDR_max", "MCC_scoreA", "TPR_scoreA", "PPV_scoreA", "FDR_scoreA")
  total <- rbind (total, MCC_total)
  
}

write.table (total, file=paste0(ResultsDir, "/ROC_total_MAX6_5prime.csv"), sep = " ", row.names = FALSE)

n = length (memory$position)

memory <- memory[order(memory$max, decreasing = TRUE),]
memory$max_rank <- 1:n
memory <- memory[order(memory$scoreA, decreasing = TRUE),]
memory$scoreA_rank <- 1:n
memory <- memory[order(memory$position, decreasing = FALSE),]
memory <- memory[order(memory$type, decreasing = FALSE),]

graphe <- memory
Nm <- graphe[grep("[AGUCP]m", graphe$identity),]
Psi <- graphe[grep ("Psi", graphe$identity),]
MetBase <- graphe[grep ("^m",graphe$identity),]
N <- graphe[grep ("No_Met", graphe$identity),]

setwd(ResultsDir)

library(ggplot2)

windows()
plotmax1 <- ggplot (graphe, aes(max)) + xlim(0,1) + geom_histogram(alpha = 0.3, fill = "green", bins = 250) + theme_bw()
ggsave (filename = paste0 ("Histogramme_global_score_max_6.pdf"), plot = plotmax1)

windows()
plotmax2 <- ggplot (Nm, aes(max)) + xlim(0,1) + geom_histogram(alpha = 0.75, fill = "deepskyblue", bins = 250) + geom_histogram(data = N, alpha = 0.5, fill = "grey", bins = 250) + geom_histogram(data = Psi, alpha = 0.5, fill = "red", bins = 250) + theme_bw()
ggsave (filename = paste0 ("Histogramme_colour_score_max_6.pdf"), plot = plotmax2)

windows()
plotA1 <- ggplot (graphe, aes(scoreA)) + xlim(0,1) + geom_histogram(alpha = 0.3, fill = "green", bins = 250) + theme_bw()
ggsave (filename = paste0 ("Histogramme_global_score_A_6.pdf"), plot = plotA1)

windows()
plotA2 <- ggplot (Nm, aes(scoreA)) + xlim(0,1) + geom_histogram(alpha = 0.75, fill = "deepskyblue", bins = 250) + geom_histogram(data = N, alpha = 0.5, fill = "grey", bins = 250) + geom_histogram(data = Psi, alpha = 0.5, fill = "red", bins = 250) + theme_bw()
ggsave (filename = paste0 ("Histogramme_colour_score_A_6.pdf"), plot = plotA2)


graphics.off()


library(ggplot2)
library(gridExtra)

ggplot (graphe, aes(max)) + 
xlim(0,1) + 
geom_density(data = Nm, alpha = 0.75, fill = "deepskyblue", colour="deepskyblue") + 
geom_density(data = N, alpha = 0.5, fill = "grey", colour="grey") + 
geom_density(data = Psi, alpha = 0.5, fill = "red", colour="red") +
theme_bw()
ggsave (filename = paste0 ("Density_score_max_6.pdf"), plot = last_plot())


ggplot (Nm, aes(scoreA)) + 
xlim(0,1) + 
geom_density(alpha = 0.75, fill = "deepskyblue", colour = "deepskyblue") + 
geom_density(data = N, alpha = 0.5, fill = "grey", colour = "grey") + 
geom_density(data = Psi, alpha = 0.5, fill = "red", colour = "red") +
theme_bw()
ggsave (filename = paste0 ("Density_score_A_6.pdf"), plot = last_plot())




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


p1 <- ggplot (Nm, aes(x = max_rank, y = scoreA_rank)) + 
geom_point(alpha = 0.75, colour = "deepskyblue") + 
geom_point(data = N, alpha = 0.5, colour = "grey", size = 0.2) + 
geom_point(data = Psi, alpha = 0.5, colour = "red") +
scale_x_continuous(expand=c(0.02,0)) +
scale_y_continuous(expand=c(0.02,0)) +
theme_bw()

p2 <- ggplot (Nm, aes(x = max_rank)) + 
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

ggsave (filename = paste0 ("Density+plot_rank_6.pdf"), plot = grid.arrange(arrangeGrob(p2,ncol=2,widths=c(3,1)),
             arrangeGrob(p1,p3,ncol=2,widths=c(3,1)),
             heights=c(1,3)), width = 7, height = 7)



p1 <- ggplot (Nm, aes(x = scoreA, y = max)) + 
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

p3 <- ggplot (graphe, aes(max)) + 
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

ggsave (filename = paste0 ("Density+plot_6_score.pdf"), plot = grid.arrange(arrangeGrob(p2,ncol=2,widths=c(3,1)),
             arrangeGrob(p1,p3,ncol=2,widths=c(3,1)),
             heights=c(1,3)), width = 7, height = 7)

graphics.off()

