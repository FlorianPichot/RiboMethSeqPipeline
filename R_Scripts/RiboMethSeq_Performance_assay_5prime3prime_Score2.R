rm(list=ls(all=TRUE))

####### Input Area #######

type <- list("5.8S", "18S", "28S")		      # rRNA list to study
extDataDir <-"Path/to/Samples/Directory"		# Directory with 'counts' files
ResultsDir <-"Path/to/Results/Directory"		# Directory where to put the results (A directory is created if needed)

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

ResultsDir <- paste0(ResultsDir,"/2DensityPlot(ScoreMEAN2)_MCC_", basename(extDataDir))
dir.create (ResultsDir, showWarnings=FALSE)

ref_Hs <- read.csv2("Hsapiens_rRNA.csv", sep=" ", blank.lines.skip= TRUE, header = FALSE, col.names=(c("type","position","identity")))
ref_Hs$identity <- paste0 (ref_Hs$identity, "_Hs") 
ref_Hs <- na.omit(ref_Hs)


#READING local DATA into separate dataframes

memory <- data.frame ()
total <- data.frame ()

for (z in DataDirs) {
  
  tot_score_mean <- data.frame()
  tot_scoreA <- data.frame()
  
  for (y in type) {
    SampleDir <- paste0(ResultsDir,"/", z)
    dir.create (SampleDir, showWarnings = FALSE)
    
    
    results <-NA
    miR_local_dir <- paste(extDataDir, "/",z, sep = "")
    #print (miR_local_dir)
    #print (z)
    miR_file <- paste (miR_local_dir,"/", "Count_hs_rRNA_",y,".csv", sep ="")
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
    
    results$counts [is.na(results$counts)] <- 1	  
    
    results$ratio <- NA	
    
    for (i in 2:n) {results$ratio[i] = (results$counts[i]/results$counts[i-1])}
    
    results$around <- NA
    for (i in 3:n) {
      j <- i-2
      k <- i-1 
      l <- i+1 
      m <- i+2
      results$around[i] = (1 - results$ratio[i]/(mean(results$ratio[j:k])+ mean (results$ratio[l:m])))}
    
    
    results$deviation <- NA
    for (i in 4:n) {
      j <- i-3 
      k <- i-1 
      l <- i+1 
      m <- i+3
      results$deviation[i] = (1 - results$around[i]/(sd(results$around[j:k])+ sd(results$around[l:m])))}
    
    results$ratio2 <- NA	
    for (i in 2:n) {results$ratio2[i] = 1/results$ratio[i]}
    
    
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
    for (i in 3:n) {results$mean[i] = results$cumul[i]}
    
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
    
    
    ref_Hs_type <- ref_Hs[grep (paste0('^',y), ref_Hs$type), ]
    score_mean <- merge (results_mean, ref_Hs_type, by="position", all = TRUE)
    score_mean$identity [is.na(score_mean$identity)] <- "No_Met"
    score_mean <- score_mean[order(score_mean$mean, decreasing = TRUE),]
    scoreA <- merge (results_scoreA, ref_Hs_type, by="position", all = TRUE)
    scoreA$identity [is.na(scoreA$identity)] <- "No_Met"
    scoreA <- scoreA[order(scoreA$scoreA, decreasing = TRUE),]
    
    write.table (score_mean, file = paste0(SampleDir, "/Score_mean_", y, ".csv"), sep = " ", row.names=FALSE)
    write.table (scoreA, file = paste0(SampleDir, "/Score_A_", y, ".csv"), sep = " ", row.names=FALSE)
    
    fusion <- merge(scoreA, score_mean, by = c("position","type","identity"))
    fusion$type <- y
    
    memory <- rbind (memory,fusion)
    tot_score_mean <- rbind(tot_score_mean, fusion)
    tot_scoreA <- rbind(tot_scoreA, fusion)
    
  }
  ###### ROC Curves calculation (with functions) ######
  
  ###### ROC curve ######
  # mean #
  tot_score_mean <- tot_score_mean[order(tot_score_mean$mean, decreasing = TRUE),]
  
  tot_score_mean$TP <- fctTP(identity = tot_score_mean$identity)
  tot_score_mean$FP <- fctFP(identity = tot_score_mean$identity)
  tot_score_mean$FN <- fctFN(identity = tot_score_mean$identity, TP = tot_score_mean$TP)
  tot_score_mean$TN <- fctTN(identity = tot_score_mean$identity, FP = tot_score_mean$FP)
  tot_score_mean$MCC <- fctMCC(TP = tot_score_mean$TP, TN = tot_score_mean$TN, FP = tot_score_mean$FP, FN = tot_score_mean$FN)
  tot_score_mean$TPR <- fctTPR(TP = tot_score_mean$TP, FN = tot_score_mean$FN)
  tot_score_mean$TNR <- fctTNR(TN = tot_score_mean$TN, FP = tot_score_mean$FP)
  tot_score_mean$PPV <- fctPPV(TP = tot_score_mean$TP, FP = tot_score_mean$FP)
  tot_score_mean$NPV <- fctNPV(TN = tot_score_mean$TN, FN = tot_score_mean$FN)
  tot_score_mean$FNR <- fctFNR(TP = tot_score_mean$TP, FN = tot_score_mean$FN)
  tot_score_mean$FPR <- fctFPR(TN = tot_score_mean$TN, FP = tot_score_mean$FP)
  tot_score_mean$FDR <- fctFDR(TP = tot_score_mean$TP, FP = tot_score_mean$FP)
  tot_score_mean$FOR <- fctFOR(TN = tot_score_mean$TN, FN = tot_score_mean$FN)
  tot_score_mean$ACC <- fctACC(TP = tot_score_mean$TP, TN = tot_score_mean$TN, FP = tot_score_mean$FP, FN = tot_score_mean$FN)
  
  write.table (tot_score_mean, file=paste0(SampleDir, "/ROC_mean.csv"), sep = " ", row.names = FALSE)
  
  # ROCurve score mean #
  dev.new()
  plot (tot_score_mean$FPR, tot_score_mean$TPR, type = "l", col = "blue", xlim = c(0,0.05), main = "ROCurve for Score mean", sub=paste0("MCC = ", max(unlist(tot_score_mean$MCC))), xlab=" ", ylab=" ")
  lines(tot_score_mean$FPR, tot_score_mean$MCC, type = "l", col = "orange")
  savePlot (filename = paste0 (SampleDir, "/ROCurve_Score_mean"), "pdf")
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
  plot (tot_scoreA$FPR, tot_scoreA$TPR, type = "l", col = "blue", xlim = c(0,0.05), main = "ROCurve for Score A", sub=paste0("MCC = ", max(unlist(tot_scoreA$MCC))), xlab=" ", ylab=" ")
  lines(tot_scoreA$FPR, tot_scoreA$MCC, type = "l", col = "orange")
  savePlot (filename = paste0 (SampleDir, "/ROCurve_Score_A"), "pdf")
  graphics.off()
  
  # Compilation of MCC results #
  
  MCC_scoreA <- tot_scoreA[grep(max(tot_scoreA$MCC, na.rm = TRUE),tot_scoreA$MCC),]
  MCC_score_mean <- tot_score_mean[grep(max(tot_score_mean$MCC, na.rm = TRUE),tot_score_mean$MCC),]
  
  MCC_total <-cbind(z, MCC_score_mean$MCC, MCC_score_mean$TPR, MCC_score_mean$PPV, MCC_score_mean$FDR, MCC_scoreA$MCC, MCC_scoreA$TPR, MCC_scoreA$PPV, MCC_scoreA$FDR)
  colnames (MCC_total) <- c("sample", "MCC_mean", "TPR_mean", "PPV_mean", "FDR_mean", "MCC_scoreA", "TPR_scoreA", "PPV_scoreA", "FDR_scoreA")
  total <- rbind (total, MCC_total)
  
}

write.table (total, file=paste0(ResultsDir, "/ROC_total_MEAN2_5prime3prime.csv"), sep = " ", row.names = FALSE)

n = length (memory$position)

memory <- memory[order(memory$mean, decreasing = TRUE),]
memory$mean_rank <- 1:n
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
plotmean1 <- ggplot (graphe, aes(mean)) + xlim(0,1) + geom_histogram(alpha = 0.3, fill = "green", bins = 250) + theme_bw()
ggsave (filename = paste0 ("Histogramme_global_score_mean_2.pdf"), plot = plotmean1)

windows()
plotmean2 <- ggplot (Nm, aes(mean)) + xlim(0,1) + geom_histogram(alpha = 0.75, fill = "deepskyblue", bins = 250) + geom_histogram(data = N, alpha = 0.5, fill = "grey", bins = 250) + geom_histogram(data = Psi, alpha = 0.5, fill = "red", bins = 250) + theme_bw()
ggsave (filename = paste0 ("Histogramme_colour_score_mean_2.pdf"), plot = plotmean2)

windows()
plotA1 <- ggplot (graphe, aes(scoreA)) + xlim(0,1) + geom_histogram(alpha = 0.3, fill = "green", bins = 250) + theme_bw()
ggsave (filename = paste0 ("Histogramme_global_score_A_2.pdf"), plot = plotA1)

windows()
plotA2 <- ggplot (Nm, aes(scoreA)) + xlim(0,1) + geom_histogram(alpha = 0.75, fill = "deepskyblue", bins = 250) + geom_histogram(data = N, alpha = 0.5, fill = "grey", bins = 250) + geom_histogram(data = Psi, alpha = 0.5, fill = "red", bins = 250) + theme_bw()
ggsave (filename = paste0 ("Histogramme_colour_score_A_2.pdf"), plot = plotA2)


graphics.off()


library(ggplot2)
library(gridExtra)

ggplot (graphe, aes(mean)) + 
  xlim(0,1) + 
  geom_density(data = Nm, alpha = 0.75, fill = "deepskyblue", colour="deepskyblue") + 
  geom_density(data = N, alpha = 0.5, fill = "grey", colour="grey") + 
  geom_density(data = Psi, alpha = 0.5, fill = "red", colour="red") +
  theme_bw()
ggsave (filename = paste0 ("Density_score_mean_2.pdf"), plot = last_plot())


ggplot (Nm, aes(scoreA)) + 
  xlim(0,1) + 
  geom_density(alpha = 0.75, fill = "deepskyblue", colour = "deepskyblue") + 
  geom_density(data = N, alpha = 0.5, fill = "grey", colour = "grey") + 
  geom_density(data = Psi, alpha = 0.5, fill = "red", colour = "red") +
  theme_bw()
ggsave (filename = paste0 ("Density_score_A_2.pdf"), plot = last_plot())




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

ggsave (filename = paste0 ("Density+plot_rank_2.pdf"), plot = grid.arrange(arrangeGrob(p2,ncol=2,widths=c(3,1)),
                                                                           arrangeGrob(p1,p3,ncol=2,widths=c(3,1)),
                                                                           heights=c(1,3)), width = 7, height = 7)



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

ggsave (filename = paste0 ("Density+plot_2_score.pdf"), plot = grid.arrange(arrangeGrob(p2,ncol=2,widths=c(3,1)),
                                                                            arrangeGrob(p1,p3,ncol=2,widths=c(3,1)),
                                                                            heights=c(1,3)), width = 7, height = 7)

graphics.off()
