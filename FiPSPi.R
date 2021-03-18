##########################################################################
#                          COPYRIGHT & SUPPORT                           #
##########################################################################
#
# FiPSPi: Random Forest-based muscle fiber peptide selection pipeline.
#
# This code has been developed by Dr. Michael Turewicz at the Ruhr University
# Bochum, Germany. It is licensed under the BSD 3-Clause License provided in
# the file 'LICENSE.txt'.
#
# For support please write an e-mail to 'michael.turewicz[at]rub.de'.
#
# This is the original version of the code written for the publication
# by Eggers et al. 'Deep proteomic characterization of skeletal muscle fiber
# types by laser microdissection and mass spectrometry (...)'. The latest
# version of FiPSPi can be found here:
# https://github.com/mpc-bioinformatics/FiPSPi
#
##########################################################################
#                            USER SETTINGS                               #
##########################################################################

### set path to current working directory
cwd <- "C:/UNI/Publikationen/2021.XX.XX_writing_Britta_Eggers_Fasertypen/FiPSPi"

### set path to data file (relative path to the working directory specified above)
data.path <- paste0(cwd, "/data/Eggers_et_al_data.txt")

### set path to metadata file (relative path to the working directory specified above)
metadata.path <- paste0(cwd, "/data/Eggers_et_al_metadata.txt")

### set output.path where all created output and graphics are saved
output.path <- paste0(cwd, "/FiPSPi_results")

### -------> Note: the above files and folders must already exist and be accessible! 

##########################################################################
#                 FURTHER SETTINGS & DEPENDENCIES                        #
##########################################################################

options(scipen=100000)
options(stringsAsFactors=FALSE)

set.seed(1234)

library(randomForest)
library(ROCR)
library(caret)
library(rpart)
library(rpart.plot)

##########################################################################
#                     FUNCTION DEFINITIONS                               #
##########################################################################

#++++++++++++++++++++++ blueyellow ++++++++++++++++++++++++++++++++++++
blueyellow<-function(n){
  colorpanel(n, "blue", "yellow")
}
#++++++++++++++++++++++ blueyellow ++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++ plotFeaturesHeatmap.2 ++++++++++++++++++++++++++++++++++

plotFeaturesHeatmap.2 <- function(features=NULL, x=NULL, n1=NULL, n2=NULL, n3=NULL, n4=NULL,
                                  output.path=NULL, file.name=FALSE, description=FALSE){
  
  if(is.null(features) || is.null(x) || is.null(n1) || is.null(n2)) {
    stop("ERROR: Not all mandatory arguments have been defined!")
  }
  
  features <- c(na.exclude(features))
  datamatrix <- x
  brc <- rownames(x)
  rows <- brc %in% features
  datamatrix <- datamatrix[rows,]
  if(description){
    rownames(datamatrix) <- elist$genes$Description[rows]
  }else{
    rownames(datamatrix) <- brc[rows]
  }
  
  #title <- ""
  if(nrow(datamatrix) > 100) {
    cexRowSize <- 0.1
  } else if(nrow(datamatrix) > 40) {
    cexRowSize <- 0.4
  } else if(nrow(datamatrix) > 10) {
    cexRowSize <- 0.8
  }else {
    cexRowSize <- 1.5
  }
  
  cexColSize <- 0.8
  
  #Col labels with default distance from heatmap (NA),
  #but aligned to middle of heatmap-cells (0.5); default is c(NA,0)
  adjCol <- c(NA,0.5) 
  
  #defining ColSideColors
  ColSideColors <- c(rep("gray25",n1), rep("gray75",n2), rep("red",n3), rep("blue",n4))
  
  #defining the correlation-based distance function
  my.dist <- function(x) as.dist((1-cor(t(x)))/2)
  
  
  #Defining the dimensions and layout of the heatmap and its color key
  lmat <- rbind(c(5,4),c(0,1),c(3,2))
  lwid <- c(1.5,4)
  lhei <- c(1.5,0.2,4)
  
  if(!is.null(output.path)){
    png(
      paste0(output.path,"/", file.name),
      width = 2000,
      height = 2000,
      pointsize = 10,
      res = 300
    )
    heatmap.2(datamatrix,
              distfun=function(x) as.dist((1-cor(t(x)))/2),
              Rowv = TRUE,                  #row and col clustering
              Colv = TRUE,
              ColSideColors=ColSideColors,
              col=blueyellow(300),          #color scheme
              scale="row",                  #scale by column
              margins=c(15, 15),            #margins around heatmap
              key=TRUE,                     #legend is present
              key.title="Color key",        #title for legend
              key.xlab="Expression values", #label for legend axis
              key.ylab="Density",           #label for legend axis
              keysize=1.1,                  #size of legend strip
              key.par=list(cex.lab=1.4),
              symkey=FALSE,                 #colors not symmetrical around 0
              density.info="density",       #no density information
              trace="none",                 
              cexRow=cexRowSize,            #row labels' size
              cexCol=cexColSize,            #column labels' size
              labRow=substr(row.names(datamatrix),1,35), #row labels
              labCol=colnames(datamatrix),   #column labels: limit to 30 chars
              adjCol=adjCol,
              lmat=lmat,
              lhei=lhei,
              lwid=lwid
    )
    #normalized device coordinates (NDC): c(xmin, xmax, ymin, ymax)
    par(fig=c(0,0.975,0,0.94), new=TRUE)
    legend("topright",      # legend location
           legend = c(unique(gsub("\\_\\d+", "", colnames(datamatrix)))),
           col = c("gray25", "gray75", "red", "blue"),  # color key
           pch = c(15,15),
           cex=1.75,
           bty="n",
           horiz=TRUE
    )
    dev.off()
  } else {
    heatmap.2(datamatrix,
              distfun=function(x) as.dist((1-cor(t(x)))/2),
              Rowv = TRUE,                  #row and col clustering
              Colv = TRUE,
              ColSideColors=ColSideColors,
              col=blueyellow(300),          #color scheme
              scale="row",                  #scale by column
              margins = c(15,15),           #margins around heatmap
              key=TRUE,                     #legend is present
              key.title="Color key",        #title for legend
              key.xlab="Expression values", #label for legend axis
              key.ylab="Density",           #label for legend axis
              keysize=1.1,                  #size of legend strip
              key.par=list(cex.lab=1.4),    #size of legend strip
              symkey=FALSE,                 #colors not symmetrical around 0
              density.info="density",       #no density information
              trace="none",                 
              cexRow=cexRowSize,            #column labels' size
              cexCol=cexColSize,            #column labels' size
              labRow=substr(row.names(datamatrix),1,35), #row labels
              labCol=colnames(datamatrix),   #column labels: limit to 30 chars
              adjCol=adjCol,
              lmat=lmat,
              lhei=lhei,
              lwid=lwid
    )
    #normalized device coordinates (NDC): c(xmin, xmax, ymin, ymax)
    par(fig=c(0,0.975,0,0.94), new=TRUE)
    legend("topright",      # legend location
           legend = c(unique(gsub("\\_\\d+", "", colnames(datamatrix)))),
           col = c("gray25", "gray75", "red", "blue"),  # color key
           pch = c(15,15),
           cex=1.75,
           bty="n",
           horiz=TRUE
    )
  }
  
}

#++++++++++++++++++++++ plotFeaturesHeatmap.2 ++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ rf.rfe +++++++++++++++++++++++++++++++++++++++
rf.rfe <- function(iteration=NULL, datamatrix=NULL, output.path=NULL,
                   label1="A", label2="B", label3="C", label4="D", n1=NULL, n2=NULL, n3=NULL,
                   n4=NULL, panel.selection.criterion="accuracy", importance.measure="MDA",
                   ntree=500, mtry=NULL, verbose=FALSE){
  if(is.null(iteration) || is.null(datamatrix) || is.null(n1) ||
     is.null(n2) || is.null(n3) || is.null(n4)) {
    stop("ERROR: Not all mandatory arguments have been defined!")
  }
  
  n.min <- min(n1,n2,n3,n4)
  print(n.min)
  classes <- factor(c(rep(label1, n1), rep(label2, n2), rep(label3, n3), rep(label4, n4)))
  print(classes)
  print(length(classes))
  acc.best.set.accuracy <- 0
  acc.best.set <- c()
  sens.best.set.sensitivity <- 0
  sens.best.set <- c()
  spec.best.set.specificity <- 0
  spec.best.set <- c()
  rf.importance.total <- matrix(0, nrow=nrow(datamatrix), ncol=2)
  colnames(rf.importance.total) <- c("MDA", "MDG")
  rownames(rf.importance.total) <- rownames(datamatrix)
  if(!is.null(output.path)){
    cat(x="set\ttotalAccuracy\ttotalSensi\ttotalSpeci\tAccuracy1\tSensitivity1\tSpecificity1\tAccuracy2\tSensitivity2\tSpecificity2\tAccuracy3\tSensitivity3\tSpecificity3\tAccuracy4\tSensitivity4\tSpecificity4\tCM11\tCM12\tCM13\tCM14\tCM21\tCM22\tCM23\tCM24\tCM31\tCM32\tCM33\tCM34\tCM41\tCM42\tCN43\tCM44\n",
        file=paste(output.path, "/rfe_", iteration, ".txt", sep=""),
        append=FALSE)
  }
  
  if(is.null(ntree)){
    ntree <- 500
  }
  
  iterations <- 0
  while({p <- nrow(datamatrix)} > 0){
    iterations <- iterations + 1
    if(is.null(mtry)){
      mtry.tmp <- sqrt(p)
    }else{
      mtry.tmp <- mtry
    }
    
    #~~~~~~~~~~~~~~~CLASSIFIER END~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dat <- t(datamatrix)
    model_randomForest <- randomForest(x=dat, y=classes, importance=TRUE, keep.forest=FALSE,
                                       ntree=ntree, mtry=mtry.tmp, sampsize=c(n.min, n.min, n.min, n.min))
    
    importance <- importance(model_randomForest)
    rf.importance <- data.frame(rownames(importance), importance[,3], importance[,4])
    names(rf.importance) <- c('Var', 'MDA', 'MDG')
    
    confusion <- model_randomForest$confusion
    print(confusion)
    #~~~~~~~~~~~~~~~CLASSIFIER END~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ### class 1:
    TP <- confusion[1,1]
    FN <- confusion[1,2]+confusion[1,3]+confusion[1,4]
    FP <- confusion[2,1]+confusion[3,1]+confusion[4,1]
    TN <- confusion[2,2]+confusion[2,3]+confusion[2,4]+confusion[3,2]+confusion[3,3]+confusion[3,4]+confusion[4,2]+confusion[4,3]+confusion[4,4]
    ACCURACY1 <- {TP+TN}/{TP+FP+TN+FN}
    SENSITIVITY1 <- TP/{TP+FN}   #=TPR
    SPECIFICITY1 <- TN/{TN+FP}  #=(1-FPR) 
    
    ### class 2:
    TP <- confusion[2,2]
    FN <- confusion[2,1]+confusion[2,3]+confusion[2,4]
    FP <- confusion[1,2]+confusion[3,2]+confusion[4,2]
    TN <- confusion[1,1]+confusion[1,3]+confusion[1,4]+confusion[3,1]+confusion[3,3]+confusion[3,4]+confusion[4,1]+confusion[4,3]+confusion[4,4]
    ACCURACY2 <- {TP+TN}/{TP+FP+TN+FN}
    SENSITIVITY2 <- TP/{TP+FN}   #=TPR
    SPECIFICITY2 <- TN/{TN+FP}  #=(1-FPR) 
    
    ### class 3:
    TP <- confusion[3,3]
    FN <- confusion[3,1]+confusion[3,2]+confusion[3,4]
    FP <- confusion[1,3]+confusion[2,3]+confusion[4,3]
    TN <- confusion[1,1]+confusion[1,2]+confusion[1,4]+confusion[2,1]+confusion[2,2]+confusion[2,4]+confusion[4,1]+confusion[4,2]+confusion[4,4]
    ACCURACY3 <- {TP+TN}/{TP+FP+TN+FN}
    SENSITIVITY3 <- TP/{TP+FN}   #=TPR
    SPECIFICITY3 <- TN/{TN+FP}  #=(1-FPR)
    
    ### class 4:
    TP <- confusion[4,4]
    FN <- confusion[4,1]+confusion[4,2]+confusion[4,3]
    FP <- confusion[1,4]+confusion[2,4]+confusion[3,4]
    TN <- confusion[1,1]+confusion[1,2]+confusion[1,3]+confusion[2,1]+confusion[2,2]+confusion[2,3]+confusion[3,1]+confusion[3,2]+confusion[3,3]
    ACCURACY4 <- {TP+TN}/{TP+FP+TN+FN}
    SENSITIVITY4 <- TP/{TP+FN}   #=TPR
    SPECIFICITY4 <- TN/{TN+FP}  #=(1-FPR)
    
    ACCURACY_total <- sum(diag(confusion))/ sum(confusion) # total Accuracy
    #ACCURACY_total <- (ACCURACY1 + ACCURACY2 + ACCURACY3 + ACCURACY4) / 4 # total Accuracy
    SENSITIVITY_total <- (SENSITIVITY1 + SENSITIVITY2 + SENSITIVITY3 + SENSITIVITY4) / 4
    SPECIFICITY_total <- (SPECIFICITY1 + SPECIFICITY2 + SPECIFICITY3 + SPECIFICITY4) / 4
    CM <- as.vector(t(confusion[,-ncol(confusion)]))  # Einträge der Confusion Matrix als Vector (zum Abspeichern)
    
    
    
    if(!is.null(output.path)){
      cat(x=paste(nrow(datamatrix), ACCURACY_total, SENSITIVITY_total, SPECIFICITY_total,
                  ACCURACY1,  SENSITIVITY1, SPECIFICITY1,
                  ACCURACY2,  SENSITIVITY2, SPECIFICITY2,
                  ACCURACY3,  SENSITIVITY3, SPECIFICITY3,
                  ACCURACY4,  SENSITIVITY4, SPECIFICITY4,
                  CM[1], CM[2],CM[3],CM[4],CM[5],CM[6],CM[7],CM[8],CM[9],CM[10],CM[11],CM[12],CM[13],CM[14],CM[15],CM[16],  sep="\t"), 
          file=paste(output.path, "/rfe_", iteration, ".txt", sep=""), append=TRUE)
      cat(x="\n", file=paste(output.path, "/rfe_", iteration, ".txt",
                             sep=""), append=TRUE)
    }
    
    
    
    if(verbose){
      message(paste0("rf.rfe - features: ", p, ", ACCURACY: ", ACCURACY_total), "\n")
    }
    
    
    
    if(acc.best.set.accuracy <= ACCURACY_total){
      acc.best.set.accuracy <- ACCURACY_total
      acc.best.set <- rownames(datamatrix)
    }
    #if(sens.best.set.sensitivity <= SENSITIVITY){
    #  sens.best.set.sensitivity <- SENSITIVITY
    #  sens.best.set <- rownames(datamatrix)
    #}
    if(spec.best.set.specificity <= SPECIFICITY_total){
      spec.best.set.specificity <- SPECIFICITY_total
      spec.best.set <- rownames(datamatrix)
    }
    
    for(j in 1:nrow(rf.importance)){
      rf.importance.total[rf.importance[j,"Var"],"MDA"] <- rf.importance.total[rf.importance[j,"Var"],"MDA"] + rf.importance[j,"MDA"]
      rf.importance.total[rf.importance[j,"Var"],"MDG"] <- rf.importance.total[rf.importance[j,"Var"],"MDG"] + rf.importance[j,"MDG"]
    }
    
    if(importance.measure == "MDA"){
      rf.importance <- rf.importance[order(-rf.importance$MDA),]
    }else if(importance.measure == "MDG"){
      rf.importance <- rf.importance[order(-rf.importance$MDG),]
    }
    
    if(nrow(rf.importance) != 1){
      next.set <- rf.importance$Var[1:(nrow(importance)-1)]
    }else{
      break
    }
    next.set <-
      rownames(datamatrix)[rownames(datamatrix) %in% next.set]
    datamatrix <- datamatrix[next.set,,drop=FALSE]
  }
  
  if(!is.null(output.path)){
    importance.results <- cbind(rownames(rf.importance.total), rf.importance.total)
    write.table(x=importance.results, file=paste0(output.path, "/rfe_importance_", iteration, ".txt"), append=TRUE, sep="\t", eol="\n", dec=".", row.names=FALSE, col.names=TRUE)
  }
  
  if(panel.selection.criterion == "accuracy"){
    if(verbose){
      message(paste("feature selection - optimal number of features: ",
                    length(acc.best.set), sep=""), "\n")
      message(paste("feature selection - best accuracy: ",
                    acc.best.set.accuracy, sep=""), "\n")
    }
    write.table(x=importance.results[acc.best.set,], file=paste0(output.path, "/rfe_importance_panel.txt"), append=TRUE, sep="\t", eol="\n", dec=".", row.names=FALSE, col.names=TRUE)
    return(acc.best.set)
  }else if(panel.selection.criterion == "specificity"){
    if(verbose){
      message(paste("feature selection - optimal number of features: ",
                    length(spec.best.set), sep=""), "\n")
      message(paste("feature selection - best specificity: ",
                    spec.best.set.specificity, sep=""), "\n")
    }
    write.table(x=importance.results[spec.best.set,], file=paste0(output.path, "/rfe_importance_panel.txt"), append=TRUE, sep="\t", eol="\n", dec=".", row.names=FALSE, col.names=TRUE)   
    return(spec.best.set)
  }
  
}
#+++++++++++++++++++++++++++ rf.rfe +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ classify.rf.evaluation +++++++++++++++++++++++++++++++++++++++
classify.rf.evaluation <- function(iteration=NULL, trainset=NULL, testset=NULL, classes_train=NULL, 
                                   classes_test=NULL, label1="A", label2="B", label3 = "C", label4 = "D", output.path=NULL, ...){
  if(is.null(trainset) || is.null(testset)) {
    stop("ERROR: Not all mandatory arguments have been defined!")
  }
  
  train.dat <- t(as.matrix(trainset))
  test.dat <- t(as.matrix(testset))
  
  
  
  model.rf <- randomForest(x=train.dat, y=classes_train, importance=TRUE, keep.forest=TRUE)
  pred.rf <- predict(object=model.rf, newdata=test.dat, type="response", norm.votes=TRUE)
  confusion <- table(observed = classes_test, predicted = pred.rf)
  
  
  
  ### class 1:
  TP <- confusion[1,1]
  FN <- confusion[1,2]+confusion[1,3]+confusion[1,4]
  FP <- confusion[2,1]+confusion[3,1]+confusion[4,1]
  TN <- confusion[2,2]+confusion[2,3]+confusion[2,4]+confusion[3,2]+confusion[3,3]+confusion[3,4]+confusion[4,2]+confusion[4,3]+confusion[4,4]
  ACCURACY1 <- {TP+TN}/{TP+FP+TN+FN}
  SENSITIVITY1 <- TP/{TP+FN}   #=TPR
  SPECIFICITY1 <- TN/{TN+FP}  #=(1-FPR) 
  
  ### class 2:
  TP <- confusion[2,2]
  FN <- confusion[2,1]+confusion[2,3]+confusion[2,4]
  FP <- confusion[1,2]+confusion[3,2]+confusion[4,2]
  TN <- confusion[1,1]+confusion[1,3]+confusion[1,4]+confusion[3,1]+confusion[3,3]+confusion[3,4]+confusion[4,1]+confusion[4,3]+confusion[4,4]
  ACCURACY2 <- {TP+TN}/{TP+FP+TN+FN}
  SENSITIVITY2 <- TP/{TP+FN}   #=TPR
  SPECIFICITY2 <- TN/{TN+FP}  #=(1-FPR) 
  
  ### class 3:
  TP <- confusion[3,3]
  FN <- confusion[3,1]+confusion[3,2]+confusion[3,4]
  FP <- confusion[1,3]+confusion[2,3]+confusion[4,3]
  TN <- confusion[1,1]+confusion[1,2]+confusion[1,4]+confusion[2,1]+confusion[2,2]+confusion[2,4]+confusion[4,1]+confusion[4,2]+confusion[4,4]
  ACCURACY3 <- {TP+TN}/{TP+FP+TN+FN}
  SENSITIVITY3 <- TP/{TP+FN}   #=TPR
  SPECIFICITY3 <- TN/{TN+FP}  #=(1-FPR)
  
  ### class 4:
  TP <- confusion[4,4]
  FN <- confusion[4,1]+confusion[4,2]+confusion[4,3]
  FP <- confusion[1,4]+confusion[2,4]+confusion[3,4]
  TN <- confusion[1,1]+confusion[1,2]+confusion[1,3]+confusion[2,1]+confusion[2,2]+confusion[2,3]+confusion[3,1]+confusion[3,2]+confusion[3,3]
  ACCURACY4 <- {TP+TN}/{TP+FP+TN+FN}
  SENSITIVITY4 <- TP/{TP+FN}   #=TPR
  SPECIFICITY4 <- TN/{TN+FP}  #=(1-FPR)
  
  
  
  ACCURACY_total <- sum(diag(confusion))/ sum(confusion) # total Accuracy
  SENSITIVITY_total <- (SENSITIVITY1 + SENSITIVITY2 + SENSITIVITY3 + SENSITIVITY4) / 4
  SPECIFICITY_total <- (SPECIFICITY1 + SPECIFICITY2 + SPECIFICITY3 + SPECIFICITY4) / 4
  
  
  
  results <- c(acc_total = ACCURACY_total, sens_total = SENSITIVITY_total, spec_total = SPECIFICITY_total, 
               acc1 = ACCURACY1,  sens1 = SENSITIVITY1, spec1 = SPECIFICITY1,
               acc2 = ACCURACY2,  sens2 = SENSITIVITY2, spec2 = SPECIFICITY2,
               acc3 = ACCURACY3,  sens3 = SENSITIVITY3, spec3 = SPECIFICITY3,
               acc4 = ACCURACY4,  sens4 = SENSITIVITY4, spec4 = SPECIFICITY4)
  
  return(results)
}
#+++++++++++++++++++++++++++ classify.rf.evaluation +++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++ combinePerf +++++++++++++++++++++++++++++++++++++++
combinePerf <- function(perf.list, label1, label2, filename, cwd){
  if(is.null(perf.list)) {
    stop("ERROR: Not all mandatory arguments have been defined!")
  }
  
  decis <- c()
  classes <- c()
  for(m in 1:length(perf.list)){
    decis <- c(decis, perf.list[[m]]$pred)
    classes <- c(classes, as.character(perf.list[[m]]$classes))
  }
  
  prediction.obj <- prediction(decis, classes)
  performance.obj <- performance(prediction.obj, measure="tpr", x.measure="fpr")
  auc <- performance(prediction.obj,"auc")@y.values[[1]]
  
  png(paste(cwd,"/roc_perflist_", filename, ".png",sep=""), width=2000, height=2000, pointsize=17, res=300)
  plot(performance.obj, col="red", lwd=5,main=paste("AUC: ", round(auc,4), sep=""))
  abline(a=0,b=1,lwd=2,lty=2,col="gray")
  dev.off()
  
}
#+++++++++++++++++++++++++++ combinePerf +++++++++++++++++++++++++++++++++++++++



##########################################################################
#                               MAIN                                     #
##########################################################################

#-----------------------------> Begin: Data Import & Pre-Processing

dat <- read.table(data.path, sep="\t", header=TRUE, na.strings = "NA")
sampledat <- read.table(metadata.path, sep="\t", header=TRUE, na.strings = "NA")
rownames(sampledat) <- sampledat$ID
colnames(dat)[1:8] <- c("IsProteinGroupSpecific", "IsProteotypic", "StrippedSequence", "ProteinAccessions", "ProteinDescriptions", "ProteinGroups", "ProteinNames", "Qvalue")
colnames(dat)[9:64] <- rownames(sampledat)


label1 <- "T1" 
label2 <- "T2a"
label3 <- "T2B"
label4 <- "T2x"
group1 <- grep("T1", rownames(sampledat), value=TRUE)
group2 <- grep("T2a", rownames(sampledat), value=TRUE)
group3 <- grep("T2B", rownames(sampledat), value=TRUE)
group4 <- grep("T2x", rownames(sampledat), value=TRUE)
n1 <- length(group1)
n2 <- length(group2)
n3 <- length(group3)
n4 <- length(group4) 



descriptive.columns <- dat[,1:8]
descriptive.columns2 <- dat[,"StrippedSequence",drop=FALSE]
rownames(descriptive.columns2) <- dat[,"StrippedSequence"]
descriptive.columns[,"ProteinAccessions"] <- gsub("sp\\|", "", descriptive.columns[,"ProteinAccessions"])
descriptive.columns[,"ProteinAccessions"] <- gsub("\\|", "\\_", descriptive.columns[,"ProteinAccessions"])
descriptive.columns[,"ProteinAccessions"] <- gsub("\\;", "\\_", descriptive.columns[,"ProteinAccessions"])
accession.list <- unique(descriptive.columns[,"ProteinAccessions"])
row.IDs <- vector(mode="character", length=length(accession.list))
for(i in 1:length(accession.list)){
  idx <- grep(accession.list[i], descriptive.columns[,"ProteinAccessions"], fixed=TRUE)
  n.peps <- length(idx)
  row.IDs[idx] <- paste0(accession.list[i], "_peptide", 1:n.peps)
}
rownames(dat) <- row.IDs
rownames(descriptive.columns) <- row.IDs
descriptive.columns <- cbind(row.IDs, descriptive.columns)
descriptive.columns2 <- cbind(row.IDs, descriptive.columns2)
write.table(x=descriptive.columns, file=paste0(output.path, "/descriptive.columns.txt"), sep="\t", col.names=TRUE, row.names=FALSE)



dat <- data.matrix(dat[,9:64])
print(paste0("original feature number: ", nrow(dat)))
dat <- na.omit(dat)
print(paste0("feature number after removing features containing NAs: ", nrow(dat)))
write.table(x=cbind(rownames(dat),dat), file=paste0(output.path, "/original_data.txt"), sep="\t", col.names=TRUE, row.names=FALSE)

dat <- dat[,c(group1,group2,group3,group4)]

#-----------------------------> End: Data Import & Pre-Processing

#-----------------------------> Begin: Diff. Analysis

classes <- gsub("\\_\\d+", "", colnames(dat))
rawdat <- 2^dat
p.values <- vector(mode="numeric", length=nrow(dat))
max.mean.ratios <- vector(mode="numeric", length=nrow(dat))
for(i in 1:nrow(dat)){
  p.values[i] <- summary(aov(dat[i,] ~ classes))[[1]][["Pr(>F)"]][1]
  
  mean.ratio1 <- max( mean(rawdat[i,group1])/mean(c(rawdat[i,group2], rawdat[i,group3], rawdat[i,group4])), mean(c(rawdat[i,group2], rawdat[i,group3], rawdat[i,group4]))/mean(rawdat[i,group1]) )
  mean.ratio2 <- max( mean(rawdat[i,group2])/mean(c(rawdat[i,group1], rawdat[i,group3], rawdat[i,group4])), mean(c(rawdat[i,group1], rawdat[i,group3], rawdat[i,group4]))/mean(rawdat[i,group2]) )
  mean.ratio3 <- max( mean(rawdat[i,group3])/mean(c(rawdat[i,group1], rawdat[i,group2], rawdat[i,group4])), mean(c(rawdat[i,group1], rawdat[i,group2], rawdat[i,group4]))/mean(rawdat[i,group3]) )
  mean.ratio4 <- max( mean(rawdat[i,group4])/mean(c(rawdat[i,group1], rawdat[i,group2], rawdat[i,group3])), mean(c(rawdat[i,group1], rawdat[i,group2], rawdat[i,group3]))/mean(rawdat[i,group4]) )
  max.mean.ratios[i] <- max(mean.ratio1, mean.ratio2, mean.ratio3, mean.ratio4)
}
p.values.adj <- p.adjust(p.values, method="fdr")
names(p.values.adj) <- rownames(dat)
diff.analysis.output <- cbind(names(p.values.adj), p.values, p.values.adj, max.mean.ratios)
colnames(diff.analysis.output) <- c("Feature", "P-value (Anova)", "Adj. p-value", "Max. mean ratio")
write.table(file=paste0(output.path, "/diff_analysis_output.txt"), x=diff.analysis.output, 
            row.names=FALSE, col.names=TRUE, sep="\t")

#-----------------------------> End: Diff. Analysis



#-----------------------------> Begin: Peptide Selection via Random Forest-based Recursive Feature Elimination

criterion <- "accuracy" #alternative criterion: specificity
imp.measure <- "MDA"
features <- rf.rfe(iteration=1, datamatrix=dat, output.path=output.path, label1="T1", 
                   label2="T2a", label3="T2B", label4="T2x",  n1=n1, n2=n2, n3=n3, n4=n4,
                   panel.selection.criterion=criterion, 
                   importance.measure=imp.measure, verbose=TRUE)
selection.results <- cbind(features, descriptive.columns[features,"StrippedSequence"])
colnames(selection.results) <- c("Feature", "Sequence")
write.table(file=paste0(output.path, "/selected_features.txt"), x=selection.results,
            sep="\t", row.names=FALSE, col.names=TRUE)

rfe.iterations <- read.table(paste0(output.path, "/rfe_1.txt"), sep="\t", header=TRUE, na.strings = "NA")

p <- nrow(dat)
png(paste(output.path,"/rf-rfe_featNumber_vs_accuracy.png",sep=""), width=4000, height=4000, pointsize=30, res=300)
plot(p-rfe.iterations$set, rfe.iterations[,"totalAccuracy"], type="l", xaxt="n", xlab="Feature number", 
     ylab="Classification accuracy", main="RF-RFE: Feature number vs. accuracy \n (total accuracy)", 
     pch=19, lwd=12, col=adjustcolor("red", alpha=0.3), cex=1)
points(p-rfe.iterations$set, rfe.iterations[,"totalAccuracy"], pch=19, col=adjustcolor("red", alpha=1.0), cex=0.5)
axis(side=1, label=p:1, at=1:p)
dev.off()

#-----------------------------> End: Peptide Selection via Random Forest-based Recursive Feature Elimination




#-----------------------------> Begin: Validation with Training Set

### 1000 times train/test split for model validation
testruns <- 1000
all.idx <- colnames(dat)
perf.final <- matrix(nrow =  testruns, ncol = 15)
colnames(perf.final) = c("acc_total", "sens_total", "spec_total",
                         "acc1",  "sens1", "spec1",
                         "acc2",  "sens2", "spec2",
                         "acc3",  "sens3", "spec3",
                         "acc4",  "sens4", "spec4")

for(i in 1:testruns){
  cat("final classification: ", i, "\r")
  test.idx1 <- sample(group1, 5)  # n1=15
  test.idx2 <- sample(group2, 3)  # n2=13
  test.idx3 <- sample(group3, 5)
  test.idx4 <- sample(group4, 3)
  test.idx <- c(test.idx1, test.idx2, test.idx3, test.idx4)
  
  train.idx <- setdiff(all.idx, test.idx)
  
  classi <- classify.rf.evaluation(iteration=i, trainset=dat[features,train.idx,drop=FALSE], 
                                   testset=dat[features,test.idx,drop=FALSE], 
                                   classes_train=as.factor(gsub("\\_\\d+", "", train.idx)), 
                                   classes_test=as.factor(gsub("\\_\\d+", "", test.idx)), 
                                   label1="T1", label2="T2a", label3="T2B", label4="T2x", output.path=NULL, 
                                   nrounds = best.nrounds, max_depth = best.max_depth)
  perf.final[i,] <- classi
}

write.table(perf.final, file = paste0(output.path, "/validation.txt"), row.names = FALSE, sep = "\t")
if(length(features) > 2){
  plotFeaturesHeatmap.2(features=features, x=rawdat, n1=n1, n2=n2, n3=n3, n4=n4, output.path=output.path, file.name="heatmap.png", description=FALSE)
}

#-----------------------------> End: Validation with Training Set


#-----------------------------> Begin: PCA using all peptides without NAs

pca.dat <- dat

group.vec <- c(
  "T1",
  "T2a",
  "T2B",
  "T2x"
)

col.vec <- c(
  adjustcolor("navy", alpha=0.3),
  adjustcolor("red", alpha=0.3),
  adjustcolor("darkorchid", alpha=0.3),
  adjustcolor("darkgreen", alpha=0.3)
)


pcdat <- prcomp(t(pca.dat), center=TRUE, scale=TRUE)
scores <- pcdat$x
for (i in 1:2){
  for (j in i:2){
    if (i<j){
      XLIM <- c(-max(abs(scores[,i])), max(abs(scores[,i])))
      XLIM <- XLIM+(XLIM*0.1)
      YLIM <- c(-max(abs(scores[,j])), max(abs(scores[,j])))
      YLIM <- YLIM+(YLIM*0.1)
      png(paste(output.path, "/01_pca_", i, "_", j, "_filteredFeatures.png", sep=""), width=3600, height=3600, pointsize=15, res=600)
      plot(scores[group1,i], scores[group1,j], xlab=paste("PC", i, sep=""), ylab=paste("PC", j, sep=""), xlim=XLIM, ylim=YLIM, pch=20, col=col.vec[1], main="PCA", cex=2)
      points(scores[group2,i], scores[group2,j], pch=20, col=col.vec[2], cex=2)
      points(scores[group3,i], scores[group3,j], pch=20, col=col.vec[3], cex=2)
      points(scores[group4,i], scores[group4,j], pch=20, col=col.vec[4], cex=2)
      legend("topleft", legend=group.vec[1:4], col=col.vec[1:4], pch=20, cex=0.75, bg="transparent")
      dev.off()
    }
  }
}

for (i in 1:2){
  for (j in i:2){
    if (i<j){
      XLIM <- c(-max(abs(scores[,i])), max(abs(scores[,i])))
      XLIM <- XLIM+(XLIM*0.1)
      YLIM <- c(-max(abs(scores[,j])), max(abs(scores[,j])))
      YLIM <- YLIM+(YLIM*0.1)
      png(paste(output.path, "/02_pca_", i, "_", j, "_filteredFeatures.png", sep=""), width=3600, height=3600, pointsize=15, res=600)
      plot(scores[group1,i], scores[group1,j], xlab=paste("PC", i, sep=""), ylab=paste("PC", j, sep=""), xlim=XLIM, ylim=YLIM, pch=20, col=col.vec[1], main="PCA", cex=2)
      points(scores[group2,i], scores[group2,j], pch=20, col=col.vec[2], cex=2)
      points(scores[group3,i], scores[group3,j], pch=20, col=col.vec[3], cex=2)
      points(scores[group4,i], scores[group4,j], pch=20, col=col.vec[4], cex=2)
      text(scores[group1,i], scores[group1,j], labels=group1, col="navy", cex=0.4)
      text(scores[group2,i], scores[group2,j], labels=group2, col="red", cex=0.4)
      text(scores[group3,i], scores[group3,j], labels=group3, col="darkorchid", cex=0.4)
      text(scores[group4,i], scores[group4,j], labels=group4, col="darkgreen", cex=0.4)
      legend("topleft", legend=group.vec[1:4], col=col.vec[1:4], pch=20, cex=0.75, bg="transparent")
      dev.off()
    }
  }
}


#-----------------------------> End: PCA using all peptides without NAs

#-----------------------------> Begin: PCA using selected peptides

pca.dat <- dat[features,]

group.vec <- c(
  "T1",
  "T2a",
  "T2B",
  "T2x"
)

col.vec <- c(
  adjustcolor("navy", alpha=0.3),
  adjustcolor("red", alpha=0.3),
  adjustcolor("darkorchid", alpha=0.3),
  adjustcolor("darkgreen", alpha=0.3)
)


pcdat <- prcomp(t(pca.dat), center=TRUE, scale=TRUE)
scores <- pcdat$x
for (i in 1:2){
  for (j in i:2){
    if (i<j){
      XLIM <- c(-max(abs(scores[,i])), max(abs(scores[,i])))
      XLIM <- XLIM+(XLIM*0.1)
      YLIM <- c(-max(abs(scores[,j])), max(abs(scores[,j])))
      YLIM <- YLIM+(YLIM*0.1)
      png(paste(output.path, "/01_pca_", i, "_", j, "_selectedFeatures.png", sep=""), width=3600, height=3600, pointsize=15, res=600)
      plot(scores[group1,i], scores[group1,j], xlab=paste("PC", i, sep=""), ylab=paste("PC", j, sep=""), xlim=XLIM, ylim=YLIM, pch=20, col=col.vec[1], main="PCA", cex=2)
      points(scores[group2,i], scores[group2,j], pch=20, col=col.vec[2], cex=2)
      points(scores[group3,i], scores[group3,j], pch=20, col=col.vec[3], cex=2)
      points(scores[group4,i], scores[group4,j], pch=20, col=col.vec[4], cex=2)
      legend("topleft", legend=group.vec[1:4], col=col.vec[1:4], pch=20, cex=0.75, bg="transparent")
      dev.off()
    }
  }
}

for (i in 1:2){
  for (j in i:2){
    if (i<j){
      XLIM <- c(-max(abs(scores[,i])), max(abs(scores[,i])))
      XLIM <- XLIM+(XLIM*0.1)
      YLIM <- c(-max(abs(scores[,j])), max(abs(scores[,j])))
      YLIM <- YLIM+(YLIM*0.1)
      png(paste(output.path, "/02_pca_", i, "_", j, "_selectedFeatures.png", sep=""), width=3600, height=3600, pointsize=15, res=600)
      plot(scores[group1,i], scores[group1,j], xlab=paste("PC", i, sep=""), ylab=paste("PC", j, sep=""), xlim=XLIM, ylim=YLIM, pch=20, col=col.vec[1], main="PCA", cex=2)
      points(scores[group2,i], scores[group2,j], pch=20, col=col.vec[2], cex=2)
      points(scores[group3,i], scores[group3,j], pch=20, col=col.vec[3], cex=2)
      points(scores[group4,i], scores[group4,j], pch=20, col=col.vec[4], cex=2)
      text(scores[group1,i], scores[group1,j], labels=group1, col="navy", cex=0.4)
      text(scores[group2,i], scores[group2,j], labels=group2, col="red", cex=0.4)
      text(scores[group3,i], scores[group3,j], labels=group3, col="darkorchid", cex=0.4)
      text(scores[group4,i], scores[group4,j], labels=group4, col="darkgreen", cex=0.4)
      legend("topleft", legend=group.vec[1:4], col=col.vec[1:4], pch=20, cex=0.75, bg="transparent")
      dev.off()
    }
  }
}

#-----------------------------> End: PCA using selected peptides

#-----------------------------> Begin: Boxplots

for(i in 1:length(features)){
  idx <- match(features[i], rownames(dat))
  seq <- descriptive.columns2[descriptive.columns2[,1] == features[i],2]
  current.p.value <- formatC(p.values.adj[idx], format = "e", digits = 2)
  png(paste0(output.path,"/boxplot_", i, "_", features[i], ".png"), width=2000, height=2000, pointsize=15, res=300)
  boxplot(dat[features[i],group1], dat[features[i],group2], dat[features[i],group3], dat[features[i],group4], main=paste0(seq, "\np-value: ", current.p.value), names=c(label1, label2, label3, label4), col=col.vec[1:4], cex.axis=1.5, cex.main=1.5)
  dev.off()
}

#-----------------------------> End: Boxplots

#-----------------------------> Begin: Decision Tree

col.list <- list(
  adjustcolor("navy", alpha=0.3),
  adjustcolor("red", alpha=0.3),
  adjustcolor("darkorchid", alpha=0.3),
  adjustcolor("darkgreen", alpha=0.3)
)
rpartdat <- data.frame(classes, t(dat[features,]), stringsAsFactors=TRUE)
for(i in 1:length(features)){
  colnames(rpartdat)[i+1] <- descriptive.columns2[descriptive.columns2[,1] == features[i],2]
}
rpart.model <- rpart(classes~., data=rpartdat, method="class")
png(paste0(output.path,"/rpart-plot.png"), width=2000, height=2000, pointsize=15, res=300)
  rpart.plot(rpart.model, roundint=FALSE,type=2,extra=101,box.palette=col.list)
dev.off()

#-----------------------------> End: Decision Tree