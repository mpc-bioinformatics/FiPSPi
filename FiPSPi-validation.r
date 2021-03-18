##########################################################################
#                          COPYRIGHT & SUPPORT                           #
##########################################################################
#
# FiPSPi-validation: R code for the validation of the results of FiPSPi
# with independent data from the publication (Murgia_et_al_data.txt):
#     Murgia M, (...), Mann M.,
#     Single muscle fiber proteomics reveals unexpected
#     mitochondrial specialization.
#     EMBO Rep. 2015 Mar;16(3):387-95.
#     doi: 10.15252/embr.201439757.
#     Epub 2015 Feb 2. PMID: 25643707.
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
# version of FiPSPi-validation can be found here:
# https://github.com/mpc-bioinformatics/FiPSPi
#
##########################################################################
#                            USER SETTINGS                               #
##########################################################################

### set path to current working directory
cwd <- "C:/UNI/Publikationen/2021.XX.XX_writing_Britta_Eggers_Fasertypen/FiPSPi"

### set path to data file (relative path to the working directory specified above)
data.path <- paste0(cwd, "/data/Murgia_et_al_data.txt")

### set path to selected features file (relative path to the working directory specified above)
selected.features.path <- paste0(cwd, "/FiPSPi_results/selected_features.txt")

### set output.path where all created output and graphics are saved
output.path <- paste0(cwd, "/FiPSPi-validation_results")

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
#                            FUNCTIONS                                   #
##########################################################################

#+++++++++++++++++++++++++++ classify.rf.evaluation +++++++++++++++++++++++++++++++++++++++
classify.rf.evaluation <- function(iteration=NULL, trainset=NULL, testset=NULL, classes_train=NULL, 
                                   classes_test=NULL, label1="A", label2="B", label3 = "C", output.path=NULL, ...){
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
  FN <- confusion[1,2]+confusion[1,3]
  FP <- confusion[2,1]+confusion[3,1]
  TN <- confusion[2,2]+confusion[2,3]+confusion[3,2]+confusion[3,3]
  ACCURACY1 <- {TP+TN}/{TP+FP+TN+FN}
  SENSITIVITY1 <- TP/{TP+FN}   #=TPR
  SPECIFICITY1 <- TN/{TN+FP}  #=(1-FPR) 
  
  ### class 2:
  TP <- confusion[2,2]
  FN <- confusion[2,1]+confusion[2,3]
  FP <- confusion[1,2]+confusion[3,2]
  TN <- confusion[1,1]+confusion[1,3]+confusion[3,1]+confusion[3,3]
  ACCURACY2 <- {TP+TN}/{TP+FP+TN+FN}
  SENSITIVITY2 <- TP/{TP+FN}   #=TPR
  SPECIFICITY2 <- TN/{TN+FP}  #=(1-FPR) 
  
  ### class 3:
  TP <- confusion[3,3]
  FN <- confusion[3,1]+confusion[3,2]
  FP <- confusion[1,3]+confusion[2,3]
  TN <- confusion[1,1]+confusion[1,2]+confusion[2,1]+confusion[2,2]
  ACCURACY3 <- {TP+TN}/{TP+FP+TN+FN}
  SENSITIVITY3 <- TP/{TP+FN}   #=TPR
  SPECIFICITY3 <- TN/{TN+FP}  #=(1-FPR)
  
  ACCURACY_total <- sum(diag(confusion))/ sum(confusion) # total Accuracy
  SENSITIVITY_total <- (SENSITIVITY1 + SENSITIVITY2 + SENSITIVITY3) / 3
  SPECIFICITY_total <- (SPECIFICITY1 + SPECIFICITY2 + SPECIFICITY3) / 3
  
  results <- c(acc_total = ACCURACY_total, sens_total = SENSITIVITY_total, spec_total = SPECIFICITY_total, 
               acc1 = ACCURACY1,  sens1 = SENSITIVITY1, spec1 = SPECIFICITY1,
               acc2 = ACCURACY2,  sens2 = SENSITIVITY2, spec2 = SPECIFICITY2,
               acc3 = ACCURACY3,  sens3 = SENSITIVITY3, spec3 = SPECIFICITY3)

  return(results)
}
#+++++++++++++++++++++++++++ classify.rf.evaluation +++++++++++++++++++++++++++++++++++++++

##########################################################################
#                               MAIN                                     #
##########################################################################

#-----------------------------> Begin: Data Import & Pre-Processing

dat <- read.table(data.path, sep="\t", header=TRUE, na.strings = "NA")

label1 <- "T1" 
#label2 <- "T2a"
label3 <- "T2b"
label4 <- "T2x"

# Exclusion of 6 samples due to following criteria:
# -	Samples without predominant Myh-isoform and/or without unambiguous fiber-type in the column 'Fiber type, renamed'
# -	Samples with at least one missing value for the peptides to be validated
# -	After the first two exclusion criteria group 'T2a' has only 3 samples --> group too small --> exclusion of complete group
# -	Samples with weakly predominant Myh-isoform and which were outliers in PCA
dat <- dat[,-match(c("T1_9", "T2b_1", "T2b_3", "T2b_6", "T2x_5", "T2x_7"), colnames(dat))]

group1 <- grep("T1", colnames(dat), value=TRUE)
#group2 <- grep("T2a", colnames(dat), value=TRUE)
group3 <- grep("T2b", colnames(dat), value=TRUE)
group4 <- grep("T2x", colnames(dat), value=TRUE)
n1 <- length(group1)
#n2 <- length(group2)
n3 <- length(group3)
n4 <- length(group4) 

descriptive.columns <- dat[,1:5]
selected.features <- read.table(selected.features.path, sep="\t", header=TRUE, na.strings = "NA")
rownames(descriptive.columns) <- selected.features[,"Sequence"]
#rownames(descriptive.columns) <- c("Q5SX40_peptide40", "O88990_peptide38")
rownames(dat) <- rownames(descriptive.columns)
write.table(x=cbind(rownames(descriptive.columns),descriptive.columns), file=paste0(output.path, "/descriptive.columns.txt"), sep="\t", col.names=TRUE, row.names=FALSE)

dat <- log2(data.matrix(dat))
dat <- dat[,c(group1,group3,group4)]
write.table(x=cbind(rownames(dat),dat), file=paste0(output.path, "/original_data.txt"), sep="\t", col.names=TRUE, row.names=FALSE)

#-----------------------------> Begin: Data Import & Pre-Processing

#-----------------------------> Begin: Diff. Analysis

classes <- gsub("\\_\\d+", "", colnames(dat))
rawdat <- 2^dat
p.values <- vector(mode="numeric", length=nrow(dat))
max.mean.ratios <- vector(mode="numeric", length=nrow(dat))
for(i in 1:nrow(dat)){
  p.values[i] <- summary(aov(dat[i,] ~ classes))[[1]][["Pr(>F)"]][1]
  
  mean.ratio1 <- max( mean(rawdat[i,group1])/mean(c(rawdat[i,group3], rawdat[i,group4])), mean(c(rawdat[i,group3], rawdat[i,group4]))/mean(rawdat[i,group1]) )
  mean.ratio3 <- max( mean(rawdat[i,group3])/mean(c(rawdat[i,group1], rawdat[i,group4])), mean(c(rawdat[i,group1], rawdat[i,group4]))/mean(rawdat[i,group3]) )
  mean.ratio4 <- max( mean(rawdat[i,group4])/mean(c(rawdat[i,group1], rawdat[i,group3])), mean(c(rawdat[i,group1], rawdat[i,group3]))/mean(rawdat[i,group4]) )
  max.mean.ratios[i] <- max(mean.ratio1, mean.ratio3, mean.ratio4)
}
p.values.adj <- p.adjust(p.values, method="fdr")
names(p.values.adj) <- rownames(dat)
diff.analysis.output <- cbind(names(p.values.adj), p.values, p.values.adj, max.mean.ratios)
colnames(diff.analysis.output) <- c("Feature", "P-value (Anova)", "Adj. p-value", "Max. mean ratio")
write.table(file=paste0(output.path, "/diff_analysis_output.txt"), x=diff.analysis.output, 
            row.names=FALSE, col.names=TRUE, sep="\t")

#-----------------------------> End: Diff. Analysis

#-----------------------------> Begin: Validation with Test Set

### 1000 times train/test split for model validation
testruns <- 1000
features <- rownames(dat)
all.idx <- colnames(dat)
#perf.final <- vector("list", testruns)
perf.final <- matrix(nrow =  testruns, ncol = 12)
colnames(perf.final) = c("acc_total", "sens_total", "spec_total",
  "acc1",  "sens1", "spec1",
  "acc2",  "sens2", "spec2",
  "acc3",  "sens3", "spec3")

for(i in 1:testruns){
  cat("final classification: ", i, "\r")
  test.idx1 <- sample(group1, 3)  # n1=10
  #test.idx2 <- sample(group2, 1)  # n2=3
  test.idx3 <- sample(group3, 5)  # n3=15
  test.idx4 <- sample(group4, 3)  # n4=10
  #test.idx <- c(test.idx1, test.idx2, test.idx3, test.idx4)
  test.idx <- c(test.idx1, test.idx3, test.idx4)
  
  train.idx <- setdiff(all.idx, test.idx)
  
  classi <- classify.rf.evaluation(iteration=i, trainset=dat[features,train.idx,drop=FALSE], 
                                   testset=dat[features,test.idx,drop=FALSE], 
                                   classes_train=as.factor(gsub("\\_\\d+", "", train.idx)), 
                                   classes_test=as.factor(gsub("\\_\\d+", "", test.idx)), 
                                   label1="T1", label2="T2b", label3="T2x", output.path=NULL, 
                                   nrounds = best.nrounds, max_depth = best.max_depth)
  perf.final[i,] <- classi
}

write.table(perf.final, file = paste0(output.path, "/validation.txt"), row.names = FALSE, sep = "\t")

#-----------------------------> End: Validation with Test Set

#-----------------------------> Begin: PCA

pca.dat <- dat[features,]

group.vec <- c(
  "T1",
  "T2b",
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
      points(scores[group3,i], scores[group3,j], pch=20, col=col.vec[3], cex=2)
      points(scores[group4,i], scores[group4,j], pch=20, col=col.vec[4], cex=2)
      legend("topleft", legend=group.vec[1:3], col=col.vec[c(1,3,4)], pch=20, cex=0.75, bg="transparent")
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
      points(scores[group3,i], scores[group3,j], pch=20, col=col.vec[3], cex=2)
      points(scores[group4,i], scores[group4,j], pch=20, col=col.vec[4], cex=2)
      text(scores[group1,i], scores[group1,j], labels=group1, col="navy", cex=0.4)
      text(scores[group3,i], scores[group3,j], labels=group3, col="darkorchid", cex=0.4)
      text(scores[group4,i], scores[group4,j], labels=group4, col="darkgreen", cex=0.4)
      legend("topleft", legend=group.vec[1:3], col=col.vec[c(1,3,4)], pch=20, cex=0.75, bg="transparent")
      dev.off()
    }
  }
}

#-----------------------------> End: PCA

#-----------------------------> Begin: Boxplots

for(i in 1:length(features)){
  idx <- match(features[i], rownames(dat))
  seq <- descriptive.columns[features[i],1]
  current.p.value <- formatC(p.values.adj[idx], format = "e", digits = 2)
  png(paste0(output.path,"/boxplot_", i, "_", features[i], ".png"), width=2000, height=2000, pointsize=15, res=300)
    boxplot(dat[features[i],group1], dat[features[i],group3], dat[features[i],group4], main=paste0(seq, "\np-value: ", current.p.value), names=c(label1, label3, label4), col=col.vec[c(1,3,4)], cex.axis=1.5, cex.main=1.5)
  dev.off()
}

#-----------------------------> End: Boxplots

#-----------------------------> Begin: Decision Tree

col.list <- list(
  adjustcolor("navy", alpha=0.3),
  adjustcolor("darkorchid", alpha=0.3),
  adjustcolor("darkgreen", alpha=0.3)
)
rpartdat <- data.frame(classes, t(dat[features,]), stringsAsFactors=TRUE)
for(i in 1:length(features)){
  colnames(rpartdat)[i+1] <- descriptive.columns[features[i],1]
}
rpart.model <- rpart(classes~., data=rpartdat, method="class")
png(paste0(output.path,"/rpart-plot.png"), width=2000, height=2000, pointsize=15, res=300)
    rpart.plot(rpart.model, roundint=FALSE,type=2,extra=101,box.palette=col.list)
dev.off()

#-----------------------------> End: Decision Tree