library(ggplot2)

if(!require(psych)){install.packages("psych")}
if(!require(RVAideMemoire)){install.packages("RVAideMemoire")}
if(!require(coin)){install.packages("coin")}
if(!require(reshape2)){install.packages("reshape2")}
if(!require(rcompanion)){install.packages("rcompanion")}

wesloc <- 'C:/Users/Xinming Zhuo/Downloads/clingen_var_2024Jan-alpha-revel-pa3d.csv'
wes <- read.csv(wesloc)


Data <- wes
#Data <- Data[which(Data$am_class != ''),]
Data$diagnosis <- 0
Data[which(Data$Class %in% c('Pathogenic','Likely_Pathogenic')),]$diagnosis <- 1  
Data[which(Data$Class %in% c('Benign','Likely_Benign')),]$diagnosis <- -1  
Data$diagnosis <- factor(Data$diagnosis,
                         levels = c(-1,0,1),
                         labels = c('B/LB','VUS','P/LP'))

Data$diagnosis2 <- -1
Data[which(Data$Class %in% c('Pathogenic','Likely_Pathogenic')),]$diagnosis2 <- 1  
Data[which(Data$Class %in% c('Benign','Likely_Benign')),]$diagnosis2 <- 0 

#gnomad 5percent
gnomad <- read.csv('C:/Users/Xinming Zhuo/Downloads/gnomad_clingen_target_popmax_mis_am_revel_pa3d.csv')
#gnomad <- gnomad[which(gnomad$am_class != ''),]
#gnomad$Variant <- paste(gnomad$chr,gnomad$pos,gnomad$ref,gnomad$alt,sep='-')
gnomad$diagnosis <- 'gnomAD'
gnomad$diagnosis <- factor(gnomad$diagnosis)
gnomad$diagnosis2 <- 0

#ambigous,pathogenic
am_default <- c(0.564)
#BM,BP,PP,PM,PS
am_threshold <- c(0,0.1227,0.1800,0.7191,0.8976,0.9784)

#revel classification
revel_default <- c(0.5)
revel_threshold <- c(0.016,0.183,0.290,0.644,0.773,0.932)

#pa3d classification
pa3d_default <- c(0.821)
#placeholder
pa3d_threshold <- c(0,0.4740,0.5816,0.7910,0.8393,0.9607)

threshold <- rbind(am_threshold,revel_threshold,pa3d_threshold)
default <- c(am_default,revel_default,pa3d_default)

performance <- function(clinvar, default, threshold, caller, output){
  #take in dataframe of score and diagnosis, violin plot and calculate confusion matrix
  colnames(clinvar) <- c('score','diagnosis')
  clinvar$diagnosis <- factor(clinvar$diagnosis,
                              levels = c(0,1),
                              labels = c('B/LB','P/LP'))
  clinvar$Variant <- rownames(clinvar)
  p <- ggplot(clinvar[,c('diagnosis','score')], aes(y=score, x=diagnosis)) +
    geom_violin() + geom_boxplot(width=0.1,outlier.shape= NA) +
    stat_summary(fun.x=mean,geom="point",shape=23,size=2, color='red') +
    geom_jitter(shape=16,position=position_jitter(0.1), size=0.5) +
    geom_hline(yintercept = default, color='gray') +
    geom_hline(yintercept = threshold[1], color='brown') +
    geom_hline(yintercept = threshold[2], color='green') + geom_hline(yintercept=threshold[3], color='cyan') +
    geom_hline(yintercept = threshold[4], color='cyan') + geom_hline(yintercept=threshold[5], color='green') +
    geom_hline(yintercept = threshold[6], color='brown') +
    labs(title=paste(caller))
  
  p
  
  pdf(paste(output,caller,'violin.pdf',sep='.'))
  print(p)
  dev.off()
  
  cmatrix <- data.frame(matrix(ncol=20,nrow=0))
  colnames(cmatrix) <- c('threshold','TP','TN','FP','FN',
                         'TPR_recall','TNR_specificity','PPV_precision','NPV',
                         'FNR','FPR','FDR','FOR',
                         'LR+','LR-','PT','Prevalence',
                         'ACC','F1','DOR')
  auc_cal <- data.frame(matrix(ncol=2,nrow=0))
  colnames(auc_cal) <- c('distance','mean_score')
  
  index <- 1
  
  for (thres in c(default,threshold,c(0:100)/100)){
    
    clinvar$alpha <- 0
    if(nrow(clinvar[which(clinvar$score >= thres),])>0){
      clinvar[which(clinvar$score >= thres),]$alpha <- 1  
    }
    clinvar$alpha <- factor(clinvar$alpha,
                            levels = c(0,1),
                            labels = c('benign','pathogenic'))
    
    ### cross tabulation
    
    Table <- xtabs( ~ diagnosis + alpha, data=clinvar)
    #Table <- Table[c(2,1),]
    
    TP <- Table[2,2]
    TN <- Table[1,1]
    FP <- Table[1,2]
    FN <- Table[2,1]
    
    TPR <- TP/(TP+FN)
    TNR <- TN/(TN+FP)
    PPV <- TP/(TP+FP)
    NPV <- TN/(TN+FN)
    
    FNR <- FN/(FN+TP)
    FPR <- FP/(FP+TN)
    FDR <- FP/(FP+TP)
    FOR <- FN/(FN+TN)
    
    LRp <- TPR/FPR
    LRn <- FNR/TNR
    PT <- sqrt(FPR)/(sqrt(TPR) + sqrt(FPR))
    Prevalence <- (TP+FN)/(TP+FN+FP+TN)
    
    ACC <- (TP+TN)/(TP+FN+FP+TN)
    F1 <- 2*PPV*TPR/(PPV+TPR)
    DOR <- LRp/LRn 
    
    cmatrix[index,] <- c(thres,TP,TN,FP,FN,
                         TPR,TNR,PPV,NPV,
                         FNR,FPR,FDR,FOR,
                         LRp,LRn,PT,Prevalence,
                         ACC,F1,DOR) 
    if (index == 1){
      auc_cal[index,] <- c(cmatrix$FPR[index],cmatrix$TPR_recall[index])
    } else if (index == (length(threshold)+length(default)+1)){
      auc_cal[index,] <- c(cmatrix$FPR[index] - 1,(cmatrix$TPR_recall[index] + 1)/2)
    } else {
      auc_cal[index,] <- c(cmatrix$FPR[index] - cmatrix$FPR[index-1],(cmatrix$TPR_recall[index] + cmatrix$TPR_recall[index-1])/2)
    }
    
    
    index <- index+1
  }

  write.csv(cmatrix, paste(output,caller,'confusion_matrix.csv',sep='.'), row.names = FALSE)
  
  auc_cal$product <- auc_cal$distance * auc_cal$mean_score
  auc_value <- sum(-auc_cal$product[(length(threshold)+length(default)+1):nrow(cmatrix)])
  
  basicplot <- ggplot(cmatrix, aes(x = FPR, y = TPR_recall)) + geom_line() +
    geom_point(data=cmatrix[1,],color='grey') +
    geom_point(data=cmatrix[4,],color='cyan') +
    geom_point(data=cmatrix[3,],color='green') +
    geom_point(data=cmatrix[2,],color='brown') +
    geom_point(data=cmatrix[5,],color='cyan') +
    geom_point(data=cmatrix[6,],color='green') +
    geom_point(data=cmatrix[7,],color='brown') +
    xlim(0,1) + ylim(0,1) + 
    labs(title=paste(caller,"ROC_AUC =",auc_value))
  
  
  basicplot 
  
  #f1plot <- ggplot(cmatrix, aes(x = threshold, y = F1)) + geom_line()
  #  gghighlight(threshold %in% c(0.564,0.7191,0.8976,0.9784),label_key = type) 
  f1plot <- ggplot(cmatrix, aes(x = threshold, y = F1)) + geom_line() +
    geom_point(data=cmatrix[1,],color='grey') +
    geom_point(data=cmatrix[4,],color='cyan') +
    geom_point(data=cmatrix[3,],color='green') +
    geom_point(data=cmatrix[2,],color='brown') +
    geom_point(data=cmatrix[5,],color='cyan') +
    geom_point(data=cmatrix[6,],color='green') +
    geom_point(data=cmatrix[7,],color='brown') +
    xlim(0,1) + ylim(0,1) +
    labs(title=paste(caller,"F1"))
  
  f1plot
  
  pdf(paste(output,caller,'confusion_matrix.pdf',sep='.'))
  print(basicplot)
  print(f1plot)
  dev.off()
  
}

test <- Data[which(Data$diagnosis2 >= 0),c('am_pathogenicity','diagnosis2')]
test1 <- gnomad[which(gnomad$diagnosis2 >= 0),c('am_pathogenicity','diagnosis2')]
test2 <- rbind(test,test1,make.row.names=FALSE)
performance(test[complete.cases(test),], am_default, am_threshold, 'ClingenVar_AM', sub('.csv','',wesloc))
performance(test2[complete.cases(test2),], am_default, am_threshold, 'ClingenVar_gnomAD_AM', sub('.csv','',wesloc))


test <- Data[which(Data$diagnosis2 >= 0),c('REVEL','diagnosis2')]
test1 <- gnomad[which(gnomad$diagnosis2 >= 0),c('REVEL','diagnosis2')]
test2 <- rbind(test,test1,make.row.names=FALSE)
performance(test[complete.cases(test),], revel_default, revel_threshold, 'ClingenVar_REVEL', sub('.csv','',wesloc))
performance(test2[complete.cases(test2),], revel_default, revel_threshold, 'ClingenVar_gnomAD_REVEL', sub('.csv','',wesloc))


test <- Data[which(Data$diagnosis2 >= 0),c('score_PAI3D','diagnosis2')]
test1 <- gnomad[which(gnomad$diagnosis2 >= 0),c('score_PAI3D','diagnosis2')]
test2 <- rbind(test,test1,make.row.names=FALSE)
performance(test[complete.cases(test),], pa3d_default, pa3d_threshold, 'ClingenVar_PAI3D', sub('.csv','',wesloc))
performance(test2[complete.cases(test2),], pa3d_default, pa3d_threshold, 'ClingenVar_gnomAD_PAI3D', sub('.csv','',wesloc))





predict <- function(score,threshold) {
  #class prediction according to calibrated PP3/BP4 score
  predict_class <- score
  predict_class <- NA
  predict_class[which(score >= 0)] <- '4_VUS_Ambiguous' 
  predict_class[which(score >= threshold[4])] <- '5_PP3_Support' 
  predict_class[which(score >= threshold[5])] <- '6_PP3_Moderate' 
  predict_class[which(score >= threshold[6])] <- '7_PP3_Strong' 
  predict_class[which(score <= threshold[3])] <- '3_BP4_Support' 
  predict_class[which(score <= threshold[2])] <- '2_BP4_Moderate' 
  predict_class[which(score <= threshold[1])] <- '1_BP4_Strong'
  return(predict_class)
}

Data$AM_predict <- predict(Data$am_pathogenicity,am_threshold)
Data$REVEL_predict <- predict(Data$REVEL,revel_threshold)
Data$PAI3D_predict <- predict(Data$score_PAI3D,pa3d_threshold)

write.csv(Data, sub('.csv','.predict_class.csv',wesloc), row.names = FALSE)

crosstab <- xtabs( ~ diagnosis + AM_predict, data=Data)
write.csv(crosstab, sub('.csv','.classtab_AM.csv',wesloc), row.names = TRUE)

crosstab <- xtabs( ~ diagnosis + REVEL_predict, data=Data)
write.csv(crosstab, sub('.csv','.classtab_REVEL.csv',wesloc), row.names = TRUE)

crosstab <- xtabs( ~ diagnosis + PAI3D_predict, data=Data)
write.csv(crosstab, sub('.csv','.classtab_PAI3D.csv',wesloc), row.names = TRUE)


individual_plot <- function(clinvar, default, threshold, caller, output){
  #take in dataframe of score and diagnosis, violin plot and calculate confusion matrix
  colnames(clinvar) <- c('score','diagnosis')
  clinvar$score <- as.double(clinvar$score)
  p <- ggplot(clinvar[which(clinvar$diagnosis %in% c('P/LP')),c('diagnosis','score')], aes(y=score, x=diagnosis)) +
    geom_violin() + geom_boxplot(width=0.1,outlier.shape= NA) +
    stat_summary(fun.x=mean,geom="point",shape=23,size=2, color='red') +
    geom_jitter(shape=16,position=position_jitter(0.1), size=0.5) +
    geom_hline(yintercept = default, color='gray') +
    geom_hline(yintercept = threshold[1], color='brown') +
    geom_hline(yintercept = threshold[2], color='green') + geom_hline(yintercept=threshold[3], color='cyan') +
    geom_hline(yintercept = threshold[4], color='cyan') + geom_hline(yintercept=threshold[5], color='green') +
    geom_hline(yintercept = threshold[6], color='brown') +
    labs(title=paste(caller))
  
  b <- ggplot(clinvar[which(clinvar$diagnosis %in% c('B/LB')),c('diagnosis','score')], aes(y=score, x=diagnosis)) +
    geom_violin() + geom_boxplot(width=0.1,outlier.shape= NA) +
    stat_summary(fun.x=mean,geom="point",shape=23,size=2, color='red') +
    geom_jitter(shape=16,position=position_jitter(0.1), size=0.5) +
    geom_hline(yintercept = default, color='gray') +
    geom_hline(yintercept = threshold[1], color='brown') +
    geom_hline(yintercept = threshold[2], color='green') + geom_hline(yintercept=threshold[3], color='cyan') +
    geom_hline(yintercept = threshold[4], color='cyan') + geom_hline(yintercept=threshold[5], color='green') +
    geom_hline(yintercept = threshold[6], color='brown') +
    labs(title=paste(caller))
  
  v <- ggplot(clinvar[which(clinvar$diagnosis %in% c('VUS')),c('diagnosis','score')], aes(y=score, x=diagnosis)) +
    geom_violin() + geom_boxplot(width=0.1,outlier.shape= NA) +
    stat_summary(fun.x=mean,geom="point",shape=23,size=2, color='red') +
    geom_jitter(shape=16,position=position_jitter(0.1), size=0.5) +
    geom_hline(yintercept = default, color='gray') +
    geom_hline(yintercept = threshold[1], color='brown') +
    geom_hline(yintercept = threshold[2], color='green') + geom_hline(yintercept=threshold[3], color='cyan') +
    geom_hline(yintercept = threshold[4], color='cyan') + geom_hline(yintercept=threshold[5], color='green') +
    geom_hline(yintercept = threshold[6], color='brown') +
    labs(title=paste(caller))
  
  g <- ggplot(clinvar[which(clinvar$diagnosis %in% c('gnomAD')),c('diagnosis','score')], aes(y=score, x=diagnosis)) +
    geom_violin() + geom_boxplot(width=0.1,outlier.shape= NA) +
    stat_summary(fun.x=mean,geom="point",shape=23,size=2, color='red') +
    geom_jitter(shape=16,position=position_jitter(0.1), size=0.5) +
    geom_hline(yintercept = default, color='gray') +
    geom_hline(yintercept = threshold[1], color='brown') +
    geom_hline(yintercept = threshold[2], color='green') + geom_hline(yintercept=threshold[3], color='cyan') +
    geom_hline(yintercept = threshold[4], color='cyan') + geom_hline(yintercept=threshold[5], color='green') +
    geom_hline(yintercept = threshold[6], color='brown') +
    labs(title=paste(caller))
  
  pdf(paste(output,caller,'ind_violin.pdf',sep='.'))
  print(p)
  print(b)
  print(v)
  print(g)
  dev.off()
}

test <- Data[,c('am_pathogenicity','diagnosis')]
test1 <- gnomad[,c('am_pathogenicity','diagnosis')]
test2 <- rbind(test,test1,make.row.names=FALSE)
individual_plot(test2[complete.cases(test2),], am_default, am_threshold, 'WES_gnomAD_AM', sub('.csv','',wesloc))

test <- Data[,c('REVEL','diagnosis')]
test1 <- gnomad[,c('REVEL','diagnosis')]
test2 <- rbind(test,test1,make.row.names=FALSE)
individual_plot(test2[complete.cases(test2),], revel_default, revel_threshold, 'WES_gnomAD_REVEL', sub('.csv','',wesloc))

test <- Data[,c('score_PAI3D','diagnosis')]
test1 <- gnomad[,c('score_PAI3D','diagnosis')]
test2 <- rbind(test,test1,make.row.names=FALSE)
individual_plot(test2[complete.cases(test2),], pa3d_default, pa3d_threshold, 'WES_gnomAD_PAI3D', sub('.csv','',wesloc))

gnomad$AM_predict <- predict(gnomad$am_pathogenicity,am_threshold)
gnomad$REVEL_predict <- predict(gnomad$REVEL,revel_threshold)
gnomad$PAI3D_predict <- predict(gnomad$score_PAI3D,pa3d_threshold)

write.csv(Data, sub('.csv','.gnomad.predict_class.csv',wesloc), row.names = FALSE)

crosstab <- xtabs( ~ diagnosis + AM_predict, data=gnomad)
write.csv(crosstab, sub('.csv','.gnomad.classtab_AM.csv',wesloc), row.names = TRUE)

crosstab <- xtabs( ~ diagnosis + REVEL_predict, data=gnomad)
write.csv(crosstab, sub('.csv','.gnomad.classtab_REVEL.csv',wesloc), row.names = TRUE)

crosstab <- xtabs( ~ diagnosis + PAI3D_predict, data=gnomad)
write.csv(crosstab, sub('.csv','.gnomad.classtab_PAI3D.csv',wesloc), row.names = TRUE)



individual_plot <- function(clinvar, default, threshold, caller, output){
  #take in dataframe of score and diagnosis, violin plot and calculate confusion matrix
  colnames(clinvar) <- c('score','diagnosis')
  clinvar$score <- as.double(clinvar$score)
  p <- ggplot(clinvar[which(clinvar$diagnosis %in% c('P/LP')),c('diagnosis','score')], aes(y=score, x=diagnosis)) +
    geom_violin() + geom_boxplot(width=0.1,outlier.shape= NA) +
    stat_summary(fun.x=mean,geom="point",shape=23,size=2, color='red') +
    geom_jitter(shape=16,position=position_jitter(0.1), size=0.5) +
    geom_hline(yintercept = default, color='gray') +
    geom_hline(yintercept = threshold[1], color='brown') +
    geom_hline(yintercept = threshold[2], color='green') + geom_hline(yintercept=threshold[3], color='cyan') +
    geom_hline(yintercept = threshold[4], color='cyan') + geom_hline(yintercept=threshold[5], color='green') +
    geom_hline(yintercept = threshold[6], color='brown') +
    labs(title=paste(caller))
  
  b <- ggplot(clinvar[which(clinvar$diagnosis %in% c('B/LB')),c('diagnosis','score')], aes(y=score, x=diagnosis)) +
    geom_violin() + geom_boxplot(width=0.1,outlier.shape= NA) +
    stat_summary(fun.x=mean,geom="point",shape=23,size=2, color='red') +
    geom_jitter(shape=16,position=position_jitter(0.1), size=0.5) +
    geom_hline(yintercept = default, color='gray') +
    geom_hline(yintercept = threshold[1], color='brown') +
    geom_hline(yintercept = threshold[2], color='green') + geom_hline(yintercept=threshold[3], color='cyan') +
    geom_hline(yintercept = threshold[4], color='cyan') + geom_hline(yintercept=threshold[5], color='green') +
    geom_hline(yintercept = threshold[6], color='brown') +
    labs(title=paste(caller))
  
  v <- ggplot(clinvar[which(clinvar$diagnosis %in% c('VUS')),c('diagnosis','score')], aes(y=score, x=diagnosis)) +
    geom_violin() + geom_boxplot(width=0.1,outlier.shape= NA) +
    stat_summary(fun.x=mean,geom="point",shape=23,size=2, color='red') +
    geom_jitter(shape=16,position=position_jitter(0.1), size=0.5) +
    geom_hline(yintercept = default, color='gray') +
    geom_hline(yintercept = threshold[1], color='brown') +
    geom_hline(yintercept = threshold[2], color='green') + geom_hline(yintercept=threshold[3], color='cyan') +
    geom_hline(yintercept = threshold[4], color='cyan') + geom_hline(yintercept=threshold[5], color='green') +
    geom_hline(yintercept = threshold[6], color='brown') +
    labs(title=paste(caller))
  
  g <- ggplot(clinvar[which(clinvar$diagnosis %in% c('gnomAD')),c('diagnosis','score')], aes(y=score, x=diagnosis)) +
    geom_violin() + geom_boxplot(width=0.1,outlier.shape= NA) +
    stat_summary(fun.x=mean,geom="point",shape=23,size=2, color='red') +
    geom_jitter(shape=16,position=position_jitter(0.1), size=0.5) +
    geom_hline(yintercept = default, color='gray') +
    geom_hline(yintercept = threshold[1], color='brown') +
    geom_hline(yintercept = threshold[2], color='green') + geom_hline(yintercept=threshold[3], color='cyan') +
    geom_hline(yintercept = threshold[4], color='cyan') + geom_hline(yintercept=threshold[5], color='green') +
    geom_hline(yintercept = threshold[6], color='brown') +
    labs(title=paste(caller))
  
  pdf(paste(output,caller,'ind_violin.pdf',sep='.'))
  print(p)
  print(b)
  print(v)
  print(g)
  dev.off()
}

test <- Data[,c('am_pathogenicity','diagnosis')]
test1 <- gnomad[,c('am_pathogenicity','diagnosis')]
test2 <- rbind(test,test1,make.row.names=FALSE)
individual_plot(test2[complete.cases(test2),], am_default, am_threshold, 'WES_gnomAD_AM', sub('.csv','',wesloc))

test <- Data[,c('REVEL','diagnosis')]
test1 <- gnomad[,c('REVEL','diagnosis')]
test2 <- rbind(test,test1,make.row.names=FALSE)
individual_plot(test2[complete.cases(test2),], revel_default, revel_threshold, 'WES_gnomAD_REVEL', sub('.csv','',wesloc))

test <- Data[,c('score_PAI3D','diagnosis')]
test1 <- gnomad[,c('score_PAI3D','diagnosis')]
test2 <- rbind(test,test1,make.row.names=FALSE)
individual_plot(test2[complete.cases(test2),], pa3d_default, pa3d_threshold, 'WES_gnomAD_PAI3D', sub('.csv','',wesloc))

