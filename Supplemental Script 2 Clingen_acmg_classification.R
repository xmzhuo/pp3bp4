#automatically classify variant with evidence

table <- 'C:/Users/Xinming Zhuo/Downloads/clingen_var_2024Jan-alpha-revel-pa3d.csv'
clingen <- read.csv(table)

am_default <- c(0.34,0.564)
#BM,BP,PP,PM,PS
am_threshold <- c(0,0.1227,0.1800,0.7191,0.8976,0.9784)

#revel classification
revel_default <- c(0.5,0.75)
revel_threshold <- c(0.016,0.183,0.290,0.644,0.773,0.932)

#pa3d classification
pa3d_default <- c(0.821)
#placeholder
pa3d_threshold <- c(0,0.4740,0.5816,0.7910,0.8393,0.9607)

predict_class <- function(score,threshold){
  #class variant according to prediction threshold
  pclass <- score
  pclass[which(score >= threshold[4])] <- 'PP3_Supporting' 
  pclass[which(score >= threshold[5])] <- 'PP3_Moderate' 
  pclass[which(score >= threshold[6])] <- 'PP3_Strong' 
  pclass[which(score <= threshold[3])] <- 'BP4_Supporting' 
  pclass[which(score <= threshold[2])] <- 'BP4_Moderate'
  pclass[which(score <= threshold[2])] <- 'BP4_Strong'
  return(pclass)
}


count_evidence <-function(current_list){
  #count current evidence by class and return a one row dataframe
  PVS <- length(grep("PVS|P_Very_Strong",current_list))
  PS <- length(grep("PS|P_Strong",current_list))
  PM <- length(grep("PM|P_Moderate",current_list))
  PP <- length(grep("PP|P_Supporting",current_list))
  BP <- length(grep("BP|B_Supporting|B_Moderate",current_list))
  BS <- length(grep("BS|B_Strong|B_Very_Strong",current_list))
  BA <- length(grep("BA",current_list))
  return(data.frame(PVS=PVS,PS=PS,PM=PM,PP=PP,BP=BP,BS=BS,BA=BA))
}

acmg_class <- function(acmg){
  
  #acmg classification logic, feed a dataframe reurn classification
  #pathogenic
  acmg$ACMG_P <- "Uncertain_Significance"
  if(acmg$PVS >= 1){
    if(acmg$PS >=1 || acmg$PM >=2 || (acmg$PM >= 1 && acmg$PP >= 1) || acmg$PP >= 2) {
      acmg$ACMG_P <- "Pathogenic"
    } else if(acmg$PM == 1){
      acmg$ACMG_P <- "Likely_Pathogenic"
    } else {
      acmg$ACMG_P <- "Uncertain_Significance"
    }
    
  } else if(acmg$PS >= 1) {
    if(acmg$PS >= 2) {
        acmg$ACMG_P <- "Pathogenic"
    } else if(acmg$PS == 1) {
      if(acmg$PM >= 3 || (acmg$PM ==2 && acmg$PP >=2) || (acmg$PM ==1 && acmg$PP >=4)){
        acmg$ACMG_P <- "Pathogenic"
      } else if(acmg$PM >= 1){
        acmg$ACMG_P <- "Likely_Pathogenic"
      } else if(acmg$PP >= 2){
        acmg$ACMG_P <- "Likely_Pathogenic"
      } else {
        acmg$ACMG_P <- "Uncertain_Significance"  
      }
    } else {
      acmg$ACMG_P <- "Uncertain_Significance"
    }
    
  } else if(acmg$PM >=1){
    if(acmg$PM >=3){
      acmg$ACMG_P <- "Likely_Pathogenic"
    } else if(acmg$PM ==2 && acmg$PP >=2){
      acmg$ACMG_P <- "Likely_Pathogenic"
    } else if(acmg$PM ==1 && acmg$PP >=4){
      acmg$ACMG_P <- "Likely_Pathogenic"
    } else {
      acmg$ACMG_P <- "Uncertain_Significance"
    }
    
  } else {
    acmg$ACMG_P <- "Uncertain_Significance"
  }
  
  
  #Benign
  acmg$ACMG_B <- "Uncertain_Significance"
  if(acmg$BA >= 1){
    acmg$ACMG_B <- "Benign"
  } else if(acmg$BS >=1){
    if(acmg$BS >=2){
      acmg$ACMG_B <- "Benign"
    }else if(acmg$BS ==1 && acmg$BP >=1){
      acmg$ACMG_B <- "Likely_Benign"
    }else{
      acmg$ACMG_B <- "Uncertain_Significance"
    }
  } else if(acmg$BP >=2){
    acmg$ACMG_B <- "Likely_Benign"
  } else {
    acmg$ACMG_B <- "Uncertain_Significance"
  }
  
  #if contradictory
  if(acmg$ACMG_P == "Pathogenic"){
    if(acmg$ACMG_B == "Benign"){
      acmg$ACMG <- "Contradictory_Uncertain_Significance"
    }else if(acmg$ACMG_B == "Likely_Benign" || acmg$BA + acmg$BS + acmg$BP > 0){
      acmg$ACMG <- "Contradictory_Likely_Pathogenic"
    }else{
      acmg$ACMG <- acmg$ACMG_P
    }
  }else if(acmg$ACMG_P == "Likely_Pathogenic"){
    if(acmg$ACMG_B == "Benign"){
      acmg$ACMG <- "Contradictory_Likely_Benign"
    }else if(acmg$ACMG_B == "Likely_Benign" || acmg$BA + acmg$BS + acmg$BP > 0){
      acmg$ACMG <- "Contradictory_Uncertain_Significance"
    }else{
      acmg$ACMG <- acmg$ACMG_P
    }
  }else if(acmg$PVS + acmg$PS + acmg$PM + acmg$PP >0){
    if(acmg$ACMG_B == "Benign"){
      acmg$ACMG <- "Contradictory_Likely_Benign"
    }else if(acmg$ACMG_B == "Likely_Benign"){
      acmg$ACMG <- "Contradictory_Uncertain_Significance"
    }else{
      acmg$ACMG <- acmg$ACMG_P
    }
  }else{
    acmg$ACMG <- acmg$ACMG_B
  }
    
  return(acmg$ACMG)
}



clingen_new <- clingen

clingen_new$am_pp3 <- predict_class(clingen$am_pathogenicity,am_threshold)
clingen_new$revel_pp3 <- predict_class(clingen$REVEL,revel_threshold)
clingen_new$pai3d_pp3 <- predict_class(clingen$score_PAI3D,pai3d_threshold)

i=1
             
for (i in 1:nrow(clingen)){
  
  #evidence <- "PS3;PM3;PP4_Moderate;PM2"
  evidence_list <- as.list(strsplit(clingen$evidence[i],";"))[[1]]
  acmg15 <- count_evidence(evidence_list)
  clingen_new$ACMG15_Class[i] <- acmg_class(acmg15)
  clingen_new$ACMG15_evidence[i] <- paste(acmg15,collapse='-')
  
  #update evidence with current classification
  current_list <- sub("^B.*[1-9]_","B_", sub("^P.*[1-9]_","P_",evidence_list))
  acmgclingen <- count_evidence(current_list)
  clingen_new$ACMGchk_Class[i] <- acmg_class(acmgclingen)
  clingen_new$ACMGchk_evidence[i] <- paste(acmgclingen,collapse='-')
  
  #update evidence without PP3 and BP4
  temp_list <- evidence_list[!grepl("PP3|BP4",evidence_list)]
  xpp3_list <- sub("^B.*[1-9]_","B_", sub("^P.*[1-9]_","P_",temp_list))
  acmgxpp3 <- count_evidence(xpp3_list)
  clingen_new$xPP3BP4_Class[i] <- acmg_class(acmgxpp3)
  clingen_new$xPP3BP4_evidence[i] <- paste(acmgxpp3,collapse='-')
  
  #add pp3/bp4 evidence AM
  addpp3_list <- sub("^B.*[1-9]_","B_", sub("^P.*[1-9]_","P_",c(temp_list,clingen_new$am_pp3[i])))
  acmgaddpp3 <- count_evidence(addpp3_list)
  clingen_new$am_PP3BP4_Class[i] <- acmg_class(acmgaddpp3)
  clingen_new$am_PP3BP4_evidence[i] <- paste(acmgaddpp3,collapse='-')
  
  addpp3_list <- sub("^B.*[1-9]_","B_", sub("^P.*[1-9]_","P_",c(temp_list,clingen_new$revel_pp3[i])))
  acmgaddpp3 <- count_evidence(addpp3_list)
  clingen_new$revel_PP3BP4_Class[i] <- acmg_class(acmgaddpp3)
  clingen_new$revel_PP3BP4_evidence[i] <- paste(acmgaddpp3,collapse='-')
  
  addpp3_list <- sub("^B.*[1-9]_","B_", sub("^P.*[1-9]_","P_",c(temp_list,clingen_new$pai3d_pp3[i])))
  acmgaddpp3 <- count_evidence(addpp3_list)
  clingen_new$pai3d_PP3BP4_Class[i] <- acmg_class(acmgaddpp3)
  clingen_new$pai3d_PP3BP4_evidence[i] <- paste(acmgaddpp3,collapse='-')
  
  i=i+1
  
}

write.csv(clingen_new, sub('.csv','.reclass.csv',table), row.names = FALSE)

