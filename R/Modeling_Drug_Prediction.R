library(sva)
library(caret)
library(glmnet)
library(pROC)
library(dplyr)
library(RSpectra)
seed = 20220511
load(file = "/Users/qingzliu/Documents/Thesis1/Data/Data_Preprocess_Outputs.RData")
load(file = "/Users/qingzliu/Documents/Thesis1/Data/fast.jive.functions.RData")

#### Traditional Models (Preparation) ####
# Drug List
drug_list = data.frame(GDSC_drug = c("Afatinib", "Bleomycin","Cetuximab","Cisplatin","Cisplatin","Cyclophosphamide",
                                     "Docetaxel","Doxorubicin","Etoposide","5-Fluorouracil", "Gemcitabine", "Irinotecan",
                                     "Oxaliplatin", "Paclitaxel", "Pemetrexed", "Temozolomide", "Vinorelbine"), 
                       TCGA_drug = c("Trastuzumab", "Bleomycin", "Cetuximab", "Cisplatin", "Carboplatin", "Cyclophosphamide", 
                                     "Docetaxel", "Doxorubicin", "Etoposide", "Fluorouracil", "Gemcitabine", "Irinotecan",
                                     "Oxaliplatin", "Paclitaxel", "Pemetrexed", "Temozolomide", "Vinorelbine"))
# Preparation for Elastic Net
fit_enet <- function(train_geno,train_pheno) {
  Data = data.frame(cbind(train_pheno, train_geno))
  colnames(Data)[1] = "pheno"
  enet_model = train(pheno ~., data = Data, method = "glmnet",
                     trControl = trainControl("cv", number = 10),
                     tuneLength = 10) 
  tune = enet_model$bestTune
  model_enet = glmnet(as.matrix(train_geno),train_pheno,
                      alpha=tune$alpha,
                      lambda=tune$lambda)
  return(model_enet)
}

read_GDSC_drug_data = function(Drug_name) {
  GDSC_response_drug = GDSC_response%>%filter(DRUG_NAME==Drug_name)%>%group_by(SANGER_MODEL_ID)%>%summarize(mean_AUC = mean(AUC, na.rm = TRUE))
  GDSC_response_drug = as.data.frame(GDSC_response_drug)
  rownames(GDSC_response_drug) = CCLE_SIDM_ACH_match$BROAD_ID[match(GDSC_response_drug$SANGER_MODEL_ID, CCLE_SIDM_ACH_match$model_id)]
  GDSC_response_drug = GDSC_response_drug[order(rownames(GDSC_response_drug)), ]
  return(GDSC_response_drug)
}

## Preparation for Combat+Elastic Net
CCLE_Tumor_exp_norm = rbind(CCLE_exp_norm, Tumor_exp_norm)
CCLE_Tumor_exp_combat = ComBat(t(CCLE_Tumor_exp_norm), batch = c(rep(1,nrow(CCLE_exp_norm)), rep(2, nrow(Tumor_exp_norm))))
CCLE_Tumor_exp_combat = t(CCLE_Tumor_exp_combat)
CCLE_exp_combat = CCLE_Tumor_exp_combat[1:nrow(CCLE_exp_norm), ]
Tumor_exp_combat = CCLE_Tumor_exp_combat[(nrow(CCLE_exp_norm)+1):nrow(CCLE_Tumor_exp_combat), ]

## Preparation for PCjoint
CCLE_exp_norm_temp = CCLE_exp_norm/norm(CCLE_exp_norm, type = "f") 
Tumor_exp_norm_temp = Tumor_exp_norm/norm(Tumor_exp_norm, type = "f")
CCLE_Tumor_Joint_svd = svd(rbind(CCLE_exp_norm_temp, Tumor_exp_norm_temp))
cv_PCjoint = function(ranks, fold_n) {
  mse_all = vector()
  for (rank in ranks) {
    CP_score_hat = CCLE_Tumor_Joint_svd$u[,1:rank]%*%diag(CCLE_Tumor_Joint_svd$d[1:rank])
    CCLE_score_hat = CP_score_hat[1:nrow(CCLE_exp_norm), ]
    rownames(CCLE_score_hat) = Cell_Line_Names_Exp
    C_score_hat = CCLE_score_hat[Cell_Line_Names_drug, ]
    mse_vec = vector()
    idx = c(1:nrow(C_score_hat))
    for (k in 1:fold_n) {
      if (k == fold_n) {
        testing_index = idx[((k-1)*(nrow(C_score_hat)%/%fold_n)+1) : nrow(C_score_hat)]
      } else {
        testing_index = idx[((k-1)*(nrow(C_score_hat)%/%fold_n)+1) : (k*(nrow(C_score_hat)%/%fold_n))]
      }
      train_X = C_score_hat[-testing_index, ]
      train_Y = C_drug_AUC[-testing_index]
      test_X = C_score_hat[testing_index, ]
      test_Y = C_drug_AUC[testing_index]
      set.seed(seed)
      enet_model = fit_enet(train_X, train_Y)
      yhat.test.PCjoint = predict(enet_model, newx = test_X)
      mse = mean((yhat.test.PCjoint - test_Y)^2)  ### Fix
      mse_vec = c(mse_vec, mse)
    }
    mean_mse_vec = mean(mse_vec)
    mse_all = c(mse_all, mean_mse_vec)
  }
  best_rank = ranks[which(mse_all == min(mse_all))]
  return(best_rank)
}

#### Implement different traditional methods across all drugs #### 
results_traditon_allDrugs1 = list()
irx_best_rank = c()
for (i in c(9,17,11,13,15,1,2,3,4,5,6,7,8,10,12,14,16)) {
  CCLE_response_drug = read_GDSC_drug_data(Drug_name = drug_list[i,1])
  TCGA_response_drug = read_TCGA_drug_data(Drug_name = drug_list[i,2])
  Cell_Line_Names_Exp = rownames(CCLE_exp_norm)
  #Cell_Line_Names_drug = intersect(rownames(CCLE_exp_norm), rownames(CCLE_response_drug))
  Cell_Line_Names_drug = intersect(rownames(CCLE_response_drug), rownames(CCLE_exp_norm))
  
  C_exp_norm = CCLE_exp_norm[Cell_Line_Names_drug, ]
  TCGA_Names_Exp = rownames(Tumor_exp_norm)
  TCGA_Uni_exp_norm = Tumor_exp_norm[panTCGA_anno$th_sampleid,]
  TCGA_Uni_exp_norm = TCGA_Uni_exp_norm[!duplicated(panTCGA_anno$site_donor_id), ]
  rownames(TCGA_Uni_exp_norm) = panTCGA_anno$site_donor_id[!duplicated(panTCGA_anno$site_donor_id)]
  P_exp_norm = TCGA_Uni_exp_norm[rownames(TCGA_response_drug),]
  
  C_drug_AUC = CCLE_response_drug[Cell_Line_Names_drug, ]%>%select(mean_AUC)
  C_drug_AUC = c(as.matrix(C_drug_AUC))
  
  P_drug_bin = TCGA_response_drug%>%select(measure_of_response)
  P_drug_bin = as.matrix(P_drug_bin)
  for (j in 1:nrow(P_drug_bin)) {
    if (P_drug_bin[j] == "Complete Response" || P_drug_bin[j] == "Partial Response") {
      P_drug_bin[j] = "Response"
    } else {P_drug_bin[j] = "Non-Response"}
  }
  
  ## Elastic Net
  set.seed(seed)
  Data = data.frame(cbind(C_drug_AUC, C_exp_norm))
  colnames(Data)[1] = "pheno"
  enet_model = train(pheno ~., data = Data, method = "glmnet",
                     trControl = trainControl("cv", number = 10),
                     tuneLength = 20) 
  tune = enet_model$bestTune
  enet_model = glmnet(C_exp_norm,C_drug_AUC,
                      alpha=tune$alpha,
                      lambda=tune$lambda)
  
  #enet_model = fit_enet(C_exp_norm, C_drug_AUC)
  yhat.test.EN = predict(enet_model, newx = P_exp_norm)
  
  ## Combat+Elastic Net
  C_exp_combat = CCLE_exp_combat[Cell_Line_Names_drug, ]
  TCGA_Uni_exp_combat = Tumor_exp_combat[panTCGA_anno$th_sampleid,]
  TCGA_Uni_exp_combat = TCGA_Uni_exp_combat[!duplicated(panTCGA_anno$site_donor_id), ]
  rownames(TCGA_Uni_exp_combat) = panTCGA_anno$site_donor_id[!duplicated(panTCGA_anno$site_donor_id)]
  P_exp_combat= TCGA_Uni_exp_combat[rownames(TCGA_response_drug),]
  Data = data.frame(cbind(C_drug_AUC, C_exp_combat))
  colnames(Data)[1] = "pheno"
  set.seed(seed)
  enet_model = train(pheno ~., data = Data, method = "glmnet",
                     trControl = trainControl("cv", number = 10),
                     tuneLength = 20) 
  tune = enet_model$bestTune
  enet_model = glmnet(C_exp_combat,C_drug_AUC,
                      alpha=tune$alpha,
                      lambda=tune$lambda)
  yhat.test.CombatEN = predict(enet_model, newx = P_exp_combat)
  
  auc(roc(c(P_drug_bin), c(yhat.test.CombatEN), direction=">"))
  # PCjoint+Elastic Net
  best_rank = cv_PCjoint(c(10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,94), fold_n = 10)
  #best_rank = cv_PCjoint(c(10,15,20,25,30,35,40,45,50,55,60), fold_n = 10)
  CP_score_hat = CCLE_Tumor_Joint_svd$u[,1:best_rank]%*%diag(CCLE_Tumor_Joint_svd$d[1:best_rank])
  CCLE_score_hat = CP_score_hat[1:nrow(CCLE_exp_norm), ]
  rownames(CCLE_score_hat) = Cell_Line_Names_Exp
  C_score_hat = CCLE_score_hat[Cell_Line_Names_drug, ] * norm(CCLE_exp_norm, type = "f") 
  Tumor_score_hat = CP_score_hat[(nrow(CCLE_exp_norm)+1):nrow(CP_score_hat), ]
  rownames(Tumor_score_hat) = Tumor_anno$th_sampleid
  TCGA_score_hat = Tumor_score_hat[panTCGA_anno$th_sampleid,]
  rownames(TCGA_score_hat) = panTCGA_anno$th_sampleid
  TCGA_Uni_score_hat = TCGA_score_hat[!duplicated(panTCGA_anno$site_donor_id), ]
  rownames(TCGA_Uni_score_hat) = panTCGA_anno$site_donor_id[!duplicated(panTCGA_anno$site_donor_id)]
  P_score_hat = TCGA_Uni_score_hat[rownames(TCGA_response_drug), ] * norm(Tumor_exp_norm, type = "f")
  
  enet_model = fit_enet(data.frame(C_score_hat), C_drug_AUC)
  yhat.test.PCjoint = predict(enet_model, newx = P_score_hat)
  auc(roc(c(P_drug_bin),c(yhat.test.PCjoint), direction=">"))
  
  ## Summary for each drug
  wilcox_test = c(wilcox.test(yhat.test.EN~P_drug_bin)$p.value,
                  wilcox.test(yhat.test.CombatEN~P_drug_bin)$p.value,
                  wilcox.test(yhat.test.PCjoint~P_drug_bin)$p.value)
  auc = c(auc(roc(c(P_drug_bin), c(yhat.test.EN), direction=">")),
          auc(roc(c(P_drug_bin), c(yhat.test.CombatEN), direction=">")),
          auc(roc(c(P_drug_bin),c(yhat.test.PCjoint), direction=">")))
  print(drug_list[i,1])
  print(best_rank)
  print(wilcox_test)
  print(auc)
  validation = data.frame(wilcox_test = wilcox_test, auc = auc)
  prediction = data.frame(true_response = P_drug_bin, yhat.EN = yhat.test.EN, yhat.CombatEN = yhat.test.CombatEN, yhat.PCjoint = yhat.test.PCjoint)
  rownames(validation) = c("EN", "CombatEN", "PCjoint")
  results_traditon_allDrugs1[[i]] = list(validation = validation, prediction = prediction)
  print(i)
}
names(results_traditon_allDrugs1) = drug_list$TCGA_drug

results_traditon_allDrugs$Vinorelbine$prediction$s0 = results_Enet[[17]]$prediction$s0
results_traditon_allDrugs$Vinorelbine$validation[1,] = results_Enet[[17]]$validation
results_traditon_allDrugs$Oxaliplatin$prediction$s0 = results_Enet[[13]]$prediction$s0
results_traditon_allDrugs$Oxaliplatin$validation[1,] = results_Enet[[13]]$validation
results_traditon_allDrugs$Etoposide$prediction$s0 = results_Enet[[9]]$prediction$s0
results_traditon_allDrugs$Etoposide$validation[1,] = results_Enet[[9]]$validation
results_traditon_allDrugs$Trastuzumab$prediction$s0 = results_Enet[[1]]$prediction$s0
results_traditon_allDrugs$Trastuzumab$validation[1,] = results_Enet[[1]]$validation
results_traditon_allDrugs$Cisplatin$prediction$s0 = results_Enet[[4]]$prediction$s0
results_traditon_allDrugs$Cisplatin$validation[1,] = results_Enet[[4]]$validation

traditon_methods_pred_5drugs = list(Trastuzumab = results_traditon_allDrugs$Trastuzumab$prediction,
                                    Etoposide = results_traditon_allDrugs$Etoposide$prediction,
                                    Vinorelbine = results_traditon_allDrugs$Vinorelbine$prediction,
                                    Cisplatin = results_traditon_allDrugs$Cisplatin$prediction,
                                    Oxaliplatin = results_traditon_allDrugs$Oxaliplatin$prediction)
save(traditon_methods_pred_5drugs, file = "traditon_methods_pred_5drugs.RData")

boxplot(results_traditon_allDrugs$Vinorelbine$prediction$s0.2[which(results_traditon_allDrugs$Vinorelbine$prediction$measure_of_response == "Response")],
        results_traditon_allDrugs$Vinorelbine$prediction$s0.2[which(results_traditon_allDrugs$Vinorelbine$prediction$measure_of_response == "Non-Response")])


results_Enet = list()
irx_best_rank = c()
for (i in c(17)) {
  CCLE_response_drug = read_GDSC_drug_data(Drug_name = drug_list[i,1])
  TCGA_response_drug = read_TCGA_drug_data(Drug_name = drug_list[i,2])
  Cell_Line_Names_Exp = rownames(CCLE_exp_norm)
  #Cell_Line_Names_drug = intersect(rownames(CCLE_exp_norm), rownames(CCLE_response_drug))
  Cell_Line_Names_drug = intersect(rownames(CCLE_response_drug), rownames(CCLE_exp_norm))
  
  C_exp_norm = CCLE_exp_norm[Cell_Line_Names_drug, ]
  TCGA_Names_Exp = rownames(Tumor_exp_norm)
  TCGA_Uni_exp_norm = Tumor_exp_norm[panTCGA_anno$th_sampleid,]
  TCGA_Uni_exp_norm = TCGA_Uni_exp_norm[!duplicated(panTCGA_anno$site_donor_id), ]
  rownames(TCGA_Uni_exp_norm) = panTCGA_anno$site_donor_id[!duplicated(panTCGA_anno$site_donor_id)]
  P_exp_norm = TCGA_Uni_exp_norm[rownames(TCGA_response_drug),]
  
  C_drug_AUC = CCLE_response_drug[Cell_Line_Names_drug, ]%>%select(mean_AUC)
  C_drug_AUC = c(as.matrix(C_drug_AUC))
  
  P_drug_bin = TCGA_response_drug%>%select(measure_of_response)
  P_drug_bin = as.matrix(P_drug_bin)
  for (j in 1:nrow(P_drug_bin)) {
    if (P_drug_bin[j] == "Complete Response" || P_drug_bin[j] == "Partial Response") {
      P_drug_bin[j] = "Response"
    } else {P_drug_bin[j] = "Non-Response"}
  }
  
  ## Elastic Net
  Data = data.frame(cbind(C_drug_AUC, C_exp_norm))
  colnames(Data)[1] = "pheno"
  enet_model = train(pheno ~., data = Data, method = "glmnet",
                     trControl = trainControl("cv", number = 10),
                     tuneLength = 20) 
  tune = enet_model$bestTune
  model_enet = glmnet(C_exp_norm,C_drug_AUC,
                      alpha=tune$alpha,
                      lambda=tune$lambda)
  
  #enet_model = fit_enet(C_exp_norm, C_drug_AUC)
  yhat.test.EN = predict(model_enet, newx = P_exp_norm)
  
  ## Summary for each drug
  wilcox_test = c(wilcox.test(yhat.test.EN~P_drug_bin)$p.value)
  auc = c(auc(roc(c(P_drug_bin), c(yhat.test.EN))))
  print(drug_list[i,1])
  print(wilcox_test)
  print(auc)
  validation = data.frame(wilcox_test = wilcox_test, auc = auc)
  prediction = data.frame(true_response = P_drug_bin, yhat.EN = yhat.test.EN)
  results_Enet[[i]] = list(validation = validation, prediction = prediction)
  print(i)
}

#### JIVE Model ####
for (rankJ in seq(30,80,by = 10)){
 for (rankC in seq(40,300,by = 20)) {
   for (rankP in seq(40,250,by = 20)) {
     assign(paste("jive_panTumor", rankJ, rankC, rankP, sep = "_"), fast.jive(list(C = CCLE_exp_norm, P = Tumor_exp_norm), rankJ = rankJ, rankA = c(rankC, rankP), scale = TRUE, center = FALSE, method = "given", conv = 15, maxiter = 300, est = TRUE, scalefloat = 10000))
     jive_model = get(paste("jive_panTumor", input[2], input[3], input[4], sep = "_"))
     CP_Joint_svd = svd(rbind(jive_model$joint[[1]], jive_model$joint[[2]]))
     
     drug_list = data.frame(GDSC_drug = c("Afatinib", "Bleomycin","Cetuximab","Cisplatin","Cisplatin","Cyclophosphamide","Docetaxel","Doxorubicin","Etoposide","5-Fluorouracil", "Gemcitabine", "Irinotecan", "Oxaliplatin", "Paclitaxel", "Pemetrexed", "Temozolomide", "Vinorelbine"), TCGA_drug = c("Trastuzumab", "Bleomycin", "Cetuximab", "Cisplatin", "Carboplatin", "Cyclophosphamide", "Docetaxel", "Doxorubicin", "Etoposide", "Fluorouracil", "Gemcitabine", "Irinotecan", "Oxaliplatin", "Paclitaxel", "Pemetrexed", "Temozolomide", "Vinorelbine"))
     results_jive_allDrugs = list()
     for (i in 1:17) {
       CCLE_response_drug = read_GDSC_drug_data(Drug_name = drug_list[i,1])
       TCGA_response_drug = read_TCGA_drug_data(Drug_name = drug_list[i,2])
       Cell_Line_Names_Exp = rownames(CCLE_exp_norm)
       Cell_Line_Names_drug = intersect(rownames(CCLE_response_drug), rownames(CCLE_exp_norm))
       C_drug_AUC = CCLE_response_drug[Cell_Line_Names_drug, ]%>%select(mean_AUC)
       C_drug_AUC = as.matrix(C_drug_AUC)
       
       P_drug_bin = TCGA_response_drug%>%select(measure_of_response)
       P_drug_bin = as.matrix(P_drug_bin)
       for (j in 1:nrow(P_drug_bin)) {
         if (P_drug_bin[j] == "Complete Response" || P_drug_bin[j] == "Partial Response") {
           P_drug_bin[j] = "Response"
         } else {P_drug_bin[j] = "Non-Response"}
       }
       #得到training的数据就够了
       #print(table(P_drug_bin))
       #print(nrow(CCLE_response_drug))
       results_jive_allDrugs[[i]] = lm.jive(get(paste("jive_panTumor", input[2], input[3], input[4], sep = "_")), CP_Joint_svd, Cell_Line_Names_Exp, panTCGA_anno, C_drug_AUC, P_drug_bin)
       print(results_jive_allDrugs[[i]]$validation)
     }
     names(results_jive_allDrugs) = drug_list$TCGA_drug
     assign(paste("results_jive_allDrugs", input[2], input[3], input[4], sep = "_"), results_jive_allDrugs)
     
   }
 }
}


library(RSpectra)
library(pROC)
library(caret)
library(glmnet)

jive_CESC_70_140_140_est = 
  fast.jive(list(C = CCLE_exp_norm, P = Tumor_exp_norm), rankJ = 70, rankA = c(140, 140), scale = TRUE, center = TRUE, method = "given", conv = 20, maxiter = 50, est = TRUE, scalefloat = 10000)
CP_Joint_svd = svd(rbind(jive_CESC_70_140_140_est$joint[[1]], jive_CESC_70_140_140_est$joint[[2]]))

drug_list = data.frame(GDSC_drug = c("Afatinib", "Bleomycin","Cetuximab","Cisplatin","Cisplatin","Cyclophosphamide","Docetaxel","Doxorubicin","Etoposide","5-Fluorouracil", "Gemcitabine", "Irinotecan", "Oxaliplatin", "Paclitaxel", "Pemetrexed", "Temozolomide", "Vinorelbine"), TCGA_drug = c("Trastuzumab", "Bleomycin", "Cetuximab", "Cisplatin", "Carboplatin", "Cyclophosphamide", "Docetaxel", "Doxorubicin", "Etoposide", "Fluorouracil", "Gemcitabine", "Irinotecan", "Oxaliplatin", "Paclitaxel", "Pemetrexed", "Temozolomide", "Vinorelbine"))
results_jive_40_140_140_crossDrugs = list()
for (i in 1:17) {
  CCLE_response_drug = read_GDSC_drug_data(Drug_name = drug_list[i,1])
  TCGA_response_drug = read_TCGA_drug_data(Drug_name = drug_list[i,2])
  Cell_Line_Names_Exp = rownames(CCLE_exp_norm)
  Cell_Line_Names_drug = intersect(rownames(CCLE_response_drug), rownames(CCLE_exp_norm))
  C_drug_AUC = CCLE_response_drug[Cell_Line_Names_drug, ]%>%select(mean_AUC)
  C_drug_AUC = as.matrix(C_drug_AUC)
  
  P_drug_bin = TCGA_response_drug%>%select(measure_of_response)
  P_drug_bin = as.matrix(P_drug_bin)
  for (j in 1:nrow(P_drug_bin)) {
    if (P_drug_bin[j] == "Complete Response" || P_drug_bin[j] == "Partial Response") {
      P_drug_bin[j] = "Response"
    } else {P_drug_bin[j] = "Non-Response"}
  }
  print(table(P_drug_bin))
  print(nrow(CCLE_response_drug))
  #results_jive_40_140_140_crossDrugs[[i]] = lm.jive(jive_CESC_70_140_140_est, CP_Joint_svd, Cell_Line_Names_Exp, panTCGA_anno, C_drug_AUC, P_drug_bin)
  print(i)
}
names(results_jive_40_140_140_crossDrugs) = drug_list$TCGA_drug


library(Rcpp)
library(rlang)
library(dplyr)
library(RSpectra)
library(ggplot2)
library(pROC)
library(generics)
library(caret)
library(glmnet)

load(file = "/home/qingzliu/JIVE/jive_array_panTumor_allDrugs/panTumor_CCLE_allDrugs.RData")

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
n <- as.numeric(slurm_arrayid)

input_array = matrix(nrow = 1320, ncol = 4)
i = 1
#for (rankJ in seq(30,80,by = 5)){
#  for (rankC in seq(80,300,by = 20)) {
#    for (rankP in seq(60,250,by = 20)) {
#      input_array[i, ] = c(i, rankJ, rankC, rankP)
#      i = i + 1
#    }
#  }
#}
for (rankJ in seq(30,60,by = 10)){
  for (rankC in c(40,60)) {
    for (rankP in c(40)) {
      input_array[i, ] = c(i, rankJ, rankC, rankP)
      i = i + 1
    }
  }
}
input = input_array[n, ]

Sys.time()
assign(paste("jive_panTumor", input[2], input[3], input[4], sep = "_"), fast.jive(list(C = CCLE_exp_norm, P = Tumor_exp_norm), rankJ = input[2], rankA = c(input[3], input[4]), scale = TRUE, center = FALSE, method = "given", conv = 15, maxiter = 300, est = TRUE, scalefloat = 10000))
#save(list = paste("jive_CESC", input[2], input[3], input[4], sep = "_"), file = paste(paste("../jive_array_panTumor_CCLE_align/jive_panTumor", input[2], input[3], input[4], sep = "_"), ".RData", sep = ""))

jive_model = get(paste("jive_panTumor", input[2], input[3], input[4], sep = "_"))
CP_Joint_svd = svd(rbind(jive_model$joint[[1]], jive_model$joint[[2]]))

drug_list = data.frame(GDSC_drug = c("Afatinib", "Bleomycin","Cetuximab","Cisplatin","Cisplatin","Cyclophosphamide","Docetaxel","Doxorubicin","Etoposide","5-Fluorouracil", "Gemcitabine", "Irinotecan", "Oxaliplatin", "Paclitaxel", "Pemetrexed", "Temozolomide", "Vinorelbine"), TCGA_drug = c("Trastuzumab", "Bleomycin", "Cetuximab", "Cisplatin", "Carboplatin", "Cyclophosphamide", "Docetaxel", "Doxorubicin", "Etoposide", "Fluorouracil", "Gemcitabine", "Irinotecan", "Oxaliplatin", "Paclitaxel", "Pemetrexed", "Temozolomide", "Vinorelbine"))
results_jive_allDrugs = list()
for (i in 1:17) {
  CCLE_response_drug = read_GDSC_drug_data(Drug_name = drug_list[i,1])
  TCGA_response_drug = read_TCGA_drug_data(Drug_name = drug_list[i,2])
  Cell_Line_Names_Exp = rownames(CCLE_exp_norm)
  Cell_Line_Names_drug = intersect(rownames(CCLE_response_drug), rownames(CCLE_exp_norm))
  C_drug_AUC = CCLE_response_drug[Cell_Line_Names_drug, ]%>%select(mean_AUC)
  C_drug_AUC = as.matrix(C_drug_AUC)
  
  P_drug_bin = TCGA_response_drug%>%select(measure_of_response)
  P_drug_bin = as.matrix(P_drug_bin)
  for (j in 1:nrow(P_drug_bin)) {
    if (P_drug_bin[j] == "Complete Response" || P_drug_bin[j] == "Partial Response") {
      P_drug_bin[j] = "Response"
    } else {P_drug_bin[j] = "Non-Response"}
  }
  #print(table(P_drug_bin))
  #print(nrow(CCLE_response_drug))
  results_jive_allDrugs[[i]] = lm.jive(get(paste("jive_panTumor", input[2], input[3], input[4], sep = "_")), CP_Joint_svd, Cell_Line_Names_Exp, panTCGA_anno, C_drug_AUC, P_drug_bin)
  print(results_jive_allDrugs[[i]]$validation)
}
names(results_jive_allDrugs) = drug_list$TCGA_drug
assign(paste("results_jive_allDrugs", input[2], input[3], input[4], sep = "_"), results_jive_allDrugs)
Sys.time()

save(list = paste("results_jive_allDrugs", input[2], input[3], input[4], sep = "_"), file = paste(paste("../jive_array_panTumor_allDrugs_Out/JIVE_panTumor_allDrugs", input[2], input[3], input[4], sep = "_"), ".RData", sep = ""))

####

function(jive_object, CP_Joint_svd, Cell_Line_Names_Exp, panTCGA_anno, Cell_Line_Drug_Response, TCGA_Drug_Response){
  library(pROC)
  library(caret)
  library(glmnet)
  Cell_Line_Names_Drug = rownames(Cell_Line_Drug_Response)
  TCGA_Names_Drug = rownames(TCGA_Drug_Response)
  #CP_Joint_svd = svd(rbind(jive_object$joint[[1]], jive_object$joint[[2]]))
  CP_score_hat = CP_Joint_svd$u[,1:jive_object$rankJ]%*%diag(CP_Joint_svd$d[1:jive_object$rankJ])
  CCLE_score_hat = CP_score_hat[1:nrow(jive_object$joint[[1]]), ]
  rownames(CCLE_score_hat) = Cell_Line_Names_Exp
  C_score_hat = CCLE_score_hat[Cell_Line_Names_Drug, ]
  Tumor_score_hat = CP_score_hat[(nrow(jive_object$joint[[1]])+1):nrow(CP_score_hat), ]
  rownames(Tumor_score_hat) = Tumor_anno$th_sampleid
  TCGA_score_hat = Tumor_score_hat[panTCGA_anno$th_sampleid,]
  rownames(TCGA_score_hat) = panTCGA_anno$th_sampleid
  TCGA_Uni_score_hat = TCGA_score_hat[!duplicated(panTCGA_anno$site_donor_id), ]
  rownames(TCGA_Uni_score_hat) = panTCGA_anno$site_donor_id[!duplicated(panTCGA_anno$site_donor_id)]
  P_score_hat = TCGA_Uni_score_hat[rownames(TCGA_Drug_Response), ]
  rank_summary = data.frame(Name = c("Joint Rank", "TCGA Indv Rank", "CCLE Indv Rank", "TCGA Dim", "CCLE Dim"), Number = c(jive_object$rankJ, jive_object$rankA[2], jive_object$rankA[1], paste(dim(jive_object$data$P)[1], dim(jive_object$data$P)[2], sep = "x"), paste(dim(jive_object$data$C)[1], dim(jive_object$data$C)[2], sep = "x")))
  
  # linear regression
  DC_model = lm(Cell_Line_Drug_Response~C_score_hat)
  beta_hat = summary(DC_model)$coefficient[,1]
  D_p_hat_lr = cbind(rep(1,nrow(P_score_hat)), P_score_hat)%*%beta_hat
  # ridge regression
  cv_ridge = cv.glmnet(C_score_hat, Cell_Line_Drug_Response, alpha = 0)
  best_lambda = cv_ridge$lambda.min
  best_ridge = glmnet(C_score_hat, Cell_Line_Drug_Response, alpha = 0, lambda = best_lambda)
  D_p_hat_ridge = predict(best_ridge, s = best_lambda, newx = P_score_hat)
  # lasso regression
  cv_lasso = cv.glmnet(C_score_hat, Cell_Line_Drug_Response, alpha = 1)
  best_lambda = cv_lasso$lambda.min
  best_lasso = glmnet(C_score_hat, Cell_Line_Drug_Response, alpha = 1, lambda = best_lambda)
  D_p_hat_lasso = predict(best_lasso, s = best_lambda, newx = P_score_hat)
  # elastic net
  Data = as.data.frame(cbind(Cell_Line_Drug_Response, C_score_hat))
  enet_model = train(mean_AUC ~., data = Data, method = "glmnet",
                     trControl = trainControl("repeatedcv", number = 10, repeats = 3),
                     tuneLength = 60)
  tune = enet_model$bestTune
  best_enet = glmnet(as.matrix(C_score_hat),Cell_Line_Drug_Response,
                     alpha=tune$alpha,
                     lambda=tune$lambda)
  D_p_hat_enet = predict(best_enet,newx = P_score_hat)
  # xgboost
  
  # output
  predicted_values = data.frame(lr = c(D_p_hat_lr), ridge = c(D_p_hat_ridge), lasso = c(D_p_hat_lasso), enet = c(D_p_hat_enet))
  wilcox_test = c(wilcox.test(D_p_hat_lr~TCGA_Drug_Response)$p.value,
                  wilcox.test(D_p_hat_ridge~TCGA_Drug_Response)$p.value,
                  wilcox.test(D_p_hat_lasso~TCGA_Drug_Response)$p.value,
                  wilcox.test(D_p_hat_enet~TCGA_Drug_Response)$p.value)
  auc = c(auc(roc(c(TCGA_Drug_Response), c(D_p_hat_lr))),
          auc(roc(c(TCGA_Drug_Response), c(D_p_hat_ridge))),
          auc(roc(c(TCGA_Drug_Response), c(D_p_hat_lasso))),
          auc(roc(c(TCGA_Drug_Response), c(D_p_hat_enet))))
  validation = data.frame(wilcox_test, auc)
  rownames(validation) = c("lr", "ridge", "lasso", "enet")
  best_pred = validation[which(auc == max(auc)), ]
  #box_plot_data = data.frame(D_p_hat_enet, TCGA_Drug_Response)
  #colnames(box_plot_data) = c("Predicted_lnIC50", "True_Response")
  #Prediction_Result_Plot = ggplot(box_plot_data, aes(x=True_Response, y=Predicted_lnIC50)) + geom_boxplot()
  return(list(predicted_values = predicted_values, validation = validation, best_pred = best_pred, C_train = list(C_score_hat, Cell_Line_Drug_Response), rank_summary = rank_summary))
}






