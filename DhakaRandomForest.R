library(tidyverse)
library(lubridate)
library(pROC)
library(fields)
library(zoo)
library(ks)
library(KernSmooth)
library(ranger)
library(broom)
library(furrr)
library(cvAUC)
library(pdp)
library(slider)
library(ranger)
library(caret)

complete_dat <- merge(finish_imputed_dat, Amoxicillin_Treatment_Failure_Sorted_27_10_2021_, by = "Study_ID")
# completes dataset by merging cleaned enrollment characteristics dataframe with treatment failure dataframe

names = colnames(complete_dat)[c(2:3, 7:11, 13:36, 38:62, 64:126, 128:134, 137:141)]
# specifices columns of variables included in analysis

resp_var = "treat_fail"
# defines response variable

out = ranger(as.formula(paste(resp_var,'~',paste(names,collapse="+"),sep="")),data=complete_dat,num.trees=2000,importance="impurity")
imps=importance(out)
df_imps=data.frame(names=names(imps),var_red=as.numeric(imps)) %>% arrange(desc(var_red))
head(df_imps,10)
# calculates overall importance using random forest

set.seed(123)
nrows_dat=nrow(complete_dat)
train=sample(1:nrows_dat,round(.80*nrows_dat))
train = 1:10 %>% purrr::map(function(x) sample(1:nrows_dat,round(.80*nrows_dat)))
cv = function(x){ 
  train=complete_dat[x,]
  test=complete_dat[-x,] 
  out=ranger(as.formula(paste(resp_var,'~',paste(names,collapse="+"),sep="")),data=train,num.trees=1000,importance="impurity")
  df_imps=data.frame(var=names(out$variable.importance),imp=as.numeric(out$variable.importance)) %>% arrange(desc(imp))
  out=ranger(as.formula(paste(resp_var,'~',paste(df_imps$var[1:5],collapse="+"),sep="")),data=train,num.trees=1000,importance="impurity")
  pred=predict(out,data=test)$predictions
  return(pred)
}                              
cv_preds=train %>% purrr::map(function(x) cv(x))
aucs=map2(cv_preds,train, function(x,y) as.numeric(roc(complete_dat[-y,]$treat_fail,x)$auc))                            
mean(unlist(aucs))
#creating test and train datasets
#repeated five fold cross validation
