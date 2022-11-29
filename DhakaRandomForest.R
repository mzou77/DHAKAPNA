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
library(party)
library(permimp)

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
  out=ranger(as.formula(paste(resp_var,'~',paste(names,collapse="+"),sep="")),data=train,num.trees=5000,importance="impurity")
  df_imps=data.frame(var=names(out$variable.importance),imp=as.numeric(out$variable.importance)) %>% arrange(desc(imp))
  out=ranger(as.formula(paste(resp_var,'~',paste(df_imps$var[1:2],collapse="+"),sep="")),data=train,num.trees=5000,importance="impurity")
  pred=predict(out,data=test)$predictions
  return(pred)
}                              
cv_preds=train %>% purrr::map(function(x) cv(x))
aucs=map2(cv_preds,train, function(x,y) as.numeric(roc(complete_dat[-y,]$treat_fail,x)$auc))                            
mean(unlist(aucs))
#creating test and train datasets
#repeated five fold cross validation

finish_imputed_dat=read.csv('finish_imputed_dat.csv') 
complete_dat=read.csv('complete_dat.csv')
names = colnames(complete_dat)[c(4, 8:12, 14:15, 17:37, 39:63, 65:123, 126:127, 129:135, 138:142)]
cdat=complete_dat %>% select(treat_fail,one_of(names))
cdat = cdat %>% mutate(treat_fail=as.factor(treat_fail))
rf=cforest(treat_fail~ .,data=data.frame(cdat),controls=cforest_unbiased(ntree=5000,mtry=floor(dim(cdat)[2]-1)/3))
rf_imp=permimp(rf,conditional=T,AUC=T) 
df_imps=data.frame(var=names(rf_imp$values),imps=as.numeric(rf_imp$values)) %>% arrange(desc(imps));df_imps

set.seed(123)
nrows_dat=nrow(cdat)

#train=sample(1:nrows_dat,round(.80*nrows_dat))
niter=100
train = 1:niter %>% purrr::map(function(x) sample(1:nrows_dat,round(.80*nrows_dat)))
#x=train[[1]]
cv = function(x,nvars){ 
  train=cdat[x,]
  test=cdat[-x,] 
  rf=cforest(treat_fail~ .,data=train, controls=cforest_unbiased(ntree=2000, mtry=floor((dim(train)[2]-1)/3)))
  rf_imp=permimp(rf,conditional=T,AUC=T)
  df_imps=data.frame(var=names(rf_imp$values), imps=as.numeric(rf_imp$values)) %>% arrange(desc(imps))
  rf=cforest(treat_fail~., data=train %>% select(treat_fail,one_of(df_imps$var[1:nvars])), controls=cforest_unbiased(ntree=500, mtry=floor(nvars/3)))
  pred=matrix(unlist(predict(rf,newdata=test,type="prob")),ncol=2,byrow=T)[,2]
  return(pred)
}      
cv_preds=train %>% purrr::map(function(x) cv(x,9))
aucs=map2(cv_preds,train, function(x,y) as.numeric(roc(cdat[-y,]$treat_fail,x)$auc))                            
mean(unlist(aucs))

#gives AUC as output for 5 fold cross validation

rf=cforest(treat_fail~ .,data=cdat, controls=cforest_unbiased(ntree=500, mtry=floor((dim(cdat)[2]-1)/3)))
rf_imp=permimp(rf,conditional=T,AUC=T)
df_imps=data.frame(var=names(rf_imp$values), imps=as.numeric(rf_imp$values)) %>% arrange(desc(imps))
df_imps
#applies model to entire dataset and extracts list of most important variables  


cdat["RSys_grunting"][cdat["RSys_grunting"] == 2] <- 0 #recoding 0 as no, 1 as yes for categorical variables
cdat["Central_Cyan_enrol"][cdat["Central_Cyan_enrol"] == 2] <- 0
cdat["sever_res_enrol"][cdat["sever_res_enrol"] == 2] <- 0
cdat["RSys_lcwi"][cdat["RSys_lcwi"] == 2] <- 0
logistic_model <- glm(treat_fail ~ RSys_grunting + SPO2_withoutO2 + Rectal_tem + RSys_lcwi + sever_res_enrol + Central_Cyan_enrol + WHZ + WAZ + Enroll_dur_poorfeed + Enroll_MUAC, data=cdat, family = "binomial")
summary(logistic_model)
exp(logistic_model$coefficients[-1])
confint(logistic_model, level = .95)
#logistic regression to calculate odd's ratios

pdp1=rf %>% pdp::partial(pred.var="RSys_grunting")
ggplot(pdp1) + geom_point(aes(x=RSys_grunting,y=1-yhat)) + ylab("Marginal Predicted Probability") + xlab("Presence of grunting on exam") + xlim("1","2")

#or y=yhat 
  
pdp2=rf %>% pdp::partial(pred.var="SPO2_withoutO2")
ggplot(pdp2,aes(x=SPO2_withoutO2,y=1-yhat)) + geom_point() + geom_line() + ylab("Marginal Predicted Probability") + xlab("Room air saturation")

pdp3=rf %>% pdp::partial(pred.var="Central_Cyan_enrol")
ggplot(pdp3) + geom_point(aes(x=Central_Cyan_enrol,y=1-yhat)) + ylab("Marginal Predicted Probability") + xlab ("Presence of central cyanosis on exam") + xlim("1","2")

pdp4=rf %>% pdp::partial(pred.var="Rectal_tem")
ggplot(pdp4,aes(x=Rectal_tem,y=1-yhat)) + geom_point() + geom_line() + ylab("Marginal Predicted Probability") + xlab("Rectal Temperature in Celsius")

pdp5=rf %>% pdp::partial(pred.var="Age_father")
ggplot(pdp5,aes(x=Age_father,y=1-yhat)) + geom_point() + geom_line()

pdp6=rf %>% pdp::partial(pred.var="WAZ")
ggplot(pdp6,aes(x=WAZ,y=1-yhat)) + geom_point() + geom_line() + ylab("Marginal Predicted Probability") + xlab("WAZ (SD)")

pdp7=rf %>% pdp::partial(pred.var="sever_res_enrol")
ggplot(pdp7) + geom_point(aes(x=sever_res_enrol,y=1-yhat)) + ylab("Marginal Predicted Probability") +  xlab("Respiratory distress at time of enrollment") + xlim("1", "2")

pdp8=rf %>% pdp::partial(pred.var="Enroll_dur_poorfeed")
ggplot(pdp8,aes(x=Enroll_dur_poorfeed,y=1-yhat)) + geom_point() + geom_line() + ylab("Marginal Predicted Probability") + xlab("Duration of Poor Feeding (Days)") + xlim(0,90)

pdp9=rf %>% pdp::partial(pred.var="RSys_lcwi")
ggplot(pdp9) + geom_point(aes(x=RSys_lcwi,y=1-yhat)) + ylab("Marginal Predicted Probability") + xlab("Lower chest wall indrawing on exam") + xlim("1", "2")

pdp10=rf %>% pdp::partial(pred.var="Antibiotic_received7days")
ggplot(pdp10) + geom_point(aes(x=Antibiotic_received7days,y=1-yhat)) 

pdp11=rf %>% pdp::partial(pred.var="WHZ")
ggplot(pdp11,aes(x=WHZ,y=1-yhat)) + geom_point() + geom_line() + ylab("Marginal Predicted Probability") + xlab("WHZ")

pdp12=rf %>% pdp::partial(pred.var="Enroll_MUAC")
ggplot(pdp12,aes(x=Enroll_MUAC,y=1-yhat)) + geom_point() + geom_line() + ylab("Marginal Predicted Probability") + xlab("MUAC (cm)")

#partial dependency plots of top 10 most important predictors
