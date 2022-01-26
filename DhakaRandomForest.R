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
split1 = sort(sample(nrow(complete_dat), nrow(complete_dat) * .7))
train <- complete_dat[split1,]
test <- complete_dat[-split1,]
out=ranger(as.formula(paste(resp_var,'~',paste(names,collapse="+"),sep="")),data=complete_dat,num.trees=1000,importance="impurity")
df_imps=data.frame(names=names(ranger::importance(out)),imps=ranger::importance(out)) %>% arrange(desc(imps))
#creating test and train datasets
