
## Random forest and Cox Model : Survival analysis

library(SurvMetrics)
library(caret)
library(randomForestSRC)
library(survival)  
library(pec)
library(ggplot2)
set.seed(123)

#simulation

#Initialization
metrics_cox = 0
metrics_rsf = 0
for (i in 1:20) {
  
  mydata = SDGM4(N = 100, p = 20, c_step = 0.2)
  index_data = createFolds(1:nrow(mydata), 2)
  train_data = mydata[index_data[[1]],]
  test_data = mydata[index_data[[2]],]

  #fit the models
  fitrsf = rfsrc(Surv(time, status)~., data = train_data, nsplit = 3, ntree = 500)
  mat_rsf = predict(fitrsf, test_data)$survival
  
  dis_time = fitrsf$time.interest
  fitcox = coxph(Surv(time, status)~., data = train_data, x = TRUE)
  mat_cox = predictSurvProb(fitcox, test_data, dis_time)
  
  #calculate the C index
  med_index = median(1:length(dis_time))
  surv_obj = Surv(test_data$time, test_data$status)
  
  #C index for Cox
  metrics_cox[i] = Cindex(surv_obj, predicted = mat_cox[, med_index])
  #C index for RSF
  metrics_rsf[i] = Cindex(surv_obj, predicted = mat_rsf[, med_index])
}
data_CI = data.frame('Cindex' = c(metrics_cox, metrics_rsf),
                     'model' = c(rep('Cox', 20), rep('RSF', 20)))

ggplot(data_CI, aes(x = model, y = Cindex, fill = model)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFBBCC", "#88CCFF"))


## ----fig.align='center', warning=FALSE, results='hide'------------------------
#Initialization

metrics_cox = 0
metrics_rsf = 0
set.seed(123)
for (i in 1:10) {
  
  mydata = SDGM4(N = 100, p = 20, c_step = 0.5)
  index_data = createFolds(1:nrow(mydata), 2)
  train_data = mydata[index_data[[1]],]
  test_data = mydata[index_data[[2]],]
  
  #fit the models
  fitrsf = rfsrc(Surv(time, status)~., data = train_data, nsplit = 3, ntree = 500)
  mat_rsf = predict(fitrsf, test_data)$survival
  
  dis_time = fitrsf$time.interest
  fitcox = coxph(Surv(time, status)~., data = train_data, x = TRUE)
  mat_cox = predictSurvProb(fitcox, test_data, dis_time)
  
  #calculate the Brier Score
  med_index = median(1:length(dis_time))
  surv_obj = Surv(test_data$time, test_data$status)
  t_star = median(fitrsf$time.interest)
  
  #Brier Score for Cox
  metrics_cox[i] = Brier(surv_obj, pre_sp = mat_cox[, med_index], t_star)
  #Brier Score for RSF
  metrics_rsf[i] = Brier(surv_obj, pre_sp = mat_rsf[, med_index], t_star)
  
}

data_BS = data.frame('BS' = c(metrics_cox, metrics_rsf),
                     'model' = c(rep('Cox', 10), rep('RSF', 10)))

ggplot(data_BS, aes(x = model, y = BS, fill = model)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFBBCC", "#88CCFF"))

## ----fig.align='center', message=FALSE, warning=FALSE, results='hide'---------
#Initialization

metrics_cox = 0
metrics_rsf = 0
set.seed(123)
for (i in 1:5) {
  
  mydata = SDGM4(N = 100, p = 20, c_step = -0.5)
  index_data = createFolds(1:nrow(mydata), 2)
  train_data = mydata[index_data[[1]],]
  test_data = mydata[index_data[[2]],]
  
  #fit the models
  fitrsf = rfsrc(Surv(time, status)~., data = train_data, nsplit = 3, ntree = 500)
  mat_rsf = predict(fitrsf, test_data)$survival
  
  dis_time = fitrsf$time.interest
  fitcox = coxph(Surv(time, status)~., data = train_data, x = TRUE)
  mat_cox = predictSurvProb(fitcox, test_data, dis_time)
  
  #calculate the IBS
  med_index = median(1:length(dis_time))
  surv_obj = Surv(test_data$time, test_data$status)
  
  
  #IBS for Cox
  metrics_cox[i] = IBS(surv_obj, sp_matrix = mat_cox, dis_time)
  #IBS for RSF
  metrics_rsf[i] = IBS(surv_obj, sp_matrix = mat_rsf, dis_time)
  
}

data_IBS = data.frame('IBS' = c(metrics_cox, metrics_rsf),
                      'model' = c(rep('Cox', 5), rep('RSF', 5)))

ggplot(data_IBS, aes(x = model, y = IBS, fill = model)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFBBCC", "#88CCFF", "#FF0000"))

## ----fig.align='center', message=FALSE, warning=FALSE, results='hide'---------
#Initialization

metrics_cox_IAE = 0
metrics_rsf_IAE = 0
metrics_cox_ISE = 0
metrics_rsf_ISE = 0
set.seed(123)
for (i in 1:10) {
  
  mydata = SDGM4(N = 100, p = 20, c_step = 0.2)
  index_data = createFolds(1:nrow(mydata), 2)
  train_data = mydata[index_data[[1]],]
  test_data = mydata[index_data[[2]],]
  
  #fit the models
  fitrsf = rfsrc(Surv(time, status)~., data = train_data, nsplit = 3, ntree = 500)
  mat_rsf = predict(fitrsf, test_data)$survival
  
  dis_time = fitrsf$time.interest
  fitcox = coxph(Surv(time, status)~., data = train_data, x = TRUE)
  mat_cox = predictSurvProb(fitcox, test_data, dis_time)
  
  #calculate the IAE and ISE
  med_index = median(1:length(dis_time))
  surv_obj = Surv(test_data$time, test_data$status)
  
  
  #IAE and ISE for Cox
  temp1 = IAEISE(surv_obj, sp_matrix = mat_cox, dis_time)
  metrics_cox_IAE[i] = temp1[1]
  metrics_cox_ISE[i] = temp1[2]
  #IAE and ISE for RSF
  temp2 = IAEISE(surv_obj, sp_matrix = mat_rsf, dis_time)
  metrics_rsf_IAE[i] = temp2[1]
  metrics_rsf_ISE[i] = temp2[2]
  
}

data_IAE = data.frame('IAE' = c(metrics_cox_IAE, metrics_rsf_IAE),
                      'model' = c(rep('Cox', 10), rep('RSF', 10)))

data_ISE = data.frame('ISE' = c(metrics_cox_ISE, metrics_rsf_ISE),
                      'model' = c(rep('Cox', 10), rep('RSF', 10)))

P1 = ggplot(data_IAE, aes(x = model, y = IAE, fill = model)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFBBCC", "#88CCFF")) #+
  #theme(legend.position = 'none')

P2 = ggplot(data_ISE, aes(x = model, y = ISE, fill = model)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFBBCC", "#88CCFF")) #+
  #theme(legend.position = 'none')

library(ggpubr)
ggarrange(P1, P2, ncol = 2)


#1. data preparation
library(survival) # to fit a Cox model
library(randomForestSRC) # to fit an RSF model
library(SurvMetrics) # to get all the metrics
library(pec) # to make predictions based on Cox model

set.seed(1)

mydata <- kidney[, -1]
train_index <- sample(1:nrow(mydata), 0.7 * nrow(mydata))
train_data <- mydata[train_index, ]
test_data <- mydata[-train_index, ]

#2. fit the RSF model and Cox model to predict the testing set
#2.1 RSF model

fit_rsf <- rfsrc(Surv(time,status)~., data = train_data) #fit the RSF model
distime <- fit_rsf$time.interest #get the survival time of events
med_index <- median(1:length(distime)) #the index of median survival time of events
mat_rsf <- predict(fit_rsf, test_data)$survival #get the survival probability matrix
vec_rsf <- mat_rsf[ ,med_index] #median survival probability of all samples

#2.2 Cox model

fit_cox <- coxph(Surv(time,status)~., data = train_data, x = TRUE) #fit the Cox model
mat_cox <- predictSurvProb(fit_cox, test_data, distime) #get the survival probability matrix
vec_cox <- mat_cox[ ,med_index]

#MAE(coxph(Surv(time, status), test_data$status, distime[med_index])
#MAE(fit_rsf, test_data, distime[med_index])

#calculate the C index
med_index = median(1:length(distime))
surv_obj = Surv(test_data$time, test_data$status)

#3. get all the metrics by SurvMetrics
#3.1 CI BS IBS IAE ISE based on RSF model: standard model input methods

Cindex_rsf <- Cindex(fit_rsf, test_data)
BS_rsf <- Brier(fit_rsf, test_data, distime[med_index])
IBS_rsf <- IBS(fit_rsf, test_data)
IAE_rsf <- IAEISE(fit_rsf, test_data)[1]
ISE_rsf <- IAEISE(fit_rsf, test_data)[2]

#MAE(Surv(time, status), pre_time)

#CI BS IBS IAE ISE based on Cox model: standard model input methods

Cindex_cox <- Cindex(fit_cox, test_data)
BS_cox <- Brier(fit_cox, test_data, distime[med_index])
IBS_cox <- IBS(fit_cox, test_data)
IAE_cox <- IAEISE(fit_cox, test_data)[1]
ISE_cox <- IAEISE(fit_cox, test_data)[2]

#MAE(Surv(time, status), pre_time)


##############################
#3.2 CI BS IBS IAE ISE based on RSF model: Non-standard model input methods

times <- test_data$time
status <- test_data$status
Cindex_rsf <- Cindex(Surv(times, status), vec_rsf)
BS_rsf <- Brier(Surv(times, status), vec_rsf, distime[med_index])
IBS_rsf <- IBS(Surv(times, status), mat_rsf, distime) # distime can be replaced by range(distime)
IAE_rsf <- IAEISE(Surv(times, status), mat_rsf, distime)[1]
ISE_rsf <- IAEISE(Surv(times, status), mat_rsf, distime)[2]

#MAE(Surv(time, status), pre_time)

#CI BS IBS IAE ISE based on Cox model: Non-standard model input methods

Cindex_cox <- Cindex(Surv(times, status), vec_cox)
BS_cox <- Brier(Surv(times, status), vec_cox, distime[med_index])
IBS_cox <- IBS(Surv(times, status), mat_cox, distime)
IAE_cox <- IAEISE(Surv(times, status), mat_cox, distime)[1]
ISE_cox <- IAEISE(Surv(times, status), mat_cox, distime)[2]

#MAE(Surv(time, status), pre_time)

#CI BS IBS IAE ISE based on RST and Cox model: Non-standard model input methods

Cindex_rsf_cox <- Cindex(Surv(times, status), vec_rsf_cox)
BS_rsf_cox <- Brier(Surv(times, status), vec_rsf_cox, distime1[med_index1])
IBS_rsf_cox <- IBS(Surv(times, status), mat_rsf_cox, distime1)
IAE_rsf_cox <- IAEISE(Surv(times, status), mat_rsf_cox, distime1)[1]
ISE_rsf_cox <- IAEISE(Surv(times, status), mat_rsf_cox, distime1)[2]


