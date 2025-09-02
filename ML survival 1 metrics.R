
## Random forest and Cox Model : Survival analysis

library(SurvMetrics)
library(caret)
library(randomForestSRC)
library(survival)  
library(pec)
library(ggplot2)

#1. data preparation

set.seed(1)

#mydata <- kidney[, -1]

mydata <-tbsurv
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

#3. get all the metrics by SurvMetrics
#3.1 CI BS IBS IAE ISE based on RSF model: standard model input methods

Cindex_rsf <- Cindex(fit_rsf, test_data)
BS_rsf <- Brier(fit_rsf, test_data, distime[med_index])
IBS_rsf <- IBS(fit_rsf, test_data)
IAE_rsf <- IAEISE(fit_rsf, test_data)[1]
ISE_rsf <- IAEISE(fit_rsf, test_data)[2]

c(Cindex_rsf, BS_rsf, IBS_rsf, IAE_rsf, ISE_rsf)

#CI BS IBS IAE ISE based on Cox model: standard model input methods

Cindex_cox <- Cindex(fit_cox, test_data)
BS_cox <- Brier(fit_cox, test_data, distime[med_index])
IBS_cox <- IBS(fit_cox, test_data)
IAE_cox <- IAEISE(fit_cox, test_data)[1]
ISE_cox <- IAEISE(fit_cox, test_data)[2]

c(Cindex_cox, BS_cox, IBS_cox, IAE_cox, ISE_cox)

#plot
dataCindex = data.frame('Cindex' = c(Cindex_cox, Cindex_rsf),
                        'model' = c(rep('Cox', 5), rep('RSF', 5)))

dataBS = data.frame('BS' = c(BS_cox, BS_rsf),
                    'model' = c(rep('Cox', 5), rep('RSF', 5)))

dataIBS = data.frame('IBS' = c(IBS_cox, IBS_rsf),
                     'model' = c(rep('Cox', 5), rep('RSF', 5)))

dataIAE = data.frame('IAE' = c(IAE_cox, IAE_rsf),
                     'model' = c(rep('Cox', 5), rep('RSF', 5)))

dataISE = data.frame('ISE' = c(ISE_cox, ISE_rsf),
                     'model' = c(rep('Cox', 5), rep('RSF', 5)))

P1 = ggplot(dataCindex, aes(x = model, y = Cindex, fill = model)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFBBCC", "#88CCFF"))

P2 = ggplot(dataBS, aes(x = model, y = BS, fill = model)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFBBCC", "#88CCFF")) 

P3 = ggplot(dataIBS, aes(x = model, y = IBS, fill = model)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFBBCC", "#88CCFF")) 

P4 = ggplot(dataIAE, aes(x = model, y = IAE, fill = model)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFBBCC", "#88CCFF")) 

P5 = ggplot(dataISE, aes(x = model, y = ISE, fill = model)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFBBCC", "#88CCFF")) 

library(ggpubr)
ggarrange(P1, P2, P3, P4, P5)

#3.2 CI BS IBS IAE ISE based on RSF model: Non-standard model input methods

times <- test_data$time
status <- test_data$status
Cindex_rsf <- Cindex(Surv(times, status), vec_rsf)
BS_rsf <- Brier(Surv(times, status), vec_rsf, distime[med_index])
IBS_rsf <- IBS(Surv(times, status), mat_rsf, distime) # distime can be replaced by range(distime)
IAE_rsf <- IAEISE(Surv(times, status), mat_rsf, distime)[1]
ISE_rsf <- IAEISE(Surv(times, status), mat_rsf, distime)[2]

c(Cindex_rsf, BS_rsf, IBS_rsf, IAE_rsf, ISE_rsf)

#CI BS IBS IAE ISE based on Cox model: Non-standard model input methods
Cindex_cox <- Cindex(Surv(times, status), vec_cox)
BS_cox <- Brier(Surv(times, status), vec_cox, distime[med_index])
IBS_cox <- IBS(Surv(times, status), mat_cox, distime)
IAE_cox <- IAEISE(Surv(times, status), mat_cox, distime)[1]
ISE_cox <- IAEISE(Surv(times, status), mat_cox, distime)[2]

c(Cindex_cox, BS_cox, IBS_cox, IAE_cox, ISE_cox)

#plot
dataCindex = data.frame('Cindex' = c(Cindex_cox, Cindex_rsf),
                        'model' = c(rep('Cox', 5), rep('RSF', 5)))

dataBS = data.frame('BS' = c(BS_cox, BS_rsf),
                    'model' = c(rep('Cox', 5), rep('RSF', 5)))

dataIBS = data.frame('IBS' = c(IBS_cox, IBS_rsf),
                     'model' = c(rep('Cox', 5), rep('RSF', 5)))

dataIAE = data.frame('IAE' = c(IAE_cox, IAE_rsf),
                     'model' = c(rep('Cox', 5), rep('RSF', 5)))

dataISE = data.frame('ISE' = c(ISE_cox, ISE_rsf),
                     'model' = c(rep('Cox', 5), rep('RSF', 5)))

P1 = ggplot(dataCindex, aes(x = model, y = Cindex, fill = model)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFBBCC", "#88CCFF"))

P2 = ggplot(dataBS, aes(x = model, y = BS, fill = model)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFBBCC", "#88CCFF")) 

P3 = ggplot(dataIBS, aes(x = model, y = IBS, fill = model)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFBBCC", "#88CCFF")) 

P4 = ggplot(dataIAE, aes(x = model, y = IAE, fill = model)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFBBCC", "#88CCFF")) 

P5 = ggplot(dataISE, aes(x = model, y = ISE, fill = model)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFBBCC", "#88CCFF")) 

library(ggpubr)
ggarrange(P1, P2, P3, P4, P5)
