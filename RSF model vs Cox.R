
--
Title: "Machine Learning Prediction of Tuberculosis Mortality: A Comparative Analysis of Random Survival Forest and Cox Regression Models"
author: "Dr. Azeez Adeboye"
date: "4/08/2025"
output: "html_document"
---

library("knitr")
opts_chunk$set(#fig.path = "regression/",
  fig.align = "center",
  fig.pos = "!htb",
  fig.show = "hold",
  # fig.height = 3,
  # fig.width = 4,
  size = "footnotesize",
  prompt = TRUE,
  highlight = FALSE,
  comment = NA,
  
  # Change echo to TRUE if you want to see all the code examples
  echo = FALSE, 
  
  results = FALSE,
  message = FALSE,
  warning = FALSE,
  error = FALSE)

# Setup the R environment
options(object.size = Inf, expressions = 100000, memory = Inf, 
        replace.assign = TRUE, width = 90)
options(mc.cores = 1, rf.cores = -1)


library("ggplot2")         # Graphics engine
library("RColorBrewer")    # Nice color palettes
library("plot3D")          # for 3d surfaces. 
library("dplyr")           # Better data manipulations
library("tidyr")           # gather variables into long format
library("parallel")        # mclapply for multicore processing
library("randomForestSRC") # random forests for survival, regression and 
# classification
library("ggRandomForests") # ggplot2 random forest figures (This!)

################ Default Settings ##################
theme_set(theme_bw())     # A ggplot2 theme with white background

set.seed(1234)
## Set open circle for censored, and x for events 
event_marks <- c(1, 13)
event_labels <- c(FALSE, TRUE)

## We want red for death events, so reorder this set.
str_col <- brewer.pal(3, "Set1")[c(2, 1, 3)]

# Set modes correctly. For binary variables: transform to logical
tbsurv$status <- as.logical(tbsurv$status)

colSums(is.na(tbsurv))
# replace missing values with median
tbpbc <- tbpbc %>% mutate(resrate = replace(resrate, is.na(resrate), 
                                        median(resrate, na.rm = T)))
tbb <- tbb %>% mutate(ART = replace(ART, is.na(ART), 
                                    median(ART, na.rm = T)))
tbb <- tbb %>% mutate(subuse = replace(subuse, is.na(subuse), 
                                       median(subuse, na.rm = T)))
colSums(is.na(tbb))

## ----cleanup, results="asis", echo=FALSE------------------------------------------------
cls <- sapply(tbsurv, class) 
 
lbls <- c("time (days)",
    "Status",
    "Sex",
    "Age",
    "Weight",
    "HIV status",
    "TB class",
    "TB type",
    "ART",
    "Diabetes",
    "Alcohol status",
    "Smoking status",
    "Substance use")

# Build a table for data description
dta_labs <- data.frame(cbind(Variable = names(cls), 
                             Description = lbls, 
                             type = cls))

# Build a named vector for labeling figures later/
st_labs <- as.character(dta_labs$Description)
names(st_labs) <- names(cls)

# Print the descriptive table.
kable(dta_labs, 
      row.names = FALSE, 
      caption  = "`TB/HIV patients Treatment status.",
      booktabs = FALSE)

## ----eda, fig.cap="EDA variable plots. Points indicate variable value against the median home value variable. Points are colored according to the chas variable.", fig.width=7, fig.height=5----
# Use tidyr::gather to transform the data into long format.
#dta <- gather(tbb, variable, value, -medv, -tbstatus)

dta <- gather(tbsurv, variable, value, -time, -status)

# plot panels for each covariate colored by the logical status variable.
ggplot(dta) +
  geom_point(alpha = 0.4, aes(x = time, y = value, color = status)) +
  geom_smooth(aes(x = time, y = value), se = FALSE) + 
  labs(y = "", x = st_labs["time"]) +
  scale_color_brewer(palette = "Set2") +
  facet_wrap(~ variable, scales = "free_y", ncol = 3)

## ----randomforest, echo=TRUE------------------------------------------------------------
# Load the data, from the call:
tbsurv$status=factor(tbsurv$status)
tbsurv$time=factor(tbsurv$time)
rfsrc_tbb <- rfsrc(status ~ ., data = as.data.frame(tbsurv), 
                   importance = TRUE, err.block = 5)
predict(rfsrc_tbb)
# print the forest summary
rfsrc_tbb

## ----error, echo=TRUE, fig.cap="Random forest generalization error. OOB error convergence along the number of trees in the forest."----
# Plot the OOB errors against the growth of the forest.
gg_e <- gg_error(rfsrc_tbb)
gg_e <- gg_e %>% filter(!is.na(error))
class(gg_e) <- c("gg_error", class(gg_e))
plot(gg_e)

## ----rfsrc, echo=TRUE, fig.cap="OOB predicted median time values. Points are jittered to help visualize predictions for each observation. Boxplot indicates the distribution of the predicted values."----
# Plot predicted median home values.
gg_dta<-gg_rfsrc(rfsrc_tbb)
plot(gg_dta)

#plot(gg_rfsrc(rfsrc_tbb), alpha = .5) +
  #coord_cartesian(ylim = c(5, 49))

## ----vimp, echo=TRUE, fig.cap="Random forest VIMP plot. Bars are colored by sign of VIMP, longer blue bars indicate more important variables.", fig.width=7, fig.height=5----
# Plot the VIMP rankings of independent variables.

plot(gg_vimp(rfsrc_tbb), lbls = st_labs)

## ----minimaldepth, echo=TRUE, fig.cap="Minimal Depth variables in rank order, most important at the top. Vertical dashed line indicates the maximal minimal depth for important variables.", fig.width=7, fig.height=5----
# Load the data, from the call:

varsel_tbb <- var.select(rfsrc_tbb)

# Save the gg_minimal_depth object for later use.
gg_md <- gg_minimal_depth(varsel_tbb)

# plot the object
plot(gg_md, lbls = st_labs)

## ----minimalvimp, echo=TRUE, fig.cap="Comparing Minimal Depth and Vimp rankings. Points on the red dashed line are ranked equivalently, points below have higher VIMP, those above have higher minimal depth ranking. Variables are colored by the sign of the VIMP measure."----
# gg_minimal_depth objects contain information about
# both minimal depth and VIMP.

plot(gg_minimal_vimp(gg_md))

## ----variable, echo=FALSE, fig.cap="Variable dependence plot. Individual case predictions are marked with points. Loess smooth curve indicates the trend as the variables increase with shaded 95\\% confidence band.", fig.width=7, fig.height=5----
# Create the variable dependence object from the random forest
gg_v <- gg_variable(rfsrc_tbb)

# We want the top ranked minimal depth variables only,
# plotted in minimal depth rank order. 
xvar <- gg_md$topvars

# plot the variable list in a single panel plot
plot(gg_v, xvar = xvar, panel = TRUE, alpha = .5) +
  labs(y = st_labs["time"], x = "")

## ----chas, echo=TRUE, fig.cap="Variable dependence for Charles River logical variable."----
plot(gg_v, xvar = "status", alpha = 0.4) +
  labs(y = st_labs["time"])

# Load the data, from the call:
partial_tbb <- plot.variable(rfsrc_tbb,
                             xvar = gg_md$topvars,
                             partial = TRUE, sorted = FALSE,
                             show.plots = FALSE)

# generate a list of gg_partial objects, one per xvar.
gg_p <- gg_partial(partial_tbb)

# plot the variable list in a single panel plot
plot(gg_p, panel = TRUE) +
  labs(y = st_labs["time"], x = "") +
  geom_smooth(method = "loess", formula = "y ~ x", se = FALSE)

## ----interactions, echo=FALSE, fig.cap="Minimal depth variable interactions. Reference variables are marked with red cross in each panel. Higher values indicate lower interactivity with reference variable.", fig.width=7, fig.height=5, cache=TRUE----
# Load the data, from the call:
interaction_tbb <- find.interaction(rfsrc_tbb)

# Plot the results in a single panel.
plot(gg_interaction(interaction_tbb), 
     xvar = gg_md$topvars, panel = TRUE)

#############################################################
library(survivalROC)
library(survival)
## Fit a Cox model
coxph1 <- coxph(formula = Surv(time, status) ~ pspline(age, df = 4) + factor(hiv) +
                  factor(sex) + factor(diabetes),
                data    = tbsurv)
## Obtain the linear predictor
tbpbc$lp <- predict(coxph1, type = "lp")
tbpbc

### Cumulative case/dynamic control ROC
## Define a helper functio nto evaluate at various t
survivalROC_helper <- function(t) {
  survivalROC(Stime        = tbsurv$time,
              status       = tbsurv$status,
              marker       = tbsurv$lp,
              predict.time = t,
              method       = "NNE",
              span = 0.25 * nrow(tbsurv)^(-0.20))
}
## Evaluate every 180 days
survivalROC_data <- data_frame(t = 60 * c(1,2,3,4,5,6,7,8,9)) %>%
  mutate(survivalROC = map(t, survivalROC_helper),
         ## Extract scalar AUC
         auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
         ## Put cut off dependent values in a data_frame
         df_survivalROC = map(survivalROC, function(obj) {
           as_tibble(obj[c("cut.values","TP","FP")])
         })) %>%
  dplyr::select(-survivalROC) %>%
  unnest() %>%
  arrange(t, FP, TP)
## Plot
survivalROC_data %>%
  ggplot(mapping = aes(x = FP, y = TP)) +
  geom_point() +
  geom_line() +
  geom_label(data = survivalROC_data %>% dplyr::select(t,auc) %>% unique,
             mapping = aes(label = sprintf("%.3f", auc)), x = 0.5, y = 0.5) +
  facet_wrap( ~ t) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.key = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank())

#The 540-day ROC looks the best. However, this is because there were very few events up until this 
#point. After the last observed event (t≥300) the AUC stabilized at 0.656 and changes later. 
#The performance does not decay because the individuals with high risk score who did die keep 
#contributing to the performance.


#### Incident case/dynamic control ROC
library(risksetROC)
## Define a helper functio nto evaluate at various t
risksetROC_helper <- function(t) {
  risksetROC(Stime        = tbsurv$time,
             status       = tbsurv$status,
             marker       = tbsurv$lp,
             predict.time = t,
             method       = "Cox",
             plot         = FALSE)
}
## Evaluate every 180 days
risksetROC_data <- data_frame(t = 60 * c(1,2,3,4,5,6,7,8,9)) %>%
  mutate(risksetROC = map(t, risksetROC_helper),
         ## Extract scalar AUC
         auc = map_dbl(risksetROC, magrittr::extract2, "AUC"),
         ## Put cut off dependent values in a data_frame
         df_risksetROC = map(risksetROC, function(obj) {
           ## marker column is too short!
           marker <- c(-Inf, obj[["marker"]], Inf)
           bind_cols(data_frame(marker = marker),
                     as_data_frame(obj[c("TP","FP")]))
         })) %>%
  dplyr::select(-risksetROC) %>%
  unnest() %>%
  arrange(t, FP, TP)
## Plot
risksetROC_data %>%
  ggplot(mapping = aes(x = FP, y = TP)) +
  geom_point() +
  geom_line() +
  geom_col() +
  geom_label(data = risksetROC_data %>% dplyr::select(t,auc) %>% unique,
             mapping = aes(label = sprintf("%.3f", auc)), x = 0.5, y = 0.5) +
  facet_wrap( ~ t) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.key = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank())

#The difference is clearer in the later period. Most notably, only individuals that are in the risk set 
#at each time point contribute data. So there are fewer data points. The decay in performance is clearer 
#perhaps because among those who survived long enough the time-zero risk score is not as relevant. 
#Once there are not events left, the ROC essentially flat-lined.

#####################################################################



