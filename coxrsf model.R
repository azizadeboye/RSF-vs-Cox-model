
# Machine learning survival model:cox and random forest

## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(survex)
library(survival)

set.seed(123)

cph <- coxph(Surv(time, status)~., data = tbsurv, model = TRUE, x = TRUE)
summary(cph)
cph_exp <- explain(cph)

rsf <- randomForestSRC::rfsrc(Surv(time, status)~., data = tbsurv)
rsf_exp <- explain(rsf)

## -----------------------------------------------------------------------------
predict(cph_exp, tbsurv[1:2,], output_type="risk")
predict(rsf_exp, tbsurv[1:2,], output_type="risk")

predict(cph_exp, tbsurv[1:2,], output_type="survival", times=seq(1, 600, 100))
predict(rsf_exp, tbsurv[1:2,], output_type="survival", times=seq(1, 600, 100))

predict(cph_exp, tbsurv[1:2,], output_type="chf", times=seq(1, 600, 100))
predict(rsf_exp, tbsurv[1:2,], output_type="chf", times=seq(1, 600, 100))

## ----warning=FALSE------------------------------------------------------------
mp_cph <- model_performance(cph_exp)
mp_rsf <- model_performance(rsf_exp)

plot(mp_cph, mp_rsf)

## ---- include=FALSE-----------------------------------------------------------
dev.off()

## -----------------------------------------------------------------------------
plot(mp_cph, mp_rsf, metrics_type="scalar")

## ---- include=FALSE-----------------------------------------------------------
dev.off()

## -----------------------------------------------------------------------------
model_parts_rsf <- model_parts(rsf_exp)
model_parts_cph <- model_parts(cph_exp)

plot(model_parts_cph,model_parts_rsf)

## ---- include=FALSE-----------------------------------------------------------
dev.off()

## -----------------------------------------------------------------------------
model_parts_rsf_auc <- model_parts(rsf_exp, loss_function=loss_one_minus_cd_auc, type="difference")
model_parts_cph_auc <- model_parts(cph_exp, loss_function=loss_one_minus_cd_auc, type="difference")

# NOTE: this may take a long time, so a progress bar is available. To enable it, use:
# progressr::with_progress(model_parts(rsf_exp, loss_function=loss_one_minus_cd_auc, type="difference"))

plot(model_parts_cph_auc,model_parts_rsf_auc)

## ---- include=FALSE-----------------------------------------------------------
dev.off()

## ---- fig.height=18-----------------------------------------------------------
model_profile_cph <- model_profile(cph_exp, categorical_variables=c("dclass", "diabetes", "sex",
                                                                     "ART", "Smoking"))
plot(model_profile_cph, facet_ncol = 1)

## ---- include=FALSE-----------------------------------------------------------
dev.off()

## ---- fig.height=18-----------------------------------------------------------
model_profile_rsf <- model_profile(rsf_exp, categorical_variables=c("dclass", "diabetes", "sex",
                                                                    "ART", "Smoking"))
plot(model_profile_rsf, facet_ncol = 1, numerical_plot_type = "contour")

## ---- include=FALSE-----------------------------------------------------------
dev.off()

## -----------------------------------------------------------------------------
predict_parts_cph_32 <- predict_parts(cph_exp, tbsurv[32,])
predict_parts_rsf_32 <- predict_parts(rsf_exp, tbsurv[32,])
plot(predict_parts_cph_32, predict_parts_rsf_32)

## ---- include=FALSE-----------------------------------------------------------
dev.off()

## -----------------------------------------------------------------------------
predict_parts_cph_12 <- predict_parts(cph_exp, tbsurv[12,])
predict_parts_rsf_12 <- predict_parts(rsf_exp, tbsurv[12,])
plot(predict_parts_cph_12, predict_parts_rsf_12)

## ---- include=FALSE-----------------------------------------------------------
dev.off()

#########
## -----------------------------------------------------------------------------
predict_parts_cph_12_lime <- predict_parts(cph_exp, tbsurv[12,], type="survlime")
predict_parts_rsf_12_lime <- predict_parts(rsf_exp, tbsurv[12,], type="survlime")
plot(predict_parts_cph_12_lime, type="local_importance")

## ---- include=FALSE-----------------------------------------------------------
dev.off()

## -----------------------------------------------------------------------------
plot(predict_parts_rsf_12_lime, type="local_importance")

## ---- include=FALSE-----------------------------------------------------------
dev.off()

## ---- fig.height=18-----------------------------------------------------------
predict_profile_cph_32 <- predict_profile(cph_exp, tbsurv[32,], categorical_variables=c("dclass", "diabetes", "sex",
                                                                                        "ART", "Smoking"))
plot(predict_profile_cph_32, facet_ncol=1)

## ---- include=FALSE-----------------------------------------------------------
dev.off()

## ---- fig.height=18-----------------------------------------------------------
predict_profile_rsf_32 <- predict_profile(rsf_exp, tbsurv[32,], categorical_variables=c("sex", "smoking"))
plot(predict_profile_rsf_32, facet_ncol=1)
