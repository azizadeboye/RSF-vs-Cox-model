
library(autoReg)
#library(moonBook) # For use of example data acs

gaze(status~.,data=tbclinic)

tbclinic$status.factor=factor(tbclinic$status,labels=c("Alive","Died"))
tbclinic$sex.factor=factor(tbclinic$sex,labels=c("Female","Male"))
tbclinic$hiv.factor=factor(tbclinic$hiv,labels=c("Positive","Negative"))
tbclinic$dclass.factor=factor(tbclinic$dclass,labels=c("EPTB","PTB"))
tbclinic$tbtype.factor=factor(tbclinic$tbtype,labels=c("DR-TB","MDR-TB"))
tbclinic$ART.factor=factor(tbclinic$ART,labels=c("No","Yes"))
tbclinic$diabetes.factor=factor(tbclinic$diabetes,labels=c("No","Yes"))
tbclinic$marital.factor=factor(tbclinic$marital,labels=c("Single","Maried", "Divorced/separated", "Widow"))
tbclinic$res.factor=factor(tbclinic$res,labels=c("Rural","Urban"))
tbclinic$alcohol.factor=factor(tbclinic$alcohol,labels=c("No","Yes"))
tbclinic$smoking.factor=factor(tbclinic$smoking,labels=c("No","Yes"))
tbclinic$subuse.factor=factor(tbclinic$subuse,labels=c("No","Yes"))

tbclinic$status.factor=setLabel(tbclinic$status.factor,"Mortality")
tbclinic$sex.factor=setLabel(tbclinic$sex, "Gender")
tbclinic$hiv.factor=setLabel(tbclinic$hiv, "HIV status")
tbclinic$dclass.factor=setLabel(tbclinic$dclass, "Disease class")
tbclinic$tbtype.factor=setLabel(tbclinic$tbtype, "TB type")
tbclinic$ART.factor=setLabel(tbclinic$ART, "ART")
tbclinic$diabetes.factor=setLabel(tbclinic$diabetes, "Diabetes")
tbclinic$alcohol.factor=setLabel(tbclinic$alcohol, "Alcohol")
tbclinic$smoking.factor=setLabel(tbclinic$smoking, "Smoking")
tbclinic$subuse.factor=setLabel(tbclinic$subuse, "SUbstance use")

library(dplyr) # for use of `%>%`
ft=gaze(status~.,data=tbclinic) %>% myft()
ft

library(rrtable)
#  
table2pptx(ft)
table2docx(ft)
