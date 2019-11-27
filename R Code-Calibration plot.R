#####################################################
# Purpose: Plot calibration plot (Figure 5)
# Author: Chrianna Bharat
# Date: 28th November 2017
# Updated: 27th November 2019
#          - Additional annotation
#####################################################

# List of libraries
library(rms)
library(survival)
library(gtools)
library(caret)
library(plyr)
library(ggplot2)
library(stringi)
library(lme4)
library(MASS)

# Step 0: 
# Load data (not provided here, but called 'dat' in code below) 
#setwd()
#load()

# Create training (datTrain) and test (datTest) datasets 
set.seed(3456)
trainIndex <- createDataPartition(dat$blockage, #vector of outcomes
                                  p = .66,      #the percentage of data that goes to training
                                  list = FALSE, #IDs allocated to training dataset
                                  times = 1)    #number of partitions to create
datTrain <- dat[ trainIndex,]
datTest  <- dat[-trainIndex,]
dim(datTrain)
dim(datTest)

# Define final model (finalTrain) using training data 
ph.full<-coxph(Surv(SurvTime,blockage)~Depth+Length+gradient+sub_index+size_band+joint_type2+
                       Land_code+Road_type+road_prox+Purpose+Install_decade+soil_group,data=datTrain, na.action=na.omit)
finalTrain<-stepAIC(ph.full,direction="backward")
summary(finalTrain)    

# step 1:
# Estimate the linear predictor for the test data from the training model
datTest$lp_finalTrainModel<-predict(finalTrain,newdata=datTest,type="lp")

# Step 2:
# Fit a Cox PHM to the Test data LP as evaluated from the final model (determined from training set) 
# At 5 years, ~1826 days

finalTest.cph <- cph(Surv(SurvTime,blockage)~lp_finalTrainModel,data=datTest,
                     surv=TRUE, x=TRUE, y=TRUE,
                     time.inc=1826)

# Step 3: 
# Use cross-validation to get bias-corrected (overfitting - corrected) estimates of predicted vs. observed
# values based on subsetting predictions into intervals
# Use B=150 in practice Validate model for accuracy of predicting survival
# B is an upper limit on the number of resamples for which information is printed about which variables were selected in each model refit. 
# at t=1 Get Kaplan-Meier estimates by divided subjects into groups of
# size 200 (for other values of u must put time.inc=u in call to cph)
# M - group predicted u-time units survival into intervals containing m subjects on the average 
# See ?calibrate for detailed description 

cal1826.cross <- calibrate(finalTest.cph, B = 150, u = 1826, m = 2000, cmethod='KM',method="crossvalidation")  
cal1826.cross[,c(1:10)]

# 'Calibrate' object contents, including those used in plotting the calibration curve
#[,1] index.orig
#[,2] training
#[,3] test
#[,4] mean.optimism
#[,5] mean.corrected
#[,6] n
#[,7] mean.predicted
#[,8] KM
#[,9] KM.corrected
#[,10] str.err

# Step 4: 
# Plot mean predicted survival probability for group interval against actual survival probability (KM)
# Including bias corrected values, and confidence intervals

#pdf("Calibration plot at 5 years for test data.pdf", width = 15, height = 6)

par(mar=c(6,7,2,2),mgp=c(4,1,0))
plot(cal1826.cross[,8]~cal1826.cross[,7],
     type="b",pch=16,cex=1.5,
     ylim=c(0.795,0.98),xlim=c(0.81,0.975),
     ylab="Fraction Blockage-Free at 5 Years",xlab="Predicted 5 Year Blockage-Free",
     cex.lab=2,cex.axis=2.1)
segments(x0=0,y0=0,x1=1,y1=1,col="grey")

#bias corrected values
points(cal1826.cross[,7],cal1826.cross[,9],col="blue",pch=4,cex=2)

# hack: we draw arrows but with very special "arrowheads"
arrows(cal1826.cross[,7], 
       cal1826.cross[,8]-1.96*cal1826.cross[,10], 
       cal1826.cross[,7], 
       cal1826.cross[,8]+1.96*cal1826.cross[,10], 
       length=0.075, angle=90, code=3)

#dev.off()
