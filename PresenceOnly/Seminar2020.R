#--------------Data Preparation---
library(dismo)
library(maptools)
library(caret)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggspatial)
library(rgeos)
data(wrld_simpl)
#loading the predictor variables as a raster stack
predictors <- stack(list.files(file.path(system.file(package="dismo"), 'ex'), pattern='grd$', full.names=TRUE ))
#loading the presence points of bradypus as coordinates
file <- file.path(system.file(package="dismo"), "ex/bradypus.csv")
bradypus <- read.table(file,  header=TRUE,  sep=',')
bradypus <- bradypus[,-1]
#loading a worldmap with borders
world <- ne_countries(scale = 'medium', returnclass = "sf")
theme_set(theme_bw())
#plot the presence points
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-107,-34),ylim = c(-40,32))+
  geom_point(data = bradypus, aes(x = lon, y = lat), color = 'dark green')


##------Find Correlation-------

covar <- layerStats(predictors,stat = 'pearson',na.rm = T)
covar
predictors <- dropLayer(predictors,c('bio1','bio16','bio6','bio12','bio7'))
#----pseudo absence------
ext <- extent(-90, -32, -33, 23)
#creating all true absence points
dummy <- predictors$bio17
values(dummy) <- NA
dummy[cellFromXY(dummy,bradypus)] <- 1
b <- raster::buffer(dummy, width=200000) 
b[is.na(b)] <- 0
b[is.na(b!=predictors$bio5)] <- NA
b[b == 1] <- 0.5
b[cellFromXY(b,bradypus)] <- 1

pa <- rasterToPoints(b)
colnames(pa) <- c("lon","lat","presence/absence")
#---------------training testing data
pres <- pa[pa[,3]==1,]
pres <- pres[,1:2]
group <- kfold(pres, 5)
pres_train <- pres[group != 1, ]
pres_test <- pres[group == 1, ]
#settting pool of all pseudo absence points to be all points in our studied area
backg <- pa[(pa[,3]==0)&(pa[,1] >= -90)&(pa[,1]<=-32)&(pa[,2]>=-33)&(pa[,2]<=23),]
backg <- backg[,1:2]
colnames(backg) <- c("lon","lat")
#plotting preence points and pool of all available pseudo absence points
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-107,-34),ylim = c(-40,32))+
  geom_point(data = data.frame(bradypus),aes(x= lon,y=lat,color = 'green'), size = .8 )+
  geom_point(data = data.frame(backg), aes(x = lon, y = lat, color = 'yellow'), size = .8)+
  scale_color_identity(guide = 'legend', name = "", breaks = c("green","yellow"),
                       labels = c("presence","pseudo absence"))
#randomly selecting 116 pseudo absence points from the pool of all true absences
backg <- backg[sample(nrow(backg),116),]
group <- kfold(backg, 5)
#dividing pseudo absence points into training and testing set
backg_train <- backg[group != 1, ]
backg_test <- backg[group == 1, ]


## ---- Visualizing Train/test Data--------------------------------------------------------------
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-107,-34),ylim = c(-40,32))+
  geom_point(data = data.frame(backg_train),aes(x= lon,y=lat,color = 'green'), size = .8 )+
  geom_point(data = data.frame(backg_test),aes(x= lon,y=lat,color = 'dark green'), size = .8)+
  geom_point(data = data.frame(pres_train),aes(x= lon,y=lat, color = 'red'), size = .8)+
  geom_point(data = data.frame(pres_test),aes(x= lon,y=lat, color = 'orange'), size = .8)+
  scale_color_identity(guide = "legend", name = "", breaks = c("green", "dark green","red","orange"),
                       labels = c("Back.Train","Back.Test","Pres.Train","Pres.Test"))

## ------------------------------------------------------------------

train <- rbind(pres_train, backg_train)
pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
#getting predictors values at training sites
envtrain <- extract(predictors, train)
envtrain <- data.frame( cbind(pa=pb_train, envtrain) )
envtrain[,'biome'] = factor(envtrain[,'biome'], levels=1:14)
#getting predictor values at test sites
testpres <- data.frame( extract(predictors, pres_test) )
testbackg <- data.frame( extract(predictors, backg_test) )
testpres[ ,'biome'] = factor(testpres[ ,'biome'], levels=1:14)
testbackg[ ,'biome'] = factor(testbackg[ ,'biome'], levels=1:14)

## ---- MaxEnt-----------------------------
maxent()
#creating MaxEnt model
xm <- maxent(predictors, pres_train, factors='biome')
#evaluating the model
e.max <- evaluate(pres_test, backg_test, xm, predictors)
e.max
#create a numerical distribution map
px <- predict(predictors, xm, ext=ext, progress='')

#plot the numerical map
ggplot(data = world) +
  layer_spatial(px)+ scale_fill_distiller(palette = "RdYlGn", name = "MaxEnt prediction")
#find threshold to convter into binary
max.tr <- threshold(e.max, 'spec_sens')
#plot binary distribution map
ggplot() + 
  layer_spatial(px > max.tr)+
  scale_fill_brewer(palette = "Greens",name = "Presence/Absence",labels=c("Absence","Presence")) +
  theme_dark()
#plot binary distribution map with presence points
ggplot() + 
  layer_spatial(px > max.tr)+
  scale_fill_brewer(palette = "Greens",name = "Presence/Absence",labels=c("Absence","Presence")) +
  geom_point(data = bradypus, aes(x=lon, y = lat))+
  theme_dark()

##------------Classification Tree--------
library(rpart)
#create CT model
ct <- rpart(pa~bio5+bio8+bio17, data = envtrain)
#create numerical distribution map
pc <- predict(predictors,ct,ext=ext)
#evaluate the model
e.ct <- evaluate(testpres,testbackg,ct)
e.ct
#find threshold to convert into binary 
c.tr <- threshold(e.ct,'spec_sens')
#plot numerical distribution map
ggplot()+
  layer_spatial(pc)+ scale_fill_distiller(palette = "RdYlGn", name = "CT prediction")
#plot binary distribution map
ggplot()+
  layer_spatial(pc > c.tr)+
  scale_fill_brewer(palette = "Greens",name = "Presence/Absence",labels=c("Absence","Presence")) +
  theme_dark()
#plot binary distribution map with presence points
ggplot()+
  layer_spatial(pc > c.tr)+
  scale_fill_brewer(palette = "Greens",name = "Presence/Absence",labels=c("Absence","Presence")) +
  geom_point(data = bradypus, aes(x=lon,y=lat))+
  theme_dark()
#compare binary pred. with MaxEnt
ggplot()+
  layer_spatial((pc > c.tr)==(px > max.tr))+
  scale_fill_brewer(palette = "Reds",name = "MaxEnt-CT",labels=c("Diff","Equal")) +
  theme_dark()
#how many cells had the same binary predictions(percentage)
abc <- as.vector((pc > c.tr)==(px > max.tr))
length(abc[(abc == T)&!(is.na(abc))])/length((px > max.tr)[!is.na(px > max.tr)])
#Diff. in numerical predictions
ggplot()+
  layer_spatial(abs(pc-px))+scale_fill_distiller(palette = "RdYlGn", name = "CT-MaxEnt Difference")
#average diff. in numerical predictions
mean(as.vector(abs(pc-px)),na.rm = T)
## ---- Random Forest-----------------------------------
library(randomForest)
#create RF model
rf1 <- randomForest(pa ~ bio5  + bio8  + bio17, data=envtrain)
#evaluate model
e.rf <- evaluate(testpres, testbackg, rf1)
e.rf
#create prediction map
pr <- predict(predictors, rf1, ext=ext)
#plot numeical distribution map
ggplot()+
  layer_spatial(pr)+scale_fill_distiller(palette = "RdYlGn", name = "RT predicition")

#find therhsold for conversion into binary
rf.tr <- threshold(e.rf, 'spec_sens')
#plot binary distribution map
ggplot()+
  layer_spatial(pr > rf.tr)+
  scale_fill_brewer(palette = "Greens",name = "Presence/Absence",labels=c("Absence","Presence")) +
  theme_dark()
#plot binary distribution map with presence points
ggplot()+
  layer_spatial(pr > rf.tr)+
  scale_fill_brewer(palette = "Greens",name = "Presence/Absence",labels=c("Absence","Presence")) +
  geom_point(data = bradypus, aes(x = lon, y = lat)) +
  theme_dark()
#compare binary predictions with MaxEnt
ggplot()+
  layer_spatial((pr > rf.tr)==(px > max.tr))+
  scale_fill_brewer(palette = "Reds",name = "MaxEnt-RF",labels=c("Diff","Equal")) +
  theme_dark()
#percentage of cells with same binary predictions
abc <- as.vector((pr > rf.tr)==(px > max.tr))
length(abc[(abc == T)&!(is.na(abc))])/length((px > max.tr)[!is.na(px > max.tr)])
#differnce in numerical predictions with MaxEnt
ggplot()+
  layer_spatial(abs(pr-px))+scale_fill_distiller(palette = "RdYlGn", name = "RF-MaxEnt Difference")
#average diff. in numerical predictions
mean(as.vector(abs(pr-px)),na.rm = T)#0.1758913

## ---- Suppport Vector Machine----------------------------------
library(kernlab)
#cereating SVM model
svm <- ksvm(pa~ bio5+bio8+bio17, data=envtrain)
#evaulate model
esv <- evaluate(testpres, testbackg, svm)
esv
#create prediction map
ps <- predict(predictors, svm, ext=ext)
#plot numerical distribution map
ggplot()+
  layer_spatial(ps)+scale_fill_distiller(palette = "RdYlGn", name = "SVM predicition")
#find threshold to convert into binary data
s.tr <- threshold(esv, 'spec_sens')
#plot binary distribution map
ggplot()+
  layer_spatial(ps > s.tr)+
  scale_fill_brewer(palette = "Greens",name = "Presence/Absence",labels=c("Absence","Presence")) +
  theme_dark()
#plot binary distribution map with presence points
ggplot()+
  layer_spatial(ps > s.tr)+
  scale_fill_brewer(palette = "Greens",name = "Presence/Absence",labels=c("Absence","Presence")) +
  geom_point(data = bradypus, aes(x=lon,y=lat))+
  theme_dark()
#compare binary predictions with MaxEnt
ggplot()+
  layer_spatial((ps > s.tr)==(px > max.tr))+
  scale_fill_brewer(palette = "Reds",name = "MaxEnt-SVM",labels=c("Diff","Equal")) +
  theme_dark()
#percentage of cells iwth same binary predictions as MaxEnt
abc <- as.vector((ps > s.tr)==(px > max.tr))
length(abc[(abc == T)&!(is.na(abc))])/length((px > max.tr)[!is.na(px > max.tr)])#0.8425171
#Diff. in numerical predictions with MAXEnt
ggplot()+
  layer_spatial(abs(ps-px))+scale_fill_distiller(palette = "RdYlGn", name = "SVM-MaxEnt Difference")
#avegreage diff. in numerical predictions
mean(as.vector(abs(ps-px)),na.rm = T)#0.1726927

# ----------ModelPerformanceEvaluation------------
predic.pres <- data.frame(extract(predictors,pres_test))
predic.backg <- data.frame(extract(predictors,backg_test))

predict1 <- cbind(rep(1,length(predic.pres)),predict(rf1, newdata = predic.pres))
predict2 <- cbind(rep(0,length(predic.backg)),predict(rf1,newdata = predic.backg))
rf.test <- rbind(predict1,predict2)
predict1 <- cbind(rep(1,length(predic.pres)),predict(ct, newdata = predic.pres))
predict2 <- cbind(rep(0,length(predic.backg)),predict(ct,newdata = predic.backg))
ct.test <- rbind(rbind(predict1,predict2))
predict1 <- cbind(rep(1,length(predic.pres)),predict(svm, newdata = predic.pres))
predict2 <- cbind(rep(0,length(predic.backg)),predict(svm,newdata = predic.backg))
svm.test <- rbind(rbind(predict1,predict2))
#calculating RMSE
#RF
RMSE(rf.test[,2],rf.test[,1])
#CT
RMSE(ct.test[,2],ct.test[,1])
#SVM
RMSE(svm.test[,2],svm.test[,1])
#Confusion matrix
#RF
rf.test.binary <- data.frame(rf.test)
rf.test.binary[rf.test.binary$X2 < rf.tr,]$X2 <- 0
rf.test.binary[rf.test.binary$X2 >= rf.tr,]$X2 <- 1
rf.cmatrix <- table(rf.test.binary$X1,rf.test.binary$X2)
rf.cmatrix
#CT
ct.test.binary <- data.frame(ct.test)
ct.test.binary[ct.test.binary$X2 < c.tr,]$X2 <- 0
ct.test.binary[ct.test.binary$X2 >= c.tr,]$X2 <- 1
ct.cmatrix <- table(ct.test.binary$X1,ct.test.binary$X2)
ct.cmatrix
#SVM
svm.test.binary <- data.frame(svm.test)
svm.test.binary[svm.test.binary$X2 < s.tr,]$X2 <- 0
svm.test.binary[svm.test.binary$X2 >= s.tr,]$X2 <- 1
svm.cmatrix <- table(svm.test.binary$X1,svm.test.binary$X2)
svm.cmatrix
# creating the ROC graph
library(pROC)
roc <- pROC::roc(response = rf.test[,1],predictor = rf.test[,2], plot = T,legacy.axes = T, percent = T,
                 xlab = "False Positive Percentage", ylab = "True Positive Percentage",lwd = 4,
                 print.auc = T, col = 'red')
pROC::roc(response = ct.test[,1],predictor = ct.test[,2],plot = T, percent = T, legacy.axes = T,
          lwd = 4, print.auc = T, add = T, col = 'blue',print.auc.y=40)
pROC::roc(response = svm.test[,1],predictor = svm.test[,2], plot = T,percent = T, legacy.axes = T,
          lwd = 4, print.auc = T, add = T, col = 'green', print.auc.y=30)
legend("bottomright", legend = c("RF","CT","SVM"),col = c('red','blue','green'), lwd = 4)

