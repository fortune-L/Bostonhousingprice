# Call R libraries

library(ROCR)
library(car)
library(mgcv)
library(dplyr)
library(tidyr)
library(spaMM)
library(ProbitSpatial)
library(arm)
library(spatialreg)
library(spdep)
library(sf)
library(spmoran)

# Call Boston data
data(boston)
data <- boston.c
dim(data)
head(data,5)

# Assign x and y variables
y <- data[, "CMEDV"]
x<- data[,c("CRIM","AGE","ZN","DIS","RAD","NOX",  "TAX","RM", "PTRATIO", "B","INDUS", "CHAS","LSTAT" )]

# creating the weight matrix
long_lat <- data[,c("LON", "LAT")]
long_lat

# Converts R data frame to a matrix 
coords <- as.matrix(long_lat) 
coords <- coordinates(coords)
coords


library(ggplot2)
library(ggmap)
library(mapproj)

BostonLL <-c(-71.30, 42.00, -70.80, 42.40)
map <- get_map(location = BostonLL, zoom = 11)
plot(map)
# map for Median value of owner-occupied homes
mapPoints <- ggmap(map) + geom_point(aes(x = LON, y = LAT,color=CMEDV),data = data, alpha = 0.7,size=3)+ scale_color_gradient(low="#9ebcda",high="#8856a7")+labs(title="Median value of owner-occupied homes in $1000's",y="Latitude",x="Longtitude",color="House Prices (in 1000$'s)" )
mapPoints

cutpoints<-quantile(data$CMEDV,seq(0,1,length=4),na.rm=TRUE)
data$CMEDVQuantiles <- cut(data$MEDV,breaks=cutpoints,include.lowest=TRUE,labels =c("Low priced","Mid priced","High priced"))

mapPoints2 <- ggmap(map) + geom_point(aes(x = LON, y = LAT,color=CMEDVQuantiles),data = data, alpha = 0.7,size=3) + ggtitle("Median value of owner-occupied homes in $1000's") +facet_grid(.~CMEDVQuantiles) +labs(title="Median House Prices Categorized",y="Latitude",x="Longtitude",color="House Prices (in 1000$'s)" )
mapPoints2

# OTHER VARIABLES

# map for NOX
mapPointsNOX <- ggmap(map) + geom_point(aes(x = LON, y = LAT,color=NOX),data = data, alpha = 1,size=3)+scale_color_gradient(low="#8856a7",high="#31a354")+labs(title="NOX",y="Latitude",x="Longtitude",color="NOX Conc. (parts per 10 million)" )
mapPointsNOX
# map for CRIM
mapPointsCRIM <- ggmap(map) + geom_point(aes(x = LON, y = LAT,color=CRIM),data = data, alpha = 1,size=3)+scale_color_gradient(low="#8856a7",high="#31a354")+labs(title="CRIM",y="Latitude",x="Longtitude",color="crime rate (per capita)")
mapPointsCRIM
#map for RM 
mapPointsRM <- ggmap(map) + geom_point(aes(x = LON, y = LAT,color=RM),data = data, alpha = 1,size=3)+scale_color_gradient(low="#8856a7",high="#31a354")+labs(title="RM",y="Latitude",x="Longtitude",color="Av no.of rooms. (per dwelling)")
mapPointsRM
# map for TAX
mapPointsTAX <- ggmap(map) + geom_point(aes(x = LON, y = LAT,color=TAX),data = data, alpha = 1,size=3)+scale_color_gradient(low="#8856a7",high="#31a354")+labs(title="TAX",y="Latitude",x="Longtitude",color="tax rate(per USD 10,000)")
mapPointsTAX
# map for DIS
mapPointsDIS <- ggmap(map) + geom_point(aes(x = LON, y = LAT,color=DIS),data = data, alpha = 1,size=3)+scale_color_gradient(low="#8856a7",high="#31a354")+labs(title="DIS",y="Latitude",x="Longtitude",color="Weighted distances")
mapPointsDIS
# map for AGE
mapPointsAGE <- ggmap(map) + geom_point(aes(x = LON, y = LAT,color=AGE),data = data, alpha = 1,size=3)+scale_color_gradient(low="#8856a7",high="#31a354")+labs(title="AGE",y="Latitude",x="Longtitude",color="owner-occupied units")
mapPointsAGE


# k-Nearest Neighbour weights,k=5
knn <- knearneigh(coords, k=5, longlat= TRUE)  
nb <- knn2nb(knn) # convert to neighbour objects
nb 

listw <- nb2listw(nb, style = "W")
plot(nb, coords, lwd=.2, col="blue", cex = .05)
W <- as(as_dgRMatrix_listw(listw),  "CsparseMatrix") 
dim(W)

### The binary model: Binary variable equals to 1 if the price >25, 0 otherwise
y_b = ifelse(y>25,1,0)  
y_b
# Using the Generalized additive models function gam(), to develop the binary regression model
# with only parametric components, there is no spatial interaction.

GAM_1 <- gam(y_b ~   -1 + x$AGE  + x$ZN + x$DIS + x$RAD +x$NOX + 
               x$TAX + x$RM + x$PTRATIO +x$B + x$INDUS + factor(x$CHAS) + 
               x$LSTAT + x$CRIM, data = data, family = binomial(link = 'probit'))

summary(GAM_1)
AIC(GAM_1)

# Semiparametric probit model, with no spatial interaction
GAM_2 <- gam(y_b ~   -1 + x$DIS + x$RAD + 
               x$TAX + x$RM + x$PTRATIO + x$INDUS + factor(x$CHAS) + 
               x$LSTAT + s(x$NOX), data = data, family = binomial(link = 'probit'))

summary(GAM_2)
AIC(GAM_2)



# Adding spatial autocorrelation(error U)

# keep smooth terms 
terms = as.data.frame(predict(GAM_2, type = 'terms'))
data$NOX.Z <- terms$`s(x$NOX)`

## fit semi-parametric probit model with smooth variables
probispatialfit3 <- ProbitSpatialFit(y_b ~   -1 + x$DIS + x$RAD + x$TAX + 
                                       x$RM + x$PTRATIO + x$INDUS + factor(x$CHAS)+ x$LSTAT + data$NOX.Z,
                                     W = W, DGP = 'SEM', method = 'conditional', varcov = 'varcov')
summary(probispatialfit3)

library(MASS)
library(mapdata)
library(maptools)
library(MapGAM)
library(cartography)
library(rgdal)
require(devtools)



#plot(x$TOWN)
COORDS_LONG <- data$LON
COORDS_LAT  <- data$LAT 
lastdata <- cbind(y_b ,COORDS_LONG, COORDS_LAT, data)

gamgrid <- predgrid(lastdata, x$TOWN)
nrow(gamgrid)

# Estimation because we are using the same data for training and test

pred <- predict(probispatialfit3, probispatialfit3$X, type = 'response', WSO = W)
colormap(list(fit=pred,grid=data.frame(X=data$LON,Y=data$LAT)),col.seq=rev(rainbow(201,start=0,end=0.66)),ptsize = 1,map=x$TOWN)

lastdata$pred=ifelse(pred>0.8,1,0) # prediction of y_b
mapPoints <- ggmap(map) + geom_point(aes(x = LON, y = LAT,color=pred),data = lastdata, alpha = 0.7,size=3)+ scale_color_gradient(low="#9ebcda",high="#8856a7")+labs(title="Prediction of owner-occupied homes in $1000's",y="Latitude",x="Longtitude",color="House Prices (in 1000$'s)" )
mapPoints


## To do prediction we divide the sample size in training and test 

library(caTools)

lastdata <- cbind(y_b,COORDS_LONG, COORDS_LAT, data)
split=sample.split(lastdata$y_b,SplitRatio = 0.5)
train=subset(lastdata,split==TRUE) # Help to build the model
test=subset(lastdata,split==FALSE) # Help to predict the model

# TRAINING


GAM_2 <- gam(y_b ~   -1 + DIS + RAD + TAX+
               RM + PTRATIO + INDUS + factor(CHAS)+ LSTAT + s(NOX),
             data = train, family = binomial(link = 'probit'))

summary(GAM_2)
AIC(GAM_2)

termplot(GAM_2, se=TRUE, col.term=1, col.se=1, data=train)# s(NOX) plot

# new spatial matrix for test 

long_lat2 <- test[,c("COORDS_LONG", "COORDS_LAT")]
long_lat2
# Converts R data frame to a matrix 
coords2 <- as.matrix(long_lat2) 
coords <- coordinates(coords2)
coords
knn <- knearneigh(coords, k=5, longlat= TRUE) # longlat=TRUE for spherical coordinates 
nb <- knn2nb(knn) # convert to neighbour objects
listw <- nb2listw(nb, style = "W")
plot(nb, coords, lwd=.2, col="blue", cex = .05)
W1 <- as(as_dgRMatrix_listw(listw),  "CsparseMatrix") 
dim(W1)

# keep smooth terms 
terms = as.data.frame(predict(GAM_2, type = 'terms'))
test$NOX.Z <- terms$`s(NOX)`
train$NOX.Z <- terms$`s(NOX)`


# training semi-parametric using the data set, see the train dataframe

probispatialfitt <- ProbitSpatialFit(train$y_b ~   -1 + train$DIS + train$RAD + train$TAX+
                                       train$RM + train$PTRATIO + train$INDUS + factor(train$CHAS)+ train$LSTAT + train$NOX.Z, data = train,
                                     W = W1, DGP = 'SEM', method = 'conditional', varcov = 'varcov')


summary(probispatialfitt)


# only to have the x of the test in the frame as the probispatialfitt 

probispatialfitt2 <- ProbitSpatialFit(test$y_b ~   -1 + test$DIS + test$RAD + test$TAX+
                                       test$RM + test$PTRATIO + test$INDUS + factor(test$CHAS)+ test$LSTAT  + test$NOX.Z, data = test,
                                     W = W1, DGP = 'SEM', method = 'conditional', varcov = 'varcov')

# TESTING 

#X=probispatialfitt2$X are used to predict y_b using the model estimated with the train

# probability prediction probabilitiers
pred1 <- predict(probispatialfitt,probispatialfitt2$X, type = 'response', WSO = W)
colormap(list(fit=pred1,grid=data.frame(X=data$LON,Y=data$LAT)),col.seq=rev(rainbow(201,start=0,end=0.66)),ptsize = 1,map=x$TOWN)

## prediction with GAM2
newdata=data.frame(test[,c("DIS","RAD", "RM","TAX","PTRATIO", "INDUS", "CHAS", "LSTAT", "NOX")])

predgam <- predict(GAM_2,newdata, type="response")


## prediction
test$pred=ifelse(pred1>0.8,1,0) ## you can change the threshold 0.8, 0.9,..
test$predgam=ifelse(predgam>0.8,1,0)
mapPoints <- ggmap(map) + geom_point(aes(x =LON, y = LAT,color=pred),data = test, alpha = 0.7,size=3)+ scale_color_gradient(low="#9ebcda",high="#8856a7")+labs(title="Prediction of owner-occupied homes",y="Latitude",x="Longtitude",color="House Prices (in 1000$'s)" )
mapPoints


mapPoints <- ggmap(map) + geom_point(aes(x =LON, y = LAT,color=predgam),data = test, alpha = 0.7,size=3)+ scale_color_gradient(low="#9ebcda",high="#8856a7")+labs(title="Prediction of owner-occupied homes",y="Latitude",x="Longtitude",color="House Prices (in 1000$'s)" )
mapPoints


#true values
mapPoints <- ggmap(map) + geom_point(aes(x =LON, y = LAT,color=y_b),data = test, alpha = 0.7,size=3)+ scale_color_gradient(low="#9ebcda",high="#8856a7")+labs(title="Owner-occupied homes",y="Latitude",x="Longtitude",color="House Prices (in 1000$'s)" )
mapPoints

#confusion table
table(test$y_b,test$pred) 

table(test$y_b,test$predgam) 
#prediction rate

sum(diag(table(test$y_b,test$pred)))/sum(table(test$y_b,test$pred))

sum(diag(table(test$y_b,test$predgam)))/sum(table(test$y_b,test$predgam))

# Plot the confusion matrices with different thresh (Sophie)


#####Compare also with the prediction done by GAM2

## GAUSSIAN SPATIAL ADDITIVE MIXED MODELS, with significant variables
samp    <- sample( dim(boston.c )[1], 405) # 80% of data for training
d       <- boston.c[ samp, ]    ## Data at observed sites
y	      <- d[,"CMEDV"]
x       <- d[,c("PTRATIO", "INDUS", "CHAS","RM","TAX","RAD","DIS","LSTAT")]
coords  <- d[,c("LON", "LAT")]


d0      <- boston.c[-samp, ]    ## Data at unobserved sites , test data
x0      <- d0[,c("PTRATIO", "INDUS", "CHAS","RM","TAX","RAD","DIS", "LSTAT")]
coords0  <- d0[,c("LON", "LAT")]
### Eigenvector spatial filtering (ESF)
meig <- meigen(coords=coords)
mod1 <- esf(y=y,x=x,meig=meig, vif=10)
mod1

### Random effects ESF (RE-ESF)
meig <- meigen(coords=coords)
mod2<- resf(y = y, x = x, meig = meig)
mod2
meig0 	<- meigen0( meig = meig, coords0 = coords0 )
pred12  <- predict0( mod = mod2, x0 = x0, meig0 = meig0 )
mean((pred12$pred[,1]-d0$CMEDV)^2) ## mean square prediction error
d0$pred12=pred12$pred[,1]

mapPoints <- ggmap(map) + geom_point(aes(x =LON, y = LAT,color=pred12),data =d0, alpha = 0.7,size=3)+ scale_color_gradient(low="#9ebcda",high="#8856a7")+labs(title="Prediction of owner-occupied homes",y="Latitude",x="Longtitude",color="House Prices (in 1000$'s)" )
mapPoints
#true values
mapPoints <- ggmap(map) + geom_point(aes(x =LON, y = LAT,color=CMEDV),data = d0, alpha = 0.7,size=3)+ scale_color_gradient(low="#9ebcda",high="#8856a7")+labs(title="Owner-occupied homes",y="Latitude",x="Longtitude",color="House Prices (in 1000$'s)" )
mapPoints

### Extended models
#Models with non-spatially varying coefficients (coefficients varying wrt covariate value)
meig <- meigen(coords=coords)
mod3 <- resf(y = y, x = x, meig = meig, nvc=TRUE)
mod3

############ Spatial prediction of y and spatially varying coefficients

samp    <- sample( dim(boston.c )[1], 405)

d       <- boston.c[ samp, ]    ## Data at observed sites
y	      <- d[, "CMEDV"]
x       <- d[,c("RM", "LSTAT")]
xconst  <- d[,c( "PTRATIO", "INDUS", "CHAS","DIS","TAX","RAD")]
coords  <- d[,c("LON", "LAT")]

d0      <- boston.c[-samp, ]    ## Data at unobserved sites , test data
x0      <- d0[,c("RM", "LSTAT")]
xconst0 <- d0[,c( "PTRATIO", "INDUS", "CHAS","DIS","TAX","RAD")]
coords0 <- d0[,c("LON", "LAT")]

############ Model estimation
meig 	  <- meigen( coords = coords )
mod4	    <- resf_vc(y=y, x=x, xconst=xconst, meig=meig )
mod4
meig0 	<- meigen0( meig = meig, coords0 = coords0 )
pred0   <- predict0_vc( mod = mod4, x0 = x0, xconst0=xconst0, meig0 = meig0 )
mean((pred0$pred[,1]-d0$CMEDV)^2)## mean square prediction error

d0$pred0=pred0$pred[,1]
mapPoints <- ggmap(map) + geom_point(aes(x =LON, y = LAT,color=pred0),data =d0, alpha = 0.7,size=3)+ scale_color_gradient(low="#9ebcda",high="#8856a7")+labs(title="Prediction of owner-occupied homes",y="Latitude",x="Longtitude",color="House Prices (in 1000$'s)" )
mapPoints
#true values
mapPoints <- ggmap(map) + geom_point(aes(x =LON, y = LAT,color=CMEDV),data = d0, alpha = 0.7,size=3)+ scale_color_gradient(low="#9ebcda",high="#8856a7")+labs(title="Owner-occupied homes",y="Latitude",x="Longtitude",color="House Prices (in 1000$'s)" )
mapPoints


# Models with spatially and non-spatially varying coefficients (mod5)
mod5 <- resf_vc(y=y,x=x,xconst=xconst,meig=meig, x_nvc =TRUE)
mod5
meig0 	<- meigen0( meig = meig, coords0 = coords0 )
pred01   <- predict0_vc( mod = mod5, x0 = x0, xconst0=xconst0, meig0 = meig0 )
mean((pred01$pred[,1]-d0$CMEDV)^2)

d0$pred01=pred01$pred[,1]
mapPoints <- ggmap(map) + geom_point(aes(x =LON, y = LAT,color=pred01),data =d0, alpha = 0.7,size=3)+ scale_color_gradient(low="#9ebcda",high="#8856a7")+labs(title="Prediction of owner-occupied homes",y="Latitude",x="Longtitude",color="House Prices (in 1000$'s)" )
mapPoints
#true values
mapPoints <- ggmap(map) + geom_point(aes(x =LON, y = LAT,color=CMEDV),data = d0, alpha = 0.7,size=3)+ scale_color_gradient(low="#9ebcda",high="#8856a7")+labs(title="Owner-occupied homes",y="Latitude",x="Longtitude",color="House Prices (in 1000$'s)" )
mapPoints

# GAM model(mod6)

GAM <- gam(CMEDV ~ CRIM+AGE+ZN+LSTAT+LON+LAT+TAX+PTRATIO+factor(CHAS)+B+RAD, data = boston.c)
summary(GAM)
AIC(GAM)

## Data at unobserved sites , test data
x0 <- d0[,c("CRIM","AGE","ZN","LSTAT","LON","LAT","TAX","PTRATIO","CHAS","B","RAD")]
predgam <- predict(GAM,x0, type="response")
predgam
mean((predgam-d0$CMEDV)^2)

mapPoints <- ggmap(map) + geom_point(aes(x =LON, y = LAT,color=predgam),data =d0, alpha = 0.7,size=3)+ scale_color_gradient(low="#9ebcda",high="#8856a7")+labs(title="Prediction of owner-occupied homes",y="Latitude",x="Longtitude",color="House Prices (in 1000$'s)" )
mapPoints
#true values
mapPoints <- ggmap(map) + geom_point(aes(x =LON, y = LAT,color=CMEDV),data = d0, alpha = 0.7,size=3)+ scale_color_gradient(low="#9ebcda",high="#8856a7")+labs(title="Owner-occupied homes",y="Latitude",x="Longtitude",color="House Prices (in 1000$'s)" )
mapPoints



pred0$pred[1:10,]  # Predicted explained variables
pred0$b_vc[1:10,]  # Predicted SVCs
pred0$bse_vc[1:10,]# Predicted standard errors of the SVCs
pred0$t_vc[1:10,]  # Predicted t-values of the SNVCs
pred0$p_vc[1:10,]  # Predicted p-values of the SNVCs

