## Load mlbench package
library(mlbench)
data("BreastCancer")
?BreastCancer
head(BreastCancer)
str(BreastCancer)

df = BreastCancer[,-1]
# convert factors to numeric
for(i in 1:9) {
  df[, i] <- as.numeric(as.character(df[, i]))
}

str(df)

#install.packages("dlookr")
#install.packages("dplyr")

library(dlookr)
library(dplyr)
diagnose(df) %>%
  filter(missing_count > 0)

breastCancer <- df %>%
  mutate(Bare.nuclei = imputate_na(df, Bare.nuclei, Class,
                                   method = "mice", no_attrs = TRUE, print_flag = FALSE))
breastCancer[c(24,41,140,146,159,165,236,250,276,293,295,298,316,322,412,618),]

?imputate_na
typeof(breastCancer$Class)
head(as.integer(breastCancer$Class))
sf =  data.frame(breastCancer[,-10], y =as.integer(breastCancer$Class)-1)
head(sf)
table(sf$y)

#pairs comparison,corralation 
pairs( sf[,1:9], col=sf[,10]+1)
?pairs

install.packages("corrplot")
library(corrplot)
sf_corr <- cor(sf %>% select(-y))
corrplot::corrplot(sf_corr, order = "hclust", tl.cex = 1, addrect = 8)
?corrplot

#install.packages("psych")
hf_corr <- psych::corr.test(sf[,1:9],)
hf_corr

?gml
sf$y
hf = sf


#install.packages("ggplot2")
install.packages("gridExtra")
library(gridExtra)
library(ggplot2)
breastCancer$Class
nrow(breastCancer)
hf$y <- ifelse(hf$y == 1, "Malignant", "Benign")
head(hf$y)
plot <- list()
box_variables <- c("Cl.thickness", "Cell.shape", "Marg.adhesion", "Bl.cromatin", "Bare.nuclei", "Marg.adhesion", "Mitoses")
for(i in box_variables) {
  plot[[i]] <- ggplot(hf, 
                      aes_string(x = "y", 
                                 y = i, 
                                 col = "y", 
                                 fill = "y")) + 
    geom_boxplot(alpha = 0.2) + 
    theme(legend.position = "none") + 
    scale_color_manual(values = c("blue", "red")) +
    scale_fill_manual(values = c("blue", "red"))
}

do.call(grid.arrange, c(plot, nrow = 1))



##Logistic regression 

(n = nrow(sf))
(p = ncol(sf)-1)

lr_fit = glm(y ~., data =sf , family = "binomial")
summary(lr_fit)
head(sf)


## Apply best subset selection logistic regression
#install.packages("bestglm")
library(bestglm)

bss_fit_AIC = bestglm(sf, family=binomial, IC="AIC")
bss_fit_BIC = bestglm(sf, family=binomial, IC="BIC")
## Examine the results
bss_fit_AIC$Subsets
bss_fit_BIC$Subsets
## Identify best-fitting models
(best_AIC = bss_fit_AIC$ModelReport$Bestk)
(best_BIC = bss_fit_BIC$ModelReport$Bestk)



##Defining functions (K-fold traning set) and random number
## Enter made-up data:
mydata = data.frame(x1=rnorm(100), y=runif(100), x2=rnorm(100))
## Make the response variable the final column:
mydata_2 = mydata[, c(1, 3, 2)]
## Check it seems okay:
head(mydata_2)

## Set the seed (say, at 5) to make the analysis reproducible
set.seed(5)
## Sample the fold-assignment index
nfolds = 10
fold_index = sample(nfolds, n, replace=TRUE)
## Print the first few fold-assignments
head(fold_index)



general_cv = function(X, y, fold_ind, fold_error_function) {
  p = ncol(X)
  Xy = cbind(X, y=y)
  nfolds = max(fold_ind)
  if(!all.equal(sort(unique(fold_ind)), 1:nfolds)) stop("Invalid fold partition.")
  fold_errors = numeric(nfolds)
  # Compute the test error for each fold
  for(fold in 1:nfolds) {
    fold_errors[fold] = fold_error_function(X, y, fold_ind==fold)
  }
  # Find the fold sizes
  fold_sizes = numeric(nfolds)
  for(fold in 1:nfolds) fold_sizes[fold] = length(which(fold_ind==fold))
  # Compute the average test error across folds
  test_error = weighted.mean(fold_errors, w=fold_sizes)
  # Return the test error
  return(test_error)
}



#logistic_reg_fold_error calculate the test error given a matrix of predictor variables
logistic_reg_fold_error = function(X, y, test_data) {
  Xy = data.frame(X, y=y)
  if(ncol(Xy)>1) tmp_fit = glm(y ~ ., data=Xy[!test_data,], family="binomial")
  else tmp_fit = glm(y ~ 1, data=Xy[!test_data,,drop=FALSE], family="binomial")
  phat = predict(tmp_fit, Xy[test_data,,drop=FALSE], type="response")
  yhat = ifelse(phat > 0.5, 1, 0) 
  yobs = y[test_data]
  test_error = 1 - mean(yobs == yhat)
  return(test_error)
}




#Using k-fold Cross-Validation in Best Subset Selection
logistic_reg_bss_cv = function(X, y, fold_ind) {
  p = ncol(X)
  Xy = data.frame(X, y=y)
  X = as.matrix(X)
  nfolds = max(fold_ind)
  if(!all.equal(sort(unique(fold_ind)), 1:nfolds)) stop("Invalid fold partition.")
  fold_errors = matrix(NA, nfolds, p+1) # p+1 because M_0 included in the comparison
  for(fold in 1:nfolds) {
    # Using all *but* the fold as training data, find the best-fitting models 
    # with 0, 1, ..., p predictors, i.e. identify the predictors in M_0, M_1, ..., M_p
    tmp_fit = bestglm(Xy[fold_ind!=fold,], family=binomial, IC="AIC")
    best_models = as.matrix(tmp_fit$Subsets[,2:(1+p)])
    # Using the fold as test data, find the test error associated with each of 
    # M_0, M_1,..., M_p
    for(k in 1:(p+1)) {
      fold_errors[fold, k] = logistic_reg_fold_error(X[,best_models[k,]], y, fold_ind==fold)
    }
  }
  # Find the fold sizes
  fold_sizes = numeric(nfolds)
  for(fold in 1:nfolds) fold_sizes[fold] = length(which(fold_ind==fold))
  # For models with 0, 1, ..., p predictors compute the average test error across folds
  test_errors = numeric(p+1)
  for(k in 1:(p+1)) {
    test_errors[k] = weighted.mean(fold_errors[,k], w=fold_sizes)
  }
  # Return the test error for models with 0, 1, ..., p predictors
  return(test_errors)
}


## Apply the cross-validation for best subset selection function
cv_errors = logistic_reg_bss_cv(sf[,1:p], sf[,p+1], fold_index)
cv_errors
## Identify the number of predictors in the model which minimises test error according to K-fold subset selection
(best_cv = which.min(cv_errors) - 1)




## Create multi-panel plotting device
par(mfrow=c(2, 2))
## Produce plots, highlighting optimal value of k
plot(0:p, bss_fit_AIC$Subsets$AIC, xlab="Number of predictors", ylab="AIC", type="b")
points(best_AIC, bss_fit_AIC$Subsets$AIC[best_AIC+1], col="red", pch=16)
plot(0:p, bss_fit_BIC$Subsets$BIC, xlab="Number of predictors", ylab="BIC", type="b")
points(best_BIC, bss_fit_BIC$Subsets$BIC[best_BIC+1], col="red", pch=16)
plot(0:p, cv_errors, xlab="Number of predictors", ylab="K-fold error", type="b")
points(best_cv, cv_errors[best_cv+1], col="red", pch=16)


#confusion testing error with full data set
lr_fit
phat = predict(lr_fit, sf, type="response")
yhat = as.numeric(ifelse(phat > 0.5, 1, 0))
(confusion = table(Observed=sf$y, Predicted=yhat) )

#Normalising the classification probabilities
normalise = function(x){
  return(x/sum(x))
}
t(apply (confusion, 1, normalise ))

#Training version of estimating of misclassification rate
1 - mean(sf$y == yhat)





##
pstar = 7
## Check which predictors are in the 7-predictor model
bss_fit_AIC$Subsets[pstar+1,]
## Construct a reduced data set containing only the 7 selected predictors
bss_fit_AIC$Subsets[pstar+1, 2:(p+1)]
(indices = which(bss_fit_AIC$Subsets[pstar+1, 2:(p+1)]==TRUE))

sf_red = sf[,c(indices, p+1)]
## Obtain regression coefficients for this model with 7-predictor model
logreg1_fit = glm(y ~ ., data=sf_red, family="binomial")
summary(logreg1_fit)


(test_error = general_cv(sf_red[,1:pstar], sf_red[,pstar+1], fold_index, logistic_reg_fold_error))




## Apply LDA:
install.packages('MASS')
library(MASS)
(lda_fit = lda(y ~., data =sf_red))
is.list(lda_fit)
lda_fit$prior
table(sf$y)/ nrow(sf)



## Compute predicted values:
lda_predict = predict(lda_fit, sf_red)
head(lda_predict$posterior)
yhat = lda_predict$class

## Calculate confusion matrix:
(confusion = table(Observed=sf_red$y, Predicted=yhat))

plot(lda_fit)


lda_fold_error = function(X, y, test_data) {
  Xy = data.frame(X, y=y)
  if(ncol(Xy)>1) tmp_fit = lda(y ~ ., data=Xy[!test_data,])
  tmp_predict = predict(tmp_fit, Xy[test_data,])
  yhat = tmp_predict$class 
  yobs = y[test_data]
  test_error = 1 - mean(yobs == yhat)
  return(test_error)
} 
(test_error = general_cv(sf_red[,1:pstar], sf_red[,pstar+1], fold_index, lda_fold_error))
summary(test_error)



(qda_fit = qda(y ~ ., data=sf_red))

## Compute predicted values:
qda_predict = predict(qda_fit, sf_red)
yhat = qda_predict$y
## Calculate confusion matrix:
(confusion = table(Observed=sf_red$y, Predicted=yhat))



qda_fold_error = function(X, y, test_data) {
  Xy = data.frame(X, y=y)
  if(ncol(Xy)>1) tmp_fit = qda(y ~ ., data=Xy[!test_data,])
  tmp_predict = predict(tmp_fit, Xy[test_data,])
  yhat = tmp_predict$class 
  yobs = y[test_data]
  test_error = 1 - mean(yobs == yhat)
  return(test_error)
}

(test_error = general_cv(sf_red[,1:pstar], sf_red[,pstar+1], fold_index, qda_fold_error))







