#-------------------------------------------------------------------------------
# EXAMPLE: CANCER DATA
#-------------------------------------------------------------------------------

# Load the data

load(file="C:\\Users\\bapti\\Downloads\\BEC_1ère_année\\Autumn_2024\\Advanced_data_analysis\\KRAS_dep.RDATA")
windows()
par(mfrow=c(2,2))

data[1:10,1:10]
y=data$KRAS_dep
X=as.matrix(data[,-1])
ncol(data[, -1])

### Do the filtering of irrelevant / sparse features 

# Calculate the proportion of 1s (mutated samples) for each gene
mutation_proportions <- colSums(data[,-1] == 1) / nrow(data[,-1])

# Filter genes where mutation proportion is >= 1%
selected_genes <- mutation_proportions[mutation_proportions >= 0.01]
length(as.vector(selected_genes))

# Subset the original data to only include selected genes
data_filtered <- data[, c("KRAS_dep", names(selected_genes))]

# Now you can proceed with the rest of your analysis
# For example, training a linear regression model on the filtered data:

X_filtered <- as.matrix(data_filtered[, -1])  # Remove KRAS_dep column
ncol(data_filtered[, -1])
y_filtered <- data_filtered$KRAS_dep
length(y_filtered)
###

hist(data$KRAS_dep,breaks = 30,main = "KRAS dependency scores",xlab = "Score",col="#FDCF76")

#-------------------------------------------------------------------------------
# Linear regression using (not) a subset of the observation (so that n<p or n~p)
#-------------------------------------------------------------------------------

lr=lm(KRAS_dep~.,data=data)
summary(lr)

names(summary(lr))
p_values <- summary(lr)$coefficients[, "Pr(>|t|)"]
range(p_values)
length(as.vector(p_values))
hist(p_values,breaks = 100,main = "p-values",xlab = "Score",col="#FDCF76")
# --> 74 predictors 
# Plot a density plot of the p-values
plot(density(p_values), 
     main = "Density Plot of P-values", 
     xlab = "P-value", 
     col = "red", 
     lwd = 2)

adjusted_pvals <- p.adjust(p_values, method = "fdr")
hist(adjusted_pvals,breaks = 10,main = "p-values",xlab = "Score",col="#FDCF76")

#-------------------------------------------------------------------------------
# Training and Test set --- TEST error
#-------------------------------------------------------------------------------

n=nrow(data_filtered)
n.train=floor(0.75*n) # 75% of the observations
set.seed(123)
train=sample(1:n, size = n.train, replace = FALSE)
X.train=X_filtered[train,]
y.train=y_filtered[train]
X.test=X_filtered[-train,]
y.test=y_filtered[-train]

nrow(X_filtered[-train,])

#-------------------------------------------------------------------------------
# Standard linear regression (OLS estimates)
#-------------------------------------------------------------------------------

lr_filtered=lm(KRAS_dep~.,data=data_filtered[train,])
summary(lr_filtered)
p_values <- summary(lr_filtered)$coefficients[, "Pr(>|t|)"]
range(p_values)
length(as.vector(p_values))
##
library(car)
range(vif(lr_filtered))
alias(lr_filtered)  # Identifie les coefficients aliasés
##
n.var.lr=length(which(summary(lr_filtered)$coeff[-1,1]!=0)) # number of variables (excluding the intercept) different from 0 (for summary table)

lr.pred=predict(lr_filtered,newdata = data_filtered[-train,])
length(as.vector(lr.pred))
test_error_ols=mean((lr.pred-y.test)^2)

#-------------------------------------------------------------------------------
# RIDGE regression 
#-------------------------------------------------------------------------------
# install.packages("glmnet")
library(glmnet)

# 1) Finding the optimal lambda
set.seed(123)
ridge.cv=cv.glmnet(x = X.train, y = y.train, family = "gaussian", type.measure = "mse", alpha = 0)
# alpha == 0 means ridge regression.
plot(ridge.cv)
bestlam=ridge.cv$lambda.min

# 2) Fit the model 
ridge.fit=glmnet(x = X.train, y = y.train, family = "gaussian",type.measure = "mse", alpha = 0)
coef(ridge.fit,s=bestlam)
coefs_df <- as.data.frame(as.matrix(coef(ridge.fit,s=bestlam)))
# Add gene names (rownames) as a column
coefs_df$gene <- rownames(coefs_df)

# Filter for coefficients with absolute value greater than 0.1 (excluding intercept)
significant_coefs <- coefs_df[abs(coefs_df$s1) > 0.1 & rownames(coefs_df) != "(Intercept)", ]

# Order by coefficient values (descending)
significant_coefs <- significant_coefs[order(abs(significant_coefs$s1), decreasing = TRUE), ]

# Display ordered significant coefficients
print(significant_coefs)


length(coef(ridge.fit,s=bestlam))
# predict(ridge.fit,type = "coefficients",s=bestlam)# same output of line above
n.var.ridge=length(which(coef(ridge.fit,s=bestlam)[-1]!=0)) # number of variables (excluding the intercept) different from 0 (for summary table)
# plot the coefficient path
# Initialize annotations with predictor names (excluding intercept)
annotations <- rownames(coef(ridge.fit, s = bestlam))[-1]
# Set predictors with absolute coefficients <= 0.1 to an empty string
annotations[abs(coef(ridge.fit, s = bestlam)[-1]) <= 0.1] <- ""

plot(ridge.fit, xvar="lambda",xlim=c(-5,5))
abline(v=log(bestlam), lty=2,col="grey")
text(-5.5, coef(ridge.fit)[-1,length(ridge.cv$lambda)], annotations, pos=4,cex=0.6)

# 3) Predictions
ridge.pred=predict(ridge.fit, s=bestlam, newx = X.test)
test_error_ridge=mean((ridge.pred-y.test)^2)

#-------------------------------------------------------------------------------
# LASSO regression 
#-------------------------------------------------------------------------------

set.seed(123)
# 1) Finding the optimal lambda
lasso.cv=cv.glmnet(X.train,y.train,family = "gaussian", type.measure = "mse",alpha = 1)
plot(lasso.cv)
bestlam.min=lasso.cv$lambda.min
bestlam.1se=lasso.cv$lambda.1se

# 2) Fit the model 
lasso.fit=glmnet(X.train,y.train,alpha = 1)
coef(lasso.fit,s=bestlam.1se)[which(coef(lasso.fit,s=bestlam.1se)!=0),,drop=FALSE]
coef(lasso.fit,s=bestlam.min)[which(coef(lasso.fit,s=bestlam.min)!=0),,drop=FALSE]
n.var.lasso=length(which(coef(lasso.fit,s=bestlam.min)[-1]!=0)) # number of variables (excluding the intercept) different from 0 (for summary table)

# plot the coefficient path
labels <- rownames(coef(lasso.fit))[-1]
labels=ifelse(coef(lasso.fit,s=bestlam.min)[-1]==0,"",labels)
plot(lasso.fit, xvar="lambda",xlim=c(-10.3,-2))
abline(v=c(log(bestlam.min),log(bestlam.1se)), lty=2,col="grey")
# Add text to the plot
text(rep(-10.5, length(labels)), 
     coef(lasso.fit)[-1, length(lasso.fit$lambda)], 
     labels, 
     pos = 4, cex = 0.6)
# 3) Predictions
lasso.pred=predict(lasso.fit, s=bestlam.min, newx = X.test)
test_error_lasso=mean((lasso.pred-y.test)^2)

#-------------------------------------------------------------------------------
# ELASTIC NET regression 
#-------------------------------------------------------------------------------

# 1) Finding the optimal alpha and lambda
alpha.grid=seq(0.1,0.9,0.1)
list.of.fits=list()
n.a=1
for (a in alpha.grid){
  set.seed(123)
  list.of.fits[[n.a]]=cv.glmnet(X.train, y.train, family = "gaussian", type.measure = "mse", alpha=a)
  names(list.of.fits)[n.a]=paste("alpha",a)
  n.a=n.a+1
}

min.cvm=c()
for (i in 1:length(list.of.fits)){
  min.cvm[i]=min(list.of.fits[[i]]$cvm)
}

bestalpha=alpha.grid[which.min(min.cvm)]
bestlam.min=list.of.fits[[which.min(min.cvm)]]$lambda.min
bestlam.1se=list.of.fits[[which.min(min.cvm)]]$lambda.1se

# 2) Fit the model 
en.fit=glmnet(X.train,y.train,family = "gaussian", type.measure = "mse",alpha=bestalpha)

coef(en.fit,s=bestlam.min)[which(coef(en.fit,s=bestlam.min)!=0),]
n.var.en=length(which(coef(en.fit,s=bestlam.min)[-1]!=0)) # number of variables (excluding the intercept) different from 0 (for summary table)

# plot the coefficient path
labels <- rownames(coef(lasso.fit))[-1]
labels=ifelse(coef(en.fit,s=bestlam.min)[-1]==0,"",labels)
plot(en.fit, xvar="lambda",xlim=c(-10.5,-2))
abline(v=c(log(bestlam.min),log(bestlam.1se)), lty=2,col="grey")
# Extract predictor labels (excluding intercept)
labels <- rownames(coef(en.fit))[-1]

# Set labels to an empty string for predictors with small coefficients
labels[abs(coef(en.fit)[-1, length(en.fit$lambda)]) <= 0.1] <- ""

# Add text to the plot
text(rep(-10.8, length(labels)), 
     coef(en.fit)[-1, length(en.fit$lambda)], 
     labels, 
     pos = 4, cex = 0.6)



# 3) Predictions
en.pred=predict(en.fit, s=bestlam.min, newx = X.test)
test_error_en=mean((en.pred-y.test)^2)

#-------------------------------------------------------------------------------
### Summary table - models comparison
#-------------------------------------------------------------------------------

mse=c(test_error_ols,test_error_ridge,test_error_lasso,test_error_en)
n.var=c(n.var.lr,n.var.ridge,n.var.lasso,n.var.en)

table=cbind(mse,n.var)
rownames(table)=c("OLS","ridge","lasso","elastic net")






