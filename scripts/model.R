


rm(list=ls())

setwd('/projects/covid-ct/imlab/users/temi/projects/TFXcan/scripts/')

library(glue)
library(GenomicRanges)
library(reticulate)
library(R.utils)
library(data.table)
library(glmnet)
library(bit64)

#as.integer64(.Machine$integer.max) + 1L


# I very much prefer to extract the bed files this way because of the way the files have been archived i.e. zip:gzip...
# and using the python code I wrote makes it easy - also a good way to learn how to use reticulate ^^
# reticulate::use_python("~/miniconda3/bin/python", required = T)
# reticulate::py_config()
#reticulate::source_python('./load-cistrome-data.py')

# np <- import("numpy")
# # data reading
# X_train <- np$load("../output/train-test-data/X_train.npy")
# X_test <- np$load("../output/train-test-data/X_test.npy")
# y_train <- np$load("../output/train-test-data/y_train.npy")
# y_test <- np$load("../output/train-test-data/y_test.npy")

dt <- read.csv('/projects/covid-ct/imlab/users/temi/projects/TFXcan/output/train-test-data/data.csv')

X_train <- as.matrix(dt[, -c(1, 2)])
y_train <- as.matrix(dt[, 2])

ENet_fit <- cv.glmnet(x=X_train[complete.cases(X_train),], y=y_train[complete.cases(X_train)], family = "binomial", type.measure = "auc", alpha = 0.5) #alpha: mixing term, lasso, 1-alpha ridge.

plot(ENet_fit)

View(X_train)
