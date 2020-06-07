library(testthat)
library(autocart)
library(tidyverse)

test_check("autocart")

# Data step
data <- read_csv(system.file("extdata", "ut2017_snow.csv", package = "autocart", mustWork = TRUE))
snow <- read_csv(system.file("extdata", "ut2017_snow.csv", package = "autocart", mustWork = TRUE))
response <- snow %>%
  select(yr50)
response <- as.vector(response)
snow <- snow %>%
  select(LONGITUDE, LATITUDE, ELEVATION, YRS, HUC, TD, FFP, MCMT,
         MWMT, PPTWT, RH, MAT)
locations <- snow %>%
  select(LONGITUDE, LATITUDE)
locations <- as.matrix(locations)

alpha = 0.5

# For the pure sake of testing, let's give all the missing values in snow just the average of the
# non-missing column values.
for (i in 1:ncol(snow)) {
  snow[is.na(snow[, i]), i] <- mean(snow[[i]], na.rm = TRUE)
}

# Autocart model
predictions <- autocart(as.matrix(response), snow, locations, alpha)
subsetRows(snow, c(2,3))
testSplit(as.matrix(response), as.matrix(snow$ELEVATION), locations, alpha)

# Test creating the tree
test_tree()

# Rpart model
library(rpart)
library(rpart.plot)
data <- data %>%
  select(yr50, LONGITUDE, LATITUDE, ELEVATION, YRS, HUC, TD, FFP, MCMT,
         MWMT, PPTWT, RH, MAT)
control_model <- rpart(yr50 ~ ., data=data)
rpart.plot(control_model)
