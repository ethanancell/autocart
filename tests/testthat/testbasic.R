context("Checking that autocart successfully creates a model and can predict with the model.")
library(autocart)

test_that("Autocart returns sensical output", {
  # Process the dataset to where it can be used by autocart
  snow <- read.csv(system.file("extdata", "ut2017_snow.csv", package = "autocart", mustWork = TRUE))
  response <- as.matrix(snow$yr50)
  snow <- data.frame(snow$LONGITUDE, snow$LATITUDE, snow$ELEVATION, snow$YRS, snow$HUC,
                snow$TD, snow$FFP, snow$MCMT, snow$MWMT, snow$PPTWT, snow$RH, snow$MAT)
  locations <- as.matrix(cbind(snow$snow.LONGITUDE, snow$snow.LATITUDE))
  alpha <- 0.85

  # Give all missing values the average of non-missing column values
  for (i in 1:ncol(snow)) {
    snow[is.na(snow[, i]), i] <- mean(snow[, i], na.rm = TRUE)
  }

  model <- autocart(response, snow, locations, alpha)

  # TESTING
  expect_is(model$prediction, "numeric")
  expect_equal(length(model$prediction), length(response))
  expect_equal(model$prediction, predictAutocart(model, snow))
})

test_that("Cross-validation with autocart is possible", {
  # Process the dataset
  snow <- read.csv(system.file("extdata", "ut2017_snow.csv", package = "autocart", mustWork = TRUE))
  response <- as.matrix(snow$yr50)
  snow <- data.frame(snow$LONGITUDE, snow$LATITUDE, snow$ELEVATION, snow$YRS, snow$HUC,
                     snow$TD, snow$FFP, snow$MCMT, snow$MWMT, snow$PPTWT, snow$RH, snow$MAT)
  locations <- as.matrix(cbind(snow$snow.LONGITUDE, snow$snow.LATITUDE))

  # Give all missing values the average of non-missing column values
  for (i in 1:ncol(snow)) {
    snow[is.na(snow[, i]), i] <- mean(snow[, i], na.rm = TRUE)
  }

  # 10-fold cross-validation. Use alpha=0 for speed of computation.
  xvs <- rep(1:10, length = nrow(snow))
  xvs <- sample(xvs)
  pred <- rep(NA, length = nrow(snow))
  for (k in 1:10) {
    train <- snow[xvs != k, ]
    test <- snow[xvs == k, ]

    train_y <- response[xvs != k]
    test_y <- response[xvs == k]

    train_loc <- locations[xvs != k, ]
    test_loc <- locations[xvs == k, ]

    # Create the model based off the test dataset
    trained_model <- autocart(train_y, train, train_loc, 0)

    pred[xvs == k] <- predictAutocart(trained_model, test)
  }

  # Testing is successful if there are no NA values to be found in pred
  expect_equal(any(is.na(pred)), FALSE)
})
