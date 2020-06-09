context("Checking that autocart returns a numeric vector.")
library(autocart)

test_that("Autocart returns a numeric vector", {
  # Process the dataset to where it can be used by autocart
  snow <- read.csv(system.file("extdata", "ut2017_snow.csv", package = "autocart", mustWork = TRUE))
  response <- as.matrix(snow$yr50)
  snow <- data.frame(snow$LONGITUDE, snow$LATITUDE, snow$ELEVATION, snow$YRS, snow$HUC,
                snow$TD, snow$FFP, snow$MCMT, snow$MWMT, snow$PPTWT, snow$RH, snow$MAT)
  locations <- as.matrix(cbind(snow$snow.LONGITUDE, snow$snow.LATITUDE))
  alpha <- 0

  # Give all missing values the average of non-missing column values
  for (i in 1:ncol(snow)) {
    snow[is.na(snow[, i]), i] <- mean(snow[, i], na.rm = TRUE)
  }

  predictions <- autocart(response, snow, locations, alpha)

  # TESTING
  expect_is(predictions, "numeric")
})
