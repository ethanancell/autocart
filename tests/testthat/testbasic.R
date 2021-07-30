context("Checking that autocart successfully creates a model and can predict with the model.")
library(autocart)

test_that("Autocart returns sensical output", {
  # Process the dataset so it can be used by autocart
  snow <- read.csv(system.file("extdata", "ut2017_snow.csv", package = "autocart", mustWork = TRUE))

  # The snow dataset does not contain factors. Here are some conjured factors that correlate with the
  # response so that the test ensures that factors can be used.
  factor_column <- rep(NA, nrow(snow))
  factor_column[snow$yr50 <= quantile(snow$yr50)[2]] <- "aa"
  factor_column[snow$yr50 > quantile(snow$yr50)[2] & snow$yr50 <= quantile(snow$yr50)[3]] <- "bb"
  factor_column[snow$yr50 > quantile(snow$yr50)[3]] <- "cc"
  factor_column <- as.factor(factor_column)
  snow <- cbind(snow, factor_column)
  snow <- na.omit(snow)

  # Cut data in half to make a quick run
  snow <- snow[1:(round(nrow(snow)) / 2), ]

  snow <- SpatialPointsDataFrame(coords = snow[, c("LONGITUDE", "LATITUDE")],
                                 data = snow,
                                 proj4string = CRS("+proj=longlat"))
  alpha <- 0.80
  beta <- 0.10

  # Set a spatial bandwidth
  myControl <- autocartControl(spatialBandwidthProportion = 0.8)

  model <- autocart(yr50 ~ ELEVATION + YRS + TD + MCMT + MWMT + factor_column,
                    alpha = alpha, beta = beta, data = snow, control = myControl)

   # TESTING
  expect_is(model$prediction, "numeric")
  expect_equal(length(model$prediction), length(snow$yr50))
  expect_equal(model$prediction, predict(model, newdata = snow))
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
    trained_model <- autocart(train_y, train, train_loc, 0, 0)

    pred[xvs == k] <- predictAutocart(trained_model, test)
  }

  # Testing is successful if there are no NA values to be found in pred
  expect_equal(any(is.na(pred)), FALSE)
})

test_that("autocartControl controls the tree correctly", {
  # Process the dataset so it can be used by autocart
  snow <- read.csv(system.file("extdata", "ut2017_snow.csv", package = "autocart", mustWork = TRUE))

  # The snow dataset does not contain factors. Here are some conjured factors that correlate with the
  # response so that the test ensures that factors can be used.
  factor_column <- rep(NA, nrow(snow))
  factor_column[snow$yr50 <= quantile(snow$yr50)[2]] <- "aa"
  factor_column[snow$yr50 > quantile(snow$yr50)[2] & snow$yr50 <= quantile(snow$yr50)[3]] <- "bb"
  factor_column[snow$yr50 > quantile(snow$yr50)[3]] <- "cc"
  factor_column <- as.factor(factor_column)
  snow <- cbind(snow, factor_column)
  snow <- na.omit(snow)

  # Cut data in half to make a quick run
  snow <- snow[1:(round(nrow(snow)) / 2), ]

  snow <- SpatialPointsDataFrame(coords = snow[, c("LONGITUDE", "LATITUDE")],
                                 data = snow,
                                 proj4string = CRS("+proj=longlat"))
  alpha <- 0.80
  beta <- 0.10

  # Arbitrary control parameters we will test against
  ms <- 35
  mb <- 9
  md <- 5
  dpower <- 2

  my_control <- autocartControl(minsplit = ms, minbucket = mb, maxdepth = md, distpower = dpower)
  model <- autocart(yr50 ~ ELEVATION + YRS + TD + MCMT + MWMT + factor_column,
                    alpha = alpha, beta = beta, data = snow, control = my_control)

  # TESTING
  splitframe <- model$splitframe
  for (i in nrow(splitframe)) {
    # For all the nodes which are not terminal, they should not have numObs < ms.
    # For all terminal nodes, they shouldn't have numobs < mb
    if (!splitframe$isterminal[i]) {
      expect_gt(splitframe$numobs[i], ms-1)
    } else {
      expect_gt(splitframe$numobs[i], mb-1)
    }
  }
})

