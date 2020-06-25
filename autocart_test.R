# Clustering
cluster_stations <- function(lon, lat, dist_adj, h,
                             method = "complete") {
  if (length(lon) != length(lat)) {
    stop("Lengths of lon, lat, and elev must be the same.")
  }
  if (!is.numeric(lon) || !is.numeric(lat)) {
    stop("Vectors lon, lat, and elev must be numeric.")
  }
  if (any(is.na(lon)) || any(is.na(lat))) {
    stop("All elev, lon, and lat not have NA values.")
  }

  # Find similarity matrix
  #=============================================================================
  # find distances between all locations in km
  dist <- fields::rdist.earth(matrix(c(lon, lat), ncol = 2), miles = FALSE)

  # find difference in elevation between all locations
  # delev <- abs(outer(elev, elev, "-"))

  # Produce "similarity matrix":
  # - distance of dist_adj km produces a score of 1
  # - elevation difference of elev_adj produces a score of 1
  similarity <- stats::as.dist((dist/dist_adj))

  # Cluster and return cluster vector
  #=============================================================================
  # Cluster all stations with the farthest neighbors in a cluster
  if (length(lon) > 1) {
    clust <- stats::hclust(similarity, method)
  } else if (length(lon) == 1) {
    return(1)
  } else return(NULL)

  # Cut heirarchical clustering tree at h and return cluster id's
  stats::cutree(clust, h = h)
}






# Cluster stations for better testing
xvs <- cluster_stations(snow$LONGITUDE, snow$LATITUDE, 60, 5)






# BEFORE RUNNING:
# Make sure you have autocart installed. Autocart can be installed with the command
# devtools::install_github("ethanancell/autocart")

library(autocart)
library(rpart)

# Import the processed data
snow <- read.csv(system.file("extdata", "ut2017_snow.csv", package = "autocart", mustWork = TRUE))

# ============================
# === EVALUATION FUNCTIONS ===
# ============================

# The functions contained here are to evaluate the results given by autocart and rpart

# Relative mean absolute error
rmae <- function(pred, obs, na.rm = TRUE) {
  if (length(pred) != length(obs)) {
    stop("Predicted and observed vectors must be same length.")
  }
  mean(abs((pred-obs)/mean(obs))*100, na.rm = na.rm)
}

# Root mean square error
rmse <- function(pred, obs, na.rm = TRUE) {
  if (length(pred) != length(obs)) {
    stop("Predicted and observed vectors must be same length.")
  }
  sqrt(mean((pred-obs)^2, na.rm = na.rm))
}

# Output both rmse and rmae in one easy paste function
test_predictions <- function(pred, obs, na.rm = TRUE) {
  if (length(pred) != length(obs)) {
    stop("Predicted and observed vectors must be same length.")
  }

  print(paste("RMAE is ", rmae(pred, obs, na.rm), ", and RMSE is ", rmse(pred, obs, na.rm), sep=""))
}

# ==========================
# ======== AUTOCART ========
# ==========================

# Process dataset to make it suitable for autocart
ac_response <- as.matrix(snow$yr50)
ac_snow <- data.frame(snow$LONGITUDE, snow$LATITUDE, snow$ELEVATION, snow$YRS, snow$HUC,
                      snow$TD, snow$FFP, snow$MCMT, snow$MWMT, snow$PPTWT, snow$RH, snow$MAT)
ac_locations <- as.matrix(cbind(snow$LONGITUDE, snow$LATITUDE))

# Give all missing values the average of non-missing column values
for (i in 1:ncol(ac_snow)) {
  ac_snow[is.na(ac_snow[, i]), i] <- mean(ac_snow[, i], na.rm = TRUE)
}

# 10-fold cross-validation with alpha = 0.85
alpha <- 0
beta <- 0

#xvs <- rep(1:10, length = nrow(ac_snow))
#xvs <- sample(xvs)
xvs <- cluster_stations(snow$LONGITUDE, snow$LATITUDE, 60, 5) # Spatial blocking
ac_pred <- rep(NA, length = nrow(ac_snow))
for (k in 1:length(unique(xvs))) {
  train <- ac_snow[xvs != k, ]
  test <- ac_snow[xvs == k, ]

  train_y <- ac_response[xvs != k]
  test_y <- ac_response[xvs == k]

  train_loc <- ac_locations[xvs != k, ]
  test_loc <- ac_locations[xvs == k, ]

  # Create the model based off the test dataset
  trained_model <- autocart(train_y, train, train_loc, alpha, beta, autocartControl(distpower=2))

  ac_pred[xvs == k] <- predictAutocart(trained_model, test)
}

# ======================================
# ====== AUTOCART + SPATIAL NODES ======
# ======================================

# Process dataset to make it suitable for autocart
acsn_response <- as.matrix(snow$yr50)
acsn_snow <- data.frame(snow$LONGITUDE, snow$LATITUDE, snow$ELEVATION, snow$YRS, snow$HUC,
                        snow$TD, snow$FFP, snow$MCMT, snow$MWMT, snow$PPTWT, snow$RH, snow$MAT)
acsn_locations <- as.matrix(cbind(snow$LONGITUDE, snow$LATITUDE))

# Give all missing values the average of non-missing column values
for (i in 1:ncol(acsn_snow)) {
  acsn_snow[is.na(acsn_snow[, i]), i] <- mean(acsn_snow[, i], na.rm = TRUE)
}

# 10-fold cross-validation with alpha = 0.85
alpha <- 0
beta <- 0

#xvs <- rep(1:10, length = nrow(acsn_snow))
#xvs <- sample(xvs)
xvs <- cluster_stations(snow$LONGITUDE, snow$LATITUDE, 60, 5) # Spatial blocking
acsn_pred <- rep(NA, length = nrow(acsn_snow))
for (k in 1:length(unique(xvs))) {
  train <- acsn_snow[xvs != k, ]
  test <- acsn_snow[xvs == k, ]

  train_y <- acsn_response[xvs != k]
  test_y <- acsn_response[xvs == k]

  train_loc <- acsn_locations[xvs != k, ]
  test_loc <- acsn_locations[xvs == k, ]

  # Create the model based off the test dataset
  trained_model <- autocart(train_y, train, train_loc, alpha, beta, autocartControl(distpower=2))

  acsn_pred[xvs == k] <- predictSpatialNodes(spatialNodes(trained_model), trained_model, test, test_loc)
}

# =========================
# ========- RPART =========
# =========================

# Process dataset to make it suitable for autocart
rp_snow <- data.frame(snow$yr50, snow$LONGITUDE, snow$LATITUDE, snow$ELEVATION, snow$YRS, snow$HUC,
                      snow$TD, snow$FFP, snow$MCMT, snow$MWMT, snow$PPTWT, snow$RH, snow$MAT)

# Give all missing values the average of non-missing column values
for (i in 1:ncol(rp_snow)) {
  rp_snow[is.na(rp_snow[, i]), i] <- mean(rp_snow[, i], na.rm = TRUE)
}

# 10-fold cross-validation with alpha = 0.85
#xvs <- rep(1:10, length = nrow(rp_snow))
#xvs <- sample
xvs <- cluster_stations(snow$LONGITUDE, snow$LATITUDE, 60, 5) # Spatial blocking
rp_pred <- rep(NA, length = nrow(rp_snow))
for (k in 1:length(unique(xvs))) {
  train <- rp_snow[xvs != k, ]
  test <- rp_snow[xvs == k, ]

  # Create the model based off the test dataset
  trained_model <- rpart(snow.yr50 ~ ., data=train, cp=0)

  rp_pred[xvs == k] <- predict(trained_model, test)
}

# =================
# ==== RESULTS ====
# =================

# Results
print("Testing autocart:")
test_predictions(ac_pred, ac_response)
print("Testing autocart + spatial nodes:")
test_predictions(acsn_pred, acsn_response)
print("Testing rpart:")
test_predictions(rp_pred, ac_response)
