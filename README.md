# autocart
A modified regression tree R package so that spatial autocorrelation can be tracked at every step of the splitting.

Existing classification and regression tree (CART) R packages allow for custom splitting functions, but not at high enough of a level to
where geographical coordinates can be tracked at every stage to include custom splitting that accounts for spatial autocorrelation.

This package is an implementation of the CART algorithm that allows for tracking geographical location in the observations at every
stage of the tree creation.
