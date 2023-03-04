# Jasjeet S. Sekhon
# Harvard University
# http://jsekhon.fas.harvard.edu
# Feb, 2005
#
# simple utility functions

get.xdata <- function(formul, datafr) {
  t1 <- terms(formul, data=datafr);
  if (length(attr(t1, "term.labels"))==0 & attr(t1, "intercept")==0) {
    m <- NULL;  # no regressors specified for the model matrix
  } else {
    m <- model.matrix(formul, data=datafr, drop.unused.levels = TRUE);
  }
  return(m);
}


# get.ydata:
# Return response vector corresponding to the formula in formul
# 
get.ydata <- function(formul, datafr) {
  t1 <- terms(formul, data=datafr);
  if (length(attr(t1, "response"))==0) {
    m <- NULL;  # no response variable specified
  }  else {
    m <- model.response(model.frame(formul, data=datafr));
  }
  return(m);
}
