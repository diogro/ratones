#' Geometric mean
#'
#' Calculates the geometric mean of a vector
#'
#' @param x input vector
#' @param na.rm logical. If true, NA are removed before calculation
#' @param zero.propagate logical. If TRUE, checks the vector for zeros and skips calculation
#' @return returns the geometric mean
#' @export
#' @examples
#' x = runif(10, 1, 5)
#' gm_mean(x)
#' gm_mean(c(1, 2, NA, 0), na.rm = FALSE)
#' gm_mean(c(1, 2, NA, 0), na.rm = FALSE, zero.propagate = TRUE)
gm_mean = function(x, na.rm = TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}
