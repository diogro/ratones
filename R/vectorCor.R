#' Vector correlation
#'
#' Calculates the cosine correlation between two vectors
#' @param x numeric vector
#' @param y numeric vector
#' @importFrom evolqg Normalize
#' @export
#' @examples
#' x = rnorm(10)
#' y = rnorm(10)
#' vectorCor(x, y)
vectorCor <- function(x, y) t(Normalize(x)) %*% Normalize(y)
