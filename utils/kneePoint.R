# https://github.com/CeMOS-Mannheim/moleculaR/tree/main/R
# https://rdrr.io/github/CeMOS-Mannheim/moleculaR/man/kneePoint.html
# https://towardsdatascience.com/detecting-knee-elbow-points-in-a-graph-d13fc517a63c
# This function calculates the knee/elbow point of a curve based on the kneedle algorithm (satopaa et al, 2011). 
# This is a simplified implementation.


kneePoint     <- function(x, y, df = 7,
                          xQuery = seq(range(x)[1], range(x)[2], 0.1),
                          plot = FALSE, sign = +1, ...) {
  
  
  
  
  # fit a smoothing spline/loess
  smoothx   <- xQuery
  smoothspl <- smooth.spline(x = x, y = y, df = df)
  smoothy <- predict(smoothspl, x = smoothx)$y
  
  
  # normalize points of the smoothing curve to unit square
  smoothnx     <- (smoothx - min(smoothx)) / (max(smoothx) - min(smoothx))
  smoothny     <- (smoothy - min(smoothy)) / (max(smoothy) - min(smoothy))
  
  
  # apply kneedle
  k <- .kneedle(x = smoothnx, y = smoothny, sign = sign)
  
  
  if(plot){
    par(mar=c(5.1, 5.1, 4.1, 2.1))
    plot(x = smoothx, y = smoothy, type = "l",
         main = "Knee-point estimation",
         xlab = "Gaussian Bandwidth",
         ylab = "Moran's I",
         ...)
    
    # defaults
    lwd <- 1; cex <- 1
    
    # extract size arg from ...
    pargs <- list(...)
    
    if(length(pargs) > 0){
      n <- names(pargs)
      
      if("lwd" %in% n){
        lwd <- pargs$lwd
      }
      
      
      if("cex.lab" %in% n){
        cex <- pargs$cex.lab
      }
    }
    
    abline(v = smoothx[k], col = "chocolate1", lty ="dashed", lwd = lwd)
    
    legend("right", bty = "n",
           legend = c(paste0("Moran's I"),
                      paste0("Knee pnt=", round(smoothx[k], 2))
           ),
           col = c("black",
                   "chocolate1"
           ),
           lty = c("solid",
                   "dashed"
           ),
           lwd = lwd, cex = cex)
  }
  
  return(smoothx[k])
  
  
}


#' Find the knee/elbow point in a vector using the Kneedle algorithm.
#'
#' This function uses the Kneedle algorithm to find the index of the knee point
#' in the provided x,y-vectors.
#'
#' @param x numeric vector, x-values.
#' @param y numeric vector, y-values.
#' @param sign +1 for increasing values and  -1 for decreasing values.
#'
#' @return The index of the knee/elbow.
#'
.kneedle <- function(x, y, sign = 1) {
  
  if(length(x) != length(y))
    stop("error in internal function .kneedle; x and y of different lengths.\n")
  
  start = c(x[1], y[1])
  end = c(x[length(x)], y[length(y)])
  
  k <- which.max(lapply(1:length(x), function(i) {
    sign * -1 * .dist2d(c(x[i], y[i]),
                        start,
                        end)
  }))
  
  
  k
  
}

.dist2d <- function(a,b,c) {
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- det(m)/sqrt(sum(v1*v1))
  d
}
