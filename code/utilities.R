
#' Convert sfc_POINT to x,y columns
#'
#' @param x
#' @param geometry
#' @param names
#'
#' @return
#' @export
#'
#' @examples
#'
sfc_as_cols <- function(x, geometry, names = c("x","y")) {
  ## from https://github.com/r-spatial/sf/issues/231#issuecomment-424929353
  if (missing(geometry)) {
    geometry <- sf::st_geometry(x)
  } else {
    geometry <- rlang::eval_tidy(enquo(geometry), x)
  }
  stopifnot(inherits(x,"sf") && inherits(geometry,"sfc_POINT"))
  ret <- sf::st_coordinates(geometry)
  ret <- tibble::as_tibble(ret)
  stopifnot(length(names) == ncol(ret))
  x <- x[ , !names(x) %in% names]
  ret <- setNames(ret,names)
  dplyr::bind_cols(x,ret)
}


#' @title Coefficient of Variation
#' @name CV
#' @description This function generates the CV and associated equi-tailed two-sided nonparametric
#'              confidence interval.
#' @param x numeric vector
#' @param R integer indicating the number of bootstrap replicates
#' @param type A vector of character strings representing the type
#'             of intervals required. The value should be any subset
#'             of the values c("norm","basic", "stud", "perc", "bca").
#' @param seed A single integer to identify the random number generator state.
#' @param conf A scalar or vector containing the confidence level(s) of the required interval(s).
#' @param na.rm a logical value indicating whether \code{NA} values should be
#'              removed before the computation proceeds.
#'
#' @return \describe{an object of type "list" which containes the estimate (est) the lower and
#'         upper intervals (lower and upper), and the bootstrap type.} \item
#' @export
#'
#' @examples
#'
#' x <- rnorm(n = 20, mean = 0, sd = 1)
#' CV(x, seed = 10)
#'
CV <- function(x,
               R = 1000,
               type = "basic",
               seed = NULL,
               conf = 0.95, na.rm = TRUE){

  if(!is.null(seed)){
    set.seed(seed)
  }

  sample_cv <- function(x, i){
    return(stats::sd(x[i], na.rm = TRUE)/mean(x[i], na.rm = TRUE))
  }

  boot_cv <- boot::boot(
    data = x,
    function(x, i) {  # coefficient of variation
      stats::sd(x[i], na.rm = na.rm)/mean(x[i], na.rm = na.rm)
    },
    R = R
  )

  boot_cv_ci <- boot::boot.ci(boot.out = boot_cv,
                           conf = conf,
                           type = type)

  lower <- ifelse(type == "normal",
                  boot_cv_ci[[type]][2],
                  boot_cv_ci[[type]][4])
  upper <- ifelse(type == "normal",
                  boot_cv_ci[[type]][3],
                  boot_cv_ci[[type]][5])

  return(list(est =  boot_cv_ci$t0,
              lower = lower,
              upper = upper,
              type = type))
}

#' Proportional Variability Index
#'
#' Calculate PV of the variable x.
#' @param x
#'
#' @return a numeric value
#' @export
#'
#' @examples
#' x <- c(1, 5, 10, 4, 3)
#' PV(x)
PV <- function (x,
                R = 1000,
                type = "basic",
                seed = NULL,
                conf = 0.95){

  if(!is.null(seed)){
    set.seed(seed)
  }

  sample_pv <- function(x, i){
    return(pv_test(x[i]))
  }

  boot_pv <- boot::boot(
    data = x,
    statistic = sample_pv,
    R = R
  )

  boot_pv_ci <- boot::boot.ci(boot.out = boot_pv,
                              conf = conf,
                              type = type)

  lower <- ifelse(type == "normal",
                  boot_pv_ci[[type]][2],
                  boot_pv_ci[[type]][4])
  upper <- ifelse(type == "normal",
                  boot_pv_ci[[type]][3],
                  boot_pv_ci[[type]][5])

  return(list(est = boot_pv_ci$t0,
              lower = lower,
              upper = upper,
              type = type))
}

cv_test <- function(x){
  return(stats::sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE))
}

pv_test <- function(x) {
  ## From: https://github.com/T-Engel/CValternatives/blob/master/R/Functions.R
  ## Citation: https://esajournals.onlinelibrary.wiley.com/doi/pdf/10.1002/ecs2.2527
  n <- length(x)
  pairs <- utils::combn(x, 2)
  min_z <- apply(pairs, 2, min)
  max_z <- apply(pairs, 2, max)
  z <- 1 - (min_z/max_z)
  return(2*sum(z)/(n*(n-1)))
}

d_test <- function(x, k = NULL){
  n = length(x)
  if(is.null(k)){
    k = mean(x, na.rm = TRUE) * .01
  }

  d_list <- NA
  for(i in (1:(n-1))){
    d_list[i]=abs(log((x[i+1] + k)/(x[i] + k)))
  }
  D=sum(d_list)/(n-1)
  return(D)
}
#' Consecutive Disparity Index
#'
#' Calculate D of the numeric vector x
#'
#' @param x numeric vector of a time-series
#' @param k constant that defaults to 1% of the time-series mean
#' @param R
#' @param type
#' @param seed
#' @param conf
#'
#' @return a numeric value
#' @export
#'
#' @examples
#' x <- c(0.2, 0.5, 1.1, 1.4, 1.8, 2.3, 2.5, 2.7, 3.5, 4.4,
#'        4.6, 5.4, 5.4, 5.7, 5.8, 5.9, 6.0, 6.6, 7.1, 7.9)
#' D(x)
D  <- function(x,
               k = NULL,
               R = 1000,
               type = "basic",
               seed = NULL,
               conf = 0.95){

  if(!is.null(seed)){
    set.seed(seed)
  }


  sample_d <- function(x, i){
    return(d_test(x[i]))
  }

  # est <- d(x)

  boot_d <-  boot::tsboot(tseries = x,
                          statistic = sample_d,
                          l = round(length(x)^(1/3)),
                          R = R,
                          sim = "geom")

  boot_d_ci <- boot::boot.ci(boot.out = boot_d,
                             conf = conf,
                             type = type)

  lower <- ifelse(type == "normal",
                  boot_d_ci[[type]][2],
                  boot_d_ci[[type]][4])
  upper <- ifelse(type == "normal",
                  boot_d_ci[[type]][3],
                  boot_d_ci[[type]][5])

  return(list(est =  boot_d_ci$t0,
              lower = lower,
              upper = upper,
              type = type))
}

#' Title
#'
#' @param x
#' @param R
#' @param type
#' @param seed
#' @param conf
#'
#' @return
#' @export
#'
#' @examples
ACF  <- function(x,
                 R = 1000,
                 type = "basic",
                 seed = NULL,
                 conf = 0.95){

  if(!is.null(seed)){
    set.seed(seed)
  }

  sample_acf <- function(x, i){
      return(acf(x[i], plot = FALSE)$acf[2])
  }


  boot_acf <-  boot::tsboot(tseries = x,
                          statistic = sample_acf,
                          l = round(length(x)^(1/3)),
                          R = R,
                          sim = "geom")

  boot_acf_ci <- boot::boot.ci(boot.out = boot_acf,
                             conf = conf,
                             type = type)

  lower <- ifelse(type == "normal",
                  boot_acf_ci[[type]][2],
                  boot_acf_ci[[type]][4])
  upper <- ifelse(type == "normal",
                  boot_acf_ci[[type]][3],
                  boot_acf_ci[[type]][5])

  return(list(est =  boot_acf_ci$t0,
              lower = lower,
              upper = upper,
              type = type))
}


# This is our example function. cg_subsize is the coarse_grain length, see
# ?coarse_grain for details
#' Title
#'
#' @param mat
#' @param cg_subsize
#'
#' @return
#' @export
#'
#' @examples
na_aware_ews <- function(mat, cg_subsize) {
  # require(moments) # required packages
  # require(raster)
  mat_coarse <- spatialwarnings::coarse_grain(mat, cg_subsize)

  c(skewness = moments::skewness(as.vector(mat_coarse), na.rm = TRUE),
    variance = var(as.vector(mat_coarse), na.rm = TRUE),
    moran    = raster::Moran(raster::raster(mat_coarse)))
}


# This is the function that given a matrix, returns a randomized matrix
randomize_matrix_no_na <- function(mat) {
  mat[!is.na(mat)] <- sample(mat[!is.na(mat)])
  return(mat)
}
