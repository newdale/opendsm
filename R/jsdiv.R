#' Jensen Shannon Divergence
#'
#' This function calculates the Jensen-Shannon Divergence
#' score between two probability distributions, P and Q.
#' In sampling design for digital soil mapping applications,
#' it can be used to determine the statistical distance between
#' a sample plan and the population from which the sample plan
#' was derived
#'
#'
#' @param p numeric
#' @param q numeric
#' @param type character Can be one of 'prob' or 'count'.
#' @param unit character Can be one of log base 2 ('log2'), log base 10 ('log10') or natural log ('log').
#'
#' @return numeric
#' @export
#'
#' @examples
#' #let us create data for testing, we need 2 probability distributions
#' prob1 <- c(0.002, 0.020, 0.127, 0.343, 0.362, 0.119, 0.025, 0.002)
#' prob2<- c(0.001, 0.019, 0.325, 0.145, 0.326, 0.028, 0.153, 0.003)
#' sum(prob1)
#' sum(prob2)
#' barplot(prob1)
#' barplot(prob2)
#' # we now run the JS Divergence
#' # it is a symmetrical test, which means Q|P == P|Q
#' # if we enter the same distribution as both P and Q, we confirm a score of 0 or no divergence
#' jsdiv(p=prob1, q=prob1, type='prob', unit='log2')
#' # P|Q
#' jsdiv(p=prob1, q=prob2, type='prob', unit='log2')
#' # Q|P
#' jsdiv(p=prob2, q=prob1, type='prob', unit='log2')
#'
jsdiv<- function(p, q, type=NULL, unit='log2'){
  # check to ensure user has specified data type as 'prob' or 'count' and stop function if not set
  if(is.null(type)) {stop('The function argument type must be either prob or count')
  }
  else if (type=='prob'){
    p<- p
    q<- q
  }
  else {
    p<- p/sum(p)
    q<- q/sum(q)
  }

  m <- 0.5 * (p + q)
  js<- 0.5 * kldiv(p, m, type, unit) + 0.5 * kldiv(q, m, type, unit)
  return(js)
}
