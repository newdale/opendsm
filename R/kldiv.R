#' Kullback-Leibler Divergence
#'
#' This function calculates the Kullback-Leibler Divergence
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
#' # we now run the KL Divergence
#' # it is a non-symmetrical test, which means you get a different score for Q|P and P|Q
#' # if we enter the same distribution as both P and Q, we confirm a score of 0 or no divergence
#' kldiv(p=prob1, q=prob1, type='prob', unit='log2')
#' # P|Q
#' kldiv(p=prob1, q=prob2, type='prob', unit='log2')
#' # Q|P
#' kldiv(p=prob2, q=prob1, type='prob', unit='log2')
#'

kldiv<- function(p, q, type=NULL, unit= 'log2'){
  # create an empty dataframe to hold results from loop
  t<- data.frame()

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

  # calculate the KL Divergence using the user-specified log units.
  # options are log2 (log base 2), log10 (log base 10) and log (natural log)
  if(unit=='log2'){
    for(i in 1:length(p)){
      r<- sum(p[i] * log2(p[i]/q[i]))
      ifelse(exists("t"), t<- rbind(t,r), t<- r)}
  } else if (unit=='log10'){
    for(i in 1:length(p)){
      r<- sum(p[i] * log10(p[i]/q[i]))
      ifelse(exists("t"), t<- rbind(t,r), t<- r)}
  } else if (unit=='log') {
    for(i in 1:length(p)){
      r<- sum(p[i] * log(p[i]/q[i]))
      ifelse(exists("t"), t<- rbind(t,r), t<- r)}
  } else {stop('Function argument unit must be either: log, log2 or log10')}

  # return the KL Divergence score
  return(sum(t))
}
