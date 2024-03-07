##############################################################################################################
#'
#' Determine minimum sample size for the clhs algorithm
#'
#' This script is based on the publication:
#' Determining minimum sample size for the conditioned Latin hypercube sampling algorithm
#' DOI: https://doi.org/10.1016/j.pedsph.2022.09.001
#'
#' @param covs dataframe*, data frame or SpacVector or SpatRaster, environmental covariates to be used for the sample size determination
#' @param lq, numeric*, lower quantile of choice for the confidence interval, default value = 0.05
#' @param uq, numeric*, upper quantile of choice for the confidence interval, default value = 0.95
#' @param vif, logical*, use variance inflation factor to reduce multicolinearity, default = FALSE
#' @param vif.thresh, numeric*, threshold to be used for VIF, typically 5 or 10, default = 5
#'
#' @return list, with 7 objects
#' Object 1. Dataframe 'Sample_Size' returning descriptive statistics for the minimum sample size
#' Object 2. Dataframe 'Transformations' provides data transformations used to normalize covariate distriburions
#' Object 3. Dataframe 'Results' provides the minimum sample size for each individual covariate
#' Object 4. Dataframe 'Bimodality' provides the results of bimodality assessment for covariates deemd to be bimodal
#' Object 5. List 'Hist_Original' containing histogtram objects for each covariate
#' Object 6. List 'Hist_Transformed' containing histogram objects for all transformed covariates, including the bimodal covariates whose distributions were split
#' Object 7. Dataframe with sample plan developed using median of minimum sample size
#' Object 8. List 'VIF_Results' containing the outputs of the variance inflation factor analysis (optional)
#'
#' @export
#'
#' @importFrom mousetrap bimodality_coefficient
#' @importFrom diptest dip.test
#' @importFrom nortest ad.test
#' @importFrom DescTools Closest
#' @importFrom bestNormalize bestNormalize
#' @importFrom mixtools normalmixEM
#' @importFrom onsoilsurvey oss.seqVIF
#' @importFrom terra crds values
#'
#' @examples
#' data(covs)
#' out<- clhs_min(covs=covs, lq=0.05,uq=0.95, vif=TRUE, vif.thresh=10)
#'
#'
#'
#'
#'
clhs_min<- function(covs,lq=0.05,uq=0.95, vif=FALSE, vif.thresh=5){

  # preprocessing of input data
  # check the format of the covs object and convert as required

  if(methods::is(covs,"SpatVector")){

    # get the data and the XY coordinates
    coords<- terra::crds(covs, df=TRUE, list=FALSE)
    covs<- terra::values(covs, dataframe=TRUE, na.rm=TRUE)

  } else if(methods::is(covs,"SpatRaster")){

    # get the data and the XY coordinates
    coords<- terra::crds(covs, df=TRUE, na.rm=TRUE)
    covs<- terra::values(covs, dataframe=TRUE, na.rm=TRUE)

  } else if(methods::is(covs,"data.frame")) {
    covs<- covs
    coords<- NA
  } else {stop('The function argument "covs" must be either "SpatVector" or "SpatRaster" or "data.frame"')}


  # Step 0: use variance inflation factor to remove multicollinearity in the predictors
  if(vif==TRUE){
    y<- onsoilsurvey::oss.seqVIF(covs,thresh=vif.thresh, trace=FALSE, show.R2.vals=FALSE)
    covs<- covs[,y$Covariates_retained]}



  # Step 1: test for bimodal distribution of covariates

  for (j in 1:ncol(covs)){
    # check for multimodality using the bimodality coefficient test and Hartigan's Dip Test
    yy<- round(mousetrap::bimodality_coefficient(covs[,j]),4)
    xx<- diptest::dip.test(covs[,j])
    xx<- xx[[2]]
    #run anderson darling test for normality
    zz<- nortest::ad.test(covs[,j])
    zz<- zz$p.value

    #bind the results
    temp<- cbind(colnames(covs[j]),yy,xx,zz)
    ifelse(exists('bidat.f'), bidat.f<- rbind(bidat.f,temp), bidat.f<- temp)
    colnames(bidat.f)<- c("Covariate","Bimodal_Coeff","DipTest","AD_Test")
  }

  #convert the matrix to dataframe and specify score as numeric
  bidat.f<-as.data.frame(bidat.f)
  bidat.f[,2]<- as.numeric(bidat.f[,2])
  bidat.f[,3]<- as.numeric(bidat.f[,3])
  bidat.f[,4]<- as.numeric(bidat.f[,4])

  # calculate the p-value to compare to the diptest p-value as per https://www.hindawi.com/journals/mpe/2019/4819475/
  bidat.f$BC_PVAL<- sqrt((0.32-0.05)^2 * bidat.f$Bimodal_Coeff) + 0.05

  # create 3 vectors of covariate names for unimodal normal, unimodal not normal and finally multimodal
  bivec<- sort(bidat.f$Covariate[bidat.f[,3]<bidat.f[,5]])
  univec<-bidat.f$Covariate[bidat.f[,3]>=bidat.f[,5] & bidat.f[,4]<0.05]
  uninorm<- bidat.f$Covariate[bidat.f[,3]>=bidat.f[,5] & bidat.f[,4]>=0.05]

  if(length(bivec)>0){
    print("The following covariates were detected as bimodal")
    cat(bivec,sep="\n")}

  #create a list of vectors from the unimodal not normal data
  if(length(univec>0)){
    dat2<- as.list(as.data.frame(covs[,univec]))
    names(dat2)<- univec
  }else{dat2<- list()}

  #create a list of vectors from the unimodal normal data
  if(length(uninorm>0)){
    dat<- as.list(as.data.frame(covs[,uninorm]))
    names(dat)<- uninorm
  }else{dat<- list()}

  # now we need to run the separation of the multimodal covariates, but only if we have some
  if(length(bivec>0)){
    bi.df<- as.data.frame(covs[,bivec])
    names(bi.df)<- bivec
    bi.df<- bi.df[,order(names(bi.df))]
    bi.df<- as.data.frame(bi.df)
    for(j in 1:ncol(bi.df)){
      z<-sort(bi.df[[j]])
      set.seed(6)
      invisible(utils::capture.output(mixmdl<- mixtools::normalmixEM(z, maxrestarts = 1000)))
      # determine lower and upper limits of 95% CI based on mean and sd of populations assumes normal distribution
      ll1<-stats::qnorm(0.025,mixmdl$mu[1], mixmdl$sigma[1])
      ul1<-stats::qnorm(0.975,mixmdl$mu[1], mixmdl$sigma[1])
      ll2<-stats::qnorm(0.025,mixmdl$mu[2], mixmdl$sigma[2])
      ul2<-stats::qnorm(0.975,mixmdl$mu[2], mixmdl$sigma[2])
      # we use the limits determined above to create the new vectors for the individual populations
      z1<- z[z<= ul1 & z>= ll1]
      z2<- z[z<= ul2 & z>= ll2]
      dat2[[length(dat2)+1]]<- z1
      names(dat2)[length(dat2)]<- paste0(bivec[j],"_",1)
      dat2[[length(dat2)+1]]<- z2
      names(dat2)[length(dat2)]<- paste0(bivec[j],"_",2)
      # and we can plot the histograms and the modeled distributions to see if it works OK
      graphics::hist(z,freq=FALSE, plot=TRUE, breaks=20,main=paste0(bivec[j]))
      graphics::lines(z,mixmdl$lambda[1]*stats::dnorm(z,mixmdl$mu[1],mixmdl$sigma[1]), col="green", lwd=2)
      graphics::lines(z,mixmdl$lambda[2]*stats::dnorm(z,mixmdl$mu[2],mixmdl$sigma[2]), col="red", lwd=2)
    }
  }

  # and finally we use the bestNormalize on the list dat2 which has the not unimodal not normal and split bimodal covariates
  # and we add them to the dat list which currently holds only the unimodal normal
  if(length(dat2)>0){
    for (i in 1:length(dat2)){
      bnTemp<- bestNormalize::bestNormalize(dat2[[i]], allow_orderNorm = FALSE, standardize=TRUE, out_of_sample = FALSE)
      dat[[length(dat)+1]]<- bnTemp$x.t
      names(dat)[length(dat)]<- names(dat2[i])
      bnName<- names(bnTemp$norm_stats)[which.min(bnTemp$norm_stats)]

      ## Bind the outputs
      temp<- cbind(names(dat2[i]),bnName)
      ifelse(exists('tname'), tname<- rbind(tname,temp), tname<- temp)
    }
    tname<- as.data.frame(tname)
    colnames(tname)<- c("Covariate","Transformation")}

  #Create empty lists to store outputs
  histo.orig<- list()
  histo.normal<- list()
  histo.order<- c(uninorm, univec,sort(rep(bivec,2)))

  # here we will store histograms for each covariate in its original form
  for (i in 1:length(covs)){
    n<- grDevices::nclass.FD(covs[[i]])
    d1<- graphics::hist(covs[[i]], breaks=seq(min(covs[i]),max(covs[i]),l=n+1), plot=FALSE)
    histo.orig[[i]] <- d1
  }

  # for each covariate we create a histogram with normalized data, and determine sample size
  for (i in 1:length(dat)){ # for each covariate in the dataframe

    ## Determine the optimal number of bins using the Freedman-Diaconis rule
    n<- grDevices::nclass.FD(dat[[i]])

    ## Create the histogram and set the number on bins based on the FD above
    graphics::par(mfrow=c(1,2))
    d1<- graphics::hist(dat[[i]], breaks=seq(min(dat[[i]]),max(dat[[i]]),l=n+1), plot=FALSE)

    ## Create an index based on bins with counts > 0 and extract mids and density associated with this
    ## I use a if...else loop so that if all bins have counts, we use the data as is, if some bins have zero counts, then we remove these
    idx<- as.logical(d1$counts)
    if(all(idx)==TRUE){
      d.mids<- d1$mids
      d.density<- d1$density
    }else{
      d.mids<- d1$mids[idx]
      d.density<- d1$density[idx]
    }

    ## Generate the quantiles for the population, default values are 5% and 95%
    a<- as.vector(stats::quantile(dat[[i]], probs = c(lq, uq)))

    ## We then identify the positions along the vector "a" from previous step where values are closest to bin midpoints
    ## We use min and max in the event the Closest function returns a list of 2 vectors instead of a vector which happens if the difference is equal
    bounds<- c(DescTools::Closest(d.mids,a[1],which=TRUE), DescTools::Closest(d.mids,a[2],which=TRUE))
    bounds.n<- c(d.mids[bounds[1]],d.mids[bounds[2]])

    ## We now use the positions determined above to extract the mean density of the 5% and 95% locations. We use mean of the two values which allows better estimates for skewed distributions
    f<- mean(c(d.density[bounds[[1]]],d.density[bounds[[2]]]))

    ## Store the histogram data to object hito.normal
    d2<- graphics::hist(dat[[i]], breaks=seq(min(dat[[i]]),max(dat[[i]]),l=n+1), plot=FALSE)
    histo.normal[[i]] <- d2

    ## Step 10: Finally we divide the density dataframe by the value from previous step, round up and add the values to get required sample size
    s<- sum(ceiling(d1$density/f))

    ## Step 11: Bind the outputs
    temp<- cbind(names(dat)[i],n,s, bounds.n[1],bounds.n[2])
    ifelse(exists('dat.f'), dat.f<- rbind(dat.f,temp), dat.f<- temp)
    colnames(dat.f)<- c("Covariate","FD_Bins","Samples","Lower","Upper")
  }

  ## Organize the outputs and return results
  dat.f<- as.data.frame(dat.f)
  dat.f[,2]<- as.numeric(dat.f[,2])
  dat.f[,3]<- as.numeric(dat.f[,3])

  sum.f<- cbind(min(dat.f[,3], na.rm=TRUE),
                round(mean(dat.f[,3], na.rm=TRUE),1),
                stats::median(dat.f[,3], na.rm=TRUE),
                max(dat.f[,3]),
                round(stats::sd(dat.f[,3], na.rm=TRUE),1))
  sum.f<- as.data.frame(sum.f)
  colnames(sum.f)<- c("Minimum", "Mean", "Median", "Maximum","StDev")
  print(sum.f)

  if(exists("tname")){tname<-tname
  }else{tname<- c("No data transformations required")}

  #generate a sample plan
  plan_idx<- clhs::clhs(x=covs, size=sum.f[1,3], iter=10000, simple=TRUE, progress=FALSE)

  if(methods::is(coords,"logical")){plan=covs[plan_idx,]

  }else{plan<- cbind(coords[plan_idx,],covs[plan_idx,])}

  if(vif==FALSE){
  clhs.size.out<- list(sum.f,tname,dat.f,bidat.f, histo.orig, histo.normal, plan)
  names(clhs.size.out)<- c("Sample_Size","Transformations","Results","Bimodality",
                           "Hist_Original","Hist_Transformed", "Sample_Plan")
  }else if(vif==TRUE){
    clhs.size.out<- list(sum.f,tname,dat.f,bidat.f, histo.orig, histo.normal,plan,y)
    names(clhs.size.out)<- c("Sample_Size","Transformations","Results","Bimodality",
                           "Hist_Original","Hist_Transformed", "Sample_Plan", "VIF_Results")}

  return(clhs.size.out)
}
