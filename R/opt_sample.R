##############################################################################################################
#'
#' Determine optimal sample size for the conditioned Latin hypercube
#' or feature (covariate) coverage sampling algorithms
#' using three different divergence metrics:
#' Kullback-Leibler Divergence, Jenson-Shannon Divergence,
#' and Jenson-Shannon Distance
#'
#' This script is based on the publication:
#' Divergence metrics for determining optimal training sample size in digital soil mapping
#' DOI: https://doi.org/10.1016/j.geoderma.2023.116553
#'
#' @param alg, character, either "clhs" or "fcs", default is "clhs"
#' @param s_min numeric, minimum sample size to be tested
#' @param s_max numeric, maximum sample size to be tested
#' @param s_step numeric, sample size step for testing, default is 20
#' @param s_reps numeric, number of repetitions of the sample design, default is 2, must be >1
#' @param bins numeric, number of bins for the divergence metric calculation, default is 30
#' @param covs data frame, environmental covariates to be used for the sample size determination
#' @param clhs_iter numeric, number of iterations to use for the clhs, default is 10000
#' @param cpus numeric, number of CPUs to be used for parallel processing, default is 1/2 of virtual CPUs
#' @param conf numeric, confidence used to select optimal sample size from the CDF, default is 0.95
#'
#' @return list with 3 dataframes, 1 list, and 2 figures are plotted
#' Dataframe 1: Dataframe with optimal sample size based on the three divergence metrics
#' Dataframe 2: List of vectors showing the index of rows for each sample plan evaluated
#' Dataframe 3: Dataframe providing mean and standard deviation of divergence metrics across the repetitions
#' Dataframe 4: Dataframe showing all raw data - divergence metrics for each covariate at each sample size
#' Figure 1: plots of exponential decay of mean divergence metrics across all covariates with increasing sample size
#' Figure 2: plots of exponential decay of normalized divergence metrics for individual covariates with increasing sample size
#'
#' @import clhs
#' @import ggplot2
#' @import stats
#' @import fields
#' @importFrom dplyr summarise_all group_by
#' @importFrom magrittr "%>%"
#' @importFrom snowfall sfInit sfLibrary sfLapply sfStop
#' @importFrom ggpubr ggarrange annotate_figure text_grob
#' @importFrom matrixStats colSds
#' @importFrom reshape2 melt
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @export
#'
#' @examples
#' # import the example data and run with 4 reps
#' data(covs)
#' out<- opt_sample(alg="clhs", s_min=10, s_max=400, s_step=20, s_reps=4, covs= covs, clhs_iter=100, cpus=NULL, conf=0.95)
#'

opt_sample<- function(alg="clhs", s_min, s_max, s_step=20, s_reps=2, bins=30, covs, clhs_iter=10000, cpus=NULL, conf=0.95){

  mean_JS <- mean_JSdist <- mean_KL <- samp_nos <- sd_JS <- sd_JSdist <- sd_KL <- value <- variable<- NULL

  ##############################################################################################################
  ##############################################################################################################
  #    Step 1. Set up variables and create replicated sample plans
  ###############################################################################################################
  ###############################################################################################################

  s_min<- s_min
  s_max<- s_max
  s_step<- s_step
  bins<- bins
  covs<- covs
  clhs_iter<- clhs_iter

  cpus<- if(is.null(cpus)){parallel::detectCores()*0.5
  }else{cpus<- cpus}

  cseq<- round(seq(s_min, s_max, (s_max-s_min)/s_step))

  # prepare a list of the dataframe from which to sample, replicate the covariates 's_reps' times
  covs.list<- rep(list(covs),s_reps) # note the index start here

  #here we use if statement and run either clhs or fcs, as requested by user
  if(is.null(alg)){stop('The function argument "type"alg" must be either "clhs" or "fcs"')

    }else if(alg=="clhs"){

    # prepare clhs package for parallel computing
      snowfall::sfInit(parallel = TRUE, cpus = cpus, type = "SOCK")
      suppressMessages(snowfall::sfLibrary(clhs))

    # create the custom function using clhs
      clhs.fun<- function(x, size=size) {
        clhs::clhs(x = x, size=size, iter=clhs_iter, simple =TRUE, progress=FALSE)}

    # loop over sample sizes and create a list of sample indexes
      for (w in 1:length(cseq)){ # for each sample size
        s.size=cseq[w]
        ## run clhs against the list using sfLapply for each value of w
        x <- snowfall::sfLapply(covs.list, fun = clhs, size=s.size)
        ifelse(exists("s.list"), s.list<- append(s.list,x), s.list<- x)
        rm(x)
        }

      snowfall::sfStop()
      rm(covs.list)

      }else if(alg=="fcs"){

        snowfall::sfInit(parallel = TRUE, cpus = cpus, type = "SOCK")
        suppressMessages(snowfall::sfLibrary(stats))
        suppressMessages(snowfall::sfLibrary(fields))

        # create the custom function using fcs

        fcs.fun<- function(x, size=size) {
          myClust<- stats::kmeans(scale(x), centers=size, iter.max=1000, nstart=20) #set the arguments in kmeans as desired
          x$clusters<- myClust$cluster
          rdist.out<- fields::rdist(x1=myClust$centers, x2=scale(x))
          ids.mindist<- apply(rdist.out,MARGIN=1,which.min)
          return(ids.mindist)
          }

    # loop over sample sizes and create a list of sample indexes
        for (w in 1:length(cseq)){ # for each sample size
          s.size=cseq[w]  # sample size

    ## run clhs against the list using sfLapply for each value of w
          x <- snowfall::sfLapply(covs.list, fun = fcs.fun, size=s.size)
          ifelse(exists("s.list"), s.list<- append(s.list,x), s.list<- x)
          rm(x)
          }

        snowfall::sfStop()
        rm(covs.list)
        }

  ##############################################################################################################
  ##############################################################################################################
  #    Step 2. Determine divergence indices and capture outputs
  ###############################################################################################################
  ###############################################################################################################

  # q.mat, cov.mat, and h.mat code chunks credited to Malone et al., 2019
  # from sample code provided in https://doi.org/10.7717/peerj.6451

  #quantile matrix of the covariate data
  q.mat<- matrix(NA, nrow=(bins+1), ncol= ncol(covs))
  j=1
  for (i in 1:ncol(covs)){ #note the index start here
    #get a quantile matrix together of the covariates
    ran1<- max(covs[,i]) - min(covs[,i])
    step1<- ran1/bins
    q.mat[,j]<- seq(min(covs[,i]), to = max(covs[,i]), by =step1)
    j<- j+1}


  cov.mat<- matrix(1, nrow=bins, ncol=ncol(covs))
  for (i in 1:nrow(covs)){ # the number of pixels
    cntj<- 1
    for (j in 1:ncol(covs)){ #for each column
      dd<- covs[i,j]
      for (k in 1:bins){  #for each quantile
        kl<- q.mat[k, cntj]
        ku<- q.mat[k+1, cntj]
        if (dd >= kl & dd <= ku){cov.mat[k, cntj]<- cov.mat[k, cntj] + 1}
      }
      cntj<- cntj+1
    }
  }

  ####################################################################
  # this loop will calculate the kl divergence, jenson divergence, and jensen distance, JS Distance,
  # for each of length(cseq) * s_reps  sample designs

  print(paste("COMPUTING DIVERGENCE METRICS FOR  ", length(cseq), "  SAMPLE SIZES * ",
              s_reps, "  REPETITIONS =  ", length(cseq)*s_reps, sep = ""))

  pb = utils::txtProgressBar(min = 0, max = length(s.list), initial = 0)

  for (nn in 1:length(s.list)){
    s.df<- covs[s.list[[nn]],]

    ####Compare whole study area covariate space with the selected sample
    h.mat<- matrix(1, nrow=bins, ncol=ncol(covs))

    for (ii in 1:nrow(s.df)){ # the number of observations
      cntj<- 1
      for (jj in 1:ncol(s.df)){ #for each column
        dd<- s.df[ii,jj]
        for (kk in 1:bins){  #for each quantile
          kl<- q.mat[kk, cntj]
          ku<- q.mat[kk+1, cntj]
          if (dd >= kl & dd <= ku){h.mat[kk, cntj]<- h.mat[kk, cntj] + 1}
        }
        cntj<- cntj+1}}


    ## First test: Kullback-Leibler (KL) divergence
    kl_loop<- matrix(NA,ncol=1, nrow=ncol(covs))
    for (m in 1:ncol(covs)){
      klo<- kldiv(c(cov.mat[,m]), c(h.mat[,m]), type='count', unit='log2')
      kl_loop[m,1]<- klo
    }
    klo<- mean(kl_loop[,1])
    #mat.f[j,1]<- klo  # value of 0 means no divergence

    ## Second test: Jensen Shannon (JS) divergence
    #Jensen-Shannon divergence (operates on densities, not counts)
    jsd_loop<- matrix(NA,ncol=1, nrow=ncol(covs))
    for (m in 1:ncol(covs)){
      jsd<- jsdiv(c(cov.mat[,m]/sum(cov.mat[,m])), c(h.mat[,m]/sum(h.mat[,m])), type= 'count', unit='log2')
      jsd_loop[m,1]<- jsd
    }
    jsd<- mean(jsd_loop[,1])
    #mat.f[j,2]<- jsd # value of 0 means no divergence

    ## Third test: Jensen Shannon Distance (JSdist)
    jsdist_loop<- matrix(NA,ncol=1, nrow=ncol(covs))
    for (m in 1:ncol(covs)){
      jsdist<- jsdist(c(cov.mat[,m]/sum(cov.mat[,m])), c(h.mat[,m]/sum(h.mat[,m])), type= 'count', unit='log2')
      jsdist_loop[m,1]<- jsdist
    }
    jsdist<- mean(jsdist_loop[,1])

    # here we will compile the mean test scores for each of j iterations
    mat.f.temp<- t(rbind(nrow(s.df),klo,jsd,jsdist))
    ifelse(exists('mat.f'),mat.f<- rbind(mat.f,mat.f.temp),mat.f<- mat.f.temp)

    # here we assign the individual scores to a matrix for further analysis
    det.f.temp<- t(rbind(nrow(s.df), kl_loop, jsd_loop, jsdist_loop))
    ifelse(exists('det.seq'),det.seq<- rbind(det.seq,det.f.temp),det.seq<- det.f.temp)

    # clean up a few temporary files
    rm(h.mat, mat.f.temp, det.f.temp, jsd_loop, jsdist_loop, kl_loop, klo, jsd, jsdist)

    utils::setTxtProgressBar(pb,nn)

    #print(paste("COMPLETED  ", nn, "  OF  ", length(s.list),sep = ""))

  } # end of the for nn loop

  close(pb)


  # and now we generate the mean results for each sample intensity
  for (n in 1:length(cseq)){
    s.size<- cseq[n]
    x<- mat.f[mat.f[,1] == s.size,]
    x<- x[,2:ncol(x)]
    a<- colMeans(x)
    aa<- matrixStats::colSds(x)
    aaa<- c(a,aa)
    ifelse(exists("mat.seq"),mat.seq<- rbind(mat.seq,aaa), mat.seq<- aaa)
    rm(x,a,aa,aaa)
  }

  # here we combine the sampling nos sequence with the mat.seq matrix
  dat.seq<- as.data.frame(cbind(cseq,mat.seq))

  colnames(dat.seq)<- c("samp_nos", "mean_KL", "mean_JS",'mean_JSdist', 'sd_KL', 'sd_JS', 'sd_JSdist')
  rownames(dat.seq)<- NULL

  # and now we can arrange the detailed outputs
  det.seq<- as.data.frame(det.seq)

  # and here we need to assign names to the columns of det.seq dataframe which holds data for all iterations
  det.seq.names1<- c("samp_nos")

  for (p in 1:ncol(covs)){
    names2<- c(paste0("kld_",colnames(covs)[p]))
    names3<- c(paste0("jsd_",colnames(covs)[p]))
    names4<- c(paste0("jsdist_",colnames(covs)[p]))

    ifelse(exists('det.seq.names2'), det.seq.names2<- c(det.seq.names2,names2), det.seq.names2<- names2)
    ifelse(exists('det.seq.names3'), det.seq.names3<- c(det.seq.names3,names3), det.seq.names3<- names3)
    ifelse(exists('det.seq.names4'), det.seq.names4<- c(det.seq.names4,names4), det.seq.names4<- names4)
  }

  det.seq.names<- c(det.seq.names1,det.seq.names2,det.seq.names3,det.seq.names4)
  colnames(det.seq)<- det.seq.names

  rm(det.seq.names1,det.seq.names2,det.seq.names3,det.seq.names4,names2,names3,names4, det.seq.names)

  ##############################################################################################################
  ##############################################################################################################
  #    Step 3. Create exponential decay functions and CDFs for determining sample size
  ###############################################################################################################
  ###############################################################################################################

  # make exponential decay functions of the KL divergence
  x<- dat.seq$samp_nos
  y1 <- dat.seq$mean_KL #PIP

  #Parametise Exponential decay function
  start1<- list(k= max(y1)*0.9, b1= 0.01, b0= 0.01)
  suppressWarnings(fit1 <- stats::nls(y1 ~ k*exp(-b1*x) + b0, start = start1, stats::nls.control(maxiter=10000)))

  ##########################################################
  # make an exponetial decay function of the JS divergence
  y2 = dat.seq$mean_JS

  #Parametise Exponential decay function
  start2 <- list(k= max(y2)*0.9, b1= 0.01, b0= 0.01)
  suppressWarnings(fit2 <- stats::nls(y2 ~ k*exp(-b1*x) + b0, start = start2, stats::nls.control(maxiter=10000)))

  ##########################################################
  # make an exponetial decay function of the JS Distance
  y3<- dat.seq$mean_JSdist

  #Parametise Exponential decay function
  start3 <- list(k= max(y3)*0.9, b1= 0.01, b0= 0.01)
  suppressWarnings(fit3 <- stats::nls(y3 ~ k*exp(-b1*x) + b0, start = start3, stats::nls.control(maxiter=10000)))

  #############################################################################
  # Now compute CDF from the exponential decay functions above
  #############################################################################

  #vector of values for the x axis for plotting
  xx<- seq(1, max(cseq),1)

  #Apply fit to KL Divergence and convert to CDF
  yy1<- stats::predict(fit1,list(x=xx))

  # determine sample number based on normalized data
  normalized = 1- (yy1-min(yy1))/(max(yy1)-min(yy1))
  z1<- which(normalized >= conf)[1]

  #############################################################################

  #############################################################################
  #Apply fit to JS Divergence and convert to CDF
  yy2<- stats::predict(fit2,list(x=xx))

  # determine sample number based on normalized data
  normalized2 = 1- (yy2-min(yy2))/(max(yy2)-min(yy2))
  z2<- which(normalized2 >= conf)[1]

  #############################################################################

  #############################################################################
  #Apply fit to JS Distance and convert to CDF
  yy3<- stats::predict(fit3,list(x=xx))

  # determine sample number based on normalized data
  normalized3 = 1- (yy3-min(yy3))/(max(yy3)-min(yy3))
  z3<- which(normalized3 >= conf)[1]

  #############################################################################
  ## and finally put the number of samples into a dataframe for each technique
  opt_sites<- data.frame(method = c("Normalized KLD","Normalized JSD","Normalized JS Distance"),
                         Sites = c(z1,z2,z3))

  rm(x,z1,z2,z3)

  ##############################################################################################################
  ##############################################################################################################
  #    Step 4. Create some plots for visualizing results
  ###############################################################################################################
  ###############################################################################################################

  # KL Divergence
  p1_fit<- as.vector(stats::fitted(fit1))
  p1_fit<- as.data.frame(cbind(cseq,p1_fit))

  p1<- ggplot2::ggplot(dat.seq, ggplot2::aes(x= samp_nos, y= mean_KL))+
    ggplot2::geom_point()+
    ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_KL-sd_KL, ymax=mean_KL+sd_KL), width=0.2) +
    ggplot2::geom_line(color="red", linewidth=0.8, data=p1_fit, ggplot2::aes(x=cseq, y=p1_fit)) +
    ggplot2::xlab("") +
    ggplot2::ylab(expression(paste("D"[KL]))) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(size =14),
                   axis.title.y = ggplot2::element_text(size =16),
                   axis.text.x= ggplot2::element_blank())

  p2_fit<- data.frame(xx,normalized)
  p2<- ggplot2::ggplot(p2_fit, ggplot2::aes(x=xx, y=normalized)) +
    ggplot2::geom_line(linewidth=1.1) +
    ggplot2::xlab("") +
    ggplot2::ylab(expression(paste("CDF of  1 - Normalized D"[KL]))) +
    ggplot2::ylim(0,1) +
    ggplot2::geom_segment(color="red",linewidth=0.8,ggplot2::aes(x=0, y=conf, xend=opt_sites[1,2]+25, yend=conf)) +
    ggplot2::geom_segment(color="red",linewidth=0.8,ggplot2::aes(x=opt_sites[1,2], y=0, xend=opt_sites[1,2], yend=1.00)) +
    ggplot2::theme_bw()+
    ggplot2::theme(axis.text = ggplot2::element_text(size =14),
                   axis.title.y = ggplot2::element_text(size =14))

  # JS Divergence
  p3_fit<- as.vector(stats::fitted(fit2))
  p3_fit<- as.data.frame(cbind(cseq,p3_fit))

  p3<- ggplot2::ggplot(dat.seq, ggplot2::aes(x= samp_nos, y= mean_JS))+
    ggplot2::geom_point()+
    ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_JS-sd_JS, ymax=mean_JS+sd_JS), width=0.2) +
    ggplot2::geom_line(color="red", linewidth=0.8, data=p3_fit, ggplot2::aes(x=cseq, y=p3_fit)) +
    ggplot2::xlab("") +
    ggplot2::ylab(expression(paste("D"[JS]))) +
    ggplot2::theme_bw()+
    ggplot2::theme(axis.text = ggplot2::element_text(size =14),
                   axis.title.y = ggplot2::element_text(size =16),
                   axis.text.x= ggplot2::element_blank())

  p4_fit<- data.frame(xx,normalized2)
  p4<- ggplot2::ggplot(p4_fit, ggplot2::aes(x=xx, y=normalized2)) +
    ggplot2::geom_line(linewidth=1.1) +
    ggplot2::xlab("") +
    ggplot2::ylab(expression(paste("CDF of  1 - Normalized D"[JS]))) +
    ggplot2::ylim(0,1) +
    ggplot2::geom_segment(color="red",linewidth=0.8,ggplot2::aes(x=0, y=conf, xend=opt_sites[2,2]+25, yend=conf)) +
    ggplot2::geom_segment(color="red",linewidth=0.8,ggplot2::aes(x=opt_sites[2,2], y=0, xend=opt_sites[2,2], yend=1.00)) +
    ggplot2::theme_bw()+
    ggplot2::theme(axis.text = ggplot2::element_text(size =14),
                   axis.title.y = ggplot2::element_text(size =14))

  # JS Distance
  p5_fit<- as.vector(stats::fitted(fit3))
  p5_fit<- as.data.frame(cbind(cseq,p5_fit))

  p5<- ggplot2::ggplot(dat.seq, ggplot2::aes(x= samp_nos, y= mean_JSdist))+
    ggplot2::geom_point()+
    ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_JSdist-sd_JSdist, ymax=mean_JSdist+sd_JSdist), width=0.2) +
    ggplot2::geom_line(color="red", linewidth=0.8, data=p5_fit, ggplot2::aes(x=cseq, y=p5_fit)) +
    ggplot2::xlab("") +
    ggplot2::ylab(expression(paste("Dist"[JS]))) +
    ggplot2::theme_bw()+
    ggplot2::theme(axis.text = ggplot2::element_text(size =14),
                   axis.title.y = ggplot2::element_text(size =16),
                   axis.text.x= ggplot2::element_blank())

  p6_fit<- data.frame(xx,normalized3)
  p6<- ggplot2::ggplot(p6_fit, ggplot2::aes(x=xx, y=normalized3)) +
    ggplot2::geom_line(linewidth=1.1) +
    ggplot2::xlab("") +
    ggplot2::ylab(expression(paste("CDF of  1 - Normalized Dist"[JS]))) +
    ggplot2::ylim(0,1) +
    ggplot2::geom_segment(color="red",linewidth=0.8,ggplot2::aes(x=0, y=conf, xend=opt_sites[3,2]+25, yend=conf)) +
    ggplot2::geom_segment(color="red",linewidth=0.8,ggplot2::aes(x=opt_sites[3,2], y=0, xend=opt_sites[3,2], yend=1.0)) +
    ggplot2::theme_bw()+
    ggplot2::theme(axis.text = ggplot2::element_text(size =14),
                   axis.title.y = ggplot2::element_text(size =14))


  fig1<- ggpubr::ggarrange(p1,p3,p5,p2,p4,p6,nrow=2,ncol=3, align="v")
  ggpubr::annotate_figure(fig1,bottom=ggpubr::text_grob("Sample Size", size=16))
  plot(fig1)

  ##here we use the detailed data to show the curves for different covariates.
  ## Since we simply calculate the arithmetic mean of the divergence scores across all covariate in the main workflow.

  # calculate means by sample numbers (we have 10 reps of each sampling intensity)
  det.melt<- det.seq %>% dplyr::group_by(samp_nos) %>% dplyr::summarise_all(mean)

  # and now normalize the values between 0 and 1 for each covariate
  det.melt2<- det.melt
  for (q in 2:(ncol(det.melt)-1)){
    diff<- det.melt[,q]
    diff_norm<- (diff - min(diff))/(max(diff)-min(diff))
    det.melt2[,q]<- diff_norm
    rm(diff, diff_norm)
  }

  # KL Divergence
  kld.n<- det.melt2[,c(1,2:11)]
  kld.n<- reshape2::melt(kld.n,id="samp_nos")

  # JS Divergence
  jsd.n<- det.melt2[,c(1,12:21)]
  jsd.n<- reshape2::melt(jsd.n,id="samp_nos")

  # JS Distance
  jsdist.n<- det.melt2[,c(1,22:31)]
  jsdist.n<- reshape2::melt(jsdist.n,id="samp_nos")

  det.plot<- rbind(cbind(metric="KLD",kld.n),
                   cbind(metric="JSD", jsd.n),
                   cbind(metric="JSDist",jsdist.n))

  rm(kld.n, jsd.n, jsdist.n, det.melt, det.melt2)

  det.plot$variable<- as.vector(det.plot$variable)
  det.plot$variable<- gsub("kld_","",det.plot$variable)
  det.plot$variable<- gsub("jsd_","",det.plot$variable)
  det.plot$variable<- gsub("jsdist_","",det.plot$variable)
  det.plot$variable<- as.factor(det.plot$variable)


  lab1<- c("Kullback-Leibler Divergence",
           "Jensen-Shannon Divergence",
           "Jensen-Shannon Distance")
  names(lab1)<- c("KLD","JSD","JSDist")

  fig2<- ggplot2::ggplot(data=det.plot, ggplot2::aes(x=samp_nos, y=value, color=variable))+
    ggplot2::geom_line(linewidth=1)+
    ggplot2::theme_bw()+
    ggplot2::guides(color=ggplot2::guide_legend(ncol=2))+
    ggplot2::facet_wrap(~metric, nrow=2, ncol=2, labeller=ggplot2::labeller(metric=lab1))+
    ggplot2::labs(y= "Normalized Divergence",x="Sample Size") +
    ggplot2::theme(axis.text = ggplot2::element_text(size =14),
                   legend.text=ggplot2::element_text(size=12),
                   legend.text.align=0,
                   axis.title.x = ggplot2::element_text(size =14, face="bold"),
                   axis.title.y = ggplot2::element_text(size =14, face="bold"),
                   strip.text.x = ggplot2::element_text(size=12, face="bold"),
                   legend.title= ggplot2::element_text(size=14, face='bold'),
                   legend.position = c(0.75,0.25))
  plot(fig2)

  print(opt_sites)

  return(out=list(optimal_sites=opt_sites, sample_index=s.list, summary=dat.seq, detailed=det.seq))

  }
