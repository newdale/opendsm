#' Perform Sequential Variable Inflation Factor Analysis
#'
#'Perform sequential Variable Inflation Factor analysis for a dataframe of predictor variables with dataframe columns representing predictor variable values
#'
#' @param cov_df data.frame* object where each column represent covariate values.
#' @param thresh numeric* object Perform sequential VIF until all covariates have VIF scores lower than this value
#' @param trace Boolean* object (default = FALSE). Should the results of each VIF be shown in the console?
#' @param show.R2.vals Boolean* object (default = FALSE) Should the VIF score be also described as an R squared value? (VIF = 1/1-R2)
#'
#' @return Returns a list of 4 objects: 1. Covariates_retained: First object is the column names of the covariates that are kept (i.e. have VIF scores lower than threshold).
#' 2. Covariates_removed: Second object is the column names of the covariates were removed (i.e., have VIF scores higher than threshold).
#' 3. VIF_removed: Third object is a summary of the VIF scores of each covariate that was removed. NA indicates scores fell below threshold and removal was not performed.
#' 4. VIF_all: Forth object is a complete report of each VIF score for each covariate for each time VIF was run.
#'
#' @noRd
#'
#' @importFrom fmsb VIF
#'
#' @examples
#' #Perform sequential VIF on an environmental raster stack
#' library(terra)
#' library(fmsb)
#'
#' #Import example data
#' data(covs)
#'
#' #Run to perform sequential VIF analyses,
#' # removing the covariate with the highest VIF value,
#' # then repeating until threshold is hit
#'
#' vif_results <- seqVIF(cov_df=covs, thresh=5, trace=FALSE, show.R2.vals=FALSE)
#'

seqVIF <- function(cov_df, thresh, trace=FALSE, show.R2.vals=FALSE){

  ###Load required library or throw error
  if(any(!'data.frame' %in% class(cov_df))){cov_df <- data.frame(cov_df)}

  ###Get initial vif value for each covariate and confirm that at least VIF is above threshold
  vif_init <- NULL
  for(covar in names(cov_df)){
    regressors <- names(cov_df)[names(cov_df) != covar]
    form <- stats::formula(paste(covar, '~', paste(regressors, collapse = '+')))
    vif_init <- rbind(vif_init, c(covar, fmsb::VIF(stats::lm(form, data = cov_df))))
  }
  max_vif <- max(as.numeric(vif_init[,2]), na.rm = TRUE)

  ###If no values are above threshold, use all covariates
  if(max_vif < thresh){
    return(list(Covariates_retained=colnames(cov_df),
                VIF_all=vif_init))
    }

  ###If at least one value is above threshold, run sequential VIF
  else{

    in_dat <- cov_df #create copy to iterate over
    #Set empty vectors to store outputs
    names_kept_vec <- NULL; names_rem_vec <- NULL; vif_dfs <- NULL
    vif_rem_vals <- data.frame(Covariate = character(), VIF_at_removal = numeric()) #Setting as numeric doesn't really matter since it gets removed later

    #Backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(max_vif >= thresh){
    if(ncol(in_dat) <= 2){print("Removed all covariates but two, stopping"); break} #If not enough covariates to do linear model, stop - shouldn't this be 1?
      vif_vals <- NULL
      for(covar in names(in_dat)){
        regressors <- names(in_dat)[names(in_dat) != covar]
        form <- stats::formula(paste(covar, '~', paste(regressors, collapse = '+')))
        vif_vals <- rbind(vif_vals, c(covar, fmsb::VIF(stats::lm(form, data = in_dat))))
      }

      #Record max value
      max_vif <- max(as.numeric(vif_vals[,2]), na.rm = TRUE)
      #if(length(which(vif_vals[,2] == max_vif)) > 1){paste0("VIF scores tied, removing alphabetically")} #This should be incredibly rare/impossible
      max_row <- which(vif_vals[,2] == max_vif)[1]

      #We need this break so it doesn't remove the next variable below the threshold
      if(max_vif<thresh) break

      #Print output of each iteration
      if(trace){print(vif_vals); cat('\n'); cat('removed: ', vif_vals[max_row,1], max_vif,'\n\n')}

      #Record names you are removing
      names_rem_vec = append(names_rem_vec, vif_vals[max_row,1])
      #Record the VIF value that is being removed
      vif_rem_vals <- rbind(vif_rem_vals, vif_vals[max_row,])

      #Remove covariate from dataframe
      in_dat <- in_dat[,!names(in_dat) == vif_vals[max_row,1]]

      #Store VIF table
      vif_dfs <- c(vif_dfs, list(vif_vals))
    }

    #Record names you are keeping
    names_kept_vec <- names(in_dat)

    #Update columns of vif_rem_vals with covariates that were not removed during VIF because either:
    #they were below threshold or because linear model breaks down with 2 covariates
    colnames(vif_rem_vals) <- c("Covariate", "VIF_score_at_removal")
    vif_rem_vals[,2] <- as.numeric(vif_rem_vals[,2])

    below_rows = data.frame(Covariate = as.character(vif_dfs[[length(vif_dfs)]][,1]),
                            VIF_score_at_removal = as.numeric(vif_dfs[[length(vif_dfs)]][,2])
                            )
    below_rows = below_rows[order(-below_rows$VIF_score_at_removal),][-1,] #remove first entry since it's already been removed
    below_rows[,2] = NA

    vif_rem_vals <- rbind(vif_rem_vals, below_rows)

    #If you want R2 and pearson r values. Not sure how to best interpret these
    if(show.R2.vals){vif_rem_vals$R2 <- 1 - 1/(vif_rem_vals$VIF_score_at_removal)
    vif_rem_vals$r <- sqrt(vif_rem_vals$R2)}

    #This list is ordered from first VIF to last VIF
    #(i.e., first variables to be removed and first VIF tables run are at the start)
    return(list(Covariates_retained=names_kept_vec,
                Covariates_removed=names_rem_vec,
                VIF_removed=vif_rem_vals,
                VIF_all=vif_dfs))

  }
}
