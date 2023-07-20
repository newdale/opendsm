#' Creating_Continuous_Covariates for digital soil mapping
#'
#' This function creates a series of 55 environmental
#' covariates derived from a digital elevation model.
#' These are commonly used covariates in digital soil mapping.
#' These are created using SAGA-GIS and WhiteboxTools through
#' R implementations in RSAGA and whitebox packages
#'
#' @param x SpatRaster
#' @param wdir character User can specify the working directory, defaults to tempdir()
#' @param cpus, numeric Used to reduce the number of cores used by RSAGA
#'
#' @return SpatRaster Returns a SPatRaster with 55 covariates pointing to the rasters in working directory
#'
#' @export
#'
#' @importFrom RSAGA rsaga.import.gdal rsaga.geoprocessor rsaga.env
#' @import whitebox
#' @importFrom terra res writeRaster rast ncol nrow
#'
#' @examples
#' #let us use the Keene digital terrain model as an example
#' data(keene)
#' keene<- terra::rast(keene, type="xyz")
#' terra::crs(keene)<- "epsg:26917"
#' covs<- dsm.covs(keene)


dsm.covs<- function(x, wdir=NULL, cpus=NULL){

  if(is.null(wdir)){temp<- tempdir()
  }else{temp<- wdir}

  # create directory for the rasters and a temp directory for intermediate files
  dir.create(paste(temp,"/covs", sep=""))
  dir.create(paste(temp,"/covs/tmp", sep=""))
  wd<- paste(temp,"/covs", sep="")
  temp_wd<- paste(temp,"/covs/tmp", sep="")
  terra::writeRaster(x,paste0(temp_wd,"/elev_in.tif", sep=""))

  DEM<- paste0(temp_wd,"/elev_in.tif", sep="")

  ####################################################################################################################
  ####################################################################################################################
  # Smooth the original DEM for topographic covariate creation
  ####################################################################################################################
  ####################################################################################################################

  whitebox::wbt_fill_missing_data(input = DEM,
                                  output = paste0(temp_wd,"/01FMD.tif", sep=""),
                                  filter = 11,
                                  weight = 2,
                                  verbose_mode = FALSE)

  whitebox::wbt_feature_preserving_smoothing(dem = paste0(temp_wd,"/01FMD.tif", sep=""),
                                             output = paste0(temp_wd,"/02SMTH.tif", sep=""),
                                             filter=11,
                                             norm_diff = 15,
                                             max_diff = 0.5,
                                             num_iter = 5,
                                             verbose_mode = FALSE)

  whitebox::wbt_fill_single_cell_pits(dem = paste0(temp_wd,"/02SMTH.tif", sep=""),
                                      output = paste0(temp_wd,"/03PITS.tif", sep=""),
                                      verbose_mode= FALSE)

  whitebox::wbt_breach_depressions_least_cost(dem = paste0(temp_wd,"/03PITS.tif", sep=""),
                                              output = paste0(temp_wd,"/04BRCH.tif", sep=""),
                                              dist = 10,
                                              max_cost = NULL,
                                              min_dist = TRUE,
                                              flat_increment = NULL,
                                              fill = TRUE ,
                                              verbose_mode = FALSE)


  ####################################################################################################################
  ####################################################################################################################
  # Now we are ready for RSAGA and derivatives
  ####################################################################################################################
  ####################################################################################################################

  DEM = paste0(temp_wd,"/03PITS.tif", sep="")
  DEMw = paste0(temp_wd,"/04BRCH.tif", sep="")
  dir.create(paste(temp_wd, "/sgrds", sep=""))
  sgrds<- paste(temp_wd, "/sgrds", sep="")

  ##Need to set a slope threshold based on the MRVBF publication by Gallant et al. using their equation based on raster resolution
  reso<-as.data.frame(terra::res(terra::rast(DEM)))[1,1]
  slopethreshold<- round(116.57*(reso^-0.62),0)

  # we create an object that holds the RSAGA environment variables that we call in when running all the tools
  work_env <- RSAGA::rsaga.env()

  # if user specifies cpus to be used, assign this to the RSAGA env file
  if(!is.null(cpus)){work_env$cores<- cpus}

  #Import DEM to SGRDs folder so it will get clipped later on
  RSAGA::rsaga.import.gdal(DEM, paste0(sgrds,'/dem.sgrd', sep=""), env = work_env, show.output.on.console = FALSE)

  ##################### Geoprocessing with RSAGA #######################

  ## Run slope, aspect, curvature tool
  ### Outputs Slope, Aspect, General Curvature, Profile Curvature, Plan Curvature, and Total Curvature
  RSAGA::rsaga.geoprocessor('ta_morphometry', 0, env = work_env,
                            param = list(ELEVATION = DEM,
                                         SLOPE = paste0(sgrds,"/sloper", sep=""),
                                         ASPECT = paste0(sgrds,"/asp", sep=""),
                                         C_GENE = paste0(sgrds,"/gcurv", sep=""),
                                         C_PROF = paste0(sgrds,"/pro", sep=""),
                                         C_PLAN = paste0(sgrds,"/plan", sep=""),
                                         C_TOTA = paste0(sgrds,"/tcurv", sep=""),
                                         UNIT_ASPECT = 0),
                            show.output.on.console = FALSE)
  ## Run MRVBF and MRRTF
  RSAGA::rsaga.geoprocessor('ta_morphometry', 8, env = work_env,
                            param = list(DEM = DEM,
                                         T_SLOPE = slopethreshold,
                                         MRVBF = paste0(sgrds,"/mrvbf", sep=""),
                                         MRRTF = paste0(sgrds,"/mrrtf", sep="")),
                            show.output.on.console = FALSE)

  ## Run Relative Heights and Slope Positions
  RSAGA::rsaga.geoprocessor('ta_morphometry', 14, env = work_env,
                            param = list(DEM = DEM,
                                         HO = paste0(sgrds,"/slopeh", sep=""),
                                         NH = paste0(sgrds,"/normh", sep=""),
                                         SH = paste0(sgrds,"/stanh", sep=""),
                                         MS = paste0(sgrds,"/msp", sep=""),
                                         HU = paste0(sgrds,"/vdepth", sep=""),
                                         W = 0.5,
                                         T = 10,
                                         E = 2),
                            show.output.on.console = FALSE)

  ## Run Terrain Roughness Index
  RSAGA::rsaga.geoprocessor('ta_morphometry', 16, env = work_env,
                            param = list(DEM = DEM,
                                         TRI = paste0(sgrds,"/tri", sep="")),
                            show.output.on.console = FALSE)

  ## Run Topographic Position Index
  RSAGA::rsaga.geoprocessor('ta_morphometry', 18, env = work_env,
                            param = list(DEM = DEM,
                                         TPI = paste0(sgrds,"/tpi", sep="")),
                            show.output.on.console = FALSE)

  ## Run Multi-Scale Topographic Position Index
  RSAGA::rsaga.geoprocessor('ta_morphometry', 28, env = work_env,
                            param = list(DEM = DEM,
                                         TPI = paste0(sgrds,"/mstpi", sep=""),
                                         SCALE_MIN = 1,
                                         SCALE_MAX = 8,
                                         SCALE_NUM = 2),
                            show.output.on.console = FALSE)

  ## Sky View Factor
  RSAGA::rsaga.geoprocessor('ta_lighting', 3, env = work_env,
                            param = list(DEM = DEM,
                                         VISIBLE = paste0(sgrds,"/vis", sep=""),
                                         SVF=paste0(sgrds,"/svf", sep="")),
                            show.output.on.console = FALSE)

  ## Run Convergence Index
  RSAGA::rsaga.geoprocessor('ta_morphometry', 1, env = work_env,
                            param = list(ELEVATION = DEM,
                                         RESULT = paste0(sgrds,"/conv", sep=""),
                                         METHOD = 0,
                                         NEIGHBOURS = 0),
                            show.output.on.console = FALSE)

  ## Run Analytical Hillshading
  RSAGA::rsaga.geoprocessor('ta_lighting', 0, env = work_env,
                            param = list(ELEVATION = DEM,
                                         SHADE = paste0(sgrds,"/hill", sep=""),
                                         METHOD = 0,
                                         AZIMUTH = 315,
                                         DECLINATION = 45,
                                         EXAGGERATION = 4,
                                         SHADOW = 1,
                                         NDIRS = 8),
                            show.output.on.console = FALSE)

  ##Create catchment areas
  RSAGA::rsaga.geoprocessor('ta_hydrology', 1, env = work_env,
                            param = list(ELEVATION = DEMw,
                                         FLOW = paste0(sgrds,"/catch", sep=""),
                                         METHOD = 1),
                            show.output.on.console = FALSE)

  ## Run Saga Wetness Index
  RSAGA::rsaga.geoprocessor('ta_hydrology', 15, env = work_env,
                            param = list(DEM = DEMw,
                                         AREA = paste0(sgrds,"/catch.sgrd", sep=""),
                                         TWI = paste0(sgrds,"/swi", sep=""),
                                         SUCTION=10,
                                         AREA_TYPE = 1,
                                         SLOPE_TYPE = 1),
                            show.output.on.console = FALSE)

  ## Run Topographic Wetness Index
  RSAGA::rsaga.geoprocessor('ta_hydrology', 20, env = work_env,
                            param = list(SLOPE = paste0(sgrds,"/sloper.sgrd", sep=""),
                                         AREA = paste0(sgrds,"/catch.sgrd", sep=""),
                                         TWI = paste0(sgrds,"/twi", sep=""),
                                         CONV = 0,
                                         METHOD = 0),
                            show.output.on.console = FALSE)

  ## Run Relative Slope Position
  RSAGA::rsaga.geoprocessor('ta_compound', 0, env = work_env,
                            param = list(ELEVATION = DEM,
                                         RSP = paste0(sgrds,"/rsp", sep="")),
                            show.output.on.console = FALSE)

  ## Run Stream Power Index
  RSAGA::rsaga.geoprocessor('ta_hydrology', 21, env = work_env,
                            param = list(SLOPE = paste0(sgrds,"/sloper.sgrd", sep=""),
                                         AREA = paste0(sgrds,"/catch.sgrd", sep=""),
                                         SPI = paste0(sgrds,"/spi", sep="")),
                            show.output.on.console = FALSE)

  ## Run Slope Length
  RSAGA::rsaga.geoprocessor('ta_hydrology', 7, env = work_env,
                            param = list(DEM = DEM,
                                         LENGTH = paste0(sgrds,"/len", sep="")),
                            show.output.on.console = FALSE)

  ## Run LS Factor
  RSAGA::rsaga.geoprocessor('ta_hydrology', 22, env = work_env,
                            param = list(SLOPE = paste0(sgrds,"/sloper.sgrd", sep=""),
                                         AREA = paste0(sgrds,"/catch.sgrd", sep=""),
                                         LS = paste0(sgrds,"/ls", sep=""),
                                         METHOD = 1,
                                         EROSIVITY = 1,
                                         STABILITY = 0),
                            show.output.on.console = FALSE)

  ######## Convert sgrds to tifs #############
  sgrd_stack <- terra::rast(list.files(path=sgrds, pattern="*.sdat$", full.names=T, recursive=FALSE))

  ## Write tiffs
  terra::writeRaster(sgrd_stack,
                     filename = paste0(wd,"/", names(sgrd_stack), '.tif'),
                     filetype="GTiff",
                     overwrite=TRUE)

  ###Remove SGRDS
  unlink(sgrds,recursive = TRUE)
  rm(sgrd_stack)

  #######################################################################################################################
  # Calculate northness and eastness from the aspect data
  #######################################################################################################################
  asp<- terra::rast(paste0(wd,"/asp.tif", sep=""))

  northness<- cos(asp)
  names(northness)<-"northness"

  terra::writeRaster(northness,filename = paste(wd,"/northness.tif", sep=""),
                                               filetype="GTiff",
                                               overwrite=TRUE)

  eastness<- sin(asp)
  names(eastness)<-"eastness"

  terra::writeRaster(eastness,filename = paste(wd,"/eastness.tif", sep=""),
                                               filetype="GTiff",
                                               overwrite=TRUE)
  rm(asp, northness, eastness)
  file.remove(paste0(wd,'/asp.tif', sep=""))

  #######################################################################################################################
  #######################################################################################################################
  # Now we can create the additional WhiteboxTools covariates
  #######################################################################################################################
  #######################################################################################################################

  # Deviation from Mean Elevation should be run 4 times at various filter sizes, starting with the default of 11
  # other filters are based on the cell width or height of the raster, whichever is larger
  size<- max(c(terra::ncol(terra::rast(DEM)), terra::nrow(terra::rast(DEM))))
  f1<- 11
  f2<- round(size*0.1, 0)
  f3<- round(size*0.3, 0)
  f4<- round(size*0.6, 0)
  filters<- c(f1,f2,f3,f4)

  # Difference from Mean Elevation should be run 4 times at various filter sizes, starting with the default of 11
  # other filters are based on the cell width or height of the raster, whichever is larger

  for(i in filters){
    whitebox::wbt_diff_from_mean_elev(dem = DEM,
                                      output = paste0(wd,"/dime",i,".tif", sep=""),
                                      filterx=i,
                                      filtery=i)
  }

  for(i in filters){
    whitebox::wbt_dev_from_mean_elev(dem = DEM,
                                     output = paste0(wd,"/deme",i,".tif", sep=""),
                                     filterx=i,
                                     filtery=i)
  }

  # Maximum Elevation Deviation should be run 3 times at various filter ranges: 3 to f2, f2 to f3 and f3 to f4

  whitebox::wbt_max_elevation_deviation(dem= DEM,
                                        out_mag= paste0(wd,"/med",f2,".tif", sep=""),
                                        out_scale= paste0(wd,"/meds",f2,".tif", sep=""),
                                        min_scale= 3,
                                        max_scale=f2,
                                        step= 2)

  whitebox::wbt_max_elevation_deviation(dem= DEM,
                                        out_mag= paste0(wd,"/med",f3,".tif", sep=""),
                                        out_scale= paste0(wd,"/meds",f3,".tif", sep=""),
                                        min_scale= f2,
                                        max_scale=f3,
                                        step= 2)

  whitebox::wbt_max_elevation_deviation(dem= DEM,
                                        out_mag= paste0(wd,"/med",f4,".tif", sep=""),
                                        out_scale= paste0(wd,"/meds",f4,".tif", sep=""),
                                        min_scale= f3,
                                        max_scale=f4,
                                        step= 2)

  # Maximum Difference from Mean Elevation should be run 3 times at various filter ranges: 3 to f2, f2 to f3 and f3 to f4

  whitebox::wbt_max_difference_from_mean(dem= DEM,
                                         out_mag= paste0(wd,"/mdm",f2,".tif", sep=""),
                                         out_scale= paste0(wd,"/mdms",f2,".tif", sep=""),
                                         min_scale= 3,
                                         max_scale=f2,
                                         step= 1)

  whitebox::wbt_max_difference_from_mean(dem= DEM,
                                         out_mag= paste0(wd,"/mdm",f3,".tif", sep=""),
                                         out_scale= paste0(wd,"/mdms",f3,".tif", sep=""),
                                         min_scale= f2,
                                         max_scale=f3,
                                         step= 1)

  whitebox::wbt_max_difference_from_mean(dem= DEM, out_mag= paste0(wd,"/mdm",f4,".tif", sep=""),
                                         out_scale= paste0(wd,"/mdms",f4,".tif", sep=""),
                                         min_scale= f3,
                                         max_scale=f4,
                                         step= 1)


  # Elevation Percentile

  for(i in filters){
    whitebox::wbt_elev_percentile(dem = DEM,
                                  output = paste0(wd,"/ep",i,".tif", sep=""),
                                  filterx=i,
                                  filtery=i,
                                  sig_digits = 2)}

  # Impoundment Size Index

  whitebox::wbt_impoundment_size_index(dem = DEM,
                                       out_mean = paste0(wd,"/isi_mean.tif", sep=""),
                                       out_max = paste0(wd,"/isi_max.tif", sep=""),
                                       damlength = 10)


  #Remove temporary files
  unlink(temp_wd, recursive = TRUE)

  covs<- terra::rast(list.files(path=wd, pattern=".tif$", full.names=T, recursive=FALSE))
  return(covs)

}
