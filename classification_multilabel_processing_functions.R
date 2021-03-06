####################################  Multilabel and fuzzy classification  #######################################
###########################################  Processing and Analyses  #######################################
#This script explores the fuzzy and multilabels concepts using classified land cover maps.
#This is the script with the functions.
#
#AUTHORS: Hichem Omrani and Benoit Parmentier                                             
#DATE CREATED: 11/03/2015 
#DATE MODIFIED: 03/14/2017
#Version: 1
#PROJECT: Multilabel and fuzzy experiment            

#
#COMMENTS: -  
#          - 
#TO DO:
# - Add function to reclassify inputs from aggregation
# - Add function to classify fuzzy/soft inputs
# - 
#
#################################################################################################
#what we need is clear step of what we are doing here..

#PART 1: Generate multi-label data from fine resolution
#1)Aggregate data at coarser resolution
#2)Reclassify pixels into multi and mono labels based on proportions?

###Loading R library and packages                                                      

library(raster)                 # loading the raster package
library(gtools)                 # loading R helper programming tools/functions
library(sp)                     # spatial objects in R
library(gplots)                 # plotting functions such as plotCI
library(rgdal)                  # gdal driver for R
library(RColorBrewer)           # color scheme, palettes used for plotting
library(gdata)                  # read different format (including .xlsx)
library(plotrix)                # plot options and functions 
library(rasterVis)              # raster visualization
library(colorRamps)             # contains matlab.like palette
library(zoo)                    # time series objects and methods
library(maptools)               #
library(rgeos)                  # spatial analysis, topological and geometric operations e.g. interesect, union, contain etc.
library(parallel)               # mclapply with cores...

###### Functions used in this script sourced from other files

#function_multilabel_fuzzy_analyses <- "classification_multilabel_processing_functions_01122016.R" #PARAM 1
#script_path <- "/home/bparmentier/Google Drive/LISER_Lux/R_scripts" #path to script #PARAM 2
#source(file.path(script_path,function_multilabel_fuzzy_analyses)) #source all functions used in this script 1.

##### Functions used in this script 

create_dir_fun <- function(outDir,out_suffix){
  #if out_suffix is not null then append out_suffix string
  if(!is.null(out_suffix)){
    out_name <- paste("output_",out_suffix,sep="")
    outDir <- file.path(outDir,out_name)
  }
  #create if does not exists
  if(!file.exists(outDir)){
    dir.create(outDir)
  }
  return(outDir)
}

#Used to load RData object saved within the functions produced.
load_obj <- function(f){
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

#This function creates a spatial polygon data frame object for the extent matching a raster input
create_polygon_from_extent<-function(reg_ref_rast,outDir=NULL,outSuffix=NULL){
  #This functions returns polygon sp from input rast
  #Input Arguments: 
  #reg_ref_rast: input ref rast
  #outDir : output directory, if NULL then the current dir in used
  #outSuffix: output suffix used for the naming of the shapefile
  #Output: 
  #reg_outline_poly: spatial polygon data.frame
  #
  if(is.null(outDir)){
    outDir=getwd()
  }
  if(is.null(outSuffix)){
    outSuffix=""
  }
  ref_e <- extent(reg_ref_rast) #extract extent from raster object
  reg_outline_poly <- as(ref_e, "SpatialPolygons") #coerce raster extent object to SpatialPolygons from sp package 
  reg_outline_poly <- as(reg_outline_poly, "SpatialPolygonsDataFrame") #promote to spdf
  proj4string(reg_outline_poly) <- projection(reg_ref_rast) #Assign projection to spdf
  infile_reg_outline <- paste("reg_out_line_",out_suffix,".shp",sep="") #name of newly crated shapefile with the extent
  writeOGR(reg_outline_poly,dsn= outDir,layer= sub(".shp","",infile_reg_outline), 
           driver="ESRI Shapefile",overwrite_layer="TRUE")
  
  return(reg_outline_poly) #return spdf
}


#Function to aggregate from fine to coarse resolution, this will change accordingly once the input raster ref is given..
#
aggregate_raster <- function(r_in, agg_fact, reg_ref_rast=NULL,agg_fun="mean",out_suffix=NULL,file_format=".tif",out_dir=NULL,out_rast_name=NULL){
  #Aggregate raster from raster input and reference file
  #INPUT arguments:
  #1) r_in: input raster layer
  #2) agg_fact: factor to aggregate
  #3) reg_ref_rast: reference raster to match in resolution, if NULL then send a message
  #4) agg_fun: default is mean
  #5) out_suffix: output suffix
  #6) file_Format: raster format used e.g. .tif
  #7) out_dir: output directory
  #8) out_rast_name: output raster name if null it is created from the input file name
  #OUTPUT:
  # out_raster_name: name of the file containing the aggregated raster
  #
  # Authors: Benoit Parmentier
  # Created: 10/15/2015
  # Modified: 05/10/2016
  # To Do: 
  # - Add option to disaggregate
  #
  ################################
  
  ##If file is provided rather than RasterLayer
  if(class(r_in)!="RasterLayer"){
    r_in <- raster(r_in)
  }
  
  if(is.null(agg_fact)){
    res_ref <- res(reg_ref_rast)[1] #assumes square cells, and decimal degrees from WGS84 for now...
    res_in <- res(r_in)[1] #input resolution, assumes decimal degrees
    agg_fact <-round(res_ref/res_in) #find the factor needed..
    #fix this to add other otpions e.g. aggregating down
  }
  
  #Default values...
  if(is.null(out_suffix)){
    out_suffix <- ""
  }
  
  if(is.null(out_dir)){
    out_dir <- "."
  }
  
  ## Create output raster name if out_rast_name is null
  if(is.null(out_rast_name)){
    raster_name <- filename(r_in)
    extension_str <- extension(raster_name)
    raster_name_tmp <- gsub(extension_str,"",basename(raster_name))
    out_rast_name <- file.path(out_dir,paste("agg_",agg_fact,"_",raster_name_tmp,out_suffix,file_format,sep="")) #output name for aster file
    #out_rast_name <- raster_name #for use in function later...
  }

  r_agg <- aggregate(r_in, fact=agg_fact,FUN=agg_fun,filename=out_rast_name,overwrite=TRUE)
  
  return(out_rast_name)
  
}

### This is a very general function to process raster to match a raster of reference,
create__m_raster_region <-function(j,list_param){
  #This processes a list of raster to match to a region of interest defined by a reference raster
  #INPUT Arguments: raster name of the file,reference file with
  # j: file to be processed with input parameters
  # raster_name: list of raster to process i.e. match the region of interest
  # reg_ref_rast: reference raster used to defined the region of interest and spatial parameters
  # out_rast_name: output raster name, if NULL then use out_suffix to the input name to generate output name
  # agg_param: aggregation parameters: this is a vector used in in the aggregate function. It has three items:
  #                                     -TRUE/FALSE: if true then aggregate
  #                                     -agg_fact: aggregation factor, if NULL compute on the fly
  #                                     -agg_fun: aggregation function to use, the default is mean
  # file_format: output format used in the raster e.g. .tif, .rst
  # NA_flag_val: flag value used for no data
  # input_proj_str: defined projection,default null in which case it is extract from the input raster
  # out_suffix : output suffix added to output names if no output raster name is given
  # out_dir:  <- list_param$out_dir
  # Output: spatial grid data frame of the subset of tiles
  #
  # Authors: Benoit Parmentier
  # Created: 10/01/2015
  # Modified: 11/01/2015
  #TODO:
  # - Add option to disaggregate...
  # - Modify agg param to be able to use different ones by file j for the mcapply function
  #
  ################################################
  ## Parse input arguments
  raster_name <- list_param$raster_name[[j]] #list of raster ot project and crop, this is a list!!
  reg_ref_rast <- list_param$reg_ref_rast #This must have a coordinate system defined!!
  out_rast_name <- list_param$out_rast_name[j] #if NULL then use out_suffix to add to output name
  agg_param <- list_param$agg_param #TRUE,agg_fact,agg_fun
  file_format <- list_param$file_format #.tif, .rst
  NA_flag_val <- list_param$NA_flag_val #flag value used for no data
  input_proj_str <- list_param$input_proj_str #default null?
  out_suffix <- list_param$out_suffix
  out_dir <- list_param$out_dir
  
  ## Start #
  
  ## Create raster object if not already present
  if(class(raster_name)!="RasterLayer"){
    layer_rast<-raster(raster_name)
  }else{
    layer_rast <- raster_name
    raster_name <- filename(layer_rast)
  }
  
  ## Create output raster name if out_rast_name is null
  if(is.null(out_rast_name)){
    extension_str <- extension(raster_name)
    raster_name_tmp <- gsub(extension_str,"",basename(raster_name))
    out_rast_name <- file.path(out_dir,paste(raster_name_tmp,"_crop_proj_reg_",out_suffix,file_format,sep="")) #for use in function later...
  }
  
  ## Get the input raster projection information if needed
  if(is.null(input_proj_str)){
    input_proj_str <-projection(layer_rast)   #Extract current coordinates reference system in PROJ4 format
  }else{
    projection(layer_rast) <- input_proj_str #assign projection info
  }
  region_temp_projected <- projectExtent(reg_ref_rast,CRS(input_proj_str))     #Project from ref to current region coord. system
  
  layer_crop_rast <- crop(layer_rast, region_temp_projected) #crop using the extent from the region tile
  #layer_projected_rast<-projectRaster(from=layer_crop_rast,crs=proj4string(reg_outline),method="ngb")
  if(agg_param[1]==TRUE){
    agg_fact <- as.numeric(agg_param[2]) #in case we have a string/char type
    agg_fun <- agg_param[3]
    #debug(aggregate_raster)
    r_agg_raster_name <- aggregate_raster(r_in=layer_crop_rast, #raster to be aggregated
                                          reg_ref_rast, #reference raster with the desired resolution
                                          agg_fact=agg_fact, #given aggregation factor
                                          agg_fun="mean", #aggregation function
                                          out_suffix=out_suffix,
                                          file_format=".tif",
                                          out_dir=out_dir,
                                          out_rast_name = NULL)
    layer_crop_rast <- raster(r_agg_raster_name)
  }
  
  layer_projected_rast <- projectRaster(from=layer_crop_rast,to=reg_ref_rast,method="ngb",
                                        filename=out_rast_name,overwrite=TRUE)
  
  #Need cleanup of tmp files here!!! building up to 19gb!
  removeTmpFiles(h=0)
  
  return(out_rast_name)
}

plot_to_file <- function(raster_name,res_pix=480,out_suffix=NULL,out_dir=NULL){
  #Quick utility function to plot raster to png file for the workflow
  #This is useful to visually check the outputs from the workflow.
  #INPUT arguments:
  #raster_name: input raster object or file name of the raster object
  #out_suffix: output suffix
  #out_dir: output directory
  #OUTPUT:
  # png_file: name of the file containing the figure/plot
  #
  # Authors: Benoit Parmentier
  # Created: 11/01/2015
  # Modified: 11/02/2015
  # To Do: 
  # - Add option to choose plot format e.g. jpeg, tif etc.
  #
  ################################
  ## Create raster object if not already present
  if(class(raster_name)!="RasterLayer"){
    layer_rast<-raster(raster_name)
  }else{
    layer_rast <- raster_name
    raster_name <- filename(layer_rast)
  }
  
  if(is.null(out_suffix)){
    out_suffix <- ""
  }
  if(is.null(out_dir)){
    out_dir <- getwd()
  }
  
  #Extract name
  extension_str <- extension(raster_name)
  raster_name_tmp <- gsub(extension_str,"",basename(raster_name))
  
  res_pix <- 480
  col_mfrow <- 1
  row_mfrow <- 1
  #Might change file format later...
  png_file <- file.path(out_dir,paste("figure_",raster_name_tmp,"_",out_suffix,".png",sep="")) #may be changed later to other format
  png(filename=png_file,
      width=col_mfrow*res_pix,height=row_mfrow*res_pix)
  plot(layer_rast)
  title(raster_name_tmp)
  dev.off()
  return(png_file)
}

generate_soft_cat_aggregated_raster_fun <- function(r,reg_ref_rast,agg_fact,agg_fun,num_cores,NA_flag_val,file_format,out_dir,out_suffix){
  ## Function to aggregate categories
  
  ##INPUTS
  #1) r: raster to aggregate and crop/project if necessary
  #2) reg_ref_rast: reference raster, it must have a coordinate system defined
  #3) agg_fact:
  #4) agg_fun:
  #5) num_cores: number of cores used in the parallel processsing
  #6) NA_flag_val: flag value used for no data
  #7) file_format: raster file format e.g. .tif, .rst
  #8) out_dir: output directory
  #9) out_suffix: output suffix added to files
  #
  #OUTPUTS
  #
  #
  #
  
  
  ##### STEP 1: Check input
  
  if(class(r)!="RasterLayer"){
    r <- raster(r)
  }
  
  NAvalue(r) <- NA_flag_val #make sure we have a flag value assigned
  
  ###### STEP 1: BREAK OUT
  ## Breakout layers: separate categories in individual layers
  
  freq_tb <- as.data.frame(freq(r)) #zero is NA?
  ## get the names
  names_val <- freq_tb$value
  names_val <- names_val[!is.na(names_val)] #remove NA

  ## Make a brick composed of multiple layers: one layer per category (in one unique multiband file)
  out_raster_name <- file.path(out_dir,paste("r_layerized_bool_",out_suffix,file_format,sep=""))
  r_layerized <- layerize(r,
                          classes=names_val,
                          filename= out_raster_name,
                          overwrite=T)
  list_out_raster_name <- file.path(out_dir,paste("r_layerized_bool_",names_val,"_",out_suffix,file_format,sep=""))
  
  ## Now write out separate layers in one unique file (one per categories)
  #notice some issue here, looping through
  #for(i in 1:nlayers){
  #  writeRaster(r_layerized,
  #              bylayer=T,
  #              suffix=paste0(names_val,"_",out_suffix),
  #              filename=paste("r_layerized_bool",
  #                             file_format,sep="")
  #             ,overwrite=T)
  #}
  
  writeRaster(r_layerized,
              bylayer=T,
              suffix=paste0(names_val,"_",out_suffix),
              filename=paste("r_layerized_bool",
                             file_format,sep="")
              ,overwrite=T)

  lf_layerized_bool <- paste("r_layerized_bool_",names_val,"_",out_suffix,file_format,sep="")
  #names(r_layerized) <- 
  #inMemory(r_date_layerized)
  #filename(r_date_layerized)
  #r_test <- raster(lf_layerized_bool[1])
  
  ### STEP 2: aggregate
  ## Aggregate
  
  ##To do
  ##Set correct names based on categories of input
  ##Set output name
  
  #r_agg_fname <-aggregate_raster(agg_fact=agg_fact,
  #                               r_in=r_date_layerized,reg_ref_rast=NULL,agg_fun="mean",out_suffix=NULL,file_format=".tif",out_dir=NULL)
  #debug(aggregate_raster)
  #r_agg_fname <-aggregate_raster(r_in=raster(lf_layerized_bool[1]),
  #                               agg_fact=agg_fact,
  #                               reg_ref_rast=NULL,
  #                               #agg_fun="mean",
  #                               agg_fun=agg_fun,
  #                               out_suffix=NULL,
  #                               file_format=file_format,
  #                               out_dir=out_dir,
  #                               out_rast_name = NULL)
  #r_var_s <- mclapply(1:length(infile_var),FUN=import_list_modis_layers_fun,list_param=list_param_import_modis,mc.preschedule=FALSE,mc.cores = num_cores) #This is the end bracket from mclapply(...) statement
  
  lf_agg <- mclapply(lf_layerized_bool,
                     FUN=aggregate_raster,
                     #r_in=raster(lf_layerized_bool[1]),
                     agg_fact=agg_fact,
                     reg_ref_rast=NULL,
                     #agg_fun="mean",
                     agg_fun=agg_fun,
                     out_suffix=NULL,
                     file_format=file_format,
                     out_dir=out_dir,
                     out_rast_name = NULL,
                     mc.preschedule=FALSE,
                     mc.cores = num_cores) 
  #r_agg <- stack(unlist(lf_agg))
  raster_outname <- unlist(lf_agg)
  #apply function by layer using lapply.
  
  ## Reclassify by labels
  return(raster_outname)
}


############################### END OF SCRIPT ########################