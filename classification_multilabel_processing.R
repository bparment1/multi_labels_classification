####################################  Multilabel and fuzzy classification  #######################################
###########################################  Processing and Analyses  #######################################
#This script explores the fuzzy and multilabels concepts using classified land cover maps.

#AUTHORS: Hichem Omrani and Benoit Parmentier                                             
#DATE CREATED: 11/03/2015 
#DATE MODIFIED: 05/11/2016
#Version: 2
#PROJECT: Multilabel and fuzzy experiment            

#
#COMMENTS: -  
#          - 
#TO DO:
# - 
# - 
# - 
#
#################################################################################################
#what we need is clear step of what we are doing here..

#PART 1: Generate multi-label data from fine resolution
#0)Break out layers into individual categories: e.g. 3 cat : 1,2,3
#1)Aggregate data at coarser resolution for each categories
#2)Reclassify pixels into multi and mono labels based on proportions:
#  
#3)Repeat for date2
#PART 2: Modeling
#Use fuzzy (continuous) values from proportion by categories: predict using a method..e.g. NN,knn etc.
#Use reclassify layers by cat using multi-label algorithm
#
#PART3
#Compare results
#
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

function_multilabel_fuzzy_analyses <- "classification_multilabel_processing_functions_05102016.R" #PARAM 1
#classification_multilabel_processing_functions_05102016.R
script_path <- "/home/bparmentier/Google Drive/LISER_Lux/R_scripts" #path to script #PARAM 2
source(file.path(script_path,function_multilabel_fuzzy_analyses)) #source all functions used in this script 1.

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
#add reclassify function..
#add aggregate function

#####  Parameters and argument set up ###########

#in_dir <- "/home/bparmentier/Google Drive/LISER_Lux/Hichem" #local bpy50
in_dir <- "/home/bparmentier/Google Drive/LISER_Lux/Data-USA-with1m-resolution" #local bpy50, MA data
#in_dir <- "//crc/profiles/RedirectFolders/hichem/Desktop/LISER/papers/Fuzzy CA/raster-USA/Hichem"#LISER
out_dir <- "/home/bparmentier/Google Drive/LISER_Lux/"
CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84 # CONST 2
proj_str<- CRS_WGS84 
CRS_reg <- CRS_WGS84 # PARAM 4

r_date1_fname <- "/home/bparmentier/Google Drive/LISER_Lux/Data-USA-with1m-resolution/lyn_1971_landuse.rst"
r_date2_fname <- "/home/bparmentier/Google Drive/LISER_Lux/Data-USA-with1m-resolution/lyn_1999_landuse.rst"
r_date3_fname <- "/home/bparmentier/Google Drive/LISER_Lux/Data-USA-with1m-resolution/lyn_2005_landuse.rst"

data_type <- "Int32"
file_format <- ".tif" #PARAM5
NA_value <- -9999 #PARAM6
NA_flag_val <- NA_value #PARAM7
out_suffix <-"multilabel_05072016" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9
num_cores <- 4 #PARAM 14
#agg_param <- 5
agg_fun <- "mean" 
  
date1 <- 1971
date2 <- 1999
date3 <- 2005

agg_fact <- 5

################# START SCRIPT ###############################

### PART I READ AND PREPARE DATA FOR REGRESSIONS #######
#set up the working directory
#Create output directory

#out_dir <- in_dir #output will be created in the input dir
out_suffix_s <- out_suffix #Can modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix_s)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

r_date1 <- raster(r_date1_fname)
r_date2 <- raster(r_date2_fname) 
r_date3 <- raster(r_date3_fname) 
  
lf <- c(r_date1_fname,r_date2_fname,r_date3_fname)

plot(r_date1,main=paste("Date 1: ",date1,sep="")) 
plot(r_date1,main=paste("Date 2: ",date2,sep="")) 
plot(r_date1,main=paste("Date 3: ",date3,sep="")) 

r_stack <- stack(r_date1,r_date2,r_date3)
names(r_stack) <- as.character(c(date1,date2,date3))
plot(r_stack)

freq_tb2 <- freq(r_date2) #zero is NA?
freq_tb2 <- freq(r_date2) #zero is NA?
## create a table of each cat per date:
freq_tb <- freq(r_stack,merge=T)
write.table(freq_tb,file = paste0("classes_freq_table_",out_suffix,".txt",sep=""))
dim(r_date2)

#### PART 1: Processing data: generate multi-label #####################

##This function assumes that NA values are set...
## Makes this function use majority rules too?
#debug(generate_soft_cat_aggregated_raster_fun)
#test <- generate_soft_cat_aggregated_raster_fun(lf[1],reg_ref_rast=NULL,agg_fact,agg_fun,num_cores,NA_flag_val=0,file_format,out_dir,out_suffix)
  
#r_1971_agg5_soft <- stack(test)
#plot(r_1971_agg5_soft)

list_lf_agg_soft <- vector("list",length=length(lf)) 

###
for(i in 1:length(lf)){

  raster_name <- lf[i]
  extension_str <- extension(raster_name)
  raster_name_tmp <- gsub(extension_str,"",basename(raster_name))
  #out_rast_name <- file.path(out_dir,paste("agg_",agg_fact,"_",raster_name_tmp,out_suffix,file_format,sep="")) #output name for aster file
  
  out_suffix_s <- paste0(raster_name_tmp,"_" ,out_suffix)
  lf_agg_soft <- generate_soft_cat_aggregated_raster_fun(raster_name,reg_ref_rast=NULL,agg_fact,agg_fun,num_cores,NA_flag_val=0,file_format,out_dir,out_suffix_s)
  list_lf_agg_soft[[i]] <- lf_agg_soft
}


#generate_soft_cat_aggregated_raster_fun <- function(r,reg_ref_rast,agg_fact,agg_fun,NA_flag_val,file_format,out_dir,out_suffix){
  
#### PART 2: Modeling ###########################

NA_flag_val_tmp <- 0
NAvalue(r_stack) <- NA_flag_val_tmp

test <- crosstab(subset(r_stack,1),subset(r_stack,2)) # 306x3
test2 <- subset(test,test$Freq > 0) # 44x3
                 
## Let's say we want to focus on transition
## First do predction using hard classified maps using neural net and logistic??

## Second do prediction using soft classifified maps,proportions?

##Needs to select validation/testing and training!!

#### PART 3: Assessment ########################


############################### END OF SCRIPT ########################