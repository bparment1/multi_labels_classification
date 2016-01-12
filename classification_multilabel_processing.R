####################################  Multilabel and fuzzy classification  #######################################
###########################################  Processing and Analyses  #######################################
#This script explores the fuzzy and multilabels concepts using classified land cover maps.

#AUTHORS: Hichem Omrani and Benoit Parmentier                                             
#DATE CREATED: 11/03/2015 
#DATE MODIFIED: 01/12/2016
#Version: 1
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

###### Functions used in this script sourced from other files

#function_multilabel_fuzzy_analyses <- "classification_multilabel_processing_function_11242015.R" #PARAM 1
#script_path <- "/home/bparmentier/Google Drive/LISER_Lux/R_scripts" #path to script #PARAM 2
#source(file.path(script_path,function_rainfall_time_series_NEST_analyses)) #source all functions used in this script 1.

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

in_dir <- "/home/bparmentier/Google Drive/LISER_Lux/Hichem" #local bpy50
#in_dir <- "//crc/profiles/RedirectFolders/hichem/Desktop/LISER/papers/Fuzzy CA/raster-USA/Hichem"#LISER
CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84 # CONST 2
proj_str<- CRS_WGS84 
CRS_reg <- CRS_WGS84 # PARAM 4

file_format <- ".rst" #PARAM5
NA_value <- -9999 #PARAM6
NA_flag_val <- NA_value #PARAM7
out_suffix <-"NEST_prism_11192015" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9
num_cores <- 11 #PARAM 14

################# START SCRIPT ###############################

### PART I READ AND PREPARE DATA FOR REGRESSIONS #######
#set up the working directory
#Create output directory

out_dir <- in_dir #output will be created in the input dir
out_suffix_s <- out_suffix #can modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix_s)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}


##########################
# Reading data from landuse layers
##########################

LU1978 <- file.path(in_dir,"landusebase.asc")
LU1998 <- file.path(in_dir, "landusefinal.asc")

##################### begin 
#########################################################################################################
raster1978 <- raster(LU1978) #There is no coordinate system!!!!
raster1998 <- raster(LU1998) 

plot(raster1978)
plot(raster1998)
plot(stack(raster1978,raster1998))
data_matrix1978 <- rasterToPoints(raster1978) #this step is not necessary...
head(data_matrix1978)

data_matrix1998 <- rasterToPoints(raster1998)
head(data_matrix1998)
hist(raster1978)
hist(raster1998)
freq(raster1998)
freq(raster1998)# 
#> freq(raster1998)
#value  count
#[1,]     0 114999
#[2,]     1  21621 #urban
#[3,]    NA  78260

#this computes the number of 1 in 5x5 pixels
r.agg_78 <- aggregate(raster1978, fact=5, fun=sum, na.rm=TRUE) # 10

r.agg_98 <- aggregate(raster1998, fact=5, fun=sum, na.rm=TRUE) #so if 25 then 100% urban, 
#could use mean!!
r_agg_98_perc <- aggregate(raster1998, fact=5, fun=mean, na.rm=TRUE)*100 #so 100% urban, 
r_agg_78_perc <- aggregate(raster1978, fact=5, fun=mean, na.rm=TRUE)*100 #so 100% urban, 

plot(r.agg_78)
plot(r.agg_98)
plot(r_agg_98_perc)
plot(r_agg_78_perc)

histogram(r.agg_98)
histogram(raster1998)

data.agg_78 <- rasterToPoints(r.agg_78) #? not sure why we are doing this right now
head(data.agg_78)

data.agg_98 <- rasterToPoints(r.agg_98) #this is matrix...
head(data.agg_98)
colnames(data.agg_78) #need comments!

r78 = data.agg_78; #make a copy?
r98 = data.agg_98;

### Reclassify data, column 3 is the landuse final
#if 25 then all the pixels are of the same class other not (ML)
#this should be done directly in the raster package classify function

#--> reclassify mixed pixels into mixed and unmixed (mono and multi) for date 78
for (i in 1:nrow(data.agg_78)){
  if ( (data.agg_78[i,3]!=0) && (data.agg_78[i,3]!=25) ){
    r78[i,3] = 2 # it means mixed class ok
  } # ML 
}

#--> reclassify mixed pixels into mixed and unmixed (mono and multi) for date 98
for (i in 1:nrow(data.agg_98)){
if ( (data.agg_98[i,3]!=0) && (data.agg_98[i,3]!=25) ) # ML 
r98[i,3] = 2
}

## if col 3 is 25 then urban ok
for (i in 1:nrow(r78)){
if (r78[i,3]==25) # Urban 
r78[i,3] = 1
}

#reclassify pixel to urban if 25 for year 98
for (i in 1:nrow(r98)){
if (r98[i,3]==25) # Urban 
r98[i,3] = 1
}

#Reclassification using raster!!
#2: urban and non urban mixed
#1: urban
#0: non urban

recmat_val <- c(-1, 0, 0,  
                0, 99, 2,  
                99, 100, 1)
rclmat <- matrix(recmat_val, ncol=3, byrow=TRUE)
r98_rec <- reclassify(r_agg_98_perc, rclmat)
freq(r98_rec)
col_pal <- c("grey","red","black")

plot(r98_rec,col=col_pal,legend=F)
cat_names <- c("non urban","ML","urban")
legend("topright",legend=cat_names,title="Categories",
       pt.cex=1,cex=1,fill=col_pal,bty="n")

### Bring it back to the raster format (use reclassify for all)
rt1 <- rasterFromXYZ(r78)
rt2 <- rasterFromXYZ(r98)

#reclassify(dd)
r_agg_98_perc <- aggregate(raster1998, fact=5, fun=mean, na.rm=TRUE)*100 #so 100% urban, 
r_agg_78_perc <- aggregate(raster1978, fact=5, fun=mean, na.rm=TRUE)*100 #so 100% urban, 


#Add legend
plot(rt1)
freq(rt1)
histogram(rt1)
writeRaster(rt1, "agg78_500meter.asc", overwrite=TRUE) #better not use asc !!! not standard file 
writeRaster(rt2, "agg98_500meter.asc", overwrite=TRUE) #note that the coordinate system is not assigned!!!

LU = cbind(r78, r98[,3]) # NU(0), U(2), ML(2), recreate a matrix object, use data.frame instead...
colnames(LU) = c("X", "Y", "LU78", "class") # class refer to LU98, maybe use year instead of class?
#LU[1:5,]
head(LU)
LU_df <- as.data.frame(LU)
table(LU_df$class)
#> table(LU_df$class)
#
#0    1    2 
#3190  168 2224 

#### GENERATE ML land use DATA

data_matrix = cbind(data_matrix1978, data_matrix1998[,3])
agg_matrix = cbind(data.agg_78, data.agg_98[,3])

head(data_matrix)
head(agg_matrix)
dim(data_matrix)
dim(data_matrix1978)

table(data_matrix[,3])
#0      1 

#117012  19608 

################Not clear what is happening here...
##### Use reclassifying function
### for the year 1978
data.agg_78_NOT_ML_equal4 = data.agg_78[data.agg_78[,3]== 0,]  
data.agg_78_NOT_ML_equal4[,3]=0

data.agg_78_NOT_ML_equal8 = data.agg_78[data.agg_78[,3]== 25,]
data.agg_78_NOT_ML_equal8[,3]=25
data.agg_78_NOT_ML = rbind(data.agg_78_NOT_ML_equal4, data.agg_78_NOT_ML_equal8) # not ML 

data.agg_78_ML_not4 = data.agg_78[data.agg_78[,3] != 0,]
data.agg_78_ML_not4and8 = data.agg_78_ML_not4[data.agg_78_ML_not4[,3] != 25,] # ML 
####

### for the year 1998
data.agg_98_NOT_ML_equal4 = data.agg_98[data.agg_98[,3]== 0,]  
data.agg_98_NOT_ML_equal4[,3]=0

data.agg_98_NOT_ML_equal8 = data.agg_98[data.agg_98[,3]== 25,]
data.agg_98_NOT_ML_equal8[,3]=25
data.agg_98_NOT_ML = rbind(data.agg_98_NOT_ML_equal4, data.agg_98_NOT_ML_equal8) # not ML 

data.agg_98_ML_not4 = data.agg_98[data.agg_98[,3] != 0,]
data.agg_98_ML_not4and8 = data.agg_98_ML_not4[data.agg_98_ML_not4[,3] != 25,] # ML 
####
 
# data_mat = cbind(data_matrix, real_change)
# colnames(data_mat) = c("X", "Y", "y78", "y98", "change01")

# urban cells are coded as 1 (red) in 1978 and 1998
# and the rest of the map coded as 0 (gray) 
real_change = rep(0,nrow(data_matrix),1)

for ( i in 1:nrow(data_matrix) ) {
 	if ( (data_matrix[i,3]==2) && (data_matrix[i,4]==1) ) { # 2== NU and 1==U
			real_change[i]=1 }
	else {
	real_change[i]=0
} }

###### BEGIN mapping 

data_mat = cbind(data_matrix, real_change)
colnames(data_mat) = c("X", "Y", "y78", "y98", "change01")

color = c("black", "gray")
colorCH = c("white", "black")
colorCH_ML = c("white", "black", "darkgreen")

########################################################################################################

# graphics.off()

# pdf("Fig_ML_USA.pdf", width=5.5, height=5.9,pointsize=12, horizontal=F)

# pdf("maps_USA.pdf")

#dat_mat is a matrix!! use standard raster...
par(mfrow=c(2,3)) #, mar=c(2,2,1,1), mgp=c(3,1,0)*0.7, lab=c(3,3,1))

plot(data_mat[,1], data_mat[,2], col=color[data_mat[,3]], cex=0.1, xlab="", ylab="", main = "Landuse base in 1978")
legend("topright", pch=16, c("NU", "U"), cex=1.2, col=c("gray", "black"))

plot(data_mat[,1], data_mat[,2], col=color[data_mat[,4]], cex=0.1, xlab="", ylab="", main = "Landuse final in 1998")
legend("topright", pch=16, c("NU", "U"), cex=1.2, col=c("gray", "black"))

# plot changes between original maps 
plot(data_mat[,1], data_mat[,2], col=colorCH[real_change+1], cex=0.1, xlab="", ylab="", main = "Real change - ml")
legend("topright", pch=16, c("New U"), cex=1.2, col=c("black"))

plot(data.agg_78_NOT_ML[,1], data.agg_78_NOT_ML[,2], col=color[data.agg_78_NOT_ML[,3]], xlab="", ylab="", main = "Landuse base - aggregated")
points(data.agg_78_ML_not4and8[,1], data.agg_78_ML_not4and8[,2], col="green", cex=0.1, type="p", xlab="", ylab="")
legend("topright", pch=16, c("NU", "U", "ML"), cex=1.2, col=c("gray", "black", "darkgreen"))

plot(data.agg_98_NOT_ML[,1], data.agg_98_NOT_ML[,2], col=color[data.agg_98_NOT_ML[,3]], pch=16, xlab="", ylab="", main = "Landuse final - aggregated")
points(data.agg_98_ML_not4and8[,1], data.agg_98_ML_not4and8[,2], col="green", cex=0.15, type="p", xlab="", ylab="")
legend("topright", pch=16, c("NU", "U", "ML"), cex=1.2, col=c("gray", "black", "darkgreen"))

# plot the changes between aggregated maps : 
plot(data_mat[,1], data_mat[,2], col=colorCH_ML[real_change_ML+1], cex=0.1, xlab="", ylab="", main = "Real change - ML")
legend("topright", pch=16, c("New U","New ML"), cex=1.2, col=c("black", "darkgreen"))

# dev.off()

###### END mapping 

# Changes in ML land use maps  

real_change_ML = rep(0,nrow(data.agg_98),1)

for ( i in 1:nrow(data.agg_98) ) {
 	if ( (data.agg_78[i,3]==0) && (data.agg_98[i,3]==25) )                                 { # 8== NU and 4==U; from NU to U: NEW U (1)
			real_change_ML[i]=1 }
	else if ( (data.agg_78[i,3] !=25) && (data.agg_78[i,3] !=0) && (data.agg_98[i,3] ==25)) { # from ML to U: NEW U(1)
	real_change_ML[i]=1 } 
	else if ( (data.agg_78[i,3] ==25) && (data.agg_98[i,3] !=0) && (data.agg_98[i,3] !=25)) { # from U to ML: NEW ML (2)
	real_change_ML[i]=2 } 
	else if ( (data.agg_78[i,3] ==0) && (data.agg_98[i,3] !=0) && (data.agg_98[i,3] !=25)) { # from NU to ML: NEW ML (2)
	real_change_ML[i]=2 } 
	else 
	real_change_ML[i]=0
}

# table(real_change_ML)
#    0     1     2 
# 31289   734  2309 

dt1ML = data.agg_78[(data.agg_78[,3]!=0) & (data.agg_78[,3]!=25),]
dt1_NU = data.agg_78[data.agg_78[,3]==0,]
dt1_U = data.agg_78[data.agg_78[,3]==25,]

dt2ML = data.agg_98[(data.agg_98[,3]!=0) & (data.agg_98[,3]!=25),]
#dt2_NU = data.     0      1 117012  19608 agg_98[data.agg_98[,3]==0,]
dt2_U = data.agg_98[data.agg_98[,3]==25,]

dt1_U[,3]=2
dt1ML[,3]=1
d1 = rbind(dt1_U, dt1_NU, dt1ML)

dt2_U[,3]=2
dt2ML[,3]=1
d2 = rbind(dt2_U, dt2_NU, dt2ML)

r1 <- rasterFromXYZ(d1)
r2 <- rasterFromXYZ(d2)

writeRaster(r1, "rrr1.asc")
writeRaster(r2, "rrr2.asc")

### small example to illustrate the concept of mono and multi-labels

mat = rbind(c(1, 1, 1, 1, 1, 2, 4, 6, 7),
	      c(1, 3, 3, 2, 5, 6, 6, 7, 8),
		c(1, 1, 3, 2, 2, 2, 4, 5, 6),
		c(1, 2, 2, 2, 2, 4, 4, 5, 6),
		c(1, 1, 1, 2, 2, 2, 4, 5, 6),
		c(1, 1, 1, 2, 2, 3, 4, 5, 6),
		c(1, 1, 1, 1, 1, 2, 3, 4, 5),
		c(0, 0, 1, 1, 1, 2, 4, 4, 5),
            c(0, 1, 1, 1, 1, 2, 3, 4, 4))

### aggregation using mono-label class assignment (by max function)

InRas <- raster(mat) 
writeRaster(InRas, file=file.path(out_dir,"InRas.rst"),overwrite=T)
plot(InRas)
freq(InRas)
OutRas_mono_label <- aggregate(InRas, fact=3, fun=max, na.rm=TRUE)
OutRas_mono_label.agg <- rasterToPoints(OutRas_mono_label)
head(OutRas_mono_label.agg)
writeRaster(OutRas_ml, "OutRas_mono_label.asc", overwrite=TRUE)

####

### aggregation using multi-label class assignment 

OutRas_multi_label <- aggregate(InRas, fact=3, fun=minmax, na.rm=TRUE)
OutRas_multi_label.agg <- rasterToPoints(OutRas_ML)
head(OutRas_multi_l#--> reclassify mixed pixels into mixed and unmixed (mono and multi)abel.agg)

writeRaster(OutRas_multi_label, "OutRas_multi_label.asc", overwrite=TRUE)

minmax <- function(x, na.rm=T) {
		if (min(x)==max(x))
		return(x)
		else
		return(table(x)) # value for ML cells 
    }

#########################################################################################################
#### drivers
#########################################################################################################
setwd("//crc/profiles/RedirectFolders/hichem/Desktop/Fuzzy CA/raster-USA/Amin/resolution_100m")  

library(sp)
library(raster)

dist_hwy = "dist_hwy.asc"
dist_lk_mi = "dist_lk_mi.asc"
dist_roads = "dist_roads.asc"
dist_rvrs = "dist_rvrs.asc"
dist_urb78 = "dist_urb78.asc"
dist_water = "dist_water.asc"

##################### begin 
raster_dist_hwy <- raster(dist_hwy) 
raster_dist_lk_mi <- raster(dist_lk_mi) 
raster_dist_roads <- raster(dist_roads) 
raster_dist_rvrs <- raster(dist_rvrs) 
raster_dist_urb78 <- raster(dist_urb78) 
raster_dist_water <- raster(dist_water) 

agg_raster_dist_hwy <- aggregate(raster_dist_hwy, fact=5, fun=mean, na.rm=TRUE) # 10
agg_raster_dist_lk_mi <- aggregate(raster_dist_lk_mi, fact=5, fun=mean, na.rm=TRUE)
agg_raster_dist_roads <- aggregate(raster_dist_roads, fact=5, fun=mean, na.rm=TRUE)
agg_raster_dist_rvrs <- aggregate(raster_dist_rvrs, fact=5, fun=mean, na.rm=TRUE)
agg_raster_dist_urb78 <- aggregate(raster_dist_urb78, fact=5, fun=mean, na.rm=TRUE)
agg_raster_dist_water <- aggregate(raster_dist_water, fact=5, fun=mean, na.rm=TRUE)

data.agg_raster_dist_hwy <- rasterToPoints(agg_raster_dist_hwy)
head(data.agg_raster_dist_hwy)

data.agg_raster_dist_lk_mi <- rasterToPoints(agg_raster_dist_lk_mi)
head(data.agg_raster_dist_lk_mi)

data.agg_raster_dist_roads <- rasterToPoints(agg_raster_dist_roads)
head(data.agg_raster_dist_roads)

data.agg_raster_dist_rvrs <- rasterToPoints(agg_raster_dist_rvrs)
head(data.agg_raster_dist_rvrs)
N Population Division
data.agg_raster_dist_urb78 <- rasterToPoints(agg_raster_dist_urb78)
head(data.agg_raster_dist_urb78)

data.agg_raster_dist_water <- rasterToPoints(agg_raster_dist_water)
head(data.agg_raster_dist_water)

Drivers = cbind(data.agg_raster_dist_hwy[,3], 
		    data.agg_raster_dist_lk_mi[,3],
		    data.agg_raster_dist_roads[,3],
		    data.agg_raster_dist_rvrs[,3],
		    data.agg_raster_dist_urb78[,3],
		    data.agg_raster_dist_water[,3])

colnames(Drivers) = c("dist_hwy", "dist_lk_mi", "dist_roads", "dist_rvrs", "dist_urb78", "dist_water")

writeRaster(agg_raster_dist_hwy, "agg_raster_dist_hwy.asc")
writeRaster(agg_raster_dist_lk_mi, "agg_raster_dist_lk_mi.asc")
writeRaster(agg_raster_dist_roads, "agg_raster_dist_roads.asc")
writeRaster(agg_raster_dist_rvrs, "agg_raster_dist_rvrs.asc")
writeRaster(agg_raster_dist_urb78, "agg_raster_dist_urb78.asc")
writeRaster(agg_raster_dist_water, "agg_raster_dist_water.asc")

LU_USA_ML = cbind(LU[,1:2], LU[,4], LU[,3], Drivers)
colnames(LU_USA_ML) = c("x", "y", "class", "LU_78", "dist_hwy", "dist_lk_mi", "dist_roads", "dist_rvrs", "dist_urb78", "dist_water")

class_98 = matrix(0, nrow(LU_USA_ML), 2)
class_78 = matrix(0, nrow(LU_USA_ML), 2)

for (i in 1:nrow(LU_USA_ML)){
	if(LU[i,3]==0)
		{ class_78[i,1]=1 }
	else if(LU[i,3]==1)
		{ class_78[i,2]=1 }
	else {
		class_78[i,1:2]=1 }
}

for (i in 1:nrow(LU_USA_ML)){
	if(LU[i,4]==0)
		{ class_98[i,1]=1 }
	else if(LU[i,4]==1)
		{ class_98[i,2]=1 }
	else {
		class_98[i,1:2]=1 }
}

LU_USA_multi_label = cbind(LU[,1:2], class_98, class_78, Drivers)

colnames(LU_USA_multi_label) = c("x", "y", "classNU", "classU", "Class_NU_78", "Class_U_78", "dist_hwy", "dist_lk_mi", "dist_roads", "dist_rvrs", "dist_urb78", "dist_water")


write(LU_USA_multi_label, file = "LU_USA_multi_label123.csv", sep = ",") 

############################### END OF SCRIPT ########################