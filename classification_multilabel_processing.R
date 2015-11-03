
##########################
# Reading data from landuse layers
##########################

setwd("//crc/profiles/RedirectFolders/hichem/Desktop/LISER/papers/Fuzzy CA/raster-USA/Hichem")  

library(sp)
library(raster)

LU1978 = "landusebase.asc"
LU1998 = "landusefinal.asc"

##################### begin 
#########################################################################################################
raster1978 <- raster(LU1978) 
raster1998 <- raster(LU1998) 

data_matrix1978 <- rasterToPoints(raster1978)
head(data_matrix1978)

data_matrix1998 <- rasterToPoints(raster1998)
head(data_matrix1998)

r.agg_78 <- aggregate(raster1978, fact=5, fun=sum, na.rm=TRUE) # 10

r.agg_98 <- aggregate(raster1998, fact=5, fun=sum, na.rm=TRUE)

data.agg_78 <- rasterToPoints(r.agg_78)
head(data.agg_78)

data.agg_98 <- rasterToPoints(r.agg_98)
head(data.agg_98)

r78 = data.agg_78;
r98 = data.agg_98;

for (i in 1:nrow(data.agg_78)){
if ( (data.agg_78[i,3]!=0) && (data.agg_78[i,3]!=25) ) # ML 
r78[i,3] = 2
}

for (i in 1:nrow(data.agg_98)){
if ( (data.agg_98[i,3]!=0) && (data.agg_98[i,3]!=25) ) # ML 
r98[i,3] = 2
}

for (i in 1:nrow(r78)){
if (r78[i,3]==25) # Urban 
r78[i,3] = 1
}

for (i in 1:nrow(r98)){
if (r98[i,3]==25) # Urban 
r98[i,3] = 1
}

rt1 <- rasterFromXYZ(r78)
rt2 <- rasterFromXYZ(r98)

writeRaster(rt1, "agg78_500meter.asc", overwrite=TRUE)
writeRaster(rt2, "agg98_500meter.asc", overwrite=TRUE)

LU = cbind(r78, r98[,3]) # NU(0), U(2), ML(2)
colnames(LU) = c("X", "Y", "LU78", "class") # class refer to LU98
LU[1:5,]

#### GENERATE ML land use DATA

data_matrix = cbind(data_matrix1978, data_matrix1998[,3])
agg_matrix = cbind(data.agg_78, data.agg_98[,3])

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
dt2_NU = data.agg_98[data.agg_98[,3]==0,]
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
writeRaster(InRas, "InRas.asc")

OutRas_mono_label <- aggregate(InRas, fact=3, fun=max, na.rm=TRUE)
OutRas_mono_label.agg <- rasterToPoints(OutRas_mono_label)
head(OutRas_mono_label.agg)
writeRaster(OutRas_ml, "OutRas_mono_label.asc", overwrite=TRUE)

####

### aggregation using multi-label class assignment 

OutRas_multi_label <- aggregate(InRas, fact=3, fun=minmax, na.rm=TRUE)
OutRas_multi_label.agg <- rasterToPoints(OutRas_ML)
head(OutRas_multi_label.agg)

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

