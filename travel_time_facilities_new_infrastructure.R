################################################################################################################################################################################
# Generating reductions in raw and population-weighted travel times based on new infrastructure in a given country
################################################################################################################################################################################

# Data used can be obtained from the following sources: 
# Country, Administrative unit shapefiles: downloaded from Malaria Atlas Project's 'getShp(country, ISO, extent)' function in the MalariaAtlas Project package
# Population Mask: WorldPop project https://www.worldpop.org/project/list
# Friction Surface Raster: https://map.ox.ac.uk/ accessed via getRaster(surface = "A global friction surface enumerating land-based travel speed for a nominal year 2015")
# Lakes Shapefile: http://193.43.36.146/map?entryId=bd8def30-88fd-11da-a88f-000d939bc5d8
# Viral Hemorrhagic Fever Environmental Suitability Maps: Pigott, DM et. al 2017 article https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(17)32092-5/fulltext

# -------------------------------------------------------------------------------------------------------------
# Load in packages for use in analysis
library(data.table)
library(raster)
library(rgdal)
library(ggplot2)
library(rasterVis)
library(RColorBrewer)
library(gdistance)
library(abind)
library(rje)
library(xml2)
library(rgeos)
library(malariaAtlas)
library(viridis)
library(seegSDM)
library(readr)
library(dplyr)
library(grid)
library(gridExtra)
# List of countries for analysis
  countries<-c("Angola", "Benin", "Botswana", "Burkina Faso", "Burundi", "Cameroon", "Central African Republic",
               "Chad", "Congo", "Democratic Republic of the Congo", "Djibouti", "Equatorial Guinea",
               "Eritrea", "Ethiopia", "Gabon", "Gambia", "Ghana", "Guinea", "Guinea Bissau", "Kenya",
               "Lesotho", "Liberia", "Madagascar", "Malawi", "Mali", "Mauritania", "Mozambique",
               "Namibia", "Niger", "Nigeria", "Rwanda", "Senegal","Sierra Leone", "Somalia", "South Africa",
               "South Sudan", "Sudan", "Tanzania", "Togo", "Uganda", "Zambia", "Zimbabwe")
  # List of countries that have any areas with populations <10pp / 5 x 5-km grid-cell
  worldpop.countries<-c('Somalia', "Ethiopia", "Mauritania", "Mali", "Niger", "Chad", "Sudan","Kenya", "Namibia", "Botswana", "Angola", "Djibouti", "Eritrea")
  #set a mask for whether the country is one of the worldpop countries:
  if(country %in% worldpop.countries){
    mask_population<-TRUE
  } else {
    mask_population<-FALSE
  }
  
  #rename Guinea-Bissau for use in MalariaAtlas functions:
  if(country=="Guinea Bissau"){country="Guinea-Bissau"}
  #Get outer border shapefile
  out_shp<-getShp(country<-country)
  
  #load in lakes / bodies of water in Africa and crop to current country
  lakes<-readOGR(paste0(dir,"/waterbodies_africa.shp"))
  lake_sub<-crop(lakes, out_shp)
  
  #crop out lakes from shapefile:
  loc_shp <-erase(out_shp, lake_sub)
  
  #Get Friction Raster from getRaster(surface = "A global friction surface enumerating land-based travel speed for a nominal year 2015") 
  #from MalariaAtlas package for the current country (without lakes) shapefile
  friction<-getRaster(surface = "A global friction surface enumerating land-based travel speed for a nominal year 2015", shp=loc_shp)
  #Convert Friction layer to a transition matrix so that we can identify the "least cost" pathway to a given point
  Tr <- transition(friction, function(x) 1/mean(x), 8)
  T.GC <- geoCorrection(Tr)
  
  #Read in WHO facility data and subset to the identified countries above
  WHO <- fread(paste0(dir,"/who_facility_list.csv"))
  healthfacility <-subset(WHO, WHO$Country %in% countries)
  #Convert Latitude and Longitude to Numeric values and remove any missing data
  healthfacility[, Lat := as.numeric(Lat)]
  healthfacility[, Long := as.numeric(Long)]
  healthfacility <- healthfacility[!is.na(Lat) & !is.na(Long)]
  setnames(healthfacility, c('Lat', 'Long'), c('Latitude', 'Longitude'))
  #Change Guinea-Bissau back again for the health facility data
  if(country=="Guinea-Bissau"){country="Guinea Bissau"}
  healthcare<- healthfacility[Country == country]
  
  dataset <- healthcare
  coordinates(dataset) <- ~ Longitude + Latitude
  proj4string(dataset) <- proj4string(loc_shp)
  points<-as.matrix(dataset@coords)
  
  # create initial access raster of time to most accessible health facility
  access.raster<-accCost(T.GC, points)
  
  # load environmental suitability maps
  ebov <- raster(paste0(dir,'/vhf_rasters/ebov.tif'))
  marburg <- raster(paste0(dir,'/vhf_rasters/marburg.tif'))
  lassa <- raster(paste0(dir,'/vhf_rasters/lassa.tif'))
  cchf <- raster(paste0(dir,'/vhf_rasters/cchf.tif'))
  
  # function to determine binary value for a raster based on a user-given threshold
  binary_threshold <- function(r, thresh) {
    return(calc(r, fun = function(x) ifelse(x >= thresh, 1, 0)))
  } 
  
  # crop and apply binary virus presence function to raster using median thresholds          
  ebov_presence <- binary_threshold(crop(ebov, extent(access.raster)), thresh = 0.3490)
  marburg_presence <- binary_threshold(crop(marburg, extent(access.raster)), thresh = 0.4267)
  lassa_presence <- binary_threshold(crop(lassa, extent(access.raster)), thresh = 0.3995)
  cchf_presence <- binary_threshold(crop(cchf, extent(access.raster)), thresh = 0.019)
  
  # sum binary rasters to get count of total number of VHFs present in grid-cell
  vhf_presence <- overlay(ebov_presence,
                          marburg_presence,
                          lassa_presence,
                          cchf_presence,
                          fun = function(r1, r2, r3, r4) {return(r1 + r2 + r3 + r4)})
  
  # create a binary raster for VHFs of whether there is at least one per grid-cell
  vhf_binary <- binary_threshold(vhf_presence, thresh = 1)
  
  # ------------------------------------------------------------------------------------------------------------------------------------------------------
  #For the countries with unpopulated areas, load in worldpop raster of places where population <10pp / 5 x 5 km grid-cell, set value=0 in those locations
  worldpop <- raster(paste0(dir,'/mask_master.tif'), vals=1)
  values(worldpop)[values(worldpop)==Inf]<-NA #removing infinite values
  values(worldpop)[values(worldpop)==-Inf]<-NA #removing infinite values
  # ---------------------------------------------------------------------------------
  #interpolate access.raster to 5km by 5km and crop to VHF map
  access_raster <- crop(projectRaster(access.raster, vhf_binary, method = 'bilinear'), extent(vhf_binary))
    if(mask_population==TRUE){
      worldpop_country<- crop(worldpop, extent(vhf_binary))
      #Convert the 1s to 0s, and then set areas with population to 1
      values(worldpop_country)[values(worldpop_country)==1]<-0
      values(worldpop_country)[is.na(values(worldpop_country))]<-1
      #mask out vhf absence here
      access_raster2<- gen_mask(access_raster, vhf_binary)
      access_raster2<-raster::mask(access_raster2, worldpop_country, maskvalue=0)
      vhf_access_raster_pop <<- only_country(access_raster2, access_raster2)
      return(vhf_access_raster_pop)
    }
    # Generic masking function
    gen_mask <- function(r1, r2) {
    r1[r2 == 0] <- NA
    return(r1)
   }
  access_raster2<- gen_mask(access_raster, vhf_binary)
  # mask outside of country
  only_country <- function(r1, r2) {
    r1[r2 == Inf] <- NA
    return(r1)
  }
  
  vhf_access_raster <<- only_country(access_raster2, access_raster2)
#convert the vhf_access_raster to points and into a dataframe, first removing infinite values - using full version here, but could use 
#version with unpopulated regions masked out
values(vhf_access_raster)[values(vhf_access_raster)==Inf]<-NA
potential_sites<-rasterToPoints(vhf_access_raster, spatial = FALSE )
potential_sites<-as.data.frame(potential_sites)
#this value below will tell you how many iterations of 1000 points you'll need to do
n <-nrow(potential_sites)
#generating "original" values based on contemprary assessment
original_mean<-cellStats(access_raster, stat='mean')
original_sum<-cellStats(access_raster, stat='sum')
#get original person-weighted estimates
#read-in a version of the worldpopulation data with population per 1x1 km grid and crop
worldpop <- raster(paste0(dir, '/population/worldpop_total_1y_2015_00_00.tif'))
worldpop<-crop(worldpop, extent(access_raster))
worldpop<-raster::mask(worldpop, access_raster)
ptt<-vhf_access_raster*worldpop
original_mean_ptt<-cellStats(ptt, stat='mean')
original_sum_ptt<-cellStats(ptt, stat='sum')
potential_sites<-as.data.frame(potential_sites)
#setting up parameters for iterations
a<-1000
n1<-floor(n/1000)
n1f<-a*n1
n2<-n-n1f
r1  <- rep(1:n1, each=a)
r2 <-rep(n1+1, each=n2)
r<-c(r1,r2)
#breaking up the potential sites into chunks of 1000 (with the last chunk containing the extras)
t<-split(potential_sites, r)
#get number of iterations needed
num.iters<-ceiling(n/a)
#iterate through all needed iterations, either using parallel computing (one command per iteration) or via a loop
#WARNING: this step can be really time and computational intensive
#This will save each iteration as an .RDS file specified in a specified folder (named ROutputs) in the directory chosen previously
for(it in 1:num.iters){
new_df<-t[[it]]
zonal_summary<-data.frame(new_df[,1:2], mean=NA, sum=NA, mean_diff=NA, ptt = NA, ptt_mean=NA)
# sparse sample potential_sites
iterations<-it
amount<-nrow(new_df)
datalist = list()
#make a personal travel time layer
start<-timestamp()
sample<-sample(1:nrow(new_df),amount, replace = F)
for (i in sample){
  # append lat long to dataset [the existing hospitals]
  new_dataset<-healthcare
  #create an NA row
  new_row<-new_dataset[1,]
  
  new_row[1,]<-NA
  new_row$Latitude<-new_df$x[i]
  new_row$Longitude<-new_df$y[i]
  new_dataset<-rbind(new_dataset, new_row)
  
  coordinates(new_dataset) <- ~ Latitude + Longitude
  proj4string(new_dataset) <- proj4string(loc_shp)
  new_points<-as.matrix(new_dataset@coords)
  
  new_access.raster<-gdistance::accCost(T.GC, new_points)
  
  #apply the vhf masks as above
  new_access_raster <- crop(projectRaster(new_access.raster, access_raster, method = 'bilinear'), extent(access_raster))
  
  # masking function
  gen_mask <- function(r1, r2) {
    r1[r2 == 0] <- NA
    return(r1)
  }
  # mask vhf absence
  new_access_raster2 <- gen_mask(new_access_raster, vhf_binary)
  #convert all Inf values to NA
  new_access_raster[values(new_access_raster)==Inf]<-NA
  new_access_raster2[values(new_access_raster2)==Inf]<-NA
  new_access<-new_access_raster
  new_ptt<-new_access_raster*worldpop
  
  
  zonal_summary$mean[i]<-cellStats(new_access, stat = 'mean')
  zonal_summary$sum[i]<-cellStats(new_access, stat = 'sum')
  zonal_summary$mean_diff[i]<-original_mean-zonal_summary$mean[i]
  zonal_summary$ptt_mean[i]<-cellStats(new_ptt, stat='mean')
  zonal_summary$ptt[i]<-original_mean_ptt-zonal_summary$ptt_mean[i]
  
  print(which(sample == i))
  
}

test<-na.omit(zonal_summary)
outfile<-paste0(dir,"/ROutputs/iteration_",it,"_", country,".rds")
saveRDS(test, file=outfile)
}
end<-timestamp()

#After all iterations are done, move to the next step:
setwd(paste0(dir,"/ROutputs/"))
temp<-list.files(pattern=paste0("*",country,".rds"))
test <- do.call(rbind, lapply(temp, function(x) readRDS(x)))

raw_travel <- rasterize(data.frame(test$x,test$y), access_raster, test$mean_diff, fun=mean)
raw_travel<-crop(raw_travel, extent(access_raster))
raw_int_travel<-raster::mask(raw_travel, access_raster)


weighted_travel <- rasterize(data.frame(test$x,test$y), access_raster, test$ptt, fun=mean)
weighted_travel<-crop(weighted_travel, extent(access_raster))
weighted_int_travel<-raster::mask(weighted_travel, access_raster)

#This function helps with plotting every third label to limit crowding of axes
label_fill <- function(orig, .offset=0, .mod=2, .fill=""){
  ## replace
  ii <- as.logical(
    ## offset==0 keeps first
    (1:length(orig)-1+.offset) %% .mod
  )
  orig[ii] <- .fill
  orig
}


#cleaner plotting versions, plot raw reduction
# Load in background shapefiles for neighboring countries: http://www.maplibrary.org/library/stacks/Africa/index.htm
africa<-readOGR(paste0(dir,"/Africa.shp"))
#use gBuffer with 0 km buffer to avoid "self-intersection" issues
africa <- gBuffer(africa, byid=TRUE, width=0)
#Get admin1 outlines from getShp:
admin1sub<-getShp(country<-country, admin_level=c("admin1"))
xlims<-c(extent(access_raster)[1,]-1.5, extent(access_raster)[2,]+1.5)
ylims<-c(extent(access_raster)[3,]-1.5, extent(access_raster)[4,]+1.5)
admin0sub<-crop(africa, extent(xlims, ylims))
Rev.Magma<-rasterTheme(region=magma(16,direction=-1), panel.background=list(col="light blue"))
vals<-seq(minValue(raw_int_travel),maxValue(raw_int_travel),length.out=16)
val<-round(seq(minValue(raw_int_travel),maxValue(raw_int_travel),length.out=16),2)
val_labs<-label_fill(val, .mod=3)
myColorkey<-list(space="bottom",at=vals, labels=val_labs, at=vals)

# set plot title
plot_title <- paste0('Change in Mean Travel Time to \n Most Accessible Health Facility \n in ', country)

legend_title<-"Reduction in mean travel time (Hours) \n "
plot_raster<-levelplot(raw_int_travel,par.settings=Rev.Magma,at=vals, margin=FALSE, xlim=xlims, ylim=ylims,ylab="", xlab=list(label=legend_title, cex=1.0, line=1),scales = list(draw = FALSE), colorkey=myColorkey,
                       main=list(plot_title,side=1,line=0.5,cex=1.2))
raw_reduction_plot<-plot_raster+layer_(sp.polygons(admin0sub, fill="white"))+layer(sp.polygons(lake_sub,lwd=0.00000001, fill="light blue", col="light blue"))
print("raw_reduction_plot3 done")
hp<-healthcare[,6:7]
coordinates(hp) <- ~ Longitude + Latitude
proj4string(hp) <- proj4string(loc_shp)
points<-as.matrix(hp@coords)
raw_reduction_plot_final<-raw_reduction_plot+layer(sp.polygons(admin0sub, lwd=1))+layer(sp.polygons(admin1sub, lwd=1, lty=2))+layer(sp.points(hp, cex=0.1, col="dark green", pch=20, alpha=0.5))
print("raw_reduction_plot_final done")

#plot the weighted reduction
values<-weighted_int_travel
values<-overlay(values, worldpop, fun=function(x, y){ x /mean(y, na.rm=T)})
vals<-seq(minValue(values),maxValue(values),length.out=16)
val<-round(seq(minValue(values),maxValue(values),length.out=16),2)
val_labs<-label_fill(val, .mod=3)
myColorkey<-list(space="bottom",at=vals, labels=val_labs, at=vals)

plot_title <- paste0('Change in Mean Population-Weighted Travel Time \n to Most Accessible Health Facility \n in ', country)
legend_title<-"Reduction in mean person travel time (Hours) \n "
plot_raster<-levelplot(values,par.settings=Rev.Magma, margin=FALSE, xlim=xlims, ylim=ylims,ylab="", xlab=list(label=legend_title, cex=1.0, line=1),scales = list(draw = FALSE), colorkey=myColorkey,  main=list(plot_title,side=1,line=0.5,cex=1.2))
weighted_reduction_plot<-plot_raster+layer_(sp.polygons(admin0sub, fill="white"))
hp<-healthcare[,6:7]
coordinates(hp) <- ~ Longitude + Latitude
proj4string(hp) <- proj4string(loc_shp)
points<-as.matrix(hp@coords)
weighted_reduction_plot_final<-weighted_reduction_plot+layer(sp.polygons(lake_sub,lwd=0.1, fill="light blue", col="light blue"))+layer(sp.polygons(admin0sub, lwd=1))+layer(sp.polygons(admin1sub, lwd=1, lty=2))+layer(sp.points(hp, cex=0.1, col="dark green", pch=20, alpha=0.5))
print("weighted_reduction_plot_final done")