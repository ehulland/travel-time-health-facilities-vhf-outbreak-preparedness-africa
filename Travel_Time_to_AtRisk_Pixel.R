################################################################################################################################################################################
# Calculating Ranked Travel Times to Nearest Location At-Risk for Viral Hemorrhagic Fevers
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

# -------------------------------------------------------------------------------------------------------------
#set country of interest as "country", and set directory with saved shapefiles and rasters as "dir"
#specify "pathogen" as Ebola, Marburg, Lassa, or CCHF, or "Any" for at least one VHF per grid-cell (default is "Any")
# -------------------------------------------------------------------------------------------------------------
# function to determine binary value for a raster based on a user-given threshold
binary_threshold <- function(r, thresh) {
  return(calc(r, fun = function(x) ifelse(x >= thresh, 1, 0)))
}
generate_at_risk_areas<-function(country, dir, pathogen="Any"){
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

#rename Guinea-Bissau for use in Malaria Atlas functions:
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
access.raster <- crop(projectRaster(access.raster, vhf_binary, method = 'bilinear'), extent(vhf_binary))

#Mask out any unpopulated areas for each VHF and the summary measure
if(pathogen=="Ebola"){
if(mask_population==TRUE){
    worldpop_country<- crop(worldpop, extent(ebov_presence))
    #Convert the 1s to 0s, and then set areas with population to 1
    values(worldpop_country)[values(worldpop_country)==1]<-0
    values(worldpop_country)[is.na(values(worldpop_country))]<-1
    #mask out here
    ebov_presenceT<-raster::mask(ebov_presence, worldpop_country, maskvalue=0)
    #convert locations with Ebola presence (e.g. any value >0) to points
    ebolarisk<- rasterToPoints(ebov_presenceT, 
                               fun = function(ebov_presenceT){
                                 ebov_presenceT > 0
                               }, spatial=TRUE)
#Retain only areas inside of the country
risk<<-gIntersection(ebolarisk, loc_shp)
return(risk)
  } else {
ebolarisk<- rasterToPoints(ebov_presence, 
                               fun = function(ebov_presence){
                                 ebov_presence > 0
                               }, spatial=TRUE)
#Retain only areas inside of the country
risk<<-gIntersection(ebolarisk, loc_shp)
return(risk)
print("Ebola Complete")
  }
}
if(pathogen=="CCHF"){
if(mask_population==TRUE){
  worldpop_country<- crop(worldpop, extent(cchf_presence))
  #Convert the 1s to 0s, and then set areas with population to 1
  values(worldpop_country)[values(worldpop_country)==1]<-0
  values(worldpop_country)[is.na(values(worldpop_country))]<-1
  #mask out here
  cchf_presenceT<-raster::mask(cchf_presence, worldpop_country, maskvalue=0)
  #convert locations with CCHF presence (e.g. any value >0) to points
  cchfrisk<- rasterToPoints(cchf_presenceT, 
                             fun = function(cchf_presenceT){
                               cchf_presenceT > 0
                             }, spatial=TRUE)
  #Retain only areas inside of the country
  risk<<-gIntersection(cchfrisk, loc_shp)
  return(risk)
}
cchfrisk<- rasterToPoints(cchf_presence, 
                             fun = function(cchf_presence){
                               cchf_presence > 0
                             }, spatial=TRUE)
#Retain only areas inside of the country
risk<<-gIntersection(cchfrisk, loc_shp)
return(risk)
print("CCHF complete")
}
if(pathogen=="Lassa"){
if(mask_population==TRUE){
  worldpop_country<- crop(worldpop, extent(lassa_presence))
  #Convert the 1s to 0s, and then set areas with population to 1
  values(worldpop_country)[values(worldpop_country)==1]<-0
  values(worldpop_country)[is.na(values(worldpop_country))]<-1
  #mask out here
  lassa_presenceT<-raster::mask(lassa_presence, worldpop_country, maskvalue=0)
  #convert locations with Lassa presence (e.g. any value >0) to points
  lassarisk<- rasterToPoints(lassa_presenceT, 
                             fun = function(lassa_presenceT){
                               lassa_presenceT > 0
                             }, spatial=TRUE)
#Retain only areas inside of the country
risk<<-gIntersection(lassarisk, loc_shp)
return(risk)
}
lassarisk<- rasterToPoints(lassa_presence, 
                             fun = function(lassa_presence){
                               lassa_presence > 0
                             }, spatial=TRUE)
#Retain only areas inside of the country
risk<<-gIntersection(lassarisk, loc_shp)
return(risk)
print("Lassa complete")
}
if(pathogen=="Marburg"){
if(mask_population==TRUE){
  worldpop_country<- crop(worldpop, extent(marburg_presence))
  #Convert the 1s to 0s, and then set areas with population to 1
  values(worldpop_country)[values(worldpop_country)==1]<-0
  values(worldpop_country)[is.na(values(worldpop_country))]<-1
  #mask out here
  marburg_presenceT<-raster::mask(marburg_presence, worldpop_country, maskvalue=0)
  #convert locations with Marburg presence (e.g. any value >0) to points
  marburgrisk<- rasterToPoints(marburg_presenceT, 
                             fun = function(marburg_presenceT){
                               marburg_presenceT > 0
                             }, spatial=TRUE)
#Retain only areas inside of the country
risk<<-gIntersection(marburgrisk, loc_shp)
return(risk)
}
marburgrisk<- rasterToPoints(marburg_presence, 
                             fun = function(marburg_presence){
                               marburg_presence> 0
                             }, spatial=TRUE)
#Retain only areas inside of the country
risk<<-gIntersection(marburgrisk, loc_shp)
return(risk)
print("Marburg Complete")
}
if(pathogen=="Any"){
if(mask_population==TRUE){
  worldpop_country<- crop(worldpop, extent(vhf_binary))
  #Convert the 1s to 0s, and then set areas with population to 1
  values(worldpop_country)[values(worldpop_country)==1]<-0
  values(worldpop_country)[is.na(values(worldpop_country))]<-1
  #mask out here
  vhf_binaryT<-raster::mask(vhf_binary, worldpop_country, maskvalue=0)
  #convert locations with VHF presence (e.g. any value >0) to points
  VHFrisk<- rasterToPoints(vhf_binaryT, 
                             fun = function(vhf_binaryT){
                               vhf_binaryT > 0
                             }, spatial=TRUE)
  #Retain only areas inside of the country
  risk<<-gIntersection(VHFrisk, loc_shp)
  return(risk)
}
VHFrisk<- rasterToPoints(vhf_binary, 
                             fun = function(vhf_binary){
                               vhf_binary > 0
                             }, spatial=TRUE)
#Retain only areas inside of the country
risk<<-gIntersection(VHFrisk, loc_shp)
return(risk)
print("Any VHF complete")
}
}

#Make sure these commands are the same as above!
plot_at_risk_map_all<-function(country,dir, pathogen="Any"){
# Generate Map (or two maps for those worldpop.countries with masked out)
# Load in background shapefiles for neighboring countries: http://www.maplibrary.org/library/stacks/Africa/index.htm
africa<-readOGR(paste0(dir,"/Africa.shp"))
#use gBuffer with 0 km buffer to avoid "self-intersection" issues
africa <- gBuffer(africa, byid=TRUE, width=0)
#Get admin1 outlines from getShp:
admin1sub<-getShp(country<-country, admin_level=c("admin1"))

legend_title <- 'Travel time (Hours) in-country \n'
# set plot title
if(pathogen=="Any"){
plot_title <- paste0('Travel time to Most Accessible \n At-Risk Grid-Cell for At Least \n One VHF in ', country)
} else{
  plot_title <- paste0('Travel time to Most Accessible \n At-Risk Grid-Cell for ', pathogen,  ' in ', country)
}
# List of countries for analysis
countries<-c("Angola", "Benin", "Botswana", "Burkina Faso", "Burundi", "Cameroon", "Central African Republic",
             "Chad", "Congo", "Democratic Republic of the Congo", "Djibouti", "Equatorial Guinea",
             "Eritrea", "Ethiopia", "Gabon", "Gambia", "Ghana", "Guinea", "Guinea Bissau", "Kenya",
             "Lesotho", "Liberia", "Madagascar", "Malawi", "Mali", "Mauritania", "Mozambique",
             "Namibia", "Niger", "Nigeria", "Rwanda", "Senegal","Sierra Leone", "Somalia", "South Africa",
             "South Sudan", "Sudan", "Tanzania", "Togo", "Uganda", "Zambia", "Zimbabwe")

#rename Guinea-Bissau for use in Malaria Atlas functions:
if(country=="Guinea Bissau"){country="Guinea-Bissau"}
#Get outer border shapefile
out_shp<-getShp(country<-country)

#load in lakes / bodies of water in Africa and crop to current country
countries.w.lakes<- c("Angola", "Benin", "Botswana", "Burkina Faso", "Burundi", "Cameroon","Chad", "Congo", 
                      "Cote d'Ivoire","Djibouti", "Democratic Republic of the Congo", "Ethiopia", "Gabon", 
                      "Ghana",  "Kenya", "Malawi", "Mali", "Mauritania", "Mozambique", "Namibia", "Niger",
                      "Nigeria", "Rwanda", "Senegal",  "Sudan", "Tanzania","Uganda", "Zambia", "Zimbabwe", "Swaziland")
if(country %in% countries.w.lakes){
  lakes<-readOGR(paste0(dir,"/waterbodies_africa.shp"))
  lake_sub<-crop(lakes, out_shp)
  
  #crop out lakes from shapefile:
  loc_shp <-erase(out_shp, lake_sub)
} else {
  lake_sub<-NA
  loc_shp<-out_shp
}

#Get Friction Raster from getRaster(surface = "A global friction surface enumerating land-based travel speed for a nominal year 2015") 
#from MalariaAtlas package for the current country (without lakes) shapefile
friction<-getRaster(surface = "A global friction surface enumerating land-based travel speed for a nominal year 2015", shp=loc_shp)
#Convert Friction layer to a transition matrix so that we can identify the "least cost" pathway to a given point
Tr <- transition(friction, function(x) 1/mean(x), 8)
T.GC <- geoCorrection(Tr)
VHFrisk_access<-accCost(T.GC, risk)
#Load in one 5x5 km VHF raster to resample to 5x5 km
ebov <- raster(paste0(dir,'/vhf_rasters/ebov.tif'))
VHFrisk.access<-crop(projectRaster(VHFrisk_access, ebov, method = 'bilinear'), extent(loc_shp))

#Remove infinte values and convert to hours
VHFrisk.access[VHFrisk.access==Inf]<-NA
values(VHFrisk.access)<-values(VHFrisk.access)/60
#Restrict upper limit to 24+ hours to standardize scale
values(VHFrisk.access)<-ifelse(values(VHFrisk.access)>24 & values(VHFrisk.access)<Inf, 24, values(VHFrisk.access))
lims=seq(0,24, length.out=16) 

#Set plot limits
xlims<-c(extent(VHFrisk.access)[1,]-1.5, extent(VHFrisk.access)[2,]+1.5)
ylims<-c(extent(VHFrisk.access)[3,]-1.5, extent(VHFrisk.access)[4,]+1.5)
admin0sub<-crop(africa, extent(xlims, ylims))
#set color scheme
Rev.Magma<-rasterTheme(region=magma(16,direction=-1), panel.background=list(col="light blue"))
np<-levelplot(VHFrisk.access,par.settings=Rev.Magma,at=lims, margin=FALSE, xlim=xlims, ylim=ylims,ylab="", xlab=list(label=legend_title, cex=1.0, line=1),scales = list(draw = FALSE), colorkey=list(space='bottom'),
                         main=list(plot_title,side=1,line=0.5,cex=1.2))
np1<-np+layer_(sp.polygons(loc_shp, lwd=0.5, col="white",fill="white"))+layer(sp.polygons(lake_sub,lwd=0.00000001, fill="light blue", col="light blue"))
np2<-np1+layer(sp.polygons(admin1sub, lwd=1, lty=2))
np3<-np2+layer_(sp.polygons(admin0sub, lwd=1, fill="white"))
#create raster data for the "at-risk" points with value of 0 for each location
pts<-rasterize(risk, VHFrisk.access, 0)
pts2<-crop(pts, loc_shp)
pts2<-raster::mask(pts2, loc_shp)
#set color
gn<-rasterTheme(region="#99CC99")
Plot_NoPop<<-np3+levelplot(pts2, par.settings=gn)+layer(sp.polygons(admin1sub, lwd=1, lty=2))+layer(sp.polygons(admin0sub, lwd=0.5))+layer(sp.polygons(loc_shp, lwd=2))
return(Plot_NoPop)
}

#only use this function for those "world pop countries"
#worldpop.countries<-c('Somalia', "Ethiopia", "Mauritania", "Mali", "Niger", "Chad", "Sudan","Kenya", "Namibia", "Botswana", "Angola", "Djibouti", "Eritrea")
plot_at_risk_map_pop<-function(country,dir, pathogen="Any"){
  if(pathogen=="Any"){
    plot_title <- paste0('Travel time to Most Accessible \n At-Risk Grid-Cell for At Least \n One VHF in ', country, 'from \n Populated Locations')
  } else {
    plot_title <- paste0('Travel time to Most Accessible \n At-Risk Grid-Cell for ', pathogen,  '\n in ', country, 'from Populated Locations')
  }
  
  # List of countries for analysis
  countries<-c("Angola", "Benin", "Botswana", "Burkina Faso", "Burundi", "Cameroon", "Central African Republic",
               "Chad", "Congo", "Democratic Republic of the Congo", "Djibouti", "Equatorial Guinea",
               "Eritrea", "Ethiopia", "Gabon", "Gambia", "Ghana", "Guinea", "Guinea Bissau", "Kenya",
               "Lesotho", "Liberia", "Madagascar", "Malawi", "Mali", "Mauritania", "Mozambique",
               "Namibia", "Niger", "Nigeria", "Rwanda", "Senegal","Sierra Leone", "Somalia", "South Africa",
               "South Sudan", "Sudan", "Tanzania", "Togo", "Uganda", "Zambia", "Zimbabwe")
  #rename Guinea-Bissau for use in Malaria Atlas functions:
  if(country=="Guinea Bissau"){country="Guinea-Bissau"}
  #Get outer border shapefile
  out_shp<-getShp(country<-country)
  
  #load in lakes / bodies of water in Africa and crop to current country
  countries.w.lakes<- c("Angola", "Benin", "Botswana", "Burkina Faso", "Burundi", "Cameroon","Chad", "Congo", 
                        "Cote d'Ivoire","Djibouti", "Democratic Republic of the Congo", "Ethiopia", "Gabon", 
                        "Ghana",  "Kenya", "Malawi", "Mali", "Mauritania", "Mozambique", "Namibia", "Niger",
                        "Nigeria", "Rwanda", "Senegal",  "Sudan", "Tanzania","Uganda", "Zambia", "Zimbabwe", "Swaziland")
  if(country %in% countries.w.lakes){
    lakes<-readOGR(paste0(dir,"/waterbodies_africa.shp"))
    lake_sub<-crop(lakes, out_shp)
    
    #crop out lakes from shapefile:
    loc_shp <-erase(out_shp, lake_sub)
  } else {
    lake_sub<-NA
    loc_shp<-out_shp
  }
  
  #Get Friction Raster from getRaster(surface = "A global friction surface enumerating land-based travel speed for a nominal year 2015") 
  #from MalariaAtlas package for the current country (without lakes) shapefile
  friction<-getRaster(surface = "A global friction surface enumerating land-based travel speed for a nominal year 2015", shp=loc_shp)
  #Convert Friction layer to a transition matrix so that we can identify the "least cost" pathway to a given point
  Tr <- transition(friction, function(x) 1/mean(x), 8)
  T.GC <- geoCorrection(Tr)
  VHFrisk_accessP<-accCost(T.GC, risk)
  #Load in one 5x5 km VHF raster to resample to 5x5 km
  ebov <- raster(paste0(dir,'/vhf_rasters/ebov.tif'))
  VHFrisk.accessP<-crop(projectRaster(VHFrisk_accessP, ebov, method = 'bilinear'), extent(loc_shp))
 
  #Remove infinte values and convert to hours
  VHFrisk.accessP[VHFrisk.accessP==Inf]<-NA
  values(VHFrisk.accessP)<-values(VHFrisk.accessP)/60
  #Restrict upper limit to 24+ hours to standardize scale
  values(VHFrisk.accessP)<-ifelse(values(VHFrisk.accessP)>24 & values(VHFrisk.accessP)<Inf, 24, values(VHFrisk.accessP))
  lims=seq(0,24, length.out=16) 
  
  #Set plot limits
  xlims<-c(extent(VHFrisk.accessP)[1,]-1.5, extent(VHFrisk.accessP)[2,]+1.5)
  ylims<-c(extent(VHFrisk.accessP)[3,]-1.5, extent(VHFrisk.accessP)[4,]+1.5)
  admin0sub<-crop(africa, extent(xlims, ylims))
  #set color scheme
  Rev.Magma<-rasterTheme(region=magma(16,direction=-1), panel.background=list(col="light blue"))
  p<-levelplot(VHFrisk.accessP,par.settings=Rev.Magma,at=lims, margin=FALSE, xlim=xlims, ylim=ylims,ylab="", xlab=list(label=legend_title, cex=1.0, line=1),scales = list(draw = FALSE), colorkey=list(space='bottom'),
               main=list(plot_title,side=1,line=0.5,cex=1.2))
  p1<-p+layer(sp.polygons(lake_sub,lwd=0.00000001, fill="light blue", col="light blue"))
  p2<-p1+layer_(sp.polygons(admin0sub, lwd=1, fill="white"))
  #create raster data for the "at-risk" points with value of 0 for each location
  pts<-rasterize(risk, VHFrisk.access, 0)
  pts2<-crop(pts, loc_shp)
  pts2<-raster::mask(pts2, loc_shp)
  #set color
  gn<-rasterTheme(region="#99CC99")
  p3<-p2+levelplot(pts2, par.settings=gn)
  worldpop <- raster(paste0(dir,'/mask_master.tif'), vals=1)
  worldpop_country<- crop(worldpop, extent(loc_shp))
  #Convert the 1s to 0s, and then set areas with population to 1
  values(worldpop_country)[values(worldpop_country)==1]<-0
  values(worldpop_country)[is.na(values(worldpop_country))]<-1
  worldPop<-worldpop_country
  worldPop[worldPop>0]<-NA
  worldPop<-raster::mask(worldPop, loc_shp)
  Plot_Pop<<-p3+levelplot(worldPop,par.settings=rasterTheme(region="light grey"))+layer(sp.polygons(admin1sub, lwd=1, lty=2))+layer(sp.polygons(admin0sub, lwd=0.5))+layer(sp.polygons(loc_shp, lwd=2))
  return(Plot_Pop)
}

