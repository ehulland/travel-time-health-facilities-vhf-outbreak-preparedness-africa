################################################################################################################################################################################
# Calculating  Travel Times to Nearest Health Facility from Areas with VHF spillover potential
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
# function to determine binary value for a raster based on a user-given threshold
binary_threshold <- function(r, thresh) {
  return(calc(r, fun = function(x) ifelse(x >= thresh, 1, 0)))
}
# Generic masking function
gen_mask <- function(r1, r2) {
  r1[r2 == 0] <- NA
  return(r1)
}
# mask outside of country
only_country <- function(r1, r2) {
  r1[r2 == Inf] <- NA
  return(r1)
}

#Set country of interest as "country", and set directory with saved shapefiles and rasters as "dir"
#Specify "pathogen" as Ebola, Marburg, Lassa, or CCHF, or "Any" for at least one VHF per grid-cell (default is "Any")
# -------------------------------------------------------------------------------------------------------------
generate_tt<-function(country, dir, pathogen="Any"){
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
  access_raster <- crop(projectRaster(access.raster, vhf_binary, method = 'bilinear'), extent(vhf_binary))
  
  #Mask out any unpopulated areas for each VHF and the summary measure
  if(pathogen=="Ebola"){
    if(mask_population==TRUE){
      worldpop_country<- crop(worldpop, extent(ebov_presence))
      #Convert the 1s to 0s, and then set areas with population to 1
      values(worldpop_country)[values(worldpop_country)==1]<-0
      values(worldpop_country)[is.na(values(worldpop_country))]<-1
      #mask out vhf absence here
      access_raster2<- gen_mask(access_raster, ebov_presence)
      access_raster2<-raster::mask(access_raster2, worldpop_country, maskvalue=0)
      vhf_access_raster <<- only_country(access_raster2, access_raster2)
      return(vhf_access_raster)

    } else {
      access_raster2<- gen_mask(access_raster, ebov_presence)
      access_raster2<-raster::mask(access_raster2, worldpop_country, maskvalue=0)
      vhf_access_raster <<- only_country(access_raster2, access_raster2)
      return(vhf_access_raster)
    }
  }
  if(pathogen=="CCHF"){
    if(mask_population==TRUE){
      worldpop_country<- crop(worldpop, extent(cchf_presence))
      #Convert the 1s to 0s, and then set areas with population to 1
      values(worldpop_country)[values(worldpop_country)==1]<-0
      values(worldpop_country)[is.na(values(worldpop_country))]<-1
      #mask out vhf absence here
      access_raster2<- gen_mask(access_raster, cchf_presence)
      access_raster2<-raster::mask(access_raster2, worldpop_country, maskvalue=0)
      vhf_access_raster <<- only_country(access_raster2, access_raster2)
      return(vhf_access_raster)
      
    } else {
      access_raster2<- gen_mask(access_raster, cchf_presence)
      access_raster2<-raster::mask(access_raster2, worldpop_country, maskvalue=0)
      vhf_access_raster <<- only_country(access_raster2, access_raster2)
      return(vhf_access_raster)
    }
  }
  if(pathogen=="Lassa"){
    if(mask_population==TRUE){
      worldpop_country<- crop(worldpop, extent(lassa_presence))
      #Convert the 1s to 0s, and then set areas with population to 1
      values(worldpop_country)[values(worldpop_country)==1]<-0
      values(worldpop_country)[is.na(values(worldpop_country))]<-1
      #mask out vhf absence here
      access_raster2<- gen_mask(access_raster, lassa_presence)
      access_raster2<-raster::mask(access_raster2, worldpop_country, maskvalue=0)
      vhf_access_raster <<- only_country(access_raster2, access_raster2)
      return(vhf_access_raster)
      
    } else {
      access_raster2<- gen_mask(access_raster, lassa_presence)
      access_raster2<-raster::mask(access_raster2, worldpop_country, maskvalue=0)
      vhf_access_raster <<- only_country(access_raster2, access_raster2)
      return(vhf_access_raster)
    }
  }
  if(pathogen=="Marburg"){
    if(mask_population==TRUE){
      worldpop_country<- crop(worldpop, extent(marburg_presence))
      #Convert the 1s to 0s, and then set areas with population to 1
      values(worldpop_country)[values(worldpop_country)==1]<-0
      values(worldpop_country)[is.na(values(worldpop_country))]<-1
      #mask out vhf absence here
      access_raster2<- gen_mask(access_raster, marburg_presence)
      access_raster2<-raster::mask(access_raster2, worldpop_country, maskvalue=0)
      vhf_access_raster <<- only_country(access_raster2, access_raster2)
      return(vhf_access_raster)
      
    } else {
      access_raster2<- gen_mask(access_raster, marburg_presence)
      access_raster2<-raster::mask(access_raster2, worldpop_country, maskvalue=0)
      vhf_access_raster <<- only_country(access_raster2, access_raster2)
      return(vhf_access_raster)
    }
  }
  if(pathogen=="Any"){
    if(mask_population==TRUE){
      worldpop_country<- crop(worldpop, extent(vhf_binary))
      #Convert the 1s to 0s, and then set areas with population to 1
      values(worldpop_country)[values(worldpop_country)==1]<-0
      values(worldpop_country)[is.na(values(worldpop_country))]<-1
      #mask out vhf absence here
      access_raster2<- gen_mask(access_raster, vhf_binary)
      access_raster2<-raster::mask(access_raster2, worldpop_country, maskvalue=0)
      vhf_access_raster <<- only_country(access_raster2, access_raster2)
      return(vhf_access_raster)
      
    } else {
      access_raster2<- gen_mask(access_raster, vhf_binary)
      access_raster2<-raster::mask(access_raster2, worldpop_country, maskvalue=0)
      vhf_access_raster <<- only_country(access_raster2, access_raster2)
      return(vhf_access_raster)
    }
  }
}

# ---------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------
#Set country of interest as "country", and set directory with saved shapefiles and rasters as "dir"
#Specify "pathogen" as Ebola, Marburg, Lassa, or CCHF grid-cell
#Specify threshold of interest as "min" for the minimum value, "p5" for 5th percentile, "median" for the median value,
#"p95" for the 95th percentile, or "max" for maximum value. Default value is "median"
# -------------------------------------------------------------------------------------------------------------
generate_tt_thresh<-function(country, dir, pathogen, threshold="median"){
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
  
  # identfy threshold values for each pathogen:
  if(pathogen=="Ebola"){
    if(threshold=="min"){value<-0.1370}
    if(threshold=="p5"){value<-0.2100}
    if(threshold=="median"){value<-0.3490}
    if(threshold=="p95"){value<-0.4580}
    if(threshold=="max"){value<-0.5620}
    vhf_presence <- binary_threshold(crop(ebov, extent(access.raster)), thresh = value)
  }
  if(pathogen=="Marburg"){
    if(threshold=="min"){value<-0.2610}
    if(threshold=="p5"){value<-0.3000}
    if(threshold=="median"){value<-0.3990}
    if(threshold=="p95"){value<-0.6080}
    if(threshold=="max"){value<-0.6510}
    vhf_presence <- binary_threshold(crop(marburg, extent(access.raster)), thresh = value)
  }
  if(pathogen=="CCHF"){
    if(threshold=="min"){value<-0.0090}
    if(threshold=="p5"){value<-0.0140}
    if(threshold=="median"){value<-0.0190}
    if(threshold=="p95"){value<-0.0250}
    if(threshold=="max"){value<-0.0350}
    vhf_presence <- binary_threshold(crop(cchf, extent(access.raster)), thresh = value)
  }
  if(pathogen=="Lassa"){
    if(threshold=="min"){value<-0.3140}
    if(threshold=="p5"){value<-0.3390}
    if(threshold=="median"){value<-0.3955}
    if(threshold=="p95"){value<-0.4430}
    if(threshold=="max"){value<-0.5090}
    vhf_presence <- binary_threshold(crop(lassa, extent(access.raster)), thresh = value)
  }
  
  # ------------------------------------------------------------------------------------------------------------------------------------------------------
  #For the countries with unpopulated areas, load in worldpop raster of places where population <10pp / 5 x 5 km grid-cell, set value=0 in those locations
  worldpop <- raster(paste0(dir,'/worldpop_10pp.tif'), vals=1)
  values(worldpop)[values(worldpop)==Inf]<-NA #removing infinite values
  values(worldpop)[values(worldpop)==-Inf]<-NA #removing infinite values
  # ---------------------------------------------------------------------------------
  #interpolate access.raster to 5km by 5km and crop to VHF map
  access_raster <- crop(projectRaster(access.raster, vhf_presence, method = 'bilinear'), extent(vhf_presence))
  
  # mask vhf absence
  access_raster2 <- gen_mask(access_raster, vhf_presence)
  
  # mask population absence
  worldpop_country<- crop(worldpop, extent(vhf_presence))
  #Convert the 1s to 0s, and then set areas with population to 1
  values(worldpop_country)[values(worldpop_country)==1]<-0
  values(worldpop_country)[is.na(values(worldpop_country))]<-1
  #mask out here
  if(mask_population == TRUE)access_raster2<-raster::mask(access_raster2, worldpop_country, maskvalue=0)
  
  vhf_access_raster_thresh <<- only_country(access_raster2, access_raster2)
  return(vhf_access_raster_thresh)
}