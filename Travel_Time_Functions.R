################################################################################################################################################################################
# Code to generate all functions
################################################################################################################################################################################
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
library(stringr)
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

#Generating raster of travel times to most accessible facility from locations with pathogen spillover potential
#Returns raster "pathogen_access_raster" with travel times to health facility at 5x5-km resolution, in minutes
# -------------------------------------------------------------------------------------------------------------
#Set country of interest as "country", and set directory with saved shapefiles and rasters as "dir"
#Specify a raster file (5 x 5km resolution) with binary classification of presence of one or more pathogens
#Can use binary_threshold function to get binary map from continuous map
# -------------------------------------------------------------------------------------------------------------
generate_tt<-function(country, dir, pathogen_map){
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
healthfacility<- fread(paste0(dir,"/who_facility_list.csv"))
#Convert Latitude and Longitude to Numeric values and remove any missing data
healthfacility[, Lat := as.numeric(Lat)]
healthfacility[, Long := as.numeric(Long)]
healthfacility <- healthfacility[!is.na(Lat) & !is.na(Long)]
setnames(healthfacility, c('Lat', 'Long'), c('Latitude', 'Longitude'))
#Subset to specified country
healthcare<- healthfacility[Country == country]
#Convert the health facilities to matrix of points for use with least cost function
dataset <- healthcare
coordinates(dataset) <- ~ Longitude + Latitude
proj4string(dataset) <- proj4string(loc_shp)
points<-as.matrix(dataset@coords)
# create initial access raster of time to most accessible health facility using least cost algorithm
access.raster<-accCost(T.GC, points)
#interpolate access.raster to 5km by 5km and crop to VHF map
access_raster <- crop(projectRaster(access.raster, pathogen_map, method = 'bilinear'), extent(pathogen_map))
#Mask out locations without VHF presence using gen_mask function
access_raster2<- gen_mask(access_raster, ebov_presence)
#Mask out any areas outside of the country
pathogen_access_raster <<- only_country(access_raster2, access_raster2)
return(pathogen_access_raster)
}

#Generating raster of travel times to most accessible at-risk grid-cell from any location
#Returns raster "pathogen.risk.access" with travel times to at-risk areas at 5x5-km resolution, in minutes
#Those locations with travel times of 0 are locations with pathogen spillover potential
# -------------------------------------------------------------------------------------------------------------
#Set country of interest as "country", and set directory with saved shapefiles and rasters as "dir"
#Specify a raster file (5 x 5km resolution) with binary classification of presence of one or more pathogens
#Can use binary_threshold function to get binary map from continuous map
generate_at_risk_areas<-function(country, dir, pathogen_map){
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
healthfacility<- fread(paste0(dir,"/who_facility_list.csv"))
#Convert Latitude and Longitude to Numeric values and remove any missing data
healthfacility[, Lat := as.numeric(Lat)]
healthfacility[, Long := as.numeric(Long)]
healthfacility <- healthfacility[!is.na(Lat) & !is.na(Long)]
setnames(healthfacility, c('Lat', 'Long'), c('Latitude', 'Longitude'))
#Subset to specified country
healthcare<- healthfacility[Country == country]
#Convert the health facilities to matrix of points for use with least cost function
dataset <- healthcare
coordinates(dataset) <- ~ Longitude + Latitude
proj4string(dataset) <- proj4string(loc_shp)
points<-as.matrix(dataset@coords)
# create initial access raster of time to most accessible health facility using least cost algorithm
access.raster<-accCost(T.GC, points)
#interpolate access.raster to 5km by 5km and crop to VHF map
access_raster <- crop(projectRaster(access.raster, pathogen_map, method = 'bilinear'), extent(pathogen_map))

#Convert locations with pathogen presence (e.g. any value >0) to points
pathogenrisk<- rasterToPoints(pathogen_map, 
                                 fun = function(pathogen_map){
                                   pathogen_map > 0
                                 }, spatial=TRUE)
#Retain only areas inside of the country
risk<-gIntersection(pathogenrisk, loc_shp)
#Get travel time to the locations at risk from any other location in country (from transition layer)
pathogen_risk_access<-accCost(T.GC, risk)
#Resample to 5x5 km
pathogen.risk.access<<-crop(projectRaster(pathogen_risk_access, pathogen_map, method = 'bilinear'), extent(loc_shp))
}

#Generating table of travel times to most accessible hospital from at-risk location
#Returns raster "hosp.table" with travel times to hospitals, admin unit, country and a hospital name / admin unit
#combination column for plotting. Returns a table of length num.rows, as specified in the function.
# -------------------------------------------------------------------------------------------------------------
#Set country of interest as "country", and set directory with saved shapefiles and rasters as "dir"
#Specify a raster file (5 x 5km resolution) of travel times to at-risk areas (pathogen.risk.access) from generate_at_risk_areas
#Specify number of hospitals and travel times to be returned with "num.rows"
generate_hospital_table<-function(country, dir, pathogen.risk.access, num.rows){
  #Get outer border shapefile
  out_shp<-getShp(country<-country)
  
  #Read in WHO facility data and subset to the identified countries above
  healthfacility<- fread(paste0(dir,"/who_facility_list.csv"))
  # Clean up the facility type using str_extract
  healthfacility$type<-stringr::str_extract(string = healthfacility$`Facility type`, 
                                            pattern = "[Hh]ospital|[Hh]ealth [Cc]entre|[Hh]ealth [Cc]enter|[Hh]ealth [Pp]ost|[Hh]ealth [Cc]linic|[Pp]oste [Dd]e [Ss]ant\xe9|[Hh]opital|[Dd]ispens|[Cc]entro de [Ss]a\xfade|[Hh]\xf4pital|[Cc]entre de [Ss]ant|[Pp]olyclin|[Cc]entre [Hh]\xf4pital|[Cc]entre [Mm]\xe9dic|Centre Medical|[Mm]edical [Cc]entre|[Cc]lini|DISPENSARY|H\xf4pital|Health Hut|Health Station|Posto de Sa\xfade|Community|Primary|Centro de Sant\xe9|Materno|Centre Medico|Centro de Sant|Referral|Unit\xe9 de Soins P\xe9riph\xe9rique|Unites de Sant\xe9 de [Vv]illage")
  healthfacility$FacilityType<-stringi::stri_enc_toutf8(healthfacility$type, is_unknown_8bit =TRUE, validate =TRUE)
  healthfacility$FacilityType<-trimws(healthfacility$FacilityType, which ="both")
  healthfacility$FacilityType[healthfacility$FacilityType %in% c("hospital", "Hospital", "H�pital", "Centre H�pital" ,"Referral")]<-"Hospital"
  healthfacility$FacilityType[healthfacility$FacilityType %in% c("DISPENSARY", "Dispens")]<-"Dispensary"
  healthfacility$FacilityType[healthfacility$FacilityType %in% c("Centre Medical", "Centre Medico", "Centre M�dic","Medical Centre")]<-"Medical Center"
  healthfacility$FacilityType[healthfacility$FacilityType %in% c("Health post","Health Post","Poste de sant�", "Poste de Sant�", "Poste De Sant�", "Posto de Sa�de")]<-"Health Post"
  healthfacility$FacilityType[healthfacility$FacilityType %in% c("Clini","Health Clinic", "Primary")]<-"Health Clinic"
  healthfacility$FacilityType[healthfacility$FacilityType %in% c("Centre de Sant","Health Centre","Centro de Sant�", "Centro de Sa�de")]<-"Health Centre"
  healthfacility$FacilityType[healthfacility$FacilityType %in% c("Unites de Sant� de village","Unites de Sant� de Village", "Unit� de Soins P�riph�rique","Community")]<-"Community Health Unit"
  healthfacility$FacilityType[healthfacility$FacilityType == "Materno"]<-"Maternity"
  healthfacility$FacilityType[healthfacility$FacilityType == "Polyclin"]<-"Polyclinic"
  #Convert Latitude and Longitude to Numeric values and remove any missing data
  healthfacility[, Lat := as.numeric(Lat)]
  healthfacility[, Long := as.numeric(Long)]
  healthfacility <- healthfacility[!is.na(Lat) & !is.na(Long)]
  setnames(healthfacility, c('Lat', 'Long'), c('Latitude', 'Longitude'))
  #Subset to specified country
  healthcare<- healthfacility[Country == country]
  
  #restrict to Hospitals only
  hospitals<-healthcare[FacilityType=="Hospital"]
  coordinates(hospitals) <- ~ Longitude + Latitude
  proj4string(hospitals) <- proj4string(out_shp)
  points<-as.matrix(hospitals@coords)
  
  #get the raster values for travel time at each of the locations of the hospitals and bind back to the locations / names; remove any missing data
  TT <- raster::extract(pathogen.risk.access, points)
  hospitals <- cbind(hospitals, TT)
  names(hospitals)[ncol(hospitals)] <- 'Travel_Time'
  hospitals <- hospitals[!is.na(hospitals$Travel_Time),]
  
  # Order hospitals with lowest n travel times by country (dictated by num.rows)
  hosp_tt <- hospitals[order(hospitals$Travel_Time),]
  hosp_tt <- hosp_tt[1:num.rows,]
  
  #clean names for printing 
  hosp_tt$name_index<-gsub("\xf4",  "o", hosp_tt$Facility.name)
  hosp_tt$name_index<-gsub("\xe9",  "e", hosp_tt$name_index)
  hosp_tt$name_index<-gsub("\U3e32393cH",  "'h", hosp_tt$name_index)
  hosp_tt$name_index<-gsub("\U3e32393c",  "'", hosp_tt$name_index)
  hosp_tt$name_index<-gsub("\U3e66653c",  "i", hosp_tt$name_index)
  hosp_tt$name_index<-gsub("\x92",  "'", hosp_tt$name_index)
  hosp_tt$name_index<-gsub("\U3e62653",  "e", hosp_tt$name_index)
  hosp_tt$name_index<-gsub("\xf3",  "o", hosp_tt$name_index)
  hosp_tt$name_index<-gsub("\xe3",  "a", hosp_tt$name_index)
  hosp_tt$name_index<-gsub("\xe1",  "a", hosp_tt$name_index)
  hosp_tt$name_index<-gsub("\xed",  "i", hosp_tt$name_index)
  hosp_tt$name_index<-gsub("\xe7",  "c", hosp_tt$name_index)
  hosp_tt$admin1<-gsub("\xf4",  "o", hosp_tt$Admin1)
  hosp_tt$admin1<-gsub("\xe9",  "e", hosp_tt$admin1)
  hosp_tt$admin1<-gsub("\U3e32393cH",  "'h", hosp_tt$admin1)
  hosp_tt$admin1<-gsub("\U3e32393c",  "'", hosp_tt$admin1)
  hosp_tt$admin1<-gsub("\U3e66653c",  "i", hosp_tt$admin1)
  hosp_tt$admin1<-gsub("\x92",  "'", hosp_tt$admin1)
  hosp_tt$admin1<-gsub("\U3e62653",  "e", hosp_tt$admin1)
  hosp_tt$admin1<-gsub("\xf3",  "o", hosp_tt$admin1)
  hosp_tt$admin1<-gsub("\xe3",  "a", hosp_tt$admin1)
  hosp_tt$admin1<-gsub("\xe1",  "a", hosp_tt$admin1)
  hosp_tt$admin1<-gsub("\xed",  "i", hosp_tt$admin1)
  hosp_tt$admin1<-gsub("\xe7",  "c", hosp_tt$admin1)
  hosp<-hosp_tt[,c("name_index", "admin1", "Country", "Travel_Time")]
  
  #Bind name of hospital and the 1st administrative district for printing
  hosp$name_admin<-paste0(hosp$name_index," (", hosp$admin1, ") ")
  #round to two decimal places for printing
  hosp$Travel_Time<-round(hosp$Travel_Time,2)
  hosp$name_admin<-reorder(hosp$name_admin, -hosp$Travel_Time)
  hosp.table<<-hosp
  return(hosp.table)
}


#This function converts the travel times to nearest health facility raster into a dataset of points - one per 5x5-km grid-cell
#and then iterates over each point by placing new infrastructure in that point and recalculating travel times. 
#Each country has it's points divided into segments of 1000, and this function iterates over each point in each segment
#WARNING: This step can be really time and computational intensive - recommend running in parallel
#Returns a series of .RDS files, saved in a subfolder (ROutputs) of "dir", of each set of 1000 points as a dataframe with long (x), lat (y), 
#original mean (mean) of travel times, original mean of population-weighted travel times (ptt) as well as the
#reduction in travel times for that location - raw (mean_diff) and population-weighted (ptt_mean)
# -------------------------------------------------------------------------------------------------------------
#Set country of interest as "country", and set directory with saved shapefiles and rasters as "dir"
#Specify a raster file (5 x 5km resolution) with binary classification of presence of one or more pathogens
#Can use binary_threshold function to get binary map from continuous map
#Requires a raster of world populations("worldpop") by 1x1-km grid-cells, cropped to country of interest
# -------------------------------------------------------------------------------------------------------------
generate_new_infrastructure_splits<-function(country, dir, pathogen_map, worldpop){
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
  healthfacility<- fread(paste0(dir,"/who_facility_list.csv"))
  #Convert Latitude and Longitude to Numeric values and remove any missing data
  healthfacility[, Lat := as.numeric(Lat)]
  healthfacility[, Long := as.numeric(Long)]
  healthfacility <- healthfacility[!is.na(Lat) & !is.na(Long)]
  setnames(healthfacility, c('Lat', 'Long'), c('Latitude', 'Longitude'))
  #Subset to specified country
  healthcare<- healthfacility[Country == country]
  #Convert the health facilities to matrix of points for use with least cost function
  dataset <- healthcare
  coordinates(dataset) <- ~ Longitude + Latitude
  proj4string(dataset) <- proj4string(loc_shp)
  points<-as.matrix(dataset@coords)
  # create initial access raster of time to most accessible health facility using least cost algorithm
  access.raster<-accCost(T.GC, points)
  #interpolate access.raster to 5km by 5km and crop to VHF map
  access_raster <- crop(projectRaster(access.raster, pathogen_map, method = 'bilinear'), extent(pathogen_map))
  #Mask out any areas outside of the country
  access.raster <<- only_country(access_raster, access_raster)
  #convert the access.raster to points and into a dataframe, first removing infinite values
  values(access.raster)[values(access.raster)==Inf]<-NA
  potential_sites<-rasterToPoints(access.raster, spatial = FALSE )
  potential_sites<-as.data.frame(potential_sites)
  #this value below will tell you how many iterations of 1000 points you'll need to do
  n <-nrow(potential_sites)
  #generating "original" values based on contemprary assessment
  original_mean<-cellStats(access.raster, stat='mean')
  #get original person-weighted estimates
  ptt<-access.raster*worldpop
  original_mean_ptt<-cellStats(ptt, stat='mean')
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
    #This will save each iteration as an .RDS file specified in a specified folder (named ROutputs) in the directory chosen previously
  for(it in 1:num.iters){
    new_df<-t[[it]]
    zonal_summary<-data.frame(new_df[,1:2], mean=NA,  mean_diff=NA, ptt = NA, ptt_mean=NA)
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
      new_access_raster <- crop(projectRaster(new_access.raster, access.raster, method = 'bilinear'), extent(access.raster))
      #convert all Inf values to NA
      new_access_raster[values(new_access_raster)==Inf]<-NA
      new_access<-new_access_raster
      new_ptt<-new_access_raster*worldpop
      
      
      zonal_summary$mean[i]<-cellStats(new_access, stat = 'mean')
      zonal_summary$mean_diff[i]<-original_mean-zonal_summary$mean[i]
      zonal_summary$ptt_mean[i]<-cellStats(new_ptt, stat='mean')
      zonal_summary$ptt[i]<-original_mean_ptt-zonal_summary$ptt_mean[i]
      print(which(sample == i))
      
    }
    out_dir<-paste0(dir, "/ROutputs/")
    if(!dir.exists(out_dir)){
      dir.create(out_dir)
    }
    test<-na.omit(zonal_summary)
    outfile<-paste0(out_dir,"iteration_",it,"_", country,".rds")
    saveRDS(test, file=outfile)
  }
  end<-timestamp()
}

#This function reads in all the .RDS files produced by "generate_new_infrastructure_split" in a given 
#subdirectory (ROutputs) in "dir" and creates the country raster of raw reduction in travel times, 
#in minutes ("raw_int_travel")
# -------------------------------------------------------------------------------------------------------------
#Set country of interest as "country", and set directory with saved .RDS files as "dir"
#Specify a raster file (5 x 5km resolution) with binary classification of presence of one or more pathogens
#Can use binary_threshold function to get binary map from continuous map
# -------------------------------------------------------------------------------------------------------------
generate_new_infrastructure_raw<-function(country, dir, pathogen_map){
  in_dir<-paste0(dir, "/ROutputs/")
  setwd(in_dir)
  temp<-list.files(pattern=paste0("*",country,".rds"))
  test <- do.call(rbind, lapply(temp, function(x) readRDS(x)))

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
  healthfacility<- fread(paste0(dir,"/who_facility_list.csv"))
  #Convert Latitude and Longitude to Numeric values and remove any missing data
  healthfacility[, Lat := as.numeric(Lat)]
  healthfacility[, Long := as.numeric(Long)]
  healthfacility <- healthfacility[!is.na(Lat) & !is.na(Long)]
  setnames(healthfacility, c('Lat', 'Long'), c('Latitude', 'Longitude'))
  #Subset to specified country
  healthcare<- healthfacility[Country == country]
  #Convert the health facilities to matrix of points for use with least cost function
  dataset <- healthcare
  coordinates(dataset) <- ~ Longitude + Latitude
  proj4string(dataset) <- proj4string(loc_shp)
  points<-as.matrix(dataset@coords)
  # create initial access raster of time to most accessible health facility using least cost algorithm
  access.raster<-accCost(T.GC, points)
  #interpolate access.raster to 5km by 5km and crop to VHF map
  access_raster <- crop(projectRaster(access.raster, pathogen_map, method = 'bilinear'), extent(pathogen_map))
  #Mask out any areas outside of the country
  access.raster <- only_country(access_raster, access_raster)


  raw_travel <- rasterize(data.frame(test$x,test$y), access.raster, test$mean_diff, fun=mean)
  raw_travel<-crop(raw_travel, extent(access.raster ))
  raw_int_travel<<-raster::mask(raw_travel, access.raster)
  return(raw_int_travel)
}

#This function reads in all the .RDS files produced by "generate_new_infrastructure_split" in a given 
#subdirectory (ROutputs) in "dir" and creates the country raster of population-weighted reduction in travel times, 
#in minutes ("weighted_int_travel")
# -------------------------------------------------------------------------------------------------------------
#Set country of interest as "country", and set directory with saved .RDS files as "dir"
#Specify a raster file (5 x 5km resolution) with binary classification of presence of one or more pathogens
#Can use binary_threshold function to get binary map from continuous map
#Requires a raster of world populations("worldpop") by 1x1-km grid-cells, cropped to country of interest
# -------------------------------------------------------------------------------------------------------------
generate_new_infrastructure_weighted<-function(country, dir, pathogen_map, worldpop){
  in_dir<-paste0(dir, "/ROutputs/")
  setwd(in_dir)
  temp<-list.files(pattern=paste0("*",country,".rds"))
  test <- do.call(rbind, lapply(temp, function(x) readRDS(x)))
  
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
  healthfacility<- fread(paste0(dir,"/who_facility_list.csv"))
  #Convert Latitude and Longitude to Numeric values and remove any missing data
  healthfacility[, Lat := as.numeric(Lat)]
  healthfacility[, Long := as.numeric(Long)]
  healthfacility <- healthfacility[!is.na(Lat) & !is.na(Long)]
  setnames(healthfacility, c('Lat', 'Long'), c('Latitude', 'Longitude'))
  #Subset to specified country
  healthcare<- healthfacility[Country == country]
  #Convert the health facilities to matrix of points for use with least cost function
  dataset <- healthcare
  coordinates(dataset) <- ~ Longitude + Latitude
  proj4string(dataset) <- proj4string(loc_shp)
  points<-as.matrix(dataset@coords)
  # create initial access raster of time to most accessible health facility using least cost algorithm
  access.raster<-accCost(T.GC, points)
  #interpolate access.raster to 5km by 5km and crop to VHF map
  access_raster <- crop(projectRaster(access.raster, pathogen_map, method = 'bilinear'), extent(pathogen_map))
  #Mask out any areas outside of the country
  access.raster <- only_country(access_raster, access_raster)
  
  
  weighted_travel <- rasterize(data.frame(test$x,test$y), access.raster, test$ptt, fun=mean)
  weighted_travel<-crop(weighted_travel, extent(access.raster))
  weighted_int_travel<<-raster::mask(weighted_travel, access.raster)
  return(weighted_int_travel)
}


