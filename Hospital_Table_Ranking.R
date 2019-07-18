################################################################################################################################################################################
# Calculating Ranked Travel Times to Hospitals closes to Ebola cases in 2018-2019 Outbreak
################################################################################################################################################################################

# Data used can be obtained from the following sources: 
# Country, Administrative unit shapefiles: downloaded from Malaria Atlas Project's 'getShp(country, ISO, extent)' function in the MalariaAtlas Project package
# Population Mask: WorldPop project https://www.worldpop.org/project/list
# Friction Surface Raster: https://map.ox.ac.uk/ accessed via getRaster(surface = "A global friction surface enumerating land-based travel speed for a nominal year 2015")
# Lakes Shapefile: http://193.43.36.146/map?entryId=bd8def30-88fd-11da-a88f-000d939bc5d8
# Viral Hemorrhagic Fever Environmental Suitability Maps: Pigott, DM et. al 2017 article https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(17)32092-5/fulltext

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
#need the results of the following file for use: 
source('Travel_Time_to_AtRisk_Pixel.R')

# List of countries for analysis
countries<-c("Angola", "Benin", "Botswana", "Burkina Faso", "Burundi", "Cameroon", "Central African Republic",
             "Chad", "Congo", "Democratic Republic of the Congo", "Djibouti", "Equatorial Guinea",
             "Eritrea", "Ethiopia", "Gabon", "Gambia", "Ghana", "Guinea", "Guinea Bissau", "Kenya",
             "Lesotho", "Liberia", "Madagascar", "Malawi", "Mali", "Mauritania", "Mozambique",
             "Namibia", "Niger", "Nigeria", "Rwanda", "Senegal","Sierra Leone", "Somalia", "South Africa",
             "South Sudan", "Sudan", "Tanzania", "Togo", "Uganda", "Zambia", "Zimbabwe")
# -------------------------------------------------------------------------------------------------------------
#set country, directory, and number of hospitals of interest for the table (num.rows)
num.rows<-25
# -------------------------------------------------------------------------------------------------------------
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
#Generate risk points from source file to get time to nearest at-risk location for pathogen of interest and country of interest
#Default Pathogen is "Any" - at least one VHF per grid-cell
#Other choices are Ebola, Lassa, Marburg or CCHF
generate_at_risk_areas(country<-country,dir<-dir,pathogen<-"Any")
VHFrisk_access<-accCost(T.GC, risk)
#set any "infinte" values (values outside of the country) to NA
VHFrisk_access[VHFrisk_access==Inf]<-NA
#Convert minutes to hours
values(VHFrisk_access)<-values(VHFrisk_access)/60

#Read in WHO facility data and subset to the identified countries above
WHO <- fread(paste0(dir, "/who_facility_list.csv"))
healthfacility <-subset(WHO, WHO$Country %in% countries)

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

#Change Guinea-Bissau back again for the health facility data
if(country=="Guinea-Bissau"){country="Guinea Bissau"}
healthcare<- healthfacility[Country == country]

#restrict to Hospitals only
hospitals<-healthcare[FacilityType=="Hospital"]
coordinates(hospitals) <- ~ Longitude + Latitude
proj4string(hospitals) <- proj4string(loc_shp)
points<-as.matrix(hospitals@coords)

#get the raster values for travel time at each of the locations of the hospitals and bind back to the locations / names; remove any missing data
TT <- raster::extract(VHFrisk_access, points)
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
#Rescale to same scale as the maps (0 to 24+ hours)
hosp$col<-ifelse(hosp$Travel_Time>24, 24, hosp$Travel_Time)
#round to two decimal places for printing
hosp$Travel_Time<-round(hosp$Travel_Time,2)
#Change the color for the travel time label based on the color of the tile (<=13 hours versus more than 13 hours)
hosp$lab<-ifelse(hosp$col>13, "#FFFFFF", "#000004FF")
hosp$name_admin<-reorder(hosp$name_admin, -hosp$Travel_Time)
list<-hosp@data

#Switch name back one last time
if(country=="Guinea Bissau"){country="Guinea-Bissau"}

#set the color for the tiles based on the number of hours as in the maps
colors<-magma(16, dir=-1)
lims=round(seq(0,24, length.out=16), 0)
#Generate table with colors
p <- ggplot(data = list, aes(x = Country, y = name_admin)) +
  geom_tile(aes(fill =col), color = 'white', size = 0.2) +
  scale_fill_gradientn(colors = colors, limits=c(0,24), breaks=lims, labels=c("0"," ", " ", "5", " ", " ", "10", " ", " ", "14", " ", " ","19",  " ", " ", "24"),name="Travel Time (Hours)") +
  geom_text(aes(col=lab, label = Travel_Time), size = 2, fontface='bold') +
  scale_colour_manual(values=c("#000004FF","#FFFFFF"), guide=FALSE)+
  ggtitle(country) +
  labs(x = ' ', y = 'Hospital')+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) 
p
