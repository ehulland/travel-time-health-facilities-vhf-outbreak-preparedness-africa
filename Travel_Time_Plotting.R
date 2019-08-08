################################################################################################################################################################################
# Code to make plots via rasterVis
################################################################################################################################################################################
# -------------------------------------------------------------------------------------------------------------
# Source functions file
source('~/Travel_Time_Functions.R')
#Generating plot of travel times in hours to most accessible facility from locations with pathogen spillover potential
#using "pathogen_access_raster" generated from "generate_tt" function
#First use mask_lake function if a shapefile of lakes exists and should be masked out for plotting
mask_lake(country,dir, lakes)
generate_tt(country,dir,pathogen_map, healthcare)
# -------------------------------------------------------------------------------------------------------------
#Set country of interest as "country", and set directory with saved shapefiles and rasters as "dir"
#Requires a database of health facilities (healthcare) with (at-minimum) geocoded latitude (Latitude) and longitude (Longitude)
# -------------------------------------------------------------------------------------------------------------
plot_tt_hrs<-function(country, dir, pathogen_access_raster, healthcare){
out_shp<-getShp(country<-country)
#Read in lakes shapefile if it exists
if(file.exists(paste0(dir, "/LakesShapefile.shp"))){
  loc_shp<-readOGR(dsn=paste0(dir, "/LakesShapefile.shp"))
} else {
  loc_shp<-out_shp
}
#Get 1st administrative unit borders for plotting from malariaAtlas "getShp"
admin1<-getShp(country=country, admin_level=c("admin1"))
#Get neighboring country outlines from malariaAtlas "getShp"
#Extend extent a bit from the shapefile outline
ext<-extent(out_shp)
ext_mat<-matrix(c(ext[1]-1,ext[3]-1,ext[2]+1,ext[4]+1), nrow=2, ncol=2, dimnames=list(c("x","y"), c("min", "max")))
admin0<-getShp(country="ALL", extent=ext_mat, admin_level="admin0")

#prep healthcare dataset to plot
dataset <- healthcare
coordinates(dataset) <- ~ Longitude + Latitude
proj4string(dataset) <- proj4string(loc_shp)

# set legend title
legend_title <- 'Travel time (Hours) \n'
plot_title <- paste0('Travel Time to Most Accessible Health \n Facility in ', country)

#convert minutes to hours for pathogen_access_raster
values(pathogen_access_raster)<-values(pathogen_access_raster)/60
#set  legend thresholds - 0 to 24+ hours
pathogen_access_raster2<-pathogen_access_raster
values(pathogen_access_raster2)<-ifelse(values(pathogen_access_raster2)>24 & values(pathogen_access_raster2)<Inf, 24, values(pathogen_access_raster2))
lims=seq(0,24, length.out=16)

#Set plotting boundaries
xlims<-c(ext[1]-1, ext[2]+1)
ylims<-c(ext[3]-1, ext[4]+1)

#crop admin0 to the limits above
admin0sub<-crop(admin0, extent(xlims, ylims))
#get standard color scheme using magma color palette and set background to light blue for lakes
Rev.Magma<-rasterTheme(region=magma(16,direction=-1), panel.background=list(col="light blue"))
#plot!
t<-levelplot(access_rasterNP2,par.settings=Rev.Magma,at=lims, margin=FALSE, xlim=xlims, ylim=ylims,ylab="", xlab=list(label=legend_title, cex=1.0, line=1),scales = list(draw = FALSE), colorkey=list(space='bottom'),
             main=list(plot_title,side=1,line=0.5,cex=1.2))
t_plot<-t+layer_(sp.polygons(loc_shp, col="white",fill="white"))+layer_(sp.polygons(out_shp, col=NA, fill="light blue"))
tt_plot<-t_plot+layer(sp.polygons(admin1, lwd=1, lty=2))+layer_(sp.polygons(admin0sub, lwd=2, fill="white"))
tt.plot<<-tt_plot+layer(sp.points(dataset, cex=0.15, col="dark green", pch=20, alpha=0.5))
return(tt.plot)
}

#Generating plot of travel times in percentiles to most accessible facility from locations with pathogen spillover potential
#using "pathogen_access_raster" generated from "generate_tt" function
# -------------------------------------------------------------------------------------------------------------
#Set country of interest as "country", and set directory with saved shapefiles and rasters as "dir"
#Requires a database of health facilities (healthcare) with (at-minimum) geocoded latitude (Latitude) and longitude (Longitude)
# -------------------------------------------------------------------------------------------------------------
plot_tt_pct<-function(country, dir, pathogen_access_raster, healthcare){
  out_shp<-getShp(country<-country)
  #Read in lakes shapefile if it exists
  if(file.exists(paste0(dir, "/LakesShapefile.shp"))){
    loc_shp<-readOGR(dsn=paste0(dir, "/LakesShapefile.shp"))
  } else {
    loc_shp<-out_shp
  }
  #Get 1st administrative unit borders for plotting from malariaAtlas "getShp"
  admin1<-getShp(country=country, admin_level=c("admin1"))
  #Get neighboring country outlines from malariaAtlas "getShp"
  #Extend extent a bit from the shapefile outline
  ext<-extent(out_shp)
  ext_mat<-matrix(c(ext[1]-1,ext[3]-1,ext[2]+1,ext[4]+1), nrow=2, ncol=2, dimnames=list(c("x","y"), c("min", "max")))
  admin0<-getShp(country="ALL", extent=ext_mat, admin_level="admin0")
  
  #prep healthcare dataset to plot
  dataset <- healthcare
  coordinates(dataset) <- ~ Longitude + Latitude
  proj4string(dataset) <- proj4string(loc_shp)
  
  # set legend title
  legend_title <- 'Travel time (Percentile) \n'
  plot_title <- paste0('Travel Time to Most Accessible Health \n Facility in ', country)

  #set  legend thresholds - 0 to 100%
  pathogen_access_raster[pathogen_access_raster==Inf]<-NA
  values(pathogen_access_raster)<-values(pathogen_access_raster)/maxValue(pathogen_access_raster)
  lims=quantile(pathogen_access_raster, probs=c(0,0.2,0.4,0.6,0.8,1))
 
   #Set plotting boundaries
  xlims<-c(ext[1]-1, ext[2]+1)
  ylims<-c(ext[3]-1, ext[4]+1)
  
  #crop admin0 to the limits above
  admin0sub<-crop(admin0, extent(xlims, ylims))
  #get standard color scheme using viridis color palette and set background to light blue for lakes
  Rev.Magma<-rasterTheme(region=viridis(5,direction=-1), panel.background=list(col="light blue"))
  myColorkey <- list(space='bottom',at=c(0,20,40,60,80,100),
                     labels=list(
                       labels=c("0-20%","20-40%","40-60%","60-80%", "80-100%"),
                       at=c(10,30,50,70,90)
                     ))
  #plot!
  t_pct<-levelplot(pathogen_access_raster,par.settings=Rev.Magma,at=lims, margin=FALSE, xlim=xlims, ylim=ylims,ylab="", xlab=list(label=legend_title, cex=1.0, line=1),scales = list(draw = FALSE), colorkey=myColorkey,
               main=list(plot_title,side=1,line=0.5,cex=1.2))
  t_plot_pct<-t_pct+layer_(sp.polygons(loc_shp, col="white",fill="white"))+layer_(sp.polygons(out_shp, col=NA, fill="light blue"))
  tt_plot_pct<-t_plot_pct+layer(sp.polygons(admin1, lwd=1, lty=2))+layer_(sp.polygons(admin0sub, lwd=2, fill="white"))
  tt.plot.pct<<-tt_plot_pct+layer(sp.points(dataset, cex=0.15, col="red", pch=20, alpha=0.3))
  return(tt.plot.pct)
}

