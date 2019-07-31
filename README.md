# Travel time to health facilities in areas of viral hemorrhagic fever outbreak potential: maps for guiding local preparedness and response - 2019

This repo contains the codes used to generate the raster files used to create the maps in the manuscript "Travel time to health facilities in areas of viral hemorrhagic fever outbreak potential: maps for guiding local preparedness and response". 

## The repo contains 2 files:

## 1) Travel_Time_Functions.R

This file generates the functions needed to produce the rasters used in the manuscript. Each of the functions is described subsequently:

`binary_threshold` : Converts a raster with continuous values into dichotomous presence (1) or absence (0) based on a user-defined threshold

`gen_mask` : Takes two rasters (or the same raster, twice) as inputs, and converts the values of the first raster to NA where the values are 0 in the second raster

`only_country` : Takes two rasters (or the same raster, twice) as inputs, and masks out values of the first raster to NA where the values are Inf in the second

`mask_lake` : Requires a country, directory, and shapefile of lakes (usually >100km) cropped to the country of interest, and returns a shapefile with the lakes cropped out in the directory provided

`generate_tt` : Requires a specified country and directory as well as a database of health facilities with latitude (named Latitude, numeric) and longitude (named Longitude, numeric), at minimum (but could contain other columns), and a raster file at 5x5-km resolution with binary classification of presence / absence of one or more pathogens. Can use `binary_threshold` to produce the binary classification, and may want to consider using `mask_lake` to produce a shapefile with lakes removed for calculating travel times (function will search for this file if it exists, but will not throw an error if it doesn't). This function returns raster *pathogen_access_raster* (5x5-km resolution) of travel times to the most accessible health facility from any grid-cell at risk for the pathogen specified. 

Examples of plots made with *pathogen_access_raster* from the manuscript include the following raw and percentile ranked travel times to health facilities from locations at-risk for at least one VHF in Central African Republic:

![Figure 1A](Maps/CAR_travel_time_raw.PNG)

![Figure 1B](Maps/CAR_travel_time_percentage.PNG)


## 2) Travel_Time_to_AtRisk_Pixel.R

This file generates the function `generate_at_risk_areas` to produce access rasters to the most accessible grid-cell (location) with environmental suitability for the pathogen of interest ("at-risk" locations) from any other location ("not-at-risk" locations) in a country.  It also includes functions to plot these maps either overall `plot_at_risk_map_all` or restricted only to populated areas `plot_at_risk_map_pop`. 

Example plots generated include travel time to locations at-risk for at least one VHF both internal to Botswana and international (using a buffer of 500km) via cross-border travel or migration: 

![Figure 2A](Maps/Botswana_travel_time_atrisk_inner.PNG)

![Figure 2B](Maps/Botswana_travel_time_atrisk_outer.PNG)


## 3) Travel_Time_Facilities_New_Infrastructure.R

This file uses many of the components of the function generate_tt to first generate the accessibility raster for a given pathogen (we used "Any" - or at least one VHF, but could be changed) and then estimates the mean reductions in travel times from this original raster if we were to place new infrastructure in any location in a country. We consider this both overall (i.e. "raw" reduction) and we also compute the population-weighted reduction in travel times to account for population affected in each location. This file culminates by producing two maps - one plotting the raw reductions and one plotting the weighted reductions.  

Example plots generated include the raw and weighted reductions in travel time for Ethiopia based on the current health facility profile in-country:

![Figure 3A](Maps/Ethiopia_travel_time_reduction_raw.PNG)

![Figure 3B](Maps/Ethiopia_travel_time_reduction_weighted.PNG)


## 4) Hosptial_Table_Ranking.R

This file uses the functions produced by generate_at_risk_areas to produce a table of hospitals in a country most accessible (i.e. with the lowest travel times) to the areas at-risk for a pathogen of choice. Users can specify the pathogen, country, and number of results displayed in the table (default n=25). The resulting table is color coded by the hours of travel time to match the maps produced in the manuscript, and provides the hospital name and first administrative unit. Here we provide a demonstration for Congo:

![Figure 4](Maps/Congo_ranked_hospital_list.PNG)

