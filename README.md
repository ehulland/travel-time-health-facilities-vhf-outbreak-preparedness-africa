# Travel time to health facilities in areas of viral hemorrhagic fever outbreak potential: maps for guiding local preparedness and response - 2019

This repo contains the codes used to generate the raster files used to create the maps in the manuscript "Travel time to health facilities in areas of viral hemorrhagic fever outbreak potential: maps for guiding local preparedness and response". 

## The repo contains 4 files:

## 1) Travel_Time_Facilities_VHFs.R

This file generates the function `generate_tt` needed to produce the health facility accessibility rasters for generating the maps, with options to specify the country, directory, and pathogen of interest (Ebola, Marburg, Lassa, CCHF, or at least one VHF ("Any")). A second function `generate_tt_thresh` produces the same raster, with the option to specify the threshold used for the pathogen of interest (min for the minimum value, p5 for the 5th percentile, median for the median value, p95 for the 95th percentile, or max for the maximum value. Other auxiliary functions for full performance are also included.

![Figure 1A](../Maps/CAR_travel_time_raw.PNG)
![Figure 1B](https://raw.githubusercontent.com/ehulland/travel-time-health-facilities-vhf-outbreak-preparedness-africa/Maps/CAR_travel_time_percentage.PNG)


## 2) Travel_Time_to_AtRisk_Pixel.R

This file generates the function `generate_at_risk_areas` to produce access rasters to the most accessible grid-cell (location) with environmental suitability for the pathogen of interest ("at-risk" locations) from any other location ("not-at-risk" locations) in a country.  It also includes functions to plot these maps either overall `plot_at_risk_map_all` or restricted only to populated areas `plot_at_risk_map_pop`. 

## 3) Travel_Time_Facilities_New_Infrastructure.R

This file uses many of the components of the function generate_tt to first generate the accessibility raster for a given pathogen (we used "Any" - or at least one VHF, but could be changed) and then estimates the mean reductions in travel times from this original raster if we were to place new infrastructure in any location in a country. We consider this both overall (i.e. "raw" reduction) and we also compute the population-weighted reduction in travel times to account for population affected in each location. This file culminates by producing two maps - one plotting the raw reductions and one plotting the weighted reductions.  

## 4) Hosptial_Table_Ranking.R

This file uses the functions produced by generate_at_risk_areas to produce a table of hospitals in a country most accessible (i.e. with the lowest travel times) to the areas at-risk for a pathogen of choice. Users can specify the pathogen, country, and number of results displayed in the table (default n=25). The resulting table is color coded by the hours of travel time to match the maps produced in the manuscript, and provides the hospital name and first administrative unit. 

