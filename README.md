# Solar-Induced Fluorescence for Modeling Urban biogenic Fluxes (SMUrF) Model

## Descriptions:
Scripts and subroutines for the SIF-based biospheric model which offer a appliable solution to NEE over urban areas around the globe and help separate biospheric fluxes from anthropogenic emissions. Hourly NEE fluxes are available at 0.05 x 0.05 deg grid spacing. 

Methodology is based on [*Wu et al*., submitted]. Please contact Dien (dienwu@caltech.edu) if you have any comments. Thank you.

## Details:
Users can start with 'main_script_*.r' for model and parameter initializations. Please refer to Figure 1 in the manuscript for required input datasets. Most data in use are stored in SMUrF/data, while one has to download MCD12Q1 data from [here](http://home.chpc.utah.edu/~u0947337/LAIR_group/OCO-2/SMUrF/MCD12Q1.zip) or [appeears](https://lpdaacsvc.cr.usgs.gov/appeears/) to SMUrF/data. 

Features:
1. Estimate 4-day mean GPP based on CSIF (*Zhang et al*., 2018) and biomes-dependent GPP-SIF slopes; with biome-specific uncertainties via model-data comparisons based on FLUXNET2015. Weighted mean GPP-SIF slopes are calculated for crop and urban areas based on the estimtated C3:C4 ratio and land fractions of urban vegetation types. 

2. Estimate daily mean Reco based on pre-trained neural network (NN) model and explantory variables including air & soil temperatures and SIF-based GPP; with biome-specific uncertainties from NN model performances. 

3. Estimate hourly mean NEE by downscaling GPP and Reco (following Fisher et al., 2016) using reanalysis and data assimilation products


## Reference:
Dien Wu, John C. Lin, Henrique F. Duarte, Vineet Yadav, Nicholas C. Parazoo, Tomohiro Oda, and Eric A. Kort: A Model for Urban Biogenic CO2 Fluxes: Solar-Induced Fluorescence for Modeling Urban biogenic Fluxes (SMUrF), submitted to *Geoscientific Model Development*. 

