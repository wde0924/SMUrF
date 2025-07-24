# Solar-Induced Fluorescence for Modeling Urban biogenic Fluxes (SMUrF) Model
# 500m resolution branch

## Descriptions:
Scripts and subroutines for the SIF-based biospheric model which offer an appliable solution to NEE over urban areas around the globe and help separate biospheric fluxes from anthropogenic emissions. Hourly NEE fluxes are available at 500 m x 500 m grid spacing. 

Methodology is based on [*Wu et al*., 2021] and [*Madsen-Colford et al*.,submitted] . Please contact Sabrina Madsen-Colford (smadsen@physics.utoronto.ca) if you have any comments. Thank you.

## Details:
Users can start with 'main_script_*.r' for model and parameter initializations. Please refer to SMUrF_instructions.pdf and [*Madsen-Colford et al*.,submitted] for required input datasets. 

Features:
1. Estimate 8-day mean GPP based on downscaled TROPOMI SIF (*Turner et al*., 2020) and biomes-dependent GPP-SIF slopes; with biome-specific uncertainties via model-data comparisons based on FLUXNET2015. Weighted mean GPP-SIF slopes are calculated for crop and urban areas based on the estimtated C3:C4 ratio and land fractions of urban vegetation types. 

2. Estimate daily mean Reco based on pre-trained neural network (NN) model and explantory variables including air & soil temperatures and SIF-based GPP; with biome-specific uncertainties from NN model performances. 

3. Estimate hourly mean NEE by downscaling GPP and Reco (following Fisher et al., 2016) using reanalysis and data assimilation products

## References:
Wu, D., Lin, J. C., Duarte, H. F., Yadav, V., Parazoo, N. C., Oda, T., and Kort, E. A.: A model for urban biogenic CO2 fluxes: Solar-Induced Fluorescence for Modeling Urban biogenic Fluxes (SMUrF v1), Geosci. Model Dev., 14, 3633–3661, https://doi.org/10.5194/gmd-14-3633-2021, 2021. 

Madsen-Colford, S., Hutyra, L., Smith, I., Wu, D., Arain, M. A., Staebler, R., Ma, W., Restrepo-Coupe, N., and Wunch, D.: Modification and comparison of two urban vegetation models over Southern Ontario, Canada, submitted to _JGR Biogeosciences_.
