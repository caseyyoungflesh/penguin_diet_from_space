# penguin_diet_from_space


## Overview

Repository containing code to quantify penguin diet from satellite sensors, how diet is responding to changing abiotic conditions, and how these trophic dynamics are linked to long-term penguin population change.


## Associated publications

Youngflesh, C, C Che-Castaldo, M Schwaller, M Polito, S Serbin, HJ Lynch. Penguin diet from space: Links between sea ice, Antarctic food webs, and population change. Submitted.


## Workflow

The workflow for this repository involves data read in as raw data at Level 0 (L0). Level 1 (L1) data represent cleaned L0 data. Level 2 (L2) data are data merged or otherwise derived from two or more L1 data, etc. Data must be organized according to the below Repository structure (i.e., in `Data/`). Scripts are designed to be run in sequence from `1-XXXX.R` -> `2-XXXX.R` -> ... The directory where the data are located (i.e., the parent directory of `Data/`) as well as the run date for each data product need to be specified at the top of associated scripts (in the `set dirs` section) before running. In each scipt, the `here` package automatically uses the Project directory (RStudio) as the root directory, though this could be specified manually.


## Repository structure

* `Scripts/`
  * `1-process-abiotic-data/`
    * `1a-process-locations.R` - get colony locations
    * `1b-process-SIC.R` - process sea ice data
    * `1c-process-bathymetry.R` - process bathymetry data
  * `2-PLSR.R` - run PLSR model, diet ~ spectra
  * `3-predict-diet.R` - predict diet based on PLSR model fit and Landsat spectra
  * `4-diet-est.R` - model to estimate diet at each site/year/day
  * `5-diet-model.R` - model to estimate for diet varies according to predictors
  * `6-gr-model.R` - model to estiamte effect of diet on pop growth rates
 
* `Data/` (ignored)
  * `L0/`
    * `Bathymetry_data/`
      * `IBCSO_v2_bed.tif` - bathymetry
    * `Landsat/`
      * `landsat_reflectance_L457-2025-07-16.rds` - reflectance values from penguin colonies
    * `SIC/` - sea ice data
    * `RSR/` - spectral response curves
      * `L4_TM_RSR.xlsx`
      * `L5_TM_RSR.xlsx`
      * `L7_ETM_RSR.xlsx`
    * `SiteLocations.csv` - site locations, produced with `1a-process-locations.R`
  * `L1/`
    * `SIA_spectra.csv` - SIA and spectra for guano samples
    * `colony_sic.csv` - colony SIC, produced with `1b-process-SIC.R`
    * `shelf_area.csv` - colony shelf area, produced with `1c-process-bathymetry.R`
  * `L2/`
    * `YYYY-MM-DD/` - run date
      * `L4_coef_en.rds` - coefficients from PLSR model, produced with `2-PLSR.R`
      * `L5_coef_en.rds` - coefficients from PLSR model, produced with `2-PLSR.R`
      * `L7_coef_en.rds` - coefficients from PLSR model, produced with `2-PLSR.R`
  * `L3/`
     * `diet_pred-YYYY-MM-DD.rds` - predicted diet, produced with `3-predict-diet.R`
* `Results/` (ignored)


## Contact Information

Casey Youngflesh - cyoungf@clemson.edu
