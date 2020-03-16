#########################################################################################################################
#####################################  RUNNING SDM ANALYSIS FOR TAXA ---- ###############################################
#########################################################################################################################



## This code runs an SDM work flow, for a subset of species (e.g. whichever you supply)
## See (Burley et al 2019) for more details on the method ::
## https://www.sciencedirect.com/science/article/pii/S0048969719323289


## When run as a package, the addresses wouldn't be the same...
## The functions run over the folders, and they will be hard to replicate


## Raster data is too big to push to github
## Changing rasters to .rda + addresses may cause errors.
## EG the plot raster doesn't work simply on the .rda.
## So



## SETUP DATA & ENVIRONMENT ===========================================================================================


## First, load the data needed to run the analysis
## This file is created by running (./R/READ_SPATIAL_DATA.R)


## Set the R library so we know all the packages and dependencies will be the same
# .libPaths("./R/win-library/3.5")


## Load only the packages needed for the analysis
package.list <- c('ff',         'things',    'raster',        'dismo',             'rJava',        'sp',                    'latticeExtra',
                  'data.table', 'devtools',  'roxygen2',
                  'rgdal',      'rgeos',     'gdalUtils',     'rmaxent',           'readr',        'plyr',                  'dplyr',
                  'tidyr',      'readr',     'rnaturalearth', 'rasterVis',         'RColorBrewer', 'latticeExtra',          'parallel',
                  'ALA4R',      'stringr',   'Taxonstand',    'CoordinateCleaner', 'gsubfn',       'PerformanceAnalytics',  'utf8',
                  'rvest',      'magrittr',  'devtools',      'ggplot2',           'reshape2',     'rmarkdown',             'flexdashboard',
                  'shiny',      'rgbif',     'ENMeval',       'tibble',            'ncdf4',        'Cairo',                 'taxonlookup',
                  'kgc',        'betareg',   'gridExtra',     'grid',              'lattice',      'ConR',                  'writexl',
                  'sf',         'ggmap',     'DataCombine',   'exactextractr',     'mgcv', 'doSNOW', 'tidyverse',  'ggpubr', 'GGally')


## Require packages
sapply(package.list, require, character.only = TRUE)


## Source functions, and set temporary directory
source('./R/SDM_GEN_PROCESSOR_FUNCTIONS.R')
source('./R/SDM_GEN_MAXENT_FUNCTIONS.R')
source('./R/SDM_GEN_MAPPING_FUNCTIONS.R')
source('./R/SDM_GEN_MODEL_LISTS.R')
# source('./R/READ_SPATIAL_DATA.R')


## Set temporary raster directory
#rasterOptions(tmpdir = './RTEMP')





## STEP 1 ======================================================================================================


## Load in all the .rda files that are save with the package ::
## this is all the data needed to run an example analysis
list.filenames <- list.files(path = './data/', pattern = ".rda$", full.names = TRUE)
list.data<-list()

## See the global environemnt for a list of objects
for (i in 1:length(list.filenames))

{
  #message('load .rda for ', i)
  list.data[[i]]<-load(list.filenames[i])
}


## Create a list of species to analyse
bat.spp <- aus.bats$Species[1:2] %>%
  str_trim()



## :: Download GBIF and ALA data
download_GBIF_all_species(species_list = bat.spp,
                          path         = "./data/GBIF/")


download_ALA_all_species(species_list = bat.spp,
                         path         = "./data/ALA/")





## STEP 2 ======================================================================================================



## :: Combine ALA data and filter to records on land taken > 1950
## The climate data is the worldclim version 1.0
## Raster :: Error in .local(.Object, ...) :
ALA.LAND = combine_ala_records(species_list      = bat.spp,
                               records_path      = "./data/ALA/",
                               records_extension = "_ALA_records.RData",
                               record_type       = "ALA",
                               keep_cols         = ALA.keep,
                               world_raster      = world.grids.current)


## :: Combine GBIF data and filter to records on land taken > 1950
GBIF.LAND = combine_gbif_records(species_list      = bat.spp,
                                 records_path      = "./data/GBIF/",
                                 records_extension = "_GBIF_records.RData",
                                 record_type       = "GBIF",
                                 keep_cols         = gbif.keep,
                                 world_raster      = world.grids.current)





## STEP 3 ======================================================================================================


## :: combine GBIF and ALA data, and extract environmental values
## Create messages for this
COMBO.RASTER.CONVERT = combine_records_extract(ala_df          = ALA.LAND,
                                               gbif_df         = GBIF.LAND,
                                               template_raster = template.raster.1km.84,
                                               world_raster    = world.grids.current,
                                               projection      = CRS.WGS.84,
                                               species_list    = bat.spp,
                                               biocl_vars      = bioclim.variables,
                                               env_vars        = env.variables,
                                               worldclim_grids = "TRUE",
                                               save_data       = "FALSE",
                                               save_run        = "TEST_BATS")





## STEP 4 ======================================================================================================


## :: Flag records as institutional or spatial outliers
COORD.CLEAN = coord_clean_records(records    = COMBO.RASTER.CONVERT,
                                  capitals   = 10000,  ## Remove records within 10km  of capitals
                                  centroids  = 5000,   ## Remove records within 5km of country centroids
                                  save_run   = "TEST_BATS",
                                  save_data  = "FALSE")


## Step 4b :: Flag spatial outliers
check_spatial_outliers(all_df      = COORD.CLEAN,
                       land_shp    = LAND,
                       clean_path  = './data/GBIF/',
                       spatial_mult = 10)


## Step 4c ::Estima ecliate niches usign pecies records
GLOB.NICHE = calc_1km_niches(coord_df     = COORD.CLEAN,
                             aus_shp      = AUS,
                             world_shp    = LAND,
                             kop_shp      = Koppen_shp,
                             ibra_shp     = IBRA,
                             species_list = bat.spp,
                             env_vars     = env.variables,
                             cell_size    = 2,
                             save_run     = "TEST_BATS",
                             save_data    = "FALSE")


## Step 4d :: plot species ranges using histograms and convex hulls
plot_range_histograms(coord_df     = COORD.CLEAN,
                      species_list = bat.spp,
                      range_path = './data/GBIF/')





## STEP 6 :: Flag spatial outliers ===============================================================================================


## This code creates the table for running SDMs in 'species with data', or SWD format.
## There is a switch in the function, that adds additional bg points from other taxa, if specified.
## In this example for bats, we'll just use
SDM.SPAT.OCC.BG = prepare_sdm_table(coord_df        = COORD.CLEAN,
                                    species_list    = unique(COORD.CLEAN$searchTaxon),
                                    sdm_table_vars  = sdm.table.vars,
                                    save_run        = "TEST_BATS",
                                    read_background = "FALSE",
                                    #BG_points       = 'SDM_SPAT_ALL_ANIMAL_BG_POINTS.rds', ## This should be just bats
                                    save_data       = "FALSE",
                                    save_shp        = "FALSE")





## STEP 7 ::: run SDMs ======================================================================================================
## Run SDMs across all taxa, but don't show extra loops, etc.


## This function takes a table of all species occurrences (rows) and environmental values (columns), and runs SDMs for a list
## of taxa. First, targetted background selection is used. Then, backwards selection is used on the occurrence and
## background points generated by the targetted selection. The function saves files in folders, it doesn't return anything.


## Here we can read in a file saved earlier
#SDM.SPAT.OCC.BG <- readRDS("./data/ANALYSIS/SDM_SPAT_OCC_BG_ALL_BATS.rds")


## Run the SDMs for a set of taxa ----
run_sdm_analysis(species_list            = bat.spp,
                 maxent_dir              = 'output/maxent/full_models',      ## Must be created in your directory
                 bs_dir                  = 'output/maxent/back_sel_models',  ## Must be created in your directory
                 sdm_df                  = SDM.SPAT.OCC.BG,
                 sdm.predictors          = bs.predictors,
                 backwards_sel           = "TRUE",        ## If TRUE, run both full models, and backwards selected models
                 cor_thr                 = 0.8,           ## The maximum allowable pairwise correlation between predictor variables
                 pct_thr                 = 5,             ## The minimum allowable percent variable contribution
                 k_thr                   = 4,             ## The minimum number of variables to be kept in the model.

                 min_n                   = 20,            ## This should be higher...
                 max_bg_size             = 70000,         ## could be 50k or lower, it just depends on the biogeography
                 background_buffer_width = 200000,        ## How far away should BG points be selected?
                 shapefiles              = TRUE,          ## Save the shapefiles?
                 features                = 'lpq',         ## Maxent features to use
                 replicates              = 5,             ## Number of replicates
                 responsecurves          = TRUE,          ## Save response curves?
                 # shp_path                = "./data/base/Contextual/", ## Path for shapefile, eg Australia's border
                 aus_shp                 = AUS,
                 Koppen                  = Koppen_1975_1km)





## STEP 8 ::: project SDMs ======================================================================================================


## Create a table of maxent results ----
## This function will just aggregate the results for models that ran successfully
MAXENT.RESULTS = compile_sdm_results(species_list = bat.spp,
                                     results_dir  = 'output/maxent/back_sel_models',
                                     save_data    = FALSE,
                                     DATA_path    = "./output/results/",
                                     save_run     = "TEST_BATS")


## First, get map_spp from the maxent results table above, change the species column, and create a list of logistic thresholds -
## These values are debateable, need a few references from the literature...
map_spp         <- MAXENT.RESULTS$searchTaxon %>% gsub(" ", "_", .,)
percent.10.log  <- MAXENT.RESULTS$Logistic_threshold
sdm.results.dir <- MAXENT.RESULTS$results_dir


## Then project the SDMS across the whole of Australia
## This function takes the output of the Maxent species distribution models using current conditions, and generates
## a prediction of habitat suitability for current and future environmental conditions.


## Use a function to convert a rasterized shapefile into a vector, based on a field from one of the shapefiles
## We are just using the Significant Urban Areas of Australia as an example
## This could be any shapefile (e.g. the IBRA regions, etc.)


## Create 2030 sdm map projections ----
## This step would need to change, to run as a package
tryCatch(
  project_maxent_grids_mess(aus_shp       = AUS,          ## Shapefile, e.g. Australian states
                            world_shp     = LAND,         ## World shapefile

                            scen_list     = scen_2030,                 ## List of climate scenarios
                            species_list  = map_spp,                   ## List of species folders with maxent models
                            maxent_path   = './output/maxent/back_sel_models/',    ## Output folder
                            climate_path  = './data/worldclim/aus/',             ## Future climate data

                            grid_names    = env.variables,             ## Names of all variables
                            time_slice    = 30,                        ## Time period
                            current_grids = aus.grids.current,         ## predictor grids
                            create_mess   = TRUE,
                            nclust        = 1),

  ## If the species fails, write a fail message to file.
  error = function(cond) {

    ## This will write the error message inside the text file, but it won't include the species
    file.create(file.path("output/maxent/back_sel_models/mapping_failed_2030.txt"))
    cat(cond$message, file = file.path("output/maxent/back_sel_models/mapping_failed_2030.txt"))
    warning(cond$message)

  })





## STEP 9 ::: project SDMs ======================================================================================================



## Rasterize a shapefile ----
areal_unit_vec <- shapefile_vector_from_raster(shp_file = SUA,
                                               prj      = ALB.CONICAL,
                                               sort_var = 'SUA_NAME16',
                                               agg_var  = 'SUA_CODE16',
                                               temp_ras = aus.grids.current[[1]],
                                               targ_ras = './data/SUA_2016_AUST.tif')


## This is a vector of all the cells that either are or aren't in the rasterized shapefile
summary(areal_unit_vec)




## Combine GCM predictions and calculate gain and loss for 2030 ----
## Then loop over the species folders and climate scenarios
tryCatch(mapply(area_cell_count,                                 ## Function aggreagating GCM predictions by spatial unit
                #unit_path     = "./data/base/Contextual/",      ## Data path for the spatial unit of analysis
                unit_shp      = SUA, #"SUA_2016_AUST.rds",             ## Spatial unit of analysis - E.G. SUAs
                unit_vec      = areal_unit_vec,                  ## Vector of rasterized unit cells
                sort_var      = "SUA_NAME16",
                code_var      = "SUA_CODE16",
                world_shp     = LAND, #"LAND_world.rds",         ## Polygon for AUS maps
                aus_shp       = AUS, #"aus_states.rds",          ## Polygon for World maps

                DIR_list      = sdm.results.dir,                 ## List of directories with rasters
                species_list  = map_spp,                         ## List of species' directories
                number_gcms   = 6,
                maxent_path   = 'output/maxent/back_sel_models/', ## Directory of maxent results
                thresholds    = percent.10.log,                   ## List of maxent thresholds
                time_slice    = 30,                               ## Time period, eg 2030
                write_rasters = TRUE),

         ## If the species fails, write a fail message to file.
         error = function(cond) {

           ## This will write the error message inside the text file,
           ## but it won't include the species
           file.create(file.path("output/maxent/back_sel_models/sua_count_failed_2030.txt"))
           cat(cond$message, file=file.path("output/maxent/back_sel_models/sua_count_failed_2030.txt"))
           warning(cond$message)

         })





#########################################################################################################################
## THINGS TO CHANGE FOR BAT MODELLING ::
#########################################################################################################################



## 1). Update all the functions for bats

##     - Test the whole code for bats
##     - Check the SUA cell count section
##     - Need to reproduce the data and code in a drive
##     -
##     -


## 2). Functionalise the workflow

##     - Iron out problems in Step 8 for mapping
##     - change SUAs to IBRA regions



## 2).  Update background points to be all bats
##
##      - Use ALA data
##      - Consider the overseas data?


## 4).  Add in taxon stand cleaing from Aniko's code

##      - Clean up files on H:



## 5). How can the code be made more portable? Convert to R package, etc. The raster step is the hard part...
##     ## The raster names from the worldclim website would need to not be hard-coded





#########################################################################################################################
## OUTSTANDING WORKFLOW TASKS:
#########################################################################################################################


## 3).  Make the code more modular, less monolithic

## 3).  Iron out some of the points which are not as reproducible - e.g the background points in step 6).

## 5).  Create a 'dashboard' for each species, using the MAXENT output (eg see "./R/PLOT_RANGE_HISTOGRAMS.R')





#########################################################################################################################
###################################################### TBC ##############################################################
#########################################################################################################################
