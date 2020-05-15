#########################################################################################################################
###################################### LISTS NEEDED TO RUN SDMS ---- ####################################################
#########################################################################################################################





## 1). CREATE OCCURRENCE LISTS ===========================================================================================


## GBIF columns to keep ----
gbif_keep <- c(## TAXONOMY
  "searchTaxon",
  "species",
  "scientificName",
  "taxonRank",
  "taxonKey",
  "genus",
  "family",

  ## CULTIVATION
  "cloc",
  "basisOfRecord",
  "locality",
  "establishmentMeans",
  "institutionCode",
  "datasetName",
  "habitat",
  "eventRemarks",

  ## RECORD ID
  "recordedBy",
  "identifiedBy",
  "gbifID",
  "catalogNumber",

  ## PLACE/TIME
  "lat",
  "lon",
  "decimalLatitude",
  "decimalLongitude",
  "country",
  "coordinateUncertaintyInMeters",
  "geodeticDatum",
  "year",
  "month",
  "day",
  "eventDate",
  "eventID")





## ALA columns to keep ----
ALA_keep <- c(## TAXONOMY
  "searchTaxon",
  "scientificName",
  "scientificNameOriginal",
  "species",
  "taxonRank",
  "rank",
  "genus",
  "family",

  ## CULTIVATION
  "occCultivatedEscapee",
  "basisOfRecord",
  "locality",
  "establishmentMeans",
  "institutionCode",
  "datasetName",
  "habitat",
  "eventRemarks",
  "taxonomicQuality",

  ## RECORD ID
  "recordedBy",
  "id",
  "catalogNumber",
  "identificationID",
  "identifiedBy",
  "occurrenceID",
  "basisOfRecord",
  "institutionCode",

  ## PLACE/TIME
  "latitude",
  "longitude",
  "lat",
  "lon",
  "coordinateUncertaintyInMetres",
  "coordinateUncertaintyInMeters",
  "zeroCoordinates",
  "country",
  "state",
  "IBRA7Regions",
  "IBRA.7.Subregions",
  "localGovernmentAreas",
  "locality",
  "geodeticDatum",
  "year",
  "month",
  "day",
  "eventDate",
  "eventID")




## 2). CREATE RASTER LISTS ===========================================================================================



## These rasters could change, but the names in the projections, etc, would also need to change


## Create the variables needed to access current environmental conditions + their names in the functions
## Names of all the worldclim variables used to extract the raster data
env_variables = c("Annual_mean_temp",
                  "Mean_diurnal_range",
                  "Isothermality",
                  "Temp_seasonality",
                  "Max_temp_warm_month",
                  "Min_temp_cold_month",
                  "Temp_annual_range",
                  "Mean_temp_wet_qu",
                  "Mean_temp_dry_qu",
                  "Mean_temp_warm_qu",
                  "Mean_temp_cold_qu",

                  "Annual_precip",
                  "Precip_wet_month",
                  "Precip_dry_month",
                  "Precip_seasonality",
                  "Precip_wet_qu",
                  "Precip_dry_qu",
                  "Precip_warm_qu",
                  "Precip_col_qu")

bioclim_variables = c('bio_01',
                      'bio_02',
                      'bio_03',
                      'bio_04',
                      'bio_05',
                      'bio_06',
                      'bio_07',
                      'bio_08',
                      'bio_09',
                      'bio_10',
                      'bio_11',

                      ## Rainfall
                      'bio_12',
                      'bio_13',
                      'bio_14',
                      'bio_15',
                      'bio_16',
                      'bio_17',
                      'bio_18',
                      'bio_19')


## Names of the sdm data table ---
sdm_table_vars <- c('searchTaxon', 'lon', 'lat', 'SOURCE',

                    'Annual_mean_temp',     'Mean_diurnal_range',  'Isothermality', 'Temp_seasonality',
                    'Max_temp_warm_month',  'Min_temp_cold_month', 'Temp_annual_range', 'Mean_temp_wet_qu',
                    'Mean_temp_dry_qu',     'Mean_temp_warm_qu',   'Mean_temp_cold_qu',

                    'Annual_precip',        'Precip_wet_month',    'Precip_dry_month',  'Precip_seasonality',
                    'Precip_wet_qu',        'Precip_dry_qu',       'Precip_warm_qu',    'Precip_col_qu')


## Names of the best 15 worldclim predictors ----
## i.e. 'backwards selected' predictors
bs_predictors <- c("Annual_mean_temp",    "Mean_diurnal_range",  "Isothermality",      "Temp_seasonality",
                   "Max_temp_warm_month", "Min_temp_cold_month", "Temp_annual_range",
                   "Mean_temp_warm_qu",   "Mean_temp_cold_qu",

                   "Annual_precip",       "Precip_wet_month",   "Precip_dry_month",    "Precip_seasonality",
                   "Precip_wet_qu",       "Precip_dry_qu")




## Just get the 6 models picked by CSIRO for Australia, for 2030, 2050 and 2070
## See the publication for why we choose this
scen_2030 = c("mc85bi30", "no85bi30", "ac85bi30", "cc85bi30", "gf85bi30", "hg85bi30")
scen_2050 = c("mc85bi50", "no85bi50", "ac85bi50", "cc85bi50", "gf85bi50", "hg85bi50")
scen_2070 = c("mc85bi70", "no85bi70", "ac85bi70", "cc85bi70", "gf85bi70", "hg85bi70")


## Make a list of SDM columns needed ----
results_columns = c("searchTaxon",        ## From the ALA/ GBIF download code
                    "Origin",             ## native/extoic : from Anthony Manea's spreadsheet, affected by taxonomy....
                    "Family",             ## From Anthony Manea's spreadsheet, will be affected by taxonomy....

                    "Maxent_records",     ## No. records used in the SDM
                    "Aus_records",        ## No. AUS records     :: from the R workflow
                    "AOO",                ## Global Area of occurrence
                    "KOP_count",          ## Number of koppen zones each species is found in...

                    "Number_var",        ## No. maxent variables :: from Maxent code
                    "Var_pcont",         ## Maxent Variable with highest permutation importance
                    "Per_cont",          ## The permutaiton importance of that variable
                    "Var_pimp",          ## Maxent Variable with highest permutation importance
                    "Perm_imp",          ## The permutaiton importance of that variable
                    "Iterations",               ## No. iterations
                    "Training_AUC",             ## training AUC
                    "Max_tss",                  ## Maximium True skill statistic
                    "Number_background_points", ## No. background points
                    "Logistic_threshold",
                    "Omission_rate"             ## Maxent threshold)
)
