sdmgen : a pacakge for rapidly estimating multiple species ranges and
habitat suitability
================
March 2020

<style>
p.caption {
font-size: 6;
}
</style>

  

The text and code below summarises a workflow in R that can be used to
relatively rapidly assess the environmental range of a species within
Australia, from downloading occurrence records, through to creating maps
of predicted climatic suitability across Australia at 1km\*1km
resolution. An example of this work is published in Science of the Total
Environment ::

  

Burley, H., Beaumont, L.J., Ossola, A., et al. (2019) Substantial
declines in urban tree habitat predicted under climate change. Science
of The Total Environment, 685, 451-462.

<https://www.sciencedirect.com/science/article/pii/S0048969719323289#f0030>

  

To install, run :

  
  
  

# Background

This code was developed at Macquarie University in Sydney, as part of
the ‘which plant where’ project. The aim was to create a pipeline to
rapidly assess the climatic suitability of large suites of horticultural
species.

  
  
  

All over the world, local government authorities are increasing their
investment in urban greening interventions, yet there is little
consideration of whether the current palette of species for these
plantings will be resilient to climate change. This pipeline was created
to assessed the distribution of climatically suitable habitat, now and
in the future, for the tree species most commonly grown by nurseries and
planted across Australia’s urban landscapes.

  
  
  

# STEP 1 :: Download species occurrence data

  
  
  

The backbone of the R workflow is a list of (taxonomically
Ridgey-Didge\!) species names that we supply. The analysis is designed
to process data for one species at a time, allowing species results to
be updated as required. We can demonstrate the workflow using a sample
of 10 plant species from the Stoten publication above.

  
  
  

``` r
## Use the first 10 plant species in the Stoten list
data("plant.spp")
analysis_spp <- plant.spp[1:5]
analysis_spp
```

  
  
  

This workfow uses four shapefiles : Australia, the World, the global
Koppen Zones, and the Significant Urban areas of Australia (or SUAs).
The SUAs are taken from the ABS :
<https://www.abs.gov.au/AUSSTATS/abs@.nsf/DetailsPage/1270.0.55.004July%202016?OpenDocument>.
The Koppen data are from CliMond, centred on 1975 :
<https://www.climond.org/Core/Authenticated/KoppenGeiger.aspx>

  
  
  

The code also uses 1km worldclim raster data. Put this folder in your
‘data’ project folder :

  
  
  

<https://drive.google.com/open?id=1T5ET5MUX3-lkqiN5nNL3SZZagoJlEOal>

  
  
  

The species list is supplied to a series of functions to calculate
enviromnetal ranges and habitat suitability. The initial functions
download all species records from the Atlas and living Australia
(<https://www.ala.org.au/>) and the Global Biodiversity Information
Facility (GBIF, <https://www.gbif.org/>). The species data are
downloaded as individual .Rdata files to the specified folders, which
must exist first, without returning anything. The functions are
separated because the ALA and GBIF columns are slightly different, but
both data sources are needed to properly quantify species ranges.

  
  
  

The package functions expect these folders (a typical R project
structure), create them if they don’t exist

  

  
  
  

Now download GBIF and ALA occurrence data for each species

  
  
  

  
  
  

# STEP 2 :: Combine species occurrence data

  
  
  

The next function combines ALA and GBIF records, filtering them to
records on land, and recorded after 1950. The climate (i.e. raster) data
used can be any worldclim layer.

  
  
  

First, get some global climate data from worldclim. You need to create a
‘data’ folder in the project directory.

  
  
  

  
  
  

Finer resolution (1km) current worldclim grids are available here
:<https://drive.google.com/open?id=1mQHVmYxSMw_cw1iGvfU9M7Pq6Kl6nz-C>

  
  
  

Then trim the occurrence records to those inside the raster boundaries
(i.e. species records in the ocean according to the Raster boundaries
will be excluded)

  
  
  

  
  
  

# STEP 3 :: extract environmental values

  
  
  

The next function in the workflow combines occurrence files from ALA and
GBIF into one table, and extracts environmental values. It assumes that
both files come from the combine\_ala\_records and
combine\_gbif\_records functions.

  
  
  

First create a template raster of 1km \* 1km cells using the downloaded
worlclim data. This raster is used to filter records to 1 per one 1km
cell. This raster needs to have the same extent and resolution of the
data used to analyse the species distributions. It should have a value
of 1 for land, and NA for the ocean. This takes ages in R…..

  
  
  

``` r
## Set the Molleweide projection
sp_epsg54009 <- '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0'


## Use gdal to create a template raster in Mollweide
template_raster_1km_mol <- gdalwarp("./data/wc2-5/bio1.bil",
                                    tempfile(fileext = '.bil'),
                                    t_srs = sp_epsg54009,
                                    output_Raster = TRUE,
                                    tr = c(1000, 1000),
                                    r = "near", dstnodata = '-9999')


## Use gdal WGS84 
template_raster_1km_WGS84 <- template_raster_1km_mol %>% 
  projectRaster(., crs = CRS("+init=epsg:4326"))


## Should be 1km*1km, It should have a value of 1 for land, and NA for the ocean
template_raster_1km_WGS84[template_raster_1km_WGS84 > 0] <- 1
template_raster_1km_WGS84[template_raster_1km_WGS84 < 0] <- 1
xres(template_raster_1km_WGS84);projection(template_raster_1km_WGS84)
```

  
  
  

A pre-prepared template raster is found on google drive :
<https://drive.google.com/open?id=1mQHVmYxSMw_cw1iGvfU9M7Pq6Kl6nz-C>

  
  
  

  
  
  

Note that the order of the raster names in ‘world\_raster’ must match
the order of names in env\_variables. In this case, it’s simply the
biolclim variables (i.e. bio1-bio19)

  
  
  

``` r
## Combine GBIF and ALA data, and extract environmental values
COMBO.RASTER.CONVERT = combine_records_extract(ala_df          = ALA.LAND,
                                               gbif_df         = GBIF.LAND,
                                               urban_df        = 'NONE',
                                               thin_records    = TRUE,
                                               template_raster = template_raster_1km_WGS84,
                                               world_raster    = worldclim_climate,
                                               prj             = CRS("+init=epsg:4326"),
                                               species_list    = analysis_spp,
                                               biocl_vars      = bioclim_variables,
                                               env_vars        = env_variables,
                                               worldclim_grids = TRUE,
                                               save_data       = FALSE,
                                               save_run        = "TEST_BATS")
```

  
  
  

# STEP 4 :: Flag institutional outliers

  
  
  

The next stage of the process runs a series of cleaning steps. This
function takes a data frame of all species records, and flag records as
institutional or spatial outliers. This function uses the
CoordinateCleaner package
<https://cran.r-project.org/web/packages/CoordinateCleaner/index.html>.
It assumes that the records data.frame is that returned by the
combine\_records\_extract function.

  
  
  

``` r
## :: Flag records as institutional or spatial outliers
COORD.CLEAN = coord_clean_records(records    = COMBO.RASTER.CONVERT,
                                  capitals   = 10000,  ## Remove records within 10km  of capitals
                                  centroids  = 5000,   ## Remove records within 5km of country centroids
                                  save_run   = "TEST_SPECIES",
                                  save_data  = FALSE)
```

  
  
  

This function takes a data frame of all species records, flags records
as spatial outliers (T/F for each record in the df), and saves images of
the checks for each. Manual cleaning of spatial outliers is very
tedious, but automated cleaning makes mistakes, so checking is handy.
This funct uses the CoordinateCleaner package
<https://cran.r-project.org/web/packages/CoordinateCleaner/index.html>.
It assumes that the input dfs are those returned by the
coord\_clean\_records function.

  
  
  

``` r
## Step 4b :: Flag spatial outliers
## This is working, but it's just not plotting the ones that are being removed
SPATIAL.CLEAN = check_spatial_outliers(all_df       = COORD.CLEAN,
                                       land_shp     = LAND,
                                       urban_df     = FALSE, 
                                       clean_path   = './data/GBIF/Check_plots/',
                                       spatial_mult = 10,
                                       prj          = CRS("+init=epsg:4326"))
```

  
  
  

This function takes a data frame of all species records, estimates
geographic and environmental ranges for each, and creates a table of
each. It uses the AOO.computing function in the ConR package :
<https://cran.r-project.org/web/packages/ConR/index.html> It assumes
that the input df is that returned by the check\_spatial\_outliers
function.

  
  
  

``` r
## Step 4c ::Estimate climate niches usign species records
GLOB.NICHE = calc_1km_niches(coord_df     = SPATIAL.CLEAN,
                             prj          = CRS("+init=epsg:4326"),
                             country_shp  = AUS,
                             world_shp    = LAND,
                             kop_shp      = Koppen_shp,
                             species_list = analysis_spp,
                             env_vars     = env_variables,
                             cell_size    = 2,
                             save_run     = "Stoten_EG",
                             data_path    = "./output/results/",
                             save_data    = TRUE)
```

  
  
  

We can also plot the environmental ranges of each species. This function
takes a data frame of all species records, and plots histograms and
convex hulls for each species in global enviromental space It assumes
that the input df is that prepared by the check\_spatial\_outliers
function

  
  
  

``` r
## Step 4d :: plot species ranges using histograms and convex hulls for rainfall and temperature distributions
plot_range_histograms(coord_df     = SPATIAL.CLEAN,
                      species_list = analysis_spp,
                      range_path   = check_dir)
```

  
  
  

# STEP 5 :: Prepare SDM table

  
  
  

The final step before modelling is to create at table we can use for
species distribution modelling. This function takes a data frame of all
species records, and prepares a table in the ‘species with data’ (swd)
format for modelling uses the Maxent algorithm. It assumes that the
input df is that returned by the coord\_clean\_records function. There
is a switch in the function, that adds additional bakground points from
other taxa, if specified. In this example for bats, we’ll just use the
species supplied

  
  
  

``` r
## The final dataset is a spatial points dataframe in the Mollwiede projection
SDM.SPAT.OCC.BG = prepare_sdm_table(coord_df        = COORD.CLEAN,
                                    species_list    = unique(COORD.CLEAN$searchTaxon),
                                    sdm_table_vars  = sdm_table_vars,
                                    save_run        = "Stoten_EG",
                                    read_background = FALSE,
                                    save_data       = FALSE,
                                    save_shp        = FALSE)
```

  
  
  

# STEP 6 :: Run Global SDMs

  
  
  

This function takes a data frame of all species records, and runs a
specialised maxent analysis for each species. It uses the rmaxent
package <https://github.com/johnbaums/rmaxent> It assumes that the input
df is that returned by the prepare\_sdm\_table function

  
  
  

It uses a rasterised version of the 1975 Koppen raster, and another
template raster of the same extent (global), resolution (1km\*1km) and
projection (mollweide) as the analysis data. This step takes ages…

  
  
  

  
  
  

A pre-prepared template raster in the Mollweide projection is found on
google drive :
<https://drive.google.com/open?id=1mQHVmYxSMw_cw1iGvfU9M7Pq6Kl6nz-C>

  
  
  

  
  
  

The sdm function runs two maxent models: a full model using all
variables, and backwards selection. Given a candidate set of predictor
variables, the function identifies a subset of variables that meets
specified multicollinearity criteria. Subsequently, backward stepwise
variable selection is used to iteratively drop the variable that
contributes least to the model, until the contribution of each variable
meets a specified minimum, or until a predetermined minimum number of
predictors remains.

  
  
  

  
  
  

# STEP 7 :: Project SDMs across Australia

  
  
  

First, extract the SDM results from the models. We need the model
results to run the maps. Each model generates a ‘threshold’ of
probability of occurrence (see), which we use to create map of habitat
suitability across Australia ().

  
  
  

  
  
  

Now we need some future climate projections. We can download raster
worldclim data using the raster package. The stoten publication uses
climate projections under six global circulation models :

  
  
  

  
  
  

Or, we can load in some 1km\*1km Worldclim rasters for current
environmental conditions :

  
  
  

``` r
## 1km*1km Worldclim rasters for current environmental conditions can be found here:
## https://drive.google.com/open?id=1B14Jdv_NK2iWnmsqCxlKcEXkKIqwmRL_
aus.grids.current <- stack(
  file.path('./data/worldclim/aus/current', 
            sprintf('bio_%02d.tif', 1:19)))
```

  
  
  

The projection function takes the maxent models created by the
‘fit\_maxent\_targ\_bg\_back\_sel’ function, and projects the models
across geographic space - currently just for Australia. It uses the
rmaxent package <https://github.com/johnbaums/rmaxent>. It assumes that
the maxent models were generated by the
’fit\_maxent\_targ\_bg\_back\_sel’function. Note that this step is
quite memory heavy, and is best run with 32GB of RAM.

  
  
  

``` r
## Create a local projection for mapping : Australian Albers
aus_albers  <- CRS('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 
                   +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')

## Create 2070 sdm map projections
tryCatch(
  project_maxent_grids_mess(country_shp   = AUS, 
                            world_shp     = LAND,   
                            country_prj   = CRS("+init=EPSG:3577"),
                            world_prj     = CRS("+init=epsg:4326"),
                            local_prj     = aus_albers,
                            
                            scen_list     = scen_2070, 
                            species_list  = map_spp,    
                            maxent_path   = './output/maxent/back_sel_models/',
                            climate_path  = './data/worldclim/aus/',               
                            
                            grid_names    = env_variables,             
                            time_slice    = 70,                       
                            current_grids = aus.grids.current,         
                            create_mess   = TRUE,
                            OSGeo_path    = 'C:/OSGeo4W64/OSGeo4W.bat', 
                            nclust        = 1),
  
  ## If the species fails, write a fail message to file.
  error = function(cond) {
    
    ## This will write the error message inside the text file, but it won't include the species
    file.create(file.path("output/maxent/back_sel_models/mapping_failed_2070.txt"))
    cat(cond$message, file = file.path("output/maxent/back_sel_models/mapping_failed_2070.txt"))
    warning(cond$message)
    
  })
```

  
  
  

![](output/maxent/back_sel_models/Acacia_dealbata/full/Acacia_dealbata_mess_panel.png)

**Figure 2.** Example of a continuous climatic suitability map for one
plant species under current conditions. Species occurrence points are
plotted in red on the left panel. The cells in the right panel are coded
from 0 : no to low suitability, to 1 : highly suitable. The shaded areas
on the right panel indicate where the maxent model is extrapolating
beyond the training data (i.e. the result of a MESS map).

  
  
  

# STEP 8 :: Aggregate SDM projections within Spatial units

  
  
  

Now that all the species models and projections have been run, we need
to aggregate them across all six global circulation models. In order to
aggregate the results, we need a shapefile to aggregate to. In this
example, we’ll use the Australian Significant Urban Areas, which were
used in the Stoten article. A geo-tif of the Significant Areas is on
Google drive :

  
  
  

``` r
## Can't use the .rda, must use the file path
areal_unit_vec <- shapefile_vector_from_raster(shp_file = SUA,
                                               prj      = CRS("+init=EPSG:3577"),
                                               agg_var  = 'SUA_CODE16',
                                               temp_ras = aus.grids.current[[1]],
                                               targ_ras = './data/SUA_2016_AUST.tif')

## This is a vector of all the cells that either are or aren't in the rasterized shapefile
summary(areal_unit_vec)
```

  
  
  

This aggregation function uses the 10th% Logistic threshold for each
species from the maxent models to threhsold the rasters of habitat
suitability (0-1) For each GCM. For each species, summ the 6 GCMS to
create a binary raster with cell values between 0-6. These cell values
represent the number of GCMs where that cell had a suitability value
above the threshold determined by maxent. We classify a cell has
suitable if it met the threshold in \> 4 GCMs, and use this combined
raster to compare current and future suitability, measuring if the
suitability of each cell is changing over time, remaining stable or was
never suitable It assumes that the maxent predictions were generated by
the ‘project\_maxent\_grids\_mess’ function. Note that this step is
quite memory heavy, and is best run with 32GB of RAM.

  
  
  

``` r
## Combine GCM predictions and calculate gain and loss for 2030 
## Then loop over the species folders and climate scenarios
tryCatch(mapply(sdm_area_cell_count,                      
                unit_shp      = './data/SUA_albers.rds',  ## This would have to change
                unit_vec      = areal_unit_vec, 
                sort_var      = "SUA_NAME16",
                agg_var       = "SUA_CODE16",
                world_shp     = './data/LAND_albers.rds', ## This would have to change
                country_shp   = './data/AUS_albers.rds',  ## This would have to change
                
                DIR_list      = sdm.results.dir,  
                species_list  = map_spp,
                number_gcms   = 6,
                maxent_path   = 'output/maxent/back_sel_models/', 
                thresholds    = percent.10.log,
                time_slice    = 30,                     
                write_rasters = TRUE),
         
         ## If the species fails, write a fail message to file.
         error = function(cond) {
           
           ## This will write the error message inside the text file,
           ## but it won't include the species
           file.create(file.path("output/maxent/back_sel_models/sua_count_failed_2030.txt"))
           cat(cond$message, file=file.path("output/maxent/back_sel_models/sua_count_failed_2030.txt"))
           warning(cond$message)
           
         })
```

  
  
  

![](output/maxent/back_sel_models/Acacia_dealbata/full/Acacia_dealbata_gain_loss_0.3799_2030.png)

**Figure 3.** Example of a combined map of change in climatic
suitability from current conditions to 2070. Species occurrence points
are plotted in red on the left panel. The cells in the right and bottom
panels are coded as either lost (orange cells - present now but not in
2070 according to 4 or more GCMs), gained (green cells - absent now, but
present in 2070), stable (blue cells - present now and in 2070), or
never suitable (white cells - never present).
