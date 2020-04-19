######################################## READ SPATIAL DATA FOR SDM MODELLING ############################################



## This code reads in the vector and raster data that is used to create species distribution models
## Shapefiles don't need locations.
## Rasters store the location as part of the file, so the actual data needs to be stored there too...
## Need to check that the updated objects have the right file location in the slot
## This can only be tested in another location :: ask Shawn to have a go


## Create a function to get the spatial data ::
#https://www.rdocumentation.org/packages/raster/versions/3.0-12/topics/getData





# str(world.grids.current[[1]]@file@name)
######################################## READ SPATIAL DATA FOR SDM MODELLING ############################################



## This code reads in the vector and raster data that is used to create species distribution models
## Shapefiles don't need locations.
## Rasters store the location as part of the file, so the actual data needs to be stored there too...
## Need to check that the updated objects have the right file location in the slot


## Change these slots to be the project paths ::
# template.raster.1km.84@file@name
# str(world.grids.current[[1]]@file@name)




## 1). READ IN SHAPEFILES ===========================================================================================


## Set coordinate system definitions
# CRS.MOL      <- CRS('+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')
# CRS.MOL.SDM  <- CRS('+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0')
# CRS.WGS.84   <- CRS("+init=epsg:4326")
# CRS.AUS.ALB  <- CRS("+init=EPSG:3577")
# ALB.CONICAL  <- CRS('+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
# sp_epsg54009 <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0"


## Need to create a function that genererates all the spatial data


## Check how the Koppen zones were calculated
## This Needs to be made transparent - need to create functions for making the koppen layers too
Koppen_zones     = unique(readOGR('data/world_koppen/WC05_1975H_Koppen_Kriticos_2012.shp')@data[, 1:2])
Koppen_shp       = readOGR('data/world_koppen/WC05_1975H_Koppen_Kriticos_2012.shp')
Koppen_1975_1km  = raster('data/world_koppen/Koppen_1000m_Mollweide54009.tif')
AUS              = readRDS("./data/AUS_albers.rds")
LAND             = readRDS("./data/LAND_albers.rds")
#IBRA             = load('./data/IBRA.rda')
SUA              = readRDS('./data/AUS_albers.rds')


# AUS_albers <- AUS %>%
#   spTransform(ALB.CONICAL)
#
# LAND_albers <- LAND %>%
#   spTransform(ALB.CONICAL)
#
# SUA_albers <- SUA %>%
#   spTransform(ALB.CONICAL)
#
#
# saveRDS(LAND_albers, 'data/LAND_albers.rds')
# saveRDS(SUA_albers,  'data/SUA_albers.rds')
# saveRDS(AUS_albers,  'data/AUS_albers.rds')


## Read in shapefiles like this
# IBRA       = readOGR('data/base/Contextual/IBRA7_subregions_states.shp')
# saveRDS(IBRA, 'data/base/Contextual/IBRA7_SUB.rds')
# SUA       = readOGR('data/base/Contextual/SUA_2016_AUST.shp')
# saveRDS(SUA, 'data/SUA_2016_AUST.rds')


length(unique(IBRA$REG_NAME_7))
length(unique(SUA$SUA_CODE16))





## 2). READ IN BACKGROUND POINTS ===========================================================================================


## If needed, read in background data
#BG.POINTS     = readRDS("./data/base/background/SDM_SPAT_ALL_ANIMAL_BG_POINTS.rds")
TI.XY   = readRDS("./data/Global_Inventory_CLEAN_21.03.2019.rds")
stoten.spp = read.csv('./data/MAXENT_SUMMARY_SUA_ANALYSIS_NATIVE_GOOD_1502_2019.csv') %>%
  .$searchTaxon %>% as.character()


## 3). CREATE RASTER COMPONENTS ===========================================================================================



## These rasters could change, but the names in the projections, etc, would also need to change
## Create a raster stack of current global environmental conditions
## Turn this into a function
world.grids.current = stack(
  file.path('./data/worldclim/world/current',
            sprintf('bio_%02d', 1:19)))


## Create a raster stack of current Australian environmental conditions, and divide the current environmental grids by 10
aus.grids.current <- stack(
  file.path('./data/worldclim/aus/current',   ## ./green_cities_sdm/data/base/worldclim/aus/1km/bio
            sprintf('bio_%02d.tif', 1:19)))

for(i in 1:11) {


  ## simple loop
  message(i)
  aus.grids.current[[i]] <- aus.grids.current[[i]]/10

}





## 4). CREATE LISTS OF NAMES FOR THE GLOBAL CIRCULATION MODELS ----


## Use GDAL to create a raster which = 1 where bio_01 has data (i.e. land), and NA where there is no data
## Also note that gdalwarp is much faster, and the trs ='+init=esri:54009' argument does not work here


## Created like this :: 1km
# template.raster.1km <- gdalwarp("data/base/worldclim/world/0.5/bio/current/bio_01",
#                                 tempfile(fileext = '.tif'),
#                                 t_srs = '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs',
#                                 output_Raster = TRUE,
#                                 tr = c(1000, 1000),
#                                 r = "near", dstnodata = '-9999')
# xres(template.raster.1km)


## Load template rasters
template.raster.1km    = raster("./data/world_koppen/template_has_data_1km.tif")
template.raster.1km.84 = raster("./data/world_koppen/template_1km_WGS84.tif")
xres(template.raster.1km)


## this website has changed :: https://www.worldclim.org/data/index.html
## Create the names for the GCMs by scraping the worldclim website
# h <- read_html('http://www.worldclim.org/cmip5_30s')
# gcms <- h %>%
#   html_node('table') %>%
#   html_table(header = TRUE) %>%
#   filter(rcp85 != '')
#
# id.50 <- h %>%
#   html_nodes(xpath = "//a[text()='bi']/@href") %>%
#   as.character %>%
#   grep('85bi50', ., value = TRUE) %>%
#   basename %>%
#   sub('\\.zip.', '', .)
#
# id.70 <- h %>%
#   html_nodes(xpath = "//a[text()='bi']/@href") %>%
#   as.character %>%
#   grep('85bi70', ., value = TRUE) %>%
#   basename %>%
#   sub('\\.zip.', '', .)
#
# h <- read_html('http://www.worldclim.org/cmip5_30s')
# gcms <- h %>%
#   html_node('table') %>%
#   html_table(header = TRUE) %>%
#   filter(rcp85 != '')
#
# id.50 <- h %>%
#   html_nodes(xpath = "//a[text()='bi']/@href") %>%
#   as.character %>%
#   grep('85bi50', ., value = TRUE) %>%
#   basename %>%
#   sub('\\.zip.', '', .)
#
# id.70 <- h %>%
#   html_nodes(xpath = "//a[text()='bi']/@href") %>%
#   as.character %>%
#   grep('85bi70', ., value = TRUE) %>%
#   basename %>%
#   sub('\\.zip.', '', .)


## Create the IDs : A work around because the 2030 data comes from the climate change in Australia website, not worldclim
#id.30 = gsub("50", "30", id.50)


## Create the IDs
# gcms.30 <- cbind(gcms, id.30)
# gcms.30$GCM = sub(" \\(#\\)", "", gcms$GCM)
#
# gcms.50 <- cbind(gcms, id.50)
# gcms.50$GCM = sub(" \\(#\\)", "", gcms$GCM)  ## sub replaces first instance in a string, gsub = global
#
# gcms.70 <- cbind(gcms, id.70)
# gcms.70$GCM = sub(" \\(#\\)", "", gcms$GCM)
# gcms.70$GCM = sub(" \\(#\\)", "", gcms$GCM)


## Now create the scenario lists across which to loop
gcms.50 ; gcms.70 ; gcms.30


## Just get the 6 models picked by CSIRO for Australia, for 2030, 2050 and 2070
## See the publication for why we choose this
scen_2030 = c("mc85bi30", "no85bi30", "ac85bi30", "cc85bi30", "gf85bi30", "hg85bi30")
scen_2050 = c("mc85bi50", "no85bi50", "ac85bi50", "cc85bi50", "gf85bi50", "hg85bi50")
scen_2070 = c("mc85bi70", "no85bi70", "ac85bi70", "cc85bi70", "gf85bi70", "hg85bi70")


## Make a list of SDM columns needed ----
results.columns = c("searchTaxon",        ## From the ALA/ GBIF download code
                    "Origin",             ## native/extoic : from Anthony Manea's spreadsheet, affected by taxonomy....
                    "Family",             ## From Anthony Manea's spreadsheet, will be affected by taxonomy....
                    "Maxent_records",     ## No. records used in the SDM
                    "Aus_records",        ## No. AUS records     :: from the R workflow
                    "AOO",                ## Global Area of occurrence
                    "KOP_count",          ## Number of koppen zones each species is found in...D

                    "Number_var",        ## No. maxent variables :: from Maxent code
                    "Var_pcont",         ## Maxent Variable with highest permutation importance
                    "Per_cont",          ## The permutaiton importance of that variable
                    "Var_pimp",          ## Maxent Variable with highest permutation importance
                    "Perm_imp",          ## The permutaiton importance of that variable
                    "Iterations",               ## No. iterations

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


##
#save.image('BAT_SDM_DATA.RData')





## 4). CREATE LIST OF OBJECTS ===========================================================================================


## Create a list of spatial objects, not functions, etc.
object.list <- c(ls()[sapply(ls(), function(i) class(get(i))) == "data.frame"],
                 ls()[sapply(ls(), function(i) class(get(i))) == "CRS"],
                 ls()[sapply(ls(), function(i) class(get(i))) == "SpatialPolygonsDataFrame"],
                 ls()[sapply(ls(), function(i) class(get(i))) == "RasterLayer"],
                 ls()[sapply(ls(), function(i) class(get(i))) == "RasterStack"])


## Clean up the creation of objects - reading in one big RData file is probably not the way to go.......................

## Find small directories (<50 MB) and delete
# cd E:/bat_sdm_repo/output/maxent/full_models
# find -mindepth 1 -maxdepth 1 -type d -exec du -ks {} + | awk '$1 <= 50' | cut -f 2- | xargs -d \\n rm -rf





#####################################################  TBC ##############################################################
