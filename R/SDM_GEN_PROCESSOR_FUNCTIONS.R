#########################################################################################################################
######################################  FUNCTIONS FOR SDM ANALYSIS ---- #################################################
#########################################################################################################################


## Below is a list of the functions used to prepare the data for SDM analysis, and run the analysis


## Get a complete df ----
#' @export
completeFun <- function(data, desiredCols) {

  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])

}





## Estiamte the niche from species records ----
#' @export
niche_estimate = function (DF,
                           colname) {

  ## R doesn't seem to have a built-in mode function
  ## This doesn't really handel multiple modes, but it doesn't matter because mode was just calculated for Renee........
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }

  ## Use ddply inside a function to create niche widths and medians for each species
  ## This syntax is tricky, maybe ask John and Stu what they think

  ## Also, need to figure out how to make the aggregating column generic (species, genus, etc.)
  summary = ddply(DF,
                  .(searchTaxon),           ## currently grouping column only works hard-wired
                  .fun = function (xx, col) {

                    ## All the different columns
                    min      = min(xx[[col]])
                    max      = max(xx[[col]])

                    q02      = quantile(xx[[col]], .02)
                    q05      = quantile(xx[[col]], .05)
                    q95      = quantile(xx[[col]], .95)
                    q98      = quantile(xx[[col]], .98)

                    median   = median(xx[[col]])
                    mean     = mean(xx[[col]])
                    mode     = Mode(xx[[col]])
                    range    = max - min
                    q95_q05  = (q95 - q05)
                    q98_q02  = (q98 - q02)

                    ## Then crunch them together
                    c(min, max, median, mode, mean, range, q05, q95,  q95_q05, q98_q02)

                  },

                  colname

  )

  ## Concatenate output
  ## Also, need to figure out how to make the aggregating column generic (species, genus, etc.)
  ## currently it only works hard-wired
  colnames(summary) = c("searchTaxon",
                        paste0(colname,  "_min"),
                        paste0(colname,  "_max"),
                        paste0(colname,  "_median"),
                        paste0(colname,  "_mode"),
                        paste0(colname,  "_mean"),
                        paste0(colname,  "_range"),
                        paste0(colname,  "_q05"),
                        paste0(colname,  "_q95"),
                        paste0(colname,  "_q95_q05"),
                        paste0(colname,  "_q98_q02"))

  ## return the summary of niche width and median
  return (summary)

}





## GBIF download ----
#' @export
download_GBIF_all_species = function (species_list, path) {

  ## create variables
  skip.spp.list       = list()
  GBIF.download.limit = 200000

  ## for every species in the list
  for(sp.n in species_list){

    ## 1). First, check if the f*&%$*# file exists
    ## data\base\HIA_LIST\GBIF\SPECIES
    file_name = paste0(path, sp.n, "_GBIF_records.RData")

    ## If it's already downloaded, skip
    if (file.exists (file_name)) {

      print(paste ("file exists for species", sp.n, "skipping"))
      next

    }
    #  create a dummy file
    dummy = data.frame()
    save (dummy, file = file_name)

    ## 2). Then check the spelling...incorrect nomenclature will return NULL result
    if (is.null(occ_search(scientificName = sp.n, limit = 1)$meta$count) == TRUE) {

      ## now append the species which had incorrect nomenclature to the skipped list
      ## this is slow, but it works for now
      print (paste ("Possible incorrect nomenclature", sp.n, "skipping"))
      nomenclature = paste ("Possible incorrect nomenclature |", sp.n)
      skip.spp.list <- c(skip.spp.list, nomenclature)
      next

    }

    ## 3). Skip species with no records
    if (occ_search(scientificName = sp.n)$meta$count <= 2) {

      ## now append the species which had no records to the skipped list
      print (paste ("No GBIF records for", sp.n, "skipping"))
      records = paste ("No GBIF records |", sp.n)
      skip.spp.list <- c(skip.spp.list, records)
      next

    }

    ## 4). Check how many records there are, and skip if there are over 200k
    if (occ_search(scientificName = sp.n, limit = 1)$meta$count > GBIF.download.limit) {

      ## now append the species which had > 200k records to the skipped list
      print (paste ("Number of records > max for GBIF download via R (200,000)", sp.n))
      max =  paste ("Number of records > 200,000 |", sp.n)

    } else {

      ## 5). Download ALL records from GBIF
      message("Downloading GBIF records for ", sp.n, " using rgbif :: occ_data")
      key <- name_backbone(name = sp.n, rank = 'species')$usageKey

      GBIF <- occ_data(taxonKey = key, limit = GBIF.download.limit)
      GBIF <- as.data.frame(GBIF$data)

      cat("Synonyms returned for :: ",  sp.n, unique(GBIF$scientificName), sep="\n")
      cat("Names returned for :: ", sp.n, unique(GBIF$name),               sep="\n")
      cat("Takonkeys returned for :: ", sp.n, unique(GBIF$taxonKey),       sep="\n")

      ## Could also only use the key searched, but that could knock out a lot of species
      message(dim(GBIF[1]), " Records returned for ", sp.n)

      ## 6). save records to .Rdata file, note that using .csv files seemed to cause problems...
      save(GBIF, file = file_name)

    }

  }

}





## ALA download ----
#' @export
download_ALA_all_species = function (species_list, path) {

  ## create variables
  skip.spp.list       = list()
  #ALA.download.limit  = 200000

  ## for every species in the list
  for(sp.n in species_list){

    ## Get the ID?
    lsid <- ALA4R::specieslist(sp.n)$taxonConceptLsid

    ## 1). First, check if the f*&%$*# file exists
    file_name = paste0(path, sp.n, "_ALA_records.RData")

    ## If it's already downloaded, skip
    if (file.exists (file_name)) {

      print (paste ("file exists for species", sp.n, "skipping"))
      next

    }
    #  create a dummy file
    dummy = data.frame()
    save (dummy, file = file_name)

    ## 2). Then check the spelling...incorrect nomenclature will return NULL result
    if (is.null(ALA4R::occurrences(taxon = paste('taxon_name:\"', sp.n, '\"',sep=""),
                                   download_reason_id = 7)$data) == TRUE) {

      ## Now, append the species which had incorrect nomenclature to the skipped list
      print (paste ("Possible incorrect nomenclature", sp.n, "skipping"))
      nomenclature = paste ("Possible incorrect nomenclature |", sp.n)
      skip.spp.list <- c(skip.spp.list, nomenclature)
      next

    }

    ## 3). Skip species with no records
    if (nrow(ALA4R::occurrences(taxon = paste('taxon_name:\"', sp.n, '\"',sep=""), download_reason_id = 7)$data) <= 2) {

      ## now append the species which had no records to the skipped list
      print (paste ("No ALA records for", sp.n, "skipping"))
      records = paste ("No ALA records |", sp.n)
      skip.spp.list <- c(skip.spp.list, records)
      next

    }

    ## Download ALL records from ALA ::
    message("Downloading ALA records for ", sp.n, " using ALA4R :: occurrences")
    ALA = ALA4R::occurrences(taxon = paste('taxon_name:\"', sp.n, '\"',sep=""), download_reason_id = 7)
    ALA = ALA[["data"]]

    cat("Synonyms returned for :: ", sp.n, unique(ALA$scientificName), sep="\n")
    message(dim(ALA[1]), " Records returned for ", sp.n)

    ## Save records to .Rdata file
    save(ALA, file = file_name)

  }

}





## Create vector of raster -----
#' @export
shapefile_vector_from_raster = function (shp_file,
                                         prj,
                                         sort_var,
                                         agg_var,
                                         temp_ras,
                                         targ_ras) {

  ## Read in shapefile
  areal_unit <- shp_file %>%
    spTransform(prj)

  ## Rasterize shapefile, insert 'sort_var'
  message('Rasterizing shapefile')
  areal_unit = areal_unit[order(areal_unit[[sort_var]]),]
  #areal_unit = areal_unit[order(areal_unit$SUA_NAME16),]
  f <- tempfile()

  ## Write a temporary raster
  writeOGR(areal_unit, tempdir(), basename(f), 'ESRI Shapefile')
  template <- raster(temp_ras)

  ## Rasterize the shapefile
  areal_unit_rast <- gdalUtils::gdal_rasterize(
    normalizePath(paste0(f, '.shp')),

    ## The missing step is the .tiff, it was created somewhere else, in R or a GIS
    ## so 'a' needs to match in this step, and the next
    targ_ras, tr = res(template),
    te = c(bbox(template)), a = agg_var, a_nodata = 0, init = 0, ot = 'UInt16', output_Raster = TRUE)

  areal_unit_vec <- c(areal_unit_rast[])
  summary(areal_unit_vec)

  ## return the vector
  return(areal_unit_vec)

}


## Combining ALA records -----
#' @export
combine_ala_records = function(species_list, records_path, records_extension, record_type, keep_cols, world_raster) {

  download = list.files(records_path, pattern = ".RData")
  length(download)

  ## Now these lists are getting too long for the combine step.
  ## Restrict them to just the strings that partially match the  species list for each run
  spp.download <- paste(species_list, records_extension, sep = "")
  download     = download[download %in% spp.download ]
  message('downloaded species ', length(download), ' analyzed species ', length(species_list))


  ## Combine all the taxa into a single dataframe at once
  ALL <- download %>%

    ## Pipe the list into lapply
    ## x = download [1]
    lapply(function(x) {

      ## Create a character string of each .RData file
      message("Reading records data for ", x)
      f <- sprintf(paste0(records_path, "%s"), x)

      ## Load each file - check if some are already dataframes
      d <- get(load(f))
      if (length(class(d)) > 1) {

        d <- d[["data"]]

      } else {

        d = d

      }

      ## Check if the dataframes have data
      if (nrow(d) <= 2) {

        ## If the species has < 2 records, escape the loop
        print (paste ("No occurrence records for ", x, " skipping "))
        return (d)

      }

      ##  type standardisation
      names(d)[names(d) == 'latitude']  <- 'lat'
      names(d)[names(d) == 'longitude'] <- 'lon'

      ##  standardi[sz]e catnum colname
      if("catalogueNumber" %in% colnames(d)) {
        message ("Renaming catalogueNumber column to catalogNumber")
        names(d)[names(d) == 'catalogueNumber'] <- 'catalogNumber'

      }

      if (!is.character(d$catalogNumber)) {
        d$catalogNumber = as.character(d$catalogNumber)

      }

      ## standardi[sz]e catnum colname
      if('coordinateUncertaintyinMetres' %in% colnames(d)) {
        message ("Renaming recordID column to id")
        names(d)[names(d) == 'coordinateUncertaintyinMetres'] <- 'coordinateUncertaintyInMetres'

      }

      ## standardi[sz]e catnum colname
      if('recordID' %in% colnames(d)) {
        message ("Renaming recordID column to id")
        names(d)[names(d) == 'recordID'] <- 'id'

      }

      ## Create the searchTaxon column
      message ('Formatting occurrence data for ', x)
      d[,"searchTaxon"] = x
      d[,"searchTaxon"] = gsub(records_extension, "", d[,"searchTaxon"])

      if(!is.character(d["id"])) {
        d["id"] <- as.character(d["id"])
      }

      ## Choose only the desired columns
      d = d %>%
        select(one_of(keep_cols))

      ## Then print warnings
      warnings()

      ## This is a list of columns in different ALA files which have weird characters
      message ('Formatting numeric occurrence data for ', x)
      d[,"coordinateUncertaintyInMetres"] = as.numeric(unlist(d["coordinateUncertaintyInMetres"]))
      d["year"]  = as.numeric(unlist(d["year"]))
      d["month"] = as.numeric(unlist(d["month"]))
      d["id"]    = as.character(unlist(d["id"]))

      return(d)

    }) %>%

    ## Finally, bind all the rows together
    bind_rows

  ## Clear the garbage
  gc()

  ## Just get the newly downloaded species
  ALL = ALL[ALL$searchTaxon %in% species_list, ]
  length(unique(ALL$searchTaxon))

  if (nrow(ALL) > 0) {

    ## What names get returned?
    sort(names(ALL))
    TRIM <- ALL%>%
      dplyr::select(dplyr::one_of(keep_cols))

    dim(TRIM)
    sort(names(TRIM))


    ## What are the unique species?
    (sum(is.na(TRIM$scientificName)) + nrow(subset(TRIM, scientificName == "")))/nrow(TRIM)*100


    ## 3). FILTER RECORDS TO THOSE WITH COORDINATES, AND AFTER 1950
    ## Now filter the ALA records using conditions which are not too restrictive
    CLEAN <- TRIM %>%

      ## Note that these filters are very forgiving...
      ## Unless we include the NAs, very few records are returned!
      filter(!is.na(lon) & !is.na(lat),
             year >= 1950 & !is.na(year))

    ## How many records were removed by filtering?
    message(nrow(TRIM) - nrow(CLEAN), " records removed")
    message(round((nrow(CLEAN))/nrow(TRIM)*100, 2),
            " % records retained using spatially valid records")

    ## Can use WORLDCIM rasters to get only records where wordlclim data is.
    message('Removing ALA points outside raster bounds for ', length(species_list),
            ' species')

    ## Now get the XY centroids of the unique 1km * 1km WORLDCLIM blocks where ALA records are found
    ## Get cell number(s) of WORLDCLIM raster from row and/or column numbers. Cell numbers start at 1 in the upper left corner,
    ## and increase from left to right, and then from top to bottom. The last cell number equals the number of raster cells
    xy <- cellFromXY(world_raster, CLEAN[c("lon", "lat")]) %>%

      ## get the unique raster cells
      unique %>%

      ## Get coordinates of the center of raster cells for a row, column, or cell number of WORLDCLIM raster
      xyFromCell(world_raster, .) %>%
      na.omit()

    ## For some reason, we need to convert the xy coords to a spatial points data frame, in order to avoid this error:
    ## 'NAs introduced by coercion to integer range'
    xy <- SpatialPointsDataFrame(coords = xy, data = as.data.frame(xy),
                                 proj4string = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))

    ## Now extract the temperature values for the unique 1km centroids which contain ALA data
    class(xy)
    z   = raster::extract(world_raster[["bio_01"]], xy)

    ## Then track which values of Z are on land or not
    onland = z %>% is.na %>%  `!` # %>% xy[.,]  cells on land or not

    ## Finally, filter the cleaned ALA data to only those points on land.
    ## This is achieved with the final [onland]
    LAND.POINTS = filter(CLEAN, cellFromXY(world_raster, CLEAN[c("lon", "lat")]) %in%
                           unique(cellFromXY(world_raster, CLEAN[c("lon", "lat")]))[onland])

    ## how many records were on land?
    records.ocean = dim(CLEAN)[1] - dim(LAND.POINTS)[1]  ## 91575 records are in the ocean

    ## Print the dataframe dimensions to screen
    dim(LAND.POINTS)
    length(unique(LAND.POINTS$searchTaxon))

    ## Add a source column
    LAND.POINTS$SOURCE = record_type
    message(round((dim(LAND.POINTS)[1])/dim(CLEAN)[1]*100, 2),
            " % records retained using spatially valid records")

    ## save data
    dim(LAND.POINTS)
    length(unique(LAND.POINTS$searchTaxon))

    ## get rid of some memory
    gc()

  } else {

    message('No ALA dat for this set of taxa, creating empty datframe to other data')
    LAND.POINTS  = setNames(data.frame(matrix(ncol = length(keep), nrow = 0)), keep)

  }

  return(LAND.POINTS)

}





## Combining GBIF records -----
#' @export
combine_gbif_records = function(species_list, records_path, records_extension, record_type, keep_cols, world_raster) {

  download = list.files(records_path, pattern = ".RData")
  length(download)

  ## Now these lists are getting too long for the combine step.
  ## Restrict them to just the strings that partially match the  species list for each run
  spp.download <- paste(species_list, records_extension, sep = "")
  download     = download[download %in% spp.download ]
  message('downloaded species ', length(download), ' analyzed species ', length(species_list))

  ALL.POINTS <- download %>%

    ## Pipe the list into lapply
    lapply(function(x) {

      ## x = download[1]
      ## Create a character string of each .RData file
      f <- sprintf(paste0(records_path, "%s"), x)

      ## Load each file
      d <- get(load(f))

      ## Now drop the columns which we don't need
      message ('Reading GBIF data for ', x)

      ## Check if the dataframes have data
      if (nrow(d) <= 2) {

        ## If the species has < 2 records, escape the loop
        print (paste ("No GBIF records for ", x, " skipping "))
        return (d)

      }

      dat <- data.frame(searchTaxon = x, d[, colnames(d) %in% keep_cols],
                        stringsAsFactors = FALSE)

      if(!is.character(dat$gbifID)) {

        dat$gbifID <- as.character(dat$gbifID)

      }

      ## Need to print the object within the loop
      names(dat)[names(dat) == 'decimalLatitude']  <- 'lat'
      names(dat)[names(dat) == 'decimalLongitude'] <- 'lon'
      dat$searchTaxon = gsub("_GBIF_records.RData", "", dat$searchTaxon)
      return(dat)

    }) %>%

    ## Finally, bind all the rows together
    bind_rows

  ## If there is GBIF data
  if (nrow(ALL.POINTS) > 0) {

    ## What proportion of the dataset has no lat/lon? Need to check this so we know the latest download is working
    formatC(dim(ALL.POINTS)[1], format = "e", digits = 2)
    (sum(is.na(ALL.POINTS$lat))            + dim(subset(ALL.POINTS, year < 1950))[1])/dim(ALL.POINTS)[1]*100

    ## Almost none of the GBIF data has no scientificName. This is the right field to use for matching taxonomy
    (sum(is.na(ALL.POINTS$scientificName)) + dim(subset(ALL.POINTS, scientificName == ""))[1])/dim(ALL.POINTS)[1]*100


    ## Now get just the columns we want to keep.
    GBIF.TRIM <- ALL.POINTS%>%
      dplyr::select(one_of(keep_cols))
    names(GBIF.TRIM)
    gc()

    ## Just get the newly downloaded species
    GBIF.TRIM = GBIF.TRIM[GBIF.TRIM$searchTaxon %in% species_list, ]
    formatC(dim(GBIF.TRIM)[1], format = "e", digits = 2)

    ## What are the unique species?
    length(unique(GBIF.TRIM$species))
    length(unique(GBIF.TRIM$searchTaxon))
    length(unique(GBIF.TRIM$scientificName))

    ## FILTER RECORDS TO THOSE WITH COORDINATES, AND AFTER 1950

    ## Now filter the GBIF records using conditions which are not too restrictive
    GBIF.CLEAN <- GBIF.TRIM %>%

      ## Note that these filters are very forgiving...
      ## Unless we include the NAs, very few records are returned!
      filter(!is.na(lon) & !is.na(lat),
             year >= 1950 & !is.na(year))

    ## How many species are there?
    names(GBIF.CLEAN)
    length(unique(GBIF.CLEAN$searchTaxon))

    ## REMOVE POINTS OUTSIDE WORLDCLIM LAYERS

    ## Can use WORLDCIM rasters to get only records where wordlclim data is.
    message('Removing GBIF points in the ocean for ', length(species_list), ' species')

    ## Now get the XY centroids of the unique 1km * 1km WORLDCLIM blocks where GBIF records are found
    ## Get cell number(s) of WORLDCLIM raster from row and/or column numbers. Cell numbers start at 1 in the upper left corner,
    ## and increase from left to right, and then from top to bottom. The last cell number equals the number of raster cells
    xy <- cellFromXY(world.grids.current, GBIF.CLEAN[c("lon", "lat")]) %>%

      ## get the unique raster cells
      unique %>%

      ## Get coordinates of the center of raster cells for a row, column, or cell number of WORLDCLIM raster
      xyFromCell(world.grids.current, .) %>%
      na.omit()

    ## For some reason, we need to convert the xy coords to a spatial points data frame, in order to avoid this error:
    ## 'NAs introduced by coercion to integer range'
    xy <- SpatialPointsDataFrame(coords = xy, data = as.data.frame(xy),
                                 proj4string = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))

    ## Now extract the temperature values for the unique 1km centroids which contain GBIF data
    message('Removing GBIF points in the ocean for ', length(species_list), ' species')
    class(xy)
    z   = raster::extract(world.grids.current[["bio_01"]], xy)

    ## Then track which values of Z are on land or not
    onland = z %>% is.na %>%  `!` # %>% xy[.,]  cells on land or not

    ## Finally, filter the cleaned GBIF data to only those points on land.
    ## This is achieved with the final [onland]
    LAND.POINTS = filter(GBIF.CLEAN, cellFromXY(world.grids.current, GBIF.CLEAN[c("lon", "lat")]) %in%
                           unique(cellFromXY(world.grids.current,    GBIF.CLEAN[c("lon", "lat")]))[onland])

    ## how many records were on land?
    records.ocean = nrow(GBIF.CLEAN) - nrow(LAND.POINTS)  ## 91575 records are in the ocean

    ## Print the dataframe dimensions to screen
    dim(LAND.POINTS)
    length(unique(LAND.POINTS$searchTaxon))

    ## Add a source column
    LAND.POINTS$SOURCE = 'GBIF'
    unique(LAND.POINTS$SOURCE)
    names(LAND.POINTS)

    ## Free some memory
    gc()

  } else {

    message('No GBIF data for this taxon, creating empty datframe to bind to GBIF data')
    LAND.POINTS  = setNames(data.frame(matrix(ncol = length(keep_cols), nrow = 0)), keep_cols)

  }

  return(LAND.POINTS)

}





## Extract environmental values for occurrence records -----
#' @export
combine_records_extract = function(ala_df,
                                   gbif_df,
                                   species_list,
                                   template_raster,
                                   world_raster,
                                   projection,
                                   biocl_vars,
                                   env_vars,
                                   worldclim_grids,
                                   save_data,
                                   save_run) {

  ## Get just the Common columns
  common.cols    = intersect(names(gbif_df), names(ala_df))
  GBIF.ALA.COMBO = bind_rows(gbif_df, ala_df)
  GBIF.ALA.COMBO = GBIF.ALA.COMBO %>%
    select(one_of(common.cols))

  message(length(unique(GBIF.ALA.COMBO$searchTaxon)))
  length(unique(GBIF.ALA.COMBO$scientificName))


  ## CHECK TAXONOMY RETURNED BY ALA USING TAXONSTAND
  GBIF.ALA.MATCH = GBIF.ALA.COMBO

  ## Create points: the 'over' function seems to need geographic coordinates for this data...
  GBIF.ALA.84   = SpatialPointsDataFrame(coords      = GBIF.ALA.MATCH[c("lon", "lat")],
                                         data        = GBIF.ALA.MATCH,
                                         proj4string = projection)

  ## The length needs to be the same
  length(unique(GBIF.ALA.84$searchTaxon))
  GBIF.ALA.84.SPLIT.ALL <- split(GBIF.ALA.84, GBIF.ALA.84$searchTaxon)
  occurrence_cells_all  <- lapply(GBIF.ALA.84.SPLIT.ALL, function(x) cellFromXY(template_raster, x))

  ## Check with a message, but could check with a fail
  message('Split prodcues ', length(occurrence_cells_all), ' data frames for ', length(species_list), ' species')

  ## Now get just one record within each 10*10km cell.
  GBIF.ALA.84.1KM <- mapply(function(x, cells) {
    x[!duplicated(cells), ]
  }, GBIF.ALA.84.SPLIT.ALL, occurrence_cells_all, SIMPLIFY = FALSE) %>% do.call(rbind, .)

  ## Check to see we have 19 variables + the species for the standard predictors, and 19 for all predictors
  message(round(nrow(GBIF.ALA.84.1KM)/nrow(GBIF.ALA.84)*100, 2), " % records retained at 1km resolution")

  ## Create points: the 'over' function seems to need geographic coordinates for this data...
  COMBO.POINTS   = GBIF.ALA.84.1KM[c("lon", "lat")]

  ## Bioclim variables
  ## Extract raster data
  message('Extracting raster values for ', length(species_list), ' species in the set ', "'", save_run, "'")
  message(projection(COMBO.POINTS));message(projection(world_raster))
  dim(COMBO.POINTS);dim(GBIF.ALA.84.1KM)

  ## Extract the raster values
  COMBO.RASTER <- raster::extract(world_raster, COMBO.POINTS) %>%
    cbind(as.data.frame(GBIF.ALA.84.1KM), .)

  ## Group rename the columns
  setnames(COMBO.RASTER, old = biocl_vars, new = env_vars)
  COMBO.RASTER <- COMBO.RASTER %>% select(-lat.1, -lon.1)

  ## Change the raster values here: See http://worldclim.org/formats1 for description of the interger conversion.
  ## All worldclim temperature variables were multiplied by 10, so then divide by 10 to reverse it.
  if (worldclim_grids == "TRUE") {

    ## Convert the worldclim grids
    message('Processing worldclim 1.0 data, divide the rasters by 10')

    COMBO.RASTER.CONVERT = as.data.table(COMBO.RASTER)
    COMBO.RASTER.CONVERT[, (env.variables [c(1:11)]) := lapply(.SD, function(x)
      x / 10 ), .SDcols = env.variables [c(1:11)]]
    COMBO.RASTER.CONVERT = as.data.frame(COMBO.RASTER.CONVERT)

  } else {

    message('Not rocessing worldclim data, dont divide the rasters')
    COMBO.RASTER.CONVERT = COMBO.RASTER
  }

  ## Print the dataframe dimensions to screen :: format to recognise millions, hundreds of thousands, etc.
  COMBO.RASTER.CONVERT = completeFun(COMBO.RASTER.CONVERT, env_vars[1])

  message(length(unique(COMBO.RASTER.CONVERT$searchTaxon)),
          ' species processed of ', length(species_list), ' original species')

  ## save data
  if(save_data == "TRUE") {

    ## save .rds file for the next session
    saveRDS(COMBO.RASTER.CONVERT, paste0(DATA_path, 'COMBO_RASTER_CONVERT_',  save_run, '.rds'))

  } else {

    return(COMBO.RASTER.CONVERT)

  }
  ## get rid of some memory
  gc()
}





## Clean coordinates of occurence records ----
#' @export
coord_clean_records = function(records,
                               capitals,
                               centroids,
                               save_run,
                               save_data) {

  ## Create a unique identifier. This is used for automated cleaing of the records, and also saving shapefiles
  ## But this will not be run for all species linearly. So, it probably needs to be a combination of species and number
  records$CC.OBS <- 1:nrow(records)
  records$CC.OBS <- paste0(records$CC.OBS, "_CC_", records$searchTaxon)
  records$CC.OBS <- gsub(" ",     "_",  records$CC.OBS, perl = TRUE)
  length(records$CC.OBS);length(unique(records$CC.OBS))
  records$species = records$searchTaxon


  ## Rename the columns to fit the CleanCoordinates format and create a tibble.
  TIB.GBIF <- records %>% dplyr::rename(coord_spp        = searchTaxon,
                                        decimallongitude = lon,
                                        decimallatitude  = lat) %>%

    ## Then create a tibble for running the spatial outlier cleaning
    timetk::tk_tbl() %>%

    ## Consider the arguments. We've already stripped out the records that fall outside
    ## the worldclim raster boundaries, so the sea test is probably not the most important
    ## Study area is the globe, but we are only projecting models onto Australia

    ## The geographic outlier detection is not working here. Try it in setp 6
    ## Also, the duplicates step is not working, flagging too many records
    clean_coordinates(.,
                      verbose         = TRUE,
                      tests = c("capitals",     "centroids", "equal", "gbif",
                                "institutions", "zeros"), ## duplicates flagged too many

                      ## remove records within 10km  of capitals
                      ## remove records within 5km of country centroids
                      capitals_rad  = capitals,
                      centroids_rad = centroids) %>%

    ## The select the relevant columns and rename
    select(., coord_spp, CC.OBS, .val,  .equ, .zer, .cap,
           .cen,    .gbf,   .inst, .summary)

  ## Then rename
  summary(TIB.GBIF)
  names(TIB.GBIF) = c("coord_spp", "CC.OBS",    "coord_val",  "coord_equ",  "coord_zer",  "coord_cap",
                      "coord_cen", "coord_gbf", "coord_inst", "coord_summary")

  ## Flagging ~ x%, excluding the spatial outliers. Seems reasonable?
  message(round(with(TIB.GBIF, table(coord_summary)/sum(table(coord_summary))*100), 2), " % records removed")

  ## Check the order still matches
  message(identical(records$CC.OBS, TIB.GBIF$CC.OBS), ' records order matches')
  message(identical(records$searchTaxon, TIB.GBIF$coord_spp))

  ## Is the species column the same as the searchTaxon column?
  COORD.CLEAN = join(records, TIB.GBIF)
  message(identical(records$searchTaxon, COORD.CLEAN$coord_spp), ' records order matches')
  identical(records$CC.OBS, COORD.CLEAN$CC.OBS)

  ## Now subset to records that are flagged as outliers
  message(table(COORD.CLEAN$coord_summary))

  ## What percentage of records are retained?
  message(length(unique(COORD.CLEAN$searchTaxon)))

  ## Remove duplicate coordinates
  drops <- c("lon.1", "lat.1")
  COORD.CLEAN <- COORD.CLEAN[ , !(names(COORD.CLEAN) %in% drops)]
  unique(COORD.CLEAN$SOURCE)

  CLEAN.TRUE <- subset(COORD.CLEAN, coord_summary == "TRUE")
  message(round(nrow(CLEAN.TRUE)/nrow(COORD.CLEAN)*100, 2), " % records retained")
  message(table(COORD.CLEAN$coord_summary))


  if(save_data == "TRUE") {

    ## save .rds file for the next session
    saveRDS(COORD.CLEAN, paste0(DATA_path, 'COORD_CLEAN_', save_run, '.rds'))

  } else {

    message('Return the cleaned occurrence data to the global environment')   ##
    return(COORD.CLEAN)

  }

}





## Save spatial outliers to file ----
#' @export
check_spatial_outliers = function(all_df,
                                  land_shp,
                                  clean_path,
                                  spatial_mult) {

  ## Try plotting the points which are outliers for a subset of spp and label them
  ALL.PLOT = SpatialPointsDataFrame(coords      = all_df[c("lon", "lat")],
                                    data        = all_df,
                                    proj4string = CRS.WGS.84)

  CLEAN.TRUE = subset(all_df, coord_summary == "TRUE")
  CLEAN.PLOT = SpatialPointsDataFrame(coords      = CLEAN.TRUE[c("lon", "lat")],
                                      data        = CLEAN.TRUE,
                                      proj4string = CRS.WGS.84)

  ## Create global and australian shapefile in the local coordinate system
  LAND.84 = land_shp %>%
    spTransform(CRS.WGS.84)
  AUS.84 = AUS %>%
    spTransform(CRS.WGS.84)

  ## spp = plot.taxa[1]
  plot.taxa <- as.character(unique(CLEAN.PLOT$searchTaxon))
  for (spp in plot.taxa) {

    ## Plot a subset of taxa
    CLEAN.PLOT.PI = CLEAN.PLOT[ which(CLEAN.PLOT$searchTaxon == spp), ]

    message("plotting occ data for ", spp, ", ",
            nrow(CLEAN.PLOT.PI), " records flagged as either ",
            unique(CLEAN.PLOT.PI$coord_summary))

    ## Plot true and false points for the world
    ## Black == FALSE
    ## Red   == TRUE
    message('Writing map of global coord clean records for ', spp)
    png(sprintf("%s%s_%s", clean_path, spp, "global_check_cc.png"),
        16, 10, units = 'in', res = 500)

    par(mfrow = c(1,2))
    plot(LAND.84, main = paste0(nrow(subset(CLEAN.PLOT.PI, coord_summary == "FALSE")),
                                " Global clean_coord 'FALSE' points for ", spp),
         lwd = 0.01, asp = 1, col = 'grey', bg = 'sky blue')

    points(CLEAN.PLOT.PI,
           pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2,
           xlab = "", ylab = "", asp = 1,
           col = factor(CLEAN.PLOT.PI$coord_summary))

    ## Plot true and false points for the world
    ## Black == FALSE
    ## Red   == TRUE
    plot(AUS.84, main = paste0(nrow(subset(CLEAN.PLOT.PI, coord_summary == "FALSE")),
                               " Global clean_coord 'FALSE' points for ", spp),
         lwd = 0.01, asp = 1, bg = 'sky blue', col = 'grey')

    points(CLEAN.PLOT.PI,
           pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2,
           xlab = "", ylab = "", asp = 1,
           col = factor(CLEAN.PLOT.PI$coord_summary))

    dev.off()

  }

  ## Create a tibble to supply to coordinate cleaner
  test.geo = SpatialPointsDataFrame(coords      = all_df[c("lon", "lat")],
                                    data        = all_df,
                                    proj4string = CRS.WGS.84)

  SDM.COORDS  <- test.geo %>%

    spTransform(., CRS.WGS.84) %>%

    as.data.frame() %>%

    select(searchTaxon, lon, lat, CC.OBS, SOURCE) %>%

    dplyr::rename(species          = searchTaxon,
                  decimallongitude = lon,
                  decimallatitude  = lat) %>%

    timetk::tk_tbl()

  ## Check
  dim(SDM.COORDS)
  head(SDM.COORDS)
  class(SDM.COORDS)
  summary(SDM.COORDS$decimallongitude)
  identical(SDM.COORDS$index, SDM.COORDS$CC.OBS)
  length(unique(SDM.COORDS$species))

  ## Check how many records each spp has...
  COMBO.LUT <- SDM.COORDS %>%
    as.data.frame() %>%
    select(species) %>%
    table() %>%
    as.data.frame()
  COMBO.LUT <- setDT(COMBO.LUT, keep.rownames = FALSE)[]
  names(COMBO.LUT) = c("species", "FREQUENCY")
  COMBO.LUT = COMBO.LUT[with(COMBO.LUT, rev(order(FREQUENCY))), ]

  ## Watch out here - this sorting could cause problems for the order of the data frame once it's stitched back together
  ## If we we use spp to join the data back together, will it preserve the order?
  LUT.100K = as.character(subset(COMBO.LUT, FREQUENCY < 100000)$species)
  LUT.100K = trimws(LUT.100K [order(LUT.100K)])
  length(LUT.100K)

  ## Create a data frame of spp name and spatial outlier
  SPAT.OUT <- LUT.100K  %>%

    ## pipe the list of spp into lapply
    lapply(function(x) {

      ## Create the spp df by subsetting by spp
      f <- subset(SDM.COORDS, species == x)

      ## Run the spatial outlier detection
      message("Running spatial outlier detection for ", x)
      message(dim(f)[1], " records for ", x)
      sp.flag <- cc_outl(f,
                         lon     = "decimallongitude",
                         lat     = "decimallatitude",
                         species = "species",
                         method  = "quantile",
                         mltpl   = spatial_mult,
                         value   = "flagged",
                         verbose = "TRUE")

      ## Now add attache column for spp, and the flag for each record
      d = cbind(searchTaxon = x,
                SPAT_OUT = sp.flag, f)[c("searchTaxon", "SPAT_OUT", "CC.OBS")]

      ## Remeber to explicitly return the df at the end of loop, so we can bind
      return(d)

    }) %>%

    ## Finally, bind all the rows together
    bind_rows

  gc()

  ## Join the data back on
  SPAT.FLAG = join(as.data.frame(test.geo), SPAT.OUT)    ## Join means the skipped spp are left out
  dim(SPAT.FLAG)

  ## Try plotting the points which are outliers for a subset of spp and label them
  SPAT.FLAG = SpatialPointsDataFrame(coords      = SPAT.FLAG[c("lon", "lat")],
                                     data        = SPAT.FLAG,
                                     proj4string = CRS.WGS.84)

  ## Get the first 10 spp
  ## spp = plot.taxa[1]
  plot.taxa <- as.character(unique(SPAT.FLAG$searchTaxon))
  for (spp in plot.taxa) {

    ## Plot a subset of taxa
    CLEAN.PLOT.PI   = subset(SPAT.FLAG, searchTaxon == spp)

    message("plotting occ data for ", spp, ", ",
            nrow(CLEAN.PLOT.PI ), " clean records")

    ## Plot true and false points for the world
    ## Black == FALSE
    ## Red   == TRUE
    message('Writing map of global coord clean records for ', spp)
    png(sprintf("%s%s_%s", clean_path, spp, "global_spatial_outlier_check.png"),
        16, 10, units = 'in', res = 500)

    par(mfrow = c(1,2))
    plot(LAND.84, main = paste0(nrow(subset(CLEAN.PLOT.PI, SPAT_OUT == "FALSE")),
                                " Spatial outlier 'FALSE' points for ", spp),
         lwd = 0.01, asp = 1, col = 'grey', bg = 'sky blue')

    points(CLEAN.PLOT.PI,
           pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2,
           xlab = "", ylab = "", asp = 1,
           col = factor(CLEAN.PLOT.PI$SPAT_OUT))

    ## Plot true and false points for the world
    ## Black == FALSE
    ## Red   == TRUE
    plot(AUS.84, main = paste0("Australian points for ", spp),
         lwd = 0.01, asp = 1, bg = 'sky blue', col = 'grey')

    points(CLEAN.PLOT.PI,
           pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2,
           xlab = "", ylab = "", asp = 1,
           col = factor(CLEAN.PLOT.PI$SPAT_OUT))

    dev.off()

  }

}




## Calculate niches using occurrence data ----
#' @export
calc_1km_niches = function(coord_df,
                           aus_shp,
                           world_shp,
                           kop_shp,
                           ibra_shp,
                           species_list,
                           env_vars,
                           cell_size,
                           save_run,
                           save_data) {

  ## Create a spatial points object
  message('Estimating global niches for ', length(species_list), ' species across ',
          length(env_vars), ' climate variables')

  NICHE.1KM    <- coord_df %>% filter(coord_summary == "TRUE")
  NICHE.1KM.84 <- SpatialPointsDataFrame(coords      = NICHE.1KM[c("lon", "lat")],
                                         data        = NICHE.1KM,
                                         proj4string = CRS.WGS.84)

  ## Use a projected, rather than geographic, coordinate system
  ## Not sure why, but this is needed for the spatial overlay step
  AUS.WGS      = spTransform(aus_shp,      CRS.WGS.84)
  LAND.WGS     = spTransform(world_shp,    CRS.WGS.84)
  KOP.WGS      = spTransform(kop_shp,      CRS.WGS.84)
  IBRA.WGS      = spTransform(ibra_shp,    CRS.WGS.84)

  ## Intersect the points with the Global koppen file
  message('Intersecting points with shapefiles for ', length(species_list), ' species')
  KOP.JOIN     = over(NICHE.1KM.84, KOP.WGS)
  IBRA.JOIN    = over(NICHE.1KM.84, IBRA.WGS)

  ## Create global niche and Australian niche for website - So we need a subset for Australia
  ## The ,] acts just like a clip in a GIS
  NICHE.AUS <-  NICHE.1KM.84[AUS.WGS, ]

  ## Aggregate the number of Koppen zones (and IBRA regions) each species is found in
  COMBO.KOP <- NICHE.1KM.84 %>%
    cbind.data.frame(., KOP.JOIN) %>%
    cbind.data.frame(., IBRA.JOIN)

  ## Aggregate the data
  KOP.AGG = tapply(COMBO.KOP$Koppen, COMBO.KOP$searchTaxon,
                   function(x) length(unique(x))) %>% ## group Koppen by species name
    as.data.frame()
  KOP.AGG =  setDT(KOP.AGG , keep.rownames = TRUE)[]
  names(KOP.AGG) = c("searchTaxon", "KOP_count")

  IBR.AGG = tapply(COMBO.KOP$REG_NAME_7, COMBO.KOP$searchTaxon,
                   function(x) length(unique(x))) %>% ## group Koppen by species name
    as.data.frame()
  IBR.AGG =  setDT(IBR.AGG , keep.rownames = TRUE)[]
  names(IBR.AGG) = c("searchTaxon", "IBRA_count")

  ## Run join between species records and spatial units :: SUA, POA and KOPPEN zones
  message('Joining occurence data to SUAs for ',
          length(species_list), ' species in the set ', "'", save_run, "'")
  projection(NICHE.1KM.84);projection(AUS.WGS)

  ## CREATE NICHES FOR PROCESSED TAXA
  ## Create niche summaries for each environmental condition like this...commit
  ## Here's what the function will produce :
  NICHE.AUS.DF = NICHE.AUS %>%
    as.data.frame() %>%
    dplyr::select(., searchTaxon, one_of(env_vars))

  NICHE.GLO.DF = NICHE.1KM.84 %>%
    as.data.frame() %>%
    dplyr::select(., searchTaxon, one_of(env_vars))

  ## Currently, the niche estimator can't handle NAs
  NICHE.AUS.DF = completeFun(NICHE.AUS.DF, env_vars[1])
  NICHE.GLO.DF = completeFun(NICHE.GLO.DF, env_vars[1])

  message('Estimating global niches for ',
          length(species_list), ' species in the set ', "'", save_run, "'")
  GLOB.NICHE <- env_vars %>%

    ## Pipe the list into lapply
    lapply(function(x) {

      ## Now, use the niche width function on each colname
      ## Also, need to figure out how to make the aggregating column generic (species, genus, etc.)
      ## currently it only works hard-wired........
      niche_estimate(DF = NICHE.GLO.DF, colname = x)

      ## Would be good to remove the duplicate columns here.....

    }) %>%

    ## finally, create one dataframe for all niches
    as.data.frame


  ## Remove duplicate Taxon columns and check the output
  GLOB.NICHE <- GLOB.NICHE %>% select(-contains("."))
  message('Calculated global niches for ', names(GLOB.NICHE), ' variables')

  message('Estimating Australian niches for ', length(species_list), ' species in the set ', "'", save_run, "'")
  AUS.NICHE <- env_vars %>%

    ## Pipe the list into lapply
    lapply(function(x) {

      ## Now, use the niche width function on each colname
      ## Also, need to figure out how to make the aggregating column generic (species, genus, etc.)
      niche_estimate(DF = NICHE.AUS.DF, colname = x)

    }) %>%

    ## finally, create one dataframe for all niches
    as.data.frame

  ## Remove duplicate Taxon columns and check the output :: would be great to skip these columns when running the function
  ## Remove duplicate Taxon columns and check the output
  ## Remove duplicate Taxon columns and check the output
  AUS.NICHE <- GLOB.NICHE %>% select(-contains("."))
  message('Calculated Australia niches for ', names(AUS.NICHE), ' variables')

  ## How are the AUS and GLOB niches related? Many species won't have both Australian and Global niches.
  ## So best to calculate the AUS niche as a separate table. Then, just use the global niche table for the rest of the code
  length(AUS.NICHE$searchTaxon); length(GLOB.NICHE$searchTaxon)

  ## Add the count of Australian records - this is not necessarily the same as maxent_records
  Aus_records = as.data.frame(table(NICHE.AUS.DF$searchTaxon))
  names(Aus_records) = c("searchTaxon", "Aus_records")
  identical(nrow(Aus_records), nrow(AUS.NICHE))

  ## Add the counts of Australian records for each species to the niche database
  GLOB.NICHE = join(Aus_records, GLOB.NICHE,  type = "right")

  ## Check the record and POA counts
  head(GLOB.NICHE$Aus_records,20)
  head(GLOB.NICHE$KOP_count,  20)

  ## Check
  message(round(nrow(AUS.NICHE)/nrow(GLOB.NICHE)*100, 2), " % species have records in Australia")
  dim(GLOB.NICHE)
  dim(AUS.NICHE)

  ## 5). CALCULATE AREA OF OCCUPANCY

  ## Given a dataframe of georeferenced occurrences of one, or more, taxa, this function provide statistics values
  ## (Extent of Occurrence, Area of Occupancy, number of locations, number of subpopulations) and provide a preliminary
  ## conservation status following Criterion B of IUCN.

  ## AREA OF OCCUPANCY (AOO).
  ## For every species in the list: calculate the AOO
  ## x = spp.geo[1]
  spp.geo = as.character(unique(COMBO.KOP$searchTaxon))

  GBIF.AOO <- spp.geo %>%

    ## Pipe the list into lapply
    lapply(function(x) {

      ## Subset the the data frame to calculate area of occupancy according the IUCN.eval
      DF   = subset(COMBO.KOP, searchTaxon == x)[, c("lat", "lon", "searchTaxon")]

      message('Calcualting geographic ranges for ', x, ', ', nrow(DF), ' records')
      AOO  = AOO.computing(XY = DF, Cell_size_AOO = cell_size)  ## Grid size in decimal degrees
      AOO  = as.data.frame(AOO)
      AOO$searchTaxon  <- rownames(AOO)
      rownames(AOO)    <- NULL

      ## Return the area
      return(AOO)

    }) %>%

    ## Finally, create one dataframe for all niches
    bind_rows

  head(GBIF.AOO)


  ## Now join on the geographic range and glasshouse data
  identical(nrow(GBIF.AOO), nrow(GLOB.NICHE))
  GLOB.NICHE <- list(GLOB.NICHE, GBIF.AOO, KOP.AGG, IBR.AGG) %>%
    reduce(left_join, by = "searchTaxon") %>%
    select(searchTaxon, Aus_records, AOO, KOP_count, IBRA_count, everything())

  if(save_data == "TRUE") {

    ## save .rds file for the next session
    message('Writing 1km resolution niche and raster data for ', length(species_list), ' species in the set ', "'", save_run, "'")
    saveRDS(GLOB.NICHE, paste0(DATA_path, 'GLOBAL_BAT_NICHE_',  save_run, '.rds'))

  } else {

    message(' Return niches to the environment for ', length(species_list), ' species analysed')
    return(GLOB.NICHE)

  }

}



## Plot histograms and convex hulls for selected taxa ----
#' @export
plot_range_histograms = function(coord_df,
                                 species_list,
                                 range_path) {

  ## Subset the occurrence data to that used by maxent
  message('Plotting global environmental ranges for ', length(species_list), ' species')
  CLEAN.TRUE <- coord_df %>% filter(coord_summary == "TRUE")

  ## Plot histograms of temperature and rainfall
  spp.plot = as.character(unique(CLEAN.TRUE$searchTaxon))
  for (spp in spp.plot) {

    ## Subset the spatial dataframe into records for each spp
    DF       <- coord_df[coord_df$searchTaxon %in% spp , ]
    DF.OCC   <- subset(coord_df, searchTaxon == spp & SOURCE != "INVENTORY")

    ## Plot convex Hull
    ## Check if the file exists
    convex.png = sprintf("./data/ANALYSIS/SPECIES_RANGES/%s_%s", spp, "_1km_convex_hull.png")
    if(!file.exists(convex.png)) {

      ## Start PNG device
      message('Writing global convex hulls for ', spp)
      png(sprintf("%s%s_%s", range_path, spp, "_1km_convex_hull.png"), 16, 10, units = 'in', res = 500)

      ## Create convex hull and plot
      find_hull <- function(DF) DF[chull(DF$Annual_mean_temp, DF$Annual_precip), ]
      hulls     <- ddply(DF, "SOURCE", find_hull)
      plot      <- ggplot(data = DF, aes(x = Annual_mean_temp,
                                         y = Annual_precip, colour = SOURCE, fill = SOURCE)) +
        geom_point() +
        geom_polygon(data = hulls, alpha = 0.5) +
        labs(x = "Annual_mean_temp", y = "Annual_precip") +

        theme(axis.title.x     = element_text(colour = "black", size = 35),
              axis.text.x      = element_text(size = 20),

              axis.title.y     = element_text(colour = "black", size = 35),
              axis.text.y      = element_text(size = 20),

              panel.background = element_blank(),
              panel.border     = element_rect(colour = "black", fill = NA, size = 1.5),
              plot.title       = element_text(size   = 40, face = "bold"),
              legend.text      = element_text(size   = 20),
              legend.title     = element_text(size   = 20),
              legend.key.size  = unit(1.5, "cm")) +

        ggtitle(paste0("Convex Hull for ", spp))
      print(plot)

      ## close device
      dev.off()

    } else {
      message('Convex hull already created for ', spp)
    }

    ## Plot temperature histograms
    temp.png = sprintf("%s%s_%s", range_path, spp, "temp_niche_histograms_1km_records.png")

    ## Check if file exists
    if(!file.exists(temp.png)) {
      message('Writing global temp histograms for ', spp)
      png(sprintf("%s%s_%s", range_path, spp, "temp_niche_histograms_1km_records.png"),
          16, 10, units = 'in', res = 500)

      ## Use the 'SOURCE' column to create a histogram for each source.
      temp.hist = ggplot(DF, aes(x = Annual_mean_temp, group = SOURCE, fill = SOURCE)) +

        geom_histogram(position = "identity", alpha = 0.5, binwidth = 0.25,
                       aes(y =..density..))  +
        geom_density(col = 4, alpha = 0.5) +

        ## Add some median lines : overall, ALA and GBIF
        geom_vline(aes(xintercept = median(DF$Annual_mean_temp)),
                   col = 'blue', size = 1) +
        geom_vline(aes(xintercept = median(DF.OCC$Annual_mean_temp)),
                   col = 'red', size = 1) +

        ggtitle(paste0("Worldclim temp niches for ", spp)) +

        ## Add themes
        theme(axis.title.x     = element_text(colour = "black", size = 35),
              axis.text.x      = element_text(size = 25),

              axis.title.y     = element_text(colour = "black", size = 35),
              axis.text.y      = element_text(size = 25),

              panel.background = element_blank(),
              panel.border     = element_rect(colour = "black", fill = NA, size = 3),
              plot.title       = element_text(size   = 40, face = "bold"),
              legend.text      = element_text(size   = 20),
              legend.title     = element_text(size   = 20),
              legend.key.size  = unit(1.5, "cm"))

      ## Print the plot and close the device
      print(temp.hist + ggtitle(paste0("Worldclim temp niches for ", spp)))
      dev.off()

    } else {
      message('Temp histogram already created for ', spp)
    }

    ## Plot rainfall histograms
    rain.png = sprintf("%s%s_%s", range_path, spp, "rain_niche_histograms_1km_records.png")

    ## Check if file exists
    if(!file.exists(rain.png)) {
      message('Writing global rain histograms for ', spp)

      message('Writing global rain histograms for ', spp)
      png(sprintf("%s%s_%s", range_path, spp, "rain_niche_histograms_1km_records.png"),
          16, 10, units = 'in', res = 500)

      ## Use the 'SOURCE' column to create a histogram for each source.
      rain.hist = ggplot(DF, aes(x = Annual_precip, group = SOURCE, fill = SOURCE)) +

        geom_histogram(position = "identity", alpha = 0.5, binwidth = 15,
                       aes(y =..density..))  +
        geom_density(col = 4, alpha = 0.5) +

        ## Add some median lines : overall, ALA and GBIF
        geom_vline(aes(xintercept = median(DF$Annual_precip)),
                   col = 'blue', size = 1) +
        geom_vline(aes(xintercept = median(DF.OCC$Annual_precip)),
                   col = 'red', size = 1) +

        ggtitle(paste0("Worldclim rain niches for ", spp)) +

        ## Add themes
        theme(axis.title.x     = element_text(colour = "black", size = 35),
              axis.text.x      = element_text(size = 25),

              axis.title.y     = element_text(colour = "black", size = 35),
              axis.text.y      = element_text(size = 25),

              panel.background = element_blank(),
              panel.border     = element_rect(colour = "black", fill = NA, size = 3),
              plot.title       = element_text(size   = 40, face = "bold"),
              legend.text      = element_text(size   = 20),
              legend.title     = element_text(size   = 20),
              legend.key.size  = unit(1.5, "cm"))

      ## Print the plot and close the device
      print(rain.hist + ggtitle(paste0("Worldclim rain niches for ", spp)))
      dev.off()

    } else {
      message('Rain histogram already created for ', spp)
    }
  }

}

