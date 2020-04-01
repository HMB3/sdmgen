#########################################################################################################################
######################################  FUNCTIONS FOR SDM ANALYSIS ---- #################################################
#########################################################################################################################


## Below is a list of the functions used to prepare the data for SDM analysis, and run the analysis
## Need to check the functions work the whole way through



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
                                   urban_df,
                                   add_urban,
                                   species_list,
                                   template_raster,
                                   thin_records,
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

  ## If Urban = TRUE
  if(urban_df != 'NONE') {

    urban_cols     <- intersect(names(GBIF.ALA.COMBO), names(urban_df))
    urban_df       <- select(urban_df, urban_cols)
    GBIF.ALA.COMBO <- bind_rows(GBIF.ALA.COMBO, urban_df)

  } else {
    message('Dont add urban data' )
  }

  ## CHECK TAXONOMY RETURNED BY ALA USING TAXONSTAND
  GBIF.ALA.MATCH = GBIF.ALA.COMBO

  ## Create points: the 'over' function seems to need geographic coordinates for this data...
  GBIF.ALA.84   = SpatialPointsDataFrame(coords      = GBIF.ALA.MATCH[c("lon", "lat")],
                                         data        = GBIF.ALA.MATCH,
                                         proj4string = projection)

  if(thin_records == TRUE) {

    ## The length needs to be the same
    length(unique(GBIF.ALA.84$searchTaxon))
    GBIF.ALA.84.SPLIT.ALL <- split(GBIF.ALA.84, GBIF.ALA.84$searchTaxon)
    occurrence_cells_all  <- lapply(GBIF.ALA.84.SPLIT.ALL, function(x) cellFromXY(template_raster, x))

    ## Check with a message, but could check with a fail
    message('Split prodcues ', length(occurrence_cells_all), ' data frames for ', length(species_list), ' species')

    ## Now get just one record within each 1*1km cell.
    GBIF.ALA.84.1KM <- mapply(function(x, cells) {
      x[!duplicated(cells), ]
    }, GBIF.ALA.84.SPLIT.ALL, occurrence_cells_all, SIMPLIFY = FALSE) %>% do.call(rbind, .)

    ## Check to see we have 19 variables + the species for the standard predictors, and 19 for all predictors
    message(round(nrow(GBIF.ALA.84.1KM)/nrow(GBIF.ALA.84)*100, 2), " % records retained at 1km resolution")

    ## Create points: the 'over' function seems to need geographic coordinates for this data...
    COMBO.POINTS   = GBIF.ALA.84.1KM[c("lon", "lat")]

  } else {
    message('dont thin the records out' )
    COMBO.POINTS = GBIF.ALA.84
  }

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




## Extract environmental values for urban occurrence records -----
#' @export
urban_records_extract = function(urban_df,
                                 species_list,
                                 thin_records,
                                 template_raster,
                                 world_raster,
                                 projection,
                                 biocl_vars,
                                 env_vars,
                                 worldclim_grids,
                                 save_data,
                                 save_run) {

  ## Get just the species list
  URBAN.XY = urban_df[urban_df$searchTaxon %in% species_list, ]
  message('Extracting raster values for ',
          length(unique(URBAN.XY$searchTaxon)), ' urban species across ',
          length(unique(URBAN.XY$INVENTORY)),   ' Councils ')

  ## Create points: the 'over' function seems to need geographic coordinates for this data...
  URBAN.XY   = SpatialPointsDataFrame(coords      = URBAN.XY[c("lon", "lat")],
                                      data        = URBAN.XY,
                                      proj4string = projection)

  if(thin_records == TRUE) {

    ## The length needs to be the same
    length(unique(URBAN.XY$searchTaxon))
    URBAN.XY.SPLIT.ALL <- split(URBAN.XY, URBAN.XY$searchTaxon)
    occurrence_cells_all  <- lapply(URBAN.XY.SPLIT.ALL, function(x) cellFromXY(template_raster, x))

    ## Check with a message, but could check with a fail
    message('Split prodcues ', length(occurrence_cells_all), ' data frames for ', length(species_list), ' species')

    ## Now get just one record within each 10*10km cell.
    URBAN.XY.1KM <- mapply(function(x, cells) {
      x[!duplicated(cells), ]
    }, URBAN.XY.SPLIT.ALL, occurrence_cells_all, SIMPLIFY = FALSE) %>% do.call(rbind, .)

    ## Check to see we have 19 variables + the species for the standard predictors, and 19 for all predictors
    message(round(nrow(URBAN.XY.1KM)/nrow(URBAN.XY)*100, 2), " % records retained at 1km resolution")

    ## Create points: the 'over' function seems to need geographic coordinates for this data...
    COMBO.POINTS = URBAN.XY.1KM[c("lon", "lat")]

  } else {
    message('dont thin the records out' )
    COMBO.POINTS = URBAN.XY
  }

  ## Bioclim variables
  ## Extract raster data
  message('Extracting raster values for ', length(species_list), ' species in the set ', "'", save_run, "'")
  message(projection(COMBO.POINTS));message(projection(world_raster))
  dim(COMBO.POINTS);dim(URBAN.XY.1KM)

  ## Extract the raster values
  COMBO.RASTER <- raster::extract(world_raster, COMBO.POINTS) %>%
    cbind(as.data.frame(URBAN.XY.1KM), .)

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
                                  urban_df,
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

      ## Now add attach column for spp, and the flag for each record
      d = cbind(searchTaxon = x,
                SPAT_OUT = sp.flag, f)[c("searchTaxon", "SPAT_OUT", "CC.OBS")]

      ## Remember to explicitly return the df at the end of loop, so we can bind
      return(d)

    }) %>%

    ## Finally, bind all the rows together
    bind_rows

  gc()

  ## Join the data back on
  SPAT.FLAG = join(as.data.frame(test.geo), SPAT.OUT)    ## Join means the skipped spp are left out
  dim(SPAT.FLAG)


  ## Try plotting the points which are outliers for a subset of spp and label them
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

  ## Could add the species in urban records in here
  if(urban_df != 'NONE') {

    ##
    message('Combinethe Spatially cleaned data with the Urban data' )
    urban_cols  <- intersect(names(SPAT.FLAG), names(urban_df)) %>% intersect(., names(all_df))
    # urban_df    <- as.data.frame(urban_df) %>%
    #   select(., urban_cols)

    SPAT.TRUE <- SPAT.FLAG %>%
      filter(SPAT_OUT == "TRUE")
    SPAT.TRUE <- SpatialPointsDataFrame(coords      = SPAT.TRUE[c("lon", "lat")],
                                        data        = SPAT.TRUE,
                                        proj4string = CRS.WGS.84)

  } else {
    message('Dont add urban data' )
  }

  return(SPAT.TRUE)

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

  NICHE.1KM    <- coord_df %>% as.data.frame () %>%
    filter(coord_summary == "TRUE")
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





## Create SDM table ----
#' @export
prepare_sdm_table = function(coord_df,
                             species_list,
                             BG_points,
                             sdm_table_vars,
                             save_run,
                             save_shp,
                             read_background,
                             save_data) {

  ## Just add clean_df to this step
  coord_df <- subset(coord_df, coord_summary == "TRUE")
  message(round(nrow(coord_df)/nrow(coord_df)*100, 2), " % records retained")
  message(table(coord_df$coord_summary))

  ## RESTRICT DATA TABLE TO ONE RECORD PER 1KM GRID CELL
  ## Create a table with all the variables needed for SDM analysis
  message('Preparing SDM table for ', length(unique(coord_df$searchTaxon)), ' species in the set ', "'", save_run, "'",
          'using ', unique(coord_df$SOURCE), ' data')

  ## Select only the columns needed. This also needs to use the variable names
  coord_df         <- coord_df[coord_df$searchTaxon %in% species_list, ]
  length(unique(coord_df$searchTaxon))

  COMBO.RASTER.ALL  <- coord_df %>%
    select(one_of(sdm_table_vars))


  ## Create a spatial points object, and change to a projected system to calculate distance more accurately
  ## This is the mollweide projection used for the SDMs
  coordinates(COMBO.RASTER.ALL)    <- ~lon+lat
  proj4string(COMBO.RASTER.ALL)    <- '+init=epsg:4326'
  COMBO.RASTER.ALL                 <- spTransform(COMBO.RASTER.ALL, CRS(sp_epsg54009))

  ## Don't filter the data again to be 1 record per 1km, that has already happened
  SDM.DATA.ALL <- COMBO.RASTER.ALL

  ## Check to see we have 19 variables + the species for the standard predictors, and 19 for all predictors
  message(round(nrow(SDM.DATA.ALL)/nrow(coord_df)*100, 2), " % records retained at 1km resolution")

  ## What resolution is the template raster at?
  xres(template.raster.1km);yres(template.raster.1km)

  ## TRY CLEANING FILTERED DATA FOR SPATIAL OUTLIERS
  ## The cc_outl function has been tweaked and sped up.

  ## Create a unique identifier for spatial cleaning. This is used for automated cleaing of the records, and also saving shapefiles
  ## But this will not be run for all species linearly. So, it probably needs to be a combination of species and number
  SDM.DATA.ALL$SPOUT.OBS <- 1:nrow(SDM.DATA.ALL)
  SDM.DATA.ALL$SPOUT.OBS <- paste0(SDM.DATA.ALL$SPOUT.OBS, "_SPOUT_", SDM.DATA.ALL$searchTaxon)
  SDM.DATA.ALL$SPOUT.OBS <- gsub(" ",     "_",  SDM.DATA.ALL$SPOUT.OBS, perl = TRUE)
  length(SDM.DATA.ALL$SPOUT.OBS);length(unique(SDM.DATA.ALL$SPOUT.OBS))

  ## Check dimensions
  dim(SDM.DATA.ALL)
  length(unique(SDM.DATA.ALL$searchTaxon))
  length(unique(SDM.DATA.ALL$SPOUT.OBS))
  unique(SDM.DATA.ALL$SOURCE)

  ## Create a tibble to supply to coordinate cleaner
  SDM.COORDS  <- SDM.DATA.ALL %>%
    spTransform(., CRS.WGS.84) %>%
    as.data.frame() %>%
    select(searchTaxon, lon, lat, SPOUT.OBS, SOURCE) %>%
    dplyr::rename(species          = searchTaxon,
                  decimallongitude = lon,
                  decimallatitude  = lat) %>%
    timetk::tk_tbl()

  ## Check
  message(identical(SDM.COORDS$index, SDM.COORDS$SPOUT.OBS))
  length(unique(SDM.COORDS$species))

  ## Check how many records each species has
  COMBO.LUT <- SDM.COORDS %>%
    as.data.frame() %>%
    select(species) %>%
    table() %>%
    as.data.frame()
  COMBO.LUT <- setDT(COMBO.LUT, keep.rownames = FALSE)[]
  names(COMBO.LUT) = c("species", "FREQUENCY")
  COMBO.LUT = COMBO.LUT[with(COMBO.LUT, rev(order(FREQUENCY))), ]

  ## Watch out here - this sorting could cause problems for the order of the data frame once it's stitched back together
  ## If we we use species to join the data back together, will it preserve the order?
  LUT.100K = as.character(subset(COMBO.LUT, FREQUENCY < 100000)$species)
  LUT.100K = trimws(LUT.100K [order(LUT.100K)])
  length(LUT.100K)

  ## Create a data frame of species name and spatial outlier
  SPAT.OUT <- LUT.100K  %>%

    ## pipe the list of species into lapply
    lapply(function(x) {

      ## Create the species df by subsetting by species
      f <- subset(SDM.COORDS, species == x)

      ## Run the spatial outlier detection
      message("Running spatial outlier detection for ", nrow(f), " records for ", x)
      sp.flag <- cc_outl(f,
                         lon     = "decimallongitude",
                         lat     = "decimallatitude",
                         species = "species",
                         method  = "quantile", #"distance",
                         mltpl   = 10,
                         #tdi     = 300,
                         value   = "flagged",
                         verbose = "TRUE")

      ## Now add attache column for species, and the flag for each record
      d = cbind(searchTaxon = x,
                SPAT_OUT = sp.flag, f)[c("searchTaxon", "SPAT_OUT", "SPOUT.OBS")]

      ## Remeber to explicitly return the df at the end of loop, so we can bind
      return(d)

    }) %>%

    ## Finally, bind all the rows together
    bind_rows

  gc()

  ## How many species are flagged as spatial outliers?
  print(table(SPAT.OUT$SPAT_OUT, exclude = NULL))
  length(unique(SPAT.OUT$searchTaxon))
  head(SPAT.OUT)

  ## FILTER DATA TO REMOVE SPATIAL OUTLIERS
  ## Join data :: Best to use the 'OBS' column here
  message('Is the order or records identical before joining?', identical(nrow(SDM.COORDS), nrow(SPAT.OUT)))
  message('Is the order or records identical after joining?', identical(SDM.DATA.ALL$searchTaxon, SPAT.OUT$searchTaxon))
  length(unique(SPAT.OUT$searchTaxon))

  ## This explicit join is required. Check the species have been analysed in exactly the same order
  SPAT.FLAG <- join(as.data.frame(SDM.DATA.ALL), SPAT.OUT,
                    by = c("SPOUT.OBS", "searchTaxon") ,
                    type = "left", match = "first")
  message('Is the order or records identical after joining?',identical(SDM.DATA.ALL$searchTaxon, SPAT.FLAG$searchTaxon))

  ## Check the join is working
  message('Checking spatial flags for ', length(unique(SPAT.FLAG$searchTaxon)), ' species in the set ', "'", save_run, "'")
  message(table(SPAT.FLAG$SPAT_OUT, exclude = NULL))

  ## Just get the records that were not spatial outliers.
  SDM.SPAT.ALL <- subset(SPAT.FLAG, SPAT_OUT == "TRUE")
  unique(SDM.SPAT.ALL$SPAT_OUT)
  unique(SDM.SPAT.ALL$SOURCE)
  length(unique(SDM.SPAT.ALL$searchTaxon))

  ## What percentage of records are retained?
  message(round(nrow(SDM.SPAT.ALL)/nrow(SPAT.FLAG)*100, 2), " % records retained after spatial outlier detection")

  ## Convert back to format for SDMs :: use Mollweide projection
  SDM.SPAT.ALL = SpatialPointsDataFrame(coords      = SDM.SPAT.ALL[c("lon", "lat")],
                                        data        = SDM.SPAT.ALL,
                                        proj4string = CRS(sp_epsg54009))
  projection(SDM.SPAT.ALL)
  message(length(unique(SDM.SPAT.ALL$searchTaxon)), ' species processed through from download to SDM table')

  ## CREATE SHAPEFILES TO CHECK OUTLIERS ARCMAP
  ## Rename the fields so that ArcMap can handle them
  SPAT.OUT.CHECK = SPAT.FLAG %>%
    select(SPOUT.OBS, searchTaxon, lat, lon, SOURCE, SPAT_OUT) %>%
    dplyr::rename(TAXON     = searchTaxon,
                  LAT       = lat,
                  LON       = lon)
  names(SPAT.OUT.CHECK)

  ## Then create a SPDF
  SPAT.OUT.SPDF    = SpatialPointsDataFrame(coords      = SPAT.OUT.CHECK[c("LON", "LAT")],
                                            data        = SPAT.OUT.CHECK,
                                            proj4string = CRS.WGS.84)

  ## Write the shapefile out
  if(save_shp == "TRUE") {

    ## save .shp for future refrence
    writeOGR(obj    = SPAT.OUT.SPDF,
             dsn    = "./data/ANALYSIS/CLEAN_GBIF",
             layer  = paste0('SPAT_OUT_CHECK_', save_run),
             driver = "ESRI Shapefile", overwrite_layer = TRUE)

  } else {

    message(' skip file saving, not many species analysed')   ##

  }

  ## CREATE BACKGROUND POINTS AND VARIBALE NAMES
  ## Use one data frame for all species analysis, to save mucking around with background points
  ## Here we are using a dataframe of mammals, reptiles and birds.
  if(read_background == "TRUE") {

    Message('Read in background data for taxa analaysed')
    background.points = readRDS(paste0(DATA_path, BG_points)) %>%
      .[, -(27:28)]
    background.points = background.points[!background.points$searchTaxon %in% species_list, ]

    ## The BG points step needs to be ironed out.
    ## For some analysis, we need to do other taxa (e.g. animals)
    ## For Bats, they just want to use other bats.
    SDM.SPAT.OCC.BG = rbind(SDM.SPAT.ALL, background.points)

  } else {

    message('Dont read in Background data, creating it in this run')
    SDM.SPAT.OCC.BG = SDM.SPAT.ALL

  }

  ## Now select the final columns needed
  message(length(unique(SDM.SPAT.OCC.BG$searchTaxon)), ' Species in occurrence and BG data')
  drops <- c('CC.OBS', 'SPOUT.OBS')
  sdm_cols <- names(select(SDM.SPAT.OCC.BG@data, searchTaxon, lon, lat, SOURCE, everything()))
  SDM.SPAT.OCC.BG <- SDM.SPAT.OCC.BG[,!(names(SDM.SPAT.OCC.BG) %in% drops)]
  #SDM.SPAT.OCC.BG <- SDM.SPAT.OCC.BG[,(names(SDM.SPAT.OCC.BG) %in% sdm_cols)]

  ## save data
  if(save_data == "TRUE") {

    ## Save .rds file of the occurrence and BG points for the next session
    saveRDS(SDM.SPAT.OCC.BG, paste0(DATA_path, 'SDM_SPAT_OCC_BG_',  save_run, '.rds'))

  } else {
    message('Return the occurrence + Background data to the global environment')   ##
    return(SDM.SPAT.OCC.BG)
  }

  ## get rid of some memory
  gc()

}




## Run the SDM analysis ----
#' @export
run_sdm_analysis = function(sdm_df,
                            species_list,
                            sdm.predictors,
                            maxent_dir,
                            #maxent_path,
                            bs_dir,
                            backwards_sel,
                            cor_thr,
                            pct_thr,
                            k_thr,
                            min_n,
                            max_bg_size,
                            background_buffer_width,
                            shapefiles,
                            features,
                            replicates,
                            responsecurves,
                            Koppen,
                            shp_path,
                            aus_shp) {

  ## Loop over all the species
  lapply(species_list, function(species){

    ## Skip the species if the directory already exists, before the loop
    outdir <- maxent_dir

    ## Could also check the folder size, so folders with no contents aren't skipped eg
    ## && sum(file.info(dir_name)$size) < 1000 (EG 1MB)
    dir_name = file.path(outdir, gsub(' ', '_', species))
    if(dir.exists(dir_name)) {
      message('Skipping ', species, ' - already run.')
      invisible(return(NULL))

    }

    ## Create the directory for the species in progress,
    ## so other parallel runs don't run the same species
    dir.create(dir_name)
    file.create(file.path(dir_name, "in_progress.txt"))

    ## Print the taxa being processed to screen
    if(species %in% SDM.SPAT.OCC.BG$searchTaxon) {
      message('Doing ', species)

      ## Subset the records to only the taxa being processed
      occurrence <- subset(SDM.SPAT.OCC.BG, searchTaxon == species)
      message('Using ', nrow(occurrence), ' occ records from ', unique(occurrence$SOURCE))

      ## Now get the background points. These can come from any species, other than the modelled species.
      ## However, they should be limited to the same SOURCE as the occ data
      background <- subset(SDM.SPAT.OCC.BG, searchTaxon != species)
      message('Using ', nrow(background), ' background records from ', unique(background$SOURCE))

      ## Finally fit the models using FIT_MAXENT_TARG_BG. Also use tryCatch to skip any exceptions
      tryCatch(
        fit_maxent_targ_bg_back_sel(occ                     = occurrence,    ## name from the .rmd CV doc
                                    bg                      = background,    ## name from the .rmd CV doc
                                    sdm.predictors          = bs.predictors,
                                    name                    = species,
                                    outdir                  = maxent_dir,
                                    bsdir                   = bs_dir,
                                    backwards_sel           = "TRUE",
                                    cor_thr                 = cor_thr,
                                    pct_thr                 = pct_thr,
                                    k_thr                   = k_thr,

                                    template.raster         = template.raster.1km,
                                    min_n                   = min_n,
                                    max_bg_size             = max_bg_size,
                                    Koppen                  = Koppen,
                                    background_buffer_width = background_buffer_width,
                                    shapefiles              = shapefiles,
                                    features                = features,
                                    replicates              = replicates,
                                    responsecurves          = responsecurves,
                                    #shp_path                = shp_path,
                                    aus_shp                 = aus_shp),


        ## Save error message
        error = function(cond) {

          ## How to put the message into the file?
          file.create(file.path(dir_name, "sdm_failed.txt"))
          message(species, ' failed')
          cat(cond$message, file = file.path(dir_name, "sdm_failed.txt"))
          warning(species, ': ', cond$message)

        })

    } else {

      message(species, ' skipped - no data.')         ## This condition ignores species which have no data...
      file.create(file.path(dir_name, "completed.txt"))

    }

    ## now add a file to the dir to denote that it has completed
    file.create(file.path(dir_name, "completed.txt"))

  })


}




## Run maxent with backwards selection
#' @export
fit_maxent_targ_bg_back_sel <- function(occ,
                                        bg, # A Spatial points data frame (SPDF) of candidate background points
                                        sdm.predictors,
                                        # sdm.predictors is a vector of enviro conditions that you want to include
                                        name,
                                        outdir,
                                        bsdir,
                                        cor_thr,
                                        pct_thr,
                                        k_thr,
                                        backwards_sel,
                                        template.raster,
                                        # template.raster is an empty raster with extent, res and projection
                                        # of final output rasters. It is used to reduce
                                        # occurrences to a single point per cell.
                                        min_n,
                                        max_bg_size,
                                        background_buffer_width,
                                        Koppen,
                                        shapefiles,
                                        features,
                                        replicates,
                                        responsecurves,
                                        rep_args,
                                        full_args,
                                        # shp_path,
                                        aus_shp) {

  #####
  ## First, stop if the outdir file exists,
  if(!file.exists(outdir)) stop('outdir does not exist :(', call. = FALSE)
  outdir_sp <- file.path(outdir, gsub(' ', '_', name))
  bsdir_sp  <- file.path(bsdir,  gsub(' ', '_', name))

  if(!missing('Koppen')) {
    if(!is(Koppen, 'RasterLayer'))
      stop('Koppen must be a RasterLayer, and should be in the same coordinate system as template.raster')
  }

  ## If the file doesn't exist, split out the features
  if(!file.exists(outdir_sp)) dir.create(outdir_sp)
  features <- unlist(strsplit(features, ''))

  ## Make sure user features are allowed: don't run the model if the
  ## features have been incorrectly specified in the main argument
  ## l: linear
  ## p: product
  ## q: quadratic
  ## h: hinge       ## disabled for this analysis
  ## t: threshold   ## disabled for this analysis
  if(length(setdiff(features, c('l', 'p', 'q', 'h', 't'))) > 1)
    stop("features must be a vector of one or more of ',
         'l', 'p', 'q', 'h', and 't'.")

  ## Create a buffer of xkm around the occurrence points
  buffer <- aggregate(gBuffer(occ, width = background_buffer_width, byid = TRUE))

  ##
  ## Get unique cell numbers for species occurrences
  cells <- cellFromXY(template.raster, occ)

  ## Clean out duplicate cells and NAs (including points outside extent of predictor data)
  ## Note this will get rid of a lot of duplicate records not filtered out by GBIF columns, etc.
  not_dupes <- which(!duplicated(cells) & !is.na(cells))
  occ       <- occ[not_dupes, ]
  cells     <- cells[not_dupes]
  message(nrow(occ), ' occurrence records (unique cells).')

  ##
  ## Skip species that have less than a minimum number of records: eg 20 species
  if(nrow(occ) < min_n) {

    print (paste ('Fewer occurrence records than the number of cross-validation ',
                  'replicates for species ', name,
                  ' Model not fit for this species'))

  } else {

    ##
    ## Subset the background records to the 200km buffered polygon
    message(name, ' creating background cells')
    system.time(o <- over(bg, buffer))
    bg <- bg[which(!is.na(o)), ]
    bg_cells <- cellFromXY(template.raster, bg)

    ## Clean out duplicates and NAs (including points outside extent of predictor data)
    bg_not_dupes <- which(!duplicated(bg_cells) & !is.na(bg_cells))
    bg           <- bg[bg_not_dupes, ]
    bg_cells     <- bg_cells[bg_not_dupes]

    ## Find which of these cells fall within the Koppen-Geiger zones that the species occupies
    ## Crop the Kopppen raster to the extent of the occurrences, and snap it
    message(name, ' intersecting background cells with Koppen zones')
    Koppen_crop <- crop(Koppen, occ, snap = 'out')

    ## Only extract and match those cells that overlap between the ::
    ## 1). cropped koppen zone,
    ## 2). occurrences and
    ## 3). background points
    message(xres(template.raster), ' metre cell size for template raster')
    message(xres(Koppen), ' metre cell size for Koppen raster')
    zones               <- raster::extract(Koppen_crop, occ)
    cells_in_zones_crop <- Which(Koppen_crop %in% zones, cells = TRUE)
    cells_in_zones      <- cellFromXY(Koppen, xyFromCell(Koppen_crop, cells_in_zones_crop))
    bg_cells            <- intersect(bg_cells, cells_in_zones)  ## this is 0 for 5km
    i                   <- cellFromXY(template.raster, bg)
    bg                  <- bg[which(i %in% bg_cells), ]

    ## For some species, we have the problem that the proportion of ALA/INV data is
    ## very different in the occurrence vs the bg records.
    ## This should be caused by the 200km / koppen restriction, etc.

    ##
    ## Reduce background sample, if it's larger than max_bg_size
    if (nrow(bg) > max_bg_size) {

      message(nrow(bg), ' target species background records for ', name,
              ', reduced to random ', max_bg_size, ' using random points from :: ', unique(bg$SOURCE))
      bg.samp <- bg[sample(nrow(bg), max_bg_size), ]

    } else {
      ## If the bg points are smaller that the max_bg_size, just get all the points
      message(nrow(bg), ' target species background records for ', name,
              ' using all points from :: ', unique(bg$SOURCE))
      bg.samp <- bg

    }

    ##
    ## Now save an image of the background points
    ## This is useful to quality control the models
    save_name = gsub(' ', '_', name)

    aus.mol = aus_shp %>%
      spTransform(projection(buffer))

    aus.kop = crop(Koppen_crop, aus.mol)
    occ.mol <- occ %>%
      spTransform(projection(buffer))

    ## Print the koppen zones, occurrences and points to screen
    plot(Koppen_crop, legend = FALSE,
         main = paste0('Occurence SDM records for ', name))

    plot(aus.mol, add = TRUE)
    plot(buffer,  add = TRUE, col = "red")
    plot(occ.mol, add = TRUE, col = "blue")

    ## Then save the occurrence points
    png(sprintf('%s/%s/%s_%s.png', outdir, save_name, save_name, "buffer_occ"),
        16, 10, units = 'in', res = 300)

    plot(Koppen_crop, legend = FALSE,
         main = paste0('Occurence SDM records for ', name))

    plot(aus.mol, add = TRUE)
    plot(buffer,  add = TRUE, col = "red")
    plot(occ.mol, add = TRUE, col = "blue")

    dev.off()

    ##
    ## Now save the buffer, the occ and bg points as shapefiles
    if(shapefiles) {

      suppressWarnings({

        message(name, ' writing occ and bg shapefiles')
        writeOGR(SpatialPolygonsDataFrame(buffer, data.frame(ID = seq_len(length(buffer)))),
                 outdir_sp, paste0(save_name, '_bg_buffer'),          'ESRI Shapefile', overwrite_layer = TRUE)
        writeOGR(bg.samp,  outdir_sp, paste0(save_name, '_bg'),       'ESRI Shapefile', overwrite_layer = TRUE)
        writeOGR(occ,           outdir_sp, paste0(save_name, '_occ'), 'ESRI Shapefile', overwrite_layer = TRUE)

      })

    }

    ## Also save the background and occurrence points as .rds files
    saveRDS(bg.samp,  file.path(outdir_sp, paste0(save_name, '_bg.rds')))
    saveRDS(occ,      file.path(outdir_sp, paste0(save_name, '_occ.rds')))

    ##
    ## SWD = species with data. Now sample the environmental
    ## variables used in the model at all the occ and bg points
    swd_occ <- occ[, sdm.predictors]
    saveRDS(swd_occ, file.path(outdir_sp, paste0(save_name,'_occ_swd.rds')))

    swd_bg <- bg.samp[, sdm.predictors]
    saveRDS(swd_bg, file.path(outdir_sp, paste0(save_name, '_bg_swd.rds')))

    ## Save the SWD tables as shapefiles
    if(shapefiles) {

      writeOGR(swd_occ, outdir_sp,  paste0(save_name, '_occ_swd'), 'ESRI Shapefile', overwrite_layer = TRUE)
      writeOGR(swd_bg,  outdir_sp,  paste0(save_name, '_bg_swd'),  'ESRI Shapefile', overwrite_layer = TRUE)

    }

    ##
    ## Now combine the occurrence and background data
    swd <- as.data.frame(rbind(swd_occ@data, swd_bg@data))
    saveRDS(swd, file.path(outdir_sp, 'swd.rds'))
    pa  <- rep(1:0, c(nrow(swd_occ), nrow(swd_bg)))

    ## Now, set the features to be used by maxent ::
    ## Linear, product and quadratic
    off <- setdiff(c('l', 'p', 'q', 't', 'h'), features)

    ## This sets threshold and hinge features to "off"
    if(length(off) > 0) {

      off <- c(l = 'linear=false',    p = 'product=false', q = 'quadratic=false',
               t = 'threshold=false', h = 'hinge=false')[off]

    }

    off <- unname(off)

    ## Run the replicates
    if(replicates > 1) {

      #if(missing(rep_args))
      rep_args <- NULL

      ## Run MAXENT for cross validation data splits of swd : so 5 replicaes, 0-4
      ## first argument is the predictors, the second is the occurrence data
      message(name, ' running xval maxent')
      me_xval <- maxent(swd, pa, path = file.path(outdir_sp, 'xval'),
                        args = c(paste0('replicates=', replicates),
                                 'responsecurves=true',
                                 'outputformat=logistic',
                                 off, paste(names(rep_args), rep_args, sep = '=')))

    }

    ## Run the full maxent model - using all the data in swd
    ## This uses DISMO to output standard files, but the names can't be altered
    #if(missing(full_args))
    full_args <- NULL
    message(name, ' running full maxent')
    me_full <- maxent(swd, pa, path = file.path(outdir_sp, 'full'),
                      args = c(off, paste(names(full_args), full_args, sep = '='),
                               'responsecurves=true',
                               'outputformat=logistic'))

    ## Save the full model. Replicate this line in the backwards selection algortithm
    ## Also worth checking that the koppen zones can be used at any resolution
    ## This is hard coded...hard to make it an argument
    saveRDS(list(me_xval = me_xval, me_full = me_full, swd = swd, pa = pa,
                 koppen_gridcode = as.character(Koppen_zones$Koppen[match(unique(zones), Koppen_zones$GRIDCODE)])),
            file.path(outdir_sp, 'full', 'maxent_fitted.rds'))

    if (backwards_sel == "TRUE") {

      ##
      ## Coerce the "species with data" (SWD) files to regular data.frames
      ## This is needed to use the simplify function
      swd_occ     <- as.data.frame(swd_occ)
      swd_occ$lon <- NULL
      swd_occ$lat <- NULL
      swd_bg      <- as.data.frame(swd_bg)
      swd_bg$lon  <- NULL
      swd_bg$lat  <- NULL

      ## Need to create a species column here
      swd_occ$searchTaxon <- name
      swd_bg$searchTaxon  <- name

      ##
      ## Run simplify rmaxent::simplify

      # Given a candidate set of predictor variables, this function identifies
      # a subset that meets specified multicollinearity criteria. Subsequently,
      # backward stepwise variable selection (VIF) is used to iteratively drop
      # the variable that contributes least to the model, until the contribution
      # of each variable meets a specified minimum, or until a predetermined
      # minimum number of predictors remains. It returns a model object for the
      # full model, rather than a list of models as does the previous function

      ## Using a modified versionof rmaxent::simplify, so that the name of the
      ## maxent model object "maxent_fitted.rds" is the same in both models.
      ## This is needed to run the mapping step over either the full or BS folder
      m <- local_simplify(

        swd_occ,
        swd_bg,
        path            = bsdir,
        species_column  = "searchTaxon",
        replicates      = replicates,  ## 5 as above
        response_curves = TRUE,
        logistic_format = TRUE,
        cor_thr         = cor_thr,
        pct_thr         = pct_thr,
        k_thr           = k_thr,
        features        = features,    ## LPQ as above
        quiet           = FALSE)

      ## Save the bg, occ and swd files into the backwards selection folder too
      saveRDS(bg.samp,  file.path(bsdir_sp, paste0(save_name, '_bg.rds')))
      saveRDS(occ,      file.path(bsdir_sp, paste0(save_name, '_occ.rds')))
      saveRDS(swd,      file.path(bsdir_sp, paste0('swd.rds')))

      ## Read the model in, because it's tricky to index
      bs.model <- readRDS(sprintf('%s/%s/full/maxent_fitted.rds', bsdir,  save_name))
      identical(length(bs.model@presence$Annual_mean_temp), nrow(occ))

      ##
      ## Save the chart corrleation file too for the training data set
      par(mar   = c(3, 3, 5, 3),
          oma   = c(1.5, 1.5, 1.5, 1.5))

      png(sprintf('%s/%s/full/%s_%s.png', bsdir,
                  save_name, save_name, "bs_predictor_correlation"),
          3236, 2000, units = 'px', res = 300)

      ## set margins
      par(mar   = c(3, 3, 5, 3),
          oma   = c(1.5, 1.5, 1.5, 1.5))

      ## Add detail to the response plot
      chart.Correlation(bs.model@presence,
                        histogram = TRUE, pch = 19,
                        title = paste0('Reduced variable correlations for ', save_name))

      dev.off()

    } else {

      message("Don't run backwards selection")

    }

  }

}




## Compile the SDM data ----
#' @export
compile_sdm_results = function(species_list,
                               results_dir,
                               save_data,
                               DATA_path,
                               save_run) {

  ## The code that adds niche info is now in './R/9_COLLATE_MAXENT_RESULTS.R'
  message('Creating summary stats for ', length(species_list), ' species in the set ', "'", save_run, "'")

  ## First, make a list of all the species with models, then restrict them to just the models on the species_list list
  map_spp_list  = gsub(" ", "_", species_list)
  map_spp_patt  = paste0(map_spp_list, collapse = "|")
  message ("map_spp_list head:")
  message (paste (head(map_spp_list), collapse=","))

  ## Now stop R from creating listing all the maxent files that have completed - this takes a long time
  #message ("DEBUGDEBUG - remember to disable next line")
  message('Compile SDM results for species in ', results_dir)
  maxent.tables = lapply (map_spp_list, FUN = function (x) {paste(results_dir , x, "full/maxent_fitted.rds", sep="/")})

  ## How many species have been modelled?
  message(paste("maxent.tables has this many entries:", length(maxent.tables)))
  message(paste(head (maxent.tables), collapse=","))
  sdm.exists = lapply(maxent.tables, FUN = function (x) {file.exists (x)}) %>% unlist()

  ## Only list the intersection between the modelled species and
  message(paste(head(sdm.exists), collapse=","))
  maxent.tables = maxent.tables[sdm.exists]

  ## Report what data the table has
  message (paste ("maxent.tables has this many entries:", length(maxent.tables)))
  maxent.tables = stringr::str_subset(maxent.tables, map_spp_patt)
  message (paste ("maxent.tables has this many entries:", length(maxent.tables)))
  message (paste (head(maxent.tables), collapse=","))

  ## Now create a table of the results
  ## x = maxent.tables[1]
  MAXENT.RESULTS <- maxent.tables %>%

    ## Pipe the list into lapply
    lapply(function(x) {

      ## We don't need to skip species that haven't been modelled
      x = gsub(paste0(results_dir, "/"), "", x)
      message (x)

      ## load the backwards selected model
      if (grepl("back", results_dir)) {
        m = readRDS(paste0(results_dir, '/',  x))

      } else {
        ## Get the background records from any source
        m = readRDS(paste0(results_dir, '/',  x))$me_full

      }

      ## Get the number of Variables
      number.var  = length(m@lambdas) - 4   ## (the last 4 slots of the lambdas file are not variables)
      mxt.records = nrow(m@presence)

      ## Get variable importance
      m.vars    = ENMeval::var.importance(m)
      var.pcont = m.vars[rev(order(m.vars[["percent.contribution"]])),][["variable"]][1]
      pcont     = m.vars[rev(order(m.vars[["percent.contribution"]])),][["percent.contribution"]][1]

      var.pimp  = m.vars[rev(order(m.vars[["permutation.importance"]])),][["variable"]][1]
      pimp      = m.vars[rev(order(m.vars[["permutation.importance"]])),][["permutation.importance"]][1]

      ## Get maxent results columns to be used for model checking
      ## Including the omission rate here
      Training_AUC             = m@results["Training.AUC",]
      Number_background_points = m@results["X.Background.points",]
      Logistic_threshold       = m@results["X10.percentile.training.presence.Logistic.threshold",]
      Omission_rate            = m@results["X10.percentile.training.presence.training.omission",]

      ## Now rename the maxent columns that we will use in the results table
      d = data.frame(searchTaxon              = x,
                     Maxent_records           = mxt.records,
                     Number_var               = number.var,
                     Var_pcont                = var.pcont,
                     Per_cont                 = pcont,
                     Var_pimp                 = var.pimp,
                     Perm_imp                 = pimp,
                     Training_AUC,
                     Number_background_points,
                     Logistic_threshold,
                     Omission_rate)

      ## Remove path gunk, and species
      d$Species     = NULL
      d$searchTaxon = gsub("/full/maxent_fitted.rds", "", d$searchTaxon)
      return(d)

    }) %>%

    ## Finally, bind all the rows together
    bind_rows

  ## Now create a list of the '10th percentile training presence Logistic threshold'. This is used in step 8 to threshold
  ## the maps to just areas above the threshold.
  message ("MAXENT.RESULTS columns")
  message (paste (colnames (MAXENT.RESULTS)))
  message (paste (nrow (MAXENT.RESULTS)))
  summary(MAXENT.RESULTS["Logistic_threshold"])
  percent.10.log <- as.list(MAXENT.RESULTS["Logistic_threshold"]) %>%
    .[["Logistic_threshold"]]

  ## Create a list of the omission files - again, don't do this for all the files, just the intersection
  omission.tables = lapply (map_spp_list, FUN = function (x) {paste(results_dir , x, "full/species_omission.csv", sep="/")})
  message (head (omission.tables))

  ## Only process the existing files
  om.exists = lapply (omission.tables, FUN = function (x) {file.exists (x)}) %>% unlist()
  # om.exists = unlist(om.exists)
  omission.tables = omission.tables[om.exists]
  message(head(omission.tables))

  ## Get the maxium TSS value using the omission data : use _training_ ommission data only
  Max_tss <- sapply(omission.tables, function(file) {

    ## For eachg species, read in the training data
    d <- read.csv(file)
    i <- which.min(d$Training.omission + d$Fractional.area)

    c(Max_tss = 1 - min(d$Training.omission + d$Fractional.area),
      thr     = d$Corresponding.logistic.value[i])

  })

  ## Add a species variable to the TSS results, so we can subset to just the species analysed
  Max_tss        = as.data.frame(Max_tss)
  MAXENT.RESULTS = cbind(MAXENT.RESULTS, Max_tss)
  summary(MAXENT.RESULTS$Max_tss)
  summary(MAXENT.RESULTS$Omission_rate)

  ## This is a summary of maxent output for current conditions
  ## All species should have AUC > 0.7
  dim(MAXENT.RESULTS)
  head(MAXENT.RESULTS, 20)[1:5]

  ## Now check the match between the species list, and the results list.
  length(intersect(map_spp_list, MAXENT.RESULTS$searchTaxon))
  MAXENT.RESULTS  =  MAXENT.RESULTS[MAXENT.RESULTS$searchTaxon %in% map_spp_list , ]
  map_spp         = unique(MAXENT.RESULTS$searchTaxon)
  length(map_spp);setdiff(sort(map_spp_list), sort(map_spp))

  ## Then make a list of all the directories containing the individual GCM rasters. This is used for combining the rasters
  SDM.RESULTS.DIR <- map_spp %>%

    ## Pipe the list into lapply
    lapply(function(species) {

      ## Create the character string...
      m <-   sprintf('%s/%s/full/', results_dir, species)                ## path.backwards.sel
      m

    }) %>%

    ## Bind the list together
    c() %>% unlist()

  length(SDM.RESULTS.DIR)

  ## Change the species column
  MAXENT.RESULTS$searchTaxon = gsub("_", " ", MAXENT.RESULTS$searchTaxon)
  MAXENT.RESULTS$results_dir = SDM.RESULTS.DIR


  ##
  if(save_data == "TRUE") {

    ## If saving, save
    saveRDS(MAXENT.RESULTS, paste0(DATA_path, 'MAXENT_RESULTS_', save_run, '.rds'))

  } else {

    ## Or return to the global environment
    message(' skip file saving, not many species analysed')
    return(MAXENT.RESULTS)

  }

}





## Project maxent files
#' @export
project_maxent_grids_mess = function(shp_path, aus_shp, world_shp, scen_list,
                                     species_list, maxent_path, climate_path, static_path,
                                     grid_names, time_slice, current_grids, create_mess, nclust) {

  ## Read in the Australian shapefile at the top
  aus_poly = readRDS(paste0(shp_path, aus_shp)) %>%
    spTransform(ALB.CONICAL)

  world_poly = readRDS(paste0(shp_path, world_shp)) %>%
    spTransform(CRS.WGS.84)

  ## First, run a loop over each scenario:
  lapply(scen_list, function(x) {

    ## Create a raster stack for each of the 6 GCMs, not for each species
    ## They need to have exactly the same extent.
    ## Could stack all the rasters, or, keep them separate
    s <- stack(c(sprintf('%s/20%s/%s/%s%s.tif', climate_path, time_slice, x, x, 1:19)))
    identical(projection(s), projection(aus_poly))

    ## Rename both the current and future environmental stack...
    ## critically important that the order of the name
    ## So this needs to include all the predicor names :: climate, soil, etc

    ## Note this step is only needed if the current grids used in the their original form, rather than being renamed
    names(s) <- names(current_grids) <- grid_names

    ## Divide the 11 temperature rasters by 10: NA values are the ocean
    message('First, divide the raster stack for ', x, ' by 10 ')
    for(i in 1:11) {
      ## Simple loop
      message(i)
      s[[i]] <- s[[ i]]/10

    }

    ## Then apply each GCM to each species.
    ## First, check if the maxent model exists
    ## Then apply each GCM to each species

    ## Define function to then send to one or multiple cores
    maxent_predict_fun <- function(species) {

      ##
      save_name = gsub(' ', '_', species)

      ## Why is this checing at the wrong level?
      if(file.exists(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, species))) {
        message('Then run maxent projections for ', species, ' under ', x, ' scenario')

        ## Then, check if the species projection has already been run...
        if(!file.exists(sprintf('%s/%s/full/%s_future_not_novel_%s.tif', maxent_path, species, species, x))) {

          ## Now read in the SDM model, calibrated on current conditions
          ## if it was run with backwards selection, just use the full model
          if (grepl("back", maxent_path)) {

            message('Read in the BS model')
            m   <- readRDS(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, save_name))

          } else {
            ## Otherwise, index the full model
            message('Read in the full model')
            m   <- readRDS(sprintf('%s/%s/full/maxent_fitted.rds', maxent_path, save_name))$me_full
          }

          ## Read in species with data and occurrence files
          swd <- as.data.frame(readRDS(sprintf('%s%s/swd.rds',    maxent_path, species, species)))
          occ <- readRDS(sprintf('%s%s/%s_occ.rds', maxent_path, species, species)) %>%
            spTransform(ALB.CONICAL)

          ## Create a file path for the current raster prediction
          f_current  <- sprintf('%s%s/full/%s_current.tif', maxent_path, species, species)

          ## If the current raster prediction has not been run, run it
          if(!file.exists(f_current)) {

            ## Report which prediction is in progress :: m$me_full, m$me_full@presence
            message('Running current prediction for ', species)

            pred.current <- rmaxent::project(
              m, current_grids[[colnames(m@presence)]])$prediction_logistic
            writeRaster(pred.current, f_current, overwrite = TRUE)

          } else {
            message('Use existing prediction for ', species)
            pred.current = raster(sprintf('%s/%s/full/%s_current.tif',
                                          maxent_path, species, species))
          }

          ## Report current mess map in progress
          ## Could work out how to the static mess once, before looping through scenarios
          MESS_dir = sprintf('%s%s/full/%s',
                             maxent_path, species, 'MESS_output')
          f_mess_current = sprintf('%s/%s%s.tif', MESS_dir, species, "_current_mess")


          ## If the current novel layer doesn't exist, create it
          if(!file.exists(sprintf('%s/%s%s.tif', MESS_dir, species, "_current_novel")))  {

            ## Set the names of the rasters to match the occ data, and subset both
            sdm_vars             = names(m@presence)
            current_grids        = subset(current_grids, sdm_vars)
            swd                  = swd [,sdm_vars]
            identical(names(swd), names(current_grids))

            ## Create a map of novel environments for current conditions.
            ## This similarity function only uses variables (e.g. n bioclim), not features
            message('Run similarity function for current condtions for ', species)
            mess_current  <- similarity(current_grids, swd, full = TRUE)
            novel_current <- mess_current$similarity_min < 0  ##   All novel environments are < 0
            novel_current[novel_current==0] <- NA             ##   0 values are NA


          } else {
            ## Otherwise, read in the current novel layer
            message(species, ' Current similarity analysis already run')
            novel_current = raster(sprintf('%s/%s%s.tif', MESS_dir, species, "_current_novel"))
          }

          ## Write out the current mess maps -
          ## create a new folder for the mess output - we are going to print it to the maps
          if(!dir.exists(MESS_dir)) {
            message('Creating MESS directory for ', species)
            dir.create(MESS_dir)

            ## Create a PNG file of MESS maps for each maxent variable
            ## raster_list  = unstack(mess_current$similarity) :: list of environmental rasters
            ## raster_names = names(mess_current$similarity)   :: names of the rasters
            message('Creating mess maps of each current environmental predictor for ', species)
            mapply(function(raster, raster_name) {

              ## Create a level plot of MESS output for each predictor variable, for each species
              p <- levelplot(raster, margin = FALSE, scales = list(draw = FALSE),
                             at = seq(minValue(raster), maxValue(raster), len = 100),
                             colorkey = list(height = 0.6),
                             main = gsub('_', ' ', sprintf('Current_mess_for_%s (%s)', raster_name, species))) +

                latticeExtra::layer(sp.polygons(aus_poly), data = list(aus_poly = aus_poly))  ## need list() for polygon

              p <- diverge0(p, 'RdBu')
              f <- sprintf('%s/%s%s%s.png', MESS_dir, species, "_current_mess_", raster_name)

              png(f, 8, 8, units = 'in', res = 300, type = 'cairo')
              print(p)
              dev.off()

            }, unstack(mess_current$similarity), names(mess_current$similarity))

          } else {
            message(species, ' MESS directory already created')
          }

          ## Write the raster of novel environments to the MESS sub-directory
          if(!file.exists(sprintf('%s/%s%s.tif', MESS_dir, species, "_current_novel")))  {

            message('Writing currently novel environments to file for ', species)
            writeRaster(novel_current, sprintf('%s/%s%s.tif', MESS_dir, species, "_current_novel"),
                        overwrite = TRUE)

          } else {
            message(species, 'Current MESS file already saved')
          }

          ## Now mask out novel environments
          ## is.na(novel_current) is a binary layer showing
          ## not novel [=1] vs novel [=0],
          ## so multiplying this with hs_current will mask out novel
          hs_current_not_novel <- pred.current * is.na(novel_current)

          ## Plot the un-novel environments
          plot(hs_current_not_novel, main = species)

          ## Write out not-novel raster :: this can go to the main directory
          message('Writing currently un-novel environments to file for ', species)
          writeRaster(hs_current_not_novel, sprintf('%s%s/full/%s%s.tif', maxent_path,
                                                    species, species, "_current_not_novel"),
                      overwrite = TRUE)

          ## Create file path for future raster doesn't exist, create it
          f_future <- sprintf('%s/%s/full/%s_future_not_novel_%s.tif',
                              maxent_path, species, species, x)

          if(!file.exists(f_future)) {

            ## Report which prediction is in progress
            message('Running future maxent prediction for ', species, ' under ', x)

            ## Create the future raster
            pred.future <- rmaxent::project(
              m, s[[colnames(m@presence)]])$prediction_logistic
            writeRaster(pred.future, f_future, overwrite = TRUE)

            ## Report future mess map in progress
            if(create_mess == "TRUE") {
              message('Running future mess map for ', species, ' under ', x)

              ## Set the names of the rasters to match the occ data, and subset both
              ## Watch the creation of objects in each run
              sdm_vars             = names(m@presence)
              future_grids         = s
              future_grids         = subset(future_grids, sdm_vars)
              swd                  = swd [,sdm_vars]
              identical(names(swd), names(future_grids))

              ## Create a map of novel environments for future conditions
              ## This similarity function only uses variables (e.g. n bioclim), not features.
              ## We don't need to repeat the static layer MESS each time - just
              ## include their results here
              mess_future  <- similarity(future_grids, swd, full = TRUE)
              novel_future <- mess_future$similarity_min < 0  ##   All novel environments are < 0
              novel_future[novel_future==0] <- NA             ##   0 values are NA


              ## Write out the future mess maps, for all variables
              writeRaster(mess_future$similarity_min, sprintf('%s/%s%s%s.tif', MESS_dir, species, "_future_mess_", x),
                          overwrite = TRUE)


              ## Create a PNG file of all the future MESS output:
              ## raster_list  = unstack(mess_current$similarity) :: list of environmental rasters
              ## raster_names = names(mess_current$similarity)   :: names of the rasters
              message('Creating mess maps of each future environmental predictor for ', species, ' under scenario ', x)
              mapply(function(raster, raster_name) {

                p <- levelplot(raster, margin = FALSE, scales = list(draw = FALSE),
                               at = seq(minValue(raster), maxValue(raster), len = 100),
                               colorkey = list(height = 0.6),
                               main = gsub('_', ' ', sprintf('Future_mess_for_%s_%s (%s)', raster_name, x, species, x))) +

                  latticeExtra::layer(sp.polygons(aus_poly), data = list(aus_poly = aus_poly))

                p <- diverge0(p, 'RdBu')
                f <- sprintf('%s/%s%s%s%s%s.png', MESS_dir, species, "_future_mess_", raster_name, "_", x)

                png(f, 8, 8, units = 'in', res = 300, type = 'cairo')
                print(p)
                dev.off()

              }, unstack(mess_future$similarity), names(mess_future$similarity))

              ## Write the raster of novel environments to the MESS maps sub-directory
              message('Writing future novel environments to file for ',    species, ' under scenario ', x)
              writeRaster(novel_future, sprintf('%s/%s%s%s.tif', MESS_dir, species, "_future_novel_",  x),
                          overwrite = TRUE)


              ## mask out future novel environments
              ## is.na(novel_future) is a binary layer showing
              ## not novel [=1] vs novel [=0],
              ## so multiplying this with hs_future will mask out novel
              hs_future_not_novel <- pred.future * is.na(novel_future)

              ## This layer of future un-novel environments can be used
              ## for the next algorithm step, where we combine the models
              plot(pred.future);plot(hs_future_not_novel)
              summary(pred.future,         maxsamp = 100000);
              summary(hs_future_not_novel, maxsamp = 100000)

              ## Write out not-novel raster
              ## Try to set up loops so different cores aren't accessing the same files............
              message('Writing un-novel environments to file under ', x, ' scenario for ', species)
              writeRaster(hs_future_not_novel, sprintf('%s%s/full/%s%s%s.tif', maxent_path,
                                                       species, species, "_future_not_novel_", x),
                          overwrite = TRUE)
            } else {
              message('Dont run future MESS maps for ', species, ' under scenario ',  x )
            }


            ## If we're on windows, use the GDAL .bat file
            novel_current_poly <- polygonizer_windows(sprintf('%s/%s%s.tif',   MESS_dir, species, "_current_novel"))
            novel_future_poly  <- polygonizer_windows(sprintf('%s/%s%s%s.tif', MESS_dir, species, "_future_novel_", x))


            ## Create the MESS path and save shapefiles
            MESS_shp_path   = sprintf('%s%s/full/%s',
                                      maxent_path, species, 'MESS_output')

            ## Check if the current MESS shapefile exists?
            novel_current_shp <- sprintf('%s/%s%s.shp',   MESS_dir, species, "_current_novel_polygon")
            if(!file.exists(novel_current_shp)) {

              ## Re-project the shapefiles
              novel_current_poly = novel_current_poly %>%
                spTransform(ALB.CONICAL)

              ## Now save the novel areas as shapefiles
              ## There is a problem with accessing the files at the same time
              message('Saving current MESS maps to polygons for ', species)
              writeOGR(obj    = novel_current_poly,
                       dsn    = sprintf('%s',  MESS_shp_path),
                       layer  = paste0(species, "_current_novel_polygon"),
                       driver = "ESRI Shapefile", overwrite_layer = TRUE)
            } else {
              message('Current MESS maps already saved to polygons for ', species)
            }

            ## Check if the future MESS shapefile exists?
            novel_future_shp <- sprintf('%s/%s%s%s.shp',   MESS_dir, species, "_future_novel_polygon_", x)
            if(!file.exists(novel_future_shp)) {

              novel_future_poly = novel_future_poly %>%
                spTransform(ALB.CONICAL)

              message('Saving current MESS maps to polygons for ', species)
              writeOGR(obj    = novel_future_poly,
                       dsn    = sprintf('%s',  MESS_shp_path),
                       layer  = paste0(species, "_future_novel_polygon_", x),
                       driver = "ESRI Shapefile", overwrite_layer = TRUE)
            } else {
              message('Future MESS maps already saved to polygons for ', species)
            }

            ## Create a SpatialLines object that indicates novel areas (this will be overlaid)
            ## Below, we create a dummy polygon as the first list element (which is the extent
            ## of the raster, expanded by 10%), to plot on panel 1). 50 = approx 50 lines across the polygon
            message('Creating polygon list under ', x, ' scenario for ', species)

            ## Cast the objects into the sf class so we avoid issues with wrong methods being called in hatch()
            novel_hatch <- list(as(extent(pred.current)*1.1, 'SpatialLines'),
                                hatch(novel_current_poly, 50),
                                hatch(novel_future_poly, 50))


            ## Now create a panel of PNG files for maxent projections and MESS maps
            ## All the projections and extents need to match
            empty_ras <- init(pred.current, function(x) NA)
            projection(novel_current_poly);projection(occ);projection(empty_ras);projection(poly)
            projection(pred.current);projection(pred.future)
            identical(extent(pred.current), extent(pred.future))

            ## Assign the scenario name (to use in the plot below)
            scen_name = eval(parse(text = sprintf('gcms.%s$GCM[gcms.%s$id == x]', time_slice, time_slice)))


            ## Use the 'levelplot' function to make a multipanel output: occurrence points, current raster and future raster
            current_mess_png = sprintf('%s/%s/full/%s_%s.png', maxent_path, species, species, "mess_panel")
            if(!file.exists(current_mess_png)) {

              ## Create level plot of current conditions including MESS
              message('Create current MESS panel maps for ', species)

              png(sprintf('%s/%s/full/%s_%s.png', maxent_path, species, species, "mess_panel"),
                  11, 4, units = 'in', res = 300)

              print(levelplot(stack(empty_ras,
                                    pred.current, quick = TRUE), margin = FALSE,

                              ## Create a colour scheme using colbrewer: 100 is to make it continuos
                              ## Also, make it a one-directional colour scheme
                              scales      = list(draw = FALSE),
                              at = seq(0, 1, length = 100),
                              col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),

                              ## Give each plot a name: the third panel is the GCM
                              names.attr = c('Australian records', 'Current'),
                              colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                              main       = list(gsub('_', ' ', species), font = 4, cex = 2)) +

                      ## Plot the Aus shapefile with the occurrence points for reference
                      ## Can the current layer be plotted on it's own?
                      ## Add the novel maps as vectors.
                      latticeExtra::layer(sp.polygons(aus_poly), data = list(aus_poly = aus_poly)) +
                      latticeExtra::layer(sp.points(occ, pch = 19, cex = 0.15,
                                                    col = c('red', 'transparent', 'transparent')[panel.number()]),
                                          data = list(occ = occ)) +
                      latticeExtra::layer(sp.lines(h[[panel.number()]]), data = list(h = novel_hatch)))
              dev.off()

            } else {
              message('Current MESS panel maps already created for ', species)
            }

          }


          ## Save the global records to PNG :: try to code the colors for ALA/GBIF/INVENTORY
          occ.world <- readRDS(sprintf('%s/%s/%s_occ.rds', maxent_path, species, species)) %>%
            spTransform(CRS.WGS.84)

          ## If the global map of occurrence points hasn't been created, create it
          global_occ_map = sprintf('%s/%s/full/%s_%s.png', maxent_path, species, species, "global_occ_records")
          if(!file.exists(global_occ_map)) {

            message('writing map of global records for ', species)
            png(sprintf('%s/%s/full/%s_%s.png', maxent_path, species, species, "global_occ_records"),
                16, 10, units = 'in', res = 500)

            ## Add land
            plot(world_poly, #add = TRUE,
                 lwd = 0.01, asp = 1, col = 'grey', bg = 'sky blue')

            ## Add points
            points(subset(occ.world, SOURCE == "GBIF"),
                   pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2,
                   xlab = "", ylab = "", asp = 1,
                   col = "orange",
                   legend(7,4.3, unique(occ.world$SOURCE), col = "orange", pch = 1))

            points(subset(occ.world, SOURCE == "ALA"),
                   pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2,
                   xlab = "", ylab = "", asp = 1,
                   col = "blue",
                   legend(7,4.3, unique(occ.world$SOURCE), col = "blue", pch = 1))

            points(subset(occ.world, SOURCE == "INVENTORY"),
                   pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2,
                   xlab = "", ylab = "", asp = 1,
                   col = "red",
                   legend(7,4.3, unique(occ.world$SOURCE), col = "red", pch = 1))

            title(main = list(paste0(gsub('_', ' ', species), ' global SDM records'), font = 4, cex = 2),
                  cex.main = 4,   font.main = 4, col.main = "black")

            dev.off()

          } else {
            message('Global occurrence maps already created for ', species)
          }


          ## Create level plot of scenario x, including MESS
          future_mess_png = sprintf('%s/%s/full/%s_%s.png', maxent_path, species, species, x)

          if(!file.exists(future_mess_png)) {
            png(sprintf('%s/%s/full/%s_%s.png', maxent_path, species, species, x),
                11, 4, units = 'in', res = 300)

            ## Create a panel of the Australian occurrences, the current layer and the future layer
            print(levelplot(stack(empty_ras,
                                  pred.current,
                                  pred.future, quick = TRUE), margin = FALSE,

                            ## Create a colour scheme using colbrewer: 100 is to make it continuos
                            ## Also, make it a one-directional colour scheme
                            scales      = list(draw = FALSE),
                            at = seq(0, 1, length = 100),
                            col.regions = colorRampPalette(rev(brewer.pal(11, 'Spectral'))),

                            ## Give each plot a name: the third panel is the GCM
                            names.attr = c('Australian records', 'Current', sprintf('%s, 20%s, RCP8.5', scen_name, time_slice)),
                            colorkey   = list(height = 0.5, width = 3), xlab = '', ylab = '',
                            main       = list(gsub('_', ' ', species), font = 4, cex = 2)) +

                    ## Plot the Aus shapefile with the occurrence points for reference
                    ## Can the current layer be plotted on it's own?
                    ## Add the novel maps as vectors.
                    latticeExtra::layer(sp.polygons(aus_poly), data = list(aus_poly = aus_poly)) +
                    latticeExtra::layer(sp.points(occ, pch = 19, cex = 0.15,
                                                  col = c('red', 'transparent', 'transparent')[panel.number()]),
                                        data = list(occ = occ)) +
                    latticeExtra::layer(sp.lines(h[[panel.number()]]), data = list(h = novel_hatch)))
            dev.off()

          } else {
            message('Future MESS maps already created for ', species)
          }

        } else {
          message(species, ' ', x, ' skipped - prediction already run')
        }

      } else {
        message(species, ' ', x, ' skipped - SDM not yet run')
      }

    }

    ## Check this is the best way to run parallel...............
    if (nclust==1) {

      lapply(species_list, maxent_predict_fun)

    } else {

      ## Export all objects from the function call
      # shp_path, aus_shp, world_shp, scen_list,
      # species_list, maxent_path, climate_path,
      # grid_names, time_slice, current_grids, create_mess, nclust

      cl <- makeCluster(nclust)
      clusterExport(cl, c(
        'shp_path',    'aus_shp',       'world_shp',   'ALB.CONICAL',  'CRS.WGS.84',
        'scen_list',   'species_list',  'maxent_path', 'climate_path', 'grid_names',
        'time_slice',  'current_grids', 'create_mess', 'hatch', 'x',
        'polygonizer', 'nclust', 'diverge0'),  envir = environment())

      clusterEvalQ(cl, {

        library(rmaxent)
        library(sp)
        library(raster)
        library(rasterVis)
        library(latticeExtra)
        library(magrittr)

      })

      message('Running project_maxent_grids_mess for ', length(species_list),
              ' species on ', nclust, ' cores for GCM ', x)
      parLapply(cl, species_list, maxent_predict_fun)
    }

  })

}





## Simplify the maxent models ----
#' @export
local_simplify = function (occ, bg, path, species_column = "species", response_curves = TRUE,
                           logistic_format = TRUE, type = "PI", cor_thr, pct_thr, k_thr,
                           features = "lpq", replicates = 1, quiet = TRUE)
{
  if (!species_column %in% colnames(occ))
    stop(species_column, " is not a column of `occ`", call. = FALSE)
  if (!species_column %in% colnames(bg))
    stop(species_column, " is not a column of `bg`", call. = FALSE)
  if (missing(path)) {
    save <- FALSE
    path <- tempdir()

  }

  else save <- TRUE
  features <- unlist(strsplit(gsub("\\s", "", features), ""))
  if (length(setdiff(features, c("l", "p", "q", "h", "t"))) >
      1)
    stop("features must be a vector of one or more of ',\n         'l', 'p', 'q', 'h', and 't'.")
  off <- setdiff(c("l", "p", "q", "t", "h"), features)
  if (length(off) > 0) {
    off <- c(l = "linear=FALSE", p = "product=FALSE", q = "quadratic=FALSE",
             t = "threshold=FALSE", h = "hinge=FALSE")[off]

  }
  off <- unname(off)
  occ_by_species <- split(occ, occ[[species_column]])
  bg_by_species <- split(bg, bg[[species_column]])
  if (!identical(sort(names(occ_by_species)), sort(names(bg_by_species)))) {
    stop("The same set of species names must exist in occ and bg")
  }
  type <- switch(type, PI = "permutation.importance", PC = "contribution",
                 stop("type must be either \"PI\" or \"PC\".", call. = FALSE))
  args <- off
  if (replicates > 1)
    args <- c(args, paste0("replicates=", replicates))
  if (isTRUE(response_curves))
    args <- c(args, "responsecurves=TRUE")
  if (isTRUE(logistic_format))
    args <- c(args, "outputformat=logistic")
  f <- function(name) {
    if (!quiet)
      message("\n\nDoing ", name)
    name_ <- gsub(" ", "_", name)
    swd <- rbind(occ_by_species[[name]], bg_by_species[[name]])
    swd <- swd[, -match(species_column, names(swd))]
    if (ncol(swd) < k_thr)
      stop("Initial number of variables < k_thr", call. = FALSE)
    pa <- rep(1:0, c(nrow(occ_by_species[[name]]), nrow(bg_by_species[[name]])))
    vc <- usdm::vifcor(swd, maxobservations = nrow(swd),
                       th = cor_thr)
    vif <- methods::slot(vc, "results")
    k <- nrow(vif)
    exclude <- methods::slot(vc, "excluded")
    if (!isTRUE(quiet) & length(exclude) > 0) {
      message("Dropped due to collinearity: ", paste0(exclude,
                                                      collapse = ", "))
    }
    if (k < k_thr)
      stop(sprintf("Number of uncorrelated variables (%s) < k_thr (%s). %s",
                   k, k_thr, "Reduce k_thr, increase cor_thr, or find alternative predictors."),
           call. = FALSE)
    swd_uncor <- swd[, as.character(vif$Variables)]
    d <- file.path(path, name_, if (replicates > 1)
      "xval"
      else "full")
    m <- dismo::maxent(swd_uncor, pa, args = args, path = d)
    if (isTRUE(save))
      saveRDS(m, file.path(d, "maxent_fitted.rds"))
    pct <- m@results[grep(type, rownames(m@results)), ,
                     drop = FALSE]
    pct <- pct[, ncol(pct)]
    pct <- sort(pct)
    names(pct) <- sub(paste0("\\.", type), "", names(pct))
    if (min(pct) >= pct_thr || length(pct) == k_thr) {
      if (replicates > 1) {
        d <- file.path(path, name_, "full")
        m <- dismo::maxent(swd_uncor, pa, args = grep("replicates",
                                                      args, value = TRUE, invert = TRUE), path = d)
      }
      if (isTRUE(save)) {
        saveRDS(m, file.path(d, "maxent_fitted.rds"))
      }
      return(m)
    }
    while (min(pct) < pct_thr && length(pct) > k_thr) {
      candidates <- vif[vif$Variables %in% names(pct)[pct ==
                                                        pct[1]], ]
      drop <- as.character(candidates$Variables[which.max(candidates$VIF)])
      message("Dropping ", drop)
      swd_uncor <- swd_uncor[, -match(drop, colnames(swd_uncor))]
      if (!quiet)
        message(sprintf("%s variables: %s", ncol(swd_uncor),
                        paste0(colnames(swd_uncor), collapse = ", ")))
      m <- dismo::maxent(swd_uncor, pa, args = args, path = d)
      pct <- m@results[grep(type, rownames(m@results)),
                       , drop = FALSE]
      pct <- pct[, ncol(pct)]
      pct <- sort(pct)
      names(pct) <- sub(paste0("\\.", type), "", names(pct))
    }
    if (replicates > 1) {
      d <- file.path(path, name_, "full")
      m <- dismo::maxent(swd_uncor, pa, args = grep("replicates",
                                                    args, value = TRUE, invert = TRUE), path = d)
    }
    if (isTRUE(save)) {
      saveRDS(m, file.path(d, "maxent_fitted.rds"))
    }
    return(m)
  }
  lapply(names(occ_by_species), f)
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



