#########################################################################################################################
###################################  FUNCTIONS FOR PREPARING SDM DATA ---- ##############################################
#########################################################################################################################


## Below are the functions used to estiamte species ranges and prepare occurrence data for SDM analysis


## GBIF download ----


#' Download species occurrence files from GBIF
#'
#' This function downloads species occurrence files from GBIF (https://www.gbif.org/).
#' It assumes that the species list supplied is taxonomically correct (haha!).
#' It downloads the species to the specified folders without returning anything
#' @param species_list   Character vector - List of species binomials to download
#' @param download_path  Character string - File path for species downloads
#' @param download_limit Numeric - How many records can be downloaded at one time? Set by server
#' @export
download_GBIF_all_species = function(species_list, download_path, download_limit) {

  ## create variables
  GBIF.download.limit = download_limit

  ## for every species in the list
  for(sp.n in species_list){

    ## First, check if the f*&%$*# file exists
    file_name = paste0(download_path, sp.n, "_GBIF_records.RData")

    ## If it's already downloaded, skip
    if (file.exists (file_name)) {

      print(paste ("file exists for species", sp.n, "skipping"))
      next

    }

    ## create a dummy file
    dummy = data.frame()
    save (dummy, file = file_name)

    ## Then check the spelling...incorrect nomenclature will return NULL result
    if (is.null(occ_search(scientificName = sp.n, limit = 1)$meta$count) == TRUE) {

      ## now append the species which had incorrect nomenclature to the skipped list
      ## this is slow, but it works for now
      print (paste ("Possible incorrect nomenclature", sp.n, "skipping"))
      nomenclature = paste ("Possible incorrect nomenclature |", sp.n)
      #skip.spp.list <- c(skip.spp.list, nomenclature)
      next

    }

    ## Skip species with no records
    if (occ_search(scientificName = sp.n)$meta$count <= 2) {

      ## now append the species which had no records to the skipped list
      print (paste ("No GBIF records for", sp.n, "skipping"))
      records = paste ("No GBIF records |", sp.n)
      #skip.spp.list <- c(skip.spp.list, records)
      next

    }

    ## Check how many records there are, and skip if there are over 200k
    if (occ_search(scientificName = sp.n, limit = 1)$meta$count > GBIF.download.limit) {

      ## now append the species which had > 200k records to the skipped list
      print (paste ("Number of records > max for GBIF download via R (200,000)", sp.n))
      max =  paste ("Number of records > 200,000 |", sp.n)

    } else {

      ## Download ALL records from GBIF
      message("Downloading GBIF records for ", sp.n, " using rgbif :: occ_data")
      key <- name_backbone(name = sp.n, rank = 'species')$usageKey

      GBIF <- occ_data(taxonKey = key, limit = GBIF.download.limit)
      GBIF <- as.data.frame(GBIF$data)

      cat("Synonyms returned for :: ",  sp.n, unique(GBIF$scientificName), sep = "\n")
      cat("Names returned for :: ",     sp.n, unique(GBIF$name),           sep = "\n")
      cat("Takonkeys returned for :: ", sp.n, unique(GBIF$taxonKey),       sep = "\n")

      ## Could also only use the key searched, but that could knock out a lot of species
      message(dim(GBIF[1]), " Records returned for ", sp.n)

      ## Save records to .Rdata file
      save(GBIF, file = file_name)

    }
  }
}





## ALA download ----


#' Download species occurrence files from the Atlas of Living Australia (ALA)
#'
#' This function downloads species occurrence files from ALA (https://www.ala.org.au/).
#' It assumes that the species list supplied is taxonomically correct.
#' It downloads the species without returning anything
#'
#' @param species_list   Character vector - List of species binomials to download
#' @param download_path  Character string - File path for species downloads
#' @param download_limit Numeric - How many records can be downloaded at one time? Set by server
#' @export
download_ALA_all_species = function (species_list, your_email, download_path, download_limit) {

  ## create variables
  download_limit  = 200000

  ## for every species in the list
  for(sp.n in species_list){

    ## Get the ID?
    lsid <- ALA4R::specieslist(sp.n)$taxonConceptLsid

    ## First, check if the f*&%$*# file exists
    file_name = paste0(download_path, sp.n, "_ALA_records.RData")

    ## If it's already downloaded, skip
    if (file.exists (file_name)) {

      print (paste ("file exists for species", sp.n, "skipping"))
      next

    }
    ## create a dummy file
    dummy = data.frame()
    save (dummy, file = file_name)

    ## Then check the spelling...incorrect nomenclature will return NULL result
    if (is.null(ALA4R::occurrences(taxon = paste('taxon_name:\"', sp.n, '\"', sep = ""),
                                   download_reason_id = 7, email = your_email)$data) == TRUE) {

      ## Now, append the species which had incorrect nomenclature to the skipped list
      print (paste ("Possible incorrect nomenclature", sp.n, "skipping"))
      nomenclature = paste ("Possible incorrect nomenclature |", sp.n)
      #skip.spp.list <- c(skip.spp.list, nomenclature)
      next

    }

    ## Skip species with no records
    if (nrow(ALA4R::occurrences(taxon = paste('taxon_name:\"', sp.n, '\"',sep=""),
                                download_reason_id = 7, email = your_email)$data) <= 2) {

      ## now append the species which had no records to the skipped list
      print (paste ("No ALA records for", sp.n, "skipping"))
      records = paste ("No ALA records |", sp.n)
      #skip.spp.list <- c(skip.spp.list, records)
      next

    }

    ## Download ALL records from ALA ::
    message("Downloading ALA records for ", sp.n, " using ALA4R :: occurrences")
    ALA = ALA4R::occurrences(taxon = paste('taxon_name:\"', sp.n, '\"',sep=""), download_reason_id = 7, email = your_email)
    ALA = ALA[["data"]]

    cat("Synonyms returned for :: ", sp.n, unique(ALA$scientificName), sep="\n")
    message(dim(ALA[1]), " Records returned for ", sp.n)

    ## Save records to .Rdata file
    save(ALA, file = file_name)

  }

}





## Combining ALA records -----


#' Combine all ala records into one file
#'
#' This function combines all occurrence files from ALA into.
#' There are slight differences between ALA and GBIF, so separate functions are useful.
#' It assumes that all the files come from the previous downloads function.
#' Although you can download all the records for in one go, this is better for
#' Doing small runs of species, or where you want to re-run them constantly
#' @param species_list       Character Vector - List of species already downloaded
#' @param records_path       File path for downloaded species
#' @param records_extension  Which R file type? RDS or RDA
#' @param record_type        Adds a column to the data frame for the data source, EG ALA
#' @param keep_cols          The columns we want to keep - a character list created by you
#' @param world_raster       An Raster file of the enviro conditions used (assumed to be global)
#' @export
combine_ala_records = function(species_list, records_path, records_extension,
                               record_type, keep_cols, world_raster) {

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
        dplyr::select(one_of(keep_cols))

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
    z   = raster::extract(world_raster, xy)

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


#' This function combines all occurrence files from GBIF into one table.
#' There are slight differences between ALA and GBIF, so separate functions are useful.
#' It assumes that all the files come from the previous GBIF downloads function.
#' Although you can download all the records for in one go, this is better for
#' Doing small runs of species, or where you want to re-run them constantly
#' @param species_list       List of species already downloaded
#' @param records_path       File path for downloaded species
#' @param records_extension  Which R file type? RDS or RDA
#' @param record_type        Adds a column to the data frame for the data source, EG ALA
#' @param keep_cols          The columns we want to keep - a character list created by you
#' @param world_raster       An Raster file of the enviro conditions used (assumed to be global)
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
    xy <- cellFromXY(world_raster, GBIF.CLEAN[c("lon", "lat")]) %>%

      ## get the unique raster cells
      unique %>%

      ## Get coordinates of the center of raster cells for a row, column, or cell number of WORLDCLIM raster
      xyFromCell(world_raster, .) %>%
      na.omit()

    ## For some reason, we need to convert the xy coords to a spatial points data frame, in order to avoid this error:
    ## 'NAs introduced by coercion to integer range'
    xy <- SpatialPointsDataFrame(coords = xy, data = as.data.frame(xy),
                                 proj4string = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))

    ## Now extract the temperature values for the unique 1km centroids which contain GBIF data
    message('Removing GBIF points in the ocean for ', length(species_list), ' species')
    class(xy)
    z   = raster::extract(world_raster, xy)

    ## Then track which values of Z are on land or not
    onland = z %>% is.na %>%  `!` # %>% xy[.,]  cells on land or not

    ## Finally, filter the cleaned GBIF data to only those points on land.
    ## This is achieved with the final [onland]
    LAND.POINTS = filter(GBIF.CLEAN, cellFromXY(world_raster, GBIF.CLEAN[c("lon", "lat")]) %in%
                           unique(cellFromXY(world_raster,    GBIF.CLEAN[c("lon", "lat")]))[onland])

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


#' This function combines occurrence files from ALA and GBIF into one table, and extracts enviro values.
#' It assumes that both files come from the previous GBIF/ALA combine function.
#' @param ala_df             Data frame of ALA records
#' @param gbif_df            Data frame of GBIF records
#' @param urban_df           Data frame of Urban records (only used if you have urban data, e.g. I-naturalist)
#' @param species_list       List of species analysed, used to cut the dataframe down
#' @param thin_records       Do you want to thin the records out? If so, it will be 1 record per 1km*1km grid cell
#' @param template_raster    A global R Raster used to thin records to 1 record per 1km grid cell
#' @param world_raster       An global R Raster of the enviro conditions used to extract values for all records
#' @param prj                The projection system used. Currently, needs to be WGS84
#' @param biocl_vars         The variables used - eg the standard bioclim names (https://www.worldclim.org/).
#' @param env_vars           The actual variable names (e.g. bio1 = rainfall, etc.) Only needed for worldlcim
#' @param worldclim_grids    Are you using worldclim stored as long intergers? If so, divide by 10.
#' @param save_data          Do you want to save the data frame?
#' @param data_path          The file path used for saving the data frame
#' @param save_run           A run name to append to the data frame (e.g. bat species, etc.). Useful for multiple runs.
#' @return                   Data frame of all ALA/GBIF records, with global enviro conditions for each record location (i.e. lat/lon)
#' @export

combine_records_extract = function(ala_df,
                                   gbif_df,
                                   urban_df,
                                   add_urban,
                                   species_list,
                                   template_raster,
                                   thin_records,
                                   world_raster,
                                   prj,
                                   biocl_vars,
                                   env_vars,
                                   worldclim_grids,
                                   save_data,
                                   data_path,
                                   save_run) {

  ## Get just the Common columns
  common.cols    = intersect(names(gbif_df), names(ala_df))
  GBIF.ALA.COMBO = bind_rows(gbif_df, ala_df)
  GBIF.ALA.COMBO = GBIF.ALA.COMBO %>%
    dplyr::select(one_of(common.cols))

  message(length(unique(GBIF.ALA.COMBO$searchTaxon)))
  length(unique(GBIF.ALA.COMBO$scientificName))

  ## If Urban = TRUE
  if(urban_df != 'NONE') {

    urban_cols     <- intersect(names(GBIF.ALA.COMBO), names(urban_df))
    urban_df       <- dplyr::select(urban_df, urban_cols)
    GBIF.ALA.COMBO <- bind_rows(GBIF.ALA.COMBO, urban_df)

  } else {
    message('Dont add urban data' )
  }

  ## CHECK TAXONOMY RETURNED BY ALA USING TAXONSTAND
  GBIF.ALA.MATCH = GBIF.ALA.COMBO

  ## Create points: the 'over' function seems to need geographic coordinates for this data...
  GBIF.ALA.84   = SpatialPointsDataFrame(coords      = GBIF.ALA.MATCH[c("lon", "lat")],
                                         data        = GBIF.ALA.MATCH,
                                         proj4string = prj)

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
  ## This relies on the bioclim order, it must be the same
  setnames(COMBO.RASTER, old = names(world_raster), new = env_vars)
  COMBO.RASTER <- COMBO.RASTER %>% dplyr::select(-lat.1, -lon.1)

  ## Change the raster values here: See http://worldclim.org/formats1 for description of the interger conversion.
  ## All worldclim temperature variables were multiplied by 10, so then divide by 10 to reverse it.
  if (worldclim_grids == TRUE) {

    ## Convert the worldclim grids
    message('Processing worldclim 1.0 data, divide the rasters by 10')

    COMBO.RASTER.CONVERT = as.data.table(COMBO.RASTER)
    COMBO.RASTER.CONVERT[, (env_variables [c(1:11)]) := lapply(.SD, function(x)
      x / 10 ), .SDcols = env_variables [c(1:11)]]
    COMBO.RASTER.CONVERT = as.data.frame(COMBO.RASTER.CONVERT)

  } else {

    message('Do not divide the rasters')
    COMBO.RASTER.CONVERT = COMBO.RASTER
  }

  ## Print the dataframe dimensions to screen :: format to recognise millions, hundreds of thousands, etc.
  COMBO.RASTER.CONVERT = completeFun(COMBO.RASTER.CONVERT, env_vars[1])

  message(length(unique(COMBO.RASTER.CONVERT$searchTaxon)),
          ' species processed of ', length(species_list), ' original species')

  ## save data
  if(save_data == TRUE) {

    ## save .rds file for the next session
    saveRDS(COMBO.RASTER.CONVERT, paste0(data_path, 'COMBO_RASTER_CONVERT_',  save_run, '.rds'))

  } else {
    return(COMBO.RASTER.CONVERT)
  }
  gc()
}




## Extract environmental values for urban occurrence records -----


#' This function takes a data frame from Ubran sources (e.g. I-naturalist), and extracts enviro values.
#' It assumes that the urban dataframe has these columns : species, lat, lon, Country, INVENTORY, SOURCE
#' @param urban_df           Data.frame of Urban records (only used if you have urban data, e.g. I-naturalist)
#' @param species_list       Character string - List of species analysed, used to cut the dataframe down
#' @param thin_records       Do you want to thin the records out? If so, it will be 1 record per 1km*1km grid cell
#' @param template_raster    A global R Raster used to thin records to 1 record per 1km grid cell
#' @param world_raster       An global R Raster of the enviro conditions used to extract values for all records
#' @param prj                The projection system used. Currently, needs to be WGS84
#' @param biocl_vars         The variables used - eg the standard bioclim names (https://www.worldclim.org/).
#' @param env_vars           The actual anmes of the variables (e.g. bio1 = rainfall, etc.) Only needed for worldlcim
#' @param worldclim_grids    Are you using worldclim stored as long intergers? If so, divide by 10.
#' @param save_data          Do you want to save the data frame?
#' @param data_path          The file path used for saving the data frame
#' @param save_run           Character string - run name to append to the data frame (e.g. bat species, etc.). Useful for multiple runs.
#' @return                   Data frame of all urban records, with global enviro conditions for each record location (i.e. lat/lon)
#' @export
urban_records_extract = function(urban_df,
                                 species_list,
                                 thin_records,
                                 template_raster,
                                 world_raster,
                                 prj,
                                 biocl_vars,
                                 env_vars,
                                 worldclim_grids,
                                 save_data,
                                 data_path,
                                 save_run) {

  ## Get just the species list
  URBAN.XY = urban_df[urban_df$searchTaxon %in% species_list, ]
  message('Extracting raster values for ',
          length(unique(URBAN.XY$searchTaxon)), ' urban species across ',
          length(unique(URBAN.XY$INVENTORY)),   ' Councils ')

  ## Create points: the 'over' function seems to need geographic coordinates for this data...
  URBAN.XY   = SpatialPointsDataFrame(coords      = URBAN.XY[c("lon", "lat")],
                                      data        = URBAN.XY,
                                      proj4string = prj)

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
  COMBO.RASTER <- COMBO.RASTER %>% dplyr::select(-lat.1, -lon.1)

  ## Change the raster values here: See http://worldclim.org/formats1 for description of the interger conversion.
  ## All worldclim temperature variables were multiplied by 10, so then divide by 10 to reverse it.
  if (worldclim_grids == TRUE) {

    ## Convert the worldclim grids
    message('Processing worldclim 1.0 data, divide the rasters by 10')

    COMBO.RASTER.CONVERT = as.data.table(COMBO.RASTER)
    COMBO.RASTER.CONVERT[, (env_variables [c(1:11)]) := lapply(.SD, function(x)
      x / 10 ), .SDcols = env_variables [c(1:11)]]
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
  if(save_data == TRUE) {

    ## save .rds file for the next session
    saveRDS(COMBO.RASTER.CONVERT, paste0(data_path, 'COMBO_RASTER_CONVERT_',  save_run, '.rds'))

  } else {

    return(COMBO.RASTER.CONVERT)

  }
  ## get rid of some memory
  gc()
}





## Clean coordinates of occurence records ----


#' This function takes a data frame of all species records, and flag records as institutional or spatial outliers.
#' It uses the CoordinateCleaner package https://cran.r-project.org/web/packages/CoordinateCleaner/index.html.
#' It assumes that the records data.frame is that returned by the combine_records_extract function
#' @param records            Data.frame. DF of all species records returned by the combine_records_extract function
#' @param capitals           Numeric. Remove records within an xkm radius of capital cites (see ?clean_coordinates)
#' @param centroids          Numeric. Remove records within an xkm radius around country centroids (see ?clean_coordinates)
#' @param save_run           Character string - run name to append to the data frame, useful for multiple runs.
#' @param save_data          Logical or character - do you want to save the data frame?
#' @param data_path          Character string - The file path used for saving the data frame
#' @return                   Data.frame of all urban records, with global enviro conditions for each record location (i.e. lat/lon)
#' @export
coord_clean_records = function(records,
                               capitals,
                               centroids,
                               save_run,
                               data_path,
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

                      ## remove records within xkm  of capitals
                      ## remove records within xkm of country centroids
                      capitals_rad  = capitals,
                      centroids_rad = centroids) %>%

    ## The select the relevant columns and rename
    dplyr::select(., coord_spp, CC.OBS, .val,  .equ, .zer, .cap,
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

  CLEAN.TRUE <- subset(COORD.CLEAN, coord_summary == TRUE)
  message(round(nrow(CLEAN.TRUE)/nrow(COORD.CLEAN)*100, 2), " % records retained")
  message(table(COORD.CLEAN$coord_summary))


  if(save_data == TRUE) {

    ## save .rds file for the next session
    saveRDS(COORD.CLEAN, paste0(data_path, 'COORD_CLEAN_', save_run, '.rds'))

  } else {

    message('Return the cleaned occurrence data to the global environment')   ##
    return(COORD.CLEAN)

  }

}





## Save spatial outliers to file ----


#' This function takes a data frame of all species records,
#' flags records as spatial outliers (T/F for each record in the df), and saves images of the checks for each.
#' Manual cleaning of spatial outliers is very tedious, but automated cleaning makes mistakes, so checking is handy
#' It uses the CoordinateCleaner package https://cran.r-project.org/web/packages/CoordinateCleaner/index.html.
#' It assumes that the input dfs are those returned by the coord_clean_records function
#' @param all_df             Data.frame. DF of all species records returned by the coord_clean_records function
#' @param urban_df           Data.frame of Urban records (only used if you have urban data, e.g. I-naturalist)
#' @param land_shp           R object. Shapefile of the worlds land (e.g. https://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-land/)
#' @param clean_path         Character string -  The file path used for saving the checks
#' @param spatial_mult       Numeric. The multiplier of the interquartile range (method == 'quantile', see ?cc_outl)
#' @return                   Data.frame of species records, with spatial outlier T/F flag for each record
#' @export
check_spatial_outliers = function(all_df,
                                  urban_df,
                                  land_shp,
                                  clean_path,
                                  spatial_mult,
                                  prj) {

  ## Try plotting the points which are outliers for a subset of spp and label them
  ALL.PLOT = SpatialPointsDataFrame(coords      = all_df[c("lon", "lat")],
                                    data        = all_df,
                                    proj4string = prj)

  CLEAN.TRUE = subset(all_df, coord_summary == TRUE)
  CLEAN.PLOT = SpatialPointsDataFrame(coords      = CLEAN.TRUE[c("lon", "lat")],
                                      data        = CLEAN.TRUE,
                                      proj4string = prj)

  ## Create global and australian shapefile in the local coordinate system
  LAND.84 = land_shp %>%
    spTransform(prj)
  AUS.84 = AUS %>%
    spTransform(prj)

  ## spp = spat.taxa[1]
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
    plot(LAND.84, main = paste0(nrow(subset(ALL.PLOT, coord_summary == FALSE)),
                                " Global clean_coord 'FALSE' points for ", spp),
         lwd = 0.01, asp = 1, col = 'grey', bg = 'sky blue')

    points(CLEAN.PLOT.PI,
           pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2,
           xlab = "", ylab = "", asp = 1,
           col = factor(ALL.PLOT$coord_summary))

    ## Plot true and false points for the world
    ## Black == FALSE
    ## Red   == TRUE
    plot(AUS.84, main = paste0(nrow(subset(ALL.PLOT, coord_summary == FALSE)),
                               " Global clean_coord 'FALSE' points for ", spp),
         lwd = 0.01, asp = 1, bg = 'sky blue', col = 'grey')

    points(ALL.PLOT,
           pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2,
           xlab = "", ylab = "", asp = 1,
           col = factor(ALL.PLOT$coord_summary))

    dev.off()

  }

  ## Create a tibble to supply to coordinate cleaner
  test.geo = SpatialPointsDataFrame(coords      = all_df[c("lon", "lat")],
                                    data        = all_df,
                                    proj4string = prj)

  SDM.COORDS  <- test.geo %>%

    spTransform(., prj) %>%

    as.data.frame() %>%

    dplyr::select(searchTaxon, lon, lat, CC.OBS, SOURCE) %>%

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
    dplyr::select(species) %>%
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
                         verbose = TRUE)

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
  SPAT.FLAG = join(as.data.frame(test.geo), SPAT.OUT) ## Join means the skipped spp are left out
  dim(SPAT.FLAG)

  ## Try plotting the points which are outliers for a subset of spp and label them
  ## Get the first 10 spp
  ## spp = spat.taxa[1]
  spat.taxa <- as.character(unique(SPAT.FLAG$searchTaxon))
  for (spp in spat.taxa) {

    ## Plot a subset of taxa
    SPAT.PLOT <- SPAT.FLAG %>% filter(searchTaxon == spp) %>%
      SpatialPointsDataFrame(coords      = .[c("lon", "lat")],
                             data        = .,
                             proj4string = prj)

    message("plotting occ data for ", spp, ", ",
            nrow(SPAT.PLOT ), " records")

    ## Plot true and false points for the world
    ## Black == FALSE
    ## Red   == TRUE
    message('Writing map of global coord clean records for ', spp)
    png(sprintf("%s%s_%s", clean_path, spp, "global_spatial_outlier_check.png"),
        16, 10, units = 'in', res = 500)

    par(mfrow = c(1,2))
    plot(LAND.84, main = paste0(nrow(subset(SPAT.PLOT, SPAT_OUT == FALSE)),
                                " Spatial outlier 'FALSE' points for ", spp),
         lwd = 0.01, asp = 1, col = 'grey', bg = 'sky blue')

    points(SPAT.PLOT,
           pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2,
           xlab = "", ylab = "", asp = 1,
           col = factor(SPAT.PLOT$SPAT_OUT))

    ## Plot true and false points for the world
    ## Black == FALSE
    ## Red   == TRUE
    plot(AUS.84, main = paste0("Australian points for ", spp),
         lwd = 0.01, asp = 1, bg = 'sky blue', col = 'grey')

    points(SPAT.PLOT,
           pch = ".", cex = 3.3, cex.lab = 3, cex.main = 4, cex.axis = 2,
           xlab = "", ylab = "", asp = 1,
           col = factor(SPAT.PLOT$SPAT_OUT))

    dev.off()

  }

  ## Could add the species in urban records in here
  if(urban_df == TRUE) {

    ##
    message('Combine the Spatially cleaned data with the Urban data' )
    SPAT.FLAG   <- SPAT.FLAG %>% filter(SPAT_OUT == TRUE)
    SPAT.FLAG   <- bind_rows(SPAT.FLAG, urban_df)
    SPAT.TRUE   <- SpatialPointsDataFrame(coords      = SPAT.FLAG[c("lon", "lat")],
                                          data        = SPAT.FLAG,
                                          proj4string = prj)
    message('Cleaned ', paste0(unique(SPAT.TRUE$SOURCE), sep = ' '), ' records')

  } else {
    message('Dont add urban data' )
    SPAT.TRUE <- SPAT.FLAG %>% filter(SPAT_OUT == TRUE)
  }

  return(SPAT.TRUE)

}





## Calculate niches using occurrence data ----


#' This function takes a data frame of all species records,
#' estimates geographic and environmental ranges for each, and creates a table of each.
#' It uses the AOO.computing function in the ConR package https://cran.r-project.org/web/packages/ConR/index.html
#' It assumes that the input df is that returned by the check_spatial_outliers function
#' @param coord_df           Data.frame. DF of all species records returned by the coord_clean_records function
#' @param aus_df             Data.frame of Urban records (only used if you have urban data, e.g. I-naturalist)
#' @param world_shp          .Rds object. Shapefile of the worlds land (e.g. https://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-land/)
#' @param kop_shp            .Rds object. Shapefile of the worlds koppen zones (e.g. https://www.climond.org/Koppen.aspx)
#' @param species_list       Character string - List of species analysed, used to cut the dataframe down
#' @param env_vars           Character string - List of environmental variables analysed
#' @param cell_size          Numeric. Value indicating the grid size in decimal degrees used for estimating Area of Occupancy (see ?AOO.computing)
#' @param save_run           Character string - run name to append to the data frame, useful for multiple runs.
#' @param save_data          Logical - do you want to save the data frame?
#' @param data_path          Character string - The file path used for saving the data frame
#' @export
calc_1km_niches = function(coord_df,
                           prj,
                           country_shp,
                           world_shp,
                           kop_shp,
                           #ibra_shp,
                           species_list,
                           env_vars,
                           cell_size,
                           save_run,
                           data_path,
                           save_data) {

  ## Create a spatial points object
  message('Estimating global niches for ', length(species_list), ' species across ',
          length(env_vars), ' climate variables')

  NICHE.1KM    <- coord_df %>% as.data.frame () %>%
    filter(coord_summary == TRUE)
  NICHE.1KM.84 <- SpatialPointsDataFrame(coords      = NICHE.1KM[c("lon", "lat")],
                                         data        = NICHE.1KM,
                                         proj4string = prj)

  ## Use a projected, rather than geographic, coordinate system
  ## Not sure why, but this is needed for the spatial overlay step
  AUS.WGS      = spTransform(country_shp,  prj)
  LAND.WGS     = spTransform(world_shp,    prj)
  KOP.WGS      = spTransform(kop_shp,      prj)

  ## Intersect the points with the Global koppen file
  message('Intersecting points with shapefiles for ', length(species_list), ' species')
  KOP.JOIN     = over(NICHE.1KM.84, KOP.WGS)

  ## Create global niche and Australian niche for website - So we need a subset for Australia
  ## The ,] acts just like a clip in a GIS
  NICHE.AUS <-  NICHE.1KM.84[AUS.WGS, ]

  ## Aggregate the number of Koppen zones (and IBRA regions) each species is found in
  COMBO.KOP <- NICHE.1KM.84 %>%
    cbind.data.frame(., KOP.JOIN)

  ## Aggregate the data
  KOP.AGG = tapply(COMBO.KOP$Koppen, COMBO.KOP$searchTaxon,
                   function(x) length(unique(x))) %>% ## group Koppen by species name
    as.data.frame()
  KOP.AGG =  setDT(KOP.AGG , keep.rownames = TRUE)[]
  names(KOP.AGG) = c("searchTaxon", "KOP_count")

  ## Run join between species records and spatial units :: SUA, POA and KOPPEN zones
  message('Joining occurence data to SUAs for ',
          length(species_list), ' species in the set ', "'", save_run, "'")

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
  GLOB.NICHE <- GLOB.NICHE %>% dplyr::select(-contains("."))
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
  AUS.NICHE <- GLOB.NICHE %>% dplyr::select(-contains("."))
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
  GLOB.NICHE <- list(GLOB.NICHE, GBIF.AOO, KOP.AGG) %>%
    reduce(left_join, by = "searchTaxon") %>%
    dplyr::select(searchTaxon, Aus_records, AOO, KOP_count, everything())

  if(save_data == TRUE) {

    ## save .rds file for the next session
    message('Writing 1km resolution niche and raster data for ',
            length(species_list), ' species in the set ', "'", save_run, "'")
    saveRDS(GLOB.NICHE, paste0(data_path, 'GLOBAL_NICHES_',  save_run, '.rds'))
    return(GLOB.NICHE)

  } else {
    message(' Return niches to the environment for ', length(species_list), ' species analysed')
    return(GLOB.NICHE)
  }

}





## Get a complete df ----

#' Remove species records with NA enviro (i.e. Raster) values - records outside the exetent of rasters
#'
#' This function downloads species occurrence files from GBIF (https://www.gbif.org/).
#' It assumes that the species list supplied is taxonomically correct.
#' It downloads the species without returning anything
#'
#' @param data            Data.frame of species records
#' @param desiredCols     Character string of columns to search for NA values EG c('rainfall', 'temp'), etc.
#' @return                A df with NA records removed
#' @export
completeFun <- function(data, desiredCols) {

  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])

}





## Estimate the niche from species records ----


#' This function takes a data frame of all species records,
#' and estimates the environmental niche for each species
#' It assumes that the input df is that prepared by the calc_1km_niches function
#' @param  DF                 Data.frame of all species records prepared by the calc_1km_niches function
#' @param  colname            Character string - the columns to estimate niches for E.G. 'rainfall', etc.
#' @return                    Data.frame of estimated environmental niches for each species
#' @export
niche_estimate = function (DF,
                           colname) {

  ## R doesn't seem to have a built-in mode function
  ## This doesn't really handel multiple modes...
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }

  ## Use ddply inside a function to create niche widths and medians for each species
  ## Also, need to figure out how to make the aggregating column generic (species, genus, etc.)
  ## Try to un-wire the grouping column using !!
  summary = ddply(DF,
                  .(searchTaxon),           ## currently grouping column only works hard-wired
                  .fun = function (xx, col) {

                    ## All the different columns
                    ## Calculate them all, and maybe supply a list of ones to keep
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
  ## Figure out how to make the aggregating column generic (species, genus, etc.)
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


#' This function takes a data frame of all species records,
#' and plots histograms and convex hulls for each species in global enviromental space
#' It assumes that the input df is that prepared by the check_spatial_outliers function
#' @param  coord_df           Data.frame of all species records prepared by the check_spatial_outliers function
#' @param  species_list       Character string - the species analysed
#' @param  range_path         Character string - file path for saving histograms and convex hulls
#' @export
plot_range_histograms = function(coord_df,
                                 species_list,
                                 range_path) {

  ## Subset the occurrence data to that used by maxent
  message('Plotting global environmental ranges for ', length(species_list), ' species')
  coord_df <- coord_df %>%
    as.data.frame()

  ## Plot histograms of temperature and rainfall
  ## spp = species_list[1]
  spp.plot = as.character(unique(coord_df$searchTaxon))
  for (spp in spp.plot) {

    ## Subset the spatial dataframe into records for each spp
    DF       <- coord_df[coord_df$searchTaxon %in% spp , ]
    DF.OCC   <- subset(coord_df, searchTaxon == spp & SOURCE != "INVENTORY")

    ## Create a new field RANGE, which is either URBAN or NATURAL
    coord_spp <- coord_df %>%
      filter(searchTaxon == spp) %>%
      mutate(RANGE = ifelse(SOURCE != "INVENTORY", "NATURAL", "URBAN"))

    ## Start PNG device
    message('Writing global convex hulls for ', spp)
    png(sprintf("%s%s_%s", range_path, spp, "1km_convex_hull.png"), 16, 10, units = 'in', res = 500)

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

    ## Check if file exists
    message('Writing global temp histograms for ', spp)
    png(sprintf("%s%s_%s", range_path, spp, "temp_niche_histograms_1km_records.png"),
        16, 10, units = 'in', res = 500)

    ## Use the 'SOURCE' column to create a histogram for each source.
    temp.hist <- ggplot(coord_spp, aes(x = Annual_mean_temp, group = RANGE, fill = RANGE)) +

      geom_density(aes(x = Annual_mean_temp, y = ..scaled.., fill = RANGE),
                   color = 'black', alpha = 1) +

      ## Change the colors
      scale_fill_manual(values = c('NATURAL' = 'coral2',
                                   'Mayor'   = 'gainsboro'), na.value = "grey") +

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

    ## Check if file exists
    png(sprintf("%s%s_%s", range_path, spp, "temp_niche_boxplots_1km_records.png"),
        10, 14, units = 'in', res = 500)

    ## Use the 'SOURCE' column to create a histogram for each source.
    temp.box = ggboxplot(coord_spp,
                         x    = 'RANGE',
                         y    = 'Annual_mean_temp',
                         fill = 'RANGE',
                         palette = c('coral2', 'gainsboro'), size = 0.6) +

      ## Use the classic theme
      theme_classic() +
      labs(y = 'Annual Mean Temp (°C)',
           x = '') +

      ## Add themes
      theme(axis.title.x     = element_text(colour = 'black', size = 25),
            axis.text.x      = element_blank(),

            axis.title.y     = element_text(colour = 'black', size = 35),
            axis.text.y      = element_text(size = 25),

            panel.border     = element_rect(colour = 'black', fill = NA, size = 1.2),
            plot.title       = element_text(size   = 25, face = 'bold'),
            legend.text      = element_text(size   = 20),
            legend.title     = element_blank(),
            legend.key.size  = unit(1.5, 'cm'))

    ## Print the plot and close the device
    print(temp.box + ggtitle(paste0("Worldclim temp niches for ", spp)))
    dev.off()


    ## Check if file exists
    message('Writing global rain histograms for ', spp)
    png(sprintf("%s%s_%s", range_path, spp, "rain_niche_histograms_1km_records.png"),
        16, 10, units = 'in', res = 500)

    ## Use the 'SOURCE' column to create a histogram for each source.
    rain.hist = ggplot(coord_spp, aes(x = Annual_precip, group = RANGE, fill = RANGE)) +

      # geom_histogram(position = "identity", alpha = 0.8, binwidth = 15,
      #                aes(y =..density..))  +
      geom_density(aes(x = Annual_precip, y = ..scaled.., fill = RANGE),
                   color = 'black', alpha = 0.8) +

      ## Change the colors
      scale_fill_manual(values = c('NATURAL' = 'coral2',
                                   'URBAN'   = 'gainsboro'), na.value = "grey") +

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

    ## Check if file exists
    message('Writing global rain boxplots for ', spp)
    png(sprintf("%s%s_%s", range_path, spp, "rain_niche_boxplots_1km_records.png"),
        10, 14, units = 'in', res = 500)

    ## Use the 'SOURCE' column to create a histogram for each source.
    rain.box = ggboxplot(coord_spp,
                         x    = 'RANGE',
                         y    = 'Annual_precip',
                         fill = 'RANGE',
                         palette = c('coral2', 'gainsboro'), size = 0.6) +

      ## Use the classic theme
      theme_classic() +
      labs(y = 'Annual Precipitation (mm)',
           x = '') +

      ## Add themes
      theme(axis.title.x     = element_text(colour = 'black', size = 25),
            axis.text.x      = element_blank(),

            axis.title.y     = element_text(colour = 'black', size = 35),
            axis.text.y      = element_text(size = 25),

            panel.border     = element_rect(colour = 'black', fill = NA, size = 1.2),
            plot.title       = element_text(size   = 25, face = 'bold'),
            legend.text      = element_text(size   = 20),
            legend.title     = element_blank(),
            legend.key.size  = unit(1.5, 'cm'))

    ## Print the plot and close the device
    print(rain.box + ggtitle(paste0("Worldclim rain niches for ", spp)))
    dev.off()

  }

}





## Create SDM table ----


#' This function takes a data frame of all species records,
#' And prepares a table in the 'species with data' (swd) format for modelling uses the Maxent algorithm.
#' It assumes that the input df is that returned by the coord_clean_records function
#' @param coord_df           Data.frame. DF of all species records returned by the coord_clean_records function
#' @param species_list       Character string - the species analysed
#' @param BG_points          Logical - Do we want to include a dataframe of background points? Otherwise, BG points taken from species not modelled
#' @param sdm_table_vars     Character string - The enviro vars to be included in species models
#' @param save_run           Character string - append a run name to the output (e.g. 'bat_species')
#' @param save_shp           Logical - Save a shapefile of the SDM table (T/F)?
#' @param read_background    Logical - Read in an additional dataframe of background points (T/F)?
#' @param save_data          Logical - do you want to save the data frame?
#' @param data_path          Character string - The file path used for saving the data frame
#' @return                   Data.frame of species records, with spatial outlier T/F flag for each record
#' @export
prepare_sdm_table = function(coord_df,
                             species_list,
                             BG_points,
                             sdm_table_vars,
                             save_run,
                             save_shp,
                             read_background,
                             save_data,
                             data_path) {

  ## Define Mollweide here. This is hard-wired, not user supplied
  sp_epsg54009 <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0"

  ## Just add clean_df to this step
  coord_df <- subset(coord_df, coord_summary == TRUE)
  message(round(nrow(coord_df)/nrow(coord_df)*100, 2), " % records retained")
  message(table(coord_df$coord_summary))

  ## RESTRICT DATA TABLE TO ONE RECORD PER 1KM GRID CELL
  ## Create a table with all the variables needed for SDM analysis
  message('Preparing SDM table for ', length(unique(coord_df$searchTaxon)),
          ' species in the set ', "'", save_run, "'",
          'using ', unique(coord_df$SOURCE), ' data')

  ## Select only the columns needed. This also needs to use the variable names
  coord_df         <- coord_df[coord_df$searchTaxon %in% species_list, ]
  length(unique(coord_df$searchTaxon))

  COMBO.RASTER.ALL  <- coord_df %>%
    dplyr::select(one_of(sdm_table_vars))


  ## Create a spatial points object, and change to a projected system to calculate distance more accurately
  ## This is the mollweide projection used for the SDMs
  coordinates(COMBO.RASTER.ALL)    <- ~lon+lat
  proj4string(COMBO.RASTER.ALL)    <- '+init=epsg:4326'
  COMBO.RASTER.ALL                 <- spTransform(COMBO.RASTER.ALL, CRS(sp_epsg54009))

  ## Don't filter the data again to be 1 record per 1km, that has already happened
  SDM.DATA.ALL <- COMBO.RASTER.ALL


  ## TRY CLEANING FILTERED DATA FOR SPATIAL OUTLIERS
  ## The cc_outl function has been tweaked and sped up.

  ## Create a unique identifier for spatial cleaning.
  ## This is used for automated cleaing of the records, and also saving shapefiles
  ## But this will not be run for all species linearly.
  ## So, it probably needs to be a combination of species and number
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
    spTransform(., CRS("+init=epsg:4326")) %>%
    as.data.frame() %>%
    dplyr::select(searchTaxon, lon, lat, SPOUT.OBS, SOURCE) %>%
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
    dplyr::select(species) %>%
    table() %>%
    as.data.frame()
  COMBO.LUT <- setDT(COMBO.LUT, keep.rownames = FALSE)[]
  names(COMBO.LUT) = c("species", "FREQUENCY")
  COMBO.LUT = COMBO.LUT[with(COMBO.LUT, rev(order(FREQUENCY))), ]

  ## Watch out here - this sorting could cause problems for the order of
  ## the data frame once it's stitched back together
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
                         method  = "quantile",
                         mltpl   = 10,
                         value   = "flagged",
                         verbose = TRUE)

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
  message('Is the order or records identical before joining?',
          identical(nrow(SDM.COORDS), nrow(SPAT.OUT)))
  message('Is the order or records identical after joining?',
          identical(SDM.DATA.ALL$searchTaxon, SPAT.OUT$searchTaxon))
  length(unique(SPAT.OUT$searchTaxon))

  ## This explicit join is required. Check the species have been analysed in exactly the same order
  SPAT.FLAG <- join(as.data.frame(SDM.DATA.ALL), SPAT.OUT,
                    by = c("SPOUT.OBS", "searchTaxon") ,
                    type = "left", match = "first")
  message('Is the order or records identical after joining?',
          identical(SDM.DATA.ALL$searchTaxon, SPAT.FLAG$searchTaxon))

  ## Check the join is working
  message('Checking spatial flags for ', length(unique(SPAT.FLAG$searchTaxon)),
          ' species in the set ', "'", save_run, "'")
  message(table(SPAT.FLAG$SPAT_OUT, exclude = NULL))

  ## Just get the records that were not spatial outliers.
  SDM.SPAT.ALL <- subset(SPAT.FLAG, SPAT_OUT == TRUE)
  unique(SDM.SPAT.ALL$SPAT_OUT)
  unique(SDM.SPAT.ALL$SOURCE)
  length(unique(SDM.SPAT.ALL$searchTaxon))

  ## What percentage of records are retained?
  message(round(nrow(SDM.SPAT.ALL)/nrow(SPAT.FLAG)*100, 2),
          " % records retained after spatial outlier detection")

  ## Convert back to format for SDMs :: use Mollweide projection
  SDM.SPAT.ALL = SpatialPointsDataFrame(coords      = SDM.SPAT.ALL[c("lon", "lat")],
                                        data        = SDM.SPAT.ALL,
                                        proj4string = CRS(sp_epsg54009))
  projection(SDM.SPAT.ALL)
  message(length(unique(SDM.SPAT.ALL$searchTaxon)),
          ' species processed through from download to SDM table')

  ## CREATE SHAPEFILES TO CHECK OUTLIERS ARCMAP
  ## Rename the fields so that ArcMap can handle them
  SPAT.OUT.CHECK = SPAT.FLAG %>%
    dplyr::select(SPOUT.OBS, searchTaxon, lat, lon, SOURCE, SPAT_OUT) %>%
    dplyr::rename(TAXON     = searchTaxon,
                  LAT       = lat,
                  LON       = lon)
  names(SPAT.OUT.CHECK)

  ## Then create a SPDF
  SPAT.OUT.SPDF    = SpatialPointsDataFrame(coords      = SPAT.OUT.CHECK[c("LON", "LAT")],
                                            data        = SPAT.OUT.CHECK,
                                            proj4string = CRS("+init=epsg:4326"))

  ## Write the shapefile out
  if(save_shp == TRUE) {

    ## save .shp for future refrence
    writeOGR(obj    = SPAT.OUT.SPDF,
             dsn    = "./data/ANALYSIS/CLEAN_GBIF",
             layer  = paste0('SPAT_OUT_CHECK_', save_run),
             driver = "ESRI Shapefile", overwrite_layer = TRUE)

  } else {
    message(' skip file saving, not many species analysed')   ##
  }

  ## CREATE BACKGROUND POINTS AND VARIBALE NAMES
  ## Use one data frame for all species analysis,
  ## to save mucking around with background points
  ## Here we are using a dataframe of mammals, reptiles and birds.
  if(read_background == TRUE) {

    Message('Read in background data for taxa analaysed')
    background.points = readRDS(paste0(data_path, BG_points)) %>%
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
  sdm_cols <- names(dplyr::select(SDM.SPAT.OCC.BG@data, searchTaxon, lon, lat, SOURCE, everything()))
  SDM.SPAT.OCC.BG <- SDM.SPAT.OCC.BG[,!(names(SDM.SPAT.OCC.BG) %in% drops)]

  ## save data
  if(save_data == TRUE) {

    ## Save .rds file of the occurrence and BG points for the next session
    saveRDS(SDM.SPAT.OCC.BG, paste0(data_path, 'SDM_SPAT_OCC_BG_',  save_run, '.rds'))

  } else {
    message('Return the occurrence + Background data to the global environment')   ##
    return(SDM.SPAT.OCC.BG)
  }

  ## get rid of some memory
  gc()

}





#########################################################################################################################
####################################################  TBC  ##############################################################
#########################################################################################################################
