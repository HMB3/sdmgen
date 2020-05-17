#########################################################################################################################
###################################  FUNCTIONS FOR RUNNING SDM ANALYSIS ---- ############################################
#########################################################################################################################


## Below are the functions used to run the SDM models


## Run the SDM analysis ----


#' This function takes a data frame of all species records,
#' and runs a specialised maxent analysis for each species.
#' It uses the rmaxent package https://github.com/johnbaums/rmaxent
#' It assumes that the input df is that returned by the prepare_sdm_table function
#' @param species_list       Character string - the species to run maxent models for
#' @param sdm_df             SpatialPointsDataFrame. Spdf of all species records returned by the 'prepare_sdm_table' function
#' @param sdm_predictors     Character string - Vector of enviro conditions that you want to include
#' @param maxent_dir         Character string - The file path used for saving the maxent output
#' @param bs_dir             Character string - The file path used for saving the backwards selection maxent output
#' @param backwards_sel      Logical - Run backwards selection using the maxent models (T/F)?
#' @param template_raster    RasterLayer - Empty raster with analysis extent (global), res (1km) and projection (Mollweide, EPSG 54009)
#' @param cor_thr            Numeric - The max allowable pairwise correlation between predictor variables
#' @param pct_thr            Numeric - The min allowable percent variable contribution
#' @param k_thr              Numeric - The min number of variables to be kept in the model
#' @param min_n              Numeric - The min number of records for running maxent
#' @param max_bg_size        Numeric - The max number of background points to keep
#' @param background_buffer_width Numeric - The max distance (km) from occ points that BG points should be selected?
#' @param shapefiles         Logical - Save shapefiles of the occ and bg data (T/F)?
#' @param features           Character string - Which features should be used? (e.g. linear, product, quadratic 'lpq')
#' @param replicates         Numeric - The number of replicates to use
#' @param responsecurves     Logical - Save response curves of the maxent models (T/F)?
#' @param country_shp             .Rds object - SpatialPolygonsDataFrame of Australia for mapping maxent points
#' @param Koppen_raster      RasterLayer of global koppen zones, in Mollweide54009 projection
#' @param Koppen_zones       Dataframe of global koppen zones, with columns : GRIDCODE, Koppen
#' @export
run_sdm_analysis = function(species_list,
                            sdm_df,
                            sdm_predictors,
                            maxent_dir,
                            bs_dir,
                            backwards_sel,
                            template_raster,
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
                            Koppen_raster,
                            Koppen_zones,
                            country_shp ) {

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
        fit_maxent_targ_bg_back_sel(occ                     = occurrence,
                                    bg                      = background,
                                    sdm_predictors          = sdm_predictors,
                                    name                    = species,
                                    outdir                  = maxent_dir,
                                    bsdir                   = bs_dir,
                                    backwards_sel           = "TRUE",
                                    cor_thr                 = cor_thr,
                                    pct_thr                 = pct_thr,
                                    k_thr                   = k_thr,

                                    template_raster         = template_raster,
                                    min_n                   = min_n,
                                    max_bg_size             = max_bg_size,
                                    Koppen_raster           = Koppen_raster,
                                    background_buffer_width = background_buffer_width,
                                    shapefiles              = shapefiles,
                                    features                = features,
                                    replicates              = replicates,
                                    responsecurves          = responsecurves,
                                    #shp_path                = shp_path,
                                    country_shp             = country_shp ),


        ## Save error message
        error = function(cond) {

          ## How to put the message into the file?
          file.create(file.path(dir_name, "sdm_failed.txt"))
          message(species, ' failed')
          cat(cond$message, file = file.path(dir_name, "sdm_failed.txt"))
          warning(species, ': ', cond$message)

        })

    } else {
      message(species, ' skipped - no data.') ## This condition ignores species which have no data...
      file.create(file.path(dir_name, "completed.txt"))
    }

    ## now add a file to the dir to denote that it has completed
    file.create(file.path(dir_name, "completed.txt"))

  })

}





## Run maxent with backwards selection ----


#' This function takes a data frame of all species records,
#' and runs a specialised maxent analysis for each species.
#' It uses the rmaxent package https://github.com/johnbaums/rmaxent
#' It assumes that the input df is that returned by the prepare_sdm_table function
#' @param occ                SpatialPointsDataFrame - Spdf of all species records returned by the 'prepare_sdm_table' function
#' @param bg                 SpatialPointsDataFrame - Spdf of of candidate background points
#' @param sdm_predictors     Character string - Vector of enviro conditions that you want to include
#' @param outdir             Character string - The file path used for saving the maxent output
#' @param bsdir              Character string - The file path used for saving the backwards selection maxent output
#' @param backwards_sel      Logical - Run backwards selection using the maxent models (T/F)?
#' @param template_raster    RasterLayer -  Empty raster with analysis extent (global), res (1km) and projection (Mollweide, EPSG 54009)
#' @param cor_thr            Numeric - The max allowable pairwise correlation between predictor variables
#' @param pct_thr            Numeric - The min allowable percent variable contribution
#' @param k_thr              Numeric - The min number of variables to be kept in the model
#' @param min_n              Numeric - The min number of records for running maxent
#' @param max_bg_size        Numeric - The max number of background points to keep
#' @param background_buffer_width Numeric - The max distance (km) from occ points that BG points should be selected?
#' @param shapefiles         Logical - Save shapefiles of the occ and bg data (T/F)?
#' @param features           Character string - Which features should be used? (e.g. linear, product, quadratic 'lpq')
#' @param replicates         Numeric - The number of replicates to use
#' @param responsecurves     Logical - Save response curves of the maxent models (T/F)?
#' @param country_shp             .Rds object - SpatialPolygonsDataFrame of Australia for mapping maxent points
#' @param rep_args             RasterLayer of global koppen zones, in Mollweide54009 projection
#' @param full_args          Dataframe of global koppen zones, with columns : GRIDCODE, Koppen
#' @export


#' @export
fit_maxent_targ_bg_back_sel <- function(occ,
                                        bg,
                                        sdm_predictors,
                                        name,
                                        outdir,
                                        bsdir,
                                        cor_thr,
                                        pct_thr,
                                        k_thr,
                                        backwards_sel,
                                        template_raster,
                                        min_n,
                                        max_bg_size,
                                        background_buffer_width,
                                        Koppen_raster,
                                        shapefiles,
                                        features,
                                        replicates,
                                        responsecurves,
                                        country_shp ) {


  ## First, stop if the outdir file exists,
  if(!file.exists(outdir)) stop('outdir does not exist :(', call. = FALSE)
  outdir_sp <- file.path(outdir, gsub(' ', '_', name))
  bsdir_sp  <- file.path(bsdir,  gsub(' ', '_', name))

  if(!missing('Koppen_raster')) {
    if(!is(Koppen_raster, 'RasterLayer'))
      stop('Koppen must be a RasterLayer, and should be in the same coordinate system as template_raster')
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

  ## Get unique cell numbers for species occurrences
  cells <- cellFromXY(template_raster, occ)

  ## Clean out duplicate cells and NAs (including points outside extent of predictor data)
  ## Note this will get rid of a lot of duplicate records not filtered out by GBIF columns, etc.
  not_dupes <- which(!duplicated(cells) & !is.na(cells))
  occ       <- occ[not_dupes, ]
  cells     <- cells[not_dupes]
  message(nrow(occ), ' occurrence records (unique cells).')

  ## Skip species that have less than a minimum number of records: eg 20 species
  if(nrow(occ) < min_n) {

    print (paste ('Fewer occurrence records than the number of cross-validation ',
                  'replicates for species ', name,
                  ' Model not fit for this species'))

  } else {

    ## Subset the background records to the 200km buffered polygon
    message(name, ' creating background cells')
    system.time(o <- over(bg, buffer))
    bg <- bg[which(!is.na(o)), ]
    bg_cells <- cellFromXY(template_raster, bg)

    ## Clean out duplicates and NAs (including points outside extent of predictor data)
    bg_not_dupes <- which(!duplicated(bg_cells) & !is.na(bg_cells))
    bg           <- bg[bg_not_dupes, ]
    bg_cells     <- bg_cells[bg_not_dupes]

    ## Find which of these cells fall within the Koppen-Geiger zones that the species occupies
    ## Crop the Kopppen raster to the extent of the occurrences, and snap it
    message(name, ' intersecting background cells with Koppen zones')
    Koppen_crop <- crop(Koppen_raster, occ, snap = 'out')

    ## Only extract and match those cells that overlap between the ::
    ## 1). cropped koppen zone,
    ## 2). occurrences and
    ## 3). background points
    message(xres(template_raster), ' metre cell size for template raster')
    message(xres(Koppen_raster), ' metre cell size for Koppen raster')
    zones               <- raster::extract(Koppen_crop, occ)
    cells_in_zones_crop <- Which(Koppen_crop %in% zones, cells = TRUE)
    cells_in_zones      <- cellFromXY(Koppen_raster, xyFromCell(Koppen_crop, cells_in_zones_crop))
    bg_cells            <- intersect(bg_cells, cells_in_zones)  ## this is 0 for 5km
    i                   <- cellFromXY(template_raster, bg)
    bg                  <- bg[which(i %in% bg_cells), ]

    ## For some species, we have the problem that the proportion of ALA/INV data is
    ## very different in the occurrence vs the bg records.
    ## This should be caused by the 200km / koppen restriction, etc.

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

    ## Now save an image of the background points
    ## This is useful to quality control the models
    save_name = gsub(' ', '_', name)

    aus.mol = country_shp  %>%
      spTransform(projection(buffer))

    aus.kop = crop(Koppen_crop, aus.mol)
    occ.mol <- occ %>%
      spTransform(projection(buffer))

    ## Print the koppen zones, occurrences and points to screen
    # plot(Koppen_crop, legend = FALSE,
    #      main = paste0('Occurence SDM records for ', name))
    #
    # plot(aus.mol, add = TRUE)
    # plot(buffer,  add = TRUE, col = "red")
    # plot(occ.mol, add = TRUE, col = "blue")

    ## Then save the occurrence points
    png(sprintf('%s/%s/%s_%s.png', outdir, save_name, save_name, "buffer_occ"),
        16, 10, units = 'in', res = 300)

    plot(Koppen_crop, legend = FALSE,
         main = paste0('Occurence SDM records for ', name))

    plot(aus.mol, add = TRUE)
    plot(buffer,  add = TRUE, col = "red")
    plot(occ.mol, add = TRUE, col = "blue")

    dev.off()

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

    ## SWD = species with data. Now sample the environmental
    ## variables used in the model at all the occ and bg points
    swd_occ <- occ[, sdm_predictors]
    saveRDS(swd_occ, file.path(outdir_sp, paste0(save_name,'_occ_swd.rds')))

    swd_bg <- bg.samp[, sdm_predictors]
    saveRDS(swd_bg, file.path(outdir_sp, paste0(save_name, '_bg_swd.rds')))

    ## Save the SWD tables as shapefiles
    if(shapefiles) {

      writeOGR(swd_occ, outdir_sp,  paste0(save_name, '_occ_swd'), 'ESRI Shapefile', overwrite_layer = TRUE)
      writeOGR(swd_bg,  outdir_sp,  paste0(save_name, '_bg_swd'),  'ESRI Shapefile', overwrite_layer = TRUE)

    }

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

      ## Only use rep_args & full_args if not using replicates
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





## Simplify the maxent models ----


#' This is a local version of rmaxent::simplify (https://github.com/johnbaums/rmaxent)
#' Given a candidate set of predictor variables, this function identifies
#' a subset that meets specified multicollinearity criteria. Subsequently,
#' backward stepwise variable selection is used to iteratively drop the variable
#' that contributes least to the model, until the contribution of each variable
#' meets a specified minimum, or until a predetermined minimum number of predictors remains.
#' It assumes that the input df is that returned by the fit_maxent_targ_bg_back_sel function
#' @param occ                SpatialPointsDataFrame - Spdf of all species records returned by the 'prepare_sdm_table' function
#' @param bg                 SpatialPointsDataFrame - Spdf of of candidate background points
#' @param path               Character string - Vector of enviro conditions that you want to include
#' @param species_column     Character string - Vector of enviro conditions that you want to include
#' @param cor_thr            Numeric - The max allowable pairwise correlation between predictor variables
#' @param pct_thr            Numeric - The min allowable percent variable contribution
#' @param k_thr              Numeric - The min number of variables to be kept in the model
#' @param features           Character string - Which features should be used? (e.g. linear, product, quadratic 'lpq')
#' @param replicates         Numeric - The number of replicates to use
#' @param type               The variable contribution metric to use when dropping variables
#' @param logistic_format    Logical value indicating whether maxentResults.csv should report logistic value thresholds
#' @param responsecurves     Logical - Save response curves of the maxent models (T/F)?
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
    stop("features must be a vector of one or more of ',\n  'l', 'p', 'q', 'h', and 't'.")
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
                                                      args, value = TRUE,
                                                      invert = TRUE), path = d)
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
                                                    args, value = TRUE,
                                                    invert = TRUE), path = d)
    }
    if (isTRUE(save)) {
      saveRDS(m, file.path(d, "maxent_fitted.rds"))
    }
    return(m)
  }
  lapply(names(occ_by_species), f)
}





## Compile the SDM data ----


#' This function extracts the SDM results from the folders.
#' It assumes that the input folders are those returned by the 'fit_maxent_targ_bg_back_sel' function
#' @param species_list      Character string - the species to run maxent models for
#' @param results_dir       Character string - The file path used for saving the maxent output
#' @param save_data         Logical or character - do you want to save the data frame?
#' @param data_path         Character string - The file path used for saving the data frame
#' @param save_run          Character string - run name to append to the data frame, useful for multiple runs.
#' @return                  Data.frame of maxent results
#' @export
compile_sdm_results = function(species_list,
                               results_dir,
                               save_data,
                               data_path,
                               save_run) {

  ## The code that adds niche info is now in './R/9_COLLATE_MAXENT_RESULTS.R'
  message('Creating summary stats for ', length(species_list),
          ' species in the set ', "'", save_run, "'")

  ## First, make a list of all the species with models, then restrict them
  ## to just the models on the species_list list
  map_spp_list  = gsub(" ", "_", species_list)
  map_spp_patt  = paste0(map_spp_list, collapse = "|")
  message ("map_spp_list head:")
  message (paste (head(map_spp_list), collapse=","))

  ## Now stop R from creating listing all the maxent files that have completed - this takes a long time
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

  ## Now create a list of the '10th percentile training presence Logistic threshold'.
  ## This is used in step 8 to threshold
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
  if(save_data == TRUE) {

    ## If saving, save
    saveRDS(MAXENT.RESULTS, paste0(data_path, 'MAXENT_RESULTS_', save_run, '.rds'))

  } else {
    ## Or return to the global environment
    message(' skip file saving, not many species analysed')
    return(MAXENT.RESULTS)
  }

}





#########################################################################################################################
####################################################  TBC  ##############################################################
#########################################################################################################################
