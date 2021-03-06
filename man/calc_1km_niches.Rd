% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SDM_GEN_PROCESSOR_FUNCTIONS.R
\name{calc_1km_niches}
\alias{calc_1km_niches}
\title{This function takes a data frame of all species records,
estimates geographic and environmental ranges for each, and creates a table of each.
It uses the AOO.computing function in the ConR package https://cran.r-project.org/web/packages/ConR/index.html
It assumes that the input df is that returned by the check_spatial_outliers function}
\usage{
calc_1km_niches(
  coord_df,
  prj,
  country_shp,
  world_shp,
  kop_shp,
  species_list,
  env_vars,
  cell_size,
  save_run,
  data_path,
  save_data
)
}
\arguments{
\item{coord_df}{Data.frame. DF of all species records returned by the coord_clean_records function}

\item{world_shp}{.Rds object. Shapefile of the worlds land (e.g. https://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-land/)}

\item{kop_shp}{.Rds object. Shapefile of the worlds koppen zones (e.g. https://www.climond.org/Koppen.aspx)}

\item{species_list}{Character string - List of species analysed, used to cut the dataframe down}

\item{env_vars}{Character string - List of environmental variables analysed}

\item{cell_size}{Numeric. Value indicating the grid size in decimal degrees used for estimating Area of Occupancy (see ?AOO.computing)}

\item{save_run}{Character string - run name to append to the data frame, useful for multiple runs.}

\item{data_path}{Character string - The file path used for saving the data frame}

\item{save_data}{Logical - do you want to save the data frame?}

\item{aus_df}{Data.frame of Urban records (only used if you have urban data, e.g. I-naturalist)}
}
\description{
This function takes a data frame of all species records,
estimates geographic and environmental ranges for each, and creates a table of each.
It uses the AOO.computing function in the ConR package https://cran.r-project.org/web/packages/ConR/index.html
It assumes that the input df is that returned by the check_spatial_outliers function
}
