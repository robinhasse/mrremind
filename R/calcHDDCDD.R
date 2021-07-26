#' @title calc HDD and CDD
#' @description heating and cooling degree days, driver for space heating and 
#'   cooling demand in buildings
#'
#' @return magpie object of heating and cooling degree days
#' @author Robin Krekeler
#' @examples
#' 
#' \dontrun{ 
#' calcHDDCDD()
#' }
#' 
#' @importFrom madrat getConfig
#' @importFrom raster brick names cellStats subset stackApply getValues ncells 
#'   xyFromCell res aggregate
#' @importFrom ncdf4 nc_open 
#' @importFrom tidyr %>% mutate
#' @importFrom rlang .data


calcHDDCDD <- function() {

  # Functions -------------------------------------------------------------

  
  # this function has to be integrated into mrcommons::readISIMIP, once 
  # madrat::readSource allows to pass raster objects
  .readISIMIP <- function(subtype) {
    if (grepl("tas", subtype)) {
      cwd <- getwd()
      setwd(file.path(getConfig("sourcefolder"), "ISIMIP", subtype))
      files <- Sys.glob("*.nc")
      r <- raster::brick(files)
      setwd(cwd)
      return(r)
    }
    if (grepl("population", subtype)) {
      cwd <- getwd()
      setwd(file.path(getConfig("sourcefolder"), "ISIMIP", subtype))
      files <- Sys.glob("*.nc4")
      r <- raster::brick(files)
      
      # rename years
      years <- strsplit(tail(strsplit(subtype, "_")[[1]], 1), "-")[[1]]
      names(r) <- paste0("y", years[1]:years[2])
      
      # filter relevant years
      r <- raster::dropLayer(
        r, as.numeric(substr(names(r), 2, 5)) < firstHistYear)
      
      # aggregate to common resolution of 0.5 deg
      if (any(raster::res(r) != 0.5)) {
        r <- raster::aggregate(r, fun = "sum",
                               fact = round(0.5 / raster::res(r), 3))
        raster::res(r)    <- 0.5
        raster::extent(r) <- round(raster::extent(r))
      }
      
      setwd(cwd)
      return(r)
    }
    if (grepl("CountryMask", subtype)) {
      cwd <- getwd()
      setwd(file.path(getConfig("sourcefolder"), "ISIMIP", subtype))
      files <- Sys.glob("*.nc")
      vars <- names(ncdf4::nc_open(files)[["var"]])
      countries <- list()
      for (var in vars) {
        countries[[var]] <- raster::brick(files, varname = var)
      }
      r <- raster::brick(countries)
      names(r) <- gsub("m_", "", vars)
      setwd(cwd)
      return(r)
    }
  }

  
  # calculate HDD/CDD in each cell from daily temperature
  calcCellHDDCDD <- function(.temp, .typeDD, .tlim) {
    # extract years
    years <- names(.temp) %>%
      substr(2, 5) %>%
      as.numeric()
    
    # calculate daily positive temperature deviation [K]
    if (.typeDD == "HDD") {
      hddcdd <- .tlim - (.temp - 273.15)
    } else if (.typeDD == "CDD") {
      hddcdd <- (.temp - 273.15) - .tlim 
    } else {
      stop(paste(".typeDD has to be either 'HDD' or 'CDD', not", .typeDD))
    }
    hddcdd[hddcdd < 0] <- 0
    
    # aggregate to yearly HDD/CDD [K.d/a]
    hddcdd <- raster::stackApply(hddcdd, years, fun = sum)
    names(hddcdd) <- gsub("index_", "y", names(hddcdd))
    
    return(hddcdd)
  }
  
  
  # aggregate cellular HDD/CDD to country-wide average (population-weighted)
  aggCells <- function(r, weight, mask) {
    years_r <- names(r)
    years_w <- names(weight)
    
    # loop: years in raster file r
    hddcdd_agg <- do.call(
      "rbind",
      lapply(
        years_r,
        function(y) {
          # filter relevant weights for year
          if (y <= lastHistYear) {
            weight_y <- weight$historical
          } else {
            weight_y <- weight[[-"historical"]]
          }
          
          # loop: SSP scenarios
          do.call(
            "rbind",
            lapply(
              names(weight_y),
              function(ssp) {
                tmp <- raster::subset(r, y) * raster::subset(weight_y[ssp], y) * mask
                tmp_tot <- raster::subset(weight_y[ssp], y) * mask
                tmp <- raster::cellStats(tmp, "sum") / 
                  raster::cellStats(tmp_tot, "sum")
                tmp <- data.frame("Region" = names(mask),
                                  "Period" = y,
                                  "Value"  = round(tmp, 1),
                                  "SSP"    = ssp)
                rownames(tmp) <- c()
                return(tmp)
              }
            )
          )
        }
      )
    )
  }
  

  # calculate all desired output from given tas file
  calcStackHDDCDD <- function(file, tlim, countries, pop) {
    # read cellular temperature
    temp <- .readISIMIP(file)
    years_t <- unique(substr(names(temp), 2, 5))
    
    # loop: typeDD
    hddcdd <- do.call(
      "rbind",
      lapply(
        c("HDD", "CDD"),
        function(typeDD) {
          
          # loop: threshold temperatures
          do.call(
            "rbind",
            lapply(
              tlim[typeDD],
              function(t) {
                hddcdd_agg <- calcCellHDDCDD(temp, typeDD, t) %>%
                  aggCells(pop, countries) %>%
                  mutate("Variable"    = typeDD,
                         "Temperature" = t)
              }
            )
          )
        }
      )
    )
  }
  
    
  plotRaster <- function(r) {
    library("ggplot2")
    r_map <- stack(as.data.frame(raster::getValues(r)))
    coords <- raster::xyFromCell(r, seq_len(raster::ncell(r)))
    names(r_map) <- c('value', 'variable')
    r_map <- cbind(coords, r_map)
    ggplot(r_map) + 
      geom_tile(aes(x, y, fill = value)) +
      facet_wrap(~ variable) +
      scale_fill_gradientn(colours = rev(terrain.colors(225))) +
      coord_equal()
  }

  

  # Setup -----------------------------------------------------------------


  # threshold temperature for heating and cooling [deg C]
  tlim <- list("HDD" = seq(17, 25), "CDD" = seq(17, 25))
  
  # historical years
  firstHistYear <- 1950
  lastHistYear  <- 2009
  
  # list of files that are processed
  files <- toolGetMapping("HDDCDDfiles.csv", type = "sectoral")

  

  # Calculation -----------------------------------------------------------
  
  
  # cells -> country  
  countries <- .readISIMIP("CountryMask")
  
  # population data as weights
  pop <- files %>%
    filter(.data$variable == "pop") %>%
    select("ssp", "file") %>%
    mutate("file" = gsub(".nc4", "", .data$file))
  pop <- split(pop$file, pop$ssp) %>%
    sapply(".readISIMIP", USE.NAMES = TRUE)
  
  # extend historical population (SSP scenarios are identical until 2009)
  pop$historical <- raster::stack(
    pop$historical,
    raster::subset(
      pop[["1"]],
      which(as.numeric(gsub("y", "", names(pop[["1"]]))) <= lastHistYear)
    )
  )
  
  # loop: GCM results for ambient temperature (RCP scenarios)
  hddcdd <- do.call(
    "rbind",
    apply(
      files[files$variable == "tas", c("file", "rcp")], 1,
      function(f, rcp) {
        hddcdd_cell <- calcStackHDDCDD(f, tlim, countries, pop) %>%
          mutate("RCP" = rcp)

      }
    )
  )

}
