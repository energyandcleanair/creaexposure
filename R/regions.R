get_all_regions <- function(){
  list(
      "NA" = c("US", "CA", "MX"),
      "IN" = "IN",
      "CN" = "CN",
      "EU" = c(
        "AD", "AT", "BA", "BE", "BG", "CH", "CY", "CZ",
        "DE", "DK", "EE", "ES", "FI", "FR", "GB", "GI", "GR",
        "HR", "HU", "IE", "IS", "IT", "LT", "LU", "LV", "ME",
        "MK", "MT", "NL", "NO", "PL", "PT", "RO", "RS", "SE", "SI", "SK", "XK"
      ),
      "TR" = "TR",
      "SEA" = "TH", # No other countries with enough measurements
      "JP" = "JP",
      "ZA" = "ZA",
      "PH" = "PH",
      "TW" = "TW",
      "CL" = "CL",
      "AU" = "AU"
    )
}



get_region_iso2s <- function(region){

  # If a list, values should be iso2s
  if(is.list(region)) return(unname(unlist(region)))

  # If a string, look into the regions list
  if(is.character(region)){
    definitions <- get_all_regions()
    if(region %in% names(definitions)){
      return(unlist(definitions[[region]]))
    } else {
      message("Region not found. Assuming this is an iso2")
      return(region)
    }
  }
  stop("region should be a character or a list of characters")
}

get_region_names <- function(regions){
  if(is.null(regions)) return(NULL)
  if(is.character(regions)) return(regions)
  if(is.list(regions)) return(names(regions))
  stop("regions should be a character or a list of characters")
}

get_bbox <- function(selected_regions = NULL){

  if(is.null(selected_regions)) return(NULL)

  iso2s <- get_region_iso2s(selected_regions)
  iso3s <- countrycode::countrycode(iso2s, "iso2c", "iso3c")

  bbox <- data.gadm0() %>% subset(GID_0 %in% iso3s) %>% sf::st_as_sf() %>% sf::st_bbox()
  return(bbox)
}


add_region_to_obs <- function(obs, regions = get_all_regions()){
  obs$region <- NA
  for (r in names(regions)) {
    iso2s <- get_region_iso2s(regions[r])
    obs[obs$country %in% regions[[r]], "region"] <- r
  }
  return(obs)
}
