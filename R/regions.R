get_regions <- function(selected_regions=NULL){
    regions <- list(
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
      "AU" = "AU",
      "BD" = "BD"
    )

    if (!is.null(selected_regions)) {
      regions <- regions[selected_regions]
    }

    return(regions)
  }


add_region_to_obs <- function(obs, regions){
  obs$region <- NA
  for (r in names(regions)) {
    obs[obs$country %in% regions[[r]], "region"] <- r
  }
  return(obs)
}

get_bb <- function(selected_regions = NULL){
  if(is.null(selected_regions)) return(NULL)

  iso2s <- unlist(get_regions(selected_regions))
  iso3s <- countrycode::countrycode(iso2s, "iso2c", "iso3c")

  bb <- data.gadm0() %>% subset(GID_0 %in% iso3s) %>% sf::st_as_sf() %>% sf::st_bbox()
  return(bb)
}
