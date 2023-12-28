get_region <- function(){
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
      "AU" = "AU"
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

