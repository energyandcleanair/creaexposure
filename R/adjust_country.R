adjust_country <- function(res, suffix, iso3, selected_regions, year){
  bb <- data.gadm0() %>% subset(GID_0==iso3) %>% sf::st_as_sf() %>% sf::st_bbox()
  pop <- data.pop(res=res, bb=bb)
  adjust_global(pop, res, suffix=suffix, selected_regions=selected_regions, year=year)
}
