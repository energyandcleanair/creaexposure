adjust_country <- function(res, suffix, iso3, selected_regions, year, poll=c('pm25', 'no2')){
  init()
  bb <- data.gadm0() %>% subset(GID_0==iso3) %>% sf::st_as_sf() %>% sf::st_bbox()
  pop <- data.pop(res=res, bb=bb)
  adjust_global(pop=pop,
                res=res,
                suffix=suffix,
                selected_regions=selected_regions,
                year=year,
                poll=poll)
}
