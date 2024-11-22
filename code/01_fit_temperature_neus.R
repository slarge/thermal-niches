# remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE)remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE)

# install.packages("Matrix")
# install.packages("TMB")
# install.packages("sdmTMB")

library(sdmTMB)
library(dplyr)
library(tidyr)
library(lubridate)
source("code/utilities.R")

## Switch off spherical geometry
sf::sf_use_s2(FALSE)

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
## Load and aggregate NEFSC BTS to functional groups ##
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

## prepare bottom depth
xlims <- c(-77, -65)
ylims <- c(35, 45)
res <- 1
bath_filename <- sprintf("marmap_coord_%s;%s;%s;%s_res_%s.csv",
                         xlims[1], ylims[1], xlims[2], ylims[2], res)

if(!bath_filename %in% list.files("grid_data/")){
  # 200 m isobath layer
  nesbath <- marmap::getNOAA.bathy(lon1 = xlims[1], lon2 = xlims[2],
                                   lat1 = ylims[1], lat2 = ylims[2],
                                   resolution = res,
                                   keep = TRUE) %>%
    marmap::as.raster()
  
  file.copy(bath_filename, "grid_data/")
  file.remove(bath_filename)
} else {
  nesbath <- marmap::read.bathy(sprintf("grid_data/%s", bath_filename), header = T) %>%
    marmap::as.raster()
}
## OK for now, but should figure out how to pull directly from Sean's survdat package
survdat <- readRDS(url("https://github.com/NOAA-EDAB/ecodata/raw/master/data-raw/survdat.rds"))$survdat

## Load EPU shapefiles as a spatialpolygonsdataframe
wg_crs <- sf::st_crs(ecodata::epu_sf)
crs <- sf::st_crs("+proj=utm +zone=19 +datum=WGS84 +units=km")

epu <- ecodata::epu_sf %>% 
  sf::st_transform(crs = wg_crs) %>% 
  select(EPU, geometry)

## Post stratify data according to EPUs
nefsc_epu <- survdat %>%
  mutate(depth = raster::extract(nesbath, y = cbind(.$LON, .$LAT)) * -1,
         depth = ifelse(is.na(DEPTH),
                               depth,
                               DEPTH),
         yday = yday(as.Date(EST_TOWDATE)),
         date = as.Date(EST_TOWDATE),
         logdepth = log(depth),
         cpue_kg_km2 = BIOMASS/0.0384) %>%
  sf::st_as_sf(coords = c("LON","LAT"), crs = wg_crs) %>%
  sf::st_join(epu) %>%
   # sf::st_join(ecodata::epu_sf) %>%
  sf::st_transform(crs = crs) %>%
  sfc_as_cols(names = c("longitude", "latitude")) %>%
  sf::st_drop_geometry()


dat <- nefsc_epu %>% 
  select(year = YEAR, survey = EPU, depth, logdepth,
         surface_temperature = SURFTEMP, temp = BOTTEMP, cpue_kg_km2, 
         latitude, longitude, yday, species = SVSPP) %>% 
  dplyr::filter(!is.na(temp), !is.na(depth), !is.na(longitude), !is.na(latitude),
                !is.na(logdepth), !is.na(yday)) 

# make mesh
spde <- try(make_mesh(dat, xy_cols = c("longitude", "latitude"),
                      cutoff = 20), silent = FALSE)

priors <- sdmTMBpriors(
  matern_s = pc_matern(
    range_gt = 5, range_prob = 0.05,
    sigma_lt = 25, sigma_prob = 0.05
  )
)
plot(spde)
#dat$fyear <- as.factor(dat$year)
mu_logdepth <- mean(dat$logdepth) #4.198351
sd_logdepth <- sd(dat$logdepth) #0.8858686
dat$logdepth <- (dat$logdepth - mu_logdepth) / sd_logdepth
dat$yday <- scale(dat$yday)[,1] # 198.3026,  sd =95.56088

fit <- sdmTMB(temp ~ s(yday) + s(logdepth),
              mesh = spde,
              time = "year",
              data = dat,
              spatial = "on",
              spatiotemporal = "ar1")
AIC(fit)
# make sure sanity checks pass
fit <- run_extra_optimization(fit, nlminb_loops = 1)
saveRDS(fit, "output/all_temp.rds")
fit <- readRDS(here::here("output/all_temp.rds"))

dat$resids <- residuals(fit)
qqnorm(dat$resids)
qqline(dat$resids)

# tidy(fit, conf.int = TRUE)

# Make predictions for each region
grid <- readRDS("grid_data/wc_grid.rds")
grid <- dplyr::rename(grid, longitude = X, latitude = Y) %>%
  dplyr::mutate(logdepth = log(-depth))

grid$lat_lon <- paste(grid$latitude, grid$longitude)

# scale the grid variables
grid$logdepth_orig <- grid$logdepth
# mu_logdepth = 5.060417, sd_logdepth = 0.689683
grid$logdepth <- (grid$logdepth - mu_logdepth) / sd_logdepth

pred_df <- expand.grid(
  lat_lon = unique(grid$lat_lon),
  year = 1990:2021
)
pred_df <- dplyr::left_join(pred_df, grid)
pred_df$yday <- (182 - 189.6727) / 35.23058 # Day 182 = July 1

# make a prediction for what this will be
pred_temp <- predict(fit, pred_df)
saveRDS(pred_temp, "output/wc_pred_temp.rds")

# Repeat for cells with depth < 250m
sub <- dplyr::filter(pred_df, abs(depth)<250)
pred_temp <- predict(fit, sub, return_tmb_object = TRUE)
index <- get_index(pred_temp, bias_correct = TRUE)
n_cells <- length(unique(sub$lat_lon))
index$est <- index$est / n_cells
index$se <- index$se / n_cells
index$lwr <- index$est - 1.96*index$se
index$upr <- index$est + 1.96*index$se
saveRDS(index, "output/temp_index_wc_250.rds")

sub <- dplyr::filter(pred_df, abs(depth)>250, abs(depth) < 500)
pred_temp <- predict(fit, sub, return_tmb_object = TRUE)
index <- get_index(pred_temp, bias_correct = TRUE)
n_cells <- length(unique(sub$lat_lon))
index$est <- index$est / n_cells
index$se <- index$se / n_cells
index$lwr <- index$est - 1.96*index$se
index$upr <- index$est + 1.96*index$se
saveRDS(index, "output/temp_index_wc_250_500.rds")

sub <- dplyr::filter(pred_df, abs(depth)>500, abs(depth) < 750)
pred_temp <- predict(fit, sub, return_tmb_object = TRUE)
index <- get_index(pred_temp, bias_correct = TRUE)
n_cells <- length(unique(sub$lat_lon))
index$est <- index$est / n_cells
index$se <- index$se / n_cells
index$lwr <- index$est - 1.96*index$se
index$upr <- index$est + 1.96*index$se
saveRDS(index, "output/temp_index_wc_500_750.rds")


# generate temp index for whole coast
pred_temp <- predict(fit, pred_df, return_tmb_object = TRUE)
index <- get_index(pred_temp, bias_correct = TRUE)
n_cells <- length(unique(pred_df$lat_lon))
index$est <- index$est / n_cells
index$se <- index$se / n_cells
index$lwr <- index$est - 1.96*index$se
index$upr <- index$est + 1.96*index$se
saveRDS(index, "output/temp_index_wc.rds")



