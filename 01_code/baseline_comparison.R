library(terra)
library(gtools)
library(pbapply)
library(exactextractr)
library(sf)
library(data.table)

setGDALconfig("GDAL_PAM_ENABLED", "FALSE")
terraOptions(memfrac = 0.85, memmax = 50)

land <- vect(rnaturalearthhires::countries10)
land <- crop(land, ext(106, 160, -50, 7))
land <- land[land$ADMIN == "Australia", ]
land <- aggregate(land)
plot(land)

plot_avg <- FALSE

# read in our downscaled data
brown <- list.files(path = "02_data/03_CHELSA_paleo/out",
                    pattern = "*1990.nc",
                    recursive = TRUE,
                    full.names = TRUE)
brown
brown <- sds(lapply(brown, function(x) {
  r <- rast(x)
  time(r) <- seq(as.Date("1900-01-06"), by = "month", l = nlyr(r))
  r
  }))
names(brown) <- c("pr", "tas", "tasmax", "tasmin")
brown

# plot the climatological average
if (plot_avg) {
  pr_avg <- (app(brown$pr, mean, na.rm = TRUE)*86400)*365 # mm/year
  tas_avg <- app(brown$tas, mean, na.rm = TRUE)-273.15
  tasmax_avg <- app(brown$tasmax, mean, na.rm = TRUE)-273.15
  tasmin_avg <- app(brown$tasmin, mean, na.rm = TRUE)-273.15
  par(mfrow = c(2,2))
  plot(pr_avg, fun = function() lines(land), range = c(0, 12500))
  plot(tas_avg, fun = function() lines(land), range = c(-10, 36))
  plot(tasmax_avg, fun = function() lines(land), range = c(-8, 38))
  plot(tasmin_avg, fun = function() lines(land), range = c(-13, 35))
}


# load in AGCD and project to area
source("01_code/00_functions/agcd_proc.R")
agcd <- lapply(c("precip", "tmin", "tmax"), agcd_proc,
               dir = "/mnt/Data/AusClim/data/raw/agcd",
               template = brown$pr[[1]], years = 1910:1989,
               proj_method = "average", type = ".nc$")

names(agcd) <- c("pr", "tasmin", "tasmax")
agcd$tas <- ((agcd$tasmax + agcd$tasmin)*0.5)*1
names(agcd$pr) <- paste0("pr_", 1:nlyr(agcd$pr))
names(agcd$tasmin) <- paste0("tasmin_", 1:nlyr(agcd$pr))
names(agcd$tasmax) <- paste0("tasmax_", 1:nlyr(agcd$pr))
names(agcd$tas) <- paste0("tas_", 1:nlyr(agcd$pr))
agcd

# Load in CHELSA baseline data
chelsa_12 <- mixedsort(list.files("/mnt/Data/CHELSA/v1.2/processed/fine",
                                  pattern = ".tif$", full.names = TRUE))
# Extract variable names for grouping
vars <- sapply(chelsa_12, function(x) sub(".*CHELSA_([a-z]+)_.*", "\\1", basename(x)))
chelsa_12 <- split(chelsa_12, vars) |>
  pblapply(function(x) {
    r <- project(rast(x), brown$pr, method = "average", use_gdal = TRUE, threads = TRUE)
    time(r) <- seq(as.Date("1980-01-16"), by = "month", l = nlyr(r))
    r*1
  })
chelsa_12

names(chelsa_12$pr) <- paste0("pr_", 1:nlyr(chelsa_12$pr))
names(chelsa_12$tasmin) <- paste0("tasmin_", 1:nlyr(chelsa_12$pr))
names(chelsa_12$tasmax) <- paste0("tasmax_", 1:nlyr(chelsa_12$pr))
names(chelsa_12$tas) <- paste0("tas_", 1:nlyr(chelsa_12$pr))

# Load in Koppen climate zones
koppen <- rast("/mnt/Data/AusClim/data/processed/koppen_zones_raster.tif")
koppen <- project(koppen, brown$pr, method = "near")
koppen <- mask(koppen, project(land, koppen), touches = TRUE)
has.colors(koppen)
kd <- data.frame(id = 1:6,
                 Koppen = c("Temperate", "Grassland", "Desert",
                            "Subtropical", "Tropical", "Equatorial"))
levels(koppen) <- kd
koppen
plot(koppen)
koppen <- as.polygons(koppen, aggregate = TRUE)
koppen
plot(koppen)

# iterate through each dataset and extract monthly averages and SD for each zone
source("01_code/00_functions/koppen_summary.R")

step_summaries <- koppen_summary(list(brown, agcd, chelsa_12))
names(step_summaries) <- c("brown", "agcd", "CHELSA")
step_summaries

plot(x = step_summaries$brown[Koppen == "Temperate" & sumstat == "mean" & climvar == "tas", ][["yd"]],
     y = step_summaries$brown[Koppen == "Temperate"& sumstat == "mean" & climvar == "tas", ][["value"]],
     type = "l", col = "#d95f02")
lines(x = step_summaries$agcd[Koppen == "Temperate" & sumstat == "mean" & climvar == "tas", ][["yd"]],
      y = step_summaries$agcd[Koppen == "Temperate"& sumstat == "mean" & climvar == "tas", ][["value"]]+273.15,
      col = "#1b9e77")
lines(x = step_summaries$CHELSA[Koppen == "Temperate" & sumstat == "mean" & climvar == "tas", ][["yd"]],
      y = step_summaries$CHELSA[Koppen == "Temperate"& sumstat == "mean" & climvar == "tas", ][["value"]]+273.15,
      col = "#7570b3")

plot(x = step_summaries$brown[Koppen == "Temperate" & sumstat == "mean" & climvar == "pr", ][["yd"]],
     y = step_summaries$brown[Koppen == "Temperate"& sumstat == "mean" & climvar == "pr", ][["value"]],
     type = "l", col = "#d95f02")
lines(x = step_summaries$agcd[Koppen == "Temperate" & sumstat == "mean" & climvar == "pr", ][["yd"]],
      y = step_summaries$agcd[Koppen == "Temperate"& sumstat == "mean" & climvar == "pr", ][["value"]],
      col = "#1b9e77")
lines(x = step_summaries$CHELSA[Koppen == "Temperate" & sumstat == "mean" & climvar == "pr", ][["yd"]],
      y = step_summaries$CHELSA[Koppen == "Temperate"& sumstat == "mean" & climvar == "pr", ][["value"]],
      col = "#7570b3")

step_summaries$brown <- copy(step_summaries$brown)[, MODEL := "BROWN"]
step_summaries$agcd <- copy(step_summaries$agcd)[, MODEL := "AGCD"]
step_summaries$CHELSA <- copy(step_summaries$CHELSA)[, MODEL := "CHELSA"]

step_summaries <- rbindlist(step_summaries)

step_summaries_roll <- copy(step_summaries)
summary(step_summaries_roll)
step_summaries_roll[, roll_value := frollmean(.SD, n = 120, align = "right",
                                              hasNA = FALSE),
                    by = c("Koppen", "sumstat", "climvar", "MODEL"),
                    .SDcols = "value"]
step_summaries_roll
# step_summaries <- pblapply(seq_along(grids), function(i) {
#   vars <- grids[[i]]
#   var_extract <- lapply(seq_along(vars), function(j,...) {
#     vj <- vars[[j]]
#     zt <- time(vj)
#     zt <- data.table(timestep = 1:length(zt), yd = zt)
#     vje <- setDT(exact_extract(vj, st_as_sf(koppen),
#                          fun = c("mean", "stdev", "coefficient_of_variation"),
#                          append_cols = "Koppen"))
#     vje <- melt(vje, id.vars = "Koppen")
#     vje[, c("sumstat", "rest") := tstrsplit(as.character(vje[["variable"]]),
#                                             "\\.(?=[^\\.]+_[0-9]+$)", perl = TRUE)]
#     vje[, c("climvar", "timestep") := tstrsplit(rest, "_")]
#     vje[, timestep := as.integer(timestep)][, rest := NULL][, variable := NULL]
#     vje <- merge(vje, zt, by = "timestep", all.x = TRUE, all.y = FALSE)
#     fout <- file.path(sprintf("scratch/var_extract_%s_%s.RDS", i, j))
#     saveRDS(vje, fout)
#     return(vje)
#   })
#   var_extract <- rbindlist(var_extract)
#   return(var_extract)
# })
# step_summaries
