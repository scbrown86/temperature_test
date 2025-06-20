library(terra)
library(gtools)
library(pbapply)
library(exactextractr)
library(sf)
library(data.table)
library(ggplot2)
library(qpdf)
library(plotrix)

setGDALconfig("GDAL_PAM_ENABLED", "FALSE")
terraOptions(memfrac = 0.85, memmax = 50)

land <- vect(rnaturalearthhires::countries10)
land <- crop(land, ext(106, 160, -50, 7))
land <- land[land$ADMIN == "Australia", ]
land <- aggregate(land)
plot(land)

plot_avg <- FALSE

# read in our downscaled and bias corrected data
brown <- list.files(path = "02_data/03_CHELSA_paleo/out",
                    pattern = ".nc",
                    recursive = FALSE,
                    full.names = TRUE)
brown
brown <- sds(lapply(brown, function(x) {
  r <- rast(x)
  r
  }))
names(brown) <- c("pr", "tas", "tasmax", "tasmin")
brown

# plot the climatological average
if (plot_avg) {
  pr_avg <- (app(brown$pr, mean, na.rm = TRUE)*86400)*365 # mm/year
  tas_avg <- app(brown$tas, mean, na.rm = TRUE)
  tasmax_avg <- app(brown$tasmax, mean, na.rm = TRUE)
  tasmin_avg <- app(brown$tasmin, mean, na.rm = TRUE)
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
chelsa_12 <- mixedsort(list.files("02_data/02_processed/CHELSA/",
                                  pattern = ".tif", full.names = TRUE))
chelsa_12 <- chelsa_12[!grepl("bspline", chelsa_12)]
# Extract variable names for grouping
vars <- sapply(chelsa_12, function(x) sub(".*CHELSA_([a-z]+)_.*", "\\1", basename(x)))
chelsa_12 <- split(chelsa_12, vars) |>
  pblapply(function(x) {
    r <- crop(rast(x), ext(brown$pr), snap = "out")
    r <- project(r, brown$pr,
                 method = "average", use_gdal = TRUE, threads = TRUE)
    time(r) <- seq(as.Date("1980-01-16"), by = "month", l = nlyr(r))
    r*1
  })
chelsa_12

names(chelsa_12$pr) <- paste0("pr_", 1:nlyr(chelsa_12$pr))
names(chelsa_12$tasmin) <- paste0("tasmin_", 1:nlyr(chelsa_12$pr))
names(chelsa_12$tasmax) <- paste0("tasmax_", 1:nlyr(chelsa_12$pr))
names(chelsa_12$tas) <- paste0("tas_", 1:nlyr(chelsa_12$pr))

# Read in the CHELSA-TraCE21k data for 1900-1989
chelsa_lut <- data.table(ID = 1:221,
                         year_ID = seq(20, -200, by = -1),
                         start_year = seq(1900, by = -100, l = 221),
                         end_year = c(1990, seq(1899, by = -100, l = 220)),
                         kBP = seq(0, by = -0.1, l = 221))
chelsa_lut[start_year >= 1850, ]
# List all .tif files in the folder
chelsa_fil <- mixedsort(list.files("/mnt/Data/CHELSA_Trace21/",
                                   pattern = ".tif$",
                                   full.names = TRUE))
head(chelsa_fil)
pattern <- "CHELSA_TraCE21k_[a-z]+_([1-9]|1[0-2])_20_V1\\.0\\.tif"
chelsa_fil <- chelsa_fil[grepl(pattern, chelsa_fil)]
chelsa_fil

chelsa_trace <- list(pr = crop(rast(chelsa_fil[1:12]), ext(brown$pr), snap = "out"),
                     tasmax = crop(rast(chelsa_fil[13:24]), ext(brown$pr), snap = "out"),
                     tasmin = crop(rast(chelsa_fil[25:36]), ext(brown$pr), snap = "out"))
chelsa_trace <- pblapply(chelsa_trace, function(i) {
  return(project(i, brown$pr[[1]], method = "average", use_gdal = TRUE, threads = TRUE))
})
sapply(chelsa_trace, function(i) plot(i[[1]]))
chelsa_trace
conv_rate <- (c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31) * 24) * 3600 # mm/month to kg/ms/2
chelsa_trace$pr <- chelsa_trace$pr / conv_rate
units(chelsa_trace$pr) <- "kg/m2/s"

# convert tasmax and tasmin to deg_C
chelsa_trace$tasmax <- (chelsa_trace$tasmax*0.1)-273.15
chelsa_trace$tasmin <- (chelsa_trace$tasmin*0.1)-273.15
units(chelsa_trace$tasmax) <- units(chelsa_trace$tasmin) <- "deg_C"
chelsa_trace$tas <- (chelsa_trace$tasmax + chelsa_trace$tasmin) * 0.5
units(chelsa_trace$tas) <- "deg_C"
time(chelsa_trace$pr) <- time(chelsa_trace$tasmax) <-
  time(chelsa_trace$tasmin) <- time(chelsa_trace$tas) <-
  seq(as.Date("1950-01-16"), by = "month", l = 12)
chelsa_trace

# Load in Koppen climate zones
koppen <- rast("02_data/01_inputs/koppen_zones_raster.tif")
koppen <- project(koppen, brown$pr, method = "near")
koppen <- mask(koppen, project(land, koppen), touches = TRUE)
has.colors(koppen)
kd <- data.frame(id = 1:6,
                 Koppen = c("Temperate", "Grassland", "Desert",
                            "Subtropical", "Tropical", "Equatorial"))
levels(koppen) <- kd
koppen
plot(koppen)
# koppen <- as.polygons(koppen, aggregate = TRUE)
koppen
plot(koppen)

# monthly climatological averages
brown_1910_1989_m <- pblapply(seq_along(brown), function(sd) {
  rsd <- brown[[sd]]
  rsd <- rsd[[which(format(time(rsd), "%Y") >= 1910)]]
  rsd <- tapp(rsd, "months", "mean")
  units(rsd) <- units(brown[[sd]])[1]
  if (units(brown[[sd]])[1] %in% c("K", "k", "kelvin")) {
    rsd <- setValues(rsd, values(rsd) - 273.15)
    units(rsd) <- "deg_C"
  }
  varnames(rsd) <- varnames(brown[[sd]])[1]
  names(rsd) <- paste0(month.abb, "_", varnames(brown[[sd]])[1])
  time(rsd, tstep = "months") <- seq(as.Date("1950-01-16"), by = "month", l = 12)
  crs(rsd) <- "EPSG:4326"
  rsd
})
brown_1910_1989_m <- rast(brown_1910_1989_m)
brown_1910_1989_m

brown_1980_1989_m <- pblapply(seq_along(brown), function(sd) {
  rsd <- brown[[sd]]
  rsd <- rsd[[which(format(time(rsd), "%Y") >= 1980)]]
  rsd <- tapp(rsd, "months", "mean")
  units(rsd) <- units(brown[[sd]])[1]
  if (units(brown[[sd]])[1] %in% c("K", "k", "kelvin")) {
    rsd <- setValues(rsd, values(rsd) - 273.15)
    units(rsd) <- "deg_C"
  }
  varnames(rsd) <- varnames(brown[[sd]])[1]
  names(rsd) <- paste0(month.abb, "_", varnames(brown[[sd]])[1])
  time(rsd, tstep = "months") <- seq(as.Date("1985-01-16"), by = "month", l = 12)
  crs(rsd) <- "EPSG:4326"
  rsd
})
brown_1980_1989_m <- rast(brown_1980_1989_m)
brown_1980_1989_m

agcd_1910_1989_m <- pblapply(seq_along(agcd), function(sd) {
  rsd <- agcd[[sd]]
  rsd <- rsd[[which(format(time(rsd), "%Y") >= 1910)]]
  rsd <- tapp(rsd, "months", "mean")
  units(rsd) <- c("kg/m2/s", "deg_C", "deg_C", "deg_C")[sd]
  varnames(rsd) <- c("pr", "tasmin", "tasmax", "tas")[sd]
  names(rsd) <- paste0(month.abb, "_", c("pr", "tasmin", "tasmax", "tas")[sd])
  time(rsd, tstep = "months") <- seq(as.Date("1950-01-16"), by = "month", l = 12)
  crs(rsd) <- "EPSG:4326"
  rsd
})
agcd_1910_1989_m <- rast(agcd_1910_1989_m)
agcd_1910_1989_m

agcd_1980_1989_m <- pblapply(seq_along(agcd), function(sd) {
  rsd <- agcd[[sd]]
  rsd <- rsd[[which(format(time(rsd), "%Y") >= 1980)]]
  rsd <- tapp(rsd, "months", "mean")
  units(rsd) <- c("kg/m2/s", "deg_C", "deg_C", "deg_C")[sd]
  varnames(rsd) <- c("pr", "tasmin", "tasmax", "tas")[sd]
  names(rsd) <- paste0(month.abb, "_", c("pr", "tasmin", "tasmax", "tas")[sd])
  time(rsd, tstep = "months") <- seq(as.Date("1985-01-16"), by = "month", l = 12)
  crs(rsd) <- "EPSG:4326"
  rsd
})
agcd_1980_1989_m <- rast(agcd_1980_1989_m)
agcd_1980_1989_m

CHELSA_1980_1989_m <- pblapply(seq_along(chelsa_12), function(sd) {
  rsd <- chelsa_12[[sd]]
  rsd <- tapp(rsd, "months", "mean")
  units(rsd) <- c("kg/m2/s", "deg_C", "deg_C", "deg_C")[sd]
  varnames(rsd) <- c("pr", "tas", "tasmax", "tasmin")[sd]
  names(rsd) <- paste0(month.abb, "_", c("pr", "tas", "tasmax", "tasmin")[sd])
  time(rsd, tstep = "months") <- seq(as.Date("1985-01-16"), by = "month", l = 12)
  crs(rsd) <- "EPSG:4326"
  rsd
})
CHELSA_1980_1989_m <- rast(CHELSA_1980_1989_m)
CHELSA_1980_1989_m

CHELSA_Trace21_1950_m <- pblapply(seq_along(chelsa_trace), function(sd) {
  rsd <- chelsa_trace[[sd]]
  rsd <- tapp(rsd, "months", "mean")
  units(rsd) <- c("kg/m2/s", "deg_C", "deg_C", "deg_C")[sd]
  varnames(rsd) <- c("pr", "tasmax", "tasmin", "tas")[sd]
  names(rsd) <- paste0(month.abb, "_", c("pr", "tasmax", "tasmin", "tas")[sd])
  time(rsd, tstep = "months") <- seq(as.Date("1950-01-16"), by = "month", l = 12)
  crs(rsd) <- "EPSG:4326"
  rsd
})
CHELSA_Trace21_1950_m <- rast(CHELSA_Trace21_1950_m)
CHELSA_Trace21_1950_m
par(mfrow = c(2, 2))
plot(crop(mask(brown_1980_1989_m[[1]],land), land), main = "brown")
plot(crop(mask(agcd_1980_1989_m[[1]], land),land), main = "AGCD")
plot(crop(mask(CHELSA_1980_1989_m[[1]],land),land), main = "CHELSA")
plot(crop(mask(CHELSA_Trace21_1950_m[[1]],land),land), main = "CHELSA-TraCE")

plot(crop(mask(brown_1980_1989_m[["Jan_tas"]],land), land), main = "brown")
plot(crop(mask(agcd_1980_1989_m[["Jan_tas"]], land),land), main = "AGCD")
plot(crop(mask(CHELSA_1980_1989_m[["Jan_tas"]],land),land), main = "CHELSA")
plot(crop(mask(CHELSA_Trace21_1950_m[["Jan_tas"]],land),land), main = "CHELSA-TraCE")

plot(crop(mask(brown_1980_1989_m[["Jan_tasmax"]],land), land), main = "brown")
plot(crop(mask(agcd_1980_1989_m[["Jan_tasmax"]], land),land), main = "AGCD")
plot(crop(mask(CHELSA_1980_1989_m[["Jan_tasmax"]],land),land), main = "CHELSA")
plot(crop(mask(CHELSA_Trace21_1950_m[["Jan_tasmax"]],land),land), main = "CHELSA-TraCE")

plot(crop(mask(brown_1980_1989_m[["Jan_tasmin"]],land), land), main = "brown")
plot(crop(mask(agcd_1980_1989_m[["Jan_tasmin"]], land),land), main = "AGCD")
plot(crop(mask(CHELSA_1980_1989_m[["Jan_tasmin"]],land),land), main = "CHELSA")
plot(crop(mask(CHELSA_Trace21_1950_m[["Jan_tasmin"]],land),land), main = "CHELSA-TraCE")

par(mfrow = c(1,1))

# save each raster
writeRaster(brown_1910_1989_m, overwrite = TRUE,
            filename = "03_comparisons/brown_1910_1989_monthly.tif")
writeRaster(brown_1980_1989_m,overwrite = TRUE,
            filename = "03_comparisons/brown_1980_1989_monthly.tif")
writeRaster(agcd_1910_1989_m,overwrite = TRUE,
            filename = "03_comparisons/agcd_1910_1989_monthly.tif")
writeRaster(agcd_1980_1989_m,overwrite = TRUE,
            filename = "03_comparisons/agcd_1980_1989_monthly.tif")
writeRaster(CHELSA_1980_1989_m,overwrite = TRUE,
            filename = "03_comparisons/CHELSA_1980_1989_monthly.tif")
writeRaster(CHELSA_Trace21_1950_m,overwrite = TRUE,
            filename = "03_comparisons/CHELSA_Trace21_1950_monthly.tif")

# iterate through each dataset and extract monthly averages and SD for each zone
source("01_code/00_functions/koppen_summary.R")

step_summaries <- koppen_summary(list(brown, agcd, chelsa_12, chelsa_trace))
names(step_summaries) <- c("brown", "agcd", "CHELSA", "CHELSA_trace")
step_summaries

saveRDS(step_summaries, "03_comparisons/step_summaries.RDS")

step_summaries <- readRDS("03_comparisons/step_summaries.RDS")

plot(x = step_summaries$brown[Koppen == "Temperate" & variable == "tas", ][["Year"]],
     y = step_summaries$brown[Koppen == "Temperate" & variable == "tas", ][["mean"]],
     type = "l", col = "#d95f02", ylim = c(12,16))
lines(x = step_summaries$agcd[Koppen == "Temperate" & variable == "tas", ][["Year"]],
      y = step_summaries$agcd[Koppen == "Temperate" & variable == "tas", ][["mean"]],
      col = "#1b9e77")
lines(x = step_summaries$CHELSA[Koppen == "Temperate" & variable == "tas", ][["Year"]],
      y = step_summaries$CHELSA[Koppen == "Temperate" & variable == "tas", ][["mean"]],
      col = "#7570b3")

plot(x = step_summaries$brown[Koppen == "Temperate" & variable == "pr", ][["Year"]],
     y = step_summaries$brown[Koppen == "Temperate" & variable == "pr", ][["mean"]],
     type = "l", col = "#d95f02")
lines(x = step_summaries$agcd[Koppen == "Temperate" & variable == "pr", ][["Year"]],
      y = step_summaries$agcd[Koppen == "Temperate" & variable == "pr", ][["mean"]],
      col = "#1b9e77")
lines(x = step_summaries$CHELSA[Koppen == "Temperate" & variable == "pr", ][["Year"]],
      y = step_summaries$CHELSA[Koppen == "Temperate" & variable == "pr", ][["mean"]],
      col = "#7570b3")

step_summaries$brown <- copy(step_summaries$brown)[, MODEL := "BROWN"]
step_summaries$agcd <- copy(step_summaries$agcd)[, MODEL := "AGCD"]
step_summaries$CHELSA <- copy(step_summaries$CHELSA)[, MODEL := "CHELSA"]
step_summaries$CHELSA_trace <- copy(step_summaries$CHELSA_trace)[, MODEL := "CHELSA_Trace"]

step_summaries <- rbindlist(step_summaries)

step_summaries_roll <- copy(step_summaries)
summary(step_summaries_roll)
step_summaries_roll

vars <- c("mean", "SD")
step_summaries_roll[, paste0("roll_", vars) := frollmean(.SD, n = 10, align = "right",
                                                         hasNA = FALSE),
                    by = c("Koppen", "variable", "MODEL"),
                    .SDcols = vars]
step_summaries_roll

step_summaries_roll[, MODEL := factor(MODEL, levels = c("BROWN", "AGCD", "CHELSA", "CHELSA_Trace"))]

set(step_summaries_roll, i = which(step_summaries_roll[["MODEL"]] == "CHELSA_Trace"),
    j = "roll_mean", value = step_summaries_roll[MODEL == "CHELSA_Trace", ][["mean"]])
set(step_summaries_roll, i = which(step_summaries_roll[["MODEL"]] == "CHELSA_Trace"),
    j = "roll_SD", value = step_summaries_roll[MODEL == "CHELSA_Trace", ][["SD"]])

# convert pr
step_summaries_roll[variable == "pr", `:=`(
  roll_mean = (roll_mean * 2629746)*12, #kg/m2/s to mm/year
  roll_SD = (roll_SD * 2629746)*12
)]
step_summaries_roll

facet_labels <- function(variable) {
  sapply(variable, function(v) {
    if (v == "pr") {
      "pr (mm/year)"
    } else {
      paste0(v, " (°C)")
    }
  })
}

# Define target years
years <- 1900:1989

# Get unique combinations of the grouping columns
unique_combos <- unique(copy(step_summaries_roll)[MODEL == "CHELSA_Trace", .(Koppen, variable, MODEL)])

# Create all combinations of those with the full year range
expanded <- CJ(Year = years, Koppen = unique_combos$Koppen,
               variable = unique_combos$variable, MODEL = unique_combos$MODEL, unique = TRUE)

# Merge the original data to the expanded data
# We'll keep only one value per group from the original and replicate it
merged <- merge(expanded, copy(step_summaries_roll)[MODEL == "CHELSA_Trace", ],
                by = c("Koppen", "variable", "MODEL"), all.x = TRUE)

# Fill in missing 'mean' and 'SD' using the first available value per group
merged[, `:=`(mean = mean[!is.na(mean)][1],
              SD = SD[!is.na(SD)][1]), by = .(Koppen, variable, MODEL)]


p1 <- ggplot(data = step_summaries_roll[MODEL %in% c("BROWN"), ],
       aes(x = Year, y = roll_mean,
           group = interaction(Koppen, MODEL),
           colour = MODEL, fill = MODEL)) +
  ggh4x::facet_grid2(factor(variable, levels = c("pr", "tas", "tasmin", "tasmax")) ~ Koppen,
                     scales = "free_y", independent = "y",
                     labeller = labeller(.rows = facet_labels)) +
  geom_ribbon(data = merged,
              inherit.aes = FALSE,
              aes(x = Year.x, y = roll_mean,
                  ymin = roll_mean - roll_SD,
                  ymax = roll_mean + roll_SD,
                  group = interaction(Koppen, MODEL),
                  fill = MODEL),
              colour = NA,
              alpha = 0.25) +
  geom_ribbon(inherit.aes = TRUE,
              colour = NA,
              aes(ymin = roll_mean - roll_SD,
                  ymax = roll_mean + roll_SD),
              alpha = 0.25) +
  geom_line(data = merged,
              inherit.aes = FALSE,
              aes(x = Year.x, y = roll_mean,
                  group = interaction(Koppen, MODEL),
                  colour = MODEL),
            linewidth = 0.5, show.legend = FALSE) +
  geom_line(linewidth = 0.5, show.legend = FALSE) +
  scale_colour_manual(values = c("BROWN" = "#1b9e77", "AGCD" = "#d95f02", "CHELSA_Trace" = "#7570b3")) +
  scale_fill_manual(values = c("BROWN" = "#1b9e77", "AGCD" = "#d95f02", "CHELSA_Trace" = "#7570b3"),
                    labels = c("BROWN" = "Brown", "AGCD" = "BoM", "CHELSA_Trace" = "CHELSA\nTraCE-21ka"),
                    guide = guide_legend(title = "Model",
                                         nrow = 1,
                                         theme = theme(
                                           legend.title.position = "top",
                                           legend.title = element_text(hjust = 0.5),
                                           legend.text.position = "bottom",
                                           legend.key.height = unit(0.5, "cm"),
                                           legend.key.width = unit(1, "cm")
                                         ))) +
  cowplot::theme_cowplot() +
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal") +
  labs(x = "Year", y = "mean")
p1

pdf(file = "03_comparisons/zonal_means_2model.pdf",
    width = 14, height = 8, bg = "white",
    pagecentre = TRUE)
print(p1)
dev.off()

# DELTA between CHELSA and Brown
## positive values mean Karger is higher
source("01_code/00_functions/raster_to_sds.R")
agcd_1910_1989 <- split_raster_by_variable(rast("03_comparisons/agcd_1910_1989_monthly.tif"))
agcd_1910_1989
agcd_1980_1989 <- split_raster_by_variable(rast("03_comparisons/agcd_1980_1989_monthly.tif"))
CHELSA_Trace21_1950 <- split_raster_by_variable(rast("03_comparisons/CHELSA_Trace21_1950_monthly.tif"))
brown_1910_1989 <- split_raster_by_variable(rast("03_comparisons/brown_1910_1989_monthly.tif"))
brown_1980_1989 <- split_raster_by_variable(rast("03_comparisons/brown_1980_1989_monthly.tif"))

pdf(file = "03_comparisons/raster_means.pdf",
    width = 9, height = 8, bg = "white",
    pagecentre = TRUE)
{p2 <- plot(crop(mask(c(app(agcd_1910_1989$pr, mean),
                        app(brown_1910_1989$pr, mean),
                        app(CHELSA_Trace21_1950$pr, mean)), land), land)*3.154e7,
            nr = 2, range = c(0, 2000),
            col = hcl.colors(100, "GnBu", rev = TRUE),
            main = c("BoM", "Brown", "CHELSA-TraCE"),
            smooth = FALSE, fill_range = TRUE,
            fun = function() lines(land),
            plg = list(title = "rainfall (mm/year)", side = 4,
                       font = 2, line = 2.5, cex = 0.8))

  p3 <- plot(crop(mask(c(app(agcd_1910_1989$tas, mean),
                         app(brown_1910_1989$tas, mean),
                         app(CHELSA_Trace21_1950$tas, mean)), land), land),
             nr = 2, range = c(5, 25),
             col = hcl.colors(100, "Batlow"),
             main = c("BoM", "Brown", "CHELSA-TraCE"),
             smooth = FALSE, fill_range = TRUE,
             fun = function() lines(land),
             plg = list(title = "tas (deg C)", side = 4,
                        font = 2, line = 2.5, cex = 0.8))

  p4 <- plot(crop(mask(c(app(agcd_1910_1989$tasmax, mean),
                         app(brown_1910_1989$tasmax, mean),
                         app(CHELSA_Trace21_1950$tasmax, mean)), land), land),
             nr = 2, range = c(10, 35),
             col = hcl.colors(100, "Batlow"),
             main = c("BoM", "Brown", "CHELSA-TraCE"),
             smooth = FALSE, fill_range = TRUE,
             fun = function() lines(land),
             plg = list(title = "tasmax (deg C)", side = 4,
                        font = 2, line = 2.5, cex = 0.8))

  p5 <- plot(crop(mask(c(app(agcd_1910_1989$tasmin, mean),
                         app(brown_1910_1989$tasmin, mean),
                         app(CHELSA_Trace21_1950$tasmin, mean)), land), land),
             nr = 2, range = c(0, 25),
             col = hcl.colors(100, "Batlow"),
             main = c("BoM", "Brown", "CHELSA-TraCE"),
             fun = function() lines(land),
             smooth = FALSE, fill_range = TRUE,
             plg = list(title = "tasmin (deg C)", side = 4,
                        font = 2, line = 2.5, cex = 0.8))

  # p6 <- plot(crop(mask(c(app(CHELSA_Trace21_1950$pr, mean), app(brown_1910_1989$pr, mean)), land), land)*3.154e7,
  #            nr = 2, range = c(0, 4000),
  #            col = hcl.colors(100, "GnBu", rev = TRUE),
  #            main = c("CHELSA TraCE", "Brown"),
  #            fun = function() lines(land),
  #            smooth = FALSE, fill_range = TRUE,
  #            plg = list(title = "rainfall (kg/m2/s)", side = 4,
  #                       font = 2, line = 2.5, cex = 0.8))
  #
  # p7 <- plot(crop(mask(c(app(CHELSA_Trace21_1950$tas, mean), app(brown_1910_1989$tas, mean)), land), land),
  #            nr = 2, range = c(3, 30),
  #            main = c("CHELSA TraCE", "Brown"),
  #            smooth = FALSE, fill_range = TRUE,
  #            fun = function() lines(land),
  #            plg = list(title = "tas (deg C)", side = 4,
  #                       font = 2, line = 2.5, cex = 0.8))
  #
  # p8 <- plot(crop(mask(c(app(CHELSA_Trace21_1950$tasmin, mean), app(brown_1910_1989$tasmin, mean)), land), land),
  #            nr = 2, range = c(-2, 29),
  #            main = c("CHELSA TraCE", "Brown"),
  #            smooth = FALSE, fill_range = TRUE,
  #            fun = function() lines(land),
  #            plg = list(title = "tasmin (deg C)", side = 4,
  #                       font = 2, line = 2.5, cex = 0.8))
  #
  # p9 <- plot(crop(mask(c(app(CHELSA_Trace21_1950$tasmax, mean), app(brown_1910_1989$tasmax, mean)), land), land),
  #            nr = 2, range = c(7, 36),
  #            main = c("CHELSA TraCE", "Brown"),
  #            smooth = FALSE, fill_range = TRUE,
  #            fun = function() lines(land),
  #            plg = list(title = "tasmax (deg C)", side = 4,
  #                       font = 2, line = 2.5, cex = 0.8))
}
dev.off()

qpdf::pdf_combine(input = c("03_comparisons/zonal_means_2model.pdf",
                            "03_comparisons/raster_means.pdf"),
                  output = "03_comparisons/combined.pdf")

regions <- unique(step_summaries_roll[["Koppen"]])

if (TRUE) {
  delta_pr <- c(
    crop(mask(app(CHELSA_Trace21_1950$pr, mean) / app(brown_1910_1989$pr, mean), land), land),
    crop(mask(app(CHELSA_Trace21_1950$pr, mean) / app(agcd_1910_1989$pr, mean), land), land),
    crop(mask(app(agcd_1980_1989$pr, mean) / app(brown_1980_1989$pr, mean), land), land))
  delta_pr # > 1 == Karger wetter
  plot(delta_pr)

  delta_tasmax <- c(
    crop(mask(app(CHELSA_Trace21_1950$tasmax, mean) - app(brown_1910_1989$tasmax, mean), land), land),
    crop(mask(app(CHELSA_Trace21_1950$tasmax, mean) - app(agcd_1910_1989$tasmax, mean), land), land),
    crop(mask(app(agcd_1980_1989$tasmax, mean) - app(brown_1980_1989$tasmax, mean), land), land))
  delta_tasmax
  hist(delta_tasmax)
  plot(delta_tasmax)

  delta_tasmin <- c(
    crop(mask(app(CHELSA_Trace21_1950$tasmin, mean) - app(brown_1910_1989$tasmin, mean), land), land),
    crop(mask(app(CHELSA_Trace21_1950$tasmin, mean) - app(agcd_1910_1989$tasmin, mean), land), land),
    crop(mask(app(agcd_1980_1989$tasmin, mean) - app(brown_1980_1989$tasmin, mean), land), land))
  hist(delta_tasmin)
  plot(delta_tasmin)

  delta_tas <- c(
    crop(mask(app(CHELSA_Trace21_1950$tas, mean) - app(brown_1910_1989$tas, mean), land), land),
    crop(mask(app(CHELSA_Trace21_1950$tas, mean) - app(agcd_1910_1989$tas, mean), land), land),
    crop(mask(app(agcd_1980_1989$tas, mean) - app(brown_1980_1989$tas, mean), land), land))
  delta_tas
  hist(delta_tas)
  plot(delta_tas)
}
pdf(file = "03_comparisons/comparisons_delta_2model.pdf",
    width = 10, height = 9, onefile = TRUE, bg = "white")
par(mfrow = c(2,2), mai = c(0.5,0.5,0.5,0.5), mar = c(0.5,0.5,0.5,0.5))
{plot(delta_pr, mar = c(5,0,2,0),
     buffer = TRUE, fill_range = TRUE,
     main = c("CHELSA-TraCE / Brown",
              "CHELSA-TraCE / BoM",
              "BoM / Brown"),
     col = hcl.colors(100, "GnBu", rev = TRUE),
     smooth = FALSE, box = TRUE, range = c(0.5, 2.5),
     plg = list(x = "bottom",
                cex = 1, bty = "n",
                size = c(0.75, 1),
                tics = "out", title = "precipitation delta"),
     fun = function() lines(land, col = "#000000"))
plot(delta_tas,
     range = c(-5, 5),
     main = c("CHELSA-TraCE - Brown",
              "CHELSA-TraCE - BoM",
              "BoM - Brown"),
     col = hcl.colors(100, "Spectral", rev = TRUE),
     fill_range = TRUE,
     buffer = TRUE,
     mar = c(5,0,2,0),
     smooth = FALSE, box = TRUE,
     plg = list(x = "bottom",
                cex = 1, bty = "n",
                size = c(0.75, 1),
                at = seq(-5, 5),
                tics = "out", title = "temperature delta"),
     fun = function() lines(land, col = "#000000"))
plot(delta_tasmax,
     range = c(-5, 5),
     main = c("CHELSA-TraCE - Brown",
              "CHELSA-TraCE - BoM",
              "BoM - Brown"),
     col = hcl.colors(100, "Spectral", rev = TRUE),
     fill_range = TRUE,
     buffer = TRUE,
     mar = c(5,0,2,0),
     smooth = FALSE, box = TRUE,
     plg = list(x = "bottom",
                cex = 1, bty = "n",
                size = c(0.75, 1),
                at = seq(-5, 5),
                tics = "out", title = "max temperature delta"),
     fun = function() lines(land, col = "#000000"))
plot(delta_tasmin,
     range = c(-5, 5),
     main = c("CHELSA-TraCE - Brown",
              "CHELSA-TraCE - BoM",
              "BoM - Brown"),
     col = hcl.colors(100, "Spectral", rev = TRUE),
     fill_range = TRUE,
     buffer = TRUE,
     mar = c(5,0,2,0),
     smooth = FALSE, box = TRUE,
     plg = list(x = "bottom",
                cex = 1, bty = "n",
                size = c(0.75, 1),
                at = seq(-5, 5),
                tics = "out", title = "min temperature delta"),
     fun = function() lines(land, col = "#000000"))
}
dev.off()

# common mask to mask all three datasets
comm <- c(agcd_1980_1989$pr[[1]], brown_1980_1989$pr[[1]], CHELSA_Trace21_1950$pr[[1]],
          agcd_1980_1989$tas[[1]], brown_1980_1989$tas[[1]], CHELSA_Trace21_1950$tas[[1]],
          agcd_1980_1989$tasmin[[1]], brown_1980_1989$tasmin[[1]], CHELSA_Trace21_1950$tasmin[[1]],
          agcd_1980_1989$tasmax[[1]], brown_1980_1989$tasmax[[1]], CHELSA_Trace21_1950$tasmax[[1]])
comm
comm_sum <- app(comm, function(i) sum(!is.na(i)))
comm_sum <- ifel(comm_sum == 12, 1, NA)
comm_mask <- mask(comm_sum, land)
comm_mask; plot(comm_mask)

koppen <- mask(koppen, comm_mask)

source("01_code/00_functions/spatrast_to_taylor.R")

#### ZONAL TAYLOR ####
{taylor_from_spatraster_zones(obs = comm$karger_tas,
                             mod = comm$brown_tas,
                             zones = koppen,
                             add = FALSE,
                             zone_names = c("Temperate", "Grassland", "Desert",
                                            "Subtropical", "Tropical", "Equatorial"),
                             col_palette = c("#1f78b4", "#ffff99",
                                             "#b15828", "#b2df8a",
                                             "#34a02c", "#cab2d6"),
                             pch = 17, pcex = 1.5,
                             sig_digits = 3,
                             main = NULL)
taylor_from_spatraster_zones(obs = comm$karger_tas,
                             mod = comm$brown_tas,
                             zones = NULL,
                             add = TRUE,
                             sig_digits = 3,
                             pch = 17, pcex = 1.5,
                             main = NULL, col = "black")
taylor_from_spatraster_zones(obs = comm$karger_tasmin,
                             mod = comm$brown_tasmin,
                             zones = koppen,
                             add = TRUE,
                             sig_digits = 3,
                             zone_names = c("Temperate", "Grassland", "Desert",
                                            "Subtropical", "Tropical", "Equatorial"),
                             col_palette = c("#1f78b4", "#ffff99",
                                             "#b15828", "#b2df8a",
                                             "#34a02c", "#cab2d6"),
                             pch = 18, pcex = 2,
                             main = NULL)
taylor_from_spatraster_zones(obs = comm$karger_tasmin,
                             mod = comm$brown_tasmin,
                             zones = NULL,
                             add = TRUE,
                             sig_digits = 3,
                             pch = 18, pcex = 2,
                             main = NULL, col = "black")
taylor_from_spatraster_zones(obs = comm$karger_tasmax,
                             mod = comm$brown_tasmax,
                             zones = koppen,
                             add = TRUE,
                             sig_digits = 3,
                             zone_names = c("Temperate", "Grassland", "Desert",
                                            "Subtropical", "Tropical", "Equatorial"),
                             col_palette = c("#1f78b4", "#ffff99",
                                             "#b15828", "#b2df8a",
                                             "#34a02c", "#cab2d6"),
                             pch = 15, pcex = 1.5,
                             main = NULL)
taylor_from_spatraster_zones(obs = comm$karger_tasmax,
                             mod = comm$brown_tasmax,
                             zones = NULL,
                             add = TRUE,
                             sig_digits = 3,
                             pch = 15, pcex = 1.5, col = "black",
                             main = NULL)
taylor_from_spatraster_zones(obs = comm$agcd_pr,
                             mod = comm$brown_pr,
                             zones = koppen,
                             add = TRUE,
                             sig_digits = 6,
                             zone_names = c("Temperate", "Grassland", "Desert",
                                            "Subtropical", "Tropical", "Equatorial"),
                             col_palette = c("#1f78b4", "#ffff99",
                                             "#b15828", "#b2df8a",
                                             "#34a02c", "#cab2d6"),
                             pch = 19, pcex = 1.5,
                             main = NULL)
taylor_from_spatraster_zones(obs = comm$agcd_pr,
                             mod = comm$brown_pr,
                             zones = NULL,
                             sig_digits = 6,
                             add = TRUE,
                             pch = 19, pcex = 1.5, col = "black",
                             main = NULL)
legend(x = 0.07, y = 1.85,
       legend = c("precipitation", "mean temperature",
                  "minimum temperature", "maximum temperature"),
       col = "black",
       pch = c(19, 17, 18, 15), box.lwd = 0,box.lty = 0,box.col = NA,
       pt.cex = c(1.5, 1.5, 2, 1.5),
       ncol = 1)
}
# Save current par settings
old_par <- par(no.readonly = TRUE)
# Add inset using par(fig = ...) and par(new = TRUE)
par(fig = c(0.55,1.0, 0.55, 1.0), new = TRUE)
# Plot the map: could be zones, obs, etc.
plot(crop(koppen, land), axes = FALSE,
     buffer = TRUE,
     legend = TRUE, box = FALSE,
     plg = list(x = 154, y = -10,
                cex = 1,
                title = "Koppen zone"))
par(old_par)


#### MONTHLY TAYLOR ####
# Define color palette
colours <- c(
  Rainfall = "#0072B2",
  AirTemp = "#D55E00",
  MinTemp = "#009E73",
  MaxTemp = "#CC79A7")
pdf(file = "03_comparisons/comparisons_taylor_updated.pdf",
    width = 8, height = 9, onefile = TRUE, bg = "white")
par(mai = c(0.5,0.5,0.5,0.5), mar = c(0.5,0.5,0.5,0.5))
{taylor_from_sds_monthly(obs_sds = agcd_1980_1989,
                        mod_sds = brown_1980_1989,
                        var_name = "tas",
                        col_palette = colours[2],
                        use_mask = comm_mask,
                        zones = NULL, add = FALSE,
                        sig_digits = 3, pch = 19,
                        main = "", ref.sd = TRUE,
                        sd.method = "population",
                        normalize = TRUE, mar = c(4,4,0,0))
taylor_from_sds_monthly(obs_sds = agcd_1980_1989,
                        mod_sds = brown_1980_1989,
                        var_name = "tasmin",
                        col_palette = colours[3],
                        use_mask = comm_mask,
                        zones = NULL, add = TRUE,
                        sig_digits = 3, pch = 19,
                        main = "", ref.sd = TRUE,
                        sd.method = "population",
                        normalize = TRUE, mar = c(4,4,0,0))
taylor_from_sds_monthly(obs_sds = agcd_1980_1989,
                        mod_sds = brown_1980_1989,
                        var_name = "tasmax",
                        col_palette = colours[4],
                        use_mask = comm_mask,
                        zones = NULL, add = TRUE,
                        sig_digits = 3, pch = 19,
                        main = "", ref.sd = TRUE,
                        sd.method = "population",
                        normalize = TRUE, mar = c(4,4,0,0))
taylor_from_sds_monthly(obs_sds = agcd_1980_1989,
                        mod_sds = brown_1980_1989,
                        var_name = "pr",
                        col_palette = colours[1],
                        use_mask = comm_mask,
                        zones = NULL, add = TRUE,
                        sig_digits = 3, pch = 19,
                        main = "", ref.sd = TRUE,
                        sd.method = "population",
                        normalize = TRUE, mar = c(4,4,0,0))

taylor_from_sds_monthly(obs_sds = CHELSA_Trace21_1950,
                        mod_sds = brown_1910_1989,
                        var_name = "tas",
                        col_palette = colours[2],
                        use_mask = comm_mask,
                        zones = NULL, add = TRUE,
                        sig_digits = 3, pch = 17,
                        main = "", ref.sd = TRUE,
                        sd.method = "population",
                        normalize = TRUE, mar = c(4,4,0,0))
taylor_from_sds_monthly(obs_sds = CHELSA_Trace21_1950,
                        mod_sds = brown_1910_1989,
                        var_name = "tasmin",
                        col_palette = colours[3],
                        use_mask = comm_mask,
                        zones = NULL, add = TRUE,
                        sig_digits = 3, pch = 17,
                        main = "", ref.sd = TRUE,
                        sd.method = "population",
                        normalize = TRUE, mar = c(4,4,0,0))
taylor_from_sds_monthly(obs_sds = CHELSA_Trace21_1950,
                        mod_sds = brown_1910_1989,
                        var_name = "tasmax",
                        col_palette = colours[4],
                        use_mask = comm_mask,
                        zones = NULL, add = TRUE,
                        sig_digits = 3, pch = 17,
                        main = "", ref.sd = TRUE,
                        sd.method = "population",
                        normalize = TRUE, mar = c(4,4,0,0))
taylor_from_sds_monthly(obs_sds = CHELSA_Trace21_1950,
                        mod_sds = brown_1910_1989,
                        var_name = "pr",
                        col_palette = colours[1],
                        use_mask = comm_mask,
                        zones = NULL, add = TRUE,
                        sig_digits = 3, pch = 17,
                        main = "", ref.sd = TRUE,
                        sd.method = "population",
                        normalize = TRUE, mar = c(4,4,0,0))
legend(x = 0.07, y = 1.2,
       legend = c("precipitation", "mean temperature",
                  "minimum temperature", "maximum temperature"),
       col = colours,
       pch = 19, box.lwd = 0,box.lty = 0,box.col = NA,
       pt.cex = 1.5,
       ncol = 1)
legend(x = 0.07, y = 1.0,
       legend = c("Aust. gridded climate data",
                  "CHELSA-TraCE21k"),
       col = "black",
       pch = c(19, 17), box.lwd = 0,box.lty = 0,box.col = NA,
       pt.cex = 1.5, ncol = 1)
mtext("*each dot represents a single month averaged over 1910-1989",
      side = 3, adj = 0.9, outer = FALSE, line = -2, cex = 0.75)
}
dev.off()
