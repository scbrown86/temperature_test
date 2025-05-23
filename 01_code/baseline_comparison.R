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

# read in our downscaled data
brown <- list.files(path = "02_data/03_CHELSA_paleo/out",
                    pattern = "*1990.nc",
                    recursive = TRUE,
                    full.names = TRUE)
brown
brown <- sds(lapply(brown, function(x) {
  r <- rast(x)
  time(r) <- seq(as.Date("1900-01-16"), by = "month", l = nlyr(r))
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

chelsa_trace <- list(pr = crop(rast(chelsa_fil[1:12]), land),
                     tasmax = crop(rast(chelsa_fil[13:24]), land),
                     tasmin = crop(rast(chelsa_fil[25:36]), land))
chelsa_trace <- pblapply(chelsa_trace, function(i) {
  return(project(i, brown$pr[[1]], method = "average", use_gdal = TRUE, threads = TRUE))
})
sapply(chelsa_trace, function(i) plot(i[[1]]))
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
# koppen <- as.polygons(koppen, aggregate = TRUE)
koppen
plot(koppen)

# iterate through each dataset and create averages for periods of overlap
brown_1910_1989 <- pblapply(seq_along(brown), function(sd) {
  rsd <- brown[[sd]]
  rsd <- rsd[[which(format(time(rsd), "%Y") >= 1910)]]
  # annual average
  rsd <- tapp(rsd, "years", "mean")
  # time avg
  rsd <- app(rsd, mean, na.rm = TRUE)
  units(rsd) <- units(brown[[sd]])[1]
  varnames(rsd) <- varnames(brown[[sd]])[1]
  names(rsd) <- varnames(brown[[sd]])[1]
  time(rsd) <- 1950
  crs(rsd) <- "EPSG:4326"
  rsd
})
brown_1910_1989 <- rast(brown_1910_1989)
brown_1910_1989[[2:4]] <- setValues(x = brown_1910_1989[[2:4]],
                                    values = values(brown_1910_1989[[2:4]]) - 273.15)
plot(crop(mask(brown_1910_1989, land), land), fun = function() lines(land))

brown_1980_1989 <- pblapply(seq_along(brown), function(sd) {
  rsd <- brown[[sd]]
  rsd <- rsd[[which(format(time(rsd), "%Y") >= 1980)]]
  # annual average
  rsd <- tapp(rsd, "years", "mean")
  # time avg
  rsd <- app(rsd, mean, na.rm = TRUE)
  units(rsd) <- units(brown[[sd]])[1]
  varnames(rsd) <- varnames(brown[[sd]])[1]
  names(rsd) <- varnames(brown[[sd]])[1]
  time(rsd) <- 1985
  crs(rsd) <- "EPSG:4326"
  rsd
})
brown_1980_1989 <- rast(brown_1980_1989)
brown_1980_1989[[2:4]] <- setValues(x = brown_1980_1989[[2:4]],
                                    values = values(brown_1980_1989[[2:4]]) - 273.15)
plot(crop(mask(brown_1980_1989, land), land), fun = function() lines(land))

agcd_1910_1989 <- pblapply(seq_along(agcd), function(sd) {
  rsd <- agcd[[sd]]
  rsd <- rsd[[which(format(time(rsd), "%Y") >= 1910)]]
  # annual average
  rsd <- tapp(rsd, "years", "mean")
  # time avg
  rsd <- app(rsd, mean, na.rm = TRUE)
  units(rsd) <- c("kg m-2 s-1", "deg_C", "deg_C", "deg_C")[sd]
  varnames(rsd) <- c("pr", "tasmin", "tasmax", "tas")[sd]
  names(rsd) <- c("pr", "tasmin", "tasmax", "tas")[sd]
  time(rsd) <- 1950
  crs(rsd) <- "EPSG:4326"
  rsd
})
agcd_1910_1989 <- rast(agcd_1910_1989)
agcd_1910_1989
plot(crop(mask(agcd_1910_1989, land), land), fun = function() lines(land))

agcd_1980_1989 <- pblapply(seq_along(agcd), function(sd) {
  rsd <- agcd[[sd]]
  rsd <- rsd[[which(format(time(rsd), "%Y") >= 1980)]]
  # annual average
  rsd <- tapp(rsd, "years", "mean")
  # time avg
  rsd <- app(rsd, mean, na.rm = TRUE)
  units(rsd) <- c("kg m-2 s-1", "deg_C", "deg_C", "deg_C")[sd]
  varnames(rsd) <- c("pr", "tasmin", "tasmax", "tas")[sd]
  names(rsd) <- c("pr", "tasmin", "tasmax", "tas")[sd]
  time(rsd) <- 1985
  crs(rsd) <- "EPSG:4326"
  rsd
})
agcd_1980_1989 <- rast(agcd_1980_1989)
agcd_1980_1989
plot(crop(mask(agcd_1980_1989, land), land), fun = function() lines(land))

plot(
  crop(mask(brown_1910_1989[[c("tas", "tasmax", "tasmin")]], land), land) -
    crop(mask(agcd_1910_1989[[c("tas", "tasmax", "tasmin")]], land), land))

plot(
  crop(mask(brown_1910_1989[[c("pr")]], land), land) /
    crop(mask(agcd_1910_1989[[c("pr")]], land), land))

CHELSA_1980_1989 <- pblapply(seq_along(chelsa_12), function(sd) {
  rsd <- chelsa_12[[sd]]
  # annual average
  rsd <- tapp(rsd, "years", "mean")
  # time avg
  rsd <- app(rsd, mean, na.rm = TRUE)
  units(rsd) <- c("kg m-2 s-1", "deg_C", "deg_C", "deg_C")[sd]
  varnames(rsd) <- c("pr", "tas", "tasmax", "tasmin")[sd]
  names(rsd) <- c("pr", "tas", "tasmax", "tasmin")[sd]
  time(rsd) <- 1985
  crs(rsd) <- "EPSG:4326"
  rsd
})
CHELSA_1980_1989 <- rast(CHELSA_1980_1989)
CHELSA_1980_1989
plot(crop(mask(CHELSA_1980_1989, land), land), fun = function() lines(land))

CHELSA_Trace21_1950 <- pblapply(seq_along(chelsa_trace), function(sd) {
  rsd <- chelsa_trace[[sd]]
  # time avg (all 1950)
  rsd <- app(rsd, mean, na.rm = TRUE)
  units(rsd) <- c("kg m-2 s-1", "deg_C", "deg_C", "deg_C")[sd]
  varnames(rsd) <- c("pr", "tasmax", "tasmin", "tas")[sd]
  names(rsd) <- c("pr", "tasmax", "tasmin", "tas")[sd]
  time(rsd) <- 1950
  crs(rsd) <- "EPSG:4326"
  rsd
})
CHELSA_Trace21_1950 <- rast(CHELSA_Trace21_1950)
CHELSA_Trace21_1950
plot(crop(mask(CHELSA_Trace21_1950, land), land), fun = function() lines(land))

# save each raster
writeRaster(brown_1910_1989,
            filename = "03_comparisons/brown_1910_1989.tif")
writeRaster(brown_1980_1989,
            filename = "03_comparisons/brown_1980_1989.tif")
writeRaster(agcd_1910_1989,
            filename = "03_comparisons/agcd_1910_1989.tif")
writeRaster(agcd_1980_1989,
            filename = "03_comparisons/agcd_1980_1989.tif")
writeRaster(CHELSA_1980_1989,
            filename = "03_comparisons/CHELSA_1980_1989.tif")
writeRaster(CHELSA_Trace21_1950,
            filename = "03_comparisons/CHELSA_Trace21_1950.tif")

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

# # Split the data.table by sumstat
# split_list <- split(step_summaries_roll, by = "sumstat", keep.by = TRUE)
#
# # Apply frollmean to decadal data and rename the column for each subset
# split_list <- lapply(split_list, function(dt) {
#   dt[, year := year(yd)]
#   # Aggregate to yearly mean
#   yearly_dt <- copy(dt)[, .(mean_value = mean(value, na.rm = TRUE)),
#                         by = .(year, Koppen, climvar, MODEL, sumstat)]
#   stat_name <- yearly_dt[["sumstat"]][1]  # assuming one stat per dt
#   # Apply frollmean over yearly data
#   yearly_dt[, paste0("roll_", stat_name) := frollmean(mean_value, n = 10, align = "right", hasNA = FALSE),
#             by = .(Koppen, climvar, MODEL)]
#   return(yearly_dt)
# })
# split_list
#
# # Combine the results back into a single data.table
# step_summaries_roll <- rbindlist(split_list, use.names = TRUE, fill = TRUE)
# step_summaries_roll
# cols <- c("Koppen", "sumstat", "climvar", "MODEL")
# step_summaries_roll[, (cols) := lapply(.SD, as.factor), .SDcols = cols]
# summary(step_summaries_roll)
# # step_summaries_roll[, roll_value := frollmean(.SD, n = 120, align = "right",
# #                                               hasNA = FALSE),
# #                     by = c("Koppen", "sumstat", "climvar", "MODEL"),
# #                     .SDcols = "value"]
# # step_summaries_roll
# levels(step_summaries_roll$MODEL)
# levels(step_summaries_roll$climvar)
# levels(step_summaries_roll$Koppen)

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


p1 <- ggplot(data = step_summaries_roll[MODEL %in% c("BROWN", "AGCD"), ],
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

pdf(file = "03_comparisons/zonal_means.pdf",
    width = 14, height = 8, bg = "white",
    pagecentre = TRUE)
print(p1)
dev.off()

pdf(file = "03_comparisons/raster_means.pdf",
    paper = "a4", bg = "white",
    pagecentre = TRUE)
{p2 <- plot(crop(mask(c(agcd_1910_1989$pr, brown_1910_1989$pr), land), land),
           nr = 2, range = c(3.4e-06, 1.3e-04),
           main = c("BoM", "Brown"),
           fun = function() lines(land),
           plg = list(title = "rainfall (kg/m2/s)", side = 4,
                      font = 2, line = 2.5, cex = 0.8))

p3 <- plot(crop(mask(c(agcd_1910_1989$tas, brown_1910_1989$tas), land), land),
           nr = 2, range = c(-2, 30),
           main = c("BoM", "Brown"),
           fun = function() lines(land),
           plg = list(title = "tas (deg C)", side = 4,
                      font = 2, line = 2.5, cex = 0.8))

p4 <- plot(crop(mask(c(agcd_1910_1989$tasmin, brown_1910_1989$tasmin), land), land),
           nr = 2, range = c(-6, 27),
           main = c("BoM", "Brown"),
           fun = function() lines(land),
           plg = list(title = "tasmin (deg C)", side = 4,
                      font = 2, line = 2.5, cex = 0.8))

p5 <- plot(crop(mask(c(agcd_1910_1989$tasmax, brown_1910_1989$tasmax), land), land),
           nr = 2, range = c(2, 36),
           main = c("BoM", "Brown"),
           fun = function() lines(land),
           plg = list(title = "tasmax (deg C)", side = 4,
                      font = 2, line = 2.5, cex = 0.8))

p6 <- plot(crop(mask(c(CHELSA_Trace21_1950$pr, brown_1910_1989$pr), land), land),
           nr = 2, range = c(4.221234e-06, 1.121215e-04),
           main = c("CHELSA TraCE", "Brown"),
           fun = function() lines(land),
           plg = list(title = "rainfall (kg/m2/s)", side = 4,
                      font = 2, line = 2.5, cex = 0.8))

p7 <- plot(crop(mask(c(CHELSA_Trace21_1950$tas, brown_1910_1989$tas), land), land),
           nr = 2, range = c(-2, 30),
           main = c("CHELSA TraCE", "Brown"),
           fun = function() lines(land),
           plg = list(title = "tas (deg C)", side = 4,
                      font = 2, line = 2.5, cex = 0.8))

p8 <- plot(crop(mask(c(CHELSA_Trace21_1950$tasmin, brown_1910_1989$tasmin), land), land),
           nr = 2, range = c(-7, 29),
           main = c("CHELSA TraCE", "Brown"),
           fun = function() lines(land),
           plg = list(title = "tasmin (deg C)", side = 4,
                      font = 2, line = 2.5, cex = 0.8))

p9 <- plot(crop(mask(c(CHELSA_Trace21_1950$tasmax, brown_1910_1989$tasmax), land), land),
           nr = 2, range = c(2, 36),
           main = c("CHELSA TraCE", "Brown"),
           fun = function() lines(land),
           plg = list(title = "tasmax (deg C)", side = 4,
                      font = 2, line = 2.5, cex = 0.8))
}
dev.off()

qpdf::pdf_combine(input = c("03_comparisons/zonal_means.pdf",
                            "03_comparisons/raster_means.pdf"),
                  output = "03_comparisons/combined.pdf")

regions <- unique(step_summaries_roll[["Koppen"]])

delta_pr <- crop(mask(agcd_1910_1989$pr / brown_1910_1989$pr, land), land)
delta_pr
plot(delta_pr)

delta_tasmax <- crop(mask(agcd_1910_1989$tasmax - brown_1980_1989$tasmax, land), land)
delta_tasmax
plot(delta_tasmax)

delta_tasmin <- crop(mask(agcd_1910_1989$tasmin - brown_1980_1989$tasmin, land), land)
delta_tasmin
plot(delta_tasmin)

delta_tas <- crop(mask(agcd_1910_1989$tas - brown_1980_1989$tas, land), land)
delta_tas
plot(delta_tas)
