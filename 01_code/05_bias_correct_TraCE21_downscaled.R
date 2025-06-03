library(terra)
setGDALconfig("GDAL_PAM_ENABLED", "FALSE") # don't write aux files!
terraOptions(memfrac = 0.85, memmax = 20)
library(gtools)
library(pbapply)
library(rnaturalearthhires)
library(qgisprocess)

# safe
safe_spline <- purrr::safely(qgisprocess::qgis_run_algorithm,
    otherwise = NULL, quiet = TRUE
)

land <- vect(rnaturalearthhires::countries10)
land <- crop(land, ext(106, 160, -45, 5))
land <- aggregate(land)
land
plot(land)

#### CHELSA ####

# load in CHELSA V1.2 data at original res
# calculate averages between 1980 and 1989

source("01_code/00_functions/chelsa_proc.R")

processed_chelsa <- lapply(c("prec", "tmin", "tmax", "tmean"), function(v, ...) {
    chelsa_proc(
        variable = v,
        mask = NULL,
        ymin = 1980, ymax = 1989,
        tras_ext = ext(105.0, 161.25, -52.5, 11.25),
        load_exist = TRUE,
        dir = "/mnt/Data/CHELSA/v1.2",
        outdir = "02_data/02_processed/CHELSA",
        cores = 10L
    )
})
names(processed_chelsa) <- c("pr", "tasmin", "tasmax", "tas")
str(processed_chelsa)

# Calculate the climatological averages
terraOptions(memfrac = 0.85, memmax = 50)
pr_avg <- rast(processed_chelsa$pr)
pr_avg
# plot(c(pr_avg[[1]], (pr_avg[[1]]*86400)*31),
#      nc = 1, main = c("kg/m2/s", "mm/month"),
#      fun = function() lines(land, col = "white", lwd = 1.5))

tmn_avg <- rast(processed_chelsa$tasmin)
tmx_avg <- rast(processed_chelsa$tasmax)
tas_avg <- rast(processed_chelsa$tas)
time(tmn_avg) <- time(tmx_avg) <- time(tas_avg) <- time(pr_avg)

idx <- format(time(pr_avg), "%b")

chelsa_climatologies <- c(
    "02_data/02_processed/CHELSA/CHELSA_pr_climatology.nc",
    "02_data/02_processed/CHELSA/CHELSA_tasmin_climatology.nc",
    "02_data/02_processed/CHELSA/CHELSA_tasmax_climatology.nc",
    "02_data/02_processed/CHELSA/CHELSA_tas_climatology.nc"
)

if (!all(file.exists(chelsa_climatologies))) {
    pr_avg <- tapp(pr_avg, idx, mean, cores = 12L)
    tmn_avg <- tapp(tmn_avg, idx, mean, cores = 12L)
    tmx_avg <- tapp(tmx_avg, idx, mean, cores = 12L)
    tas_avg <- tapp(tas_avg, idx, mean, cores = 12L)
    tas_avg
    # set time for variables
    time(pr_avg) <- time(tmn_avg) <- time(tas_avg) <- time(tmx_avg) <- seq(as.Date("1984-01-16"), by = "month", l = 12)
    # precip units
    units(pr_avg) <- "kg/m2/s"
    varnames(pr_avg) <- "precip"
    terra::longnames(pr_avg) <- "precipitation"
    # temperature units
    units(tmn_avg) <- units(tmx_avg) <- units(tas_avg) <- "deg_c"
    writeCDF(pr_avg, "02_data/02_processed/CHELSA/CHELSA_pr_climatology.nc",
        varname = "precip", longname = "rainfall", # compression = 6L,
        unit = "kg/m2/s", zname = "time", prec = "float",
        overwrite = TRUE
    )
    writeCDF(tmn_avg, "02_data/02_processed/CHELSA/CHELSA_tasmin_climatology.nc",
        varname = "tasmin",
        longname = "minimum near surface air temperature",
        # compression = 6L,
        unit = "deg_c", zname = "time", prec = "float",
        overwrite = TRUE
    )
    writeCDF(tmx_avg, "02_data/02_processed/CHELSA/CHELSA_tasmax_climatology.nc",
        varname = "tasmax",
        longname = "maximum near surface air temperature",
        # compression = 6L,
        unit = "deg_c", zname = "time", prec = "float",
        overwrite = TRUE
    )
    writeCDF(tas_avg, "02_data/02_processed/CHELSA/CHELSA_tas_climatology.nc",
        varname = "tas",
        longname = "mean near surface air temperature",
        # compression = 6L,
        unit = "deg_c", zname = "time", prec = "float",
        overwrite = TRUE
    )
}
rm(list = ls())

# safe
safe_spline <- purrr::safely(qgisprocess::qgis_run_algorithm,
    otherwise = NULL, quiet = TRUE
)

land <- vect(rnaturalearthhires::countries10)
land <- crop(land, ext(106, 160, -45, 5))
land <- aggregate(land)
land

pr_avg <- rast("02_data/02_processed/CHELSA/CHELSA_pr_climatology.nc")
pr_avg
pr_avg[[1]] * 1

chelsa_coarse_clim <- c(
    "02_data/02_processed/CHELSA/CHELSA_coarse_pr_climatology.nc",
    "02_data/02_processed/CHELSA/CHELSA_coarse_tas_climatology.nc",
    "02_data/02_processed/CHELSA/CHELSA_coarse_tasmax_climatology.nc",
    "02_data/02_processed/CHELSA/CHELSA_coarse_tasmin_climatology.nc"
)

if (!all(file.exists(chelsa_coarse_clim))) {
    #### PR_BIL ####
    # function for downscaling climatologies to common 0.5° grid
    pr_avg <- wrap(pr_avg)
    pr_bil <- pblapply(seq_len(12), function(i) {
        if (file.exists(sprintf("02_data/02_processed/CHELSA/pr_bspline_%02d.tif", i))) {
            return(terra::wrap(terra::rast(sprintf("02_data/02_processed/CHELSA/pr_bspline_%02d.tif", i))))
        }
        # precip constant to avoid div/0 errors later down the line
        prec_con <- 0.0001 / (86400 * c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))[i]
        tmp_r <- tempfile(
            pattern = sprintf("bspline_%02d", i), fileext = ".tif",
            tmpdir = tempdir()
        )
        bspline <- safe_spline(
            "sagang:multilevelbsplinefromgridpoints",
            GRID = unwrap(pr_avg)[[i]] * 1,
            "TARGET_USER_XMIN TARGET_USER_XMAX TARGET_USER_YMIN TARGET_USER_YMAX" = ext(105, 160, -50, 7),
            TARGET_USER_SIZE = 0.5,
            TARGET_USER_FITS = 1,
            METHOD = 0,
            DATATYPE = 0,
            EPSILON = 0.000100,
            LEVEL_MAX = 14,
            TARGET_OUT_GRID = tmp_r
        )
        if (is.null(bspline$error)) {
            b <- qgis_as_terra(bspline$result)
            # add a very small constant if interp rain <= 0
            b <- ifel(b <= 0, prec_con, b)
            writeRaster(b,
                filename = sprintf("02_data/02_processed/CHELSA/pr_bspline_%02d.tif", i),
                gdal = c("COMPRESS=LZW", "TFW=NO", "PREDICTOR=3")
            )
            b <- terra::wrap(b)
            return(b)
        } else {
            return(NULL)
        }
    }, cl = 12)
    pr_bil <- rast(lapply(pr_bil, function(f) unwrap(f)))
    time(pr_bil) <- seq(as.Date("1984-01-16"), by = "month", l = 12)
    units(pr_bil) <- "kg/m2/s"
    varnames(pr_bil) <- "precip"
    longnames(pr_bil) <- "precipitation"
    pr_bil * 1

    writeCDF(pr_bil,
        "02_data/02_processed/CHELSA/CHELSA_coarse_pr_climatology.nc",
        varname = "precip", longname = "rainfall", # compression = 6L,
        unit = "kg/m2/s", zname = "time", prec = "float",
        overwrite = TRUE
    )

    #### TSMX_BIL ####
    # function for downscaling climatologies to common 0.5° grid
    tmx_avg <- rast("02_data/02_processed/CHELSA/CHELSA_tasmax_climatology.nc")
    tmx_avg
    tmx_avg <- wrap(tmx_avg)

    tmx_bil <- pblapply(seq_len(12), function(i) {
        if (file.exists(sprintf("02_data/02_processed/CHELSA/tmx_bspline_%02d.tif", i))) {
            return(terra::wrap(terra::rast(sprintf("02_data/02_processed/CHELSA/tmx_bspline_%02d.tif", i))))
        }
        tmp_r <- tempfile(
            pattern = sprintf("bspline_tmx_%02d", i), fileext = ".tif",
            tmpdir = tempdir()
        )
        bspline <- safe_spline(
            "sagang:multilevelbsplinefromgridpoints",
            GRID = unwrap(tmx_avg)[[i]] * 1,
            "TARGET_USER_XMIN TARGET_USER_XMAX TARGET_USER_YMIN TARGET_USER_YMAX" = ext(105, 160, -50, 7),
            TARGET_USER_SIZE = 0.5,
            TARGET_USER_FITS = 1,
            METHOD = 0,
            DATATYPE = 0,
            EPSILON = 0.000100,
            LEVEL_MAX = 14,
            TARGET_OUT_GRID = tmp_r
        )
        if (is.null(bspline$error)) {
            b <- terra::wrap(qgis_as_terra(bspline$result))
            writeRaster(unwrap(b),
                filename = sprintf("02_data/02_processed/CHELSA/tmx_bspline_%02d.tif", i),
                gdal = c("COMPRESS=LZW", "TFW=NO", "PREDICTOR=3")
            )
            return(b)
        } else {
            return(NULL)
        }
    }, cl = 12)
    tmx_bil <- rast(lapply(tmx_bil, function(f) unwrap(f)))
    time(tmx_bil) <- seq(as.Date("1984-01-16"), by = "month", l = 12)
    units(tmx_bil) <- "deg_c"
    varnames(tmx_bil) <- "tasmax"
    longnames(tmx_bil) <- "maximum temperature at surface"
    tmx_bil * 1
    writeCDF(tmx_bil,
        "02_data/02_processed/CHELSA/CHELSA_coarse_tasmax_climatology.nc",
        varname = "tasmax", longname = "maximum temperature at surface", # compression = 6L,
        unit = "deg_c", zname = "time", prec = "float",
        overwrite = TRUE
    )

    #### TSMN_BIL ####
    # function for downscaling climatologies to common 0.5° grid
    tmn_avg <- rast("02_data/02_processed/CHELSA/CHELSA_tasmin_climatology.nc")
    tmn_avg <- wrap(tmn_avg)

    tmn_bil <- pblapply(seq_len(12), function(i) {
        if (file.exists(sprintf("02_data/02_processed/CHELSA/tmn_bspline_%02d.tif", i))) {
            return(terra::wrap(terra::rast(sprintf("02_data/02_processed/CHELSA/tmn_bspline_%02d.tif", i))))
        }
        tmp_r <- tempfile(
            pattern = sprintf("bspline_tmn_%02d_", i), fileext = ".tif",
            tmpdir = tempdir()
        )
        bspline <- safe_spline(
            "sagang:multilevelbsplinefromgridpoints",
            GRID = unwrap(tmn_avg)[[i]] * 1,
            "TARGET_USER_XMIN TARGET_USER_XMAX TARGET_USER_YMIN TARGET_USER_YMAX" = ext(105, 160, -50, 7),
            TARGET_USER_SIZE = 0.5,
            TARGET_USER_FITS = 1,
            METHOD = 0,
            DATATYPE = 0,
            EPSILON = 0.000100,
            LEVEL_MAX = 14,
            TARGET_OUT_GRID = tmp_r
        )
        if (is.null(bspline$error)) {
            b <- terra::wrap(qgis_as_terra(bspline$result))
            writeRaster(unwrap(b),
                filename = sprintf("02_data/02_processed/CHELSA/tmn_bspline_%02d.tif", i),
                gdal = c("COMPRESS=LZW", "TFW=NO", "PREDICTOR=3")
            )
            return(b)
        } else {
            return(NULL)
        }
    }, cl = 12)
    tmn_bil <- rast(lapply(tmn_bil, function(f) unwrap(f)))
    time(tmn_bil) <- seq(as.Date("1984-01-16"), by = "month", l = 12)
    units(tmn_bil) <- "deg_c"
    varnames(tmn_bil) <- "tasmin"
    longnames(tmn_bil) <- "minimum temperature at surface"
    tmn_bil * 1
    plot(tmn_bil)
    writeCDF(tmn_bil,
        "02_data/02_processed/CHELSA/CHELSA_coarse_tasmin_climatology.nc",
        varname = "tasmin", longname = "minimum temperature at surface", # compression = 6L,
        unit = "deg_c", zname = "time", prec = "float",
        overwrite = TRUE
    )

    #### TAS_BIL ####
    # function for downscaling climatologies to common 0.5° grid
    tas_avg <- rast("02_data/02_processed/CHELSA/CHELSA_tas_climatology.nc")
    tas_avg <- wrap(tas_avg)

    tas_bil <- pblapply(seq_len(12), function(i) {
        if (file.exists(sprintf("02_data/02_processed/CHELSA/tas_bspline_%02d.tif", i))) {
            return(terra::wrap(terra::rast(sprintf("02_data/02_processed/CHELSA/tas_bspline_%02d.tif", i))))
        }
        tmp_r <- tempfile(
            pattern = sprintf("bspline_tas_%02d", i), fileext = ".tif",
            tmpdir = tempdir()
        )
        bspline <- safe_spline(
            "sagang:multilevelbsplinefromgridpoints",
            GRID = unwrap(tas_avg)[[i]] * 1,
            "TARGET_USER_XMIN TARGET_USER_XMAX TARGET_USER_YMIN TARGET_USER_YMAX" = ext(105, 160, -50, 7),
            TARGET_USER_SIZE = 0.5,
            TARGET_USER_FITS = 1,
            METHOD = 0,
            DATATYPE = 0,
            EPSILON = 0.000100,
            LEVEL_MAX = 14,
            TARGET_OUT_GRID = tmp_r
        )
        if (is.null(bspline$error)) {
            b <- terra::wrap(qgis_as_terra(bspline$result))
            writeRaster(unwrap(b),
                filename = sprintf("02_data/02_processed/CHELSA/tas_bspline_%02d.tif", i),
                gdal = c("COMPRESS=LZW", "TFW=NO", "PREDICTOR=3")
            )
            return(b)
        } else {
            return(NULL)
        }
    }, cl = 12)
    tas_bil <- rast(lapply(tas_bil, function(f) unwrap(f)))
    time(tas_bil) <- seq(as.Date("1984-01-16"), by = "month", l = 12)
    units(tas_bil) <- "deg_c"
    varnames(tas_bil) <- "tas"
    longnames(tas_bil) <- "mean temperature at surface"
    tas_bil * 1
    plot(tas_bil)
    writeCDF(tas_bil,
        "02_data/02_processed/CHELSA/CHELSA_coarse_tas_climatology.nc",
        varname = "tas", longname = "mean temperature at surface", # compression = 6L,
        unit = "deg_c", zname = "time", prec = "float",
        overwrite = TRUE
    )
} else {
    tas_bil <- terra::rast("02_data/02_processed/CHELSA/CHELSA_coarse_tas_climatology.nc")
    tmn_bil <- terra::rast("02_data/02_processed/CHELSA/CHELSA_coarse_tasmin_climatology.nc")
    tmx_bil <- terra::rast("02_data/02_processed/CHELSA/CHELSA_coarse_tasmax_climatology.nc")
    pr_bil <- terra::rast("02_data/02_processed/CHELSA/CHELSA_coarse_pr_climatology.nc")
}

#### TRACE CLIMATOLOGIES ####
trace_pr <- rast("02_data/02_processed/pr.nc")
trace_pr <- trace_pr[[which(time(trace_pr) >= as.Date("1980-01-01"))]]
trace_pr
plot(trace_pr[[1]], fun = function() lines(land, col = "#FFFFFF"))

trace_tmn <- rast("02_data/02_processed/tasmin.nc")
trace_tmn <- trace_tmn[[which(time(trace_tmn) >= as.Date("1980-01-01"))]]
trace_tmn <- trace_tmn - 273.15
trace_tmn
plot(trace_tmn[[1]], fun = function() lines(land, col = "#FFFFFF"))

trace_tmx <- rast("02_data/02_processed/tasmax.nc")
trace_tmx <- trace_tmx[[which(time(trace_tmx) >= as.Date("1980-01-01"))]]
trace_tmx <- trace_tmx - 273.15
trace_tmx
plot(trace_tmx[[1]], fun = function() lines(land, col = "#FFFFFF"))

trace_tas <- rast("02_data/02_processed/tas.nc")
trace_tas <- trace_tas[[which(time(trace_tas) >= as.Date("1980-01-01"))]]
trace_tas <- trace_tas - 273.15
trace_tas
plot(trace_tas[[1]], fun = function() lines(land, col = "#FFFFFF"))

idx <- format(time(trace_pr), "%b")

trace_pr <- tapp(trace_pr, idx, mean)
trace_pr

trace_tmn <- tapp(trace_tmn, idx, mean)
trace_tmx <- tapp(trace_tmx, idx, mean)
trace_tas <- tapp(trace_tas, idx, mean)
trace_tas
trace_tmn
trace_tmx

time(trace_pr) <- time(trace_tmn) <- time(trace_tmx) <- time(trace_tas) <- seq(as.Date("1984-01-16"), by = "month", l = 12)
units(trace_pr) <- "kg/m2/s"
varnames(trace_pr) <- "precip"
terra::longnames(trace_pr) <- "precipitation"

units(trace_tmn) <- units(trace_tmx) <- units(trace_tas) <- "degC"
varnames(trace_tmn) <- varnames(trace_tmx) <- varnames(trace_tas) <- "tas"

trace_pr
trace_tas
trace_tmn
trace_tmx

writeCDF(trace_pr, "02_data/02_processed/TRACE/TRACE_pr_climatology.nc",
    varname = "pr", longname = "rainfall", # compression = 6L,
    unit = "kg/m2/s", zname = "time", prec = "float",
    overwrite = TRUE
)
writeCDF(trace_tmn, "02_data/02_processed/TRACE/TRACE_tasmin_climatology.nc",
    varname = "tas_min",
    longname = "minimum near surface air temperature",
    # compression = 6L,
    unit = "deg_c", zname = "time", prec = "float",
    overwrite = TRUE
)
writeCDF(trace_tmx, "02_data/02_processed/TRACE/TRACE_tasmax_climatology.nc",
    varname = "tas_max",
    longname = "maximum near surface air temperature",
    # compression = 6L,
    unit = "deg_c", zname = "time", prec = "float",
    overwrite = TRUE
)
writeCDF(trace_tas, "02_data/02_processed/TRACE/TRACE_tas_climatology.nc",
    varname = "tas",
    longname = "mean near surface air temperature",
    # compression = 6L,
    unit = "deg_c", zname = "time", prec = "float",
    overwrite = TRUE
)

#### TRACE DOWNSCALE ####

# functions for downscaling climatologies to common 0.5° grid
pr_trace <- rast("02_data/02_processed/TRACE/TRACE_pr_climatology.nc")
pr_trace
plot(pr_trace[[1]])
par(mfrow = c(2, 1))
plot(pr_trace[[1]] * 86400 * 31,
    fun = function() lines(land), range = c(0, 1100),
    col = hcl.colors(100, "Batlow"), main = "TraCE21"
)
plot(rast("02_data/02_processed/CHELSA/CHELSA_coarse_pr_climatology.nc")[[1]] * 86400 * 31,
    fun = function() lines(land), range = c(0, 1100), main = "CHELSA",
    col = hcl.colors(100, "Batlow")
)
graphics.off()

pr_trace <- wrap(pr_trace)

pr_Tbil <- pblapply(seq_len(12), function(i) {
    if (file.exists(sprintf("/02_data/02_processed/TRACE/pr_bspline_%02d.tif", i))) {
        return(terra::wrap(terra::rast(sprintf("02_data/02_processed/TRACE/pr_bspline_%02d.tif", i))))
    }
    prec_con <- 0.1 / (86400 * c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))[i]
    tmp_r <- tempfile(
        pattern = sprintf("bspline_%02d_", i), fileext = ".tif",
        tmpdir = tempdir()
    )
    ingrid <- unwrap(pr_trace)[[i]] * 1
    # convert to mm/month
    dmon <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)[i]
    ingrid <- ingrid * (86400 * dmon)
    ingrid <- ifel(ingrid < 5, 5, ingrid)
    bspline <- safe_spline(
        "sagang:multilevelbsplinefromgridpoints",
        GRID = ingrid,
        "TARGET_USER_XMIN TARGET_USER_XMAX TARGET_USER_YMIN TARGET_USER_YMAX" = ext(105, 160, -50, 7),
        TARGET_USER_SIZE = 0.5,
        TARGET_USER_FITS = 1,
        METHOD = 0,
        DATATYPE = 0,
        EPSILON = 0.000100,
        LEVEL_MAX = 14,
        TARGET_OUT_GRID = tmp_r
    )
    if (is.null(bspline$error)) {
        b <- qgis_as_terra(bspline$result)
        # b <- ifel(b <= 0, prec_con*(86400*dmon), b) # add a very small constant if interp rain < 0
        b <- ifel(b < 1, 1, b)
        # convert to back to kg/m2/2
        b <- b / (86400 * dmon)
        writeRaster(b,
            filename = sprintf("02_data/02_processed/TRACE/pr_bspline_%02d.tif", i),
            gdal = c("COMPRESS=LZW", "TFW=NO", "PREDICTOR=3")
        )
        b <- terra::wrap(b)
        return(b)
    } else {
        return(NULL)
    }
}, cl = 12)
pr_Tbil <- rast(lapply(pr_Tbil, function(f) unwrap(f)))
time(pr_Tbil) <- seq(as.Date("1984-01-16"), by = "month", l = 12)
units(pr_Tbil) <- "kg/m2/s"
varnames(pr_Tbil) <- "pr"
longnames(pr_Tbil) <- "precipitation"
pr_Tbil * 1
plot(pr_Tbil)


c(
    pr_Tbil[[1]] * 86400 * 31,
    rast("02_data/02_processed/CHELSA/CHELSA_coarse_pr_climatology.nc")[[1]] * 86400 * 31
)

plot(
    c(
        pr_Tbil[[1]] * 86400 * 31,
        rast("02_data/02_processed/CHELSA/CHELSA_coarse_pr_climatology.nc")[[1]] * 86400 * 31
    ),
    main = c("TraCE downscaled", "CHELSA"),
    fun = function() lines(land), range = c(0, 1000),
    col = hcl.colors(100, "Batlow")
)

minmax(mask(
    rast("02_data/02_processed/CHELSA/CHELSA_coarse_pr_climatology.nc") / pr_Tbil,
    land
))
minmax(rast("02_data/02_processed/CHELSA/CHELSA_coarse_pr_climatology.nc") / pr_Tbil)

plot(rast("02_data/02_processed/CHELSA/CHELSA_coarse_pr_climatology.nc") / pr_Tbil,
    main = "delta", fun = function() lines(land),
    range = c(0, 20),
    col = hcl.colors(100, "Batlow")
)

writeCDF(pr_Tbil,
    "02_data/02_processed/TRACE/TRACE_fine_pr_climatology.nc",
    varname = "pr", longname = "rainfall", # compression = 6L,
    unit = "kg/m2/s", zname = "time", prec = "float",
    overwrite = TRUE
)

## TMX TRACE ##
tmx_trace <- rast("02_data/02_processed/TRACE/TRACE_tasmax_climatology.nc") * 1
tmx_trace
plot(tmx_trace[[1]])
par(mfrow = c(2, 1))
plot(tmx_trace[[1]],
    fun = function() lines(land), range = c(0, 60),
    col = hcl.colors(100, "Batlow")
)
plot(tmx_bil[[1]],
    fun = function() lines(land), range = c(0, 60),
    col = hcl.colors(100, "Batlow")
)
dev.off()

tmx_trace <- wrap(tmx_trace)

tmx_Tbil <- pblapply(seq_len(12), function(i) {
    if (file.exists(sprintf("02_data/02_processed/TRACE/tmx_bspline_%02d.tif", i))) {
        return(terra::wrap(terra::rast(sprintf("02_data/02_processed/TRACE/tmx_bspline_%02d.tif", i))))
    }
    tmp_r <- tempfile(
        pattern = sprintf("bspline_tmx_%02d_", i), fileext = ".tif",
        tmpdir = tempdir()
    )
    ingrid <- unwrap(tmx_trace)[[i]] * 1
    bspline <- safe_spline(
        "sagang:multilevelbsplinefromgridpoints",
        GRID = ingrid,
        "TARGET_USER_XMIN TARGET_USER_XMAX TARGET_USER_YMIN TARGET_USER_YMAX" = ext(105, 160, -50, 7),
        TARGET_USER_SIZE = 0.5,
        TARGET_USER_FITS = 1,
        METHOD = 0,
        DATATYPE = 0,
        EPSILON = 0.000100,
        LEVEL_MAX = 14,
        TARGET_OUT_GRID = tmp_r
    )
    if (is.null(bspline$error)) {
        # plot(qgis_as_terra(bspline$result), fun = function() lines(land))
        b <- terra::wrap(qgis_as_terra(bspline$result))
        writeRaster(unwrap(b),
            filename = sprintf("02_data/02_processed/TRACE/tmx_bspline_%02d.tif", i),
            gdal = c("COMPRESS=LZW", "TFW=NO", "PREDICTOR=3")
        )
        return(b)
    } else {
        return(NULL)
    }
}, cl = 12)
tmx_Tbil <- rast(lapply(tmx_Tbil, function(f) unwrap(f)))
time(tmx_Tbil) <- seq(as.Date("1984-01-16"), by = "month", l = 12)
units(tmx_Tbil) <- "degC"
varnames(tmx_Tbil) <- "tasmax"
longnames(tmx_Tbil) <- "maximum temperature at surface"
tmx_Tbil * 1
plot(tmx_Tbil[[1]], fun = function() lines(land))
c(app(tmx_Tbil, mean), app(tmx_bil, mean))
plot(c(app(tmx_Tbil, mean), app(tmx_bil, mean)),
    fun = function() lines(land), range = c(0, 46),
    col = hcl.colors(100, "Batlow")
)

writeCDF(tmx_Tbil,
    "02_data/02_processed/TRACE/TRACE_fine_tasmax_climatology.nc",
    varname = "tasmax", longname = "maximum temperature at surface", # compression = 6L,
    unit = "degC", zname = "time", prec = "float",
    overwrite = TRUE
)

## TMN TRACE ##
tmn_trace <- rast("02_data/02_processed/TRACE/TRACE_tasmin_climatology.nc") * 1
tmn_trace
plot(tmn_trace[[1]])
plot(tmx_trace[[1]] - tmn_trace[[1]])

par(mfrow = c(2, 1))
plot(tmn_trace[[1]],
    fun = function() lines(land), range = c(0, 30),
    col = hcl.colors(100, "Batlow")
)
plot(tmn_bil[[1]],
    fun = function() lines(land), range = c(0, 30),
    col = hcl.colors(100, "Batlow")
)
dev.off()

tmn_trace <- wrap(tmn_trace)

tmn_Tbil <- pblapply(seq_len(12), function(i) {
    if (file.exists(sprintf("02_data/02_processed/TRACE/tmn_bspline_%02d.tif", i))) {
        return(terra::wrap(terra::rast(sprintf("02_data/02_processed/TRACE/tmn_bspline_%02d.tif", i))))
    }
    tmp_r <- tempfile(
        pattern = sprintf("bspline_tmn_%02d_", i), fileext = ".tif",
        tmpdir = tempdir()
    )
    ingrid <- unwrap(tmn_trace)[[i]] * 1
    bspline <- safe_spline(
        "sagang:multilevelbsplinefromgridpoints",
        GRID = ingrid,
        "TARGET_USER_XMIN TARGET_USER_XMAX TARGET_USER_YMIN TARGET_USER_YMAX" = ext(105, 160, -50, 7),
        TARGET_USER_SIZE = 0.5,
        TARGET_USER_FITS = 1,
        METHOD = 0,
        DATATYPE = 0,
        EPSILON = 0.000100,
        LEVEL_MAX = 14,
        TARGET_OUT_GRID = tmp_r
    )
    if (is.null(bspline$error)) {
        b <- terra::wrap(qgis_as_terra(bspline$result))
        writeRaster(unwrap(b),
            filename = sprintf("02_data/02_processed/TRACE/tmn_bspline_%02d.tif", i),
            gdal = c("COMPRESS=LZW", "TFW=NO", "PREDICTOR=3")
        )
        return(b)
    } else {
        return(NULL)
    }
}, cl = 12)
tmn_Tbil <- rast(lapply(tmn_Tbil, function(f) unwrap(f)))
time(tmn_Tbil) <- seq(as.Date("1984-01-16"), by = "month", l = 12)
units(tmn_Tbil) <- "degC"
varnames(tmn_Tbil) <- "tasmin"
longnames(tmx_Tbil) <- "minimum temperature at surface"
tmn_Tbil * 1
plot(tmn_Tbil[[1]], fun = function() lines(land))
c(app(tmn_Tbil, mean), app(tmn_bil, mean))
plot(c(app(tmn_Tbil, mean), app(tmn_bil, mean)),
    fun = function() lines(land), range = c(0, 30),
    col = hcl.colors(100, "Batlow")
)

writeCDF(tmn_Tbil,
    "02_data/02_processed/TRACE/TRACE_fine_tasmin_climatology.nc",
    varname = "tasmin", longname = "minimum temperature at surface", # compression = 6L,
    unit = "degC", zname = "time", prec = "float",
    overwrite = TRUE
)

## TAS TRACE ##
tas_trace <- rast("02_data/02_processed/TRACE/TRACE_tas_climatology.nc") * 1
tas_trace
plot(tas_trace[[1]])
par(mfrow = c(2, 1))
plot(tas_trace[[1]],
    fun = function() lines(land), range = c(0, 35),
    col = hcl.colors(100, "Batlow")
)
plot(tas_bil[[1]],
    fun = function() lines(land), range = c(0, 35),
    col = hcl.colors(100, "Batlow")
)
dev.off()

tas_trace <- wrap(tas_trace)

tas_Tbil <- pblapply(seq_len(12), function(i) {
    if (file.exists(sprintf("02_data/02_processed/TRACE/tas_bspline_%02d.tif", i))) {
        return(terra::wrap(terra::rast(sprintf("02_data/02_processed/TRACE/tas_bspline_%02d.tif", i))))
    }
    tmp_r <- tempfile(
        pattern = sprintf("bspline_tas_%02d", i), fileext = ".tif",
        tmpdir = tempdir()
    )
    ingrid <- unwrap(tas_trace)[[i]] * 1
    bspline <- safe_spline(
        "sagang:multilevelbsplinefromgridpoints",
        GRID = ingrid,
        "TARGET_USER_XMIN TARGET_USER_XMAX TARGET_USER_YMIN TARGET_USER_YMAX" = ext(105, 160, -50, 7),
        TARGET_USER_SIZE = 0.5,
        TARGET_USER_FITS = 1,
        METHOD = 0,
        DATATYPE = 0,
        EPSILON = 0.000100,
        LEVEL_MAX = 14,
        TARGET_OUT_GRID = tmp_r
    )
    if (is.null(bspline$error)) {
        b <- terra::wrap(qgis_as_terra(bspline$result))
        writeRaster(unwrap(b),
            filename = sprintf("02_data/02_processed/TRACE/tas_bspline_%02d.tif", i),
            gdal = c("COMPRESS=LZW", "TFW=NO", "PREDICTOR=3")
        )
        return(b)
    } else {
        return(NULL)
    }
}, cl = 12)
tas_Tbil <- rast(lapply(tas_Tbil, function(f) unwrap(f)))
time(tas_Tbil) <- seq(as.Date("1984-01-16"), by = "month", l = 12)
units(tas_Tbil) <- "degC"
varnames(tas_Tbil) <- "tas"
longnames(tas_Tbil) <- "mean temperature at surface"
tas_Tbil * 1
plot(tas_Tbil[[1]], fun = function() lines(land))
c(app(tas_Tbil, mean), app(tas_bil, mean))
plot(c(app(tas_Tbil, mean), app(tas_bil, mean)),
    fun = function() lines(land), range = c(0, 30),
    col = hcl.colors(100, "Batlow")
)

writeCDF(tas_Tbil,
    "02_data/02_processed/TRACE/TRACE_fine_tas_climatology.nc",
    varname = "tas", longname = "mean temperature at surface", # compression = 6L,
    unit = "deg_c", zname = "time", prec = "float",
    overwrite = TRUE
)

#### DELTA ####
pr_bil <- rast("02_data/02_processed/CHELSA/CHELSA_coarse_pr_climatology.nc") * 1
tmx_bil <- rast("02_data/02_processed/CHELSA/CHELSA_coarse_tasmax_climatology.nc") * 1
tmn_bil <- rast("02_data/02_processed/CHELSA/CHELSA_coarse_tasmin_climatology.nc") * 1
tas_bil <- rast("02_data/02_processed/CHELSA/CHELSA_coarse_tas_climatology.nc") * 1

pr_bil
tmx_bil
tmn_bil
tas_bil

plot(pr_bil / pr_Tbil)
plot((tmx_bil + 273.15) - (tmx_Tbil + 273.15))
plot((tmn_bil + 273.15) - (tmn_Tbil + 273.15))
plot((tas_bil + 273.15) - (tas_Tbil + 273.15))

prec_delta <- pr_bil / pr_Tbil
minmax(prec_delta)
names(prec_delta) <- month.abb
prec_delta
plot(prec_delta, fun = function() lines(land, col = "white"))
plot(mask(prec_delta, land), fun = function() lines(land, col = "black"))

c(pr_Tbil[[1]] * 1, pr_bil[[1]] * 1, pr_Tbil[[1]] * prec_delta[[1]]) * 86400 * 365
plot(c(pr_Tbil[[1]], pr_bil[[1]], pr_Tbil[[1]] * prec_delta[[1]]) * 86400 * 365,
    fun = function() lines(land),
    range = c(0, 12000),
    main = c("Trace fine", "Chelsa coarse", "Trace corrected"),
    col = hcl.colors(100, "Batlow")
)

tsmx_delta <- (tmx_bil + 273.15) - (tmx_Tbil + 273.15)
tsmx_delta
plot(tsmx_delta, fun = function() lines(land, col = "white"))
plot(c(tmx_Tbil[[1]], tmx_bil[[1]], tmx_Tbil[[1]] + tsmx_delta[[1]]) + 273.15,
    fun = function() lines(land),
    range = c(0, 40) + 273.15,
    main = c("Trace fine", "Chelsa coarse", "Trace corrected"),
    col = hcl.colors(100, "Batlow")
)

tsmn_delta <- (tmn_bil + 273.15) - (tmn_Tbil + 273.15)
tsmn_delta
plot(tsmn_delta, range = c(-25, 15))
plot(c(tmn_Tbil[[1]], tmn_bil[[1]], tmn_Tbil[[1]] + tsmn_delta[[1]]) + 273.15,
    fun = function() lines(land),
    range = c(0, 30) + 273.15,
    main = c("Trace fine", "Chelsa coarse", "Trace corrected"),
    col = hcl.colors(100, "Batlow")
)

tas_delta <- (tas_bil + 273.15) - (tas_Tbil + 273.15)
tas_delta
plot(tas_delta, range = c(-22, 7))
plot(c(tas_Tbil[[1]], tas_bil[[1]], tas_Tbil[[1]] + tas_delta[[1]]) + 273.15,
    fun = function() lines(land),
    range = c(5, 35) + 273.15,
    main = c("Trace fine", "Chelsa coarse", "Trace corrected"),
    col = hcl.colors(100, "Batlow")
)

time(prec_delta, tstep = "months") <- seq(as.Date("1984-01-16"), by = "months", l = 12)
names(prec_delta) <- month.abb
prec_delta
writeCDF(prec_delta,
    "02_data/02_processed/TRACE/TRACE_pr_delta.nc",
    varname = "precip_delta", longname = "correction factor for precipitation", # compression = 6L,
    unit = "", zname = "year", prec = "float",
    overwrite = TRUE
)

time(tsmn_delta, tstep = "months") <- seq(as.Date("1984-01-16"), by = "months", l = 12)
names(tsmn_delta) <- month.abb
tsmn_delta
writeCDF(tsmn_delta,
    "02_data/02_processed/TRACE/TRACE_tsmn_delta.nc",
    varname = "tasmin_delta", longname = "correction factor for minimum temperature", # compression = 6L,
    unit = "", zname = "year", prec = "float",
    overwrite = TRUE
)

time(tsmx_delta, tstep = "months") <- seq(as.Date("1984-01-16"), by = "months", l = 12)
names(tsmx_delta) <- month.abb
tsmx_delta
writeCDF(tsmx_delta,
    "02_data/02_processed/TRACE/TRACE_tsmx_delta.nc",
    varname = "tasmax_delta", longname = "correction factor for maximum temperature", # compression = 6L,
    unit = "", zname = "year", prec = "float",
    overwrite = TRUE
)

time(tas_delta, tstep = "months") <- seq(as.Date("1984-01-16"), by = "months", l = 12)
names(tas_delta) <- month.abb
tas_delta
writeCDF(tas_delta,
    "02_data/02_processed/TRACE/TRACE_tas_delta.nc",
    varname = "tas_delta", longname = "correction factor for average temperature", # compression = 6L,
    unit = "", zname = "month", prec = "float",
    overwrite = TRUE
)

tas_delta <- rast("02_data/02_processed/TRACE/TRACE_tas_delta.nc") * 1
tsmx_delta <- rast("02_data/02_processed/TRACE/TRACE_tsmx_delta.nc") * 1
tsmn_delta <- rast("02_data/02_processed/TRACE/TRACE_tsmn_delta.nc") * 1
prec_delta <- rast("02_data/02_processed/TRACE/TRACE_pr_delta.nc") * 1
time(tas_delta, tstep = "months") <- time(tsmx_delta, tstep = "months") <-
    time(tsmn_delta, tstep = "months") <- time(prec_delta, tstep = "months") <- seq(as.Date("1984-01-16"), by = "months", l = 12)

#### BIAS CORRECTION ####
# downscale and then bias correct pr, tasmin, tasmax, and tas
in_files <- list.files("02_data/02_processed",
    pattern = "\\.nc$",
    full.names = TRUE, recursive = FALSE
)
in_files <- in_files[grepl("tas|pr.nc", in_files)]
in_files

delta_sds <- sds(prec_delta, tsmn_delta, tsmx_delta, tas_delta)
names(delta_sds) <- c("pr", "tasmin", "tasmax", "tas")
delta_sds

source("01_code/00_functions/downscale_bias_correct.r")

downscaled_and_corrected <- lapply(in_files, downscale_bias_correct,
    delta_sds = delta_sds,
    convert_to_annual = FALSE,
    cores = 72L
)
names(downscaled_and_corrected) <- c("pr", "tas", "tasmax", "tasmin")
downscaled_and_corrected

sapply(seq_along(downscaled_and_corrected), function(i) {
    r <- downscaled_and_corrected[[i]]
    writeCDF(
        x = r,
        filename = sprintf("02_data/02_processed/TRACE/%s_downscaled_and_bc.nc", varnames(r[[1]])),
        varname = varnames(r[[1]]), longname = varnames(r[[1]]),
        unit = units(r[[1]]), zname = "time", prec = "float",
        overwrite = TRUE
    )
})

par(mfrow = c(2, 2))
plot(app(downscaled_and_corrected[["pr"]][[961:1080]], mean), main = "pr", fun = function() lines(land))
plot(app(downscaled_and_corrected[["tas"]][[961:1080]], mean) - 273.15, main = "tas", fun = function() lines(land))
plot(app(downscaled_and_corrected[["tasmax"]][[961:1080]], mean) - 273.15, main = "tasmax", fun = function() lines(land))
plot(app(downscaled_and_corrected[["tasmin"]][[961:1080]], mean) - 273.15, main = "tasmin", fun = function() lines(land))

#### SUBSET FOR TEST RUN ####
library(terra)

## SUBSET IS CLIMATOLOGICAL MONTHLY AVERAGES OVER 1900-1990 CENTERED ON 1950
in_files <- list.files("02_data/02_processed",
    pattern = ".nc$", full.names = TRUE,
    recursive = FALSE
)
in_files <- in_files[!grepl("tas|pr.nc", in_files)]

# downscaled temp and precip
in_files_bc <- list.files("02_data/02_processed/TRACE",
    pattern = ".nc$", full.names = TRUE,
    recursive = FALSE
)
in_files_bc <- in_files_bc[grepl("downscaled", in_files_bc)]
in_files_bc

in_files <- c(in_files_bc, in_files)
in_files

template <- rast(
    res = 0.5, crs = "EPSG:4326",
    extent = ext(105, 160, -50, 7)
)

sapply(paste0("02_data/03_CHELSA_paleo/", c("orog", "static", "clim", "out/tas", "out/tasmax", "out/tasmin", "out/pr")), dir.create, recursive = TRUE)

# move files to chelsa_paleo folder.
# Manually move into subfolders.
pbsapply(in_files, function(f) {
    in_r <- rast(f)
    if (varnames(in_r) %in% c("pr", "tas", "tasmin", "tasmax")) {
        f_out <- sapply(strsplit(f, "/"), "[", 1:4)
        v_out <- strsplit(f_out[4], "_")[[1]][1]
        f_out <- file.path(paste(c(f_out[1], "03_CHELSA_paleo"), collapse = "/"), paste0(v_out, ".nc"))
        writeCDF(in_r,
            filename = f_out,
            zname = "time", prec = "float",
            varname = terra::varnames(in_r),
            longname = terra::longnames(in_r),
            unit = terra::units(in_r),
            overwrite = TRUE
        )
        return(f_out)
    } else {
        f_out <- sapply(strsplit(f, "/"), "[", 1:4)
        f_out <- file.path(paste(c(f_out[1], "03_CHELSA_paleo"), collapse = "/"), basename(f))
        writeCDF(in_r,
            filename = f_out,
            zname = "time", prec = "float",
            varname = terra::varnames(in_r),
            longname = terra::longnames(in_r),
            unit = terra::units(in_r),
            overwrite = TRUE
        )
        return("f_out")
    }
})
