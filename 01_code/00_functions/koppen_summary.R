koppen_summary <- function(grids) {
  require(data.table); require(pbapply); require(exactextractr)
  ext <- pblapply(seq_along(grids), function(i) {
    vars <- grids[[i]]
    var_extract <- lapply(seq_along(vars), function(j,...) {
      vj <- vars[[j]]
      zt <- time(vj)
      zt <- data.table(timestep = 1:length(zt), yd = zt)
      vje <- setDT(exact_extract(vj, st_as_sf(koppen),
                                 fun = c("mean", "stdev", "coefficient_of_variation"),
                                 append_cols = "Koppen"))
      vje <- melt(vje, id.vars = "Koppen")
      vje[, c("sumstat", "rest") := tstrsplit(as.character(vje[["variable"]]),
                                              "\\.(?=[^\\.]+_[0-9]+$)", perl = TRUE)]
      vje[, c("climvar", "timestep") := tstrsplit(rest, "_")]
      vje[, timestep := as.integer(timestep)][, rest := NULL][, variable := NULL]
      vje <- merge(vje, zt, by = "timestep", all.x = TRUE, all.y = FALSE)
      fout <- file.path(sprintf("scratch/var_extract_%s_%s.RDS", i, j))
      saveRDS(vje, fout)
      return(vje)
    })
    var_extract <- rbindlist(var_extract)
    return(var_extract)
  })
  return(ext)
}
