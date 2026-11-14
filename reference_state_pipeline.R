# ---------------------------------------------------------------------------
# Global availability of contemporary protected reference states
# R / terra pipeline -- NOT YET RUN.
#
# Every step caches its output and is skipped if that output already exists,
# so the script can be stopped and restarted. Downloads are large (on the order
# of 150-250 GB in transit, mostly WorldCover and gHM); WorldCover tiles are
# deleted once aggregated unless cfg$keep_wc_tiles is TRUE.
#
# Packages: terra, data.table, arrow, dplyr, RANN, nabor
#   renv::install(c("terra", "data.table", "arrow", "dplyr", "RANN", "nabor", "curl"))
#   renv::snapshot()
#
# WDPA / WD-OECM terms of use: non-commercial use only, with citation. Read
# https://www.protectedplanet.net/en/legal before downloading.
# ---------------------------------------------------------------------------
pkgs <- c("terra", "data.table", "arrow", "dplyr", "RANN", "nabor", "curl")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) {
  stop("Install first: ", paste(miss, collapse = ", "))
}
library(terra)
library(data.table)
library(arrow)
library(dplyr)
library(RANN)

cfg <- list(
  data_dir = "data",
  out_dir = "out",
  ref_year = 2021, # land cover year; PAs filtered to STATUS_YR <= ref_year
  ghm_year = 2020, # nearest year in the gHM v3 consistent time series
  ymin = -56, # EarthEnv topography covers 56 S to 84 N
  ymax = 84,
  nat_min = 0.90, # X: minimum natural fraction for a natural pixel
  ghm_max = 0.10, # hard filter on the reference pool (Inf = no filter)
  pa_cover = 0.50, # union PA area fraction for a cell to count as protected
  pa_set = "strict", # "strict" = IUCN Ia-IV + OECMs; "all" = every category + OECMs
  iucn_strict = c("Ia", "Ib", "II", "III", "IV"),
  pc_var = 0.95, # variance retained in the PCA
  cal_q = 0.90, # caliper = this quantile of pool-to-pool (other PA) NN distances
  cal_n = 2e5, # pool pixels sampled for calibration
  cal_k = 200, # neighbours searched before the per-PA fallback
  k = 100, # max candidates returned per focal pixel
  chunk = 5e4, # focal pixels per nn2 call
  pa_tile = 10, # degrees; tile size for PA rasterization
  near_km = c(1, 5, 10, 25, 50), # bins for distance from focal pixel to nearest PA
  keep_wc_tiles = FALSE
)
for (d in c(cfg$data_dir, cfg$out_dir, file.path(cfg$out_dir, "tmp"))) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
terraOptions(
  tempdir = file.path(cfg$out_dir, "tmp"),
  memfrac = 0.6,
  progress = 1
)
options(timeout = 6 * 3600)

dpath <- function(...) file.path(cfg$data_dir, ...)
opath <- function(...) file.path(cfg$out_dir, ...)

## Helpers -------------------------------------------------------------------
# Download to dest unless it already exists; writes to .part first so an
# interrupted download is not mistaken for a complete one.
# Data go to dest.part with a progress bar. If interrupted, the next call asks the
# server for the remaining bytes (HTTP Range) and appends, rather than starting over.
# The .part file is renamed to dest only once complete. Retries a few times.
dl <- function(url, dest, required = TRUE, tries = 5) {
  if (file.exists(dest) && file.size(dest) > 0) {
    return(invisible(TRUE))
  }
  dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(dest, ".part")
  for (i in seq_len(tries)) {
    res <- tryCatch(
      curl::multi_download(url, tmp, resume = TRUE, progress = interactive()),
      error = function(e) data.frame(success = FALSE, status_code = NA, error = conditionMessage(e))
    )
    # 200 = full download, 206 = resumed, 416 = .part was already complete
    if (isTRUE(res$success) && res$status_code %in% c(200, 206, 416)) {
      return(invisible(file.rename(tmp, dest)))
    }
    if (res$status_code %in% c(403, 404)) break   # not there; retrying won't help
    message(sprintf("Download attempt %d/%d failed (%s): %s", i, tries,
                    res$status_code, basename(dest)))
    Sys.sleep(5 * i)
  }
  if (res$status_code %in% c(403, 404)) unlink(tmp)   # keep partial data otherwise
  if (required) {
    stop("Download failed: ", url)
  }
  invisible(FALSE)
}

unzip_once <- function(zip, exdir) {
  if (!dir.exists(exdir) || !length(list.files(exdir, recursive = TRUE))) {
    unzip(zip, exdir = exdir)
  }
  exdir
}

# Evaluate expr (a SpatRaster) and write it to path, unless path already exists.
cached_rast <- function(path, expr, datatype = "FLT4S") {
  if (!file.exists(path)) {
    writeRaster(expr, path, datatype = datatype, overwrite = TRUE)
  }
  rast(path)
}

# First of several candidate attribute names present in a table (field names
# differ between WDPA schema versions).
fld <- function(x, cands, required = TRUE) {
  hit <- intersect(cands, names(x))
  if (length(hit)) {
    return(hit[1])
  }
  if (required) {
    stop(
      "None of [",
      paste(cands, collapse = ", "),
      "] found. Check the schema."
    )
  }
  NA_character_
}

chord_km <- function(chord) 2 * 6371 * asin(pmin(chord / 2, 1))

## 0. Downloads ------------------------------------------------------------------
# Climate: WorldClim 2.1, 30 arc-sec bioclim (1970-2000 normals), ~10 GB zip
dl(
  "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_bio.zip",
  dpath("worldclim", "wc2.1_30s_bio.zip")
)
unzip_once(dpath("worldclim", "wc2.1_30s_bio.zip"), dpath("worldclim", "bio"))

# Topography: EarthEnv (Amatulli et al. 2018), 1 km, median aggregation of GMTED2010 median
topo_vars <- c("elevation", "slope", "roughness", "tpi")
for (v in topo_vars) {
  dl(
    sprintf("https://data.earthenv.org/topography/%s_1KMmd_GMTEDmd.tif", v),
    dpath("earthenv", sprintf("%s_1KMmd_GMTEDmd.tif", v))
  )
}

# Human modification: gHM v3 (Theobald et al. 2025, Sci Data), 300 m, all threats
# combined (AA), consistent time series. ~9.5 GB. Zenodo record 14449495.
ghm_name <- sprintf("HMv20240801_%sc_AA.tif", cfg$ghm_year)
dl(
  sprintf("https://zenodo.org/api/records/14449495/files/%s/content", ghm_name),
  dpath("ghm", ghm_name)
)

# Ecoregions: WWF / RESOLVE Ecoregions 2017 (Dinerstein et al. 2017)
dl(
  "https://storage.googleapis.com/teow2016/Ecoregions2017.zip",
  dpath("ecoregions", "Ecoregions2017.zip")
)
unzip_once(
  dpath("ecoregions", "Ecoregions2017.zip"),
  dpath("ecoregions", "shp")
)

# Countries: Natural Earth 1:10m admin 0
dl(
  "https://naciscdn.org/naturalearth/10m/cultural/ne_10m_admin_0_countries.zip",
  dpath("countries", "ne_10m_admin_0_countries.zip")
)
unzip_once(
  dpath("countries", "ne_10m_admin_0_countries.zip"),
  dpath("countries", "shp")
)

# Protected areas: current WDPA + WD-OECM release. The link redirects to the
# month's file (e.g. WDPA_WDOECM_Nov2026_Public.zip). Record which release was used.
wdpa_zip <- dpath("wdpa", "WDPA_WDOECM_current_Public.zip")
if (!file.exists(wdpa_zip)) {
  dl("https://wcmc.io/wdpa_current_release", wdpa_zip)
  writeLines(
    c(
      paste("downloaded:", Sys.Date()),
      "source: https://wcmc.io/wdpa_current_release"
    ),
    dpath("wdpa", "RELEASE.txt")
  )
}
unzip_once(wdpa_zip, dpath("wdpa", "unz"))
# Shapefile releases nest further zips; expand them too.
for (z in list.files(
  dpath("wdpa", "unz"),
  "\\.zip$",
  full.names = TRUE,
  recursive = TRUE
)) {
  unzip_once(z, sub("\\.zip$", "", z))
}

# Land cover: ESA WorldCover 2021 v200, 3 x 3 degree tiles at 10 m. Fetched in
# step 3 one tile at a time, aggregated, then deleted.
dl(
  "https://esa-worldcover.s3.eu-central-1.amazonaws.com/esa_worldcover_grid.geojson",
  dpath("worldcover", "esa_worldcover_grid.geojson")
)

## 1. Template grid and layer alignment -------------------------------------------
template <- rast(
  xmin = -180,
  xmax = 180,
  ymin = cfg$ymin,
  ymax = cfg$ymax,
  resolution = 1 / 120,
  crs = "EPSG:4326"
) # 30 arc-sec

align <- function(r, method = "bilinear") {
  r <- crop(r, template, snap = "near")
  if (isTRUE(compareGeom(r, template, stopOnError = FALSE))) {
    return(r)
  }
  if (!same.crs(r, template)) {
    return(project(r, template, method = method))
  }
  resample(r, template, method = method)
}

bio_files <- list.files(
  dpath("worldclim", "bio"),
  "bio_[0-9]+\\.tif$",
  full.names = TRUE,
  recursive = TRUE
)
bio_files <- bio_files[order(as.integer(sub(
  ".*bio_([0-9]+)\\.tif$",
  "\\1",
  bio_files
)))]
stopifnot(length(bio_files) == 19)
clim <- cached_rast(opath("clim.tif"), align(rast(bio_files)))
names(clim) <- paste0("bio", 1:19)
topo <- cached_rast(
  opath("topo.tif"),
  align(rast(dpath("earthenv", sprintf("%s_1KMmd_GMTEDmd.tif", topo_vars))))
)
names(topo) <- topo_vars
env <- c(clim, topo)

# gHM: 300 m to 1 km by area average
ghm <- cached_rast(
  opath("ghm.tif"),
  align(rast(dpath("ghm", ghm_name)), "average")
)
local({
  s <- spatSample(ghm, 1e5, na.rm = TRUE)[[1]]
  if (min(s) < 0 || max(s) > 1) {
    stop("gHM values outside 0-1; check scaling of ", ghm_name)
  }
})

# Ecoregions (exact-match stratum)
eco_v <- vect(list.files(
  dpath("ecoregions", "shp"),
  "\\.shp$",
  full.names = TRUE,
  recursive = TRUE
)[1])
eco_lookup <- as.data.table(as.data.frame(eco_v))[, .(
  eco = get(fld(eco_v, "ECO_ID")),
  ECO_NAME = get(fld(eco_v, "ECO_NAME")),
  BIOME_NUM = get(fld(eco_v, "BIOME_NUM")),
  BIOME_NAME = get(fld(eco_v, "BIOME_NAME")),
  REALM = get(fld(eco_v, "REALM"))
)]
eco_lookup <- unique(eco_lookup, by = "eco")
saveRDS(eco_lookup, opath("eco_lookup.rds"))
eco <- cached_rast(
  opath("eco.tif"),
  rasterize(eco_v, template, field = fld(eco_v, "ECO_ID")),
  datatype = "INT4S"
)

# Countries. ADM0_A3 rather than ISO_A3, which is "-99" for France, Norway and others.
cty_v <- vect(list.files(
  dpath("countries", "shp"),
  "\\.shp$",
  full.names = TRUE,
  recursive = TRUE
)[1])
cty_v$cid <- seq_len(nrow(cty_v))
cty_lookup <- as.data.table(as.data.frame(cty_v))[, .(
  cid,
  ADM0_A3 = get(fld(cty_v, "ADM0_A3")),
  NAME = get(fld(cty_v, c("NAME", "ADMIN")))
)]
saveRDS(cty_lookup, opath("cty_lookup.rds"))
cty <- cached_rast(
  opath("cty.tif"),
  rasterize(cty_v, template, field = "cid"),
  datatype = "INT4S"
)

## 2. Protected areas ------------------------------------------------------------
# Polygons only (point records dropped). Kept: designated / inscribed / established,
# not marine, not UNESCO MAB, STATUS_YR <= ref_year (0 = year unknown, kept).
# "strict" set: OECMs plus PAs in IUCN Ia-IV. "all" set: every category plus OECMs.
pa_files <- c(
  opath("pa_frac_s.tif"),
  opath("pa_id_s.tif"),
  opath("pa_frac_a.tif"),
  opath("pa_id_a.tif")
)
if (!all(file.exists(pa_files))) {
  gdb <- list.files(
    dpath("wdpa", "unz"),
    "\\.gdb$",
    full.names = TRUE,
    recursive = TRUE,
    include.dirs = TRUE
  )
  if (length(gdb)) {
    lyr <- grep(
      "poly",
      vector_layers(gdb[1]),
      ignore.case = TRUE,
      value = TRUE
    )[1]
    pa_v <- vect(gdb[1], layer = lyr)
  } else {
    shp <- list.files(
      dpath("wdpa", "unz"),
      "polygons\\.shp$",
      full.names = TRUE,
      recursive = TRUE
    )
    pa_v <- do.call(rbind, lapply(shp, vect))
  }
  a <- as.data.frame(pa_v)
  status <- a[[fld(a, "STATUS")]]
  yr <- suppressWarnings(as.integer(a[[fld(a, "STATUS_YR")]]))
  desig <- a[[fld(a, c("DESIG_ENG", "DESIG"))]]
  iucn <- a[[fld(a, "IUCN_CAT")]]
  # Marine flag: MARINE ("0" terrestrial, "1" coastal, "2" marine) or REALM in newer schemas
  marine <- if (!is.na(fld(a, "MARINE", FALSE))) {
    as.character(a$MARINE) %in% "2"
  } else {
    grepl("^marine$", a[[fld(a, "REALM")]], ignore.case = TRUE)
  }
  # OECM flag: PA_DEF ("1" PA, "0" OECM) or SITE_TYPE in newer schemas
  oecm <- if (!is.na(fld(a, "PA_DEF", FALSE))) {
    as.character(a$PA_DEF) %in% "0"
  } else {
    grepl("OECM", a[[fld(a, "SITE_TYPE")]], ignore.case = TRUE)
  }
  keep <- status %in%
    c("Designated", "Inscribed", "Established") &
    !grepl("Biosphere", desig, ignore.case = TRUE) &
    !marine &
    !is.na(yr) &
    yr <= cfg$ref_year
  pa_v <- pa_v[keep, ]
  a <- a[keep, ]
  pa_v$pid <- seq_len(nrow(pa_v))
  pa_v$strict <- oecm[keep] | iucn[keep] %in% cfg$iucn_strict
  pa_v$area_km2 <- a[[fld(a, c("GIS_AREA", "REP_AREA"))]]
  pa_lookup <- data.table(
    pid = pa_v$pid,
    WDPA_PID = a[[fld(a, c("WDPA_PID", "SITE_PID"))]],
    NAME = a[[fld(a, c("NAME", "NAME_ENG"))]],
    ISO3 = a[[fld(a, c("ISO3", "PARENT_ISO3", "PARENT_ISO"))]],
    IUCN_CAT = iucn[keep],
    OECM = oecm[keep],
    STATUS_YR = yr[keep],
    strict = pa_v$strict,
    area_km2 = pa_v$area_km2
  )
  saveRDS(pa_lookup, opath("pa_lookup.rds"))
  pa_v <- makeValid(pa_v[, c("pid", "strict", "area_km2")])

  # Tile by tile:
  #   pa_frac = area fraction of the cell covered by the union of PAs (dissolved first,
  #             so overlapping polygons are not double counted)
  #   pa_id   = the PA covering the largest share of the cell; ties go to the smaller PA
  tdir <- opath("pa_tiles")
  dir.create(tdir, showWarnings = FALSE)
  xs <- seq(-180, 180 - cfg$pa_tile, by = cfg$pa_tile)
  ys <- seq(
    floor(cfg$ymin / cfg$pa_tile) * cfg$pa_tile,
    cfg$ymax - 1e-9,
    by = cfg$pa_tile
  )
  for (x0 in xs) {
    for (y0 in ys) {
      dst <- file.path(tdir, sprintf("pa_%04d_%04d.tif", x0, y0))
      if (file.exists(dst)) {
        next
      }
      te <- intersect(
        ext(x0, x0 + cfg$pa_tile, y0, y0 + cfg$pa_tile),
        ext(template)
      )
      if (is.null(te)) {
        next
      }
      pv <- crop(pa_v, te)
      if (nrow(pv) == 0L) {
        next
      }
      tt <- init(crop(template, te), 1)
      ex <- as.data.table(extract(tt, pv, exact = TRUE, cells = TRUE))
      ex[, `:=`(
        pid = pv$pid[ID],
        strict = pv$strict[ID],
        area_km2 = pv$area_km2[ID]
      )]
      layer_pair <- function(sel) {
        fr <- rast(tt, vals = 0)
        id <- rast(tt, vals = NA_real_)
        if (any(sel)) {
          u <- tryCatch(aggregate(pv[sel, ]), error = function(e) NULL)
          if (!is.null(u)) {
            fr <- rasterize(u, tt, cover = TRUE, background = 0)
          } else {
            # dissolve failed: sum of fractions capped at 1 (overstates overlaps)
            warning(
              "PA dissolve failed in tile ",
              basename(dst),
              "; using capped sum"
            )
            s <- ex[
              pid %in% pv$pid[sel],
              .(f = min(1, sum(fraction))),
              by = cell
            ]
            fr[s$cell] <- s$f
          }
          b <- ex[pid %in% pv$pid[sel]][order(cell, -fraction, area_km2)]
          b <- unique(b, by = "cell")
          id[b$cell] <- b$pid
        }
        c(fr, id)
      }
      s_pair <- layer_pair(pv$strict)
      a_pair <- layer_pair(rep(TRUE, nrow(pv)))
      out <- c(s_pair, a_pair)
      names(out) <- c("pa_frac_s", "pa_id_s", "pa_frac_a", "pa_id_a")
      writeRaster(out, dst, datatype = "FLT4S", overwrite = TRUE)
    }
  }
  pa_all <- extend(
    vrt(list.files(tdir, "\\.tif$", full.names = TRUE)),
    template
  )
  pa_all <- crop(pa_all, template)
  for (i in 1:4) {
    lyr <- pa_all[[i]]
    if (i %in% c(1, 3)) {
      lyr <- ifel(is.na(lyr), 0, lyr)
    } # cells outside every tile: unprotected
    writeRaster(lyr, pa_files[i], datatype = "FLT4S", overwrite = TRUE)
  }
}
pa_frac_s <- rast(pa_files[1])
pa_id_s <- rast(pa_files[2])
pa_frac_a <- rast(pa_files[3])
pa_id_a <- rast(pa_files[4])

## 3. WorldCover class fractions at 1 km -----------------------------------------------
# Each 3 x 3 degree tile is 36000 x 36000 cells at 1/12000 degree, so 100 x 100
# blocks land exactly on the 1/120 grid. Fractions are of valid (non-nodata) pixels.
wc_classes <- c(
  tree = 10,
  shrub = 20,
  grass = 30,
  crop = 40,
  built = 50,
  sparse = 60,
  ice = 70,
  water = 80,
  herb_wet = 90,
  mangrove = 95,
  moss = 100
)

wc_tile_fractions <- function(src, dst) {
  r <- rast(src)
  stopifnot(nrow(r) %% 100 == 0, ncol(r) %% 100 == 0)
  codes <- c(0, wc_classes) # 0 = no data
  nk <- length(codes)
  lut <- rep(1L, 256)
  lut[codes + 1] <- seq_len(nk) # value -> 1..nk; unknown -> no data
  nbc <- ncol(r) / 100
  nbr <- nrow(r) / 100
  colblk <- rep((seq_len(ncol(r)) - 1L) %/% 100L, 100L) # block column, row-major order
  out <- array(NA_real_, c(nbr, nbc, nk - 1))
  readStart(r)
  on.exit(readStop(r))
  for (b in seq_len(nbr)) {
    v <- readValues(r, row = (b - 1) * 100 + 1, nrows = 100)
    v[is.na(v)] <- 0
    cnt <- matrix(
      tabulate(colblk * nk + lut[v + 1], nbins = nbc * nk),
      nrow = nk
    )
    valid <- colSums(cnt[-1, , drop = FALSE])
    out[b, , ] <- t(cnt[-1, , drop = FALSE]) / ifelse(valid > 0, valid, NA)
  }
  fr <- rast(out, extent = ext(r), crs = crs(r))
  names(fr) <- names(wc_classes)
  writeRaster(fr, dst, datatype = "FLT4S", overwrite = TRUE)
}

frac_dir <- opath("wc_frac")
dir.create(frac_dir, showWarnings = FALSE)
wc_grid <- vect(dpath("worldcover", "esa_worldcover_grid.geojson"))
tiles <- wc_grid$ll_tile
tile_lat <- as.integer(substr(tiles, 2, 3)) *
  ifelse(substr(tiles, 1, 1) == "S", -1, 1)
tiles <- tiles[tile_lat + 3 > cfg$ymin & tile_lat < cfg$ymax]
for (t in tiles) {
  dst <- file.path(frac_dir, sprintf("wc_frac_%s.tif", t))
  if (file.exists(dst)) {
    next
  }
  fn <- sprintf("ESA_WorldCover_10m_2021_v200_%s_Map.tif", t)
  src <- dpath("worldcover", "tiles", fn)
  if (
    !dl(
      paste0(
        "https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/",
        fn
      ),
      src,
      required = FALSE
    )
  ) {
    warning("Missing WorldCover tile ", t)
    next
  }
  wc_tile_fractions(src, dst)
  if (!cfg$keep_wc_tiles) unlink(src)
}
wc <- crop(
  extend(vrt(list.files(frac_dir, "\\.tif$", full.names = TRUE)), template),
  template
)
names(wc) <- names(wc_classes)

## 4. Natural fraction and modal natural class ------------------------------------------
# Collapsed classes, in tie-break order: when two classes have equal shares,
# which.max keeps the first, so wetland > tree > shrub > grass > sparse.
nat_cls <- c("wet", "tree", "shrub", "grass", "sparse")
nat <- cached_rast(
  opath("nat.tif"),
  c(
    wc$herb_wet + wc$mangrove, # wetland (herbaceous wetland + mangroves)
    wc$tree,
    wc$shrub,
    wc$grass,
    wc$sparse + wc$moss
  )
) # sparse / bare (+ moss and lichen)
names(nat) <- nat_cls
nat_frac <- cached_rast(opath("nat_frac.tif"), sum(nat))
dom_idx <- cached_rast(opath("dom_idx.tif"), which.max(nat), datatype = "INT1U") # 1..5
excl <- cached_rast(
  opath("excl.tif"),
  (wc$water + wc$ice) > 0.5,
  datatype = "INT1U"
)

## 5. Environmental space: PCA and whitening ---------------------------------------------
# All 19 bioclim variables plus four topography variables. The PCA absorbs their
# correlation; components are then scaled to unit variance (whitened).
pca_file <- opath("pca.rds")
if (!file.exists(pca_file)) {
  set.seed(1)
  smp <- spatSample(env, 2e5, method = "regular", na.rm = TRUE, as.df = TRUE)
  pca <- prcomp(smp, center = TRUE, scale. = TRUE)
  npc <- which(cumsum(pca$sdev^2) / sum(pca$sdev^2) >= cfg$pc_var)[1]
  saveRDS(list(pca = pca, npc = npc), pca_file)
}
pca_obj <- readRDS(pca_file)
npc <- pca_obj$npc
pcs <- cached_rast(
  opath("pcs.tif"),
  predict(env, pca_obj$pca, index = seq_len(npc)) /
    pca_obj$pca$sdev[seq_len(npc)]
)
names(pcs) <- paste0("PC", seq_len(npc))
pcs_cols <- names(pcs)

## 6. Independent protected patches of the baseline pool -------------------------------------
# Connected components (8-neighbour) of the baseline reference pool. Patch counts
# per stratum use this baseline definition in every run.
patch <- cached_rast(
  opath("pool_patch.tif"),
  {
    pool_r <- (pa_frac_s >= cfg$pa_cover) &
      (nat_frac >= cfg$nat_min) &
      (ghm <= cfg$ghm_max) &
      (excl == 0)
    patches(pool_r, directions = 8, zeroAsNA = TRUE)
  },
  datatype = "INT4S"
)

## 7. Tabulate cells to disk, partitioned by ecoregion -----------------------------------------
# Raw layers are stored so thresholds and PA set can change at matching time
# without re-tabulating. ECO_ID 0 ("Rock and Ice") is dropped; check this in the
# ecoregion attribute table.
cells_dir <- opath("cells")
if (!file.exists(file.path(cells_dir, "DONE"))) {
  unlink(cells_dir, recursive = TRUE)
  stk <- c(
    eco,
    cty,
    dom_idx,
    nat_frac,
    excl,
    ghm,
    pa_frac_s,
    pa_id_s,
    pa_frac_a,
    pa_id_a,
    patch,
    cellSize(template, unit = "km"),
    pcs
  )
  names(stk) <- c(
    "eco",
    "cty",
    "dom_idx",
    "nat_frac",
    "excl",
    "ghm",
    "pa_frac_s",
    "pa_id_s",
    "pa_frac_a",
    "pa_id_a",
    "patch",
    "area",
    pcs_cols
  )
  rows_per <- 200
  starts <- seq(1, nrow(stk), by = rows_per)
  readStart(stk)
  for (i in seq_along(starts)) {
    nr <- min(rows_per, nrow(stk) - starts[i] + 1)
    v <- as.data.table(readValues(
      stk,
      row = starts[i],
      nrows = nr,
      dataframe = TRUE
    ))
    v[, cell := (starts[i] - 1) * ncol(stk) + seq_len(.N)]
    v <- v[!is.na(eco) & eco > 0 & !is.na(PC1) & !is.na(nat_frac) & excl == 0]
    if (nrow(v) == 0L) {
      next
    }
    v[is.na(pa_frac_s), pa_frac_s := 0]
    v[is.na(pa_frac_a), pa_frac_a := 0]
    v[, dom_idx := as.integer(dom_idx)]
    xy <- xyFromCell(stk, v$cell)
    lon <- xy[, 1] * pi / 180
    lat <- xy[, 2] * pi / 180
    v[, `:=`(X = cos(lat) * cos(lon), Y = cos(lat) * sin(lon), Z = sin(lat))] # unit sphere
    v[, excl := NULL]
    write_dataset(
      v,
      cells_dir,
      format = "parquet",
      partitioning = "eco",
      basename_template = paste0("chunk", i, "-{i}.parquet"),
      existing_data_behavior = "overwrite"
    )
  }
  readStop(stk)
  file.create(file.path(cells_dir, "DONE"))
}
ecos <- open_dataset(cells_dir) |>
  select(eco) |>
  distinct() |>
  collect() |>
  pull(eco)
read_eco <- function(e) {
  as.data.table(collect(filter(open_dataset(cells_dir), eco == !!e)))
}

## 8. Distance from each cell to the nearest protected cell (baseline PA set) -----------------
# For reporting how much "suitable reference" area sits next to a PA. One kd-tree
# over all protected cell centres (nabor keeps the tree between queries).
dist_dir <- opath("dist_pa")
if (!file.exists(file.path(dist_dir, "DONE"))) {
  dir.create(dist_dir, showWarnings = FALSE)
  pc <- open_dataset(cells_dir) |>
    filter(pa_frac_s >= !!cfg$pa_cover) |>
    select(X, Y, Z) |>
    collect() |>
    as.matrix()
  tree <- nabor::WKNNF(pc)
  rm(pc)
  gc()
  for (e in ecos) {
    d <- read_eco(e)[, .(cell, X, Y, Z)]
    q <- tree$query(as.matrix(d[, .(X, Y, Z)]), k = 1, eps = 0, radius = 0)
    write_parquet(
      data.table(cell = d$cell, dist_pa_km = chord_km(q$nn.dists[, 1])),
      file.path(dist_dir, sprintf("eco_%s.parquet", e))
    )
  }
  file.create(file.path(dist_dir, "DONE"))
}

## 9. Matching --------------------------------------------------------------------------
prep_flags <- function(d, cfg) {
  sfx <- if (cfg$pa_set == "strict") "s" else "a"
  d[, pa_frac := get(paste0("pa_frac_", sfx))]
  d[, pa_id := get(paste0("pa_id_", sfx))]
  d[, is_nat := nat_frac >= cfg$nat_min]
  d[, grp := fifelse(is_nat, 1L, 2L)] # 1 natural, 2 modified
  d[, dom := fifelse(is_nat, dom_idx, 0L)] # modal natural class, 0 if modified
  d[, prot := pa_frac >= cfg$pa_cover]
  ghm_ok <- if (is.finite(cfg$ghm_max)) {
    !is.na(d$ghm) & d$ghm <= cfg$ghm_max
  } else {
    rep(TRUE, nrow(d))
  }
  # hard filter; pool cells also need a PA id so that reliance can be attributed
  d[, pool := prot & is_nat & ghm_ok & !is.na(pa_id)]
  d[, focal := !prot]
  d
}

### 9a. Caliper calibration (one global rule) ---------------------------------------------
# For a sample of baseline pool pixels, find the environmental distance to the
# nearest pool pixel in the same ecoregion that lies in a *different* PA (pixels
# in the same PA are spatially autocorrelated and would make the caliper tiny).
# The caliper is the cal_q quantile of these distances, applied everywhere.
cal_file <- opath("caliper_calibration.rds")
if (!file.exists(cal_file)) {
  n_pool <- open_dataset(cells_dir) |>
    filter(
      pa_frac_s >= !!cfg$pa_cover,
      nat_frac >= !!cfg$nat_min,
      !is.na(ghm),
      ghm <= !!cfg$ghm_max,
      !is.na(pa_id_s)
    ) |>
    summarise(n = n()) |>
    collect() |>
    pull(n)
  frac <- min(1, cfg$cal_n / n_pool)
  set.seed(2)
  cal <- rbindlist(lapply(ecos, function(e) {
    p <- prep_flags(read_eco(e), cfg)[pool == TRUE]
    if (nrow(p) < 2L) {
      return(NULL)
    }
    ns <- rbinom(1, nrow(p), frac)
    if (ns == 0L) {
      return(NULL)
    }
    s <- sample(nrow(p), ns)
    P <- as.matrix(p[, ..pcs_cols])
    nn <- nn2(P, P[s, , drop = FALSE], k = min(cfg$cal_k, nrow(p)))
    own <- p$pa_id[s]
    other <- matrix(p$pa_id[nn$nn.idx], nrow = ns) != own
    has <- rowSums(other) > 0
    d <- rep(NA_real_, ns)
    d[has] <- nn$nn.dists[cbind(
      which(has),
      max.col(other, ties.method = "first")[has]
    )]
    # Pixels whose cal_k neighbours all sit in their own PA (large PAs): search
    # the pool outside that PA directly, so these larger distances are not lost.
    for (a in unique(own[!has])) {
      q <- which(!has & own == a)
      o <- which(p$pa_id != a)
      if (!length(o)) {
        next
      } # single-PA ecoregion: no other PA to compare
      d[q] <- nn2(
        P[o, , drop = FALSE],
        P[s[q], , drop = FALSE],
        k = 1
      )$nn.dists[, 1]
    }
    data.table(eco = e, pa = own, d = d)
  }))
  saveRDS(cal, cal_file)
}
cal <- readRDS(cal_file)
caliper_at <- function(q) unname(quantile(cal$d, q, na.rm = TRUE))
cfg$caliper <- caliper_at(cfg$cal_q)
message(sprintf(
  "Caliper = %.3f SD (q = %.2f; %d of %d sampled pixels had no other PA in their ecoregion)",
  cfg$caliper,
  cfg$cal_q,
  sum(is.na(cal$d)),
  nrow(cal)
))

### 9b. Matching within a stratum ----------------------------------------------------------
# Returns list(res = one row per focal pixel, pa = per-PA tallies for this stratum).
match_stratum <- function(f, p, cfg) {
  n <- nrow(f)
  base <- f[, .(cell, eco, cty, grp, dom, area)]
  if (nrow(p) == 0L) {
    res <- cbind(
      base,
      ref_cell = NA_real_,
      ref_cty = NA_real_,
      ref_pa = NA_real_,
      env_d = NA_real_,
      geo_km = NA_real_,
      n_cand = 0L,
      n_pa = 0L,
      fail = "no_pool_in_stratum"
    )
    return(list(res = res, pa = NULL))
  }
  P <- as.matrix(p[, ..pcs_cols])
  k <- min(cfg$k, nrow(p))
  complete_list <- k == nrow(p) # every pool pixel was considered
  out <- vector("list", ceiling(n / cfg$chunk))
  tallies <- vector("list", length(out))
  for (j in seq_along(out)) {
    ii <- ((j - 1) * cfg$chunk + 1):min(j * cfg$chunk, n)
    Q <- as.matrix(f[ii, ..pcs_cols])
    nn <- nn2(P, Q, k = k, searchtype = "radius", radius = cfg$caliper)
    idx <- nn$nn.idx
    ok <- idx > 0
    gd <- matrix(Inf, nrow(idx), ncol(idx)) # chord distance on the unit sphere
    cpa <- matrix(NA_real_, nrow(idx), ncol(idx)) # PA id of each candidate
    for (m in seq_len(ncol(idx))) {
      r <- which(ok[, m])
      q <- idx[r, m]
      gd[r, m] <- sqrt(
        (f$X[ii][r] - p$X[q])^2 +
          (f$Y[ii][r] - p$Y[q])^2 +
          (f$Z[ii][r] - p$Z[q])^2
      )
      cpa[r, m] <- p$pa_id[q]
    }
    fa <- f$area[ii]
    n_cand <- rowSums(ok)
    has <- n_cand > 0

    # Distinct PAs among the candidates, for shared and sole-supplier importance
    ok_v <- as.vector(ok)
    u <- unique(data.table(
      row = rep(seq_len(nrow(idx)), ncol(idx))[ok_v],
      pa = as.vector(cpa)[ok_v]
    ))
    u[, npa := .N, by = row]
    n_pa <- integer(nrow(idx))
    n_pa[u$row] <- u$npa
    u[, w := fa[row] / npa] # focal area split equally across its PAs
    sole_row <- n_pa == 1L & (n_cand < k | complete_list) # only if the list was not capped
    shared <- u[, .(shared_area = sum(w)), by = pa]
    sole <- u[npa == 1L & sole_row[row], .(sole_area = sum(fa[row])), by = pa]
    tallies[[j]] <- merge(shared, sole, by = "pa", all.x = TRUE)

    best <- max.col(-gd, ties.method = "first") # closest in geographic space
    sel <- cbind(seq_along(ii), best)
    pick <- pmax(idx[sel], 1L)
    out[[j]] <- data.table(
      ref_cell = ifelse(has, p$cell[pick], NA_real_),
      ref_cty = ifelse(has, p$cty[pick], NA_real_),
      ref_pa = ifelse(has, p$pa_id[pick], NA_real_),
      env_d = ifelse(has, nn$nn.dists[sel], NA_real_),
      geo_km = ifelse(has, chord_km(gd[sel]), NA_real_),
      n_cand = n_cand,
      n_pa = n_pa,
      fail = ifelse(has, NA_character_, "none_in_caliper")
    )
  }
  tally <- rbindlist(tallies)
  if (nrow(tally)) {
    tally[is.na(sole_area), sole_area := 0]
    tally <- tally[, lapply(.SD, sum), by = pa]
  }
  list(res = cbind(base, rbindlist(out)), pa = tally)
}

### 9c. A matching run over all ecoregions -------------------------------------------------
# mode "lfl": natural focal pixels, strata = ecoregion x modal natural class
# mode "env": all focal pixels,     strata = ecoregion
# by_country: strata additionally split by country, so only domestic references count
run_matching <- function(cfg, mode, by_country, run_name) {
  odir <- opath("runs", run_name)
  if (file.exists(file.path(odir, "DONE"))) {
    return(invisible(odir))
  }
  for (s in c("focal", "pa", "pool", "strata")) {
    dir.create(file.path(odir, s), showWarnings = FALSE, recursive = TRUE)
  }
  keys <- c(if (mode == "lfl") "dom", if (by_country) "cty")
  for (e in ecos) {
    d <- prep_flags(read_eco(e), cfg)
    fl <- d[focal == TRUE & (mode == "env" | grp == 1L)]
    pl <- d[pool == TRUE]
    if (by_country) {
      fl <- fl[!is.na(cty)]
    }
    if (nrow(fl) == 0L) {
      next
    }
    if (length(keys)) {
      fs <- split(fl, by = keys)
      ps <- split(pl, by = keys)
    } else {
      fs <- list(all = fl)
      ps <- list(all = pl)
    }
    parts <- lapply(names(fs), function(nm) {
      match_stratum(fs[[nm]], if (is.null(ps[[nm]])) pl[0] else ps[[nm]], cfg)
    })
    write_parquet(
      rbindlist(lapply(parts, `[[`, "res")),
      file.path(odir, "focal", sprintf("eco_%s.parquet", e))
    )
    pa_t <- rbindlist(lapply(parts, `[[`, "pa"))
    if (nrow(pa_t)) {
      write_parquet(
        pa_t[, lapply(.SD, sum), by = pa],
        file.path(odir, "pa", sprintf("eco_%s.parquet", e))
      )
    }
    if (nrow(pl)) {
      write_parquet(
        pl[, .(pool_area = sum(area)), by = .(pa = pa_id)],
        file.path(odir, "pool", sprintf("eco_%s.parquet", e))
      )
    }
    # Per-stratum counts, so a stratum with a handful of reference cells is not
    # read like one with thousands. Patches are from the baseline pool (step 6).
    st <- rbindlist(lapply(names(fs), function(nm) {
      p <- ps[[nm]]
      data.table(
        eco = e,
        stratum = nm,
        n_focal = nrow(fs[[nm]]),
        focal_area = sum(fs[[nm]]$area),
        n_pool = if (is.null(p)) 0L else nrow(p),
        pool_area = if (is.null(p)) 0 else sum(p$area),
        n_pa = if (is.null(p)) 0L else uniqueN(p$pa_id),
        n_patch = if (is.null(p)) 0L else uniqueN(na.omit(p$patch))
      )
    }))
    write_parquet(st, file.path(odir, "strata", sprintf("eco_%s.parquet", e)))
  }
  file.create(file.path(odir, "DONE"))
  invisible(odir)
}

## 10. Summaries --------------------------------------------------------------------------
# main_run and its country-restricted twin are read eco by eco and joined on cell:
#   matched        a reference exists within the caliper anywhere in the ecoregion
#   domestic       a reference exists within the caliper in the pixel's own country
#   only_foreign   matched, but not domestically (strict dependence on other countries)
#   chosen_foreign the chosen (nearest) reference is in another country (overstates dependence)
summarise_run <- function(main_run, cty_run, cfg) {
  sdir <- opath("summaries", main_run)
  if (file.exists(file.path(sdir, "DONE"))) {
    return(invisible(sdir))
  }
  dir.create(sdir, showWarnings = FALSE, recursive = TRUE)
  rd <- function(run, sub, e) {
    f <- opath("runs", run, sub, sprintf("eco_%s.parquet", e))
    if (file.exists(f)) as.data.table(read_parquet(f)) else NULL
  }
  dist_bins <- c(-Inf, cfg$near_km, Inf)
  acc <- setNames(
    rep(list(list()), 8),
    c("cty", "eco", "fail", "fdist", "envd", "near", "reuse", "chosen")
  )
  for (e in ecos) {
    x <- rd(main_run, "focal", e)
    if (is.null(x)) {
      next
    }
    y <- rd(cty_run, "focal", e)
    x[, dom_ok := FALSE]
    if (!is.null(y)) {
      x[y[!is.na(ref_cell)], dom_ok := TRUE, on = "cell"]
    }
    dp <- as.data.table(read_parquet(file.path(
      dist_dir,
      sprintf("eco_%s.parquet", e)
    )))
    x[dp, dist_pa_km := i.dist_pa_km, on = "cell"]
    x[, `:=`(
      a_match = fifelse(is.na(ref_cell), 0, area),
      a_cf = fifelse(!is.na(ref_cty) & !is.na(cty) & ref_cty != cty, area, 0),
      a_dom = fifelse(dom_ok, area, 0),
      a_of = fifelse(!is.na(ref_cell) & !is.na(cty) & !dom_ok, area, 0)
    )]
    acc$cty[[as.character(e)]] <- x[
      !is.na(cty),
      .(
        area = sum(area),
        matched = sum(a_match),
        domestic = sum(a_dom),
        chosen_foreign = sum(a_cf),
        only_foreign = sum(a_of)
      ),
      by = .(cty, grp)
    ]
    acc$eco[[as.character(e)]] <- x[,
      .(area = sum(area), matched = sum(a_match)),
      by = .(eco, grp)
    ]
    acc$fail[[as.character(e)]] <- x[, .(area = sum(area)), by = .(grp, fail)]
    acc$fdist[[as.character(e)]] <- x[
      a_cf > 0,
      .(area = sum(area)),
      by = .(km_bin = floor(geo_km / 50) * 50)
    ]
    acc$envd[[as.character(e)]] <- x[
      !is.na(env_d),
      .(area = sum(area)),
      by = .(grp, rel_bin = floor(env_d / cfg$caliper * 10) / 10)
    ]
    acc$near[[as.character(e)]] <- x[,
      .(area = sum(area), matched = sum(a_match)),
      by = .(grp, dist_pa = as.character(cut(dist_pa_km, dist_bins)))
    ]
    acc$reuse[[as.character(e)]] <- x[
      !is.na(ref_cell),
      .(n_focal = .N, n_ref = uniqueN(ref_cell))
    ]
    acc$chosen[[as.character(e)]] <- x[
      !is.na(ref_pa),
      .(chosen_area = sum(area), chosen_foreign_area = sum(a_cf)),
      by = .(pa = ref_pa)
    ]
  }
  sumby <- function(l, by) rbindlist(l)[, lapply(.SD, sum), by = by]

  by_cty <- sumby(acc$cty, c("cty", "grp"))
  by_cty[, `:=`(
    prop_matched = matched / area,
    prop_domestic = domestic / area,
    prop_chosen_foreign = chosen_foreign / area,
    prop_only_foreign = only_foreign / area
  )]
  by_cty <- merge(
    by_cty,
    readRDS(opath("cty_lookup.rds")),
    by.x = "cty",
    by.y = "cid",
    all.x = TRUE
  )

  by_eco <- merge(
    sumby(acc$eco, c("eco", "grp")),
    readRDS(opath("eco_lookup.rds")),
    by = "eco",
    all.x = TRUE
  )
  by_eco[, prop_matched := matched / area]
  by_biome <- by_eco[,
    .(area = sum(area), matched = sum(matched)),
    by = .(BIOME_NUM, BIOME_NAME, grp)
  ][, prop_matched := matched / area]

  strata <- rbindlist(lapply(
    list.files(opath("runs", main_run, "strata"), full.names = TRUE),
    read_parquet
  ))
  reuse <- rbindlist(acc$reuse)

  # Importance per protected area:
  #   chosen_*   focal area whose chosen (nearest) reference lies in the PA
  #   shared_*   focal area split equally across all PAs among its candidates
  #   sole_*     focal area for which the PA is the only supplier within the caliper
  #   *_per_pool each of the above divided by the PA's qualifying pool area
  rd_all <- function(sub) {
    rbindlist(lapply(
      list.files(opath("runs", main_run, sub), full.names = TRUE),
      read_parquet
    ))
  }
  pa_imp <- Reduce(
    function(a, b) merge(a, b, by = "pa", all = TRUE),
    list(
      sumby(acc$chosen, "pa"),
      sumby(list(rd_all("pa")), "pa"),
      sumby(list(rd_all("pool")), "pa")
    )
  )
  pa_imp[, `:=`(
    chosen_per_pool = chosen_area / pool_area,
    shared_per_pool = shared_area / pool_area,
    sole_per_pool = sole_area / pool_area
  )]
  pa_imp <- merge(
    pa_imp,
    readRDS(opath("pa_lookup.rds")),
    by.x = "pa",
    by.y = "pid",
    all.x = TRUE
  )

  res <- list(
    by_cty = by_cty,
    by_eco = by_eco,
    by_biome = by_biome,
    fail = sumby(acc$fail, c("grp", "fail")),
    foreign_dist = sumby(acc$fdist, "km_bin"),
    env_dist = sumby(acc$envd, c("grp", "rel_bin")),
    near_pa = sumby(acc$near, c("grp", "dist_pa")),
    strata = strata,
    pa_importance = pa_imp,
    reuse_ratio = sum(reuse$n_focal) / sum(reuse$n_ref),
    caliper = cfg$caliper
  )
  for (nm in names(res)) {
    if (is.data.frame(res[[nm]])) {
      fwrite(res[[nm]], file.path(sdir, paste0(nm, ".csv")))
    }
  }
  saveRDS(res, file.path(sdir, "summary.rds"))
  file.create(file.path(sdir, "DONE"))
  invisible(sdir)
}

## 11. Runs: baseline and one-at-a-time sensitivity --------------------------------------
# Each variant changes one setting from the baseline. The caliper value is held
# at its baseline calibration except in the cal_q variants, so "environmentally
# similar" means the same thing when the pool definition changes.
# Patch counts (step 6) and distance to PA (step 8) always use the baseline pool / PA set.
variants <- list(
  base = list(),
  ghm_0.4 = list(ghm_max = 0.4),
  ghm_none = list(ghm_max = Inf), # unfiltered; pixels without gHM are kept
  calq_0.75 = list(cal_q = 0.75),
  calq_0.95 = list(cal_q = 0.95),
  nat_0.80 = list(nat_min = 0.80),
  nat_0.95 = list(nat_min = 0.95),
  pa_all = list(pa_set = "all")
)
for (vn in names(variants)) {
  cfg_v <- modifyList(cfg, variants[[vn]])
  if (!is.null(variants[[vn]]$cal_q)) {
    cfg_v$caliper <- caliper_at(cfg_v$cal_q)
  }
  for (mode in c("lfl", "env")) {
    main <- paste(vn, mode, sep = "_")
    run_matching(cfg_v, mode, by_country = FALSE, run_name = main)
    run_matching(
      cfg_v,
      mode,
      by_country = TRUE,
      run_name = paste0(main, "_cty")
    )
    summarise_run(main, paste0(main, "_cty"), cfg_v)
  }
}
# Compare base_lfl vs base_env for grp 1 to see how much availability depends on cover type.
# Shared and sole-supplier tallies use candidate lists capped at k; raise k to firm them up.
