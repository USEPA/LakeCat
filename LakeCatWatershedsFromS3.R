library(DBI)
library(duckdb)
library(sf)

.configure_s3_spatial <- function(con, anonymous = TRUE, region = "us-east-1") {
  # S3 access
  dbExecute(con, "INSTALL httpfs; LOAD httpfs;")
  dbExecute(con, sprintf("SET s3_region='%s';", region))
  dbExecute(con, "SET s3_use_ssl=1;")
  dbExecute(con, "SET s3_url_style='path';")
  if (isTRUE(anonymous)) {
    dbExecute(con, "SET s3_access_key_id='';")
    dbExecute(con, "SET s3_secret_access_key='';")
    dbExecute(con, "SET s3_session_token='';")
  }
  # Geometry support
  dbExecute(con, "INSTALL spatial; LOAD spatial;")
}

.as_sf <- function(df, crs = 4326) {
  if (!nrow(df)) {
    return(sf::st_sf(sf::st_sfc(), crs = crs)[, 0])
  }
  geom_wkb <- structure(df$geometry, class = "WKB")
  sfc <- tryCatch(sf::st_as_sfc(geom_wkb, EWKB = TRUE, crs = crs),
                  error = function(e) sf::st_as_sfc(geom_wkb, EWKB = FALSE, crs = crs))
  sf::st_sf(df[setdiff(names(df), "geometry")], geometry = sfc, crs = crs)
}

# Resolve requested names against actual schema (case-insensitive), ignore extras not found
.resolve_cols <- function(con, pattern, id_col, geometry_col, extras) {
  probe <- dbGetQuery(con, sprintf("SELECT * FROM parquet_scan('%s') LIMIT 0", pattern))
  avail <- names(probe)
  to_map <- function(x) {
    idx <- match(tolower(x), tolower(avail))
    if (is.na(idx)) NA_character_ else avail[idx]
  }
  id_real   <- to_map(id_col)
  geom_real <- to_map(geometry_col)
  extras_real <- vapply(extras, to_map, character(1))
  missing_extras <- extras[is.na(extras_real)]
  extras_real <- extras_real[!is.na(extras_real)]
  
  if (is.na(id_real)) stop(sprintf('id_col "%s" not found. Available: %s', id_col, paste(avail, collapse=", ")), call. = FALSE)
  if (is.na(geom_real)) {
    # Try common fallbacks
    fb <- c("wkb_geometry", "geom")
    fb_idx <- match(tolower(fb), tolower(avail))
    if (any(!is.na(fb_idx))) {
      geom_real <- avail[fb_idx[which(!is.na(fb_idx))[1]]]
      warning(sprintf('geometry_col not found as "%s". Using "%s".', geometry_col, geom_real), call. = FALSE)
    } else {
      stop(sprintf('geometry_col "%s" not found. Available: %s', geometry_col, paste(avail, collapse=", ")), call. = FALSE)
    }
  }
  if (length(missing_extras)) warning(sprintf('Ignoring extras not present: %s', paste(missing_extras, collapse=", ")), call. = FALSE)
  list(id = id_real, geom = geom_real, extras = extras_real)
}

# Fast path when you KNOW the HUC2
read_one_feature_s3_huc2_duckdb <- function(
    huc2,
    comid,
    bucket_root = "s3://dmap-data-commons-ow/data/streamcat/LakeCatWatersheds",
    extras = character(),
    id_col = "COMID",
    geometry_col = "geometry",
    limit = 1,
    anonymous = TRUE,
    region = "us-east-1",
    crs = 4326
) {
  con <- dbConnect(duckdb::duckdb())
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  .configure_s3_spatial(con, anonymous = anonymous, region = region)
  
  pattern <- sprintf("%s/HUC2=%s/*.parquet", bucket_root, huc2)
  cols <- .resolve_cols(con, pattern, id_col, geometry_col, extras)
  
  # Build SELECT list: ST_AsWKB(geom) AS geometry, id, extras
  geom_expr <- sprintf('ST_AsWKB("%s") AS geometry', cols$geom)
  other_cols <- unique(c(cols$id, cols$extras))
  other_expr <- if (length(other_cols)) paste(sprintf('"%s"', other_cols), collapse = ", ") else ""
  select_list <- paste(c(geom_expr, other_expr), collapse = if (nzchar(other_expr)) ", " else "")
  
  value_sql <- if (is.numeric(comid)) as.character(comid) else dbQuoteString(con, as.character(comid))
  
  sql <- sprintf(
    "SELECT %s
     FROM parquet_scan('%s')
     WHERE \"%s\" = %s
     %s",
    select_list, pattern, cols$id, value_sql,
    if (!is.null(limit) && limit > 0) sprintf("LIMIT %d", limit) else ""
  )
  
  df <- dbGetQuery(con, sql)
  g <- .as_sf(df, crs = crs)
  if (nrow(g)) g$HUC2 <- huc2
  g
}

# If you DO NOT KNOW the HUC2: try each partition in order and stop at first hit
read_one_feature_s3_find_duckdb <- function(
    comid,
    huc2_candidates = sprintf("%02d", 1:21),
    bucket_root = "s3://dmap-data-commons-ow/data/streamcat/LakeCatWatersheds",
    extras = character(),
    id_col = "COMID",
    geometry_col = "geometry",
    limit = 1,
    anonymous = TRUE,
    region = "us-east-1",
    crs = 4326
) {
  for (huc2 in huc2_candidates) {
    g <- tryCatch(
      read_one_feature_s3_huc2_duckdb(
        huc2 = huc2, comid = comid, bucket_root = bucket_root, extras = extras,
        id_col = id_col, geometry_col = geometry_col, limit = limit,
        anonymous = anonymous, region = region, crs = crs
      ),
      error = function(e) NULL
    )
    if (!is.null(g) && nrow(g)) return(g)
  }
  # Empty result
  .as_sf(data.frame(geometry = I(list())), crs = crs)[, 0]
}

g <- read_one_feature_s3_huc2_duckdb(
  huc2 = "01",
  comid = 19334077,
  extras = c("GNIS_NAME"),  # use actual column names from that partition's schema
  anonymous = TRUE,
  region = "us-east-1",
  crs = 4326
)
print(g)

g <- read_one_feature_s3_find_duckdb(
  comid = 6779567,
  extras = c("GNIS_NAME"),
  anonymous = TRUE,
  region = "us-east-1",
  crs = 4326
)
print(g)
plot(g$geometry)

