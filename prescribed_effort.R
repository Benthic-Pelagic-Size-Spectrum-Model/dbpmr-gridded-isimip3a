# prescribed_effort.R -- per-cell PRESCRIBED fishing effort for region L (replaces spatial gravity).
#
# Instead of re-allocating the region's effort by an ideal-free "gravity" rule (which co-locates
# effort with model biomass and inflates gridded catch ~3.3x), the per-cell effort SHARES are read
# directly from the ISIMIP3a gridded effort product (industrial + artisanal, 0.5deg), aggregated to
# the 1deg mask cells and normalised to sum to 1 within the region. The regional pelagic/benthic
# (U/V) effort MAGNITUDE comes from effort_split_lme.csv (Rousseau regional effort). Spin-up years
# (<1961, no gridded effort) hold the earliest gridded pattern (1961), scaled by each year's regional
# total. Cells are then UNCOUPLED (effort no longer depends on biomass).
#
# COORDINATE CONVENTION: the 1deg FAO-LME mask labels lon by the cell's WESTERN EDGE (integer) and
# lat by the cell CENTRE (x.5). floor(lon),floor(lat) maps BOTH the 0.5deg effort subcells and the
# 1deg mask cell to the same integer SW-corner key, so the aggregation is exact (median 4 subcells
# per 1deg cell). Do NOT use round() here -- it splits a cell's two lon subcells across two keys.
#
# build_prescribed_effort(reg, yrs, ind_csv, art_csv, split_csv, L) -> list(S, EU, EV)
#   reg : data.frame with Lon, Lat (mask cell labels), in the ORDER the caller runs cells
#   yrs : integer model years (e.g. 1841:2010)
#   S   : matrix [nrow(reg) x length(yrs)] spatial shares -- each column sums to 1 (or all-0 if the
#         region has no gridded effort that year)
#   EU,EV: numeric[length(yrs)] normalised regional pelagic(U)/benthic(V) effort (max(EU+EV)=1;
#          absolute level is absorbed by the region-level q rescale)
# Caller sets per-cell fishing mortality:  F_U_i(y)=qpel*N*S[i,y]*EU[y] ; F_V_i(y)=qben*N*S[i,y]*EV[y]
#   (mean over cells of N*S = 1, so mean F_U = qpel*EU -- same scale the 0-D q was fit on.)

build_prescribed_effort <- function(reg, yrs, ind_csv, art_csv, split_csv, L) {
  stopifnot(all(c("Lon","Lat") %in% names(reg)))
  fk <- function(lon, lat) paste(floor(lon), floor(lat))        # shared 1deg SW-corner key (see header)
  rk <- fk(reg$Lon, reg$Lat); N <- nrow(reg)

  ## --- spatial shares from gridded effort (industrial + artisanal) ---
  # 6.8M-row files -> use a fast reader (data.table, else arrow), keep only the 4 needed columns
  rd <- function(f) {
    if (requireNamespace("data.table", quietly = TRUE))
      as.data.frame(data.table::fread(f, select = c("Lat","Lon","Year","NomActive")))
    else { d <- as.data.frame(arrow::read_csv_arrow(f)); d[, c("Lat","Lon","Year","NomActive")] }
  }
  e  <- rbind(rd(ind_csv), rd(art_csv))
  e$fk <- fk(e$Lon, e$Lat)
  e  <- e[e$fk %in% rk, ]                                        # this region's cells only
  gy <- sort(unique(e$Year))                                     # gridded-effort years (1961-2010)
  S  <- matrix(0, N, length(yrs), dimnames = list(NULL, yrs))
  for (y in gy) {
    sub <- e[e$Year == y, ]
    ag  <- tapply(sub$NomActive, sub$fk, sum, na.rm = TRUE)
    v   <- as.numeric(ag[rk]); v[!is.finite(v)] <- 0
    col <- match(y, yrs)
    if (!is.na(col) && sum(v) > 0) S[, col] <- v / sum(v)
  }
  ## spin-up / any pre-gridded year: hold the earliest gridded pattern (1961)
  y0 <- min(gy); pat <- S[, match(y0, yrs)]
  for (iy in seq_along(yrs)) if (yrs[iy] < y0) S[, iy] <- pat

  ## --- regional U/V effort magnitude (Rousseau; effort_split_lme) ---
  sp  <- read.csv(split_csv); sp <- sp[sp$LME == L, ]
  EUy <- tapply(sp$effort[sp$class == "U"], sp$Year[sp$class == "U"], sum)
  EVy <- tapply(sp$effort[sp$class == "V"], sp$Year[sp$class == "V"], sum)
  EU  <- as.numeric(EUy[as.character(yrs)]); EU[!is.finite(EU)] <- 0
  EV  <- as.numeric(EVy[as.character(yrs)]); EV[!is.finite(EV)] <- 0
  denom <- max(EU + EV); if (!is.finite(denom) || denom <= 0) denom <- 1
  list(S = S, EU = EU / denom, EV = EV / denom)
}
