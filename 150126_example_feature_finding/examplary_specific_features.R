# exemplary script to look for specific features #
# examples used here are: a) bacterial features from a known list
                        # b) lipids with with masses of x.4, x.5, x.6, x.7

# Set the tolerance to which matches are looked for. Opposite to alignment this tolerance is the half window. Therefore the standard would be alignment_tolerance/2
matching_tolerance <- as.numeric(alignment_tolerance/2)
# Set method for scaling
scale_method <- 'quantile'
# set TRUE for mode that experiment was measured in.
pos_mode <- FALSE
neg_mode <- TRUE
#### a) bacterial features from a known list ####
# 1. define list: The example list is from an excel document

# read in an excel sheet where all bacterial lipids can be found. If the package has not been installed, use: install.packages('readxl')
# adjust starting column  and end accordingly
library(readxl)
backies <- read_excel('C:/Users/COMPUTER/Documents/Detectability_bactfeat_MPIMM_584.xlsx')
mz_backies <- as.numeric(backies$mz[1: (nrow(backies)-6)]) # this column has the name 'mz' and has non-important information in the last 6 rows.

# 2. find matches from the list inside the experiment
find_mz_matches <- function(mz_library, mz_find, tolerance = 5L) {
  # Compute lower and upper tolerance bounds
  tols <- (tolerance / 1e6) * mz_library
  lows <- mz_library - tols
  ups  <- mz_library + tols
  # Use findInterval to get indices in mz_find
  start_idx <- findInterval(lows, mz_find, left.open = TRUE) + 1 # returns value that is <  lows before values are higher
  end_idx   <- findInterval(ups, mz_find)                        # returns value that is <= ups before values are higher
  # Count occurrences
  ocs <- pmax(end_idx - start_idx + 1, 0)
  matches <- c(ocs > 0)
  return(matches)
}

mzs <- mz(mse) # list of mz values of the experiment
bac_matches <- find_mz_matches(mzs, mz_backies, matching_tolerance)
cat('sum of matches found: ')
sum(bac_matches)
cat('\nmz values of the matches are saved as bac_mz. \n')
bac_mz <- mzs[bac_matches]
no_matches <- mz_backies[!find_mz_matches(mz_backies, bac_mz, matching_tolerance)]
cat('no match was found for:\n', no_matches)

# how many pixels do the matching features occupy
inten <- intensity(mse)
bac_inten <- inten[bac_matches,]
pix1 <- rowSums(bac_inten >0)
pix1 <- unlist(pix1, use.names = FALSE)
cat('The matched masses occur in:\n', pix1, '\npixel\n')

# This function only works efficiently for smaller subsets.
# Otherwise scale_image can be used for every image individually.
.scale_row <- function(intensities, method, percent) {

  if (method == 'log') {
    logs <- log1p(intensities)
    return(logs)
  }

  if (method == '%') {
    max_int <- max(intensities)
    quant_int <- max_int * (percent / 100)
    intensities[intensities >= max_int] <- quant_int
    return(intensities)
  }

  if (method == 'quantile') {
    sorted <- sort(intensities)
    all_int <- sorted[sorted > 0]
    upper_q <- ceiling(0.75 * length(all_int))
    quan <- all_int[upper_q]
    intensities[intensities > quan] <- quan
    return(intensities)
  }
}
scale_all <- function(experiment, method = c('log', '%', 'quantile'), percent = NULL, tolerance = 5) {

  method <- match.arg(method)

  tics <- pData(experiment)$'TIC'

  if (method == 'log') {
    cat('selected method: log \n')
    pData(experiment)$TIC <- log1p(tics)
  }

  if (method == '%') {
    if (is.null(percent)) stop('select a percentage: percent= x')
    cat('selected method:', percent, '% \n')
    max_int <- max(tics)
    quant_int <- max_int * (percent/100)
    tics[tics >= max_int] <- quant_int
    pData(experiment)$TIC <- tics
  }

  if (method == 'quantile') {
    cat('selected method: quantile \n')
    sorted <- sort(tics)
    upper_q <- round(0.75 * length(tics))
    quan <- tics[upper_q]
    tics[upper_q:length(tics)] <- quan
    pData(experiment)$TIC <- tics
  }

  spec <- spectra(experiment)[]

  for (i in seq_len(nrow(spec))) {
    spec[i, ] <- .scale_row(spec[i, ], method, percent)
  }

  experiment@spectraData$'intensity' <- spec
  return(experiment)
}
mse_bac <- subsetFeatures(mse, bac_matches)

# scale the bacteria to enhance visibility
mse_scale <- scale_all(mse_bac, method=scale_method)

### imaging: create a pdf file with all images of bacterial features. ###
date <- Sys.Date()
file_name <- paste0(date,'_bac_lipids_nt_', round(noise_threshhold), '.pdf')
pdf(file = file_name, width =6, height = 4)
for (uz in bac_mz) {
  p <-image(mse_scale, mz=uz, style='dark', ylab="intensity", xlab= 'm/z')
  print(p) # allows for the MSI to be "printed" onto the PDF
}
dev.off()
rm(p)

#### b) finding lipids
look_for <- function(mz_list, what) {
  mz_character <- as.character(mz_list)
  where <- vapply(mz_character, function(mz) grepl(what, mz), logical(1))
  return(unname(where))
} # function to look for specific patterns
mzs <- mz(mse) #'\d*\.98\d*'
lipids <- look_for(mzs, '\\d*\\.[456]\\d*' )

lipidsTRUE <- which(unname(lipids))
mse_lipids <- subsetFeatures(mse, lipidsTRUE)
mse_lipids <- summarizePixels(mse_lipids, stat=c(TIC='sum'))
mzlip <- mz(mse_lipids)

### network analysis ###
# https://www.acdlabs.com/blog/common-adduct-and-fragment-ions-in-mass-spectrometry/
H <- 1.0078250319
C <- 12.00
Na= 	22.98976928
K =  38.96370648
C2H4 <- C+C+H+H+H+H
CH2<- C+H+H
Cl <- 34.969402
Cl_iso <- 36.969402
CHO2 <- 44.998201
CH3O2 <- 59.013851
Br <- 78.918885
Br_iso <- 80.918885
CF3CO2 <- 112.985586

# mzlip_r4 <- round(mzlip, 4)
# total_diffs <-outer(mzlip, mzlip, "-")
# colnames(total_diffs) <- mzlip_r4
# rownames(total_diffs) <- mzlip_r4
# diff_mat <- round(total_diffs, 2)

positive_ions <- list('[M+H]'= 1.0078250319,
                 '[M+Na]'=22.98976928,
                 '[M+K]'=38.96370648)
negative_ions <- list('[M-H]'= 1.0078250319,
                 '[M+Cl]'= 34.968852682,
                 '[M+Cl_iso]'= 36.965902602,
                 '[M+CHO2]'=44.998201,
                 '[M+CH3CO2]' =59.013851,
                 '[M+Br]'=78.9183376,
                 '[M+Br_iso]'=80.9162897)


adducts <- function(subset, ppm, pos_ions, neg_ions) {
  if (class(subset) == 'MSImagingExperiment') {
    mz_list <- round(mz(subset),4)
  }
  else {mz_list <- round(subset,4)}

  diff_mat <- outer(mz_list, mz_list, "-")
  diff_mat <- round(diff_mat, 2)
  result_mat <- NULL
  cnames <- c()
  rnames <- round(mz_list,4)
  # test ionization
  if (pos_mode){
    for (g in 1:(length(pos_ions)-1)) {
      for (i in (g+1):length(pos_ions)) {
        vec_ind <- which((-1*diff_mat) == round((pos_ions[[i]]-pos_ions[[g]]),2), arr.ind = TRUE)
        if (nrow(vec_ind) == 0) next
        vec <- rep(NA_real_, length(mz_list))
        vec[vec_ind[, 1]] <- mz_list[vec_ind[,2]]
        toler <- rep(NA_real_, length(mz_list))
        tv <- mz_list[vec_ind[,1]] + pos_ions[[i]]-pos_ions[[g]]
        toler[vec_ind[, 1]] <- ((mz_list[vec_ind[,2]]-tv)/tv)*1e6
        vec[abs(toler)>ppm] <- NA_real_
        toler[(abs(toler)>ppm)] <- NA_real_             # filters out values with deviation above tolerance.
        result_mat <- cbind(result_mat, vec, toler)
        cnames <- c(cnames, names(pos_ions)[i], 'tolerance [ppm]')
      }
    }
  }
  if (neg_mode){
    ### for [M-H] adducts
    for (q in 2:length(neg_ions)) {
      vec_ind <- which((-1*diff_mat) == round((neg_ions[[q]]+ neg_ions[[1]]),2), arr.ind = TRUE)
      if (nrow(vec_ind) == 0) next
      vec <- rep(NA_real_, length(mz_list)) # creates a list the same length as mz values in experiment with NA.
      vec[vec_ind[, 1]] <- mz_list[vec_ind[,2]] #replaces NA with value for matches
      toler <- rep(NA_real_, length(mz_list))   # same process but for tolerances
      tv <- mz_list[vec_ind[,1]] + neg_ions[[q]] + neg_ions[[1]]
      toler[vec_ind[, 1]] <- ((mz_list[vec_ind[,2]]-tv)/tv)*1e6 # calculation of tolerance in ppm
      vec[abs(toler)>ppm] <- NA_real_
      toler[(abs(toler)>ppm)] <- NA_real_             # filters out values with deviation above tolerance.
      result_mat <- cbind(result_mat, vec, toler)
      name <- paste0(names(neg_ions)[q], '_', names(neg_ions[1]))
      cnames <- c(cnames, name, 'tolerance [ppm]')
      # rnames[vec_ind[,1]] <- paste0(rnames[vec_ind[,1]], names(neg_ions[1]))
    }
    ## for all added adducts [M+]
    for (g in 2:(length(neg_ions)-1)) {
      for (i in (g+1):length(neg_ions)) {
        vec_ind <- which((-1*diff_mat) == round((neg_ions[[i]]-neg_ions[[g]]),2), arr.ind = TRUE)
        if (nrow(vec_ind) == 0) next
        vec <- rep(NA_real_, length(mz_list)) # creates a list the same length as mz values in experiment with NA.
        vec[vec_ind[, 1]] <- mz_list[vec_ind[,2]] #replaces NA with value for matches
        toler <- rep(NA_real_, length(mz_list))   # same process but for tolerances
        tv <- mz_list[vec_ind[,1]] -neg_ions[[g]] + neg_ions[[i]]
        toler[vec_ind[, 1]] <- ((mz_list[vec_ind[,2]]-tv)/tv)*1e6 # calculation of tolerance in ppm
        vec[abs(toler)>ppm] <- NA_real_
        toler[(abs(toler)>ppm)] <- NA_real_             # filters out values with deviation above tolerance.
        result_mat <- cbind(result_mat, vec, toler)
        name <- paste0(names(neg_ions)[i], '_', names(neg_ions[g]))
        cnames <- c(cnames, name, 'tolerance [ppm]')
      }
    }
  }
  colnames(result_mat) <- cnames
  rownames(result_mat) <- rnames
  return(result_mat)
}
net <- adducts(mzlip,5, positive_ions, negative_ions)

# occurences <- function(mz_list, tab) {
#   mz_vals <- round(mz_list, 4)
#   net_r  <- round(tab, 4)
#   row_hits_mat <- sapply(mz_vals, function(mz) {
#     rowSums(net_r == mz, na.rm = TRUE) > 0
#   })
#
#   dim(row_hits_mat)
#   # nrow(net) x length(mz_vals)
#   return(row_hits_mat)
# }
# occs <- occurences(mzlip, net)

add_occured<- function(tab, mz_list, occs) {
  mz_list <- round(mz_list, 4)

  # create empty list to store rows for each m/z
  res <- vector("list", length(mz_list))
  names(res) <- mz_list

  for (col in seq_len(ncol(occs))) {
    # find rows where this mz occurs
    rows <- which(occs[, col])
    if (length(rows) > 0) {
      # store those rows as a sub-matrix
      res[[col]] <- tab[rows, , drop = FALSE]
    }
  }

  return(res)
}

occurence_list <- add_occured(net, mzlip,occs)


.diag_list <- function(n) {
  ma <- matrix(NA, n,n)
  idx <- which(upper.tri(ma, diag = TRUE), arr.ind = TRUE)
  # Order row-wise: first by row, then by column
  ord <- order(idx[, 1], idx[, 2])
  ma[idx[ord, ]] <- seq(1, n*(n+1)/2)
  return(as.data.frame(ma))
}


adduct_prob <- function(tab, occ_list, mzlist, pos_ions, neg_ions, pos_mode, neg_mode) {
  result <- data.frame(mz = numeric(), ion = character(), prob= numeric(), relations = I(list()))
  mzlist <- round(mzlist,4)
  rem <- seq(1, ncol(tab), by=2)
  tab_no_tol <- tab[, rem]
  if (neg_mode){
    ions <- neg_ions
  }
  else if (pos_mode){
    ions <- pos_ions
  }
  slices <- cumsum((length(ions)-1):1)
  ma <- .diag_list(length(ions)-1)
  for (pos in seq_len(nrow(tab))) {
    if (sum(!is.na(tab[pos, ])) != 0) {
      # tab occurences
      nona <- which(!is.na(unname(tab_no_tol[pos, ])))
      most <- findInterval(nona, slices, left.open = TRUE)+1
      # in other rows occurences as
      if (!is.null(occ_list[[pos]])) {
        ox <- (unname(which(occ_list[[pos]]== mzlist[pos], arr.ind = TRUE)[,2])+1)/2
        occin <- sapply(ox, function(x) unname(which(ma == x, arr.ind = TRUE)[,2]))
      }
      else {occin <- c()}

      if (length(most)+length(occin) == 0) next
      together <- prop.table(table(c(most, occin)))
      max_p <- max(together)
      posion <- as.integer(names(together)[together==max_p])

      for (pi in posion) {
        io <- names(ions)[pi]
        if (pi == 1) {
          nona_filt <- nona[nona >0 & nona <= slices[pi]]
        }
        else {
        nona_filt <- nona[nona > slices[pi-1] & nona <= slices[pi]]
        }

      related_mzs <- c(unname(tab_no_tol[pos, nona_filt]), as.numeric(rownames(try[[pos]])))
      result <- rbind(result, data.frame(mz= mzlist[pos], ion = io, prob = round(max(together),2), relations = I(list(related_mzs))))
      }
    }
  }
  return(result)
}

adduct_probabilities <- adduct_prob(net, occurence_list, mzlip, positive_ions, negative_ions, pos_mode, neg_mode)


