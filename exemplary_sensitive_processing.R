# sensitive processing pipeline.

#### set paths to define in which directory data is stored. ####
  # windows directory separations like '\' instead of '/' will be converted automatically when copying the path into the brackets of the raw string.

# 1. path to the working directory in which created data will be stored.
wd_path <- r'(C:\Users\COMPUTER\Documents\MALDI_Daten)'
# 2. path to the MALDI data file (.imzML). If the file is stored in the working directory the specifying ending is enough.
msa_path <- r'(20251023_MST_MS161_mY9colon_slide430_DAN_15um_A38_neg_400-1200_461x436.imzml)'

#### set processing parameter ####
# 1. set the alignment_tolerance in ppm. Standard for MALDI experiment is 10 ppm
alignment_tolerance <- 10L
# 2. Set the maximum number of pixels in which a feature may occur to be removed by strict filtering.
strict_pix <- 3L
# 3. Set the maximum number of pixels in which a feature may occur to be removed by neighbor filtering.
nb_pixel <- 100L
# 4. Set the proximity for pixels to be considered neighbors.
nb_proxim <- 2L

#### From this point on the script won't need adjustment but can be adjusted for parameters. After selecting the noise threshold in the opened window and closing it, no supervison needed. ####
#### load needed packages ####
# If certain packeges need to be installed, use the following argument: install.packages(c('ggplot2', 'plotly', 'shiny', 'Cardinal'))
library(Cardinal)
library(shiny)
library(ggplot2)
library(plotly)

#### functions ####
path_convert <- function(file_path) {
  gsub("\\\\", "/", file_path)
}

wd_path <- path_convert(wd_path)
setwd(wd_path)

path_msa <- path_convert(msa_path)
msa<- readMSIData(path_msa)

#### converting the msa data into usable data frame ####
combine_msa <- function(experiment, method= c('mean', 'median'), round_by = 3L){
  cat("Combining spectra .")
  # extracts intensity data and mz data from MSI Array
  intensities <- experiment@spectraData@data$intensity
  mz <- experiment@spectraData@data$mz
  # creates lists with all m/z values or all intensity values of all spectra.
  all_mz <- unlist(unname(mz[]))
  cat(" .")
  all_intensities <- unlist(unname(intensities[]))
  # sort intensity list in the order of mz values so thant integers will match after sorting mz values.
  inten_sort <- all_intensities[order(all_mz)]
  mz_sort <- sort(all_mz)
  mz_sort_rounded <- round(mz_sort, round_by)
  cat(" .")
  # create list where integers stay the same when a value is repeated. E.g. c(10,10,10,23,24,24,48,152,152) -> c(1,1,1,2,3,3,4,5,5)
  grp <- with(rle(mz_sort_rounded), rep(seq_along(values), lengths))
  # create mean or median for intensities belonging to same rounded m/z value. listed as duplicates for each m/z value
  if (length(method)>1) {
    method <- method[1]
  }
  inten_sort     <- ave(inten_sort, grp, FUN = method)
  cat(" .")
  # remove consecutive duplicates from m/z list and intensity list
  keep <- c(TRUE, mz_sort_rounded[-1] != mz_sort_rounded[-length(mz_sort_rounded)])
  inten_sort_short      <- inten_sort[keep]
  mz_sort_rounded_short <- mz_sort_rounded[keep]

  msa_df <- data.frame('mz'= mz_sort_rounded_short, "intensity"=inten_sort_short)
  return(msa_df)
}
msa_df <- combine_msa(msa, round_by = 3L)

# use function to select noise threshold in the mean/median msa plot
spectraViewer <- function(experiment, i=NULL, xcol = "mz", ycol = "intensity") {
  #### transcribing MSI Experiment data into data frame ####
  if (type(experiment) != 'double'){
    y_intensities <- c()
    if (is.null(i)) {
      i = 1
    }
    if (is.numeric(i)){
      y_intensities <- as.numeric(spectra(experiment)[,i])
    }
    else if (is.character(i) && i %in% colnames(fData(experiment))) { # <<< fix: check column exists
      y_intensities <- as.numeric(fData(experiment)[[i]])
    } else {
      stop("i must be numeric index of a spectrum or 'mean' ")
    }
    df <- setNames(
      data.frame(
        as.numeric(mz(experiment)),
        y_intensities
      ),
      c(xcol, ycol)
    )
  }
  else if (type(experiment) == 'double') {
    df <- experiment
  }

  plot_label <- if (is.numeric(i)) {
    paste("pixel", i)
  } else if (is.character(i)) {
    paste("mean spectrum", i)
  } else {
    ""
  }

  #### layout of window that opens for interactive plot ####
  ui <- fluidPage(
    tags$head(
      tags$script(HTML("
      document.addEventListener('keydown', function(e) {
        if (e.key === 'Enter') {
          Shiny.setInputValue('enter_pressed', Math.random());
        }
      });
    "))
    ),
    style = "padding:0; margin:0;",
    plotlyOutput("msaPlot", height = "80vh"),
    verbatimTextOutput("hoverText"),
    h5("Saved Value (Click):"),
    verbatimTextOutput("clickedText")
  )

  #### creation of plot and programming cursor functions ####
  server <- function(input, output, session) {
    # Stores last clicked y-value inside Shiny session
    saved_val <- reactiveVal(NULL)
    ####  Render plot ####
    output$msaPlot <- renderPlotly({
      p <- ggplot(df, aes(.data[[xcol]], .data[[ycol]])) +
        geom_line(linewidth = 0.5) +
        theme_minimal() +
        labs(x = xcol, y = ycol, title = plot_label)
      # interactive through plotly
      ggplotly(p, tooltip = NULL) %>%
        layout(
          hovermode = "closest",

          xaxis = list(showspikes = FALSE),
          # horizontal red cursor line
          yaxis = list( showspikes      = TRUE,
                        spikemode       = "across",
                        spikesnap       = "cursor",
                        spikethickness  =  2,
                        spikecolor      = 'violet',
                        spikedash       = 'solid'
          ),
          # Clean compact margins
          margin = list(l = 40, r = 10, b = 40, t = 10)
        ) %>%
        config(displayModeBar = TRUE) %>%
        event_register("plotly_hover") %>%
        event_register("plotly_click")
    })
    ####  Hover text below plot ####
    output$hoverText <- renderPrint({
      d <- event_data("plotly_hover")
      if (is.null(d)) {
        cat("Move cursor over plot to see", xcol, "&", ycol)
      } else {
        cat(xcol, ": ", d$x, "\n", sep = "")
        cat(ycol, ": ", d$y, sep = "")
      }
    })
    ####  Save clicked Y-value ####
    observeEvent(event_data("plotly_click"), {
      click <- event_data("plotly_click")
      if (!is.null(click)) {
        saved_val(click$y)           # Save inside Shiny
        noise_threshhold <<- click$y # Save globally (persistent)
      }
    })
    ####  Show saved value ####
    output$clickedText <- renderPrint({
      if (is.null(saved_val())) {
        cat("Click on the plot to save a y-value.")
      } else {
        cat("noise_threshhold: ", saved_val(), "\n")
        cat("Press Enter or close window to save chosen noise threshhold")
      }
    })

    #### Close window on Enter ####
    observeEvent(input$enter_pressed, {
      stopApp()
    })
  }
  shinyApp(ui, server)
}
spectraViewer(msa_df)

# remove m/z values and their respective intensity that is below the noise threshold.
# When not using spectraViewer for  selecting the noise threshold instead use: noise_threshhold <- X
remove_noise <- function(experiment, noise) {
  cat("removing noise signals ...")
  filtered_mz <- mapply(
    function(mz, int) mz[int >= noise],
    experiment@spectraData@data$mz,
    experiment@spectraData@data$intensity,
    SIMPLIFY = FALSE
  )
  filtered_intensity <- mapply(
    function(int, thr) int[int >= thr],
    experiment@spectraData@data$intensity,
    MoreArgs = list(thr = noise),
    SIMPLIFY = FALSE
  )
  msa_no_noise <- experiment
  msa_no_noise@spectraData@data$intensity <- filtered_intensity
  msa_no_noise@spectraData@data$mz <- filtered_mz
  return(msa_no_noise)
}
msa_no_noise <- remove_noise(msa, noise_threshhold)

## clean up R environment ##
rm(msa, path_msa, msa_df, spectraViewer, combine_msa ,remove_noise)
gc()

# align msa spectra with set tolerance
mse <- normalize(msa_no_noise, method='tic')
mse <- peakAlign(mse, tolerance= alignment_tolerance, units= 'ppm')

#### filtering sparse signals ####
# strict filter to reduce noise befor summarizing:
pixel_filter <- function(experiment, min_count_pix) {
  cat('removing pixel')
  spectra_num <- experiment@spectraData@data@listData$intensity@dim[2]
  count <-experiment@featureData$count
  count_pix <- c(count <= min_count_pix)
  if (!any(count_pix)) {
    return('there are no features appearing in less or equal amount of given minimal pixel')
  }
  # creating an experiment that only holds features with count >= selected pixel
  mse_count <- subsetFeatures(experiment, count_pix)
  gc()
  non_empty_pix <- colSums(mse_count@spectraData$intensity[])> 0
  mse_count <- subsetPixels(mse_count, non_empty_pix)
  cat(' .')
  int_mat <- round(mse_count@spectraData$intensity[],4)
  cat(' .')
  # calculating max intensity value of a feature across all spectra.
  row_sums <- rowSums(int_mat>0)
  proof <- row_sums <= min_count_pix
  # comparing whether max value divided by all spectra is equal to the mean. If so, that feature only occurs in one pixel.
  cat(' .')
  mse_count <- subsetFeatures(mse_count, proof) # only contains features that occur in one pixel
  tossed <- !(mz(experiment) %in% mz(mse_count))
  mse_1pp <- subsetFeatures(experiment, tossed)
  return(mse_1pp)
}
mse <- pixel_filter(mse, strict_pix)
gc()

#### remove pixel without neighbors ####
chunks <- function(experiment, num) {
  feat_num <- experiment@spectraData@data@listData$intensity@dim[1]
  spec_num <- experiment@spectraData@data@listData$intensity@dim[2]
  chunk_id <- cut(seq_len(feat_num), num, labels = FALSE) # cut creates integer list 1 to 4 with equal amounts of features in seq_len from feat_num
  starts <- tapply(seq_len(feat_num), chunk_id, min) # looks for smallest index in each chunk
  ends   <- tapply(seq_len(feat_num), chunk_id, max) # looks for largest index in each chunk
  exp_list <- list()
  for (i in 1:num) {
    subsetf <- subsetFeatures(experiment, (starts[[i]]:ends[[i]]))
    exp_list[length(exp_list) + 1] <- list(subsetf)
  }
  return(exp_list)
} # helper function for pixel_filter_nb
pixel_filter_nb <- function(experiment, min_pixel, proxim_range, chunked=FALSE, chunk_amount = 4L){
  FUN <- function(experiment, sp) {
    row_counts <- rowSums(sp > 0)
    small_10 <- row_counts <= min_pixel
    mse_s10 <- subsetFeatures(experiment, small_10) # generates a data set only holding those features which occur 10 times or less
    empty_pixel <- c(colSums(mse_s10@spectraData$intensity[]) > 0)
    mse_s10 <- subsetPixels(mse_s10, empty_pixel)
    ### Now features that were previosly discarded are checked if any of these "cluster" by checking the distance
    # to all pixels holding any given feature

    #generates a list of the length of the features to check holding all the non 0 positions
    spec <- spectra(mse_s10)
    nro_features <- nrow(mse_s10@featureData)
    n0_features_position <- apply(spec[1:nro_features, ], 1, function(row) which(row > 0))

    # Checking if any of the non 0'S for each feature are neighbors by a defined radius
    xy_cord <- data.frame(x = pData(mse_s10)$x, y = pData(mse_s10)$y) # saves xy cords separate to reduce computational time
    distance_test_bool <-  vapply(seq_along(n0_features_position), function(i) {
      idx <- n0_features_position[[i]]
      return(any(dist(xy_cord[idx, ]) <= proxim_range))
    },
    logical(1) )

    idx_small_10_true <- which(small_10) # holds the indexes for all values smaller n pixels in the small_10 list
    idx_small_change <- idx_small_10_true[distance_test_bool] # only holds the indexes which are smaller n but also have neighbors
    small_10co <- small_10 #generates a new list holding the same information as small_10
    small_10co[idx_small_change] <- FALSE # changes the TRUE value of all features with neighbors but smaller n to FALSE
    big_10co <- !small_10co #inverts the whole frame for removal of said feature
    return(big_10co)
  }
  # To chunk or not to chunk
  feat_num <- experiment@spectraData@data@listData$intensity@dim[1]
  spec_num <- experiment@spectraData@data@listData$intensity@dim[2]
  if (isFALSE(chunked)) {
    if (((spec_num/1e4) * (feat_num/1e4)) > 45) {chunked <- TRUE}  # threshold decided by trying different sizes of spectra.
  }
  if (chunked) {
    exp_list <- chunks(experiment, chunk_amount)
    vec <- c()

    for (i in 1:chunk_amount) {
      cat('chunk ', i, '/', chunk_amount, '\n')
      mse1 <- exp_list[[i]]
      sp  <- mse1@spectraData$intensity[]
      vece <- FUN(mse1, sp)
      vec <- c(vec, vece)
    }
  }
  else {
    sp <- experiment@spectraData$intensity[]
    vec <- FUN(experiment, sp)
  }
  experiment <- subsetFeatures(experiment, vec) # removes all pixel smaller 10 which do not have neighbors
  return(experiment)
}
mse <- pixel_filter_nb(mse, nb_pixel, nb_proxim, chunk_amount = 4L)
gc()

# summarizing experiment data for mean and TIC
mse <- summarizeFeatures(mse, stat = 'mean')
gc()
mse<- summarizePixels(mse, stat=c(TIC='sum'))
gc()

