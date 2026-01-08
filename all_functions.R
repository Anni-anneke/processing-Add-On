# function script

# This function creates an dark image of the TIC in a new window in.
pic <- function(experiment) {
  dev.new()
  image(experiment, 'TIC', style= 'dark')}

# This function creates two lists containing m/z values or intensity values of all spectra combined
# It will then order intensities and m/z lists so that m/z values are inclining and intensities are sorted to match their respective m/z value.
# All m/z values are rounded to a decimal number e.g 3 (0.001).
# Duplicates of m/z values are deleted and their intensities are combined either by mean or median.
# The output is a data frame containing a list of the sorted m/z values and of the intensities.
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

# This function opens the shiny app to plot the previously created data frame.
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

# old, needs mean. less effective!
one_pix_remover <- function(experiment) {
  cat('removing pixel')
  count <-experiment@featureData$count
  count_1 <- c(count ==1)
  if (sum(count_1)==0) {
    return('there are no features appearing in only one pixel')
  }
  # creating an experiment that only holds features with count=1
  mse_count1 <- subsetFeatures(experiment, count_1)

  spectra_num <-  mse_count1@spectraData@data@listData$intensity@dim[2]
  means <- round(mse_count1@featureData$mean, 8)
  int_mat <- mse_count1@spectraData$intensity[]
  cat(' .')
  # calculating max intensity value of a feature across all spectra.
  row_max <- apply(int_mat, 1, max)
  cat(' .')
  # comparing wether max value devided by all spectra is equal to the mean. If so, that feature only occurs in one pixel.
  ratio <- round((row_max / spectra_num),8)
  proof <- c(ratio == means)
  cat(' .')
  mse_count1 <- subsetFeatures(mse_count1, proof) # only contains features that occur in one pixel
  tossed <- !(round(mz(experiment),4) %in% round(mz(mse_count1),4))
  mse_1pp <- subsetFeatures(experiment, tossed)
  return(mse_1pp)
}

select_matrix <- function(experiment, corner_size_tl, corner_size_tr, corner_size_bl,corner_size_br) {
  x_list <- pData(experiment)$x
  y_list <- pData(experiment)$y

  xmax <- max(x_list)
  ymax <- max(y_list)

  selection_vec <- (
    # top-left triangle
    (x_list <= corner_size_tl & y_list <= (corner_size_tl - x_list + 1)) |
      # top-right triangle
      (x_list >= (xmax - corner_size_tr +1) & y_list <= (corner_size_tr-(xmax - x_list) + 1)) |
      # bottom-left triangle
      (x_list <= corner_size_bl & y_list >= (ymax - corner_size_bl + x_list)) |
      # bottom-right triangle
      (x_list >= (xmax - corner_size_br +1) & y_list >= (ymax -(corner_size_br-(xmax - x_list)))))

  return(selection_vec)
}

subtract_matrix <- function(experiment, mse_sample, mse_matrix) {
  matrix_mean <- mse_matrix@featureData$mean
  sample_mean <- mse_sample@featureData$mean

  ratio_matrix_sample <- mapply("/", matrix_mean, sample_mean) # calculates the ratio between the mean intensities of the regions selected
  ratio_matrix_sample[is.infinite(ratio_matrix_sample)] <- 10000 # sets all values to 10000 which were calculated as infinite (therefor all values with x/0)
  ratio_matrix_sample[is.na(ratio_matrix_sample)] <- 0 # sets all values to 0 which are calculated as NaN (therefor all values calculated by 0/0)

  exp_vec <- ratio_matrix_sample < 1 # generates a logic vector holding TRUE for all ratios < 1 (up more matrix signal remains)
  exp_total <- subsetFeatures(experiment, exp_vec) # generates a experiments of all features with a small ratio
  exp_total  <-  summarizePixels(exp_total, stat=c(TIC="sum")) # summarizes the TIC for the new experiment
  return(exp_total)
}

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

.interpret <- function(experiment, mz=NULL, i=NULL, tolerance = 5L) {
  mz_list <- mz(experiment)
  if (is.null(mz) & is.null(i))
    stop("You must provide either mz, i, or mz='TIC'")
  if (!is.null(mz) && mz == "TIC" | !is.null(i) && i == "TIC") {
    return(FALSE) }
  if (!is.null(mz) | !is.null(i)) {
    if (!is.null(i) && i > length(mz_list))
      stop("Index is out of range")
    if (!is.null(mz)) {
      i <- which(find_mz_matches(mz_list, as.numeric(mz), tolerance))
      if (length(i) == 0)
        stop("No m/z matched within tolerance. Closest m/z values are: ", mz_list[findInterval(mz, mz_list)], 'and ', mz_list[findInterval(mz, mz_list)+1])
      cat("selected m/z value matched with ", mz_list[i], '\n')
      if (length(i)>1) {
        cat('identifying closest match... \n')
        comp_mz <- mz_list[i]
        diffs <- comp_mz - mz
        mina <- which(diffs== min(diffs))
        if (length(mina)>1) {
          cat(comp_mz[mina], 'all match with same distance. Will work with', comp_mz[mina[1]], '\n')
          i <- mina[1]
        } else {
            i <- i[mina]
            cat('closest match is: ', comp_mz[mina], '\n')
        }
      }
    }
    return(as.integer(i))
  }
}

scale <- function(experiment, mz=NULL, i=NULL, method = c('log', '%', 'quantile'), percent=NULL, tolerance = 5) {
  method <- match.arg(method)
  i <- .interpret(experiment, mz, i)
  if (isFALSE(i)) {
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
      tics[tics >= quant_int] <- quant_int
      pData(experiment)$TIC <- tics
    }
    if (method == 'quantile') {
      cat('selected method: quantile \n')
      sorted <- sort(tics)
      upper_q <- round(0.75*length(tics))
      quan <- tics[upper_q]
      tics[upper_q:length(tics)] <- quan
      pData(experiment)$TIC <- tics
    }
    return(experiment)
  }
  if (is.integer(i)){
    spec <- spectra(experiment)[]
    intensities <- spec[i, ]
    if (method == 'log') {
      cat('selected method: log \n')
      logs <- log1p(intensities)
      spec[i,] <- logs
      experiment@spectraData$'intensity' <- spec
    }
    if (method == '%') {
      if (is.null(percent)) stop('select a percentage: percent= x')
      cat('selected method:', percent, '% \n')
      max_int <- max(intensities)
      quant_int <- max_int * (percent/100)
      intensities[intensities >= quant_int] <- quant_int
      spec[i, ] <- intensities
      experiment@spectraData$'intensity'<- spec
    }
    if (method == 'quantile') {
      cat('selected method: quantile \n')
      sorted <- sort(intensities)
      all_int <- sorted[sorted>0]
      upper_q <- ceiling(0.75*length(all_int))
      quan <- all_int[upper_q]
      intensities[intensities>quan] <- quan
      spec[i, ] <- intensities
      experiment@spectraData$'intensity'<- spec
    }
  return(experiment)
  }
}

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

scale_image <- function(experiment, mz=NULL, i=NULL, method=c('log', '%', 'quantile'), percent=NULL, tolerance = 5) {
  method <- match.arg(method)
  i <- .interpret(experiment, mz, i)
  if (isFALSE(i)) {
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
      tics[tics >= qunat_int] <- quant_int
      pData(experiment)$TIC <- tics
    }
    if (method == 'quantile') {
      cat('selected method: quantile \n')
      sorted <- sort(tics)
      upper_q <- round(0.75*length(tics))
      quan <- tics[upper_q]
      tics[upper_q:length(tics)] <- quan
      pData(experiment)$TIC <- tics
    }
    return(experiment)
  }
  if (is.integer(i)){
    experiment <- subsetFeatures(experiment, i)
    spec <- spectra(experiment)[]
    if (method == 'log') {
      cat('selected method: log \n')
      logs <- log1p(spec)
      experiment@spectraData$'intensity' <- logs
    }
    if (method == '%') {
      if (is.null(percent)) stop('select a percentage: percent= x')
      cat('selected method:', percent, '% \n')
      max_int <- max(spec)
      quant_int <- max_int * (percent/100)
      spec[spec >= quant_int] <- quant_int
      experiment@spectraData$'intensity'<- spec
    }
    if (method == 'quantile') {
      cat('selected method: quantile \n')
      sorted <- sort(spec)
      all_int <- sorted[sorted>0]
      upper_q <- ceiling(0.75*length(all_int))
      quan <- all_int[upper_q]
      spec[spec>quan] <- quan
      experiment@spectraData$'intensity'<- spec
    }
    return(experiment)
  }
}

### R can not read \ as string. It will use the according regular expression. This is why one needs to put r (stands for raw) infront of the string.
# e.g.  r'(C:\Users\MALDI_Data\perfect-experiment.IBD)'
path_convert <- function(file_path) {
  gsub("\\\\", "/", file_path)
}

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

.chunks <- function(experiment, num) {
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
}

## usesmse## uses chunks. Don't know which is better.
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
    exp_list <- .chunks(experiment, chunk_amount)
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

PixelViewer <- function(experiment, mz=NULL, i=NULL, tolerance= 5L) {
  #### creating a data frame from experiment data ####
  i <- .interpret(experiment, mz, i, tolerance)
  if (is.integer(i)) {
    intensities <- spectra(experiment)[i, ]
  }
  if (isFALSE(i)) {
    intensities <- pData(experiment)$TIC
  }
  df <- data.frame(x_coord = pData(experiment)$x, y_coord = pData(experiment)$y, value = intensities)
  plot_label <- if (is.numeric(mz)) {
    paste("i = ", i, '/ m/z = ', round(mz(experiment)[i]),4)
  } else if (is.character(mz)) {
    paste("TIC")
  } else if (is.null(mz) && is.numeric(i)){
    paste("i = ", i)
  } else {
    ""
  }
  #### set up for interactive shiny window ####
  ui <- fluidPage(
    tags$head(
      tags$script(HTML("
      document.addEventListener('keydown', function(e) {
        if (e.key === 'Enter') {
          Shiny.setInputValue('enter_pressed', Math.random());
        }
        if (e.key === 'z') {
          Shiny.setInputValue('z_pressed',Math.random());
        }
      });
    "))
    ),
    style = "padding: 0; margin:0;", # uses full window without boarders
    plotlyOutput("heat_im", height = "80vh"), # plot takes up 80 percent of height.
    # Hover information
    verbatimTextOutput("hoverText"),
    # Saved click value (display only, stored globally)
    h5("Selected pixels (last click):"),
    verbatimTextOutput("clickedText"),
    h5('All selected pixels:'),
    verbatimTextOutput('allClicks')
  )

  #### creation of plot and programming cursor functions ####
  server <- function(input, output, session) {

    # Stores last clicked y-value inside Shiny session
    selected_df <- reactiveVal(data.frame(x=numeric(0), y=numeric(0)))

    ####  Render Plot ####
    output$heat_im <- renderPlotly({
      # ggplot2 heat map
      p <- ggplot(df, aes(x=x_coord, y=y_coord, fill= value)) +
        geom_tile(color= 'black') +
        scale_fill_gradientn(colors = c("#4a0057", '#00598d', '#009b95', '#5dcd63','#eee42d')) +
        scale_y_reverse() +
        coord_fixed()

      # Convert to interactive plotly
      ggplotly(p, tooltip = NULL) %>%
        layout(
          hovermode = "closest",
          xaxis = list(showspikes = FALSE),
          yaxis = list(showspikes = FALSE),
          # Clean compact margins
          margin = list(l = 40, r = 10, b = 40, t = 30),
          title = list(text = plot_label, x = 0.5)
        ) %>%
        config(displayModeBar = TRUE) %>%
        event_register("plotly_hover") %>%
        event_register("plotly_click")
    })

    ####  Hover Text Below Plot ####
    output$hoverText <- renderPrint({
      d <- event_data("plotly_hover")
      if (is.null(d)) {
        cat("Move cursor over plot to see coordinates and integer of pixel")
      } else {
        cat('x = ', d$x, "\n", sep = "")
        cat('y = ', (-1*d$y), "\n", sep = "")
      }
    })

    ####  Save Clicked X- and Y-coordinates ####
    observeEvent(event_data("plotly_click"), {
      click <- event_data("plotly_click")
      if (!is.null(click$x) && !is.null(click$y)) {
        new_df <- rbind(selected_df(), data.frame(x=click$x, y= -1*click$y))
        selected_df(new_df)
      }
    })

    #### shows list of all selected pixel ####
    output$allClicks <- renderPrint({
      df_sel <- selected_df()
      if (nrow(df_sel)==0) {
        cat("No pixels selected.")
      } else {
        print(df_sel)
      }
    })

    #### delete pixel that was selected last ####
    observeEvent(input$z_pressed, {
      df_sel <- selected_df()
      if (nrow(df_sel)>0){
        df_sel <- df_sel[-nrow(df_sel), , drop=FALSE]
        selected_df(df_sel)
      }
    })
    ####  Show Saved Value ####
    output$clickedText <- renderPrint({
      df_sel <- selected_df()
      if (nrow(df_sel)==0) {
        cat("Click to save coordinates of selected pixel.")
      } else {
        cat("pixel coord:", df_sel$x[nrow(df_sel)], df_sel$y[nrow(df_sel)])
      }
    })
    #### Close window on Enter ####
    observeEvent(input$enter_pressed, {
      final_df <- selected_df()
      indices <-  match(paste(final_df$x, final_df$y),
                        paste(df$x, df$y))
      final_df <- cbind(final_df, indices)
      assign('roi_pixel_df', final_df, envir=.GlobalEnv)
      stopApp()
    })
  }
  shinyApp(ui, server)
}

Pixel2Plot <- function(experiment, mz=NULL, i=NULL, tolerance=5L) {
  #### crreating the data frame ####
  i <- .interpret(experiment, mz, i, tolerance)
  if (is.integer(i)) {
    intensities <- spectra(experiment)[i, ]
  } else if (isFALSE(i)){
    intensities <- pData(experiment)$TIC
  }
  x_coord <- pData(experiment)$x
  y_coord <- pData(experiment)$y
  df <- data.frame(x = x_coord, y= y_coord, value= intensities)
  plot_label <- if (is.numeric(mz)) {
    paste("i = ", i, '/ m/z = ', round(mz(experiment)[i]),4)
  } else if (is.character(mz)) {
    paste("TIC")
  } else if (is.null(mz) && is.numeric(i)){
    paste("i = ", i)
  } else {
    ""
  }
  #### set up for interactive shiny window ####
  ui <- fluidPage(
    tags$head(
      tags$script(HTML("
      document.addEventListener('keydown', function(e) {
        if (e.key === 'Enter') {
          Shiny.setInputValue('enter_pressed', Math.random());
        }
        if (e.key === 'p') {
          Shiny.setInputValue('p_pressed', Math.random());
        }
      });
    "))
    ),
    style = "padding: 0; margin:0;", # uses full window without boarders
    plotlyOutput("heat_im", height = "80vh"), # plots takes up 80 percent of height.
    # Hover information
    verbatimTextOutput("hoverText"),
    # Saved click value (display only, stored globally)
    h5("Selected pixels (last click):"),
    verbatimTextOutput("clickedText"),
    plotlyOutput("create_plot")
  )
  #### creation of plot and programming cursor functions ####
  server <- function(input, output, session) {
    # Stores last clicked y-value inside Shiny session
    saved_x <- reactiveVal(NULL)
    saved_y <- reactiveVal(NULL)
    selected_df <- reactiveVal(data.frame(mz=NULL, intensity=NULL))
    selec_label <- reactiveVal(NULL)
    ####  Render Plot ####
    output$heat_im <- renderPlotly({
      # ggplot2 heat map
      p <- ggplot(df, aes(x=x_coord, y=y_coord, fill= value)) +
        geom_tile(color= 'black') +
        scale_fill_gradientn(colors = c("#4a0057", '#00598d', '#009b95', '#5dcd63','#eee42d')) +
        scale_y_reverse() +
        coord_fixed()

      # Convert to interactive plotly
      ggplotly(p, tooltip = NULL) %>%
        layout(
          hovermode = "closest",
          xaxis = list(showspikes = FALSE),
          yaxis = list(showspikes = FALSE),
          # Clean compact margins
          margin = list(l = 40, r = 10, b = 40, t = 30),
          title = list(text = plot_label, x = 0.5)
        ) %>%
        config(displayModeBar = TRUE) %>%
        event_register("plotly_hover") %>%
        event_register("plotly_click")
    })

    ####  Hover Text Below Plot ####
    output$hoverText <- renderPrint({
      d <- event_data("plotly_hover")
      if (is.null(d)) {
        cat("Move cursor over plot to see coordinates of pixel", '\n', 'Click to select pixel.')
      } else {
        cat('x = ', d$x, "\n", sep = "")
        cat('y = ', (-1*d$y), "\n", sep = "")
      }
    })
    ####  Save Clicked X- and Y-coordinates ####
    observeEvent(event_data("plotly_click"), {
      click <- event_data("plotly_click")
      if (!is.null(click)) {
        saved_x(click$x)
        saved_y(-1*click$y)
      }
    })
    ####  Show Saved Value ####
    output$clickedText <- renderPrint({
      e <- event_data('plotly_click')
      if (is.null(e)) {
        cat("Click to select pixel.")
      } else {
        cat("pixel coord:", saved_x(), saved_y(), '\n')
        cat('To plot the spectrum press "p" on the keyboard.')
      }
    })
    #### plot selected pixel ####
    observeEvent(input$p_pressed, {
      xx <- saved_x()
      yy <- saved_y()
      if (!is.null(saved_x) && !is.null(saved_y)) {
        integer <- which(df$x == xx & df$y == yy)
        plot_df <- data.frame(mzs = mz(experiment), intensities = spectra(experiment)[ ,integer])
        p_label <- paste('pixel ', integer)
        selected_df(plot_df)
        selec_label(p_label)
      }
      #### create data frame for plot ####
      output$create_plot <-renderPlotly({
        plot_df <- selected_df()
        if (!is.null(nrow(plot_df)))
        # ggplot2 line plot
        pl <- ggplot(plot_df, aes(x= mzs, y= intensities)) +
          geom_line(linewidth = 0.5) +
          theme_minimal() +
          labs(x = 'm/z', y = 'intensity', title = p_label)
        ggplotly(pl)
      })
    })
    #### Close window on Enter ####
    observeEvent(input$enter_pressed, {
      stopApp()
    })
  }
  shinyApp(ui, server)
}

PixelViewer(mse_try, 618.4507)
Pixel2Plot(mse_try, 618.4507)

mse_try <- scale(mse, mz=618.4507, method = "75%")

find_mz_matches_LEO <- function(mz_result, mz_searched, tolerance) {
  # Compute lower and upper tolerance bounds
  find_match_integer <- function(ocs, ox, mz_searched, mz_result, end_idx, start_idx) {
    results <- logical(length(ocs))
    for(n in 1:length(ocs)){
      if(ocs[n] < 1){
        results[n] <- NA
        next
      } else if (ocs[n] >= 1) {

        win_idx <- (start_idx[n]:end_idx[n])[which.min(mz_searched[start_idx[n]:end_idx[n]] - mz_result[n])]
        results[n] <- win_idx
      }
    }
    return(results)
  }

  tols <- (tolerance / 1e6) * mz_result
  lows <- mz_result - tols
  ups  <- mz_result + tols
  # Use findInterval to get indices in mz_minus
  start_idx <- findInterval(lows, mz_searched, left.open = TRUE) + 1 # returns value that is <=lows before values are higher
  end_idx   <- findInterval(ups, mz_searched)                        # returns value that is < ups before values are higher
  # Count occurrences
  ocs <- pmax(end_idx - start_idx + 1, 0)
  ocs_inverse <- c(ocs==0)
  ox <- c(ocs_inverse==FALSE)

  mse_int <- find_match_integer(ocs, ox, mz_searched, mz_result, end_idx, start_idx)

  #print(ox_int <- ocs > 1)
  #print(ocs)
  return(list(matches = ox, inverse=ocs_inverse, ocs = ocs, int = mse_int))
}

