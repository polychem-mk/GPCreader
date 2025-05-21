plotMWD <-
function(dataGPC = NULL,
                    dataMWD = NULL,
                    dataPDI = NULL,
                    add.fitted = FALSE,
                    text.size = 1){

  # (1) Check input   ----------------------

  # (1.1) dataGPC
  # This is the output of find_peaks(), containing the full GPC chromatogram and
  # the intensity after baseline correction (ybl).
  # If peaks were found, dataGPC has yi and peak_id; if peak tails are missing and
  # fitted values are added, y_fit and fitted variables are also present in dataGPC.

  # dataGPC must be a list of data frames or a data frame to be converted to
  # a list.
  if(any(class(dataGPC) %in% c("tbl_df","tbl","data.frame"))){
    dataGPC = list(as.data.frame(dataGPC))

  }else if(any(!"list"  %in% class(dataGPC) )){
    dataGPC = NULL
  }

  n_GPCpeaks = ifelse(is.null(dataGPC), 0, length(dataGPC))

  # (1.2) dataMWD
  # This is the output of the calculateMWD() function, containing the Mi and
  # dWi_over_dLogMi. These data will be used to plot Molecular Weight
  # Distribution (MWD) graphs.
  # If peak tails are missing and fitted values are added,
  # dWi_over_dLogMi_fit variable will be used instead of dWi_over_dLogMi.

  # dataMWD must be a list of data frames or a data frame to be converted to
  # a list.
  if(any(!"list"  %in% class(dataMWD) )){
    dataMWD = NULL
  }

  if("list" %in% class(dataMWD[1])){
    dataMWD = dataMWD[[1]]
  }

  n_dataMWD = ifelse(is.null(dataMWD), 0, length(dataMWD))

  # If there is no data or the length of dataGPC and dataMWD are not the same,
  # plotMWD() returns NULL.
  if(all(is.null(dataGPC), is.null(dataMWD) ) ){
    return(NULL)
  }

  if(all(!is.null(dataGPC), !is.null(dataMWD), n_GPCpeaks != n_dataMWD ) ){
    return(NULL)
  }

  # (1.3) dataPDI
  # dataPDI must be a data frame; it is an output of the  calculatePDI() function.
  # and must have FileName column.
  if(all(!class(dataPDI) %in% c("tbl_df","tbl","data.frame") )){
    dataPDI = NULL
  }

  if(!"FileName" %in% colnames(dataPDI)){
    dataPDI = NULL
  }

  n_dataPDI = ifelse(is.null(dataPDI), 0, length(unique(dataPDI$FileName)))

  # (1.4) add.fitted
  # add.fitted must be a logical and have a length of 1 or a value that can
  # be coerced to a logical, otherwise add.fitted will be automatically set
  # to the default value.
  if(any(!is.vector(add.fitted, mode = "logical"),
         length(add.fitted) != 1,
         is.na(add.fitted))){

    if(any(!is.vector(add.fitted),
           length(add.fitted) != 1,
           is.na(add.fitted))){
      add.fitted = FALSE  # set to the default value

    }else{
      add.fitted = as.logical(add.fitted)
      if(is.na(add.fitted)){
        add.fitted = FALSE  # set to the default value
      }
    }
  }

  # (2) Set variables -------------

  # (2.1) Plot colors and text sizes
  col_peaks = c("#0364cc",   "#c90259", "#357503", "#E69F00","#9223f5",
                "#0703f6", "#CC79A7","#11c902",  "#920103", "#6e03f6",
                "#02a4c1", "#d903f6","#009E73", "#d00416", "#490182")

  M_color = "#34495e"

  if(any(!is.vector(text.size, mode = "numeric"), length(text.size) != 1)){
    text.size = 1
  }else{
    if(any(text.size < 0.3, text.size > 3 )){
      text.size = 1
    }
  }

  title_cex = text.size
  peak_cex = text.size*0.95
  pdi_cex = text.size*0.95
  M_cex = text.size*0.7

  # (2.2) dataGPC and dataMWD
  # n_plot is the number of chromatograms and MWD plots to be plotted.
  if(is.null(dataMWD)){
    # If only chromatogram data are provided (no dataMWD)
    n_plot = n_GPCpeaks

  }else{
    n_plot = n_dataMWD
  }

  # (2.3) dataPDI
  # If dataPDI is not specified or FileName is missing, the graphs will not
  # display calculated data and the graph names will be integers starting from 1.

  if(all(!is.null(dataPDI), n_dataPDI == n_plot, "FileName" %in% colnames(dataPDI))){
    plot_names = unique(dataPDI$FileName)
  }else{
    # If dataPDI is NULL or file name (FileName) is missing
    plot_names = 1:n_plot
    dataPDI = NULL
  }

  # (3) By file  ----------------
  for(i in 1:n_plot) {

    # (3.1) Plot1. Chromatogram using dataGPC_i data

    # (3.1.1) Data for the chromatogram:
    dataGPC_i = dataGPC[[i]]

    # (3.1.2) The names of the graphs are file names or numbers starting from 1
    plot1_title_i <- plot_names[i]

    # (3.1.3) Check the data (dataGPC_i)
    if(!is.null(dataGPC_i)){
      if(any(!c("x", "ybl") %in% colnames(dataGPC_i))){

        dataGPC_i = NULL
      }
    }

    # (3.1.4) Plot 1. Chromatogram
    #       Make a plot with the chromatogram
    if(!is.null(dataGPC_i)){

      if(any(!c("yi", "peak_id") %in% colnames(dataGPC_i))){
        # Plot dataGPC_i if it is not null, but there no peaks
        # (yi column is missing), then no peaks are detected and only
        # the chromatogram after baseline correction can be shown.
        plot(dataGPC_i$x,
             dataGPC_i$ybl,
             type = "l",
             col = "#b3b6b7",
             lwd = 3,
             main = plot1_title_i,
             xlab = "retention time",
             ylab = "intensity",
             cex = title_cex)

        if(!is.null(dataMWD)){

          # And no MWD data
          plot.new()
        }
        next

      }else{
        # If dataGPC_i has data for plotting the chromatogram.

        # maximum intensity for the chromatogram
        max_y = max(dataGPC_i$yi, na.rm = TRUE)

        x_lim = dataGPC_i %>%
          filter(peak_id != 0)

        if(nrow(x_lim) > 0){
          x_lim = x_lim %>%
          summarise(min_x = first(x),
                    max_x = last(x))
        }else if(nrow(dataGPC_i) != 0){
          x_lim = dataGPC_i %>%
            summarise(min_x = first(x),
                      max_x = last(x))
        }else{
          x_lim = data.frame(min_x = 0,
                             max_x = 100)
        }

        # Plot1 chromatogram.
        plot(dataGPC_i$x,
             dataGPC_i$ybl,
             type = "l",
             lwd = 3,
             col = "#b3b6b7",
             xlab = "retention time",
             ylab = "intensity",
             ylim = c(-0.05*max_y, 1.15*max_y),
             xlim = c(0.5*x_lim$min_x, 1.5*x_lim$max_x),
             main = plot1_title_i,
             cex = title_cex)

        abline(0, 0, lty = 2, lwd =0.3 )
      }
    }

    # (3.2) Set the variables that will be used to add peaks and plot MWD graphs.

    # (3.2.1)  Calculated PDI, Mp, Mn and Mw from dataPDI for file i.
    # dataPDI is the output of calculatePDI() function, which must have file names,
    # otherwise, the calculation results will not be displayed on the graphs.
    if(!is.null(dataPDI)){
      dataPDI_i = dataPDI %>% filter(FileName == plot_names[i])
    }else{
      dataPDI_i = NULL
    }

    # (3.2.2) Molecular weight distribution data for file i.
    dataMWD_i = dataMWD[[i]]

    if(!is.null(dataMWD_i)){

      if( any(!c("x", "yi", "peak_id", "Mi","dWi_over_dLogMi") %in%
          colnames(dataMWD_i))){
        dataMWD_i = NULL

      }else if(all(add.fitted, c("y_fit", "dWi_over_dLogMi_fit") %in% colnames(dataMWD_i))){

        dataMWD_i = dataMWD_i %>%
          mutate(dWi_over_dLogMi = ifelse(!is.na(dWi_over_dLogMi_fit),
                                          dWi_over_dLogMi_fit,
                                          dWi_over_dLogMi)) %>%
          select(x, yi, y_fit, peak_id, peak, Mi, dWi_over_dLogMi)

      }else{
        dataMWD_i = dataMWD_i %>%
            select(x, yi, peak_id, peak, Mi, dWi_over_dLogMi)
      }
    }

    # If there is not enough data to add peaks and molecular weight distribution plots:
    if(any(is.null(dataPDI_i), is.null(dataMWD_i))){
      next
    }

    # (3.2.3) Peak names
    peaks = unique(dataPDI_i$peak)

    # (3.3) Plot 1. Add peaks
    # If there are dataMWD and dataPDI inputs and a chromatogram plot

    if(!is.null(dataGPC_i)){

      for(p in peaks){
        # Set peak boundaries for peak p

        ## PDI is calculated within time1 and time2
        time1_p = dataPDI_i$time1[dataPDI_i$peak == p]
        time2_p = dataPDI_i$time2[dataPDI_i$peak == p]

        # Intensity and time data for the peak p from dataMWD_i
        dataMWD_p = dataMWD_i %>%
          filter(peak == p) %>%
          filter(x >= time1_p & x <= time2_p)

        ## If peak tails have been added, the peak boundaries from
        # the original chromatogram will differ from those found after fitting
        # and adding tails.
        peak_ind = which(dataMWD_p$peak_id != 0)

        # Add peak p to the chromatogram
        lines(dataMWD_p$x[peak_ind],
              dataMWD_p$yi[peak_ind],
              col = col_peaks[p],
              lwd = 3)

        # Add peak name
        text(x = dataPDI_i$time_max[dataPDI_i$peak == p],
             y = max(dataMWD_p$yi[peak_ind]) + max_y*0.1 ,
             labels = p, col = col_peaks[p], cex = peak_cex)

        # Add peak boundaries
        points(x = c(time1_p, time2_p),
               y =c(0, 0),
               col = col_peaks[p],
               pch = "|", cex = 2 )

        if(all(add.fitted, "y_fit" %in% colnames(dataMWD_i)) ){

          # Add peak with fitted tails
          if(sum(!is.na(dataMWD_p$y_fit)) > 10 ){
            lines(dataMWD_p$x,
                  dataMWD_p$y_fit,
                  col = "black",
                  # lty = 3,
                  lwd = 1.5)
          }
        }
      }
    }

    # (3.4)  Plot2. Molecular weight distribution

    if(any(sum(is.na(dataMWD_i$Mi)) > 0,
           sum(dataMWD_i$Mi  <= 0) > 0)){
      plot.new()
      next
    }

   # Maximum y-axis limit for the plot
    max_dWi_over_dLogMi = dataMWD_i %>%
      filter(peak_id != 0) %>%
      summarise(max(dWi_over_dLogMi, na.rm = TRUE)) %>%
      pull()

    plot(log10(dataMWD_i$Mi),
         dataMWD_i$dWi_over_dLogMi,
         main = "Molecular weight distribution",
         xlab = "log(Mi)",
         ylab = "dWi/dlog(Mi)",
         type = "l", lwd = 3, col = "white", cex = title_cex,
         ylim = c(-0.25*max_dWi_over_dLogMi, max_dWi_over_dLogMi*1.15))

    abline(0, 0, lty = 2, lwd =0.3 )

    # (3.5) Add molecular weight distribution by peak

    for(p in peaks){

      # Set peak boundaries for peak p
      time1_p = dataPDI_i$time1[dataPDI_i$peak == p]
      time2_p = dataPDI_i$time2[dataPDI_i$peak == p]

      # Molecular weight distribution data for the peak p
      dataMWD_p = dataMWD_i %>%
        filter(peak == p) %>%
        filter(x >= time1_p & x <= time2_p)

      # Add dataMWD_p to the plot2
      lines(log10(dataMWD_p$Mi),
            dataMWD_p$dWi_over_dLogMi,
            col = col_peaks[p],
            lwd = 3)

      # Add PDI to the plot:
      PDI_i = dataPDI_i$PDI[dataPDI_i$peak == p]

      text(x = median(log10(dataMWD_p$Mi), na.rm = TRUE),
           y = max(dataMWD_p$dWi_over_dLogMi) + max_dWi_over_dLogMi*0.05 ,
           labels = paste("peak ", p,
                          ", PDI = ", PDI_i,
                          sep = ''),
           col = col_peaks[p], cex = pdi_cex)

      # Add Mp to the plot
      x_Mp = log10(dataPDI_i$Mp[dataPDI_i$peak == p])
      points(x =  x_Mp,
             y = 0,
             col = M_color,  pch = "|", cex = 0.8)

      text(x = x_Mp,
           y = -0.25*max_dWi_over_dLogMi,
           labels = paste("Mp =", dataPDI_i$Mp[dataPDI_i$peak == p]),
           col = M_color, cex = M_cex)

      # Add Mn to the plot
      x_Mn = log10(dataPDI_i$Mn[dataPDI_i$peak == p])
      points(x =  x_Mn,
             y = 0,
             col = M_color, pch = "|", cex = 0.8)

      text(x = x_Mn,
           y = -0.08*max_dWi_over_dLogMi,
           labels = paste("Mn =", dataPDI_i$Mn[dataPDI_i$peak == p]),
           col = M_color, cex = M_cex)

      # Add Mw to the plot
      x_Mw = log10(dataPDI_i$Mw[dataPDI_i$peak == p])
      points(x =  x_Mw,
             y = 0,
             col = M_color, pch = "|", cex = 0.8)

      text(x = x_Mw,
           y = -0.16*max_dWi_over_dLogMi,
           labels = paste("Mw =", dataPDI_i$Mw[dataPDI_i$peak == p]),
           col = M_color, cex = M_cex)
    }
  }
}
