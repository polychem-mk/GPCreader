# GPC reader, server

server <- function(input, output, session) {

  #  1.Calibration files                  -------------------------------------
  ## 1.1.  Reactive values                    ---------------------------------

  ### Table 1. Calibration standards             ------------------------------

  # 1.1.1 Data frame with information about calibration files: file name (FileName),
  #       sample name (SampleName), number of numeric columns (DataColumns)
  #       first row with numeric data (DataStart) and file path (FilePath).
  calib_data_info <- reactive( {
    file_info(input$calib_file$datapath, input$calib_file$name)
  })

  # 1.1.2 Calibration data files.
  # 'calib_data' is a list of data frames with GPC chromatograms for calibration
  # standards.

  calib_data <- reactive( {

    # System peak for calibration standards. If added, find_peaks() will trim
    # the chromatogram to this value
    s_peak <- input$syst_peak_calib

    withProgress(
      message = 'Calculation in progress',
      {
        # 'read_calib_files' is a list of data frames with the original data, where
        # x is the time (or volume) and y is the detector signal intensity
        read_calib_files = readGPC(calib_data_info())

        # 'calib_peaks' is a list of data frames with data after intensity correction
        # and peak detection. If the chromatogram is noisy, smoothing is also applied
        # to the intensity. shoulder.sens = 1 is setting for calibration files; if
        # the chromatogram is a mixture of standards, then merged peaks (shoulder peaks)
        # will be detected if their intensity exceeds 25% of the maximum intensity
        calib_peaks = find_peaks(read_calib_files,
                                 syst.peak = s_peak,
                                 shoulder.sens = 1 )

        # If the input 'first_peak' is checked and more than one file is loaded,
        # then for files with calibration standards, only the first peak from each
        # chromatogram is detected.
        if(input$first_peak & length(calib_peaks) > 1){
          filter_peaks(calib_peaks,
                       first.peak = TRUE,
                       tails = FALSE)
        }else{

          # Otherwise, all peaks are detected, including shoulder peaks.
          filter_peaks(calib_peaks,
                       first.peak = FALSE,
                       tails = FALSE)
        }
      })
  })

  # 1.1.3 Calibration files information is in the data frame with columns:
  #       fileName (file name), Mp (peak molecular weight), RT (retention time).
  #       'calib_DT' is displayed as Table 1 in the Calibration tab.
  calib_DT <- reactiveValues(data =  {
    data.frame(fileName = character(0), Mp = numeric(0), RT = numeric(0))
  })

  # 'calib_Mp_RT' is created from the input data and assigned to 'calib_DT'
  calib_Mp_RT <- reactive( {
    # calibration_info() function finds RT (retention time) from calib_data()
    # and Mp (peak molecular weight) values from calib_data_info().
    calibration_info(calib_data_info(), calib_data())
  })

  ### Plot 1. GPC chromatogram of standards      ------------------------------

  # Plot 1 is displayed in the calibration tab and shows all the chromatograms for
  # the calibration standards. It is built from calib_data(). Plot 1 can be enlarged
  # by brushing over the area and reduced by clicking the "Zoom Out" button.

  # (1.1.4) Axis limits for the Plot 1.
  xy_lims_calib_plot <- reactiveValues(data = {
    data.frame(x_min = numeric(0),
               x_max = numeric(0),
               y_min = numeric(0),
               y_max = numeric(0) )
  })

  # 'xy_lims_from_data' is x and y axis limits for Plot 1 from data
  xy_lims_from_data <- reactive( {

    if(length(input$calib_file$datapath) == 1){

      # If all calibration standards are in one file, then the x and y limits
      # for plot 1 are the minimum and maximum values from that file.
      calibration_data = calib_data()[[1]]

      if(all(c("x", "yi") %in% colnames(calibration_data))){
        calibration_data = calibration_data %>%
          summarise(x_min = min(x, na.rm = TRUE),
                    x_max = max(x, na.rm = TRUE),
                    y_min = min(yi, na.rm = TRUE),
                    y_max = max(yi, na.rm = TRUE)*1.1 )

      }else if(all(c("x", "ybl") %in% colnames(calibration_data))){
        calibration_data = calibration_data %>%
          summarise(x_min = min(x, na.rm = TRUE),
                    x_max = max(x, na.rm = TRUE),
                    y_min = min(ybl, na.rm = TRUE),
                    y_max = max(ybl, na.rm = TRUE)*1.1 )
      }else{
        calibration_data = data.frame(x_min = 0,
                                      x_max = 100,
                                      y_min = 0,
                                      y_max = 100)
      }
      calibration_data

    }else{

      # If each calibration chromatogram is in its own file, then the x and y
      # limits for graph 1 are the minimum and maximum values across all files.
      map(calib_data(), function(dat){

        if(all(c("x", "yi") %in% colnames(dat))){
          calibration_data = dat %>%
            summarise(x_min = min(x, na.rm = TRUE),
                      x_max = max(x, na.rm = TRUE),
                      y_min = min(yi, na.rm = TRUE),
                      y_max = max(yi, na.rm = TRUE)*1.1 )

        }else if(all(c("x", "ybl") %in% colnames(dat))){
          calibration_data = dat %>%
            summarise(x_min = min(x, na.rm = TRUE),
                      x_max = max(x, na.rm = TRUE),
                      y_min = min(ybl, na.rm = TRUE),
                      y_max = max(ybl, na.rm = TRUE)*1.1 )
        }else{
          calibration_data = data.frame(x_min = 0,
                                        x_max = 10,
                                        y_min = 0,
                                        y_max = 5 )
        }
        calibration_data

      }) %>%
        bind_rows() %>%
        summarise(x_min = min(x_min, na.rm = TRUE),
                  x_max = max(x_max, na.rm = TRUE),
                  y_min = min(y_min, na.rm = TRUE),
                  y_max = max(y_max, na.rm = TRUE)*1.1 )
    }
  })

  # 1.1.5 Make Plot 1.
  calibration_files_plot <- reactive({

    par(mar = c(5, 4, 2, 2) + 0.1, mgp = c(2, 1, 0))

    if(length(input$calib_file$name) == 1){

      # If all calibration standards are in one file with the one chromatogram.
      calib_df = calib_data()[[1]]

      if("yi" %in% colnames(calib_df)){

        calib_chromatogram = calib_df %>%
          distinct(x, .keep_all = TRUE) %>%
          select(x, yi)
        plot(calib_chromatogram$x, calib_chromatogram$yi,
             type = "l", lwd = 2,  frame.plot = FALSE,
             xlab = "retention time",  ylab = "intensity",
             col = "#b3b6b7",
             ylim = c(xy_lims_calib_plot$data$y_min, xy_lims_calib_plot$data$y_max)  ,
             xlim = c(xy_lims_calib_plot$data$x_min, xy_lims_calib_plot$data$x_max)   )

        indx = which(calib_df$peak_id != 0)

        lines(calib_df$x[indx], calib_df$yi[indx],  lwd = 2)

      }else{
        plot(1,  frame.plot = FALSE, xlab = "retention time", ylab = "intensity", col = "white"  )
      }

    }else{

      # If calibration samples are in multiple files.
      plot(1, frame.plot = FALSE, xlab = "retention time", ylab = "intensity",
           col = "white",
           ylim = c(xy_lims_calib_plot$data$y_min, xy_lims_calib_plot$data$y_max)  ,
           xlim = c(xy_lims_calib_plot$data$x_min, xy_lims_calib_plot$data$x_max)     )

      for (calib_file_i in calib_data()) {

        if("yi" %in% colnames(calib_file_i)){
          calib_file_i = calib_file_i %>%
            distinct(x, .keep_all = TRUE)

          indx = which(calib_file_i$peak_id != 0)

          lines(calib_file_i$x[indx], calib_file_i$yi[indx],  lwd = 2)

        }else{
          next
        }
      }
    }
    abline(v = xy_from_calib_plot$data$x)
    abline(h = xy_from_calib_plot$data$y)
  })

  # 1.1.6 The x and y values obtained when clicking on Plot 1 appear in the text
  #       box below the graph.
  xy_from_calib_plot <- reactiveValues(data = {
    data.frame(x = 0, y = 0)
  })

  ### Table 2. Calibration coefficients          ------------------------------

  # The calibration coefficients are summarized in Table 2 ('calib_coeffs') on
  # the calibration tab. 'calib_coeffs' can have values calculated from data or
  # input from the sidebar. It also depends on the input 'calib_fit_order'.

  # 1.1.7 Calibration coefficients that are set at the beginning

  calib_coeffs <- reactiveValues(data = {
    data.frame(intercept = numeric(0),
               coef1 = numeric(0),
               coef2 = numeric(0),
               coef3 = numeric(0) )
  })

  # calibration_coeffs() function takes the output of calibration_info()
  # function (calib_Mp_RT() reactive) and returns a two-row data frame with
  # the calibration coefficients and the order for the calibration curve.

  # 1.1.7.1 Calibration coefficients calculated from data

  calib_coeffs_data <- reactive( {

    if(length(input$calib_file$datapath) != 0){

      suppressWarnings({
        coeffs =  calibration_coeffs(calib_Mp_RT())
      })

      if(any(nrow(coeffs) == 0, is.null(coeffs))){
        coeffs = data.frame(intercept = numeric(0),
                            coef1 = numeric(0),
                            coef2 = numeric(0),
                            coef3 = numeric(0) )
      }else{
        coeffs = coeffs %>%
          filter( fit_order == as.numeric(input$calib_fit_order ) ) %>%
          select(intercept, coef1, coef2, coef3)
      }
    }
  })

  # 1.1.7.2 Calibration coefficients calculated from edited data

  calib_coeffs_update <- reactive( {
    if(length(input$calib_file$datapath) != 0 ){

      suppressWarnings({
        coeffs =  calibration_coeffs(calib_DT$data)
      })

      if(any(nrow(coeffs) == 0, is.null(coeffs))){
        coeffs = data.frame(intercept = numeric(0),
                            coef1 = numeric(0),
                            coef2 = numeric(0),
                            coef3 = numeric(0) )
      }else{
        coeffs = coeffs %>%
          filter( fit_order == as.numeric(input$calib_fit_order ) ) %>%
          select(intercept, coef1, coef2, coef3)
      }
    }
  })

  # 1.1.7.3 Calibration coefficients from sidebar

  calib_coeffs_sidebar <- reactive( {
    intercept <- input$intercept

    if( input$calib_fit_order == "1"){
      coef1 <-  input$coef1
      coef2 <-  0
      coef3 <-  0
    }else{
      coef1 <-  input$coef1
      coef2 <-  input$coef2
      coef3 <-  input$coef3
    }

    data.frame( intercept = intercept,
                coef1 = coef1,
                coef2 = coef2,
                coef3  = coef3)
  })

  ### Plot 2. GPC calibration curve              ------------------------------

  # 1.1.9 Create a calibration curve, which will be displayed as Plot2 in
  #       the Calibration tab if 'use_input_coefs' input is FALSE.
  calibration <- reactive({

    # Actual data are shown as points, where each point represents the logarithm
    # of peak molecular weight (log Mp) versus retention time (RT) for each
    # calibration standard
    actual = calib_DT$data %>%
      drop_na() %>% filter(RT != 0) %>%  filter( Mp != 0 )

    if(nrow(actual) >= 2){
      coefs <- calib_coeffs$data

      if(all(nrow(coefs) != 0,
             sum(c(coefs$intercept, coefs$coef1), na.rm = TRUE) != 0 ,
             sum(is.na(coefs)) == 0)){

        # The purple line shows the fit: linear or polynomial, depending on
        # 'calib_fit_order' input.
        fitted_vals = data.frame(rt = seq(min(actual$RT),
                                          max(actual$RT),
                                          length.out = 20) ) %>%
          mutate(log_mp = coefs$coef3*(rt^3) +
                   coefs$coef2*(rt^2) +
                   coefs$coef1*rt + coefs$intercept)
        par(mar = c(5, 4, 2, 2) + 0.1, mgp = c(2, 1, 0))
        plot( actual$RT, log10(actual$Mp ),
              xlab = "retention time",  ylab = "log(Mp)"  )
        lines( fitted_vals, col = "purple" , lwd = 2)

      }
    }
  })

  ### Help text on the side panel                ----------------------------

  # Text is displayed if the entered coefficients are missing or invalid.
  # For a linear calibration curve, the intercept is expected to be greater
  # than 0 and the slope (coef1) to be negative, while for a polynomial curve,
  # all 4 coefficients are not equal to NA.
  calib_help_text1 <- reactive({

    use_coefs = input$use_input_coefs
    intercept_is_0 = calib_coeffs$data$intercept == 0
    intercept_is_na = is.na(input$intercept)
    coef1_is_na = is.na(input$coef1)

    if(all(use_coefs, input$calib_fit_order == "1",
           any(intercept_is_0, intercept_is_na, coef1_is_na )) ){
      "One or more calibration coefficients are missing or set to an invalid
      value. The Intercept must be greater than 0."

    }else if(all(use_coefs, input$calib_fit_order == "3",
                 any(intercept_is_0,
                     intercept_is_na,
                     coef1_is_na,
                     is.na(input$coef2),
                     is.na(input$coef3) ) )){
      "One or more calibration coefficients are missing or set to an invalid value. "

    }else{
      ""
    }
  })

  ### Help text on the main panel                ------------------------------

  # Help text is displayed if values in Table 1 are missing or equal
  # to 0, or if there are less than 6 standard points left in Table 1
  # after all corrections
  calib_help_text2 <- reactiveValues(data = {
      'Select calibration files or enter "intercept" and "slope". '
  })

  calib_help_text2_update <- reactive( {

    if( 0 %in%  calib_DT$data$Mp  ){
       ' Missing Mp values for calibration samples (zeros or empty fields in Table 1)
        will be excluded from the calculations. To change Mp or retention time (RT),
      double-click the corresponding cell in the Table 1 and enter a new value.'

    }else if( any(c(0, NA) %in%  calib_DT$data$RT)){
      'Missing retention time (RT) values or values equal to zero in Table 1 will
       be excluded from the calculations. To enter this values manually double-click
       the corresponding cell in the Table 1 and enter a new value.'
    }else if(nrow(calib_DT$data) < 1){
      'Select calibration files or enter "intercept" and "slope". '
    }else if(nrow(calib_DT$data) < 6){
      'To construct a calibration curve, at least 6 standard points are required.
      It is recommended to use 10 to 12 calibration standards. '
    }else{
      ""
    }
  })

  ## 1.2 Events                               ---------------------------------

  #                      Events for the calibration files.

  ### calib_file, first_peak, syst_peak_calib    ------------------------------

  events_1 <- reactive({
    list(
      input$calib_file,      # loading calibration files
      input$first_peak,      # changing the checkbox 'one sample per file'
      input$syst_peak_calib  # adding a system peak for the calibration chromatogram(s)
      )
  })

  # events_1  affects reactive: calib_data_info() and calib_data();
  #           these two are inputs to the reactive calib_Mp_RT()

  observeEvent(events_1(), {

    validate( need(input$calib_file != "" ,    ' ') )

    # Get peak retention times from calibration files and Mp values from sample
    # names. This will create Table 1, which is rendered from 'calib_DT'
    calib_DT$data = calib_Mp_RT()

    # Set axis limits for the Plot 1.
    xy_lims_calib_plot$data <- xy_lims_from_data()

    # calibration coefficients:
    if(input$use_input_coefs){
      # Get calibration coefficients from sidebar,
      # if input 'use_input_coefs' is TRUE
      calib_coeffs$data <- calib_coeffs_sidebar()

    }else{

      # Find calibration coefficients from data
      # if calibration standard files are loaded
      if(!is.null(calib_coeffs_data())){
        calib_coeffs$data <- calib_coeffs_data()
      }
    }

    # Update a help message in the Calibration tab
    calib_help_text2$data <- calib_help_text2_update()

  })

  ### calib_samples_table_cell_edit              ------------------------------

  # Event2: edit and update Table 1

  # This adds new values to Table 1 ('calib_DT') and causes the calibration
  # coefficients to be recalculated and the 'calib_coeffs' table to be updated.
  observeEvent(input$calib_samples_table_cell_edit, {

    # Edit cells of the calib_DT table
    info = input$calib_samples_table_cell_edit
    i = as.numeric(info$row)
    j = as.numeric(info$col)
    val = as.numeric(info$value)

    if(is.na(val)){     # set to 0 if NA
      val <- 0
    }
    if(val < 0){        # convert to positive if negative
      val <- val * -1
    }
    calib_DT$data[i, j] <- val

    # update 'calib_coeffs'
    calib_coeffs$data <- calib_coeffs_update()

    # update help text
    calib_help_text2$data <- calib_help_text2_update()

  }   )

  ### plot_click                                 ------------------------------

  # Event3:  click on the Plot 1 and show x and y values in the text box below

  observeEvent(input$plot_click, {
    validate( need(input$calib_file != "" ,    ' ') )

    xy_from_calib_plot$data$x <- input$plot_click$x
    xy_from_calib_plot$data$y <- input$plot_click$y
  })

  ### plot_brush                                 ------------------------------

  # Event4: input$plot_brush; brush the plotting area, and
  # zoom in Plot 1 by setting new axis limits xy_lims_calib_plot
  observeEvent(input$plot_brush, {
    validate( need(input$calib_file != "" ,    ' ') )

    xy_lims_calib_plot$data <-
      data.frame(x_min = input$plot_brush$xmin,
                 x_max = input$plot_brush$xmax,
                 y_min = input$plot_brush$ymin,
                 y_max = input$plot_brush$ymax   )
  })

  ### add_rt_from_click                          ------------------------------

  # Event5: click on 'add RT' and add new row to the Table 1 ('calib_DT')
  # with RT = x and Mp = 0 (then the missing Mp need to be added manually)
  observeEvent(input$add_rt_from_click, {

    calib_DT$data <- bind_rows(calib_DT$data ,
                               data.frame(FileName = "",
                                          Mp = 0,
                                          RT = round(as.numeric(xy_from_calib_plot$data$x), 4)))

    calib_help_text2$data <- calib_help_text2_update()
  }   )

  ### zoom_out                                   ------------------------------

  # Event6: Click the 'zoom_out' button; zoom out from plot 1, setting the axis
  # limits to the initial values stored in xy_lims_from_data()
  observeEvent(input$zoom_out, {
    validate( need(input$calib_file != "" ,    ' ') )

    plot_zoom_out <- xy_lims_from_data() %>%
      mutate(x_min = x_min*0.5,
             x_max = x_max*1.5,
             y_min = y_min - y_min*0.3,
             y_max = y_max*1.3 )

    xy_lims_calib_plot$data <- plot_zoom_out
  })

  ### calib_coeffs                               ------------------------------

  # Recalculate calibration coefficients stored in 'calib_coeffs' or switch to
  # manually entered calibration coefficients, even if calibration files are loaded

  events_7 <- reactive({
    list(
      # If 'use_input_coefs' is checked, then manually entered calibration
      # coefficients will be used.
      input$use_input_coefs,

      # Selecting between polynomial and linear calibration curve
      input$calib_fit_order,

      # If 'first_peak' is checked, then for files with calibration standards,
      # only the first peak from each chromatogram is detected
      input$first_peak,

      # Adding a system peak value to trim the calibration chromatogram(s)
      input$syst_peak_calib
      )
  })

  observeEvent(events_7(), {

    validate( need(input$calib_file != "" ,    ' ') )

    if(input$use_input_coefs){
      # using manually entered coefficients from the sidebar
      calib_coeffs$data <- calib_coeffs_sidebar()

    }else{
      # using coefficients calculated from data
      calib_coeffs$data <- calib_coeffs_update()
    }

    # using coefficients calculated from data
       calib_help_text2$data <- calib_help_text2_update()
  })

  ### use_input_coefs                            ------------------------------

  # Update 'calib_coeffs' or switch to entering calibration coefficients manually,
  # regardless of whether calibration files are loaded or not.

  events_8 <- reactive({
    list(
      # 'use_input_coefs'  toggles between manually entered calibration
      #  coefficients and calculating these coefficients from the data.
      input$use_input_coefs,

      # manually entered calibration coefficients
      input$intercept,
      input$coef1,
      input$coef2,
      input$coef3,

      # 'calib_fit_order'  toggles between linear calibration and polynomial
      input$calib_fit_order)
  })

  # Event8 result in setting new calibration coefficients
  observeEvent(events_8(), {

    if(input$use_input_coefs){
      # using manually entered coefficients from the sidebar
      calib_coeffs$data <- calib_coeffs_sidebar()

    }else if(length(input$calib_file$name) != 0){
      # If there are calibration files
      calib_coeffs$data <- calib_coeffs_update()

    }else{
      # If there are no calibration files and 'use_input_coefs' is FALSE,
      # no calculations can be performed.
      calib_coeffs_NULL <- data.frame(intercept = numeric(0),
                                     coef1 = numeric(0),
                                     coef2 = numeric(0),
                                     coef3 = numeric(0) )

      calib_coeffs$data <- calib_coeffs_NULL
    }
  }, ignoreInit = TRUE)

  #  2. GPC files                         -------------------------------------

   ## 2.1.  Reactive values                   ---------------------------------

  ### GPC data                                   ------------------------------

  # 'gpc_data' is a list of data frames with data after intensity correction and
  # peak detection. If the chromatogram is noisy, smoothing is also applied to the
  # intensity. shoulder.sens = 2 is setting for GPC samples: merged peaks (shoulder
  # peaks) are detected if their intensity exceeds 5% of the maximum intensity.

  # If the input 'syst_peak' is added, find_peaks() will trim the chromatogram(s)
  # to this value.

  gpc_data  <- reactive( {
    s_peak = input$syst_peak

    withProgress(
      message = 'Calculation in progress',
      {
        file_info_GPC_files = file_info(path = input$gpc_file$datapath,
                                        name = input$gpc_file$name)
        read_GPC_files = readGPC(file_info_GPC_files)

        find_peaks(dataGPC = read_GPC_files,
                   syst.peak = s_peak,
                   shoulder.sens = 2)
      })
  })

  ### GPC peaks                                  ------------------------------

  # 'gpc_filter' is a list of data frames with part of the chromatogram that has
  # only detected peaks (and a small part of the chromatogram before and after
  # each peak); it is used for calculations and to create graphs on the Data tab.

  # If the peak tails are absent from the merger with other peaks or has a shoulder,
  # and the input 'fit_missing_tail' is checked, then filter_peaks() tries to fit
  # the bell curve and add missed values.

  gpc_filter <- reactive( {

    filter_peaks(gpc_data(),
                 first.peak = FALSE,
                 above = 0,
                 tails = 0.02)
  })

  ### MWD data                                   ------------------------------

  # Molecular weight distribution (MWD).

  # 'gpc_calc' is a list of length 2:
  # the first is a list of data frames with calculated data necessary to
  # make MWD graph and to find PDI and average molecular weight (Mw and Mn);
  # the second is a vector with calibration coefficients.
  gpc_calc <- reactive( {

    suppressWarnings({
      calculateMWD( gpc_filter(),
                    intercept = calib_coeffs$data$intercept,
                    coef1 = calib_coeffs$data$coef1,
                    coef2 = calib_coeffs$data$coef2,
                    coef3 = calib_coeffs$data$coef3 )
    })
  })

  ### Peak start/end                             ------------------------------

  # The area included in the each peak is detected by functions, but the boundaries
  # can be moved using the input values ('t1', 't2', 'use_t1' and 'use_t2')
  peak_boundaries_shift <- reactiveValues(data =  {
    data.frame(start = 0,
                 end = 0)
  })

  # 't1' and 't2' are numeric input values (percentage of the peak width),
  # and the use_t1 and use_t2 checkbox are used to apply changes.
  peak_boundaries_shift_input <- reactive( {

    req(input$use_t1 == TRUE | input$use_t2 == TRUE)

    data.frame(start = input$t1,
               end = input$t2)
  })

  ### PDI results                                ------------------------------

  # Data frame with the results for each GPC file is shown on the Data tab as
  # 'Table 3. Results' (from 'results' reactive)

  # This table has the following columns:
  # FileName     GPC file name;
  # peak         peak number;
  # Mp           peak molecular weight;
  # Mn and Mw    average molecular weights;
  # PDI          polydispersity index;
  # time1 and time2 are the start and end times of the peak, respectively;
  # time_max     the retention time at the peak's maximum;
  # Af           asymmetry factor;
  # Tf           tailing factor;
  # intercept, coef1, coef2 and are calibration coefficients from input.

  results <- reactive( {

    # The values 'start'/'end' are used to change the boundaries of the peak
    # before calculations, if not set for 0
    start = peak_boundaries_shift$data$start
    end   = peak_boundaries_shift$data$end

    # 'calib' indicates whether it is necessary to add calibration coefficients
    # to the results table
    calib = input$add_calib_info

    # 'fit' indicates whether the values for missing peak tails were added and
    # the calculations for peaks will include these values
    fit = input$fit_missing_tail

    suppressWarnings({
      calculatePDI(dataMWD = gpc_calc(),
                   name = input$gpc_file$name,
                   shift.start = start,
                   shift.end = end,
                   use.fitted = fit,
                   add.calib.info = calib)
    })
  })

  ### Help text                                  ------------------------------

  # Help text is displayed if calibration coefficients are not specified or
  # are incorrect
  help_text3 <- reactive( {

    if("PDI" %in% colnames(results())){
      if(nrow(results()) > 0){
        PDI_na = sum(!is.na(results()$PDI))
        PDI_max = ifelse(PDI_na > 0,
                         max(results()$PDI[which(results()$peak == 1)], na.rm = TRUE),
                         10)

        if(any(sum(is.na(results()$PDI)) == nrow(results()), PDI_max > 10)){
          "Check calibration coefficients... "
        }else{
          ""
        }

      }else{
        "Calibration coefficients are required to calculate PDI"
      }

    }else{
      "Calibration coefficients are required to calculate PDI"
    }
  })

  ## 2.2.  Events                            ----------------------------------

  # Changing peak boundaries using inputs:
  # the 'use_t1' checkbox and the numeric input 't1' are intended to shift
  # the start of the peak by the value t1,
  # the 'use_t2' checkbox and the numeric input 't2' are intended to shift
  # the end of the peak by the value t2

  events_9 <- reactive({
    list(input$use_t1,
         input$use_t2,
         input$t1,
         input$t2)
  })

  observeEvent(events_9(), {

    validate( need(input$gpc_file != "" ,    ' ') )

    # The change in peak start/end values is stored in "peak_boundaries_shift"
    # which is updated if events_9() occurs
    if(input$use_t1){
      peak_boundaries_shift$data$start = peak_boundaries_shift_input()$start
    }else{
      peak_boundaries_shift$data$start = 0
    }

    if(input$use_t2){
      peak_boundaries_shift$data$end = peak_boundaries_shift_input()$end
    }else{
      peak_boundaries_shift$data$end = 0
    }
  })

  # 3.  Events for the 'Help' buttons       -----------------------------------

  # Each tab has 'Help' buttons that display a module with help text (help_text4
  # and help_text5).

  observeEvent(input$help1, {
    showModal(  modalDialog(help_text4,
                            size = "l",
                            easyClose = TRUE )    )
  })
  observeEvent(input$help2, {
    showModal( modalDialog(help_text5,
                           size = "l" ,
                           easyClose = TRUE)    )
  })

  # 4. Render                             -------------------------------------

  ## 4.1.  Sidebar                           ----------------------------------

  #  If incorrect/missing values are entered in the fields,
  # 'intercept', 'coef1', 'coef2' or 'coef3', the help text 1 is displayed.
  output$calib_help1 <- renderText({
    validate( need( input$use_input_coefs == TRUE,  '') )
    calib_help_text1()
  })

  ## 4.2  Calibration TAB                    ----------------------------------

  # 4.2.1.  Table 1  shows the names of the calibration files,
  #                  Mp values and peak retention times for each standard.
  output$calib_samples_table <- renderDT({

    validate( need(input$calib_file != "" , '') )

    DT::datatable(calib_DT$data,
                  class = "row-border",
                  selection = "none",
                  options = list(searching = FALSE,
                                 filter = "none",
                                 lengthChange = FALSE),
                  editable = list(target = "cell",
                                  numeric = "all",
                                  disable = list(columns = 1)))
  } )

  # 4.2.2. Plot 1. GPC chromatogram of standards: intensity  vs time (or volume)
  output$calibration_files <- renderPlot( {
    validate( need(input$calib_file != ""  ,  '') )
    calibration_files_plot()
  })

  # 4.2.3. Text box, shows x and y values obtained by clicking on the Plot 1 area.
  output$xy_calib <- renderText({
    validate( need(input$calib_file != "" & input$use_input_coefs == FALSE,    ' ') )
    paste0(" x = ", round(xy_from_calib_plot$data$x, 4),
           ",  y = ", round(xy_from_calib_plot$data$y,4)  )
  })

  # 4.2.4. If Table 1 contains zeros or missing values, help text 2 is displayed.
  output$calib_help2 <- renderText({

    calib_help_text2$data
  })

  # 4.2.5. Table 2 shows the calibration plot coefficients calculated from the
  # calibration files (in 'calib_coeffs') or entered in the sidebar.
  output$calib_coef <- renderTable({

    withProgress(
      message = 'Calculation in progress',
      {
        if(nrow(calib_coeffs$data) != 0){
          data.frame( coefficient = colnames(calib_coeffs$data ),
                      value = unlist(calib_coeffs$data ) )
        }
      })
  }, digits = 6)

  # 4.2.6. Plot 2. Calibration curve. The logarithm of the molecular weight of
  # polymer standards (column Mp in the Table 1), vs the retention time (column
  # RT in the Table 1, time at the calibration standard peak maximum).
  output$calibration_plot <- renderPlot( {
    validate( need( input$calib_file != "" & input$use_input_coefs == FALSE, '') )
    calibration()
  })

  ## 4.3. Data TAB                           ----------------------------------

  # 4.3.1  Table 3 shows results for each GPC file.
  output$data_table <- DT::renderDataTable( {

    validate( need( length(input$gpc_file) != 0, '') )

    DT::datatable(results(),
                  class = "row-border",
                  options = list(searching = FALSE,
                                 filter = "none",
                                 lengthChange = FALSE))
  })

  # 4.3.2  Handle saving results
  output$save_results <- downloadHandler(
    filename = function() {
      paste("results", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(results() , file)
    }
  )

  # 4.3.3 Help text 3 is displayed if the PDI values in Table 3 for peak 1 are
  # greater than 10 or all are missing, or if calibration coefficients are not
  # provided.
  output$help3 <- renderText({
    validate( need(  input$gpc_file != "",  '') )
    help_text3()
  })

  ### Plot 3                                    -------------------------------

  # Plots displayed in the Data tab.
  # The plotMWD() function creates two plots for each GPC file: a chromatogram
  # and a molecular weight distribution (MWD)

   output$plot1 <- renderPlot( {
    validate( need( input$gpc_file != "", 'Select GPC files.') )

    # 'fit' indicates whether the values for missing peak tails were added.
    # If TRUE, then the fitted peaks are also plotted.
    fit = input$fit_missing_tail

    par(mfrow = c(length(input$gpc_file$name) ,2),
        cex = 0.9, mgp = c(2, 1, 0),
        mar = c(4, 3, 2, 1) + 0.1)

    withProgress(
      message = 'Making plots',

      {
        plotMWD(dataGPC = gpc_data(),
                dataMWD = gpc_calc(),
                dataPDI = results(),
                add.fitted = fit,
                text.size = 1)
      })
  } ,
  height = function(){50 + round(length(input$gpc_file$datapath)*1.3)*250}
  )
}
