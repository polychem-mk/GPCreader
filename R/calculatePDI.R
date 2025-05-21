calculatePDI <-
function(dataMWD,
                           name = NULL,
                           shift.start = 0,
                           shift.end = 0,
                           use.fitted = FALSE,
                           add.calib.info = TRUE
  ){

  # (1)  input   --------------------------------------------------

  # (1.1) shift.start
  # shift.start must be a numeric value between -20 and 35 or NA;
  # otherwise shift.start will be automatically set to the default value.

  shift_start_min = -20
  shift_start_max = 35

  if(all(is.vector(shift.start),
         length(shift.start) == 1)){

    if(is.numeric(shift.start)){
      if( shift.start < shift_start_min | shift.start > shift_start_max){
        shift.start = 0 # set to the default value
      }

    }else{
      shift.start = 0 # set to the default value
    }

  }else{
    shift.start = 0 # set to the default value
  }

  # (1.2) shift.end
  # shift.end must be a numeric value between -20 and 35 or NA;
  # otherwise shift.end will be automatically set to the default value.

  shift_end_min = -20
  shift_end_max = 35

  if(all(is.vector(shift.end),
         length(shift.end) == 1)){

    if(is.numeric(shift.end)){
      if( shift.end < shift_end_min | shift.end > shift_end_max){
        shift.end = 0 # set to the default value
      }

    }else{
      shift.end = 0 # set to the default value
    }

  }else{
    shift.end = 0 # set to the default value
  }

  # (1.3) use.fitted
  # use.fitted must be a logical and have a length of 1 or a value that can
  # be coerced to a logical, otherwise use.fitted will be automatically set
  # to the default value.
  if(any(!is.vector(use.fitted, mode = "logical"),
         length(use.fitted) != 1,
         is.na(use.fitted))){

    if(any(!is.vector(use.fitted),
           length(use.fitted) != 1,
           is.na(use.fitted))){
      use.fitted = FALSE  # set to the default value

    }else{
      use.fitted = as.logical(use.fitted)
      if(is.na(use.fitted)){
        use.fitted = FALSE  # set to the default value
      }
    }
  }

  # (1.4) add.calib.info
  # add.calib.info must be a logical and have a length of 1 or a value that can
  # be coerced to a logical, otherwise add.calib.info will be automatically set
  # to the default value.
  if(any(!is.vector(add.calib.info, mode = "logical"),
         length(add.calib.info) != 1,
         is.na(add.calib.info))){

    if(any(!is.vector(add.calib.info),
           length(add.calib.info) != 1,
           is.na(add.calib.info))){
      add.calib.info = TRUE  # set to the default value

    }else{
      add.calib.info = as.logical(add.calib.info)
      if(is.na(add.calib.info)){
        add.calib.info = TRUE  # set to the default value
      }
    }
  }

  # (1.5) dataMWD
  # dataMWD must be a list of 2; if it is not, filter_peaks() returns NULL.
  if(any(!"list"  %in% class(dataMWD) , length(dataMWD) != 2)){

    warning(" 'dataMWD' must be the output of calculateMWD() function ... NULL is returned ...")
    return(NULL)
  }

  # (1.5.1) first is the list of the data frames.
  if(!"list"  %in% class(dataMWD[1])){
    warning(" 'dataMWD' must be the output of calculateMWD() function ... NULL is returned ...")
    return(NULL)
  }

  # (1.5.2) the second is a numeric vector of length 4
  coeffs = unlist(dataMWD[2])

  if(any(!is.vector( coeffs, mode = "numeric"), length(coeffs) != 4)){
    warning(" 'dataMWD' must be the output of calculateMWD() function ... NULL is returned ...")
    return(NULL)
  }

  Intercept = coeffs[1]
  Coef1 = coeffs[2]
  Coef2 = coeffs[3]
  Coef3 = coeffs[4]

  # (1.6) name
  # name must be a character vector of the same length as dataMWD[[1]];
  # otherwise, name will be automatically set to integers starting from 1
  if(any(!is.vector(name, mode = "character"),
         length(name) != length(dataMWD[[1]]))){
    name = 1:length(dataMWD[[1]])
  }

  # (2) map by file -------------------------------------

  dataMWD_1 = dataMWD[[1]]
  rersults = map(1:length(dataMWD_1), function(j) {

    dat = dataMWD_1[[j]]
    file_name = name[j]

    if(any(!class(dat) %in% c("tbl_df","tbl","data.frame"))){
      rersults_i <- data.frame(FileName = file_name,
                               peak = NA,
                               PDI = NA,
                               Mn = NA,
                               Mw = NA,
                               Mp = NA,
                               time1 = NA,
                               time2 = NA,
                               time_max = NA,
                               Af = NA,
                               Tf = NA)
      return(rersults_i)
    }

    if( any(nrow(dat) < 10, !"yi" %in% colnames(dat), is.null(dat))){
      # If no peaks were detected:
      rersults_i <- data.frame(FileName = file_name,
                               peak = NA,
                               PDI = NA,
                               Mn = NA,
                               Mw = NA,
                               Mp = NA,
                               time1 = NA,
                               time2 = NA,
                               time_max = NA,
                               Af = NA,
                               Tf = NA)
      return(rersults_i)
    }

    # define the peak region to be used for calculations (peak ID is not 0) and
    # intensity; if peak tails were added, y_fit will be used for further
    # calculations:
    if(all("y_fit" %in% colnames(dat), use.fitted)){
      dat = dat %>%
        group_by(peak) %>%
        mutate(peak_id = ifelse(fitted == 1, peak, peak_id)) %>%
        mutate(max_y = ifelse(peak_id != 0, max(yi), 0 )) %>%
        ungroup() %>%
        mutate(yi = ifelse(!is.na(y_fit), y_fit, yi))

    }else{
      dat = dat %>%
        group_by(peak_id) %>%
        mutate(max_y = ifelse(peak_id != 0, max(yi), 0 )) %>%
        ungroup()
    }

    # The start and the end of the peaks:
    peaks_summ = dat %>%
      filter(peak_id != 0) %>%
      group_by(peak_id) %>%
      summarise(peak = first(peak),
                start_x = first(x),
                end_x = last(x))

    # Time (x) and intensity (yi) at the peak maximum:
    max_xy = dat %>%
      filter(peak_id != 0) %>%
      group_by(peak_id) %>%
      filter(yi == max(yi)) %>%
      summarise(max_y = first(yi), max_x = first(x))

    # assymetry factor = (width0.1 - f0.1)/f0.1
    assym_factor = dat %>%
      filter(peak_id != 0 ) %>%
      left_join(max_xy, by = "peak_id") %>%
      select(x, yi, max_x, peak_id) %>%
      group_by(peak_id) %>%
      mutate(left_tail = ifelse(x < max_x, 1, 0)) %>%
      filter(yi >= max(yi)*0.1) %>%
      group_by(peak_id, left_tail) %>%
      arrange(yi) %>%
      filter(x == first(x)) %>%
      group_by(peak_id) %>%
      arrange(left_tail) %>%
      summarise(af =  round((first(x)- first(max_x)) /(first(max_x) - last(x)),2 ) )

    # Tailing factor =  width0.05 / 2*f0.05.
    tail_factor = dat %>%
      filter(peak_id != 0 ) %>%
      left_join(max_xy, by = "peak_id") %>%
      select(x, yi, max_x, peak_id) %>%
      group_by(peak_id) %>%
      mutate(left_tail = ifelse(x < max_x, 1, 0)) %>%
      filter(yi >= max(yi)*0.05) %>%
      group_by(peak_id, left_tail) %>%
      arrange(yi) %>%
      filter(x == first(x)) %>%
      group_by(peak_id) %>%
      arrange(left_tail) %>%
      summarise(tf=  round((first(x)- last(x)) /( (first(max_x) - last(x))*2 ),2 ) )

    # (3) map by peak.
    rersults_i = map(1:nrow(peaks_summ), function(i) {

      Time1 = peaks_summ$start_x[i]    # start of the peak
      Time2 = peaks_summ$end_x[i]      # end of the peak
      time_max = max_xy$max_x[i]       # peak maximum

      # if there is the user input to change start/end time of the peak
      if(shift.start != 0){
        Time1 = Time1 + (Time2-Time1)*(shift.start/100) # change peak start time
        Time1 = ifelse(Time1 < 0, 0, Time1 )
      }

      if(shift.end != 0){
        Time2 <- Time2 - (Time2-Time1)*(shift.end/100) # change peak end time
      }

      # If intensity goes below a baseline it will be dropped from the calculations
      peak_i = dat %>% filter(peak == peaks_summ$peak[i]) %>%
        filter(x >= Time1 & x <= Time2) %>% filter(yi > 0 )

      Time1 = min(peak_i$x)   # start of the peak
      Time2 = max(peak_i$x)

      is_Mi = "Mi" %in% colnames(dat)
      if(is_Mi){
        is_Mi = sum(!is.na(dat$Mi)) > 0
      }

      if( all(!is.na(Intercept),
              !is.na(Coef1),
              is_Mi )){

        peak_i <- peak_i %>%
          summarise(Mn = round( sum(yi)/sum(yi/Mi ) ),      # number-average molecular weight
                    Mw = round( sum(yi*Mi)/sum(yi) ) ) %>%  # weight-average molecular weight
          mutate(peak = i,
                 Mp = round( 10^(Coef3*(time_max^3)+ Coef2*(time_max^2)+ Coef1*time_max+ Intercept )),
                 PDI = ifelse(Mw != 0 | Mn != 0  ,round(Mw/Mn, 2), NA), # polydispersity index
                 time1 = round(Time1, 2),
                 time2 = round(Time2, 2),
                 time_max = round(time_max, 2) )
      }else{
        peak_i <- data.frame(Mn = NA,         # number-average molecular weight
                             Mw = NA ) %>%    # weight-average molecular weight
          mutate(peak = i,
                 Mp = NA,
                 PDI = NA,   # polydispersity index
                 time1 = round(Time1, 2),
                 time2 = round(Time2, 2),
                 time_max = round(time_max, 2) )
      }
      peak_i

    } ) %>% bind_rows()

    rersults_i = rersults_i %>%
      mutate(Af = assym_factor$af,
             Tf = tail_factor$tf,
             FileName = file_name) %>%
      select(FileName, peak, PDI, Mn,  Mw, Mp,  time1, time2, time_max,
              Af, Tf)

  }) %>% bind_rows()

  if(add.calib.info){
    rersults$intercept = round(Intercept, 5)
    rersults$coef1 = round(Coef1, 5)
    rersults$coef2 = round(Coef2, 5)
    rersults$coef3 = round(Coef3, 5)
  }

  #       Output:
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

  rersults
}
