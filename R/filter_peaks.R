filter_peaks <-
function(dataGPC,
                         first.peak = FALSE,
                         above = 0,
                         tails = 0.02){

  # (1) Check input   -----------------------

  # (1.1) tails
  # 'tails' is the proportion of the maximum intensity.
  # If more than 5% (tails*100%) of a peak's tail is missing, filter_peaks()
  # will try to fit and add the missing values.

  # 'tails' must be a numeric value between 0.02 and 0.3 (from 2 to 30%) or
  # FALSE; otherwise 'tails' will be automatically set to the default value.

  fit_tails_min = 0.02
  fit_tails_max = 0.3

  if(all(is.vector(tails),
         length(tails) == 1)){

    if(is.numeric(tails)){
      if( tails < fit_tails_min | tails > fit_tails_max){
        tails = 0.05    # set to the default value
      }

    }else if(tails != FALSE){
      tails = 0.05      # set to the default value
    }

  }else{
    tails = 0.05        # set to the default value
  }

  # (1.2) first.peak
  # By default filter_peaks() keeps all detected peaks, if the peak of interest
  # is the first peak in the chromatogram, then first.peak can be set to TRUE.
  # In this case, all other peaks will be excluded from further calculations.

  # first.peak must be a logical and have a length of 1 or a value that can be coerced
  # to a logical, otherwise first.peak will be automatically set to the default value.
  if(any(!is.vector(first.peak, mode = "logical"),
         length(first.peak) != 1)){

    if(any(!is.vector(first.peak),
           length(first.peak) != 1)){
      first.peak = FALSE

    }else{
      first.peak = as.logical(first.peak)
      if(is.na(first.peak)){
        first.peak = FALSE
      }
    }
  }

  # (1.3) above
  # To keep only peaks above a certain proportion of the largest peak, set
  # above to that value. Thus, if above is set to 0.1, then all
  # peaks below 10% of the maximum intensity will be discarded.

  # 'above' must be a numeric value between 0 and 0.99;
  # otherwise 'above' will be automatically set to the default value.

  peaks_above_min = 0
  peaks_above_max = 0.99

  if(all(is.vector(above, mode = "numeric"),
         length(above) == 1)){

    if( above < peaks_above_min | above > peaks_above_max){
      above = 0      # set to the default value
    }

  }else{
    above = 0        # set to the default value
  }

  # (1.4) dataGPC
  # dataGPC must be a list of data frames or a data frame to be converted to a list.
  # If it is not, filter_peaks() returns NULL.

  if(any(class(dataGPC) %in% c("tbl_df","tbl","data.frame"))){

    if(any(!c("x","y","ybl","yi","peak_id") %in% colnames(dataGPC),
           nrow(dataGPC) < 50)){
      warning(" 'list_df' must be a list of data frames ... NULL is returned ...")
      return(NULL)

    }else{
      dataGPC = list(as.data.frame(dataGPC))
    }
  }

  if(!is.list(dataGPC)){
    warning(" 'list_df' must be a list of data frames ... NULL is returned ...")
    return(NULL)
  }

  # (2) Functions  ------------------

    # (2.1) parameters()
    # finds the slope, intercept, and x (time) at half y (half intensity)

    # arguments:
    # - df is a data frame with x and yi for one peak;
    # - tail_missing_prop is proportion y of maximum y at the beginning/end of the peak;
    # - x1 for the left side of the peak is the peak start time or peak shoulder time,
    #      for the right side of the peak this is the time of the peak's maximum;
    # - x2 is the maximum peak time for the left side, and for the right side of
    #      the peak this value is the peak end time or peak shoulder time;
    # - max_intensity  maximum y;
    # - up is 1 for the left side of the peak and up = -1 for the right side.

    parameters = function(df, tail_missing_prop, x1, x2, max_intensity, dif_x, up = 1){
      # At half peak intensity, the slope and intercept have their maximum
      # absolute value. We find the slopes for each yi between time x1 and x2
      # and compute the intercept for the maximum slope.
      slope = df %>%
        mutate(slope_local = (lead(yi) - yi)/dif_x ) %>%
        drop_na() %>%
        filter( x >= x1  & x <= x2) %>%
        arrange(slope_local*up) %>%
        last() %>%
        mutate(intercept_local = yi - x*slope_local)

      if(tail_missing_prop < 0.45 & nrow(slope) == 1){
        # If less than 45% of the tail is missing, we set the left slope to
        # the maximum of the "local slope" (and the intercept to the minimum
        # of the "local intercept"); and vice versa for the right tail.
        intercept = slope$intercept_local
        slope = slope$slope_local

      }else{
        # If more than 45% of the tail is missing, the peak merges with another peak
        # and we cannot determine the slope exactly, but we can make some approximation
        # by calculating the slope between the point with maximum local slope and
        # the point with zero intensity, which is between the peak center and the
        # start/end of the peak, which is calculated from the maximum slope and the
        # intercept (maximum absolute values).
        # 'approx_narrow' is a factor, the higher this value, the closer
        # the peak tails are to the center.
        approx_narrow = 0.9

        approx_narrow = ifelse(tail_missing_prop > 0.8, tail_missing_prop+0.2, approx_narrow)

        slope = slope %>%
          mutate(y_prop = yi/max_intensity) %>%
          mutate( x_at_y0 = (-intercept_local/slope_local) + up*approx_narrow*(y_prop-0.5)/(1/y_prop) ) %>%
          mutate(slope = yi/(x - x_at_y0)) %>%
          mutate(intercept = yi - x*slope)

        intercept = slope$intercept
        slope = slope$slope
      }

      # x at the half intensity
      half_x = (max_intensity/2  - intercept)/slope

      data.frame(half_x = half_x,
                 slope = slope,
                 intercept = intercept)
    }

    # (2.2) find_coeffs()  fits a bell curve:
    # y =  max_y/(1 + exp(slope_left*(x_50_left - x) ) ) -
    #                             max_y/(1 + exp(slope_right*(x - x_50_right)) )

    # arguments:
    # - df is a data frame with x and yi for one peak;
    # - coeffs is a vector with
    #             - maximum intensity,
    #             - left slope,
    #             - left half x,
    #             - right slope,
    #             - right half x;
    # - pars is a vector with initial parameters for the optim() function.

    find_coeffs =  function(df, pars = c(1, 1, 1, 1), coeffs){
      # func() calculates the sum of squares of the residuals
      func = function( par)   {
        fitted = par[1]*coeffs[1]/(1 + par[2]*exp(par[3]*coeffs[2]*(coeffs[3] - x) ) ) -
          par[1]*coeffs[1]/(1 + par[2]*exp(par[4]*coeffs[4]*(x - coeffs[5])) )

        sum((y - fitted)^2)
      }

      # Set x and y
      df = df %>% filter(y_full >= 0) %>%
        slice_head(prop = 0.85)

      x = df %>%  pull(x)
      y = df %>%  pull(y_full)

      # Try to find coefficients
      try_fit =   tryCatch(
        expr = {
          fit_coeffs = optim(pars, func, method = "BFGS", control=list(reltol = 1e-10))

        },
        error =  function(e){
          return(NA)
        }
      )

      if( any(!is.na(try_fit)) ){
        return(fit_coeffs$par)

      }else{
        return(NULL)
      }
    }

  # (3) Map list with GPC files (dataGPC)  -----------------

  map(dataGPC, function(dat){

    # (3.1) Check dat
    if(any(is.null(dat), !is.data.frame(dat))){
      return(NULL)
    }

    # If the data does not contain x and y columns, filter_peaks() returns NULL
    if(any(ncol(dat) < 2,
           nrow(dat) < 12,
           !c("x", "y")  %in%  colnames(dat) )){
      return(NULL)
    }

    # If there are no peaks, filter_peaks() returns the input data
   if(any(!c("yi", "peak_id", "type")  %in%  colnames(dat) )){
      return(dat)
    }

    # The time difference between two measurements
    dif_x  =  mean(diff(dat$x))

    # (3.2) All peaks summary
    peak_summary = dat %>% filter(peak_id != 0) %>%
      group_by(peak_id) %>%
      mutate(peak_width = n()) %>%
      filter(type %in% c("start", "end", "max", "shoulder")) %>%
      mutate(trend = ifelse(x <= x[which(type == "max")], "up", "down")) %>%
      ungroup() %>%
      group_by(peak_id, trend, type) %>%
      filter(yi == max(yi)) %>%
      distinct(type, .keep_all = TRUE) %>%
      ungroup() %>%
      group_by(peak_id) %>%
      # tail_rel_y is the intensity proportion of the maximum intensity at
      # the beginning/end (or at the shoulder) of the peak
      mutate(tail_rel_y = yi/max(yi)) %>%
      # The lim_x variable will be used to add the area before and after
      # the peak where the fitted y values need to be added.
      mutate(lim_x = x +
               ifelse(trend == "up", -1, 1)*100*dif_x*tail_rel_y +
               ifelse(trend == "up", -1, 1)*20*dif_x ) %>%

      ungroup() %>%
      mutate(rel_y = yi/max(yi))

    if(above > 0){

      # If 'above' is not 0, then all peak_ids corresponding to peaks
      # below the 'above' proportion of maximum intensity will be removed
      # from the calculations.
      peak_ID_rm = peak_summary %>%
        filter(type == "max") %>%
        filter(rel_y < above) %>%
        pull(peak_id)

      peak_summary = peak_summary %>%
        filter(!peak_id %in% peak_ID_rm)
    }

    # All peak_id(s) for the chromatogram
    peak_IDs = unique(peak_summary$peak_id)

    peak_summary = peak_summary %>%
      filter(type!= "max")

    if(first.peak == TRUE){
      # if only the first peak is the peak of interest:
      peak_IDs = peak_IDs[1]
    }

    #   (4) Map by peaks (peak_id)
    peaks_stack = map(1:length(peak_IDs), function(i){

      ### Summary for a peak i:

      # peak_id     peak ID;
      # peak_width  is the number of measurements for one peak;
      # type        "start" and "end" are the start and the end of the peak
      #             (where peak_id is not 0); if there is a "shoulder", then
      #             the x of the start or end is moved to the x of the shoulder;
      # tail_rel_y  is relative intensity at the start/end of the peak
      #             (or at the shoulder), if it is high then the tails
      #             are above the baseline and can be added;
      # lim_x       adds x to the right and left of the peak, these
      #             values will be used if the peak tails are missing
      #             and the fitted values need to be added.

      summary_i = peak_summary %>% filter(peak_id == peak_IDs[i])

      # the number of points that will be added to correct the peak region
      n_left = summary_i %>% filter(type == "start") %>%
        summarise(round(tail_rel_y*peak_width)) %>% pull()

      n_right = summary_i %>% filter(type == "end") %>%
        summarise(round(tail_rel_y*peak_width)) %>% pull()


      # x limits for right and left sides of the peak
      time_lim_start = summary_i %>% filter(trend == "up") %>%
        summarise(min(lim_x)) %>% pull()
      time_lim_start = ifelse(time_lim_start <= 0, 0.1, time_lim_start)

      time_lim_end = summary_i %>% filter(trend == "down") %>%
        summarise(max(lim_x)) %>% pull()


      # The peak i, including the area before and after the peak
      peak_i = dat %>%
        filter( x > time_lim_start  & x < time_lim_end) %>%
        mutate(other_peak = ifelse(peak_id != peak_IDs[i] & peak_id != 0, 1, 0) ) %>%

        mutate(peak_id = ifelse(peak_id == peak_IDs[i], peak_id, 0) ) %>%
        # Add a variable "peak" that will have the same value as peak_id but
        # includes the regions before and after the original peak. If missing
        # tails need to be added, the x values from this region of the
        # chromatogram will be added to find the fitted intensity (yi values).
        mutate(peak = i) %>%
        # Change peak_id so that it is in the numeric order and the same as 'peak'
        mutate(peak_id = ifelse(peak_id != 0, i, 0))   %>%
        # To have only one "max"/"end"/"start" type for the peak:
        mutate(type = ifelse(peak_id == 0, "other", type))   %>%
        # Add a 'trend' variable: 1 indicates the left side of the peak
        # (before the maximum), and -1 indicates the right side of the peak
        # (after the maximum)
        mutate(trend = ifelse(x < x[which(type == "max")], 1, -1)) %>%
        group_by(trend) %>%
        # The 'other_peak' variable ensures that we don't add part of another peak.
        mutate(other_peak = ifelse(trend == -1 ,
                                   cummax(other_peak),
                                   other_peak)) %>%
        arrange(desc(x)) %>%
        mutate(other_peak = ifelse(trend == 1 ,
                                   cummax(other_peak),
                                   other_peak)) %>%
        arrange(x) %>%
        ungroup() %>%
        # Add a relative intensity ('rel_y')
        mutate(rel_y = yi/yi[which(type == "max")])  %>%
        group_by(peak_id, trend) %>%
        # 'n' indicates how far the x values are outside the peak region
        mutate(n = case_when(
          trend == 1 & peak_id == 0 ~ row_number(desc(x)),
          trend == -1 & peak_id == 0 ~ row_number(x),
          .default = 0
        )) %>%
        ungroup() %>%
        # To correct the peak region (where peak_id is not 0), we add
        # a variable 'add_peak_id' which has a value of 1 if the measurement
        # is close to the start/end of the peak, has an intensity above
        # the baseline and the same trend as the side of the peak.
        mutate(add_peak_id = case_when(
          trend == 1 &
              n <= n_left + 5 &
              rel_y > 0.01 &
              rel_y < lead(rel_y) &
              rel_y < summary_i$yi[which(summary_i$type == "start")] &
              peak_id == 0 &
              other_peak == 0     ~ 1,
          trend == -1 &
              n <= n_right + 5 &
              rel_y > 0.01 &
              rel_y < lag(rel_y) &
              rel_y < summary_i$yi[which(summary_i$type == "end")] &
              peak_id == 0 &
              other_peak == 0      ~ 1,
          .default = 0
        )) %>%
        mutate(peak_id = ifelse(add_peak_id == 1, peak, peak_id))

      # add_tails specifies whether tails should be added: if add_tails is
      # not FALSE and the peak start/end intensity (or shoulder intensity) is
      # greater than the add_tails value
      if(tails == FALSE){
        add_tails = FALSE
      }else{
        add_tails = any(summary_i$tail_rel_y > tails)
      }

      if(!add_tails){
      # If there is no need to add tails, the function returns a data frame with peaks.
        peak_i = peak_i %>%
        mutate(y_fit = NA, fitted = 0) %>%
        select(x, yi, peak_id, peak, y_fit, fitted, type)

        return(peak_i)
      }

      #                                 Add tails
      max_y = peak_i %>%
        filter(peak_id != 0) %>%
        filter(yi == max(yi)) %>%
        first()
      max_x = max_y %>% pull(x)    # time of maximum peak intensity i (x)
      max_y =  max_y %>% pull(yi)  # maximum intensity of peak i (y)

      # proportions of the missing parts
      left_tail_missing_prop = summary_i %>% filter(trend == "up") %>%
        summarise(max(tail_rel_y)) %>% pull()

      right_tail_missing_prop = summary_i %>% filter(trend == "down") %>%
        summarise(max(tail_rel_y)) %>% pull()

      # If more than 95% of one side of a peak is missing, no further
      # calculations will be performed for that peak.
      if(any(left_tail_missing_prop > 0.95, right_tail_missing_prop >0.95)){
        peak_i = peak_i %>%
          mutate(y_fit = NA,  fitted = -1) %>%
          select(x, yi, peak_id, peak, y_fit, fitted, type)

        return(peak_i)
      }

      #                  Find slope, intercept and x at the half y

      # Unlike HPLC chromatograms, GPC chromatograms are not necessarily
      # symmetrical, so to fit the equation we need parameters for each side
      # of the peak (left and right): slope and time (x) at half intensity (y).

      # LEFT part ___________
      # 'x_start' is the start of the peak, but if there is a shoulder, then
      # 'x_start' is x (time) of the left shoulder
      x_start = summary_i %>% filter(trend == "up") %>%
        summarise(max(x)) %>% pull()

      # Find the slope and the intersecept for the left side of the peak
      # between the start (or left shoulder) and the maximum of the peak.
      parameters_left = parameters(df = peak_i,
                                   x1 = x_start,
                                   x2 = max_x,
                                   dif_x = dif_x,
                                   max_intensity = max_y,
                                   tail_missing_prop = left_tail_missing_prop,
                                   up = 1)
      # RIGHT part __________
      # 'x_end' is the end of the peak, but if there is a shoulder, then
      # 'x_end' is x (time) of the right shoulder
      x_end = summary_i %>% filter(trend == "down") %>%
        summarise(min(x)) %>% pull()

      # Find the slope and the intersecept for the right side of the peak
      # between the maximum and the end (or right shoulder) of the peak.
      parameters_right = parameters(df = peak_i,
                                    x1 = max_x,
                                    x2 = x_end,
                                    dif_x = dif_x,
                                    max_intensity = max_y,
                                    tail_missing_prop = right_tail_missing_prop,
                                    up = -1)

      if(any( is.na(parameters_left$slope),
              is.na(parameters_right$slope),
              parameters_left$slope <= 0,
              parameters_right$slope >= 0)){

        # If the slopes are incorrect, no further calculations will be
        # performed for that peak.
        peak_i = peak_i %>%
          mutate(y_fit = NA,  fitted = -1) %>%
          select(x, yi, peak_id, peak, y_fit, fitted, type)

        return(peak_i)

      }else{

        # AB line
        # If more than 10% of the tail is missing, we add abline to replace
        # the missing tail, otherwise optim() may fail to fit the bell curve.
        peak_i = peak_i %>%
          mutate(y_full = ifelse(x <= max_x,
                                 x*parameters_left$slope +parameters_left$intercept,
                                 x*parameters_right$slope + parameters_right$intercept))

        # LEFT tail ____________
        if(left_tail_missing_prop > 0.1){

          x_mix1 =  peak_i %>%
            filter(x <= max_x) %>%
            filter(y_full < yi & lead(y_full) >= lead(yi))

          if(nrow(x_mix1) == 0){
            x_mix1 = x_start

          }else if(nrow(x_mix1) == 1){
            x_mix1 = x_mix1 %>% pull(x)

          }else{
            x_mix1 = x_mix1 %>% last() %>% pull(x)
          }

        }else{
          x_mix1 = x_start
        }

        # RIGHT tail ____________
        if(right_tail_missing_prop > 0.1){

          x_mix2 =  peak_i %>%
            filter(x >= max_x) %>%
            filter(y_full < yi & lag(y_full) >= lag(yi))

          if(nrow(x_mix2) == 0){
            x_mix2 = x_end

          }else if(nrow(x_mix2) == 1){
            x_mix2 = x_mix2 %>% pull(x)

          }else{
            x_mix2 = x_mix2 %>% first() %>% pull(x)
          }

        }else{
          x_mix2 = x_end
        }

        peak_i = peak_i %>%
          mutate(y_full = ifelse(x >= x_mix1 & x <= x_mix2, yi, y_full))

        #                       Fit bell curve
        # try to fit:
        # y =  max_y/(1 + exp(slope_left*(x_50_left - x) ) ) -
        #                     max_y/(1 + exp(slope_right*(x - x_50_right)))

        bell_coeffs2 = find_coeffs(peak_i,
                                   coeffs = c(
                                     max_y ,
                                     parameters_left$slope,
                                     parameters_left$half_x ,
                                     parameters_right$slope,
                                     parameters_right$half_x
                                   ),
                                   pars = c(1, 1, 1, 1))

        if(all(!is.null(bell_coeffs2))){
          peak_i = peak_i %>%
            mutate(
              y_bell2 = bell_coeffs2[1]*max_y/
                (1 + bell_coeffs2[2]*exp(bell_coeffs2[3]* parameters_left$slope*
                                           (parameters_left$half_x - x)))  -
                bell_coeffs2[1]*max_y/
                (1 + bell_coeffs2[2]*exp(bell_coeffs2[4]*parameters_right$slope*
                                           (x - parameters_right$half_x))) )

          # check if the fitted curve is bell shaped:
          err2  = peak_i %>%
            filter(y_bell2 == max(y_bell2, na.rm = T)) %>%
            select(x)

          if(any(nrow(err2) > 5, err2$x < x_start, err2$x > x_end)){
            err2 = 1
          }else{
            # the difference between the fitted values and y:
            err2  = peak_i %>% filter(y_full > max_y*0.1 ) %>%
              summarise( mean(abs(y_full - y_bell2))/max_y) %>%
              pull()

            err2  = ifelse(is.na(err2), 1, err2)
          }
        }else{
          err2 = 1
        }

        if(err2 < 0.021){
          # if err2 is less than 2%, add y_bell2 as fitted values
          peak_i = peak_i %>%
            mutate(y_bell = y_bell2,
                   fitted = 1)
        }else{
          # If the error exceeds 7%, the fitted values will not be used:
          err2 = ifelse(err2 > 0.07, 1, err2)

          # If the error exceeds 2%, try to fit:
          # y =  max_y/(1 + exp(log(abs(slope_left))*(x_50_left - x) ) ) -
          #          max_y/(1 + exp(-log(abs(slope_right))*(x - x_50_right)))

          bell_coeffs1 = find_coeffs(peak_i,
                                     coeffs = c(
                                       max_y ,
                                       log(abs( parameters_left$slope)),
                                       parameters_left$half_x ,
                                       (-log(abs(parameters_right$slope))),
                                       parameters_right$half_x
                                     ),
                                     pars = c(1, 1, 1, 1))

          if(all(!is.null(bell_coeffs1))){
            peak_i = peak_i %>%
              mutate(
                y_bell1 = bell_coeffs1[1]*max_y/
                  (1 + bell_coeffs1[2]*exp(bell_coeffs1[3]* log(abs( parameters_left$slope))*
                                             (parameters_left$half_x - x)))  -
                  bell_coeffs1[1]*max_y/
                  (1 + bell_coeffs1[2]*exp(bell_coeffs1[4]*(-log(abs(parameters_right$slope)))*
                                             (x - parameters_right$half_x)))
              )

            err1  = peak_i %>%
              filter(y_bell1 == max(y_bell1)) %>%
              select(x)

            if(any(nrow(err1) > 5, err1$x < x_start, err1$x > x_end)){
              err1 = 1
            }else{
              err1  = peak_i %>% filter(y_full > max_y*0.1 ) %>%
                summarise( mean(abs(y_full - y_bell1))/max_y) %>%
                pull()

              err1  = ifelse(is.na(err1), 1, err1)
            }

          }else{
            err1 = 1
          }

          if(err1 < 0.021){
            # if err1 is less than 2%, add y_bell1 as fitted values
            peak_i = peak_i %>%
              mutate(y_bell = y_bell1,
                     fitted = 1)

          }else{
            # If the error exceeds 7%, the fitted values will not be used:
            err1 = ifelse(err1 > 0.07, 1, err1)

            # If both errors exceed 2%, try to fit:
            # y =  max_y/(1 + exp((slope*4/max_y)*(x - x_50_right) -
            #              max_y/(1 + exp((slope*4/max_y)*(x - x_50_right)))

            bell_coeffs3 = find_coeffs(peak_i,
                                       coeffs = c(
                                         max_y ,
                                         abs(parameters_left$slope*4/max_y),
                                         parameters_left$half_x ,
                                         abs(parameters_right$slope*4/max_y),
                                         parameters_right$half_x
                                       ),
                                       pars = c(1, 1, 1, 1))

            if(all(!is.null(bell_coeffs3))){
              peak_i = peak_i %>%
                mutate(
                  y_bell3 = bell_coeffs3[1]*max_y/
                    (1 + bell_coeffs3[2]*exp(bell_coeffs3[3]* abs(parameters_left$slope*4/max_y)*
                                               (parameters_left$half_x - x)))  -
                    bell_coeffs3[1]*max_y/
                    (1 + bell_coeffs3[2]*exp(bell_coeffs3[4]*(abs(parameters_right$slope*4/max_y))*
                                               (x - parameters_right$half_x)))
                )

              err3  = peak_i %>%
                filter(y_bell3 == max(y_bell3)) %>%
                select(x)

              if(any(nrow(err3) > 5, err3$x < x_start, err3$x > x_end)){
                err3 = 1
              }else{
                err3  = peak_i %>% filter(y_full > max_y*0.1 ) %>%
                  summarise( mean(abs(y_full - y_bell3))/max_y) %>%
                  pull()

                err3  = ifelse(is.na(err3), 1, err3)
              }

            }else{
              err3 = 1
            }

            # If the error exceeds 7%, the fitted values will not be used:
            err3 = ifelse(err3 > 0.07, 1, err3)

            # Errors for all three fitted bell curves:
            errs = c(err1, err2, err3)

            if(sum(errs) == 3){
              # If all errors exceed 7% fitted values will not be added:
              peak_i = peak_i %>%
                mutate(y_fit = NA,
                       fitted = -1) %>%
                select(x, yi, peak_id, peak, y_fit, fitted, type)

              return(peak_i)

            }else{

              # add the fitted values that give the minimum error:
              if(err2 == min(errs)){
                peak_i = peak_i %>%
                  mutate(y_bell = y_bell2,
                         fitted = 1)

              }else if(err1 == min(errs)){
                peak_i = peak_i %>%
                  mutate(y_bell = y_bell1,
                         fitted = 1)

              }else{
                peak_i = peak_i %>%
                  mutate(y_bell = y_bell3,
                         fitted = 1)
              }
            }
          }
        }

        # Add fitted values.
        peak_i = peak_i %>% mutate(y_fit = yi)

        # Mp is calculated from the time (x) at the peak maximum, therefor
        # we only add the missing tails to preserve the peak shape (optim()
        # does not necessarily keeps the maximum position of the bell curve).

        # LEFT tail ____________
        if(left_tail_missing_prop > tails){

          x_mix1 =  peak_i %>%
            filter(x < max_x) %>%
            filter(peak_id != 0) %>%
            mutate(fit_dif = abs(y_fit - yi))

          x_mix1_alt = x_mix1 %>%
            filter(y_bell > yi & lead(y_bell) <= lead(yi)) %>%
            last() %>% pull(x)

          x_mix1 = x_mix1 %>%
            filter(y_bell < yi & lead(y_bell) >= lead(yi))

          if(nrow(x_mix1) == 0){
            x_mix1 = ifelse(!is.na(x_mix1_alt), x_mix1_alt,  x_start)

          }else {
            x_mix1 = x_mix1 %>%
              filter(fit_dif == min(fit_dif)) %>%
              last() %>% pull(x)
          }

          peak_i = peak_i %>%
            mutate(y_fit = ifelse(x <= x_mix1, y_bell, y_fit))
        }

        # RIGHT tail ____________
        if(right_tail_missing_prop > tails){

          x_mix2 =  peak_i %>%
            filter(x > max_x) %>%
            filter(peak_id != 0) %>%
            mutate(fit_dif = abs(y_fit - yi))

          x_mix2_alt = x_mix2 %>%
            filter(y_bell > yi & lag(y_bell) <= lag(yi)) %>%
            first() %>% pull(x)

          x_mix2 = x_mix2 %>%
            filter(y_bell < yi & lag(y_bell) >= lag(yi))

          if(nrow(x_mix2) == 0){
            x_mix2 = ifelse(!is.na(x_mix2_alt), x_mix2_alt,  x_end)

          }else {
            x_mix2 = x_mix2 %>%
              filter(fit_dif == min(fit_dif)) %>%
              first() %>% pull(x)
          }

          peak_i = peak_i %>%
            mutate(y_fit = ifelse(x >= x_mix2, y_bell, y_fit))
        }

        # variable 'fitted' is
        #   1 if the peak is fitted
        #  -1 if fitting is not possible (error greater than 7%, more than
        #     95% of peak missing on one side)
        #   0 otherwise

        if(left_tail_missing_prop > tails & right_tail_missing_prop > tails){
          peak_i = peak_i %>%
            mutate(fitted = ifelse(y_fit <= max_y*0.001, 0, fitted))

        }else if(left_tail_missing_prop > tails ){
          peak_i = peak_i %>%
            mutate(fitted = ifelse(y_fit <= max_y*0.001 | x > x_end, 0, fitted))

        }else{
          peak_i = peak_i %>%
            mutate(fitted = ifelse(y_fit <= max_y*0.001 | x < x_start, 0, fitted))
        }

        # Check if there is noise due to merging yi and y_fit.
        noise =  peak_i %>%
          filter(  x < x_mix1 | x > x_mix2 ) %>%
          group_by(trend) %>%
          filter(y_fit == max(y_fit)) %>%
          mutate(noise = abs(y_fit - yi)/max_y)

        trend_up = which(noise$trend == 1)
        trend_down = which(noise$trend == -1)
        noise_l = ifelse(length(trend_up) != 0, max(noise$noise[trend_up] ), 0)
        noise_r = ifelse(length(trend_down) != 0, max(noise$noise[trend_down] ), 0)

        fix_merge = 0.01

        if(any(noise_l > fix_merge, noise_r > fix_merge)){
          # If the difference between yi and y_fit at the merge point
          # exceeds 1%, loess smoothing is applied.
          y_loess =  smooth(peak_i$y_fit, kind = "3RS3R")

          suppressWarnings({
            y_loess = loess(y_loess ~ peak_i$x, span = 0.15)
          })

          y_loess = y_loess$fitted

          if(nrow(peak_i) == length(y_loess)){
            peak_i = peak_i %>%
              mutate(y_fit = y_loess)
          }else{
            y_loess = data.frame(y_fit = y_loess,
                                 x = peak_i$x[1:length(y_loess)])

            peak_i = peak_i %>% select(x, yi, peak_id, peak, fitted, type) %>%
              right_join(y_loess, by = "x")
          }
        }
        peak_i = peak_i %>% select(x, yi, peak_id, peak, y_fit, fitted, type)
      }

    }) %>% bind_rows()

    # (5)   Output:
    # - x  is time;
    # - yi: if the peak is noisy yi is smoothed ybl, otherwise it is ybl;
    #       it will be used for the calculations if there are no fitted y values;
    # - peak_id is an integer, it is the peak id for the peak region above
    #           zero intensity, and peak_id is 0 for all other regions;
    # - peak is integer, same as peak_id, specifies the peak plus area before and
    #        after the peak (will be used for graphs)
    # - y_fit are the fitted y values that will be used for calculations;
    # - fitted:
    #              1 if the peak is fitted;
    #             -1 if fitting is not possible (error greater than 7%, more
    #                than 95% of peak missing on one side);
    #              0 otherwise.

    peaks_stack
  })
}
