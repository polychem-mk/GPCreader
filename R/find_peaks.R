find_peaks <-
function(dataGPC, syst.peak = 0, shoulder.sens = 2){

  # (1) input, variables   --------------------------

  # (1.1) check_x
  # check_x, if TRUE check first column x (time or volume) for time difference
  # between two measurements, is it constant (or close to constant)
  check_x = TRUE

  # (1.2) syst.peak
  # syst.peak must be a positive numeric value;
  # otherwise syst.peak will be automatically set to the default value.

  if(all(is.vector(syst.peak, mode = "numeric"),
         length(syst.peak) == 1)){

    if( syst.peak < 0){
      syst.peak = 0  # set to the default value
    }

  }else{
    syst.peak = 0    # set to the default value
  }

  # (1.3) noise_red
  # If the chromatogram is noisy and the noise exceeds the noise_red level, then
  # smoothing methods will be applied
  noise_red = 0.0001

  # (1.4) noise_lim
  # If the chromatogram is too noisy and the noise exceeds the noise_lim level,
  # then the find_peaks() function returns the original data (x and y)
  noise_lim = 0.01

  # (1.5) shoulder.sens
  # Limits to detect merged peaks (shoulder peaks)

  if(all(is.vector(shoulder.sens, mode = "numeric"),
         length(shoulder.sens) == 1)){

    if(shoulder.sens == 1){  # settings for calibration files
      shoulder.sens = c(0.20, 0.95, 1.10, 0.25)

    }else{  # settings for samples
      shoulder.sens = c(0.15, 0.85, 1.15, 0.05)
    }
  }else{
    shoulder.sens = c(0.15, 0.85, 1.15, 0.05)
  }

  # (1.6) dataGPC
  # dataGPC must be a list of data frames or a data frame to be converted to a list.
  # If it is not, find_peaks() returns NULL.

  if(any(class(dataGPC) %in% c("tbl_df","tbl","data.frame"))){

    if(any(!c("x", "y") %in% colnames(dataGPC),
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

  # (2) Map by file -------------------------
  map(dataGPC, function(dat){

    # (2.1) Check the format of the input data

    ### dat
    # check if dat is a data frame.
    # If there are less than 50 lines, then the resolution of the chromatogram
    # in dat is too low for further calculations.
    if(any(is.null(dat), !is.data.frame(dat),
           !c("x", "y") %in% colnames(dat),
           nrow(dat) < 50)){
      return(NULL)
    }

    ### Check x
    # The first column (x) is time (or  volume) and the time difference between
    # the two measurements should be approximately the same: the difference
    # between the maximum and minimum differences should be be less than 1% of
    # the average difference between the two measurements.

    dif_x = mean(diff(dat$x)) # x difference between the two measurements

    if(check_x){
      dif_x_pct  =  100*( max(diff(dat$x)) -  min(diff(dat$x)))/dif_x

      if(dif_x_pct > 1){
        return(dat)
      }
    }

    # (2.2)  Baseline

    # To find the baseline, we will look for areas of the chromatogram with the
    # smallest changes in intensity and add the variable bl, which indicates the
    # baseline (1) or another part, such as a peak (negative peak also) or noise (0).

    row_count = nrow(dat)

    dat = dat %>%
      # divide the chromatogram into 5 parts by retention time (variable x);
      # and look for the baseline in the first 4 parts of the chromatogram
      mutate(part = cut(x, breaks = 5, labels = FALSE)) %>%
      # divide the chromatogram by intensity (y variable) so that each
      # interval is 0.5% of the intensity range
      mutate(y_interval = cut(y, breaks = 200, labels = FALSE)) %>%
      group_by(y_interval) %>%
      # baseline y values are the most common in GPC chromatograms;
      # find proportion of y values that fall in each interval:
      mutate(interval_prop = n()/row_count) %>%
      ungroup()

    ### x and y for maximum intensity in the chromatogram:
    peak_y = max(dat$y)    # Intensity of the largest peak in the chromatogram (y)
    peak_x = dat$x[which.max(dat$y)[1]] # and the corresponding time (x)

    ### Keep 10 intervals for y values that have the largest number of y values;
    # this covers 5% of the maximum intensity. If the baseline has significant
    # drift, it may fall outside this range, but in this case the quality of the
    # chromatogram is too low for correct calculations.
    y_interval_top_prop = dat %>%
      distinct(y_interval, .keep_all = TRUE) %>%
      arrange(desc(interval_prop)) %>%
      slice_max(interval_prop, n = 10) %>%
      # a 'ratio' variable shows how the proportion of y values changes in each interval:
      mutate(ratio = interval_prop/max(interval_prop))

    ### We define bl as 1 if the y values fall into the intervals with the largest
    # number of y values. To find out how many intervals include the baseline,
    # we use the variable 'ratio'.
    if(y_interval_top_prop$ratio[2] < 0.2){

      # If the 'ratio'  drops sharply from interval 1 to interval 2, then most of
      # the y values fall in the same interval, and the baseline is flat.
      baseline_intervals = y_interval_top_prop$y_interval[1]

    }else if(y_interval_top_prop$ratio[2] < 0.5){
      # If most y values fall within 1% of the range (two intervals), the baseline
      # is close to flat.
      baseline_intervals = y_interval_top_prop$y_interval[1:2]

    }else {
      # If most y values fall within a few intervals, then baseline drift occurs,
      # and we use the number of intervals that cover half of the y values
      baseline_intervals = y_interval_top_prop %>%
        filter(ratio >= 0.5)  %>%
        pull(y_interval)
    }

    # Add a variable bl which indicates the baseline (1) or any other areas (0):
    dat = dat %>%
      mutate(bl = ifelse(y_interval %in% baseline_intervals, 1, 0)) %>%
      group_by(part) %>%
      mutate(n = sum(bl ==1)) %>%
      ungroup()

    ### Only baseline.
    baseline = dat %>%
      filter(part != 5) %>%
      filter(bl == 1) %>%
      filter(n > (row_count/5)*0.1)   # remove small groups of points from bl:

    if(nrow(baseline) == 0){
      # If the baseline cannot be found, the function returns NULL.
      dat = dat %>% select(x, y)
      return(dat)
    }

    # Filter out the region after this largest peak , if the maximum intensity
    # is not in the first half of the chromatogram.
    if( all(!dat$part[which(dat$x == peak_x)] %in% c(1, 2),
       length(which(baseline$x < peak_x)) > 30)){

      baseline = baseline %>%
        filter(x < peak_x)
    }

    if(nrow(baseline) == 0){
      dat = dat %>% select(x, y)
      return(dat)
    }

    # additionally filter the values at the beginning (close to the beginning
    # of the chromatogram) and at the end (close to the beginning of the peak)
    if(nrow(baseline) > 30){

      # If all baseline values fall within one interval, the baseline is flat
      # and we remove 5% before the first peak and 15% closer to the beginning,
      # since the baseline can be unstable in this part
      if(length(baseline_intervals) == 1){   #
        baseline = baseline %>% slice_head(prop = 0.95) %>%
          slice_tail(prop = 0.85)

        # If baseline drift is observed, we remove more of the values from the
        # beginning of the chromatogram to correct the baseline.
      }else if(length(baseline_intervals) == 2){   # baseline is close to flat
        baseline = baseline %>% slice_head(prop = 0.95) %>%
          slice_tail(prop = 0.8)

      }else{       # baseline drift occurs.
        baseline = baseline %>% slice_head(prop = 0.95) %>%
          slice_tail(prop = 0.65)
      }
    }

    ### Intensity correction.

    # A chromatogram consists of two main regions: the baseline and the peak.
    # At the baseline: the y values (intensity) are close to 0.
    # However, GPC chromatograms can have:
    #  - baseline shift (the chromatogram moves below or above the baseline),
    #  - baseline drift (when the baseline is not flat), or
    #  - other problems (irregularities with too high or too low intensity)

    # Add a variable ybl which is the intensity (y) after baseline correction.

    if( nrow(baseline) < 25 ){
      # Set ybl to be the difference between y and the mean (y) in the 'baseline'
      # if there are less than 25 baseline points.
      bl_mean_y = dat %>%
        filter(bl == 1) %>%
        summarise(mean(y)) %>%
        pull()

      if(nrow(bl_mean_y) == 0){
        y_min = min(dat$y)
        dat = dat %>% mutate( ybl = y - y_min)
      }else{
        dat = dat %>% mutate( ybl = y - bl_mean_y)
      }

    }else{
      # Determine whether the baseline is sloped or flat

      # Find the correlation between y and x at the 'baseline' if the sd of
      # the variable y is not 0
      if( sd(baseline$y) == 0){
        cor_bl = 0
      }else{
        cor_bl = cor(baseline$y, baseline$x , method = "spearman")
      }

      if( abs(cor_bl) < 0.85 ){
        # If correlation is less than 0.85, we set ybl as the difference
        # between y and the mean value y in the 'baseline'.
        baseline = baseline %>% summarise(mean(y)) %>% pull()
        dat = dat %>% mutate( ybl = y - baseline)

      }else{
        # Otherwise we calculate the ab line and subtract it from y.
        fit_bl = lm(y ~ x , data = baseline)$coefficients
        baseline = dat$x*fit_bl[2] + fit_bl[1]
        dat = dat %>% mutate( ybl = y - baseline)
      }
    }

    # (2.3) Find Peaks

    min_x = min(dat$x)  # start of chromatogram, time (x)
    max_x = max(dat$x)  # end of chromatogram, time (x)

    ###                   System peak
    # A system peak occurs after the main peak and may be large compared to
    # the peak of interest. The part of the chromatogram containing the system
    # peak may be removed to improve peak detection.

    # There is a limit for the system peak as 'not less than 20% of the total time (x)';
    # this limit is set for the shiny app.
    # If there is an input for the system peak time, filter rows with x greater than syst.peak:
    if(all(syst.peak != 0 , !is.na(syst.peak) , syst.peak/max_x > 0.2, syst.peak < max_x) ){
      dat = dat %>% filter(x < syst.peak)
    }

    # Intensity (y) at the peak (or maximum intensity) after the baseline correction
    peak_ybl = max(dat$ybl)
    peak_ybl = ifelse(peak_ybl <= 0, 0.000001, peak_ybl)

    # Average intensity (y)  at the baseline
    bl_mean_y = dat %>%
      filter(bl == 1) %>%
      summarise(mean(ybl)) %>%
      pull()
    bl_mean_y = ifelse(any(bl_mean_y <= 0, is.na(bl_mean_y)), 0.000001, bl_mean_y)

    # Intensity range (y) at baseline
    intensity_range = diff(range(dat$ybl))

    intensity_negative = dat %>%  filter(ybl < 0)

    if(nrow(intensity_negative) > 0){
      # Check if there is a significant negative peak
      intensity_negative  = intensity_negative %>%
        summarise(y_negative = min(ybl), width_negative = sum(ybl < min(ybl)/2))
    }else{
      intensity_negative  = data.frame(y_negative = 0, width_negative = 0)
    }

    ###        Compare peak intensity with overall intensity changes.
    # If the maximum intensity (ybl) is too low compared to the negative system
    # peak or not sufficiently above the baseline, then only a solvent
    # is present in the sample and peaks cannot be detected.
    if(any( peak_ybl/intensity_range < 0.03,
            bl_mean_y/peak_ybl > 0.7,
            intensity_range == 0,
            abs(intensity_negative$y_negative)/intensity_range > 0.95 &
            intensity_negative$width_negative > 50)){
      dat = dat %>% select(x, y, ybl)
      return(dat)
    }

    ###      Check if the chromatogram is noisy.
    # Compare y values and smoothed y values.

    if(all(syst.peak != 0 , !is.na(syst.peak) , syst.peak/max_x > 0.2, syst.peak < max_x) ){
      # Remove first 10% of chromatogram if system peak value is specified.
      ind = round(row_count*0.1):nrow(dat)
      ybl_sub = dat$ybl[ind]

    }else{
      # Remove the first and last 10% of the chromatogram that may contain abnormalities.
      ind = round(row_count*0.1):round(row_count*0.9)
      ybl_sub = dat$ybl[ind ]
    }

    y_smooth = smooth(ybl_sub, kind = "3RS3R")

    # 'ybl_sub' is vector with intensity values after baseline correction,
    # 'y_smooth' are the same values after applying the smooth() function;
    # 'noise' is the root mean square value divided by the maximum intensity
    # of the chromatogram.

    noise = mean(sqrt((y_smooth - ybl_sub)^2))/peak_ybl
    noise_top = noise
    if(length(ind) > 20){
      ind_top = dat[ind, ]$y_interval > 150
      noise_top = mean(sqrt((y_smooth[ind_top] - ybl_sub[ind_top])^2))/peak_ybl
      noise_top = ifelse(is.na(noise_top), noise, noise_top)
    }

    # If the nose is above 'noise_lim' then no feather calculations will be done
    if(noise > noise_lim  ){
      dat = dat %>% select(x, y, ybl)
      return(dat)
    }

    # If nose is above 'noise_red', the smoothed y will be used as  intensity
    # 'yi' for the calculations.
    if(noise > noise_red | noise_top > noise_red*3){
      k1 = round( mean(dat$ybl > peak_ybl*0.5)*row_count/4)
      k1 = ifelse(k1%%2 == 0, k1-1, k1)
      k1 = ifelse(k1 < 3, 3, k1)
      k1 = ifelse(k1 > 9, 9, k1)

      k2 = round( mean(dat$ybl < peak_ybl*0.2)*row_count/20)
      k2 = ifelse(k2%%2 == 0, k2-1, k2)
      k2 = ifelse(k2 < 15, 15, k2)
      k2 = ifelse(k2 > 59, 59, k2)

      bl_min_y = dat %>%
        filter(bl == 1)
      if(nrow(bl_min_y) == 0){
        bl_min_y =  -0.000001
      }else{
        bl_min_y = bl_min_y %>%
          summarise(min(y)) %>%
          pull()
      }

      # The baseline is smoothed more and the peak area less to preserve the
      # features and shape of the peaks (including negative):
      dat = dat %>%
        mutate(runmedian1 = runmed(ybl, k = k1),
               runmedian2 = runmed(ybl, k = k2)) %>%
        mutate(yi = ifelse(ybl > peak_ybl*0.05 | ybl < bl_min_y, runmedian1, runmedian2 ))

      # Additionally added loess smoothing
      s = 0.02
      suppressWarnings({
        y_loess_smooth = loess( dat$yi ~ dat$x,
                                span =  s)
      })

      if(length(y_loess_smooth$fitted)== nrow(dat)){
        dat = dat %>% mutate(yl = y_loess_smooth$fitted)

      }else{
        dat = dat %>% mutate(yl = predict(y_loess_smooth, x))
      }

      dat = dat %>%
        mutate(top = ifelse(bl == 0 & yl > (peak_ybl - bl_mean_y)*0.02, 1, 0)) %>%
        mutate(yi = yl)

    }else{
      # peaks with lower than 2% of the maximum intensity will not be detected
      dat = dat %>%
        mutate(top = ifelse(bl == 0 & ybl > (peak_ybl - bl_mean_y)*0.02, 1, 0),
               yi = ybl)
    }

    # The variable dat$top has a value of 1 for all positive peaks (above the baseline)

    ###    Find values above the baseline.
    # Add 10 values before and 10 values after each potential peak.
    top = dat %>%
      filter(top == 1 |
               top == 0 & lag(top, n = 10) == 1 |
               top == 0 & lead(top, n = 10) == 1) %>%
      filter(yi >=0)

    if(nrow(top) == 0){
      dat = dat %>% select(x, y, ybl)
      return(dat)
    }

    top = top %>%
      # again divide the y values (intensity) into 0.5% intervals
      mutate(y_interval = cut(yi, breaks = 200, labels = FALSE)) %>%
      # add a trend variable where the left side of the peak has a value of 1
      # and the right side has a value of -1
      mutate(trend = sign(y_interval - lag(y_interval)),
             dif_x_i = x - lag(x)) %>%
      fill(trend, .direction = "up") %>%
      fill(dif_x_i, .direction = "up") %>%
      mutate(trend = ifelse(trend == 0, NA, trend)) %>%
      mutate(trend = ifelse(is.na(trend) & dif_x_i > 5*dif_x, 1, trend)) %>%  #!!!!!
      fill(trend, .direction = "down") %>%
      fill(trend, .direction = "up") %>%
      mutate(trend = ifelse(dif_x_i > dif_x*5, lead(trend), trend)) %>%
      # add variable peak_id: integer, identifier for each peak and
      # peak_id = 0 for baseline and negative peaks
      mutate(peak_id = ifelse(trend == 1 & lag(trend) == -1 | dif_x_i > 20*dif_x, 1, 0)) %>%
      fill(peak_id, .direction = "up") %>%
      mutate(peak_id = cumsum(peak_id)+1)

    ### Peak maxima correspond to values with maximum intensity (y) within a single peak.
    maxs = top %>%
      group_by(peak_id) %>%
      mutate(start_y = first(y_interval)) %>%
      filter(yi == max(yi)) %>%
      mutate(n_max = row_number()) %>%
      filter(n_max == ceiling(max(n_max)/2)) %>%
      ungroup() %>%
      distinct(peak_id, .keep_all = TRUE) %>%
      mutate(type = "max") %>%
      mutate(dif_x_i = x - lag(x)) %>%
      # If the maxima of the peaks are too close (in x and in intensity),
      # we merge these peaks.
      # 'merge'  indicates whether the current peak is part of the previous one.
      mutate(merge = ifelse(dif_x_i < dif_x*10 &
                              y_interval - lag(y_interval) < 15 |
                              dif_x_i < dif_x*20 &
                              y_interval - lag(y_interval) < 10 &
                              abs(start_y - y_interval) < 10,
                            1, 0) ) %>%
      mutate(new = ifelse(merge == 1,
                            0, 1) ) %>%
      mutate(merge = ifelse(is.na(merge), 0, merge)) %>%
      # There may be more than one peak to merge, so the variable 'new' indicates
      # whether the current peak is part of the next peak and convert it to a new peak ID.
      mutate(new = ifelse(is.na(new), 1, new)) %>%
      mutate(new = cumsum(new)) %>%
      select(x, type, new, merge)

    # The top$type variable specifies the extreme points of the peak, where
    # 'start' is the beginning of the peak, 'end' is the end, 'max' is the
    # maximum intensity of the peak, all other points will be 'other'.
      top =  top %>%
        left_join(maxs, by = "x") %>%
        group_by(peak_id) %>%
        fill(new, .direction = "updown" ) %>%
        fill(merge, .direction = "updown" ) %>%
        ungroup() %>%
        mutate(peak_id = ifelse(merge == 0, peak_id, new)) %>%
        group_by(peak_id) %>%
        mutate(dif_x_i = lead(x) - x) %>% # dif_x_i was changed in 'maxs', so reset it
        fill(dif_x_i, .direction = "down") %>%
        # 'start' and 'end' variables specify whether to adjust the start and end
        # of the peak if part of the peak tail falls below the baseline
        mutate(start = ifelse(dif_x_i > dif_x*3 & trend == 1, 1, 0)) %>%
        arrange(desc(x)) %>%
        mutate(start = cumsum(start)) %>%
        arrange(x) %>%
        mutate(dif_x_i = x - lag(x)) %>%
        fill(dif_x_i, .direction = "up") %>%
        mutate(end = ifelse(dif_x_i > dif_x*3 & trend == -1, 1, 0)) %>%
        mutate(end = cumsum(end)) %>%
        ungroup() %>%
        mutate(peak_id = ifelse(start != 0 | end != 0, 0, peak_id)) %>%
        group_by(peak_id) %>%
        mutate(type = case_when(
          x == first(x) & peak_id != 0 ~ "start",
          x == last(x) & peak_id != 0 ~ "end",
          .default = type
          )) %>%
        mutate(type = ifelse(type == "max" & merge != 0, "other",type)) %>%
        ungroup()

      ###   Filter "bad" peaks

      peak_ids = top %>%
        filter(peak_id > 0) %>%
        group_by(peak_id) %>%
        # variable 'sme' indicates whether all 3 required parts of the peak are
        # present (start, end, and maximum of the peak)
        mutate(sme = sum(c("start", "end", "max") %in% type)) %>%
        filter(sme == 3) %>%
        ungroup()

      if(nrow(peak_ids) == 0){
        # if all potential peaks have been filtered, only data after baseline
        # correction is returned
        dat = dat %>% select(x, y, ybl)
        return(dat)
      }

      # Add variables that will be used to exclude "bad" peaks.
      peak_ids =  peak_ids %>%
        group_by(peak_id) %>%
        mutate(peak_width = n(),  # number of points in a single peak;
               # 'from_start' and 'from_end' indicate whether the peak's maximum
               # is too close to the start or end
               from_start = row_number(x),
               from_end = row_number(desc(x)),
               peak_height = diff(range(yi)),
               y_start = y_interval[which(type == "start")], # intensity at the start;
               y_end = y_interval[which(type == "end")],     # intensity at the end.
               # 'noise_top' shows the proportion of points that change trend
               noise_top = ifelse(top == 1 & lead(top) == 1 |
                                   top == 0 & lead(top) == 0, 0, 1)) %>%
        fill(noise_top, .direction = "down") %>%
        # 'unique_prop' shows the proportion of unique y values within a single peak,
        # ideally it should be 1 or close to 1, but if a peak is just part of a
        # baseline it may have many y values with the same intensity (y)
        mutate(unique_prop = length(unique(yi))/peak_width,
        # 'max_prop' shows the proportion of y values equal to the maximum intensity,
        # ideally there should only be one maximum, but for a wider peak there may be more.
        # If the peak is due to baseline noise, unique_prop is increased.
               max_prop = sum(yi == max(yi))/peak_width,
               noise_top = sum(noise_top)/peak_width) %>%
        ungroup() %>%
        filter(type == "max") %>%
        filter(bl == 0 )

      if(nrow(peak_ids) == 0){
        dat = dat %>% select(x, y, ybl)
        return(dat)
      }

      if(nrow(peak_ids) == 1 ){

          if(peak_ids$yi > bl_mean_y){
            peak_ids = peak_ids %>% pull(peak_id)
          }else{
            return(dat)
          }

        }else{
          #  different limits for noisy peaks:
          if(noise > noise_red*1.5 ){

            peak_ids = peak_ids %>%
              mutate(from_start = ifelse(from_start/peak_width < 0.2 & y_interval < 30 |
                                           from_start < 4 & yi/intensity_range < 0.03 , 0, 1),
                     from_end = ifelse(from_end/peak_width < 0.2 & y_interval < 30 |
                                         from_end < 4 & yi/intensity_range < 0.03 , 0, 1),
                     y_start = ifelse(abs(y_start - y_interval) < 5 & y_interval < 20, 0, 1),
                     y_end = ifelse(abs(y_end - y_interval) < 5 & y_interval < 20, 0, 1),
                     peak_height = ifelse(noise > 0.001 & peak_height/max(peak_height) > 0.1 |
                                            noise <= 0.001 & peak_height > bl_mean_y, 1, 0))
          }else{
            peak_ids = peak_ids %>%
              mutate(from_start = ifelse(from_start/peak_width < 0.15 & y_interval < 20 |
                                           from_start < 3 & yi/intensity_range < 0.03, 0, 1),
                     from_end = ifelse(from_end/peak_width < 0.15 & y_interval < 20 |
                                         from_end < 4 & yi/intensity_range < 0.03, 0, 1),
                     y_start = ifelse(y_start >= y_interval & y_interval < 10, 0, 1),
                     y_end = ifelse(y_end >= y_interval & y_interval < 10, 0, 1),
                     peak_height = ifelse(peak_height > bl_mean_y, 1, 0))
          }

          # Set peak_id to 0 for peaks that have
          # - few points or low intensity (baseline bumps),
          # - maxima close to the beginning or end of the peak,
          # - a large proportion of identical y values,
          # - low peak height
          peak_ids = peak_ids %>%
            filter(x > (max_x - min_x)*0.05) %>%
            filter(x < (max_x - min_x)*0.95 ) %>%
            filter(peak_width > 8) %>%
            filter(y_interval > 6 | y_interval > 2 & peak_width > 20) %>%
            filter(max_prop < 0.2 ) %>%
            filter(unique_prop > 0.2) %>%
            filter(from_start == 1 ) %>%
            filter(from_end == 1) %>%
            filter(y_start == 1 ) %>%
            filter(y_end == 1) %>%
            filter(peak_height == 1) %>%
            filter(noise_top < 0.2) %>%
            pull(peak_id)
        }

        top = top %>%
          filter(peak_id  %in% peak_ids)

      if(nrow(top) == 0){
        # if all potential peaks have been filtered, only data after baseline
        # correction is returned
        dat = dat %>% select(x, y, ybl)
        return(dat)
      }

        # (2.4)   Merged peaks (shoulder peaks)
        # A slowing of the intensity change (difference in y values) indicates
        # the presence of a merged peak or shoulder peak.

        # If the change in intensity is 15% of the maximum change in intensity
        # (argument shoulder.sens[1]) and occurs at a level less than 85% of
        # the peak height (argument shoulder.sens[2]) and greater than 15% of
        # the initial/final intensity (argument shoulder.sens[3]), it is defined
        # as a 'shoulder' in the 'type' variable.
        # shoulder.sens[4] is intensity above the baseline.

      shoulders = top %>%
        mutate(y_dif = abs(yi - lag(yi)))  %>%  # y_dif is the intensity change
        fill(y_dif, .direction = "up") %>%
        group_by(peak_id, trend) %>%
        # y_dif_rate is the relative rate of change of intensity;
        # only the left side of the peak is used to find the maximum y_dif
        mutate(y_dif_rate = max(y_dif)*shoulder.sens[1]) %>%
        ungroup() %>%
        mutate(y_dif_rate = ifelse(trend == -1, NA, y_dif_rate)) %>%
        fill(y_dif_rate, .direction = "down") %>%
        # 'flat' indicates whether the intensity changes slower than it should
        # for this peak.
        mutate(flat = as.numeric(y_dif < y_dif_rate)) %>%
        group_by(peak_id) %>%
        filter(yi > first(yi)*shoulder.sens[3] & yi < max(yi)*shoulder.sens[2] & trend == 1 |
                 yi > last(yi)*shoulder.sens[3] & yi < max(yi)*shoulder.sens[2] & trend == -1) %>%
        ungroup() %>%
        # do not look for the shoulders below 10% (shoulder.sens[4]) of the maximum intensity
        filter(yi > peak_ybl*shoulder.sens[4] & flat == 1) %>%
        mutate(shoulder = 1) %>%
        select(x, shoulder)

      if(nrow(shoulders) > 0){
        top =  top %>%
          left_join(shoulders, by = "x") %>%
          mutate(type = ifelse(!is.na(shoulder), "shoulder", type))
      }

      top =  top %>%
         select(x, type, peak_id)

      # (3)   Output
      # - x and y are original data;
      # - ybl is y after the baseline correction;
      # - yi: if the peak is noisy yi is smoothed ybl, otherwise it is ybl;
      # - peak_id is an integer, it is the peak id for the peak region above
      #            zero intensity, and peak_id is 0 for all other regions;
      # - type shows the peak part:
      #         'start' is the beginning ,
      #           'end' is the end,
      #           'max' is the maximum intensity of the peak,
      #      'shoulder' is the merged peak;
      #         'other'  all other points.

      dat = dat %>% left_join(top, by = "x") %>%
        replace_na(list(type = "other", peak_id = 0) )  %>%
        group_by(peak_id) %>%
        mutate(peak_id_rm = sum(c("max", "start", "end") %in% type)) %>%
        ungroup() %>%
        mutate(peak_id = ifelse(peak_id > 0 & peak_id_rm != 3, 0, peak_id))

      if(length(unique(dat$peak_id)) < 2){
        # if all potential peaks have been filtered, only data after baseline
        # correction is returned
        dat = dat %>% select(x, y, ybl)
        return(dat)
      }

      dat %>%
        select(x, y, ybl, yi, peak_id, type)
  })
}
