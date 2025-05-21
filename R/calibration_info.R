calibration_info <-
function(filesGPC, dataGPC){

  # (1) Check input ---------------------------

  # (1.1) dataGPC
  # dataGPC must be a list of data frames or a data frame to be converted
  # to a list. If it is not, calibration() returns NULL.
  if(any(class(dataGPC) %in% c("tbl_df","tbl","data.frame"))){
    dataGPC = list(as.data.frame(dataGPC))

  }else if(any(!"list"  %in% class(dataGPC) )){
    warning(" 'dataGPC' must be a list of data frames ... NULL is returned ...")

    return(NULL)
  }

  # (1.2) filesGPC
  # filesGPC must be a data frame with the same number of rows as the length
  # of the dataGPC list. It also must have columns "FileName" and "SampleName".
  # Otherwise FileName and Mp will be NA.
  n_files = length(dataGPC)

  if(any(class(filesGPC) %in% c("tbl_df","tbl","data.frame"))){

    if(nrow(filesGPC) != n_files){
      filesGPC = data.frame(FileName = NA, SampleName = NA)
      warning("the number of samples in 'filesGPC' must be the same as in 'dataGPC'")

    }else if(any(!c("FileName", "SampleName") %in% colnames(filesGPC))){

      if(all(!c("FileName", "SampleName") %in% colnames(filesGPC))){
        filesGPC = data.frame(FileName = NA, SampleName = NA)

      }else if(!"FileName" %in% colnames(filesGPC)){
        filesGPC$FileName = NA

      }else{
        filesGPC$SampleName = NA
      }
      warning("'filesGPC' must be a data frame with columns FileName, SampleName")
    }

  }else{
    filesGPC = data.frame(FileName = NA, SampleName = NA)
    warning("'filesGPC' must be a data frame with columns FileName, SampleName")
  }

  # (2) Find RT and Mp  --------------------

  if(length(dataGPC) == 1){
    # (2.1) If all calibration standards are in one file with the one chromatogram.
    # (2.1.1) find the retention time (RT) as x that corresponds to the maximum
    #         intensity (y) or shoulder of all the peaks.
    RT = dataGPC[[1]]

    if(any(!c("x", "peak_id", "peak", "yi", "type") %in% colnames(RT))){
      data.frame(FileName = filesGPC$FileName[i], Mp= NA, RT = NA)

    }else{
      RT = RT  %>%
        group_by(peak) %>%
        mutate(up = ifelse(x < x[which(type == "max" & peak_id != 0)], 1, 0)) %>%
        filter(type == "max" | type == "shoulder") %>%
        group_by(peak, up, type) %>%
        slice_min(yi, with_ties = FALSE) %>%
        pull(x)

      # (2.1.2) Data frame with one file name and Mp set to 0 for all samples
      #         (Mp needs to be added manually)
      data.frame(FileName = filesGPC$FileName,
                 Mp = rep(0, length(RT)),
                 RT = RT) %>%
        arrange(RT)
    }

  }else{

    # (2.2) If calibration samples are in multiple files.

    RT = map(1:length(dataGPC), function(i){
      rt = dataGPC[[i]]

      if(any(!c("x", "peak_id", "peak", "yi", "type") %in% colnames(rt))){
        # If the input data does not contain all required columns, the retention
        # time cannot be determined.
        data.frame(FileName = filesGPC$FileName[i], RT = NA)

      }else{
        rt = rt  %>% filter(peak_id != 0)     # data with peaks
        n_peaks = length(unique(rt$peak_id))  # the number of peaks

        if(n_peaks == 1){
          # (2.2.1) If there is one sample per file, find the retention time (RT)
          #         as x that corresponds to the maximum intensity (y) of the
          #         first peak in each calibration file.

          rt = rt %>%
            filter(type == "max") %>%
            pull(x)

          data.frame(FileName = filesGPC$FileName[i], RT = rt)

        }else{

          # (2.2.3) If the mixture of calibration standards is in one file, find
          #         the retention time (RT) as x that corresponds to the maximum
          #         intensity (y) or shoulder of the peak.
          rt = dataGPC[[i]] %>%
            group_by(peak) %>%
            mutate(up = ifelse(x < x[which(type == "max" & peak_id != 0)], 1, 0)) %>%
            filter(type == "max" | type == "shoulder") %>%
            group_by(peak, up, type) %>%
            slice_min(yi, with_ties = FALSE) %>%
            pull(x)

          data.frame(FileName = rep(filesGPC$FileName[i], length(rt)), RT = rt)
        }
      }
    }) %>% bind_rows()

    # (2.2.2) Mp, molecular weight (in Da), is extracted from the input "Mp"
    #         (in the shiny app), or "SampleName", or FileName,

    # All numeric characters in the sample name are taken as the molecular weight
    # of the calibration sample, all subsequent characters such as k, K or kDa
    # are replaced by 000.

    if(any(c("SampleName", "FileName") %in% colnames(filesGPC))){
      suppressWarnings({
        Mp = filesGPC %>%
          mutate(Mp_sn = str_extract(SampleName, "\\d+\\.?_?\\d*[Kk]?"),
                 Mp_fn = str_extract(FileName, "\\d+\\.?_?\\d*[Kk]?")) %>%
          mutate(Mp_sn = str_replace(Mp_sn, "[\\._](\\d)[Kk]", "\\100"),
                 Mp_fn = str_replace(Mp_fn, "[\\._](\\d)[Kk]", "\\100") )%>%
          mutate(Mp_sn = str_replace(Mp_sn, "[\\._](\\d{2})[Kk]", "\\10"),
                 Mp_fn = str_replace(Mp_fn, "[\\._](\\d{2})[Kk]", "\\10") )%>%
          mutate(Mp_sn = str_replace(Mp_sn, "[Kk]", "000"),
                 Mp_fn = str_replace(Mp_fn, "[Kk]", "000"))  %>%
          mutate(Mp_fn = str_replace(Mp_fn, "\\.$", "") ) %>%
          mutate( across( c(Mp_sn, Mp_fn) , as.numeric )) %>%
          mutate(Mp = ifelse( Mp_sn < 100 | is.na(Mp_sn),  Mp_fn, Mp_sn)) %>%
          mutate(Mp = ifelse( is.na(Mp), 0, Mp) ) %>%
          mutate(Mp = ifelse( Mp < 200, 0, Mp) ) %>%
          pull(Mp)
      })

      # If the number of Mp values does not match the number of RT values, then
      # all Mp values are set to 0 and must be added manually.
      if(length(Mp) != nrow(RT)){
        Mp = rep(0, nrow(RT))
      }

    }else{
      Mp = rep(0, nrow(RT))
    }

    # (2.2.3) Data frame with file names, Mps and RTs
    bind_cols(RT, Mp = Mp) %>%
      select(FileName, Mp, RT) %>%
      arrange(RT)
  }
}
