readGPC <-
function(filesGPC){

  # Check input
  # filesGPC must be a data frame with columns FilePath, DataStart, DataColumns
  if(all(!class(filesGPC) %in% c("tbl_df","tbl","data.frame") )){
    warning("'filesGPC' must be a dataframe with columns FilePath, DataStart and DataColumns")
    return(NULL)

  }else{
    if( any(!c("FilePath", "DataStart", "DataColumns") %in% colnames(filesGPC))){
      warning("'filesGPC' must be a dataframe with columns FilePath, DataStart and DataColumns")
      return(NULL)
    }
  }

  # Files with an invalid data format or invalid file format have a value of NA
  # in the DataStart column and will not be used further.

  filesGPC = filesGPC %>%
    select(FilePath, DataStart, DataColumns)

  # Read files.
  map(1:nrow(filesGPC), function(i){

    if(is.na(filesGPC$DataStart[i])){
      return(NULL)
    }

    ### Read file i by line.
    suppressWarnings({
      file_i =  readLines(filesGPC$FilePath[i])
    })

    ### Remove metadata.
    file_i = file_i[ filesGPC$DataStart[i]:length(file_i)]

    ### Set column names.
    # If there two numeric columns we set column names as x and y
    if(filesGPC$DataColumns[i] == 2){
      data_col_names = c("x" , "y")

    }else{
      # If there more than two numeric columns we set column names as x, y, y1, y2 ... etc
      data_col_names = c("x", "y", paste("y", 1:(filesGPC$DataColumns[i]-2), sep = "") )
    }

    ### Find a pattern to split columns.
    sep_columns = str_extract(file_i[1], "[^\\d|\\.|-]")

    ### Split columns.

    # Note:   tidyr::separate() is superseded (tidyr_1.3.1) in favor of
    # tidyr::separate_wider_delim(), which is in experimental stage.
    # separate_wider_delim() is almost twice as fast as separate()

    file_i = data.frame(xy = file_i) %>%
      separate_wider_delim(cols = xy,
                           names = data_col_names,
                           delim = sep_columns,
                           too_many = "drop",
                           too_few = "align_start")

    # If there is any non-numeric data at the end of the file, it will be ignored.
    suppressWarnings({
      file_i = file_i %>%
        as.data.frame() %>%
        mutate( across(everything(),  as.numeric ) ) %>%
        drop_na()
    })

    # If there is no data, the readGPC() function returns NULL.
    if(nrow(file_i) < 10){
      return(NULL)
    }

    file_i = file_i %>%
     select(x, y) %>%
     filter(x > 0)

   if(nrow(file_i) < 10){
     return(NULL)
   }

    # Output:
    #          x     time (or volume)
    #          y     intensity (of the first detector, if there are more than one)
    file_i
  })
}
