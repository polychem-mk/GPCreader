file_info <-
function(path, name = NULL) {

  # (1) Check  input

  ## (1.1) 'path' is file path(s) or path to folder with files.
  if(any(length(path) == 0,
         !is.vector(path, mode = "character")) ){
    # If 'path' is not a character vector, file_info() returns NULL with a warning
    warning("'path' must be a valid directory path or path to file(s)")
    return(NULL)
  }

  # 'is_dir' indicates whether path is a directory name or a path to file(s)
  is_dir = dir.exists(path)

  if(any(is_dir)){
    # If 'path' is a directory name,  change it to the full file names:
    file_full_name = list.files(path, recursive = TRUE, full.names = TRUE)
  }else{
    file_full_name = path
  }

  ## (1.2) If 'name' is a file name (base file name) and there are more files
  #        in 'path' ,then the output of file_info() will only contain the files
  #        listed in the 'name' vector.
  #        If 'name' is any other name (character or numeric vector), it must
  #        have the same length as the number of files in 'path', otherwise it
  #        will be replaced by the file names from 'path'

  if(any(is.null(name), !is.vector(name), length(name) > length(file_full_name))){
    name = basename(file_full_name)

  }else if(length(name) < length(file_full_name)){

    if(is.vector(name, mode = "character")){
      name_ind = which( basename(file_full_name) %in% name)

     if(length(name_ind) > 0 ){
       file_full_name = file_full_name[name_ind]
       name = basename(file_full_name)

     }else{
       name = basename(file_full_name)
     }
    }
  }

  ## (2) Create a data frame sample_info with columns:
  # SampleName     the name of the sample, extracted from the metadata at
  #                the beginning of the file;
  # DataColumns    the number of columns with numeric values;
  # DataStart      the first line containing numeric data.

  sample_info =  data.frame( )

  for(path_i in file_full_name){

    ### (2.1) Check file extension.
    file_ext = tools::file_ext(path_i) %in%
      c("txt", "arw", "xls", "xlsx", "xlsm", "xltx",  "xltm", "csv", "tsv")

    if(any(!file.exists(path_i), !file_ext)){
      # If the file does not exist or has the wrong extension, the output of
      # file_info() will only contain the FileName and FilePath for that file
      # and will show a warning message
      warn_text = paste("cannot open file '", basename(path_i), sep = "" )
      warning(warn_text)

      file_md = data.frame(  SampleName = "",
                             DataColumns = NA,
                             DataStart = NA)
    }else{
      # If the file type is correct:

      ### (2.2) Find the first numeric data line.
      # Check the first n lines and find lines that have only numeric values.
      # Typical GPC files have non numeric data at the beginning of the file
      # (about 2 - 72 lines, or more)

      suppressWarnings({
        first_data_line = readLines(path_i, n = 35 )
        first_data_line = which(!is.na(as.numeric(str_remove_all(first_data_line,
                                                                 "\\s|\\.|,|-" ))))
      })

      if(sum(!is.na(first_data_line)) <= 20){

        suppressWarnings({
          first_data_line = readLines(path_i, n = 100 )
          first_data_line = which(!is.na(as.numeric(str_remove_all(first_data_line,
                                                                   "\\s|\\.|,|-" ))))
        })
      }

      if(sum(!is.na(first_data_line))<= 20){

        suppressWarnings({
          first_data_line = readLines(path_i, n = 220 )
          first_data_line = which(!is.na(as.numeric(str_remove_all(first_data_line,
                                                                   "\\s|\\.|,|-" ))))
        })
      }

      # The GPC file must have at least 20 lines of numeric data within the first
      # 35 to 220 lines of the file.
      if(sum(!is.na(first_data_line))<= 20){
        first_data_line = NA

      }else{
        # Find the first line, which has only numeric values, and not part of
        # the metadata: the difference in the row indexes should be 1
        first_data_line_diff = diff(first_data_line)

        ind = ifelse(any(first_data_line_diff != 1) ,
                      last(which(first_data_line_diff >1 )) + 1,
                      1)
        first_data_line = ifelse(length(first_data_line) != ind, first_data_line[ind], NA)
      }

      if(is.na(first_data_line ) ){
        # If there are no numeric lines among the first 35 to 220 lines, then
        # the file format is incorrect.
          file_md = data.frame(  SampleName = "",
                                 DataColumns = NA,
                                 DataStart = NA)
      }else{
        # If there are numeric lines in the file.

        ### (2.3) Find the number of the numeric columns.
        # The chromatogram file must contain at least 2 data columns.
        data_columns = readLines(path_i)[first_data_line] %>%
          str_split("\\t") %>%  unlist() %>% as.numeric()
        data_columns = length(data_columns[!is.na(data_columns)])

        if(data_columns < 2){
          # If there are less than two numeric columns, the data format is incorrect.
          file_md = data.frame(  SampleName = "",
                                 DataColumns = data_columns,
                                 DataStart = NA)
        }else{
          # If the data format is correct.

           ### (2.4) Find Sample Name.
          suppressWarnings({
            file_i = readLines(path_i, n = first_data_line-1) %>%
              str_split("\\t")

            # Sample Name is at the line:
            sample_name_line = which(str_detect( tolower( file_i) ,
                                                 "sample\\s?name") )[1]
          })

          if(is.na(sample_name_line) ){
            # There is no Sample Name:
            file_md = data.frame(  SampleName = "",
                                   DataColumns = data_columns,
                                   DataStart = first_data_line)
          }else {
            # There is Sample Name:

            # If there are two entries in each line, then sample name is the second
            # value (metadata is in long format, two columns: name value).
            if(length(file_i[[sample_name_line]]) == 2){
              file_md = data.frame(SampleName =  str_remove_all(file_i[[sample_name_line]][2],"\""),
                                   DataColumns = data_columns,
                                   DataStart = first_data_line)
            }else{
              # Otherwise, the sample name is on the next line
              # (metadata in wide format, two lines: names in the first line,
              # values in the second line).
              sample_name_cols = which(str_detect( tolower( file_i[[sample_name_line]]) ,
                                                   "sample\\s?name") )

              file_md = data.frame(
                SampleName =  str_remove_all(file_i[[sample_name_line+1]][sample_name_cols],"\""),
                DataColumns = data_columns,
                DataStart = first_data_line)
            }
          }
        }
      }
    }

    sample_info =  bind_rows(sample_info , file_md)
  }

  # (3) Add a column with the file name
  sample_info$FileName = name

  # (4) Add a column with the file path
  sample_info$FilePath = file_full_name

  # (5) Output:
  # FileName     file name (input or from file path);
  # SampleName   sample name, from metadata;
  # DataColumns  number of numeric columns;
  # DataStart    first numeric row of data;
  # FilePath     file path.

  sample_info %>% select(FileName,  SampleName,  DataColumns, DataStart, FilePath)

}
