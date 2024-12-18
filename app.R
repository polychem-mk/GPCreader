# GPC reader 
#                     Libraries          #######

#if(!require(shiny)) install.packages("shiny")
#if(!require(tidyverse)) install.packages("tidyr")
#if(!require(stringr)) install.packages("stringr")
#if(!require(purrr)) install.packages("purrr")
#if(!require(purrr)) install.packages("dplyr")
#if(!require(purrr)) install.packages("DT")
library(shiny)
library(tidyr)
library(stringr)
library(dplyr)
library(purrr)
library(DT)

#                     Help texts         #########
help_text1 <- HTML("<h4> Calibration</h4> 
There are two ways to set up calibration. <br>
1. If you know the Intercept and Slope coefficients, you can enter these values
in the corresponding fields in the sidebar. <br>
Then check the <strong>Use these coefficients</strong> checkbox. <br>
For a linear calibration curve: <br>
&nbsp;&nbsp;&nbsp;&nbsp;∘ set the  Intercept (it is greater than 0);<br>
&nbsp;&nbsp;&nbsp;&nbsp;∘  the Slope (it is negative number);<br>
&nbsp;&nbsp;&nbsp;&nbsp;∘ leave the second and the third coefficients equal to 0. <br>

For a polynomial calibration curve:<br>
&nbsp;&nbsp;&nbsp;&nbsp;∘ switch the <strong>Calibration curve</strong> to <em>polynomial 3</em>;<br>
&nbsp;&nbsp;&nbsp;&nbsp;∘ enter Intercept  and other coefficients;<br>
&nbsp;&nbsp;&nbsp;&nbsp;∘ the coefficients must be very precise; 
do not round them for polynomial  fit.<br><br>

2.  Another way to set up calibration is to download GPS data file(s) with calibration standards. <br>
2.1. If   each calibration sample is in separate file, upload all calibration files at once 
(<span style='color:navy;'>each calibration standard has been injected into the GPC
separately</span>). <br>
&nbsp;&nbsp;&nbsp;&nbsp;This application extracts <em>Mp</em> values for each calibration sample 
from the SampleName or FileName (<strong>Table1</strong>, the second column).
If these values are incorrect or missing, double-click the appropriate cells 
in the <strong>Table1</strong>  and enter the new <em>Mp</em>(s) values.
If necessary, do the same for the peak <em>retention time</em> values 
(<strong>Table1</strong>, the third column) and
hit  <strong>Recalculate </strong> button.<br><br>

2.2.Upload a single file, if <span style='color:navy;'> a mixture of calibration 
standards has been injected into the GPC</span>. 
Add the <em>Mp</em> values to the 2nd column of <strong>Table1</strong>
(double-click the corresponding cells) and hit  <strong>Recalculate </strong> button. <br><br>

• The <em>retention time</em> values can also be added from <strong>Plot 1</strong>.
If you click on this plot, the x and y values will appear in the text box below.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;To add this RT value to the <strong>Table1</strong>, 
click the <strong>Add RT</strong> button.
<br>
• Missing values or values equal to zero
in the <strong>Table1</strong>  will be excluded from  calculations.
<br>
• You can also enter the time for the <em>system peak</em> in the sidebar. 
In this case all chromatograms will be trimmed by this value. "  )

help_text2 <- HTML("<h4>Data tab</h4>
1. Upload GPC file(s)  <br><br>
2. The calculations are summarized in the  <strong>Table 3. Results</strong>. 
Check the <strong>Add calibration information</strong> box if you want this
information to be added to the  table.
To save the results click on <strong>Save</strong> button. <br><br>
3. Plots. This tab displays the chromatogram and MWD plots for each peak.
<br><br>
If the intercept and slope coefficients are missing or set to incorrect values,
only the chromatogram is displayed and the PDIs cannot be calculated and are 
not displayed in the results table.
<br><br>
<h4>Sidebar</h4>
• You can switch between <strong>linear</strong> and <strong>polynomial 3</strong> 
 calibration; the <strong>Results</strong> table will be automatically updated. <br><br>
• To cut off the right part of the chromatogram enter the value for the <strong>System peak </strong> time.<br><br>
• If the peak is incomplete, to add the  missing tails check 
<strong>Add tails</strong> box.<br><br>
•  To shift the <em>start</em> and the <em>end</em> of the peak enter
<strong>t1</strong> and <strong>t2</strong> values (percentage of peak width); 
check the boxes below to apply the changes.          ")


###                                 1.  FUNCTIONS          #########
# (1.1) file_info #####
# Returns data frame
# Extracts:     file name    sample name               
#              <FileName>  <SampleName>              
file_info <- function( file_path, file_name ) {  # 
  
  md =   data.frame( ) 
  for(path in file_path){
    # Find the first numeric data line
    suppressWarnings({
      first_data_line = readLines(path, n = 100 ) 
      first_data_line = which(!is.na(as.numeric(str_remove_all(first_data_line, "\\s|\\.|,|-" ))))[1]
    })
    
    if(is.na(first_data_line ) ){ # If there are no numeric lines among the first 100 lines, 
      suppressWarnings({          # then we have an incorrect file format.
        num_y = readLines(path, n = 1) %>% str_split("\\t") %>% unlist() %>% length()
        file_md = data.frame(  SampleName = "no data",
                               # Date = date_col, 
                               DataColumns = num_y) 
      })
      
    }else if(first_data_line >1 & !is.na(first_data_line) ){
      suppressWarnings({    # find Sample Name line
        num_y = readLines(path)[first_data_line] %>% str_split("\\t") %>% 
          unlist() %>% length()
        file_i = readLines(path, n = first_data_line-1) %>% 
          str_split("\\t") 
        sample_name_lines = which(str_detect( tolower( file_i) , "sample\\s?name") )[1] 
      })
      
      if(is.na(sample_name_lines) ){  # if there is no Sample Name 
        file_md = data.frame(  SampleName = "no sample name",
                               # Date = date_col, 
                               DataColumns = num_y) 
      }else {      # if there is Sample Name 
        # If there are two entries in each line, then sample name is the second value
        if(length(file_i[[sample_name_lines]]) == 2){          
          file_md = data.frame(SampleName =  str_remove_all(file_i[[sample_name_lines]][2],"\""),
                               # Date = date_col, 
                               DataColumns = num_y) 
        }else{
          # Otherwise, the sample name is on the next line 
          sample_name_cols = which(str_detect( tolower( file_i[[sample_name_lines]]) ,
                                               "sample\\s?name") )
          
          file_md = data.frame(SampleName =  str_remove_all(file_i[[sample_name_lines+1]][sample_name_cols],"\""),
                               # Date = date_col, 
                               DataColumns = num_y) 
        }
      }
    }else{   # if first_data_line is 1, then there is no sample information at the beginning of the file
      file_md = data.frame(  SampleName = "no sample name",
                             #   Date = date_col, 
                             DataColumns = NA) 
    }
    
    md =  bind_rows(md , file_md)
  }  
  
  md$FileName <- file_name   # adding file name
  md %>% select(FileName,  SampleName,
                #  Date,
                DataColumns)
}

# (1.2) readGPC #######
# Reads numeric data from the file; returns list of the data frames;
# it keeps only two numeric columns (x and y)
#         <x>    retention time
#         <y>    intensity 
# If the file is not in the correct format, the output data frame will 
# contain only one row with x = 0 and y = 0
readGPC <- function(file_path, x_col = 1, y_col = 2, sep = "\\t|\\s|,", header = 100 ){
  map(file_path, function(path){     
    
    # check file extension
    file_extension = str_extract(path, "\\.[aclmrstvwx]{3,4}$") %in% 
      c(".txt", ".arw", ".xls", ".xlsx", ".xlsm", ".xltx",  ".xltm", ".csv", ".tsv")
    
    suppressWarnings({
      first_data_line = readLines(path, n = header )  # Find the first numeric data line
      first_data_line = which(!is.na(as.numeric(str_remove_all(first_data_line, "\\s|\\.|,|-" ))))[1]
      num_y = readLines(path)[first_data_line] %>% str_split(sep) %>% unlist() %>% length()
    })
    
    if(num_y <= 1 | is.na(first_data_line) | !file_extension ){     # If there is only one numeric column or 
      file_i = data.frame(x = 0, y = 0)           # no numeric data at all, the file is not in the correct format
      
    }else{               # If there two or more numeric columns.
      if(num_y == 2){    # If there two numeric columns we set column names as x and y
        data_col_names = c("x" , "y")
      }else{             # If there more than two numeric columns we set column names as x, y, y1, y2 ... etc
        data_col_names = c("x", "y", str_remove(paste("y", 1:(num_y-2) ), " ") )
      }  
      # filter out non-numeric rows, drop rows with NAs
      suppressWarnings({
        file_i = data.frame(xy = readLines(path)[-seq(1,first_data_line)]) %>%   
          separate(xy, into = data_col_names ,sep = sep) %>%  drop_na() %>%
          mutate( across(everything(),  as.numeric ) ) %>%  drop_na()
      })
      
      # If 'x' column (time) or 'y' column (intensity) are not the 1st and the 2nd columns
      if( x_col != 1 | y_col != 2){
        file_i = file_i[ , c(x_col, y_col) ]
        colnames(file_i) = c("x", "y")
      }
    }
    file_i %>% select(x, y) %>% filter(x > 0)
  })
}

# (1.3) baselineGPC #######
# makes baseline correction; returns list of the data frames;
# columns:   <x>    retention time
#            <y>    intensity 
#            <ybl> intensity after baseline correction
baselineGPC <- function(GPCdata,
                        signif_dif = 0.02,  # percentage of maximum intensity that is significant
                        # minimum proportion of approximately equal y values to be defined as baseline values                        
                        prop_eql_y = 0.25,
                        # check the data  is chromatogram                        
                        check_interval = TRUE){
  map(GPCdata, function(GPCdata_i){
    file_i = GPCdata_i  
    nr = nrow(file_i)   
    dif_x  =  mean(diff(file_i$x))
    range_dif_x  =  mean(range(diff(file_i$x)))
    diff_x_lgl = TRUE
    if(check_interval){
      diff_x_lgl = abs(dif_x - range_dif_x)/dif_x < 0.01
    }
    
    # Check data format. x suppose to have the same interval between measurements; 
    # if it is not then no calculations will be performed and the output will have only original data
    if(nr != 1 | diff_x_lgl){
      # to find y values that are not significantly different from each other, we set an argument for the round()
      #  function that depends on the range of y values (divided by 2, since negative peaks are common)
      range_y = abs( diff(range(file_i$y, rm.na = TRUE)))
      signif_dif_rnd = -round(log10(signif_dif*range_y/2))  
      if(signif_dif_rnd == nchar(max(file_i$y)) - nchar(round(max(file_i$y)))-1){
        signif_dif_rnd = signif_dif_rnd-1
      }
      # Find the proportion of these values. The baseline values are the most common in a typical chromatogram.
      file_i = file_i %>% mutate(yr = round(y, signif_dif_rnd)) %>%
        group_by(yr) %>% mutate(yr_prop = n()/nr) %>% ungroup() %>%
        mutate(bl = ifelse(yr_prop == max(yr_prop) | yr_prop > prop_eql_y, 1, 0))
      
      # If there are less than 25% approximately the same y values then there is significant baseline drift. 
      # In this case, we add more baseline values, rounding y to less significant digits.
      if(mean(file_i$bl != 0) < 0.25){
        file_i = file_i %>% mutate(yr = round(y, signif_dif_rnd-1)) %>%
          group_by(yr) %>% mutate(yr_prop = n()/nr) %>% ungroup() %>%
          mutate(bl = ifelse(yr_prop == max(yr_prop) | yr_prop >prop_eql_y, 1, 0))
      }
      
      # Filter baseline part of the chromatogram. If there are irregularities with high 
      # intensity at the beginning of the chromatogram, this part should be filtered out
      cutoff1_x = file_i %>% filter(bl == 1) %>% slice_tail(prop = 0.95) %>%
        first() %>% pull(x)
      # also filter out the area after the system peak (often negative peaks)
      cutoff2_x = file_i %>% filter(x > cutoff1_x) %>% filter(y == max(y)) %>% first() %>% pull(x)
      bl_mean_y = file_i %>% filter(bl == 1) %>% summarise(mean(y)) %>% pull()
      cutoff3_x = file_i %>% 
        filter(y < bl_mean_y - (max(file_i$y) - bl_mean_y)*0.2 ) %>% first() %>% pull(x)
      
      cutoff2_x = min(c(cutoff2_x, cutoff3_x ), na.rm = TRUE)
      
      bl_correctn = file_i %>%  filter(x > cutoff1_x) %>%
        filter(x < cutoff2_x) %>% filter(bl == 1 ) %>% select(x, y)
      
      if(nrow(bl_correctn) > 30){
        bl_correctn = bl_correctn %>% slice_head(prop = 0.95) %>%
          slice_tail(prop = 0.85)
      }
      
      if( nrow(bl_correctn) < 25 ){  # if we can not find baseline
        file_i = file_i %>% mutate( ybl = y - bl_mean_y)  
        
      }else{
        # if base line is flat we find minimum y within bl_correctn
        if( sd(bl_correctn$y) == 0){
          cor_bl = 0
        }else{
          cor_bl = cor(bl_correctn$y, bl_correctn$x , method = "spearman")  
        }
        if( abs(cor_bl) < 0.85 ){
          bl_correctn = bl_correctn %>% summarise(mean(y)) %>% pull()
          file_i = file_i %>% mutate( ybl = y - bl_correctn) 
          
        }else{                                                        
          # If the base line is not flat we calculate the baseline slope.
          fit_bl = lm(y ~ x , data = bl_correctn)$coefficients
          bl_correctn = file_i$x*fit_bl[2] + fit_bl[1]
          file_i = file_i %>% mutate( ybl = y - bl_correctn) 
        }
      } 
    }  
    file_i %>% select(x, y, bl, ybl)
  })  
}

# (1.4) findPeaksGPC #######
# adds peak IDs; returns list of the data frames;
# columns:   <x>    retention time
#            <y>    intensity 
#           <ybl>  intensity after baseline correction
#           <peak_id> peak ID 
findPeaksGPC <- function(GPCbaseline,  # a list of data frames
                         signif_dif = 0.02, # percent of the max intensity that is significant
                         # cut off for the noise reduction 
                         noise_red = 0.003,
                         noise_lim = 0.03,
                         syst_peak = 0,
                         # y difference, from max, from min:
                         shoulder_sens = c(0.3, 0.9, 1.1)  ){
  map(GPCbaseline, function(GPCbaseline_i){
    file_i = GPCbaseline_i  
    nr = nrow(file_i)   
    if(nr != 1 & "ybl" %in% colnames(file_i)){  
      # There is a limit for the system peak as 'not less than 20% of the total time (x)'
      last_x = max(file_i$x)
      # if there is input for the system peak time, filter rows with x greater than syst_peak
      if(syst_peak != 0 & !is.na(syst_peak) & syst_peak/last_x > 0.2 ){ 
        file_i = file_i %>% filter(x < syst_peak)  
      }
      
      # find maximum possible intensity (y) value at the baseline 'max_intens_bl'
      max_intens_bl = file_i %>% filter(bl == 1) 
      if(nrow(max_intens_bl) > 0){
        max_intens_bl = max_intens_bl %>% summarise(max(ybl)) %>% pull()
        max_intens_bl = ifelse(max_intens_bl<= 0, 0.000001,max_intens_bl)
      }else{
        max_intens_bl = 0.000001
      }
      # maximum intensity (y) after a baseline correction (ybl)
      max_intens = max(file_i$ybl) 
      max_intens2 = ifelse(max_intens <= 0, 0.000001, max_intens)
      
      if(max_intens <=0 | max_intens_bl/max_intens2 > 0.7){ # If nothing above 0
        file_i = file_i %>% select(x, y)
      }else{    # If there are peaks above the baseline.
        dif_x  =  mean(diff(file_i$x)) # The time difference between two measurements
        
        signif_dif_rnd = -round(log10(signif_dif*max_intens)) # parameter for the round() function
        if(signif_dif_rnd == nchar(max(file_i$y)) - nchar(round(max(file_i$y)))-1){
          signif_dif_rnd = signif_dif_rnd-1        }
        
        ###                  Check if the chromatogram is noisy.
        # Compare y values and loess values.
        y_loess = loess(ybl ~ x, data = file_i, span = 0.02, control=loess.control(surface = "direct"))
        
        if(length(y_loess$fitted)== nrow(file_i)){
          file_i = file_i %>% mutate(yl = y_loess$fitted)
        }else{
          file_i = file_i %>% mutate(yl = predict(y_loess, x))
        }
        
        ###         Find intensity values ('y') above the baseline ('top == 1').
        # We find the maximum intensity for the baseline 'bl_max_yr' and y values (yi) above 
        # will be considered as detector signal. To account for the small noise in the chromatogram  
        # we round yi values to find 'top' y values.  If the chromatogram is noisy, 
        # we use smoothed intensity to find  points above the baseline and maximums.
        noisy = mean(sqrt(y_loess$residuals^2))/max_intens2
        
        if(noisy > noise_lim){  
          # If the nose is above 'noise_lim' then no feather calculations will be done
          maxs = data.frame()
        }else{
          if(noisy > noise_red){ 
            # If nose is above 'noise_red', the smoothed y will be used as  intensity 'yi'
            file_i = file_i %>% mutate(yr = round(yl, signif_dif_rnd), yi = yl)
            bl_max_yr = file_i %>% filter(bl == 1) %>% summarise(max(yr)) %>% pull()
            file_i = file_i %>% mutate(top = ifelse(yr > bl_max_yr, 1, 0 ) )
            
          }else{ # otherwise we use intensity after baseline correction as 'yi'
            file_i = file_i %>% mutate(yr = round(ybl, signif_dif_rnd),  yi = ybl)
            bl_max_yr = file_i %>% filter(bl == 1) %>% summarise(max(yr)) %>% pull()
            file_i = file_i %>% mutate(top = ifelse(yr > bl_max_yr, 1, 0 ) )
          }
          
          # The lowest intensity above baseline.
          min_intens_top = file_i %>% filter(top == 1) %>% summarise(min(ybl)/max_intens2 ) %>% pull()
          
          if(sum(file_i$top == 1) < 20 | min_intens_top > 0.05){
            ## If minimal intensity above baseline is too high or too many 
            # values were filtered out in the previous step, we low 
            # the the cutoff to 5 % of the maximum and 15 % if the chromatogram is noisy
            if(noisy > noise_red ){
              bl_max_yr = max_intens2*0.15
            }else{
              bl_max_yr = max_intens2*0.05
            }
            # to not include noise after system peak, set right limit for the x
            x_lim_right = file_i %>% filter(top == 1) %>% last() %>% pull(x)
            file_i = file_i %>% 
              mutate(top = ifelse(yi > bl_max_yr & x < x_lim_right + dif_x*35 , 1, 0 ) )
          }
          
          # If there are too many points were detected above the baseline then there is a significant
          # drift in the chromatogram; so we set a slightly larger cut off for the 'yr' 
          if( sum(file_i$top == 1)/nr >0.5 ){ 
            bl_max_yr = file_i %>% filter(bl == 0) %>% 
              filter(yr > bl_max_yr) %>% first() %>% pull(yr)
            file_i = file_i %>% mutate(top = ifelse(yr > bl_max_yr, 1, 0 ) )
          }
          
          # find groups of values above the baseline and filter out groups with 
          # low count and low intensity (y)
          tops = file_i %>% filter(top == 1) %>% 
            mutate(top_group = ifelse( (lead(x) - x) > dif_x*2, x, NA)) %>%
            fill(top_group, .direction = "up") %>%
            mutate(top_group = ifelse(is.na(top_group), max(x), top_group)) %>%
            group_by(top_group) %>% 
            mutate(n_top = n(),
                   # to filter out bumps and baseline drift at the beginning and the end
                   max_y_top = max(yi), 
                   at_beginning = ifelse(first(x) < last_x*0.1, 1, 0 ),
                   at_end = ifelse(last(x) > last_x*0.9, 1, 0 ),
                   # the rate of y changing between two measurements  'y_dif2' will be used
                   # to find shoulder peaks
                   y_dif2 = abs(lead(yi) - yi ) ) %>% 
            fill(y_dif2) %>% mutate( y_dif2_max = max(y_dif2)) %>%
            ungroup() %>% mutate(n_top =n_top/max(n_top)) %>%
            mutate(top = ifelse(  n_top < 0.2 & max_y_top < max_intens*0.05 | 
                                    y_dif2_max == 0 |
                                    max_y_top <= 0 |
                                    n_top < 0.05 |
                                    at_beginning == 1|
                                    at_end == 1, 0, top)) %>%
            filter( top ==1) 
          
          top_col = tops %>% select(x, top)   # update 'top' in "file_i"
          file_i = file_i %>% select(-top) %>% left_join(top_col, by = "x") %>%
            mutate(top = ifelse(is.na(top), 0, top))
          
          # intensity may change in different rate for different parts of the chromatogram, 
          # especially if there is large system peak. To account for this we group
          #   'top' y values and find how intensity is changing whithin these groups
          res = tops %>%  group_by(top_group)  %>%    
            #  at list this number of points between the peak maximums and baseline
            summarise(m = max(yi)/max(y_dif2, na.rm = TRUE)  ) 
          tg = unique(tops$top_group)
          
          ###              Find maximums (peaks)
          maxs =  map( 1:length(tg), function(i){
              # 'dif_x_group' is difference in time (x) that will be used to find 
              # if there is a significant gap between groups of values
              dif_x_group = dif_x*res$m[i]/2
              if(res$m[i]/2 >15){
                dif_x_group = dif_x*res$m[i]/3
              }
              if(res$m[i]/2 < 3 ){
                dif_x_group = dif_x*3
              }
              
              maxs_i = tops %>% filter( top_group == tg[i]) %>%
                # define decreasing y values as 'down ==1' and increasing intensity as 'down ==0'
                mutate(down = ifelse(yi >= lead(yi, ceiling(res$m[i]/3) ), 1, 0 ) )%>%
                fill(down) %>%   group_by(down) %>%  # group these values
                mutate(group_down = ifelse( (lead(x) - x) > dif_x_group, x, NA)) %>%
                fill(group_down, .direction = "up") %>%
                mutate(group_down = ifelse(is.na(group_down), max(x), group_down)) %>%
                ungroup()  %>% group_by(group_down ) %>%
                # set 'type' variable, where    1 is maximum (peak)
                #                               4 is minimum
                #                               3 is shoulder (the rate of y changing is lower than within a peak slope)
                #                               2 is shoulder peak
                mutate(type = ifelse( y_dif2/max(y_dif2) < shoulder_sens[1] & 
                                        yi < max(yi)*shoulder_sens[2] & 
                                        yi > min(yi)*shoulder_sens[3], 3, 10)) %>% 
                mutate(type = ifelse(yi == max(yi),1, type) ) %>% 
                ungroup() %>% group_by(type) %>%   # remove duplicated maximums
                mutate(dif_x_type = as.numeric(lead(x) - x > dif_x_group )) %>% 
                mutate(dif_x_type = ifelse( is.na(dif_x_type ), 1, dif_x_type)  ) %>%
                ungroup() %>%  # combine maximums
                mutate(type = ifelse( type == 1 & dif_x_type ==0 , 10, type )) %>%
                mutate(between_max = ifelse(type == 1, x, NA)) %>%
                fill(between_max, .direction = "up") %>%
                mutate(between_max = ifelse(is.na(between_max), max(x), between_max)) %>%
                group_by(between_max) %>%  # find minimums between maximums (type = 4)
                mutate(type = ifelse(yi == min(yi), 4, type)) %>% 
                ungroup() %>%  # if max is close to the min then it is a part of the shoulder
                mutate(type = ifelse(type == 1 & lead(type)== 4 |
                                       type == 1 & lead(type, 2)== 4|
                                       type == 1 & lead(type, 3)== 4|
                                       type == 1 & lag(type)== 4 |
                                       type == 1 & lag(type, 2)== 4|
                                       type == 1 & lag(type, 3)== 4, 3, type))
              # find shoulder peaks
              shoulders = maxs_i %>% filter( type == 3) %>%
                mutate(dif_x_type = ifelse(lead(x) - x > dif_x_group, x, NA)) %>%
                fill(dif_x_type, .direction = "up") %>%
                mutate(dif_x_type = ifelse(is.na(dif_x_type), max(x), dif_x_type)) %>%
                group_by(dif_x_type) %>% mutate(n = n()) %>% filter(n > 1) %>%
                arrange(y_dif2) %>% mutate(type = ifelse( row_number() == 1 , 2, type)) %>%
                ungroup() %>%  select(x,  type)
              
              maxs_i = maxs_i %>% filter( type == 1) %>% select(x,  type ) %>%
                mutate(dif_x_type = 0) %>% bind_rows(shoulders ) %>%  arrange(x)   
              
            } ) %>% bind_rows()
        }   
        
        if(nrow(maxs) == 0 | !1 %in% maxs$type){   # if we filter out all the maximums
          file_i = file_i %>% select(x, y)
        }else{
          file_i = file_i %>% left_join(maxs, by = "x") %>%
            mutate(type = ifelse(is.na(type), 10, type))%>% 
            # find minimums between maximums (above or equal to 0 intensity)
            mutate(between_max = ifelse(type == 1, x, NA)) %>%
            fill(between_max, .direction = "up") %>%
            mutate(between_max = ifelse(is.na(between_max), max(x), between_max))# %>%
          # keep only first match
          divs =  file_i %>%  group_by(between_max) %>%
            mutate(divs = as.numeric(yi == min(yi)))  %>% filter( divs == 1) %>% 
            summarise(x = first(x) , divs= first(divs)) %>% select(x, divs)
          
          file_i = file_i %>% left_join(divs, by = "x") %>%
            # type == 4 is minimum between two peaks
            mutate(type = ifelse(!is.na(divs), 4, type)) %>% 
            mutate(between_min = ifelse(type == 4, x, NA)) %>%
            fill(between_min, .direction = "up") %>%
            # set to 0 the first and the last between_min groups
            mutate(between_min = ifelse(is.na(between_min), 0, between_min)) 
          
          first_between_min = file_i$between_min[1]
          file_i = file_i %>%
            mutate(between_min = ifelse(between_min == first_between_min, 0, between_min)) 
          
          ###                 Find Starts and ends of the peaks
          starts_ends = file_i %>% filter(between_min != 0) %>%
            group_by(between_min) %>%
            # find where top points start/end and add 10 more points to the left/right
            mutate(starts = ifelse(  top == 0 &  lead(top, 10) == 1 , 1, 0),
                   ends = ifelse(top == 1 & lead(top, 10) == 0 , 1, 0) )  %>%
            # also top points on the left from maximum for the start and on the right for the ends
            mutate(starts = ifelse( top == 1 & between_max != max(between_max, na.rm = TRUE) , 1,starts),
                   ends =  ifelse( top == 1 & between_max != min(between_max, na.rm = TRUE) , 1, ends))  %>%
            # if min between two peaks is above a baseline, it is start/end of the peak
            mutate(starts = ifelse( top == 1 & x == min(x, na.rm = TRUE), 1, starts),
                   ends = ifelse( top == 1 & x == max(x, na.rm = TRUE), 1, ends))  %>%
            filter(starts == 1 | ends == 1) %>% filter(yi >= 0) %>%
            filter(x  == first(x) | x == last(x)) %>% 
            ungroup() %>% select(x, starts, ends)
          
          # check if we have the same number of the starts and the ends
          s = starts_ends %>% filter( starts == 1 ) %>% pull(x)
          e = starts_ends %>% filter( ends == 1 ) %>% pull(x)
          if(length(s) != length(e)){
            if( length(s) < length(e) ){
              starts_ends2 = maxs %>% filter(type == 1) %>%
                full_join(starts_ends, by = "x") %>%
                full_join(divs, by = "x" ) %>% arrange(x) %>%
                mutate(starts = ifelse(lead(type) == 1 & is.na(starts) ,1 , starts) )
            }else{
              starts_ends2 = maxs %>% filter(type == 1) %>%
                full_join(starts_ends, by = "x") %>%
                full_join(divs, by = "x" ) %>% arrange(x) %>%
                mutate(ends = ifelse(lag(type) == 1 & is.na(ends) ,1 , ends) )
            }
            starts_ends = starts_ends2 %>% filter(starts == 1 | ends == 1) %>%
              select(x, starts, ends)
          }
          
          file_i = file_i %>% left_join(starts_ends, by = "x") %>%
            mutate(starts = ifelse(is.na(starts), 0, starts),
                   ends = ifelse(is.na(ends), 0, ends)) %>%
            mutate(peak_id = cumsum(starts+ ends))  %>%
            mutate(peak_id = ifelse(ends == 1, lag(peak_id), peak_id)) %>%
            group_by(peak_id) %>% 
            mutate(n_id = n(),
                   max_y_prop = max(yi)/max_intens2,
                   maxs = as.numeric(min(type) == 1) ) %>% ungroup() %>%
            mutate(peak_id = ifelse(n_id < 8 |
                                      max_y_prop < 0.03 |
                                      max(starts, na.rm = TRUE) == 0 |
                                      maxs == 0,   0, peak_id)) %>%
            ungroup() %>%
            mutate(type = ifelse(peak_id == 0, 10, type))  %>%
            mutate(type = ifelse(type ==10  , 0, type) ) %>%
            select(x, y, ybl, yi, type, peak_id)
        } 
      }
    }
    file_i
  }) 
} 

# (1.5) filterPeaksGPC #######
# adds peak IDs; returns list of the data frames
filterPeaksGPC <- function(GPCpeakID,  first_peak = TRUE){
  map(GPCpeakID, function(GPCpeakID_i){
    file_i = GPCpeakID_i  
    if(ncol(file_i)>2 & nrow(file_i) >1 ){
      peak_IDs = data.frame(  table(file_i$peak_id)   )  
      colnames(peak_IDs) = c("peak_ID", "n")
      peak_IDs = peak_IDs %>% 
        mutate(peak_ID = as.numeric( as.character(peak_ID) )) %>%
        filter(peak_ID != 0) %>%  filter(n > 8)
      
      # If we filter out all the peaks, then only the original x and y will remain at the output.
      if(nrow(peak_IDs) < 1 ){
        file_i = file_i %>% select(x, y)
      }else{
        dif_x  =  mean(diff(file_i$x)) # The time difference between two measurements
        
        peaks_stack = data.frame()
        for(i in 1:nrow(peak_IDs)){
          time_lims = file_i %>% filter( peak_id == peak_IDs$peak_ID[i]) %>%
            filter(x == min(x) | x == max(x)  | yi == max(yi)) %>%
            mutate(prop_y = yi/max(yi) ) %>% filter(x == first(x) | x == last(x))
          
          time_lim1 = time_lims$x[1]- peak_IDs$n[i]*dif_x*time_lims$prop_y[1]*2.5
          
          # the right tail is often longer than the left one
          time_lim2 = time_lims$x[2] + peak_IDs$n[i]*dif_x*time_lims$prop_y[2]*2.5
          
          peak_i = file_i %>% filter( x > time_lim1 & x < time_lim2) %>%
            mutate(peak_id = ifelse(peak_id == peak_IDs$peak_ID[i], peak_id, 0) ) %>%
            mutate(peak = i) %>% mutate(peak_id = ifelse(peak_id != 0, i, peak_id))
          
          peaks_stack = rbind(peaks_stack, peak_i)
        }
      }
      
    }else{
      peaks_stack = file_i %>% select(x, y)
    }
    peaks_stack 
  })
}


# (1.6) fitPeaksGPC   ##########
fitPeaksGPC <- function(GPCfilterPeaks,
                        fit_tails = 0.05,
                        calc_coef = 0.45,
                        add_ab = 0.1,
                        err_min = 0.007){
  map(GPCfilterPeaks, function(GPCpeakID_i){  
    file_i = GPCpeakID_i            
    if(ncol(file_i)>2 & nrow(file_i) >1 ){
      
      peaks_summ = file_i %>% filter(peak_id != 0) 
      should = file_i %>% filter(type %in% c( 2,3)) 
      
      if(nrow(should) != 0){
        should = peaks_summ %>% group_by(peak_id) %>%
          arrange(desc(yi)) %>% mutate(max_x = first(x)) %>%
          filter(type %in% c( 2,3)) %>%
          mutate(is_left = x < max_x) %>% group_by(peak_id, is_left) %>%
          filter(row_number() == 1) %>%
          
          mutate(start_y = ifelse( is_left, yi, NA),
                 end_y = ifelse( !is_left, yi, NA) ,
                 start_x =ifelse( is_left, x, NA),
                 end_x =ifelse( !is_left, x, NA) ) %>%
          group_by(peak_id) %>% fill(start_y , .direction = "updown" ) %>%
          fill(end_y , .direction = "updown" ) %>%  fill(start_x , .direction = "updown" ) %>%
          fill(end_x , .direction = "updown" ) %>%  filter(row_number() == 1) %>%
          ungroup() %>% arrange(x) %>%  select(peak_id, peak,start_x,  start_y, end_x, end_y)
      }
      
      peaks_summ = peaks_summ %>%
        group_by(peak_id) %>%
        summarise( peak = first(peak),
                   start_x = first(x),
                   start_y = first(yi),
                   end_x = last(x),
                   end_y = last(yi)) %>%
        mutate(start_y = ifelse(start_y <= 0, 0.0001, start_y),
               end_y = ifelse(end_y <= 0, 0.0001, end_y)) 
      
      if(nrow(should) != 0){
        suppressMessages({
          peaks_summ = peaks_summ %>% full_join(should) %>% group_by(peak_id) %>%
            summarise( peak = first(peak),
                       start_x = max(start_x, na.rm = TRUE),
                       start_y = max(start_y, na.rm = TRUE),
                       end_x = min(end_x, na.rm = TRUE),
                       end_y = max(end_y, na.rm = TRUE)) %>%
            mutate(start_y = ifelse(start_y <= 0, 0.0001, start_y),
                   end_y = ifelse(end_y <= 0, 0.0001, end_y)) 
        })
      }
      
      dif_x  =  mean(diff(file_i$x))
      
      peaks_stack = data.frame()  
      for(i in 1:nrow(peaks_summ)){ 
        #                        "df_fitted" is a data frame for each peak.
        ###  Change 'peak_id' if we cut off shoulders
        df_fitted = file_i %>% filter(peak == peaks_summ$peak[i]) %>%
          mutate(peak_id = ifelse(x < peaks_summ$start_x[i] |
                                    x > peaks_summ$end_x[i]  ,0, peak_id))
        ###  Define maximum intensity and the corresponding time
        max_y = df_fitted %>% filter(peak_id == peaks_summ$peak_id[i]) %>% 
          summarise(max(yi)) %>% pull()
        max_x = df_fitted %>% filter(peak_id == peaks_summ$peak_id[i]) %>% 
          filter(yi == max(yi) ) %>% first() %>% pull(x)
        ###  Check intensity at the left and right tails
        left = peaks_summ$start_y[i]/max_y
        right = peaks_summ$end_y[i]/max_y
        
        if( left < fit_tails & right < fit_tails){ 
          # If both tails go to the baseline and no shoulders, we do not need any fit
          df_fitted = df_fitted %>% mutate( fit_type = "none")
        }else{ # If we need to fit peak   
          
          ###                   Find slopes  
          df_fitted = df_fitted %>%
            mutate(left_tail = as.numeric(x < max_x), # Find local slopes   (between two points) 
                   slope_loc = ifelse(peak_id != 0 , (lead(yi) - yi)/(lead(x) - x), NA ) ) 
          
          # left tail supposed to have positive slopes and the right tail negative ones. The  maximum of 
          # the 'local slope' for the left tail and minimum for the right are at the half way to the maximum intensity
          # however, if chromatogram is noisy then there is no relationship between yi and 'int_loc' or 'int_slope'
          niosy = df_fitted %>% filter( peak_id != 0 ) %>% 
            slice_head(prop = 0.95) %>% slice_tail(prop = 0.95) %>%
            mutate(nois = ifelse(left_tail == 1 & slope_loc <= 0 |
                                   left_tail == 0 & slope_loc >= 0  , 1, 0)) %>%
            summarise(n = sum(nois), avg = mean(nois)) 
          
          if( is.na(niosy$avg)){
            niosy = data.frame( n = 100 , avg =1)
          }
          
          if(niosy$avg < 0.2){   
            #                                 Find slopes
            if(niosy$n < 3){   # Not noisy chromatogram
              # If there is more than 70% of the tail. 
              # - the slope will be the one at the minimum of the 'slope_loc' for 
              #   the right tail and maximum for the left tail.
              slopes = df_fitted %>% filter( !is.na(slope_loc)) %>% 
                group_by(left_tail) %>% arrange(abs(slope_loc)) %>%
                summarise( slope_loc = last(slope_loc)) %>%  arrange(left_tail) %>%
                summarise( slope_l = last(slope_loc),        
                           slope_r = first(slope_loc)  )  %>% mutate(  ab_l = 0, ab_r =0) 
              
              if(left > add_ab |  right > add_ab){  
                # If  more than 15% (add_ab) of the tail is missing.
                # - we add abline to the missing tail
                # - slope and intercept are calculated if left or right tails are greater than 45% (calc_coef). 
                
                df_fitted = df_fitted %>%
                  mutate(int_loc = ifelse(peak_id != 0 ,lead(yi) - lead(x)*slope_loc, NA))
                
                slopes = df_fitted %>% filter( !is.na(slope_loc)) %>% 
                  group_by(left_tail) %>% arrange(abs(slope_loc)) %>%
                  summarise( slope_loc = last(slope_loc), int_loc = last(int_loc), x = last(x)) %>%  # 
                  arrange(left_tail) %>%
                  summarise( slope_l = last(slope_loc), int_l = last(int_loc), x_l = last(x),        
                             slope_r = first(slope_loc) , int_r = first(int_loc) ,x_r = first(x)) %>%
                  mutate(ab_l = ifelse(left < add_ab, 0, 1 ),
                         ab_r = ifelse(right < add_ab, 0, 1 ))
                
                # add abline
                df_fitted = df_fitted %>%
                  mutate(y_full = case_when(
                    x < slopes$x_l & slopes$ab_l == 1 ~ x*slopes$slope_l + slopes$int_l,
                    x > slopes$x_r & slopes$ab_r == 1 ~ x*slopes$slope_r + slopes$int_r,
                    .default = yi) ) %>%
                  mutate(peak_id = ifelse(
                    x <= slopes$x_l & slopes$ab_l == 1 & y_full >0  |
                      x >= slopes$x_r & slopes$ab_r == 1 & y_full >0 , peaks_summ$peak_id[i], peak_id) )
                
              }else{  # if we do not need to add abline (not noisy chromatogram)
                df_fitted = df_fitted %>% mutate(y_full = yi)
              }   
              
            }else{ # if chromatogram is noisy we use all available values to find slopes with lm()
              if(left > add_ab |  right > add_ab){  # if we need to add abline
                left_half = df_fitted %>%filter( peak_id != 0) %>% filter(left_tail == 1 )
                slopes_l = lm( yi~ x, data = left_half)$coefficients 
                
                right_half = df_fitted %>%filter( peak_id != 0) %>% filter(left_tail == 0 )
                slopes_r = lm( yi~ x, data = right_half)$coefficients 
                
                slopes = data.frame( slope_l = slopes_l[2], int_l = slopes_l[1], x_l = first(left_half$x),
                                     slope_r = slopes_r[2], int_r = slopes_r[1], x_r = last(right_half$x)) %>%
                  mutate(int_l = ifelse(left < add_ab, NA, int_l ),
                         int_r = ifelse(right < add_ab, NA, int_r ))
                
                df_fitted = df_fitted %>%   # add abline
                  mutate(y_full = case_when(
                    x < slopes$x_l & !is.na(slopes$int_l) ~ x*slopes$slope_l + slopes$int_l,
                    x > slopes$x_r & !is.na(slopes$int_r) ~ x*slopes$slope_r + slopes$int_r,
                    .default = yi) ) %>%
                  mutate(peak_id = ifelse(
                    x <= slopes$x_l & !is.na(slopes$int_l) & y_full >0  |
                      x >= slopes$x_r & !is.na(slopes$int_r) & y_full >0 , peaks_summ$peak_id[i], peak_id) )
              
              }else{   # if we do not need to add abline
                left_half = df_fitted %>%filter( peak_id != 0) %>% filter(left_tail == 1 )
                slopes_l = lm( yi~ x, data = left_half)$coefficients 
                
                right_half = df_fitted %>%filter( peak_id != 0) %>% filter(left_tail == 0 )
                slopes_r = lm( yi~ x, data = right_half)$coefficients 
                
                slopes = data.frame( slope_l = slopes_l[2], 
                                     slope_r = slopes_r[2])   
                df_fitted = df_fitted %>%     mutate(y_full = yi)
              } 
            }
            
            #                          x at the half intensity
            half_x = df_fitted %>% filter(peak_id != 0 & y_full > max_y*0.5) %>%
              group_by(left_tail) %>% arrange(y_full) %>% filter(x == first(x)  ) %>% 
              ungroup() %>% arrange(left_tail) %>% 
              summarise(half_x_l = last(x), half_x_r = first(x))
            # __________________________________
            #                      set x, y,  parameters   and function
            x = df_fitted %>% filter(peak_id == peaks_summ$peak_id[i])  %>% pull(x)
            y = df_fitted %>% filter(peak_id == peaks_summ$peak_id[i])  %>% pull(y_full)
            
            pars = c(max_y,             #  set parameters:    top
                     half_x$half_x_l,   #  half x left tial 
                     slopes$slope_l,    #  left slope
                     half_x$half_x_r,   #  half x right tial 
                     slopes$slope_r  )  #  right slope
            
            # max_y/(1+10^(slope_l*(half_l - x) ) ) - max_y/(1+10^(slope_r*( x - half_r) ) )
            fit_bell_func = function( par)   {  # function: calculate sum of squares for the left part of the peak
              fitted = par[1]/(1+exp(par[3]*(par[2] - x) ) ) - par[1]/(1+exp(par[5]*( x - par[4]) ) )
              sum((y - fitted)^2)
            }
            
            try_bell =   tryCatch(     # Try to find coefficients
              expr = {
                bell_coeffs = optim(pars, fit_bell_func, method="BFGS", control=list(reltol=1e-9))
              },
              error =  function(e){
                return(NA)
              }
            )
            
            if(any(!is.na(try_bell))   ){  # estimate RMSE
              err_bell = sqrt(bell_coeffs$value)/length(y)/max_y
            }else{
              err_bell = 1
            }
            
            #                                        Add fitted values
            if(any(!is.na(try_bell))  & err_bell < err_min ){  # 
              
              if(niosy$n < 3){
                df_fitted = df_fitted %>% mutate( 
                  yf = bell_coeffs$par[1]/(1+exp(bell_coeffs$par[3]*(bell_coeffs$par[2] - x) ) ) -
                    bell_coeffs$par[1]/(1+exp(bell_coeffs$par[5]*( x - bell_coeffs$par[4]) ) ),
                  fit_type = "bell", dif_y = abs(y_full -  yf))
                
                # find x where to add fitted values
                if(left > fit_tails ){                 #  left tial
                  x_mix_l = df_fitted %>% filter(peak_id == peaks_summ$peak_id[i]) %>%
                    filter(x <= max_x) %>%  filter(y_full < max_y*0.9)
                  
                  if(nrow(x_mix_l) <3 ){
                    x_mix_l = df_fitted %>% filter(peak_id == peaks_summ$peak_id[i]) %>%
                      filter(x <= max_x) %>% filter(dif_y == min(dif_y, na.rm = TRUE)) %>% 
                      last() %>% pull(x)
                  }else{
                    x_mix_l = x_mix_l %>% filter(dif_y == min(dif_y, na.rm = TRUE)) %>%
                      last() %>% pull(x)
                  }
                }else{
                  x_mix_l = df_fitted %>% first() %>% pull(x)
                }
                #                           right tail
                if(right > fit_tails ){
                  x_mix_r = df_fitted %>% filter(peak_id == peaks_summ$peak_id[i]) %>%
                    filter(x > max_x) %>% filter(y_full < max_y*0.9)
                  
                  if(nrow(x_mix_r) <3){
                    x_mix_r = df_fitted %>% filter(peak_id == peaks_summ$peak_id[i]) %>%
                      filter(x > max_x) %>%
                      filter(dif_y == min(dif_y, na.rm = TRUE)) %>%  first() %>% pull(x)
                  }else{
                    x_mix_r = x_mix_r  %>%
                      filter(dif_y == min(dif_y, na.rm = TRUE)) %>%  first() %>% pull(x)
                  }
                }else{
                  x_mix_r = df_fitted %>% last() %>% pull()
                }
                df_fitted = df_fitted %>%
                  mutate(yi = ifelse(x < x_mix_l ,  yf, y_full )) %>%
                  mutate(peak_id = ifelse(x < x_mix_l & left > fit_tails & 
                                            yi > max_y*0.01, peaks_summ$peak_id[i], peak_id)) %>%
                  mutate(yi = ifelse(x > x_mix_r ,  yf, yi )) %>%
                  mutate(peak_id = ifelse(x > x_mix_r & right > fit_tails & 
                                            yi > max_y*0.01, peaks_summ$peak_id[i], peak_id)) %>%
                  mutate(fit_type = "bell") 
              }else{
                df_fitted = df_fitted %>%
                  mutate(yf = bell_coeffs$par[1]/(1+exp(bell_coeffs$par[3]*(bell_coeffs$par[2] - x) ) ) -
                           bell_coeffs$par[1]/(1+exp(bell_coeffs$par[5]*( x - bell_coeffs$par[4]) ) )) %>%
                  mutate(yi = yf) %>%
                  mutate( peak_id = ifelse(yi > max_y*0.01, peaks_summ$peak_id[i], peak_id))  %>%
                  mutate(fit_type = "nois_red") 
              }
              
            }else{ # if we fail to fit bell curve try to fit tail 
              df_fitted = df_fitted %>% 
                mutate( yf = NA, fit_type = "bell fail" )
              
              #  fit left tail/the Boltzmann function
              if(left > fit_tails ){    
                x = df_fitted %>% filter(peak_id == peaks_summ$peak_id[i]) %>%
                  filter(x <= max_x) %>% pull(x)
                y = df_fitted %>% filter(peak_id == peaks_summ$peak_id[i]) %>%
                  filter(x <= max_x) %>% pull(y_full) 
                
                pars = c(max_y,    # set parameters:     top
                         half_x$half_x_l                #  half_x
                         , abs(slopes$slope_l*4/max_y)  #  slope 
                )  
                
                fit_left_func = function( par)   {  # function: calculate sum of squares for the left part of the peak
                  fitted = par[1]  - par[1] / (1 + exp(( par[2] - x  )/par[3]))
                  sum((y - fitted)^2)
                }
                try_left =   tryCatch(     # Try to find coefficients
                  expr = {
                    left_coeffs = optim(pars, fit_left_func , method="BFGS", control=list(reltol=1e-9))
                  },
                  error =  function(e){
                    return(NA)
                  }
                )
                
                if(any(!is.na(try_left))   ){  # estimate RMSE
                  err_left = sqrt(left_coeffs$value)/length(y)/max_y
                }else{
                  err_left = 1
                }
                
                if(any(!is.na(try_left)) & err_left < 0.02 ){# Add fitted values
                  df_fitted = df_fitted %>% 
                    mutate(yf = ifelse(x <= max_x,
                                       left_coeffs$par[1] - 
                                         left_coeffs$par[1] / 
                                         (1 + exp((left_coeffs$par[2] - x)/left_coeffs$par[3])  ), yf ) ) %>%
                    mutate( fit_type = ifelse(x <= max_x, "tail", fit_type),
                            dif_y = abs(y_full -  yf))
                  x_mix_l = df_fitted %>% filter(peak_id == peaks_summ$peak_id[i]) %>%
                    filter(x <= max_x) %>%
                    filter(dif_y == min(dif_y, na.rm = TRUE)) %>%  last() %>% pull(x)
                  df_fitted = df_fitted %>%                               
                    mutate(yi = ifelse(x <= x_mix_l & !is.na(yf),  yf, y_full )) %>%
                    mutate(peak_id = case_when(x < x_mix_l  & yi > 0 ~ peaks_summ$peak_id[i],
                                               x < x_mix_l  & yi <= 0 ~ 0,
                                               .default = peak_id))
                }else{
                  df_fitted = df_fitted %>%
                    mutate( fit_type = "fit fail")
                }
              }
              #    fit right tail/the Boltzmann function 
              if(right > fit_tails ){
                pars = c(max_y,      # set parameters:     top
                         half_x$half_x_r                #  half_x                         
                         , abs(slopes$slope_r*4/max_y)  #  slope 
                )                     
                x = df_fitted %>% filter(peak_id == peaks_summ$peak_id[i]) %>%
                  filter(x > max_x) %>% pull(x)
                y = df_fitted %>% filter(peak_id == peaks_summ$peak_id[i]) %>%
                  filter(x > max_x) %>% pull(y_full) 
                
                fit_right_func = function( par)   {  # function: calculate sum of squares for the right part of the peak
                  fitted = par[1]  - par[1] / (1 + exp(( par[2] - x  )/par[3]))
                  sum((y - fitted)^2)
                }
                try_right =   tryCatch(     # Try to find coefficients
                  expr = {
                    right_coeffs = optim(pars, fit_right_func, method="BFGS", control=list(reltol=1e-9))
                  },
                  error =  function(e){
                    return(NA)
                  }
                )
                
                if(any(!is.na(try_right))   ){  # estimate RMSE
                  err_right = sqrt(right_coeffs$value)/length(y)/max_y
                }else{
                  err_right = 1
                }
                
                # Add fitted values
                if(any(!is.na(try_right)) & err_right < 0.02){
                  df_fitted = df_fitted %>% 
                    mutate(yf = ifelse(x > max_x , 
                                       right_coeffs$par[1] - 
                                         right_coeffs$par[1] / 
                                         (1 + exp((right_coeffs$par[2] - x)/right_coeffs$par[3])  ), yf) ) %>%
                    mutate( fit_type = ifelse(x > max_x, "tail", fit_type),
                            dif_y = abs(y_full -  yf))
                  
                  x_mix_r = df_fitted %>% filter(peak_id == peaks_summ$peak_id[i]) %>%
                    filter(x > max_x) %>%
                    filter(dif_y == min(dif_y, na.rm = TRUE)) %>% first() %>% pull(x)
                  df_fitted = df_fitted %>%   
                    mutate(yi = ifelse(x >= x_mix_r & !is.na(yf),  yf, y_full )) %>%
                    mutate(peak_id = case_when(x > x_mix_r  & yi > 0 ~ peaks_summ$peak_id[i],
                                               x > x_mix_r  & yi <= 0 ~ 0,
                                               .default = peak_id))
                }else{
                  df_fitted = df_fitted %>%
                    mutate( fit_type = "fit fail")
                }
              } 
            }
          }else{
            df_fitted = df_fitted %>% 
              mutate( fit_type = "noise only") 
          }
        }  
        df_fitted = df_fitted %>% select(x, yi, peak_id, peak,  fit_type)
        
        peaks_stack = rbind(peaks_stack, df_fitted)
      }  
      file_i = peaks_stack
    } 
    file_i
  })  
}

# _________________________  UI _____________________________#######
ui <- fluidPage(   
  h4(" GPC reader "),    # Header
  sidebarLayout(         # Sidebar layout                                                    
    ######   sidebar Panel #################################
    sidebarPanel( width = 2,
                  tags$style(".well {background-color: #FAFAFA;}"),
                  radioButtons("calib_fit_order", "Calibration:", selected = "1",
                               choices = c("linear" = "1", "polynomial 3" = "3")),
              div(HTML('<hr style="border-color: #5d6d7e;">')),
                  numericInput("intercept", "Enter the Intercept", value =  10.99),
                  numericInput("coef1",  "and the Slope", value = -0.65),
                  numericInput("coef2", 
                               "for a polynomial calibration curve, enter the 2nd ",
                               value = 0),
                  numericInput("coef3", 
                               "and 3rd coefficients ", value = 0),
                  checkboxInput("use_input_coefs", "use these coefficients", 
                                value = FALSE),
                  span(textOutput("calib_help1"), style="color:firebrick"),
              div(HTML('<hr style="border-color: #5d6d7e;">')),
                  numericInput("syst_peak",  "System peak",
                               value = 0, min = 0, step = 0.1),
              div(HTML('<hr style="border-color: #5d6d7e;">')),
                  checkboxInput("fit_missing_tail", strong("Add tails"), 
                                value = FALSE),
              div(HTML('<hr style="border-color: #5d6d7e;">')),
                  paste("Change the peak boundaries by t1 and t2 values:"),
                  numericInput("t1", "t1", value = 1,  min = -19,  max = 35), 
                  checkboxInput("use_t1", "start time", value = FALSE),
                  numericInput("t2", "t2", value = 1, min = -19, max = 35), 
                  checkboxInput("use_t2", "end time", value = FALSE)
              ), 
    ################################## Main panel ###################################
    mainPanel( width = 10, tags$style(".shiny-file-input-progress {display: none}"), 
               tabsetPanel(  
                 ##  Calibration TAB ######################################
                 tabPanel( 'Calibration',                        
                           div(HTML("<br>")),
                           fluidRow(                              
                             column(width = 6,
                                    fileInput("calib_file", "Select calibration file(s)",
                                              multiple = TRUE,
                                              accept = c(".txt", ".arw", ".xls", ".xlsx", ".xlsm", ".xltx",  ".xltm", ".csv", ".tsv"))
                             ),
                             column(width = 2,
                                    actionButton("zoom_out", "Zoom Out", width = '100%')
                             ),
                             column(width = 2),
                             column(width = 2,
                                    actionButton("help1", "Help", 
                                                 style="color: black; background-color: #EAEBEC; border-color: darkgray",
                                                 width = '100%') )
                           ), 
                           fluidRow(                
                             column(width = 6,
                                    strong(em("Table 1. Calibration samples") ),
                                    DTOutput('calib_samples_table')
                             ),
                             column(width = 6,                           
                                    strong(em("Plot 1. GPC chromatogram of standards")),
                                    plotOutput('calibration_files',
                                               click = "plot_click",
                                               brush = brushOpts("plot_brush", 
                                                                 fill = "white",
                                                                 stroke = "white",
                                                                 opacity = 0,
                                                                 resetOnNew = TRUE,
                                                                 clip = FALSE) ),  
                             )
                           ),
                           fluidRow(         
                             column(width = 6,
                                    span(textOutput("calib_help2"), style="color: #565374"),
                                    div(HTML("<br>")),
                                    actionButton("recalcutate_calib", "Recalculate",
                                                 style="color: black; background-color: #DCF0FC; border-color: gray")
                             ),
                             column(width = 6,
                                    column(width = 8,
                                           verbatimTextOutput("xy_calib", placeholder = TRUE)
                                    ),
                                    column(width = 4,
                                           actionButton("add_rt_from_click", "Add RT", width = '100%')
                                    )
                             ),
                           ),              
                           fluidRow(      
                             column(width = 6,
                                    div(HTML("<br>")),
                                    strong(em("Table 2. Calibration coefficients")),
                                    tableOutput("calib_coef")
                             ),
                             column(width = 6,
                                    div(HTML("<br>")),
                                    strong(em("Plot 2. GPC calibration curve")),
                                    plotOutput('calibration_plot', height = "380px")
                             )
                           )          
                 ),         
                 ##  DATA TAB #################################
                 tabPanel( 'Data',                               # 2ndst TAB
                           div(HTML("<br>")),
                           fluidRow(
                             column(width = 10,
                                    fileInput("gpc_file", "Select GPC file(s)", 
                                              accept = c(".txt", ".arw", ".xls", ".xlsx", ".xlsm", ".xltx",  ".xltm", ".csv", ".tsv"),
                                              multiple = TRUE)),
                             column(width = 2, 
                                    actionButton("help2", "Help", 
                                                 style="color: black; background-color: #EAEBEC; border-color: darkgray",
                                                 width = '100%'))  
                           ),
                           div(HTML("<br>")),
                       strong(em("Table 3. Results") ),
                       DT::DTOutput('data_table'), 
                           fluidRow(
                             column(width = 4,
                                    checkboxInput("add_calib_info",
                                                  "Add calibration information",
                                                  value = FALSE)
                             ),
                             column(width = 6,
                                    downloadButton("save_results", "Save",
                                                   style="color: black; background-color: #F0F8FF; border-color: darkgray",
                                                   width = '100%')
                             )
                           ),
                           div(HTML('<hr style="border-color: #CBD1E6;">')),
                       plotOutput('plot1')
                 )
               )
    )
  )                                                       
)

#_________________________  server  ________________________###########
server <- function(input, output, session) {
  
  # (1.7) calculations #####
  # adds columns Mi, dWi_over_dVi, dVi_over_dLogMi  and dWi_over_dLogMi
  # returns list of the data frames.
  MWDcalculateGPC <- function(GPCpeaks, 
                              intercept ,coef1, coef2, coef3 
  ){
    map(GPCpeaks, function(table_i) {  
      nr = nrow(table_i)
      
      if(nr <= 1 | !"yi" %in%  colnames(table_i)  ){   
        # If there is no data, the third graph will be empty and no calculations will be performed.
        table_i =  data.frame(x = NA, yi = NA)
      }else{       
        
        dif_x  =  mean(diff(table_i$x))  # time difference between two measurements
        
        #  Calculations.
        table_i = table_i %>% mutate(yi_above0 = ifelse(yi < 0, 0, yi ) ) %>%
          group_by(peak) %>%
          mutate(
            # Mi is molecular weight at the i-th point:
            Mi = 10^(coef3*(x^3)+ coef2*(x^2)+ coef1*x+ intercept ) ,
            # Values for the molecular weight distribution:
            dWi_over_dVi = yi/sum(yi_above0)/dif_x,
            dVi_over_dLogMi = abs(1/ (coef1 + coef2*2*x + coef3*3*x^2) ),
            dWi_over_dLogMi = dWi_over_dVi*dVi_over_dLogMi
          ) %>% ungroup()
      }      
    } )
  }
  
  
  
  #                   2.  REACTIVES for the Calibration files     #####
  ### For Table 1.
  # (2.1) Table 1. Calibration samples. Data table with columns: FileName, Mp, retention_time
  calib_DT <- reactiveValues(data =  {
    data.frame(FileName = character(0), Mp = numeric(0), RT = numeric(0))
  })
  
  # (2.2) Calibration data files: read, adjust baseline, and add peak IDs
  calib_data <- reactive( {  
    s_peak <- input$syst_peak 
    read_calib_files = readGPC(file_path = input$calib_file$datapath)
    baseline_calib = baselineGPC(read_calib_files)
    findPeaksGPC(baseline_calib, syst_peak = s_peak )
  })  
  
  # (2.3) Peak retention time values for calibration samples 
  calib_peak_time <- reactiveValues(data = {
    numeric(0)
  })
  
  calib_peak_time_file <- reactive({   
    data_list = calib_data()
    
    if(length(data_list) == 1){  
      # If all calibration samples are in one file, find peaks and shoulder peaks
      data_list = data_list[[1]]  
      if( "type" %in% colnames(data_list)){
        calib_peak_xy = data_list %>% filter(type %in% c(1,2)) %>% 
          group_by(peak_id) %>%
          mutate(type = ifelse(type == 2 & yi/max(yi) < 0.5, 0, type)) %>% 
          ungroup() %>%   filter( type %in% c(1,2)) %>% 
          mutate(xMax = max(x), yMin = min(yi)) %>% select(x, yi, xMax, yMin)
      }else if( all(c("x", "y") %in% colnames(data_list) )) {
        calib_peak_xy = data.frame(x = 0, yi = max(y) , xMax = max(x), yMin = min(y))
      }else{
        calib_peak_xy = data.frame(x = 0, yi = 1 , xMax = 1, yMin = 0)
      }
      
    }else{   # If calibration samples are in different files, find first peak in each file
      calib_peak_xy = data.frame()
      for(i in 1:length(data_list)){ 
        peak_xy_i = data_list[[i]] 
        if( "type" %in% colnames(peak_xy_i)){
          peak_xy_i = peak_xy_i %>%
            filter(type == 1) %>% first() %>% select(x, yi)
        }else {
          peak_xy_i = data.frame(x = 0, yi = 0)
        }
        calib_peak_xy = rbind(calib_peak_xy, peak_xy_i)
      }
      calib_peak_xy = calib_peak_xy %>% 
        mutate(xMax = max(x, na.rm = TRUE), yMin = min(yi, na.rm = TRUE))
    }
    calib_peak_xy
  })  
  
  # (2.4) Mp values for calibration samples
  calib_Mp <- reactiveValues(data = {
    numeric(0)
  }) 
  
  calib_Mp_file <- reactive({
    # Molecular weight (in Da) is extracted from
    #                   1. input "Mp"
    #                   2. "SampleName"
    #                   3. File Name
    data_path <- input$calib_file$datapath
    
    if( length(data_path) == 1 ){
      calib_Mp_in <- rep(NA ,length( calib_peak_time$data ))
    }else{
      file_n <- input$calib_file$name
      calib_Mp_in <- file_info(data_path , file_n  )
      
      suppressWarnings({
        calib_Mp_in <- calib_Mp_in %>% 
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
          mutate(Mp = ifelse( Mp < 200, 0, Mp) ) %>%  pull(Mp)
      })
    }
    calib_Mp_in
  }) 
  
  ### For Plot 1.
  # (2.5) Axis limits for the Plot 1.
  xy_lims_calib_plot <- reactiveValues(data = {
    data.frame(x_min = numeric(0),
               x_max = numeric(0),
               y_min = numeric(0),
               y_max = numeric(0)      )
  })
  
  # (2.6) Plot 1. GPC chromatogram of standards.
  calibration_files_plot <- reactive({
    
    par(mar = c(5, 4, 2, 2) + 0.1, mgp = c(2, 1, 0))
    if(length(input$calib_file$name) == 0){
      
    }else if(length(input$calib_file$name) == 1){
      calib_df = calib_data()[[1]]
      
      if("ybl" %in% colnames(calib_df)){
        plot(calib_df$x, calib_df$ybl,
             type = "l", lwd = 2,  frame.plot = FALSE,
             xlab = "retention time",  ylab = "intensity",
             ylim = c(xy_lims_calib_plot$data$y_min, xy_lims_calib_plot$data$y_max)  ,
             xlim = c(xy_lims_calib_plot$data$x_min, xy_lims_calib_plot$data$x_max)   )
      }else{
        plot(1,  frame.plot = FALSE, xlab = "retention time", ylab = "intensity", col = "white"  )
      }
      
    }else{
      plot(1, frame.plot = FALSE, xlab = "retention time", ylab = "intensity", col = "white",
           ylim = c(xy_lims_calib_plot$data$y_min, xy_lims_calib_plot$data$y_max)  ,
           xlim = c(xy_lims_calib_plot$data$x_min, xy_lims_calib_plot$data$x_max)     )
      
      for (calib_file_i in calib_data()[1:length(calib_data())]) {
        
        if("ybl" %in% colnames(calib_file_i)){
          lines(calib_file_i$x, calib_file_i$ybl, lwd = 2)
        }else{
          next
        }
      }
    }
    abline(v = xy_from_calib_plot$data$x)
    abline(h = xy_from_calib_plot$data$y)
  })
  
  # (2.7) x and y values from the click on Plot 1. 
  xy_from_calib_plot <- reactiveValues(data = {
    data.frame(x = 0, y = 0)
  })
  
  ### For Plot 2.
  # (2.8) Calibration curve.
  calibration <- reactive({
    actual <- data.frame( calib_time = calib_peak_time$data,
                          calib_mp = calib_Mp$data  ) %>%
      drop_na() %>% filter(calib_time > 0) %>%  filter( calib_mp > 0 )
    
    if( nrow(actual) >= 2){
      coefs <- coeffs() 
      fitted_vals = data.frame(rt = seq(min(actual$calib_time), 
                                        max(actual$calib_time), 
                                        length.out = 20) ) %>%
        mutate(log_mp = coefs$coef3*(rt^3) + 
                 coefs$coef2*(rt^2) + 
                 coefs$coef1*rt + coefs$intercept)
      par(mar = c(5, 4, 2, 2) + 0.1, mgp = c(2, 1, 0))
      plot( actual$calib_time, log10(actual$calib_mp ),
            xlab = "retention time",  ylab = "log(Mp)"  )
      lines( fitted_vals, col = "purple" , lwd = 2)
    }
  })
  
  ### for help texts
  # (2.9) help text on the side panel.
  calib_help_text1 <- reactive({
    if(input$use_input_coefs == TRUE & input$calib_fit_order == "1" &
       any(coeffs()$intercept == 0, is.na(input$intercept), is.na(input$coef1) ) ){
      "One or more calibration coefficients are missing or set to an invalid value. 
      The Intercept must be greater than 0. 
      To fit a linear calibration curve, leave the 2nd and the 3rd coefficients equal to 0."
    }else if(input$use_input_coefs == TRUE & input$calib_fit_order == "3" &
             any(coeffs()$intercept == 0, is.na(input$intercept), is.na(input$coef1), is.na(input$coef2), is.na(input$coef3) ) ){
      "One or more calibration coefficients are missing or set to an invalid value. "
    }else{
      ""
    }
  })
  
  # (2.10) help text on the main panel.
  calib_help_text2 <- reactive( {
    if( 0 %in%  calib_DT$data$Mp  ){
      HTML(  ' Missing Mp values for calibration samples (zeros or empty fields in Table 1)
        will be excluded from the calculations. To change Mp or retention time (RT),  
      double-click the corresponding cells in the Table 1, enter the new  values
      and click  "Recalculate"'       )   
      
    }else if(!identical( calib_DT$data$Mp, calib_Mp$data) |
             !identical( calib_DT$data$RT, calib_peak_time$data)){
      ' To change Mp or retention time (RT), double-click the corresponding cells in the Table 1,
      enter the new values and click  "Recalculate"  '
    }else if( 0 %in%  calib_DT$data$RT){  # calib_Mp$data <- calib_DT$data$Mp
      ' Missing retention time values (RT) are set to zero in the Table 1 and will be excluded 
      from the  calculations. To enter this values manually double-click
       the corresponding cells in the Table 1, enter new values 
      and click  "Recalculate"'
    }else{
      ""
    }
  })
  
  ### calculations
  # (2.11) fit calibration curve and calculate coefficients
  fit <- reactive({
    
    df <- data.frame(RT = calib_peak_time$data, Mp =  calib_Mp$data) %>%
      drop_na() %>% filter(RT > 0 ) %>%  filter(Mp > 0 )
    
    if( nrow(df)<2  ){
      coefs <- data.frame(fit_order = c(3, 1), Intercept = c(0, 0),
                          coef1 = c(0, 0), coef2 = c(0, 0), coef3 = c(0, 0) )
    }else{
      fit_calib <- lm(log10(df$Mp) ~ poly(df$RT, 3, raw=TRUE))
      coefs <- data.frame(fit_order = 3,
                          Intercept = fit_calib$coefficients[1],
                          coef1 = fit_calib$coefficients[2],
                          coef2 = fit_calib$coefficients[3],
                          coef3 = fit_calib$coefficients[4])
      rownames(coefs) <- NULL
      fit_calib <- c( 1, lm(log10(df$Mp) ~ df$RT)$coefficients, 0, 0)
      names(fit_calib ) <- NULL
      coefs <- rbind( coefs, fit_calib )
      coefs <- coefs %>% mutate(across(everything() , as.numeric))
    }
    coefs
  })  
  
  # (2.12) Calibration coefficients (for Table 2)
  coeffs <- reactive({
    # setting an intercept and a slope/ coefficients 
    if(input$use_input_coefs == TRUE ){  
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
    }else{
      fit_all <- fit() %>%  filter( fit_order == as.numeric(input$calib_fit_order ) )
      intercept <- fit_all$Intercept 
      coef1 <- fit_all$coef1
      coef2 <- fit_all$coef2
      coef3 <- fit_all$coef3
    }
    data.frame( intercept = intercept,
                coef1 = coef1,
                coef2 = coef2,
                coef3  = coef3)
  })
  
  #                   3.  Events for the Calibration      #####
  # (3.1) Event1: input$calib_file; loading calibration file(s)
  observeEvent(input$calib_file, {
    # Get peak retention time values from calibration files and assign to calib_peak_time 
    calib_peak_time$data = calib_peak_time_file()$x
    
    # Get peak retention time values from calibration files and assign to calib_Mp 
    calib_Mp$data = calib_Mp_file()
    
    # Update Table 1 (calib_DT )
    calib_DT$data =
      data.frame(FileName = input$calib_file$name,
                 Mp = calib_Mp$data,
                 RT = round(calib_peak_time$data, 4) )
    
    # Set axis limits for the Plot 1.
    if(length(input$calib_file$datapath) == 1){
      xy_lims_calib_plot$data <-
        data.frame(x_min = 0,
                   x_max = max(calib_data()[[1]]$x, na.rm = TRUE),
                   y_min = min(calib_data()[[1]]$ybl, na.rm = TRUE),
                   y_max = max(calib_data()[[1]]$ybl, na.rm = TRUE)*1.1  )
    }else{
      xy_lims_calib_plot$data <-
        data.frame(x_min = 0,
                   x_max = first(calib_peak_time_file()$xMax),
                   y_min = first(calib_peak_time_file()$yMin),
                   y_max = max(calib_peak_time_file()$yi, na.rm = TRUE)*1.2  )
    }
  })
  
  # (3.2) Event2: input$calib_samples_table_cell_edit; editing the Table 1
  observeEvent(input$calib_samples_table_cell_edit, {
    info = input$calib_samples_table_cell_edit  # Update Table 1 (calib_DT )
    i = as.numeric(info$row)
    j = as.numeric(info$col)
    val = as.numeric(info$value)
    if(is.na(val)){ # set to 0 if NA
      val <- 0
    }
    if(val < 0){ # convert to positive if negative
      val <- val * -1
    }
    calib_DT$data[i, j] <- val
  }   )
  
  # (3.3) Event3: input$plot_click; click on the Plot 1 and show x and y values in the text box below
  observeEvent(input$plot_click, {
    validate( need(input$calib_file != "" ,    ' ') ) 
    
    xy_from_calib_plot$data$x <- input$plot_click$x
    xy_from_calib_plot$data$y <- input$plot_click$y
  })
  
  # (3.4) Event4: input$plot_brush; brush the plotting area, and
  # zoom in Plot 1 by setting new axis limits xy_lims_calib_plot 
  observeEvent(input$plot_brush, {
    validate( need(input$calib_file != "" ,    ' ') ) 
    
    xy_lims_calib_plot$data <-
      data.frame(x_min = input$plot_brush$xmin,
                 x_max = input$plot_brush$xmax,
                 y_min = input$plot_brush$ymin,
                 y_max = input$plot_brush$ymax   )
  })
  
  # (3.5) Event5: input$add_rt_from_click; click on 'add retention time' and
  # add new row to the Table 1 with RT = x
  observeEvent(input$add_rt_from_click, {
    calib_DT$data <- 
      data.frame(FileName = c( calib_DT$data$FileName, ""  ),
                 Mp = c( calib_DT$data$Mp, 0 ),
                 RT = c(calib_DT$data$RT, 
                        round(as.numeric(xy_from_calib_plot$data$x), 4) )  )
    
  }   )
  
  # (3.6)  Event6: Zoom out Plot 1 by setting axis limits
  observeEvent(input$zoom_out, {
    validate( need(input$calib_file != "" ,    ' ') )  
    
    xy_lims_calib_plot$data <-
      data.frame(x_min = min(calib_data()[[1]]$x, na.rm = TRUE),
                 x_max = max(calib_data()[[1]]$x, na.rm = TRUE),
                 y_min = min(calib_data()[[1]]$ybl, na.rm = TRUE),
                 y_max = max(calib_data()[[1]]$ybl, na.rm = TRUE)*4  )
  })
  
  #  (3.7)  Event7: Enter peak retention time value for calibration samples
  observeEvent(input$syst_peak, { # syst_peak_c
    req(input$syst_peak != 0 & !is.na(input$syst_peak) & input$calib_file != "" )   
    
    # Get peak retention time values from calibration files and assign to calib_peak_time 
    calib_peak_time$data = calib_peak_time_file()$x
    
    # Get peak retention time values from calibration files and assign to calib_Mp 
    calib_Mp$data = calib_Mp_file()
    
    # Update Table 1 (calib_DT )
    calib_DT$data =
      data.frame(FileName = input$calib_file$name,
                 Mp = calib_Mp_file(),
                 RT = calib_peak_time_file()$x )
    
    # Set axis limits for the Plot 1.
    if(length(input$calib_file$datapath) == 1){
      xy_lims_calib_plot$data <-
        data.frame(x_min = 0,
                   x_max = max(calib_data()[[1]]$x, na.rm = TRUE),
                   y_min = min(calib_data()[[1]]$ybl, na.rm = TRUE),
                   y_max = max(calib_data()[[1]]$ybl, na.rm = TRUE)*1.1  )
    }else{
      xy_lims_calib_plot$data <-
        data.frame(x_min = 0,
                   x_max = max(calib_peak_time_file()$xMax, na.rm = TRUE),
                   y_min = min(calib_peak_time_file()$yMin, na.rm = TRUE),
                   y_max = max(calib_peak_time_file()$yi, na.rm = TRUE)*1.2  )
    }
  })
  
  #  (3.8)  Event8: Click Recalculate button
  observeEvent(input$recalcutate_calib, {
    calib_Mp$data <- calib_DT$data$Mp
    calib_peak_time$data <- calib_DT$data$RT
  })
  
  #                   4.  REACTIVES for the GPC files      #####
  
  # (4.2) A list of the data frames for each file. 
  # It contains the original data and the y values after baseline correction and peak IDs.
  gpc_data  <- reactive( {  
    s_peak = input$syst_peak
    read_GPC_files = readGPC(file_path = input$gpc_file$datapath)
    baseline_GPC = baselineGPC(read_GPC_files)
    findPeaksGPC(baseline_GPC, syst_peak = s_peak )
  })   
  
  # (4.3) Data for all peaks for each GPC file
  gpc_filter <- reactive( { 
    if(input$fit_missing_tail == FALSE){
      filterPeaksGPC(gpc_data() )
    }else{
      filterPeaks_GPC = filterPeaksGPC(gpc_data() )
      fitPeaksGPC(filterPeaks_GPC,  fit_tails = 0.05, calc_coef = 0.45, add_ab = 0.1)
    }
  })
  
  # (4.4) Calculations for GPC files; calibration coefficients are used to find Mp, Mn, Mw and PDI	
  gpc_calc <- reactive( {   
    MWDcalculateGPC( gpc_filter(),
                  intercept = coeffs()$intercept ,
                  coef1 = coeffs()$coef1,
                  coef2 = coeffs()$coef2,
                  coef3 = coeffs()$coef3    )
  })
  
  # (4.5) A data frame with the results for each GPC file
  results <- reactive( {   
    use_t1 = input$use_t1 
    use_t2 = input$use_t2
    t1 = input$t1
    t2 = input$t2
    # setting an intercept and a slope/ coefficients 
    int <- coeffs()$intercept 
    c1 <- coeffs()$coef1
    c2 <- coeffs()$coef2
    c3 <- coeffs()$coef3
    
    n_files = length(gpc_calc())
    results <- map(1:n_files, function(j) {  
      file_i = gpc_calc()[[j]]   
      file_name = input$gpc_file$name[j]                 
      if( nrow(file_i) ==1 | ncol(file_i) == 2 ){      # If no data
        table_i <- data.frame(FileName = file_name, peak = NA, PDI = NA, Mn = NA,
                              Mw = NA, Mp = NA,  time1 = NA, time2 = NA,
                              time_max = NA, AssymFactor = NA, TailFactor = NA)
      }else{     
        file_i = file_i %>% group_by(peak_id) %>%
          mutate(max_y = ifelse(peak_id != 0, max(yi), 0 )) %>% 
          mutate(start_y_prop = ifelse(peak_id != 0, first(yi)/max_y, 0 ),
                 end_y_prop = ifelse(peak_id != 0, last(yi)/max_y, 0 ) ) %>%
          ungroup() %>%  group_by(peak) %>%
          mutate(peak_id = ifelse(lead(start_y_prop, 5) > 0.005 & yi > 0 &
                                    yi < lead(yi, 5) & peak_id == 0,
                                  lead(peak_id, 5) ,
                                  peak_id)     ) %>%
          mutate(peak_id = ifelse(lag(end_y_prop, 5) > 0.005 & yi > 0 &
                                    yi < lag(yi, 5) & peak_id == 0,
                                  lag(peak_id, 5) ,
                                  peak_id)     ) %>%  ungroup()
        
        peaks_summ = file_i %>% filter(peak_id != 0) %>%  group_by(peak_id) %>%
          summarise( peak = first(peak),
                     start_x = first(x),
                     start_y = first(yi),
                     end_x = last(x),
                     end_y = last(yi),
                     left = first(yi)/max(yi),
                     right = last(yi)/max(yi)) 
        
        max_xy = file_i %>% filter(peak_id != 0) %>%
          group_by(peak_id) %>% filter(yi == max(yi)) %>%
          summarise(max_y = first(yi), max_x = first(x))
        
        # assymetry factor = (width0.1 - f0.1)/f0.1
        assym_factor = file_i %>% filter(peak_id != 0 ) %>%
          left_join(max_xy, by = "peak_id") %>% select(x, yi, max_x, peak_id) %>%
          group_by(peak_id) %>% mutate(left_tail = ifelse(x < max_x, 1, 0)) %>%
          filter(yi >= max(yi)*0.1) %>%  group_by(peak_id, left_tail) %>%
          arrange(yi) %>%  filter(x == first(x)) %>%  group_by(peak_id) %>%
          arrange(left_tail) %>%
          summarise(af =  round((first(x)- first(max_x)) /(first(max_x) - last(x)),2 ) )
        
        # Tailing factor =  width0.05 / 2*f0.05. 
        tail_factor = file_i %>% filter(peak_id != 0 ) %>%
          left_join(max_xy, by = "peak_id") %>% select(x, yi, max_x, peak_id) %>%
          group_by(peak_id) %>% mutate(left_tail = ifelse(x < max_x, 1, 0)) %>%
          filter(yi >= max(yi)*0.05) %>%  group_by(peak_id, left_tail) %>% arrange(yi) %>%
          filter(x == first(x)) %>%  group_by(peak_id) %>% arrange(left_tail) %>%
          summarise(tf=  round((first(x)- last(x)) /( (first(max_x) - last(x))*2 ),2 ) )
        
        table_i = data.frame()
        for(i in 1:nrow(peaks_summ)){ 
          
          Time1 = peaks_summ$start_x[i]   # start of the peak 
          Time2 = peaks_summ$end_x[i]   # end of the peak 
          time_max = max_xy$max_x[i] # peak maximum
          
          # # if there is the user input to change start/end time of the peak
          if(input$use_t1 == TRUE){
            Time1 = Time1 + (Time2-Time1)*(input$t1/100) # change peak start time
            Time1 = ifelse(Time1 < 0, 0, Time1 )
          }
          if(input$use_t2 == TRUE){
            Time2 <- Time2 - (Time2-Time1)*(input$t2/100) # change peak end time
          }
          
          # If intensity goes below a baseline it will be dropped from the calculations
          peak_i = file_i %>% filter(peak == peaks_summ$peak[i]) %>%
            filter(x >= Time1 & x <= Time2) %>% filter(yi > 0 )  %>%
            mutate(shape = ifelse(first(yi)/max(yi) > 0.1 | last(yi)/max(yi) > 0.1, "shoulder", "peak"))
           
          Time1 = min(peak_i$x)   # start of the peak 
          Time2 = max(peak_i$x)  
          
          if( all(coeffs()$intercept != 0, !is.na(input$intercept), !is.na(input$coef1), !is.na(input$coef2), !is.na(input$coef3)) ){
            peak_i <- peak_i %>%
              summarise(Mn = round( sum(yi)/sum(yi/Mi ) ),    # number-average molecular weight
                        Mw = round( sum(yi*Mi)/sum(yi) ) ,    # weight-average molecular weight
                        shape = first(shape)  ) %>%
              mutate(peak = i,
                     Mp = round( 10^(  c3*(time_max^3)+ c2*(time_max^2)+ c1*time_max+ int )),
                     PDI = ifelse(Mw != 0 | Mn != 0  ,round(Mw/Mn, 2), NA),   # polydispersity index
                     time1 = round(Time1, 2),
                     time2 = round(Time2, 2),
                     time_max = round(time_max, 2) ) 
          }else{
            peak_i <- peak_i %>%
              summarise(Mn = NA,    # number-average molecular weight
                        Mw = NA,    # weight-average molecular weight
                        shape = first(shape)  ) %>%
              mutate(peak = i,
                     Mp = NA,
                     PDI = NA,   # polydispersity index
                     time1 = round(Time1, 2),
                     time2 = round(Time2, 2),
                     time_max = round(time_max, 2) ) 
          }
          table_i = rbind(table_i, peak_i)
        }
        table_i = table_i %>%
          mutate(AssymFactor = assym_factor$af, 
                 TailFactor = tail_factor$tf,
                 FileName = file_name) %>%
          select(FileName, peak, PDI, Mn,  Mw, Mp,  time1, time2, time_max, 
                 shape, AssymFactor, TailFactor)
      }
    }) %>% bind_rows()
    
    if(input$add_calib_info == TRUE){
      # If we want to add calibration information
      if(input$use_input_coefs == FALSE){
        # shows only 1st calibration file name
        results$CalibFile <- input$calib_file$name[1]
      }
      # adding calibration coefficients
      results$intercept <- round(int, 6)    
      results$coef1 <- round(c1, 6)
      
      if(input$calib_fit_order == "3"){
        results$coef2 <- round(c2, 6)
        results$coef3 <- round(c3, 6)
      }
    }
    results
  })
  
  #                   5.  Events for the help buttons     #####
  # Each tab has “help” buttons showing help_text1, help_text2, and help_text3 
  # describing what to do on each tab.
  observeEvent(input$help1, {
    showModal(  modalDialog(help_text1, size = "l" )    )
  })
  observeEvent(input$help2, {
    showModal( modalDialog(help_text2 )    )
  })
  
  #                               6. Render sidebar  #######
  # (6.1) If incorrect/missing values are entered in the fields, 
  # 'intercept', 'coef1', 'coef2' or 'coef3', the help text 1 is displayed.
  output$calib_help1 <- renderText({
    validate( need( input$use_input_coefs == TRUE,  '') )
    calib_help_text1()
  })
  #                               7. Render 1st TAB   "Calibration"  #######
  
  # (7.1) Table 1  shows the names of the calibration files, 
  # Mp values and peak retention times for each standard.
  output$calib_samples_table <- renderDT({
    validate( need(input$calib_file != "" , 'Select calibration files or enter "intercept" and "slope". ') )
    
    DT::datatable(calib_DT$data,
                  class = "row-border",   selection = "none",
                  options = list(searching = FALSE, filter = "none", lengthChange = FALSE),
                  editable = list(target = "cell",  disable = list(columns = 1)))
  } )
  
  # (7.2) Plot 1. GPC chromatogram of standards. ybl (intensity  after base line correction) vs x (time)
  output$calibration_files <- renderPlot( {       
    validate( need(input$calib_file != ""  ,  '') )
    calibration_files_plot()
  })
  
  # (7.3) Text box, shows x and ybl values obtained by clicking on the Plot 1 area.
  output$xy_calib <- renderText({
    validate( need(input$calib_file != "" & input$use_input_coefs == FALSE,    ' ') )
    paste0(" x = ", round(xy_from_calib_plot$data$x, 4),
           ",  y = ", round(xy_from_calib_plot$data$y,4)  )
  })
  
  # (7.4) Missing values in Table 1 are set to 0. If Table 1 contains zeros,
  # Help text 2 will be displayed.
  output$calib_help2 <- renderText({
    validate( need( input$calib_file != "" & input$use_input_coefs == FALSE, '') )
    calib_help_text2()
  })
  
  # (7.5) Table 2 shows the calibration plot coefficients calculated from 
  # calibration files (in fit()) or entered in the sidebar.
  output$calib_coef <- renderTable({
    validate( need(input$calib_file != "" ,    ' ') )
    data.frame( coefficient = colnames(coeffs() ),
                value = unlist(coeffs() ) )
  }, digits = 6)
  
  # (7.6) Plot 2. calibration curve. The logarithm of the molecular weight of
  # polymer standards (column Mp in the Table 1), vs 
  # the retention time (column RT in the Table 1, time at the calibration standard peak maximum).
  output$calibration_plot <- renderPlot( {       # Calibration curve
    validate( need( input$calib_file != "" & input$use_input_coefs == FALSE, '') )
    calibration()
  })   
  
  #                               8. Render 2nd TAB   "Data"   #######
  
  # (8.1)  Table 4 shows results for each GPC file.
  output$data_table <- DT::renderDataTable( {    # 
    validate( need( length(input$gpc_file) != 0, '  ') )
    
    DT::datatable(results(), filter = "none", class = "row-border",
                  options = list(searching = FALSE, filter = "none", lengthChange = FALSE))
  })
  
  # (8.2)  Handle saving results
  output$save_results <- downloadHandler(
    filename = function() {
      paste("results", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(results() , file)
    }
  )
  
  # (8.3) Plot. Chromatogram and molecular weight distribution (MWD)
  output$plot1 <- renderPlot( {       # Plot gpc files, original data
    validate( need( input$gpc_file != "", 'Select GPC files.') )
    
    file_name = unique(results()$FileName)
    
    col_peaks = c("#0364cc",   "#c90259", "#357503", "#E69F00","#9223f5",
                  "#0703f6", "#CC79A7","#11c902",  "#920103", "#6e03f6",
                  "#02a4c1", "#d903f6","#009E73", "#d00416", "#490182")
    if(length(input$calib_file$name) != 0 &
       all(coeffs()$intercept != 0, !is.na(input$intercept), !is.na(input$coef1), !is.na(input$coef2), !is.na(input$coef3)) |
       input$use_input_coefs == TRUE &
       all(coeffs()$intercept != 0, !is.na(input$intercept), !is.na(input$coef1), !is.na(input$coef2), !is.na(input$coef3)) ){ # chromatograms + MWD
      n_plot_rows = round(length(file_name) *2) 
      par(mfrow = c(n_plot_rows ,3),
          cex = 0.9, mgp = c(2, 1, 0),
          mar = c(4, 3, 2, 1) + 0.1)
      
      for(i in 1:length(input$gpc_file$datapath ) ) { # for each file
        df1 = gpc_data()[[i]] 
        plot1_title <- paste( input$gpc_file$name[i] )  
        
        if(nrow(df1) <= 1 ){    # if no proper data
          next
          
        }else if( ncol(df1) == 2 ){
          plot(df1, type = "l", col = "gray",
               xlab = "retention time",   ylab = "intensity",
               main = plot1_title )
          abline(0, 0, lty = 2, lwd =0.3 )
          plot.new()
          plot.new()
          
        }else{                 # if there are peaks
          res = results() %>% filter(FileName == file_name[i])
          df2 = gpc_filter()[[i]]
          df4 = gpc_calc()[[i]]  
          
          # 1. plot chromatogram
          plot(df1$x, df1$yi, type = "l", lwd = 3, col = "gray",
               xlab = "retention time",   ylab = "intensity",
               ylim = c(-max(df1$yi)*0.05, max(df1$yi)*1.15),
               xlim = c(min(df2$x)*0.7, max(df2$x)*1.3),
               main = plot1_title )
          abline(0, 0, lty = 2, lwd =0.3 )
          
          # 2. for each peak: a. add lines(); b. lable peak; c. add end/start 
          for(k in 1:nrow(res)){
            df3 = df2 %>% filter(peak == res$peak[k]) %>%
              filter(x >= res$time1[k] & x <= res$time2[k]) 
            
            if("fit_type" %in% colnames(df3)){
              lines(df3$x, df3$yi  ,lwd = 2, col = "red")
            }else{
              lines(df3$x, df3$yi  ,lwd = 2) #, col = "steelblue"
            }
            
            text(x = res$time_max[k], y = max(df3$yi)*1.1 , 
                 labels = res$peak[k], col = col_peaks[res$peak[k] ], cex = 1.2)
            points(x = c(res$time1[k], res$time2[k]), y =c(0, 0),
                   col = col_peaks[res$peak[k] ], pch = "|", cex = 2 )
          }
          
          # 3. Plot MWD for each peak
          if("Mi" %in%  colnames(df4)){
            for(k in 1:nrow(res)){   # loop for each peak
              df5 = df4 %>% filter(peak == res$peak[k]) %>%
                filter(x >= res$time1[k] & x <= res$time2[k]) 
              
              xlim_max = max(log10(df5$Mi)) + (max(log10(df5$Mi)) - min(log10(df5$Mi)))/2
              
              plot3_title =  paste("peak", res$peak[k] , ",   PDI =" , res$PDI[k])
              plot( log10(df5$Mi),  df5$dWi_over_dLogMi,    
                    xlab = "log M",   ylab = "dW/d(log M)",
                    type = "l", lwd = 3,
                    frame.plot = FALSE,
                    xlim = c(min(log10(df5$Mi)) , xlim_max),
                    main = plot3_title
              )
              abline(0, 0, lty = 2 )
              abline(v =  log10(res$Mp[k] ), col = "gray", lwd = 2)
              abline(v = log10(res$Mn[k] ), col = "blue", lwd = 2)
              abline(v = log10(res$Mw[k] ), col = "maroon", lwd = 2)
              text(x = max(log10(df5$Mi)),
                   y = max(df5$dWi_over_dLogMi)*0.9,
                   labels = paste("Mp =", res$Mp[k]), 
                   col = "gray35")
              text(x = max(log10(df5$Mi)),
                   y = max(df5$dWi_over_dLogMi)*0.7,
                   labels = paste("Mn =",res$Mn[k]) , 
                   col = "blue")
              text(x = max(log10(df5$Mi)),
                   y = max(df5$dWi_over_dLogMi)*0.8,
                   labels = paste("Mw =",res$Mw[k]) , 
                   col = "maroon")
            }
          }else{
            plot.new()
            plot.new()
          }
          
          # 4. Add empty plots to fill the row
          if( (nrow(res)+1)%%3 != 0  ){
            if(nrow(res)%%3 == 0  ){
              plot.new()
              plot.new()
            }else{
              plot.new()
            }
          }
        }
      }
      
    }else{   # only chromatogtams, if no calculations done
      par(mfrow = c(length(input$gpc_file$datapath)*1.5, 1 ) ,
          cex = 0.9, mgp = c(2, 1, 0),
          mar = c(4, 3, 2, 1) + 0.1)
      
      for(i in 1:length(input$gpc_file$datapath ) ) {
        df1 = gpc_data()[[i]] 
        
        if(nrow(df1) <= 1 ){
          next
        }
        
        plot1_title <- paste( input$gpc_file$name[i] )   
        
        if( ncol(df1) == 2 ){
          plot(df1, type = "l", col = "gray",
               xlab = "retention time",   ylab = "intensity",
               main = plot1_title )
          abline(0, 0, lty = 2, lwd =0.3 )
        }else{
          res = results() %>% filter(FileName == file_name[i])
          df2 = gpc_filter()[[i]]
          
          plot(df1$x, df1$yi, type = "l", lwd = 3, col = "gray",
               xlab = "retention time",   ylab = "intensity",
               ylim = c(-max(df1$yi)*0.05, max(df1$yi)*1.15),
               main = plot1_title )
          abline(0, 0, lty = 2, lwd =0.3 )
          
          for(k in 1:nrow(res)){
            df3 = df2 %>% filter(peak == res$peak[k]) %>%
              filter(x >= res$time1[k] & x <= res$time2[k]) 
            
            if("fit_type" %in% colnames(df3)){
              lines(df3$x, df3$yi  ,lwd = 2, col = "red")
            }else{
              lines(df3$x, df3$yi  ,lwd = 2) 
            }
            
            text(x = res$time_max[k], y = max(df3$yi)*1.1 , cex = 1.5,
                 labels = res$peak[k], col = col_peaks[res$peak[k] ])
            points(x = c(res$time1[k], res$time2[k]), y =c(0, 0),
                   col = col_peaks[res$peak[k] ], pch = "|", cex = 2 )
          }
        }
      }
    }
    
  } , 
  height = function(){50 + round(length(input$gpc_file$datapath)*1.3)*450} 
  ) 
}

# Run the application ########################
shinyServer(function(input, output, session){
  session$onSessionEnded(function() {
    stopApp()
  })
})

shinyApp(ui = ui, server = server)







