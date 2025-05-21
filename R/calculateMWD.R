calculateMWD <-
function(dataGPC,
                           intercept = 0 ,
                           coef1 = 0,
                           coef2 = 0,
                           coef3 = 0
){

  # (1) Check input   -------------------------

  # (1.1) intercept, coef1, coef2, coef3
  # intercept must be a numeric value; otherwise calculateMWD() returns dataGPC
  # with warning. If coef2 and coef3 are NA, they are set to 0. But intercept
  # and coef1 are required

  int_coef1 = TRUE

  if(any(!is.vector(intercept, mode = "numeric"),
         !length(intercept) == 1,
         !is.vector(coef1, mode = "numeric"),
         !length(coef1) == 1)){

    warning(" coefficients must be numeric... calculations were not performed... ")

    int_coef1 = FALSE
  }

  if(any(is.na(coef2),
         !is.vector(coef2, mode = "numeric"),
         !length(coef2) == 1)){
    coef2 = 0
  }

  if(any(is.na(coef3),
         !is.vector(coef3, mode = "numeric"),
         !length(coef3) == 1)){
    coef3 = 0
  }

  # Intercept must be greater than 0, for linear approximation coef1 must be negative.
   if(coef2 == 0 & coef3 == 0 & int_coef1){
     if(any(intercept < 0,coef1 > 0 )){

       int_coef1 = FALSE
     }
   }

  # (1.2) dataGPC
  # dataGPC must be a list of data frames or a data frame to be converted to a list.
  # If it is not, filter_peaks() returns NULL.
  if(any(class(dataGPC) %in% c("tbl_df","tbl","data.frame"))){
    dataGPC = list(as.data.frame(dataGPC))

  }else if(any(!"list"  %in% class(dataGPC) )){
    warning(" 'dataGPC' must be a list of data frames ... NULL is returned ...")
    return(NULL)
  }

  # (2) Map by file -----------------------

  # (2.1)
  data_MWD =
  map(dataGPC, function(dat) {

    if(int_coef1){
      # If there is intercept and coeff1 (int_coef1 is TRUE),
      # then calculations are possible.
      nr = nrow(dat)

      if(any(is.null(dat) ,
             nr <= 1,
             !c("x", "yi", "peak_id", "peak") %in%  colnames(dat) ) ){
        # If there is no data, the MWD graph will be empty and no calculations
        # will be performed.
        return(NULL)
      }

      dif_x  =  mean(diff(dat$x))  # time difference between two measurements

      #  Calculations.
      dat = dat %>%
        # set 'yi_above0' to remove all negative y values
        mutate(yi_above0 = ifelse(yi < 0 | peak_id == 0, 0, yi ) ) %>%
        group_by(peak) %>%
        mutate(p = coef3*(x^3)+ coef2*(x^2)+ coef1*x+ intercept) %>%
        ungroup()

      p_max = dat  %>% filter(peak_id != 0 ) %>%
        summarise( max(p, na.rm = TRUE)) %>%
        pull()

      if(p_max > 10){
        # The calculated molecular weight is 10^p. Polymers can have millions of
        # carbon atoms, so setting p_max to 10 is more than enough to cover the
        # maximum possible molecular weight and discard results if the calibration
        # coefficients are incorrect and result in p values that are too large.
        dat = dat %>%
          group_by(peak) %>%
          mutate(
            Mi = NA,       # Mi is molecular weight at the i-th point;
            dWi_over_dVi = NA,  # Values for the molecular weight distribution:
            dVi_over_dLogMi = NA,
            dWi_over_dLogMi = NA
          ) %>% ungroup()

      }else{

        dat = dat %>%
          group_by(peak) %>%
          mutate(
            Mi = 10^p,   # Mi is molecular weight at the i-th point;
            # Values for the molecular weight distribution:
            dWi_over_dVi = yi/sum(yi_above0)/dif_x,
            dVi_over_dLogMi = abs(1/ (coef1 + coef2*2*x + coef3*3*x^2) ),
            dWi_over_dLogMi = dWi_over_dVi*dVi_over_dLogMi
          ) %>% ungroup()

        if("y_fit" %in% colnames(dat)){

          # If there is a y_fit column, two more columns will be added:
          # dWi_over_dVi_fit and dWi_over_dLogMi_fit.

          # Mi and dVi_over_dLogMi are functions of x, so there is no need
          # to recalculate them.

          dat = dat %>%
            mutate(yi_above0 = ifelse(y_fit < 0 | fitted != 1, 0, y_fit ) ) %>%
            group_by(peak) %>%
            mutate(
              # Values for the molecular weight distribution for the fitted values:
              dWi_over_dVi_fit = y_fit/sum(yi_above0)/dif_x,
              dWi_over_dLogMi_fit = dWi_over_dVi_fit*dVi_over_dLogMi
            ) %>% ungroup()
        }
      }

      dat %>%
        select(-yi_above0)
    }else{
      dat
    }
  } )

  # (2.2) calibration coefficients (input):
  coeffs = c(intercept = intercept,
             coef1 = coef1,
             coef2 = coef2,
             coef3 = coef3)

  # (3) The output is a list of length 2.

  # (3.1) The first is the list of the data frames:
  # The columns x, yi, peak_id, peak, y_fit, fit, type are kept from the
  # previous functions since the output of the calculateMWD() can be used
  # to plot more than just MWD

  # The new columns are
  # Mi  molecular weight at the i-th point;
  # dWi_over_dVi (or dWi_over_dVi_fit) is  molecular weight fraction over
  #              a constant volume interval dV (the intensity (y) divided by
  #              the total intensity above the baseline over a dif_x);
  # dVi_over_dLogMi is the weight fraction per unit log Mi;
  # dWi_over_dLogMi and dWi_over_dLogMi_fit" is differential molecular weight.
  # [S. Mori and H. G. Barth, Size Exclusion Chromatography, 1st ed.
  #                                        Springer Berlin, Heidelberg, 1999.]

  # (3.2) The second is the vector 'coeffs' with calibration coefficients.

  return(list(data_MWD, coeffs))

}
