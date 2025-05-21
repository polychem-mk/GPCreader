calibration_coeffs <-
function(rt.mp = NULL, rt = NULL, mp = NULL){

  #  Mp (peak molecular weight) and RT (retention time)

  # To find the intercept and other coefficients, calibration_coeffs() requires
  # either a data frame 'rt.mp' with two columns (RT and Mp), or two vectors mp
  # and rt. If Mps were not in the SampleName, they can be listed in mp.

  # Thus calibration_coeffs() looks for retention time in rt.mp, then in rt
  # and it looks for molecular weight in mp if it is not NULL, even if
  # there is a Mp column in rt.mp

  if(any(class(rt.mp) %in% c("tbl_df","tbl","data.frame"))){

    if("RT" %in% colnames(rt.mp)){
      retention_time = rt.mp$RT
    }else{
      retention_time = rt
    }

    if(all(!is.null(mp), length(mp) == length(retention_time))){
      molecular_weight = mp

    }else if("Mp" %in% colnames(rt.mp)){
      molecular_weight = rt.mp$Mp

      if(all(molecular_weight == 0)){
        molecular_weight = NULL
      }

    }else{
      molecular_weight = NULL
    }

    if(any(!c("RT", "Mp") %in% colnames(rt.mp))){
      warning("'rt.mp' must be a data frame with columns RT and Mp")
    }

  }else{
    retention_time = rt
    molecular_weight = mp

    if(any(!is.null(rt.mp))){
      warning("'rt.mp' must be a data frame with columns RT and Mp")
    }
  }

  if(any(is.null(retention_time), is.null(molecular_weight),
         length(molecular_weight) != length(retention_time))){
    warning("'mp' or 'rt' are missing or has different number of values")
    return(NULL)
  }

  # Filter Mp and RT
  RT_Mp = data.frame(RT = retention_time, Mp = molecular_weight) %>%
    filter(Mp != 0) %>%
    filter(RT != 0) %>%
    filter(!is.na(RT))

  # If retention time or molecular weight values are missing or have different
  # number of values, calibration_coeffs() returns NULL with a warning.
  if(nrow(RT_Mp) < 2){
    # If there are less than 2 data points, then further calculations are
    # not possible.
    warning("not enough data points to find calibration coefficients")
    return(NULL)
  }

  # To construct a calibration curve, at least 6 standard points are required.
  # It is recommended to use 10 to 12 calibration standards. However, the
  # calibration_coeffs() function attempts to find calibration coefficients if
  # there are at least two data points.

  # Find the coefficients for the calibration curve as:

  # linear:         log10(Mp) = coef1*RT + intercept
  fit_calib_linear = lm(log10(RT_Mp$Mp) ~ RT_Mp$RT)

  coefs = data.frame(fit_order = 1,
                     intercept = fit_calib_linear$coefficients[1],
                     coef1 = fit_calib_linear$coefficients[2],
                     coef2 = 0,
                     coef3 = 0)
  rownames(coefs) = NULL

  # polynomial:     log10(Mp) = coef3*RT^3 + coef2*RT^2 + coef1*RT + intercept

  fit_calib = lm(log10(RT_Mp$Mp) ~ poly(RT_Mp$RT, 3, raw=TRUE))

  fit_calib <- c( 3, fit_calib$coefficients)
  names(fit_calib ) = NULL

  # Output is a data frame with columns fit_order, intercept, coef1, coef2 and
  # coef3. Where fit_order 1 is for a linear calibration curve and fit_order 3
  # is a polynomial
  rbind( coefs, fit_calib ) %>%
    mutate(across(everything() , as.numeric))
}
