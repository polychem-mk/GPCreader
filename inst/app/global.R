# Libraries --------------------------------------

#  R version 4.4.3 (2025-02-28)

library(tidyr)    # tidyr_1.3.1
library(stringr)  # stringr_1.5.1
library(dplyr)    # dplyr_1.1.4.9000
library(purrr)    # purrr_1.0.2

#library(here)    # here_1.0.1
#library(tools)   # base package

library(shiny)    # shiny_1.9.1
library(DT)       # DT_0.33


# Functions  ------------------------------------------

# (1) file_info() reads the first n lines of a chromatogram file, extracts the
# sample name, finds the line where the numeric data starts, and the number of
# numeric columns. This is summarized in a data frame that also contains the
# file name and file path

# (2) readGPC() reads a chromatogram file, extracts numerical data as x
# (retention time or retention volume) and y (detector signal intensity).

# (3) find_peaks() performs intensity correction and finds peaks.

# (4) filter_peaks() filters peaks and adds peak tails if they are missing.

# (5) calibration_info() function finds RT (retention time) from chromatograms of
# calibration standards and Mp (peak molecular weight) values from sample names.

# (6) calibration_coeffs() function calculates coefficients for linear and
# polynomial calibration curve from RT and Mp values.

# (7) calculate_MWD() function calculates the molecular weight distribution.

# (8) calculatePDI() function calculates Mp (peak molecular weight), Mn, Mw
# (average molecular weight) and PDI (polydispersity index). It also estimates
# the asymmetry factor and the tailing factor.

# (9) plotMWD() function plots two graphs: a chromatogram and a molecular weight
# distribution

path_functions <- list.files(here::here("R"), full.names = TRUE)

for(i in path_functions){
  source(i)
}

# Help texts  -----------------------------------------------

help_text4 <- HTML("<h4> Calibration</h4>
There are two ways to set up calibration. <br>
1. If you know the Intercept and Slope coefficients, you can enter these values
in the corresponding fields in the sidebar. <br>
Then check the <strong>use these coefficients</strong> checkbox. <br>
For a linear calibration curve: <br>
&nbsp;&nbsp;&nbsp;&nbsp;∘ set the  Intercept (it is greater than 0);<br>
&nbsp;&nbsp;&nbsp;&nbsp;∘  the Slope (it is negative number);<br>
&nbsp;&nbsp;&nbsp;&nbsp;∘ leave the second and the third coefficients equal to 0. <br>

For a polynomial calibration curve:<br>
&nbsp;&nbsp;&nbsp;&nbsp;∘ switch the <strong>Calibration</strong> to <em>polynomial 3</em>;<br>
&nbsp;&nbsp;&nbsp;&nbsp;∘ enter Intercept  and other coefficients;<br>
&nbsp;&nbsp;&nbsp;&nbsp;∘ the coefficients must be very precise;
do not round them for polynomial  fit.<br><br>

2.  Another way to set up calibration is to download GPS data file(s) with calibration standards. <br>
2.1. If   each calibration sample is in separate file, upload all calibration files at once
(<span style='color:navy;'>each calibration standard has been injected into the GPC
separately</span>). The <strong>one sample per file</strong> checkbox is enabled by default. <br>

&nbsp;&nbsp;&nbsp;&nbsp;This application extracts <em>Mp</em> values for each calibration sample
from the SampleName or FileName (<strong>Table1</strong>, the second column).
If these values are incorrect or missing, double-click the appropriate cells
in the <strong>Table1</strong>  and enter the new <em>Mp</em>(s) values.
If necessary, do the same for the peak <em>retention time</em> values
(<strong>Table1</strong>, the third column).<br><br>

2.2.Upload one or more files, if <span style='color:navy;'> a mixture of calibration
standards has been injected into the GPC</span>. Uncheck <strong>one sample per file</strong> checkbox.
Add the <em>Mp</em> values to the 2nd column of <strong>Table1</strong>
(double-click the corresponding cells). <br><br>

If only one calibration file is provided, the application will automatically assume
that all calibration samples are in this file and will try to find all peaks and shoulder peaks.

<br><br>

• The <em>retention time</em> values can also be added from <strong>Plot 1</strong>.
If you click on this plot, the x and y values will appear in the text box below.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;To add this RT value to the <strong>Table1</strong>,
click the <strong>Add RT</strong> button.
<br>
• Missing values or values equal to zero
in the <strong>Table1</strong>  will be excluded from  calculations.
<br>
• You can also trim the calibration chromatograms by entering a time in
                   the sidebar in the <strong>System Peak, Calibration</strong> field."  )

help_text5 <- HTML("<h4>Data tab</h4>
1. Upload GPC file(s)  <br><br>
2. The calculations are summarized in the  <strong>Table 3. Results</strong>.
Check the <strong>Add calibration coefficients</strong> box if you want this
information to be added to the table.
To save the results click on <strong>Save</strong> button. <br><br>
3. Plots. This tab displays the chromatogram and MWD plots.
<br><br>
If the intercept and slope coefficients are missing or set to incorrect values,
only the chromatogram is displayed and the PDIs cannot be calculated and are
not displayed in the results table.
<br><br>
<h4>Sidebar</h4>
• You can switch between <strong>linear</strong> and <strong>polynomial 3</strong>
 calibration; the <strong>Results</strong> table will be automatically updated. <br><br>
• To cut off the right part of the chromatogram enter the value for the <strong>System peak, samples</strong> time.<br><br>
• If the peak is incomplete, to add the  missing tails check
<strong>Add tails</strong> box.<br><br>
•  To shift the <em>start</em> and <em>end</em> of the peak, enter values in the sidebar above the
<strong>shift start time</strong> / <strong>shift end time</strong> checkboxes.
These values are percentages of the peak width. Check the boxes below to apply the changes.
                   ")



