# GPC reader,  UI

ui <- fluidPage(
  h4(" GPC reader "),    # Header

  sidebarLayout(         # Sidebar layout

    # Sidebar Panel -----------------------------------------------------------

    sidebarPanel(
      width = 2,
      tags$style(".well {background-color: #FAFAFA;}"),

      #                            Calibration information

      # Select the equation for the calibration curve. For a linear calibration
      # equation, the input 'calib_fit_order' is 1, for a polynomial, it is 3.
      radioButtons("calib_fit_order",
                   "Calibration:",
                   selected = "1",
                   choices = c("linear" = "1", "polynomial 3" = "3")),
      div(HTML('<hr style="border-color: #5d6d7e;">')),

      # Calibration coefficients (numeric input):
      numericInput("intercept",
                   "Enter the Intercept",
                   value =  10.99),

      numericInput("coef1",
                   "and the Slope",
                   value = -0.65),

      numericInput("coef2",
                   "for a polynomial calibration curve, enter the 2nd ",
                    value = 0),

      numericInput("coef3",
                   "and 3rd coefficients ",
                   value = 0),

      # Checkbox input. If 'use_input_coefs' is checked, then manually entered
      # calibration coefficients will be used.
      checkboxInput("use_input_coefs",
                    "use these coefficients",
                    value = FALSE),

      # Help text is displayed if the entered coefficients are missing or invalid
      span(textOutput("calib_help1"),
           style="color:firebrick"),
      div(HTML('<hr style="border-color: #5d6d7e;">')),

      # System peak for calibration standards (numeric input).
      # The calibration chromatogram can be trimmed by 'syst_peak_calib' value.
      # (Plot1 in the Calibration tab shows these chromatograms)
      numericInput("syst_peak_calib",
                   "System peak, calibration",
                   value = 0,
                   min = 0,
                   step = 0.01),

      # System peak for samples (numeric input). The chromatogram(s) can be
      # trimmed by 'syst_peak' value (graphs in the Data tab)
      numericInput("syst_peak",
                   "System peak, samples",
                   value = 0,
                   min = 0,
                   step = 0.01),
      div(HTML('<hr style="border-color: #5d6d7e;">')),

      # Checkbox input. If 'fit_missing_tail' is checked, missing tales of the
      # peaks will be added.
      checkboxInput("fit_missing_tail",
                    strong("Add tails"),
                    value = FALSE),
      div(HTML('<hr style="border-color: #5d6d7e;">')),

      # Peak boundaries may be changed.
      strong("Change the peak boundaries:"),

      # the 'use_t1' checkbox and the numeric input 't1' are intended to shift
      # the start of the peak by the value t1
      numericInput(inputId = "t1",
                   "",
                   step = 1,
                   value = 0,
                   min = -20,
                   max = 35),

      checkboxInput("use_t1",
                    "shift start time",
                    value = FALSE),

      # the 'use_t2' checkbox and the numeric input 't2' are intended to shift
      # the end of the peak by the value t2
      numericInput(inputId = "t2",
                   "",
                   step = 1,
                   value = 0,
                   min = -20,
                   max = 35),

      checkboxInput("use_t2",
                    "shift end time",
                    value = FALSE)
    ),

    # Main Panel -------------------------------------------------------------

    mainPanel(
      width = 10,
      tags$style(".shiny-file-input-progress {display: none}"),
      tabsetPanel(

        ##  Calibration TAB ------------------------------
        tabPanel(
          'Calibration',
          div(HTML("<br>")),

          fluidRow(

            # File input 'calib_file' is for loading files with calibration
            # chromatograms.
            column(
              width = 6,
              fileInput("calib_file", "Select calibration file(s)",
              multiple = TRUE,
              accept = c(".txt", ".arw", ".xls", ".xlsx", ".xlsm", ".xltx",
                         ".xltm", ".csv", ".tsv"))
              ),

            # Zoom out button for Plot1
            column(
              width = 2,
              actionButton("zoom_out", "Zoom Out",
                           style="color: black; background-color: #F0F8FF; border-color: #808b96",
                           width = '100%')
              ),
            column(width = 2),

            # Help button for the Calibration tab
            column(
              width = 2,
              actionButton("help1", "Help",
              style="color: black; background-color: #EAEBEC; border-color:  #808b96",
              width = '100%') )
            ),

          fluidRow(

            # Table 1 output
            column(
              width = 6,
              strong(em("Table 1. Calibration standards") ),
              DTOutput('calib_samples_table'),

             # The input 'first_peak' checkbox is displayed as 'one sample per
             # file'. If checked, then for files with calibration standards,
             # only the first peak from each chromatogram is detected.
              checkboxInput("first_peak",
                            "one sample per file",
                            value = TRUE)
              ),

            # Plot1 output
            column(
              width = 6,
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

            # Help text is displayed if values in Table 1 are missing or equal
            # to 0, or if there are less than 6 standard points left in Table 1
            # after all corrections
            column(
              width = 6,
              span(textOutput("calib_help2"), style="color: #565374"),
              div(HTML("<br>"))
              ),

            # A text output box that shows the x and y axis values after
            # clicking on Plot1
            column(
              width = 6,

              column(
                width = 8,
                verbatimTextOutput("xy_calib", placeholder = TRUE)
                ),

              # The 'add_rt_from_click' button is displayed as 'Add RT'. It adds
              # the x value from the 'xy_calib' text field to the RT column of
              # Table1.
              column(
                width = 4,
                actionButton("add_rt_from_click",
                             "Add RT",
                             style="color: black; background-color: #F0F8FF; border-color: #808b96",
                             width = '100%')
                )
              ),
            ),

          fluidRow(

            # Table 2 output
            column(
              width = 6,
              div(HTML("<br>")),
              strong(em("Table 2. Calibration coefficients")),
              tableOutput("calib_coef")
              ),

            # Plot 2 output
            column(
              width = 6,
              div(HTML("<br>")),
              strong(em("Plot 2. GPC calibration curve")),
              plotOutput('calibration_plot', height = "380px")
              )
            )
          ),

        ##  Data TAB -----------------------------------

        tabPanel(
          'Data',
          div(HTML("<br>")),

          fluidRow(

            # File input 'gpc_file' is for loading files with sample chromatograms.
            column(
              width = 10,
              fileInput("gpc_file",
                        "Select GPC file(s)",
                        accept = c(".txt", ".arw", ".xls", ".xlsx", ".xlsm",
                                   ".xltx",  ".xltm", ".csv", ".tsv"),
                        multiple = TRUE)),

            # Help button for the Data tab
            column(
              width = 2,
              actionButton("help2",
                           "Help",
                           style="color: black; background-color: #EAEBEC; border-color: #808b96",
                           width = '100%'))
            ),
          div(HTML("<br>")),

          # Table 3 output
          strong(em("Table 3. Results") ),
          DT::DTOutput('data_table'),

          fluidRow(

            # Checkbox input. If the input "add_calib_info" is checked, then 4 more
            # columns with calibration coefficients are added to the results table.
            column(
              width = 2,
              checkboxInput("add_calib_info",
                            "Add calibration coefficients",
                            value = FALSE)
              ),

            # Save results button
            column(
              width = 2,
              downloadButton("save_results",
                             "Save",
                             style="color: black; background-color: #F0F8FF; border-color: #808b96",
                             width = '100%')
              ),

            # Help text is displayed if the PDI values in Table 3 for peak 1 are
            # greater than 10 or all are missing, or if calibration coefficients
            # are not provided (these coefficients should be shown in Table 2).
            column(
              width = 8,
              span(textOutput("help3"), style="color:firebrick"))
            ),

          div(HTML('<hr style="border-color: #CBD1E6;">')),

          # Chromatograms and MWD Plots
          plotOutput('plot1')
        )
      )
    )
  )
)











