# Run the application

shinyServer(function(input, output, session){

  session$onSessionEnded(function() {
    stopApp()
  })

})

shinyAppDir("inst/app")  # , options = list()

# runApp("inst/app")
# shinyApp(ui = ui, server = server)



