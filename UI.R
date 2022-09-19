if (interactive()) {
  options(device.ask.default = FALSE)
shinyApp(
  ui = fluidPage(
    sidebarPanel(
            fileInput(
            "GFA_input","GFA data input",multiple=FALSE, buttonLabel="Browse",placeholder="No file selected"
            ),
            fileInput(
              "GAF_input2","GAF data input",multiple=FALSE, buttonLabel="Browse",placeholder="No file selected"
            ),
            fileInput(
              "Annotation_Input","Add Annotation",multiple=FALSE, buttonLabel="Browse",placeholder="No file selected"
            ),
#Not so sure about this part yet
            selectInput(
              "select_graph","Graph",selected = NULL, multiple = FALSE,list("outputfile1","outputfile2")),
            
            downloadButton("Fasta_download", "FASTA file Download",icon = shiny::icon("download")),
            downloadButton("BED_download", "BED file Download",icon = shiny::icon("download"))
      ),
    mainPanel(
      "Graph Visualization Window",
      plotOutput(
      "visualization","Graph Visualization Window",width="100%",height="400px",brush=brushOpts(id="plot_brush")
      ),
      "Linear Visualization Window",
      plotOutput(
        "visualization",width="100%",height="400px",brush=brushOpts(id="plot_brush")
      )
    )),
    
  server = function(input, output,session) {

  }

  )}