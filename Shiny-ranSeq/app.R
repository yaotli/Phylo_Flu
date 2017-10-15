require(shiny)
require(shinyjs)
require(seqinr)
require(ape)
require(stringr)

# function needed 

selectseq <- function( 
                       in.fas,
                       in.tre,
                       s.n, 
                       min.lth, 
                       max.amb, 
                       seed,
                       type, 
                       key )
{
  fasfile <- read.fasta( in.fas )
  rawid   <- attributes( fasfile )$names
  rawseq  <- getSequence( fasfile )
  
  
  if(type == "s")
  {
    
    test_lthamb <- sapply( rawseq, function( x )
    {
      y   <- gsub("-|~", "", c2s(x) ) 
      x.l <- length( s2c(y) )
      x.a <- length( grep( pattern = "a|t|c|g", x = s2c(y), invert = TRUE) )
      
      return( x.l > min.lth & x.a < max.amb )  
    })
    test_key   <- grep( pattern = key, rawid, ignore.case = TRUE )
    test_out   <- intersect( which(test_lthamb), test_key)
    
    if( grepl("\\%", s.n ) ){ s.n = floor(length(test_out)*as.numeric( str_match(s.n, "[0-9.]+")[1] )/100)  }else
    { s.n = as.numeric(s.n) }
    
    set.seed( seed )
    s = test_out[ sample( length(test_out), s.n) ] 
    
  }else{
    
    trefile <- read.nexus( in.tre )
    treid   <- gsub("'", "", trefile$tip.label)
    
    test_lthamb <- sapply( rawseq, function( x )
    {
      y   <- gsub("-|~", "", c2s(x) ) 
      x.l <- length( s2c(y) )
      x.a <- length( grep( pattern = "a|t|c|g", x = s2c(y), invert = TRUE) )
      
      return( x.l > min.lth & x.a < max.amb )  
    })
    
    test_key   <- match( grep( pattern = key, treid, ignore.case = TRUE, value = TRUE), rawid )
    test_out   <- intersect( which(test_lthamb), test_key)
    
    if( grepl("\\%", s.n ) ){ s.n = floor(length(test_out)*as.numeric( str_match(s.n, "[0-9.]+")[1] )/100)  }else
    { s.n = as.numeric(s.n) }
    
    set.seed( seed )
    s = test_out[ sample( length(test_out), s.n) ]
    
  }
  
  rawout <- list( rawid[s], rawseq[s] )
  return(rawout)
  
}

### ui --------------------------------

ui <- 
  navbarPage( "ranSeq", 
              
              # Main panel
              tabPanel( "Main", 
                        fluidRow( 
                          column(3, 
                                 wellPanel(
                                   
                                   #  input_fas
                                   fileInput( "input_fas", "Your FAS file", multiple = FALSE),
                                   # seed
                                   sliderInput( "seed", "Set seed", min = 1, max = 1000, value = 666),
                                   # s.n
                                   textInput( "s.n", label = "Sample number or %", value = "10%"),
                                   # min.lth
                                   textInput( "min.lth", label = "Min length", value = "10"),
                                   # max.amb
                                   textInput( "max.amb", label = "Max ambiguous nt.", value = "1000"),
                                   # download botton
                                   actionButton( "run", label = "Run")
                                   
                                 ),
                                 
                                 wellPanel(
                                   h4("More options"),
                                   # key
                                   textInput( "key", label = "Keyword", value = ""),
                                   # type
                                   radioButtons( "type", "Type", c("sequence only" = "s", "with tree" = "t") ),
                                   # input_tre
                                   fileInput( "input_tre", "A tree may be helpful (NEX)", multiple = FALSE)
                                   
                                  )
                          ),
                         column(4, 
                                wellPanel(
                                h4("Download when you see DONE"),  
                                # state
                                verbatimTextOutput("state", TRUE),
                                # download botton
                                downloadButton( "download", "Result")
                                ), tableOutput("table"))
                                
                        )
                        
                        ),
              
              tabPanel( "Readme") )

### server --------------------------------

server  <- function( input, output )
{
  options(shiny.maxRequestSize=50*1024^2)
  # table
  inputvalues <- reactive(
    {
      data.frame(
      Criteria = c( "Sampling", "Min.lth", "Max.amb", "useTree",  "useKeyword", "Seed"),
      Values   = c( input$s.n, 
                    input$min.lth,
                    input$max.amb,
                    ifelse( input$type == "s", "No", "Yes" ),
                    ifelse( input$key == "", "No", "Yes" ),
                    input$seed ), stringsAsFactors = FALSE )  
    })
  
  output$table <- renderTable( { inputvalues() } )
  
  # uploaded file
  
  
  # run the file

  show <- eventReactive(input$run, 
                        {
                        
                          fasfile <- read.fasta( input$input_fas$datapath )
                          
                          temp0 = 
                          selectseq( 
                                     in.fas  = input$input_fas$datapath,
                                     in.tre  = input$input_tre$datapath,
                                     s.n     = input$s.n,
                                     min.lth = as.numeric(input$min.lth),
                                     max.amb = as.numeric(input$max.amb),
                                     seed    = as.numeric(input$seed),
                                     type    = input$type,
                                     key     = as.character(input$key) )
                          
                          return( paste0("Done with n = ", length(temp0[[1]]) ) )
                        })
  
  
  out <- eventReactive(input$run, 
                        {
                          
                          fasfile <- read.fasta( input$input_fas$datapath )
                          
                          temp = 
                            selectseq( 
                              in.fas  = input$input_fas$datapath,
                              in.tre  = input$input_tre$datapath,
                              s.n     = input$s.n,
                              min.lth = as.numeric(input$min.lth),
                              max.amb = as.numeric(input$max.amb),
                              seed    = as.numeric(input$seed),
                              type    = "s",
                              key     = as.character(input$key) )
                          
                          return( temp )
                          
                        })
  
  output$state <- renderText(  show() )
  
  output$download <- downloadHandler( filename = "subSeq", 
                                      content = function( file )
                                        {
                                        write.fasta( sequences = out()[[2]], 
                                                     names     = out()[[1]], 
                                                     file.out  = file )
                                        
                                      } )
  
}


shinyApp(ui = ui, server = server)