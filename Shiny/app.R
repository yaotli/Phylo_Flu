library(shiny)
library(shinyjs)

library(ape)
library(ggtree)
library(stringr)

library(XML)
library(RCurl)

source("helpers.R") 

# function needed
# BLAST related functions were derived from annotate::blastSequences

 parsexmltoDF <- function(xml){
  
  iter <- xml["//Iteration"]
  iterlen <- sapply(iter, xpathSApply, "count(.//Hsp)")
  iterdf <- xmlToDataFrame(iter, stringsAsFactors=FALSE)
  
  hit <- xml["//Hit"]
  hitlen <- sapply(hit, xpathSApply, "count(.//Hsp)")
  hitdf <- xmlToDataFrame(hit, stringsAsFactors=FALSE)
  hitdf <- hitdf[, names(hitdf) != "Hit_hsps", drop=FALSE]
  
  hsp <- xmlToDataFrame(xml["//Hsp"] , stringsAsFactors=FALSE)
  
  df <- cbind(
    iterdf[rep(seq_len(nrow(iterdf)), iterlen),, drop=FALSE],
    hitdf[rep(seq_len(nrow(hitdf)), hitlen),, drop=FALSE],
    hsp)
  rownames(df) <- NULL
  df
}
  parseResult <- function(baseUrl, rid, rtoe){
  
  timeout = 9999
  
  start <- Sys.time()
  end <- Sys.time() + timeout
  url <- sprintf("%s?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=%s",
                 baseUrl, rid)
  
  Sys.sleep(min(rtoe, timeout))
  
  repeat {
    elapsed <- as.double(Sys.time() - start, units="secs")
    
    result <- as(htmlParse(getURL(url, followlocation=TRUE),
                           error = xmlErrorCumulator(immediate=FALSE)),
                 "character")
    
    if (grepl("Status=FAILED", result))
      
      stop("BLAST search failed")
    
    else if  (grepl("Status=UNKNOWN", result))
      
      stop("BLAST search expired")
    
    else if (grepl("Status=READY", result)) {
      url <- sprintf("%s?RID=%s&FORMAT_TYPE=XML&CMD=Get", baseUrl, rid)
      
      result <- xmlParse(getURL(url, followlocation=TRUE),
                         error = xmlErrorCumulator(immediate=FALSE))
      return(result)
      
    } else if (grepl("Status=WAITING", result)) {
      
      if (Sys.time() > end && interactive()) {
        msg <- sprintf("wait another %d seconds? [y/n] ", timeout)
        
        repeat {
          ans <- substr(trimws(tolower(readline(msg))), 1, 1)
          
          if (ans %in% c("y", "n"))
            break
        }
        if (ans == "n")
          break
        
        end <- Sys.time() + timeout
      }
      Sys.sleep(10)
      
    } else
      stop("BLAST search unknown response") 
  }
  msg <- sprintf("'blastSequences' timeout after %.0f seconds",
                 elapsed)
  stop(msg, call.=FALSE)
}

     blastseq <- function(x){
  
  t1 = Sys.time()
  library(XML)
  library(RCurl)
  
  database = "nr"
  hitListSize = "10"
  filter = "L"
  expect = "2"
  program = "blastn"
  
  baseUrl <- "https://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
  
  query <- paste("QUERY=", URLencode(as.character(x)), "&DATABASE=", database,
                 "&HITLIST_SIZE=",hitListSize,"&FILTER=",filter,
                 "&EXPECT=",expect,"&PROGRAM=",program, sep="")
  
  url0 <- sprintf("%s?%s&CMD=Put", baseUrl, query)
  post <- htmlParse(getURL(url0, followlocation=TRUE))  
  
  x <- post[['string(//comment()[contains(., "QBlastInfoBegin")])']]
  rid <- sub(".*RID = ([[:alnum:]]+).*", "\\1", x)
  rtoe <- as.integer(sub(".*RTOE = ([[:digit:]]+).*", "\\1", x))
  
  result <- parseResult(baseUrl, rid, rtoe)
  y = parsexmltoDF(result)
  z = y[[9]]
  td = Sys.time() - t1
  
  print(td)
  
  return(z)
}

    blastinfo <- function(x){
      
          x <- gsub(" ", "_", x)      
          nn <- "A/([A-Za-z0-9-_]+)/([A-Za-z0-9-_]+)/([A-Za-z0-9-_]+)"
      hit_id <- str_match(x, nn)[,1]
      
      hitlabel = c()
      for(i in 1: 10){
        
        add = grep(hit_id[i], treedata_subTree[,10])
        hitlabel = c(hitlabel, add)
        
      }
      
      hitlabel = unique(hitlabel)
      
      return(hitlabel)
    }

# read tree file

            subTreefile <- read.tree("www/sub_7225")
       treedata_subTree <- fortify(subTreefile)
  treedata_subTree[,10] <- gsub("'", "", treedata_subTree$label)
  treedata_subTree[,11] <- "grey"
  
  colnames(treedata_subTree)[11] = "hit_color"
  
  H5_tip <- data.frame(node= c(1:length(subTreefile$tip.label)), shape = NA)  
  
# ui ####
  
ui <- fluidPage(
  
  useShinyjs(),
  tags$style(appCSS),

  # textinput
  textInput(inputId = "seq", 
            label = "Example: Gs/GD/1/96_H5N1", 
            value = "gatcagatttgcattggttaccatgcaaacaactcgacagagcaggttgacacaataatggaaaagaacgttactgttacacatgcccaagacatactggaaaagacacacaatgggaagctctgcgatctaaatggagtgaagcctctcattttgagagattgtagtgtagctggatggctcctcggaaaccctatgtgtgacgaattcatcaatgtgccggaatggtcttacatagtggagaaggccagtccagccaatgacctctgttacccaggggatttcaacgactatgaagaactgaaacacctattgagcagaacaaaccattttgagaaaattcagatcatccccaaaagttcttggtccaatcatgatgcctcatcaggggtgagctcagcatgtccataccatgggaggtcctcctttttcagaaatgtggtatggcttatcaaaaagaacagtgcatacccaacaataaagaggagctacaataataccaaccaagaagatcttttagtactgtgggggattcaccatcctaatgatgcggcagagcagacaaagctctatcaaaacccaaccacttacatttccgttggaacatcaacactgaaccagagattggttccagaaatagctactagacccaaagtaaacgggcaaagtggaagaatggagttcttctggacaattttaaagccgaatgatgccatcaatttcgagagtaatggaaatttcattgctccagaatatgcatacaaaattgtcaagaaaggggactcagcaattatgaaaagtgaattggaatatggtaactgcaacaccaagtgtcaaactccaatgggggcgataaactctagtatgccattccacaacatacaccccctcaccatcggggaatgccccaaatatgtgaaatcaaacagattagtccttgcgactggactcagaaatacccctggactatttggagctatagcaggttttatagagggaggatggcagggaatggtagatggttggtatgggtaccaccatagcaatgagcaggggagtggatacgctgcagacaaagaatccactcaaaaggcaatagatggagtcaccaataaggtcaactcgatcattgacaaaatgaacactcagtttgaggccgttggaagggaatttaataacttggaaaggaggatagagaatttaaacaagcagatggaagacggattcctagatgtctggacttataatgctgaacttctggttctcatggaaaatgagagaactctagactttcatgactcaaatgtcaagaacctttatgacaaggtccgactacagcttagggataatgcaaaggagctgggtaatggttgtttcgagttctatcacaaatgtgataatgaatgtatggaaagtgtaaaaaacggaacgtatgactacccgcagtattcagaagaagcaagactaaacagagaggaaataagtggagtaaaattggaatcaatgggaacttaccaaatactgtcaatttattcaacagtggcgagttccctagcactggcaatcatggtagct"),
  
  # button
  
  withBusyIndicatorUI(
    actionButton(inputId = "click", 
               label = "Enter")
  ),
  
  # textoutput
  tableOutput("blastinfo"),
  
  # plotoutput
  plotOutput("tree")
  
  
  )



# server ####
  
server <- function(input, output){
  
  # data reactive to click
  c_data <- eventReactive(input$click, {
    
    withBusyIndicatorServer("click", {
      
      blastinfo(blastseq(input$seq))
      
    })
    
     }) 

  # rendertable  
  output$blastinfo <- renderTable({
    
    refSequence = treedata_subTree[,10][c_data()] })
    
  # renderplot
  output$tree <- renderPlot({
    
    treedata_subTree[,11][c_data()] <- "orange"
    H5_tip$shape[c_data()] <- 1
    
    p1 <- ggtree(subTreefile) %<+% treedata_subTree + 
      aes(color = I(hit_color), alpha = 0.5)
    
    p1 %<+% H5_tip + geom_tippoint(aes(shape = factor(shape)), size = 2.5)
    
  })
  
          
  
}


shinyApp(ui = ui, server = server)

