library(seqinr) 

# ID cleaning ####

g_20170112_8135 <- cleanID()

# remove _A_/
         file = read.fasta(file.choose())
    seq_name0 = attributes(file)$names
         seq0 = getSequence(file)
         
    seq_name = sub("_A_/", "",seq_name0)
    
    write.fasta(seq0, 
                file.out = "~/Desktop/cleanID2.fasta", 
                names = seq_name)

n_20170112_3964 <- cleanID()


# seq curation 

# n = 6840
g_20170112_8135 <- curateSeq(5, 1500, 1)

# n = 3413
n_20170112_3964 <- curateSeq(5, 1500, 1)


pool_10653 <- curateSeq(5, 1500, 0)


# clade define file (from Dr. Gavin Smith)
# n = 238 

         file = read.fasta(file.choose())
    seq_name0 = attributes(file)$names
         seq0 = getSequence(file)

    seq_name = sub("\\{new\\}", "",seq_name0)
    
    seq_name = gsub(" ", "_", seq_name)
    seq_name = gsub("\\(", "-", seq_name)
    seq_name = gsub("\\)", "-", seq_name)
    seq_name = gsub("\\[", "-", seq_name)
    seq_name = gsub("\\]", "-", seq_name)
    seq_name = gsub("\\'", "", seq_name)
    seq_name = gsub(">", "", seq_name)  
    seq_name = gsub("\\.", "-", seq_name) 
    seq_name = gsub("\\|", "_", seq_name)  
    seq_name = gsub("\\{", "_", seq_name)  
    seq_name = gsub("\\}", "", seq_name)  
    
    write.fasta(seq0, 
                file.out = "~/Desktop/smallref_cleanID.fasta", 
                names = seq_name)
    
# align and trim file: pool_7982ref = pool_10653_curateSeq + smallref_cleanID
# another round of cleanID and curate
# resulting n = 7229    
    
    pool_7982ref_trim = cleanID()
    pool_7982ref_trim = curateSeq(99, 99, 0)
    

# read-in tree ####
    
    library(ape)
    library(ggtree)
    
    sub_7225 <- read.tree(file.choose())
    ggtree(sub_7225)
    
    
# function needed ####
    
    parsexmltoDF <- function(xml) {
      
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
    parseResult <- function(baseUrl, rid, rtoe) {
      
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
    
    x = "tgtcaaatcagataaactggtccttgcaacaggactgaggaacgtgcctggtctgtttggagcaatagcaggattcatagaaggggggtggcaaggaatggtagatggatggtatggttaccatcatagcaacgagcagggaagtg"
    y = blastseq(x)

# tree annotation ####
    
    
    nn <- "A/([A-Za-z0-9-_]+)/([A-Za-z0-9-_]+)/([A-Za-z0-9-_]+)"
    hit_id <- str_match(y, nn)[,1]
    
    treedata_sub7225 <- fortify(sub_7225)
    
    treedata_sub7225[,10] <- gsub("'", "", treedata_sub7225$label)
    
    # find out all hit on tree
    hitlabel = c()
    for(i in 1: 10){
      
      add = grep(hit_id[i], treedata_sub7225[,10])
      hitlabel = c(hitlabel, add)
      
    }

    hitlabel = unique(hitlabel)
    
    # coloring the the tip on tree file
              treedata_sub7225[,11] <- "grey"
    treedata_sub7225[,11][hitlabel] <- "orange"
    
      colnames(treedata_sub7225)[11] = "hit_color"
    
    # coloring the tip on tip file
      
    H5_tip = data.frame(node= c(1:length(sub_7225$tip.label)), shape = NA)  
    H5_tip$shape[hitlabel] = 1
    
    # tree
    p1 <- ggtree(sub_7225) %<+% treedata_sub7225 + aes(color = I(hit_color), alpha = 0.5)
    p1 %<+% H5_tip + geom_tippoint(aes(shape = factor(shape)), size = 2.5)
    
    
    
    
    
    
    
    