### fastaEx -------------------------------- 

fastaEx <- function(filedir = file.choose())
{
  require(seqinr)
  
  file     <- read.fasta(filedir)
  file_seq <- getSequence(file)
  file_id  <- attributes(file)$names
  
  return( list(seq = file_seq, 
               id  = file_id ) )
  #v201706
}





### idInfo --------------------------------
idInfo <- function( rawid, 
                    datasource = "n",
                    g.csv      = "")
{
  # format: 
  # N:  >{accession}_{strain}_{serotype}_|{country}|_{year}-{month}-{day}
  # G:  Isolate name Type Collection date Isolate ID
  # Both need replace the blank with underline 
  
  library(seqinr)
  library(stringr)
  
  # g
  a.string.g  <- "EPI_ISL_([0-9]+)"
  s.string.g  <- "_A_/_(H[5,9]N[0-9xX]{1,2})_"
  y.string.g  <- "_[0-9]{4}[-0-9]{6}|_[0-9]{4}-[0-9]{2}_\\(Day_unknown\\)|_[0-9]{4}_\\(Month_and_day_unknown\\)" 
  n.Nstring.g <- "_EPI_ISL_([0-9]+)|_A_/_H[5,9]N[0-9xX]{1,2}|_[0-9]{4}[-0-9]{6}|_[0-9]{4}-[0-9]{2}_\\(Day_unknown\\)|_[0-9]{4}_\\(Month_and_day_unknown\\)"
  
  # n 
  a.string.n   <- "[A-Z]{1,2}[0-9]{5,6}"
  s.g.string.n <- "_(H[5,9][N0-9xX]{0,2})_\\|([a-zA-Z_\\']+)\\|"
  y.string.n   <- "_[0-9]{4}-[0-9]{2}-[0-9]{2}|_[0-9]{4}-[0-9]{2}-|_[0-9]{4}--|_--"
  n.Nstring.n  <- "[A-Z]{1,2}[0-9]{5,6}_|_(H[5,9][NxX0-9]{0,2})_\\|([a-zA-Z_\\']+)\\||_([0-9]{4}[-0-9]{2,6})$|_--$"
  
  if( datasource == "g")
  {
    id.a <- gsub("_ISL_", "", str_match( rawid, a.string.g )[, 1] )
    
    id.s <- str_match( rawid, s.string.g )[,2] 
    id.s[ which(id.s == "H5"|
                  id.s == "H5N"|
                  id.s == "H5Nx"|
                  id.s == "H5NX")  ] = "H5N0"
    id.s[ which(id.s == "H9"|
                  id.s == "H9Nx"|
                  id.s == "H9NX"|
                  id.s == "H9N")  ] = "H9N0"
    
    id.y <- str_match( rawid, y.string.g )
    id.y <- gsub( "^_", "", x = id.y)[,1]
    
    id.n <- gsub( n.Nstring.g, rawid, replacement = "")
    
    id.n[ which( startsWith(id.n, "A/") == FALSE) ] <- gsub( "_A/", "A/", id.n[ which( startsWith(id.n, "A/") == FALSE) ] )
    
    id.n  <- gsub("\\(|\\)|\\[|\\]|\\.|:|-|/", "_", id.n)
    id.n  <- gsub("\\'|\\?|>", "", id.n)
    id.n  <- gsub("A_", "", id.n)
    id.n  <- gsub("_$", "", id.n)
    id.n  <- gsub("__", "_", id.n)
    
    g <- gsub( " ", "_", read.csv( g.csv, header = TRUE, stringsAsFactors = FALSE)$Location )
    g <- gsub("_$", "",  str_match( g, "([A-Za-z_]+)_/_([A-Za-z_]+)" )[,3] ) 
    
    g[ which( is.na(g) == TRUE ) ] = "Unknown"
    g[ which(g == "Russian_Federation") ] = "Russia"
    g[ which(g == "United_States") ] = "USA"
    g[ which(g == "Korea") ] = "South_Korea"
    
    id.g <- g[ match( id.a, gsub("_ISL_", "", read.csv( g.csv, header = TRUE, stringsAsFactors = FALSE)$Isolate_Id ) ) ]
    
    
  }else
  {
    id.a <- str_match( rawid, a.string.n)[,1]
    
    id.s <- str_match( rawid, s.g.string.n)[,2]
    id.s[ which(id.s == "H5"|
                id.s == "H5N"|
                id.s == "H5Nx"|
                id.s == "H5NX")  ] = "H5N0"
    id.s[ which(id.s == "H9"|
                id.s == "H9Nx"|
                id.s == "H9NX"|
                id.s == "H9N")  ] = "H9N0"
    
    id.g <- str_match( rawid, s.g.string.n)[,3]
    id.g[ which( id.g == "Viet_Nam") ] = "Vietnam"
    id.g[ which( id.g == "Cote_d'Ivoire") ] = "Cote_dIvoire"
    
    id.y <- str_match( string = rawid, y.string.n)
    id.y <- gsub( "_--", "1900-01-01", id.y)
    id.y <- gsub( "^_", "", id.y)
    
    id.n <- gsub( n.Nstring.n, "", rawid)
    
    id.n[ which( startsWith(id.n, "A/") == FALSE) ] <- 
      paste0("A/", id.n[ which( startsWith(id.n, "A/") == FALSE) ])
    
    id.n  <- gsub("\\(|\\)|\\[|\\]|\\.|:|-|/|__", "_", id.n)
    id.n  <- gsub("\\'|\\?|>", "", id.n)
    id.n  <- gsub("A_", "", id.n)
    id.n  <- gsub("_$", "", id.n)
    id.n  <- gsub("__", "_", id.n)
    
  }
  
  infolist = list(id.a, id.s, id.g, id.y, id.n)
  
  e = 
    which(
      sapply( infolist, 
              function(x)
              {
                TRUE %in% is.na(x)
                
              })  == TRUE )
  
  
  print( paste("ERROR in ", c("ac", "sero", "geo", "year", "name")[e] )  )
  
  return(infolist)
  
  #v20170920b
}

### strainSelect --------------------------------

strainSelect <- function( infolist )
{
  
  infolist.n <- infolist[[ length(infolist) - 1 ]]
  infolist.q <- infolist[[ length(infolist) ]]
  infolist.y <- infolist[[ length(infolist) - 2 ]]
  infolist.a <- infolist[[1]]
  
  toberemove <- c() 
  dup        <- which( duplicated( infolist.n ) )
  
  for(i in 1: length(dup) )
  {
    id_dup_ii <- which( infolist.n %in% infolist.n[ dup[i] ] == TRUE )
    lth_ii    <- sapply( infolist.q[id_dup_ii],
                         
                         function(x)
                         {
                           y = c2s(x)
                           z = gsub("-|~", "", y)
                           z = grep( pattern = "a|t|c|g", 
                                     x = y, 
                                     ignore.case = TRUE, value = TRUE )
                           
                           l = length( s2c(z) )
                           
                           return(l)
                           
                         } )
    
    SeqL      <- which.max( lth_ii ) 
    
    if ( length(  which( lth_ii == max(lth_ii) )  ) > 1 )
    {
      
      id_dup_jj  <- id_dup_ii[ which( lth_ii == max(lth_ii) ) ] 
      
      nchar_jj   <- nchar( gsub( pattern     = "[-\\(\\)A-Za-z]+", 
                                 replacement = "",
                                 x           = infolist.y[id_dup_jj] ) )
      
      id_dup_j   <- which.max( nchar_jj )
      SeqL       <- which( lth_ii == max(lth_ii) )[id_dup_j]
      
      
      if (  length( which( nchar_jj == max(nchar_jj) ) ) > 1  )
      {
        
        id_dup_kk <- id_dup_jj[ which(nchar_jj == max(nchar_jj) ) ]
        
        ac        <- infolist.a[id_dup_kk]
        ac.a      <- nchar( gsub( pattern = "[0-9]+", replacement = "", x = ac) )
        ac.d      <- as.numeric( gsub( pattern = "[a-zA-Z]+", replacement = "", x = ac ) )
        ac.df     <- data.frame( id_dup_kk, ac.a, ac.d )
        
        SeqL      <- which( id_dup_ii == ac.df[order( ac.df[,2], ac.df[,3] ),][1,1] )
        
      }
    }
    
    toberemove = c( toberemove, id_dup_ii[-SeqL] )
    
  }
  
  remain <- seq( 1, length(infolist.q) )[- toberemove]
  
  newlist = list()
  for(l in 1 : length(infolist) )
  {
    newlist[[l]] <- infolist[[l]][remain]
    
  }
  
  newlist[[ length(newlist) + 1 ]] <- ifelse( grepl( pattern = "1900-01-01|--$|Month", x = newlist[[4]] ), 1, 0)
  
  if( TRUE %in% is.na( unlist(newlist) ) ){ print("ERROR") }
  
  return( newlist )
  
  #v20170920b
}


### seqDate ----------------------------------

seqDate <- function( rawdata )
{
  library(stringr)
  
  # gisaid
  
  rawdata.1 <- gsub( "_\\(Day_unknown\\)", "-15", 
                     gsub( "_\\(Month_and_day_unknown\\)", "-07-01", rawdata ) )
  
  # ncbi
  
  rawdata.2 <- gsub( "-$", "-15", 
                     gsub( "--$", "-07-01", rawdata.1) )
  
  # parse into numeric
  
  d = "([0-9]{4})-([0-9]{2})-([0-9]{2})"
  
  yr   <- as.numeric( str_match(rawdata.2, d)[,2] )
  yr.0 <- paste0(yr, "-01-01")
  
  daydifference <- as.numeric( difftime( strptime( rawdata.2, "%Y-%m-%d"),
                                         strptime( yr.0, "%Y-%m-%d"), 
                                         units = "days") 
  )/365
  
  # bug?
  if ( TRUE %in% is.na(daydifference) )
  {
    
    rawdata.2[ which(is.na(daydifference)) ] <- 
      sub("01$", "02", rawdata.2[ which(is.na(daydifference)) ] )
    
    
    daydifference <- as.numeric( difftime( strptime( rawdata.2, "%Y-%m-%d"),
                                           strptime( yr.0, "%Y-%m-%d"), 
                                           units = "days") 
    )/365
  }
  
  yr.daydifference <- yr + daydifference
  yr.daydifference <- format( round( yr.daydifference, 3 ), nsmall = 3)
  
  return(yr.daydifference)
  
  #v20170921b
}


### seqSelect --------------------------------

seqSelect <- function( minlth  = 1000, 
                       maxamb  = 1,
                       rmdup   = TRUE,
                       seqlist )
{
  df   <- data.frame( a = seqlist[[1]], 
                      s = seqlist[[2]],
                      y = seqlist[[4]],
                      i = seqlist[[7]], stringsAsFactors = FALSE) 
  
  df.s <- df[ order(df$i, df$y, df$a, df$s), ]
  
  idx  <- as.numeric( rownames(df.s) )
  q    <- seqlist[[6]][ idx ]
  
  
  # length and ambiguous nucleotide
  lth_amb <- which( sapply( q, 
                            
                            function(x)
                            {
                              ATCG  <-  c("a", "t", "c", "g")
                              
                              x.s   <- gsub( "-|~", "", c2s( x ) )
                              x.l   <- length( s2c(x.s) )
                              
                              x.c.s <- grep( "a|t|c|g", s2c(x.s) )[1]
                              x.c.e <- grep( "a|t|c|g", s2c(x.s) )[ length( grep( "a|t|c|g", s2c(x.s) ) ) ]
                              
                              x.a   <- length( which(! s2c(x.s)[x.c.s: x.c.e] %in% ATCG ) )
                              
                              return( x.l < minlth | x.a > (maxamb/100)*x.l )
                              
                            } ) ) 
  # duplicated sequence
  if (rmdup)
  {
    dup     <- which( duplicated( sapply( q,  
                                          function(x)
                                          {
                                            x.s <- gsub( "~|-", "", c2s(x) )
                                            return(x.s)
                                          }) 
    ) )
    
  }else{
    dup = c()
  }
  
  if ( ( length(dup) + length(lth_amb) ) > 0 )
  {
    remain  <- seq(1, length( seqlist[[6]] ) )[ - unique( sort( c(dup, lth_amb) )) ]
    
  }else
  {
    remain  <- seq(1, length( seqlist[[6]] ) )  
  }
  
  
  newlist = list()
  for(l in 1 : length(seqlist) )
  {
    newlist[[l]] <- seqlist[[l]][idx][remain]
  }
  
  print( paste0("n = ", length( newlist[[1]] )) )
  
  return(newlist)
  
  #v20170921b
}


### trimtool --------------------------------

trimtool <- function( propblank = 0.8, 
                      filedir   = file.choose()){
  
  library(stringr)
  library(seqinr)
  
  file       = read.fasta(filedir)
  seq_name0  = attributes(file)$names
  seq0       = getSequence(file)
  seq_matrix = do.call(rbind, seq0)
  
  
  coltoberemove = apply(seq_matrix, 2, 
                        function(x)
                        {
                          blank = ( length( which( x == "-") ) + length( which( x == "~") ) ) 
                          fl    = length(x)
                          
                          if ( fl*propblank < blank ){ return(1) }else{ return(0) }    
                        }
  )
  
  cut_matrix = seq_matrix[ ,-which(coltoberemove == 1) ]
  
  seq_cut    = as.list( data.frame(t(cut_matrix), stringsAsFactors = FALSE) )
  
  filename <- str_match(filedir, "([a-zA-Z0-9_-]+)(\\.)(fas)" )[,2]
  
  write.fasta( seq_cut, 
               file.out = paste0("./trim_", filename, ".fasta"),
               names    = seq_name0)
  print("DONE")
  
  #20170924  
}


### rmDup --------------------------------

rmDup <- function( fasfile = file.choose(), 
                   year    = c(1000,3000),
                   geo     = c(),
                   sero    = "",
                   rmdup   = TRUE)
{
  require(seqinr)
  require(stringr)
  
  readin <- read.fasta( fasfile )
  seq    <- getSequence( readin )
  id     <- attributes( readin )$names
  
  if( rmdup )
  {
    
    # order: time ( data completeness ), accession number ( data source )
    
    id.y  <- as.numeric( str_match( id, "_([0-9]{4}\\.[0-9]+$)")[,2] )
    id.a  <- str_match( id, "EPI[0-9]+|[A-Z]{1,2}[0-9]{5,6}" )[,1]
    
    id.d  <- 
      ifelse( ( endsWith( str_match( id, "_([0-9]{4}\\.[0-9]+$)")[,2], ".496" )|
                  endsWith( str_match( id, "_([0-9]{4}\\.[0-9]+$)")[,2], ".499" ) )
              , 1, 0)
    
    id.a.c <- nchar( gsub("[0-9]+", "", id.a) )
    id.a.d <- as.numeric( gsub("[A-Za-z]+", "", id.a) )
    
    
    df <- data.frame( id.y, id.d, id.a.c, id.a.d )  
    df <- df[ order( df[,2], df[,1], df[,3], df[,4] ), ]
    
    seq <- seq[ as.numeric( rownames(df) ) ]
    id  <- id[ as.numeric( rownames(df) ) ]
    
    
    dup <- which( duplicated( sapply( seq, 
                                      function(x)
                                      {
                                        x.s <- gsub( "~|-", "", c2s(x) )
                                        return(x.s)
                                      } ) 
    ) )
    
  }else
  {
    dup = NA
  } 
  
  if( length(dup) > 1 )
  {
    remain = seq( 1, length(seq) )[ - sort( unique(dup) ) ]
    
  }else
  {
    remain = seq( 1, length(seq) )
    
  }
  
  # year
  y = which( (str_match( id, "_([0-9]{4}\\.[0-9]+$)")[,2] > year[1] & 
                str_match( id, "_([0-9]{4}\\.[0-9]+$)")[,2] < year[2]) )
  
  # geo
  geo.p <- paste0( geo, collapse = "|" )
  g     <- grep( geo.p, str_match( id, "\\|([A-Za-z_]+)\\|" )[,2] )
  
  # sero
  
  s     <- grep( sero, str_match( id, "_(H[5,9]N[0-9]{1,2})_" )[,2] )
  
  if ( TRUE %in% is.na( 
    c(str_match( id, "_([0-9]{4}\\.[0-9]+$)")[,2], 
      str_match( id, "\\|([A-Za-z_]+)\\|" )[,2], 
      str_match( id, "_(H[5,9]N[0-9]{1,2})_" )[,2]) ) )
  {
    stop( 
      c("Year", "Geo", "Serotype")[ ceiling( 
        which( is.na( 
          c(str_match( id, "_([0-9]{4}\\.[0-9]+$)")[,2], 
            str_match( id, "\\|([A-Za-z_]+)\\|" )[,2], 
            str_match( id, "_(H[5,9]N[0-9]{1,2})_" )[,2]) ))/3) ] 
    )
  }  
  
  remain <- sort( Reduce( intersect, list( remain, y, g, s ) ) )
  
  write.fasta( seq[ remain ], 
               id[remain], 
               file.out = sub(".fasta", "_cr.fasta", fasfile) )
  
  print( length( remain ) )
  #v20170927v
}


### taxainfo --------------------------------

taxaInfo <- function( file    = file.choose(), 
                      useTree = FALSE, 
                      makecsv = FALSE )
{
  # input: 
  # 1 colored .tre file
  # 2 .fas file with clean id 
  
  require(seqinr)
  require(stringr)
  
  if ( useTree )
  {
    anno.tre <- read.csv( file, stringsAsFactors = FALSE)
    taxa.s   <- grep( "taxlabels", anno.tre[,1] ) + 1
    
    ntax     <- as.numeric( str_match( grep( "ntax", anno.tre[,1],  value = TRUE ), 
                                       "(ntax=)([0-9]+)" )[,3] )
    taxa.e   <- taxa.s + ntax - 1
    
    id  <- str_match( anno.tre[, 1][taxa.s: taxa.e], "\'([0-9A-Za-z_\\|.]+)\'" )[,2]
    tag <- str_match( string = anno.tre[, 1][taxa.s: taxa.e], 
                      pattern = "color=#([a-z0-9]{6})")[, 2]
    
    id.a <- str_match( id, "[A-Z]{1,2}[0-9]{5,6}|EPI[0-9]+" )[,1]
    id.g <- str_match( id, "\\|([A-Za-z_]+)\\|")[,2]
    id.s <- str_match( id, "_(H[5,9]N[0-9]{1,2})_")[,2]
    id.y <- as.numeric( str_match( id, "_([0-9]{4}.[0-9]{3})$")[,2] )  
    id.n <- gsub("^[A-Z]{1,2}[0-9]{5,6}_|^EPI[0-9]+_|_\\|[A-Za-z_]+\\|_|H5N[0-9]{1,2}_[0-9]{4}.[0-9]{3}$", "", id)
    
    
    ls <- list( id.a, id.g, id.s, id.y, id.n, id, tag)
    df <- data.frame( id.a, id.g, id.s, id.y, id.n, id, tag )
    
    
    if( TRUE %in% is.na(unlist( ls[c(1:6)] ) ) ){ stop() }
    
    
  }else
  {
    fas <- read.fasta( file )
    id  <- attributes( fas )$names
    seq <- getSequence( fas )
    
    id.a <- str_match( id, "[A-Z]{1,2}[0-9]{5,6}|EPI[0-9]+" )[,1]
    id.g <- str_match( id, "\\|([A-Za-z_]+)\\|")[,2]
    id.s <- str_match( id, "_(H[5,9]N[0-9]{1,2})_")[,2]
    id.y <- as.numeric( str_match( id, "_([0-9]{4}.[0-9]{3})$")[,2] )  
    id.n <- gsub("^[A-Z]{1,2}[0-9]{5,6}_|^EPI[0-9]+_|_\\|[A-Za-z_]+\\|_|H5N[0-9]{1,2}_[0-9]{4}.[0-9]{3}$", "", id)
    
    seq.l <- sapply( seq,
                     function(x)
                     {
                       y = c2s(x)
                       z = gsub("-|~", "", y)
                       z = grep( pattern = "a|t|c|g", 
                                 x = y, 
                                 ignore.case = TRUE, value = TRUE )
                       l = length( s2c(z) )
                       
                       return(l)
                     } )
    
    ls <- list( id.a, id.g, id.s, id.y, id.n, id, seq.l)
    df <- data.frame( id.a, id.g, id.s, id.y, id.n, id, seq.l )
    
    if( TRUE %in% is.na(unlist( ls[c(1:6)] ) ) ){ stop() }
    
  }
  
  if ( makecsv ){ write.csv(df, file = sub( ".fasta", "_info.csv", file) , row.names = FALSE) }
  
  return( ls )
  
  #v20170928v
}


### findtaxa --------------------------------

findtaxa <- function(type, 
                     tree, 
                     targetid, 
                     target)
{
  
  # type 1 = branch coloring 
  # type 0 = tip shape
  # default branch color = black
  
  library(ape)
  library(ggtree)
  
  # extract tree data
  
  tree.d                         <- fortify(tree)
  tree.d[, ncol(tree.d) + 1]     <- gsub("'", "", tree.d$label)
  colnames(tree.d)[ncol(tree.d)] <- "taxaid"
  
  # for tip shape
  
  if (type == 0)
  {
    
    shapetaxa <- data.frame(node = c(1:length(tree.d$isTip)), shapee = NA)
    
    for (i in 1: length(targetid))
    {
      shapetaxa$shapee[  grep(  tolower( targetid[i] ), tolower(tree.d$taxaid) ) ] <- target[i]
      
    }
    
    return(shapetaxa)
    
  }else {
    
    # for branch colorring  
    
    # new column
    
    tree.d[, ncol(tree.d) + 1]     <- "black"
    colnames(tree.d)[ncol(tree.d)] <- "colorr"
    
    # for branch extension
    
    edgematrix <- as.matrix(tree.d[,c(2,1)])
    
    # color grouping 
    
    group_color <- unique(target)
    
    for (i in 1: length(group_color) )
    {
      
      # color as group to combine key word to targetno
      
      sub_color <- which(target == group_color[i] )
      targetno  <- c()
      
      for (t in 1: length(sub_color) )
      {
        
        targetno <- 
          unique( c(targetno, grep( tolower( targetid[ sub_color[t] ] ), tolower(tree.d$taxaid) )) )
        
      }
      
      tobecolor     <- c()
      pre_targetno  <- length(targetno)
      post_targetno = 0
      
      # while loop 
      
      while( pre_targetno != post_targetno )
      {
        
        pre_targetno = length(targetno)
        
        for(k in 1:length(targetno))
        {
          
          # all sibiling 
          sibs <- edgematrix[
            which(edgematrix[,1] == 
                    edgematrix[which(edgematrix[,2] == targetno[k]),][1]),][,2]
          
          if (length(sibs) == 1)
          {
            
            targetno = c(targetno, edgematrix[which(edgematrix[,2] == targetno[k]),][1])
            
          }else{
            
            if (length(which(sibs %in% targetno == "FALSE")) == 0){
              
              tobecolor = c(edgematrix[which(edgematrix[,2] == targetno[k]),][1], tobecolor)
              targetno  = c(targetno, edgematrix[which(edgematrix[,2] == targetno[k]),][1])
            }
            
          }
          targetno  = unique(targetno)
          tobecolor = unique(c(targetno, tobecolor))
          
        }
        
        post_targetno = length(targetno)
        
      }
      
      # coloring
      
      tree.d$colorr[tobecolor] <- group_color[i]
      
    }
    return(tree.d)    
    
  }
  #v201706
}


### cladesampling --------------------------------




cladeSampling <- function( trefile      = file.choose(),
                           fasfile      = file.choose(),
                           listinput    = list(),
                           seed         = 666,
                           grid         = 1,
                           minBranchlth = TRUE, 
                           showTree     = FALSE, 
                           saveFasta    = FALSE,  
                           suppList     = FALSE,
                           list.x       = c("id", "y", "geo") )
{
  require( ape )
  require( seqinr )
  require( ggtree )
  require( stringr )
  
  # getdescendants
  getDes <- function( node, curr = NULL )
  {
    if( is.null(curr) ){ curr <- vector() }
    
    edgemax   <- tre.d[ c(2,1) ]
    daughters <- edgemax[which( edgemax[,1] == node ), 2]
    
    curr <- c(curr, daughters)
    nd   <- which( daughters >= length( which( tre.d$isTip )) )
    
    if( length(nd) > 0)
    {
      for(i in 1:length(nd) ){ curr <- getDes( daughters[ nd[i] ], curr ) }
    }
    return(curr)
  }
  
  
  nex <- read.nexus( trefile )
  fas <- read.fasta( fasfile )
  seq <- getSequence( fas )
  id  <- attributes( fas )$names
  
  tre.d   <- fortify( nex )
  N.tip   <- length( which( tre.d$isTip ) )
  N.node  <- nex$Nnode
  edgemax <- tre.d[ c(2,1) ]
  
  
  t.id <- gsub("'", "", tre.d$label)
  
  if( suppList )
  {
    
    tem.m  <- match( t.id[ 1:  N.tip], listinput[[ list.x[1] ]])
    
    if( TRUE %in% is.na(tem.m) ){stop()}
    
    i.id.y <- listinput[[ list.x[2] ]][ tem.m ]
    i.id.g <- listinput[[ list.x[3] ]][ tem.m ]
    
  }else
  {
    i.id.y <- as.numeric( str_match( t.id, "_([0-9]{4}\\.[0-9]+)$")[,2] )
    i.id.g <- str_match( t.id, "\\|([A-Za-z_]+)\\|")[,2]  
  }
  
  # 1st search for node with homogeneous descendants 
  inner.node <- seq( 1, dim(tre.d )[1])[ - seq(1, N.tip+1) ]
  c1.node    <- inner.node[ which( sapply( as.list(inner.node), 
                                           function(x)
                                           {
                                             all <- getDes(x)[ getDes(x) <= N.tip ]
                                             r   <- range( i.id.y[all])[2] - range(i.id.y[all] )[1]
                                             g   <- unique( i.id.g[all] )
                                             return( (r <= grid) & ( length( g ) == 1 ) )
                                           } )) 
                            ]
  # reduce redndant nodes
  c2.node    <- c1.node[ which( sapply( as.list(c1.node), 
                                        function(x)
                                        {
                                          if( edgemax[,1][ which( edgemax[,2] == x) ] %in% c1.node )
                                          {
                                            return( FALSE )
                                            
                                          }else
                                          {
                                            return( TRUE )
                                          }
                                          
                                        } )
  )  
  ] 
  
  c2.node_des <- c( c2.node, unlist( sapply( as.list(c2.node), getDes) ) )
  c2.node_tip <- c2.node_des[ c2.node_des <= N.tip ]
  nogroup_tip <- seq(1, N.tip)[ - c2.node_tip ]
  
  # sample within a group
  if( minBranchlth )
  {
    selected_tip <- sapply( as.list(c2.node), 
                            function(x)
                            {
                              alltip <- getDes(x)[ getDes(x) <= N.tip ]
                              m      <- which.min( tre.d$x[ alltip ] )
                              
                              if( length( which( tre.d$x[ alltip ] == tre.d$x[ alltip ][m] ) )  > 1 )
                              {
                                set.seed( seed ) 
                                s = sample( which( tre.d$x[ alltip ] == tre.d$x[ alltip ][m] ), 1 )
                                
                              }else
                              {
                                s = m
                              }
                              return( alltip[s] )
                            })
    
  }else
  {
    selected_tip <- sapply( as.list(c2.node), 
                            function(x)
                            {
                              set.seed( seed ) 
                              s = sample( getDes(x)[ getDes(x) <= N.tip ], 1)
                              return(s)
                              
                            } 
    )
  }  
  
  if( showTree )
  {
    # view the result
    tre.d[, ncol(tre.d) + 1 ] = "gray"
    colnames(tre.d)[ ncol(tre.d) ] = "colorr"
    tre.d$colorr[c2.node_des] = "red"
    
    tre.d[, ncol(tre.d) + 1 ] = NA
    colnames(tre.d)[ ncol(tre.d) ] = "shapee"
    tre.d$shapee[selected_tip] = 16
    
    g1 <- ggtree( nwk ) %<+% tre.d + aes(color = I(colorr)) + geom_tiplab(size = 1)
    g1 + geom_tippoint(aes( shape = factor(shapee) ), size = 2)
  }
  
  remain <- c( nogroup_tip, selected_tip )
  seq.o  <- seq[ match( t.id[ sort(remain) ], id) ]
  id.o   <- id[ match( t.id[ sort(remain) ], id) ]
  
  if( saveFasta )
  {
    write.fasta( seq.o, id.o, 
                 file.out = gsub( ".fasta", "_s.fasta", fasfile) )
  }
  
  print( paste0("sampled n = ", length(remain), " from ", length(id) ) )
  print( table( floor( i.id.y[remain] ), i.id.g[remain]) )
  
  
  #v20171005b
}
