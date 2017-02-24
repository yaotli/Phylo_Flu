# getDescendants (from http://blog.phytools.org; 20120126) #####


getDescendants<-function(tree,node,curr=NULL){
  
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  curr<-c(curr,daughters)
  w<-which(daughters>=length(tree$tip))
  if(length(w)>0) for(i in 1:length(w)) 
    curr<-getDescendants(tree,daughters[w[i]],curr)
  ll = length(tree$tip.label)
  curr<-curr[which(curr <= ll)]
  
  return(curr) }



# Clean tree IDs ####

# should follow strict ID arrangement before download .fasta file
# require package seqinr and stringr
# deal with duplicated, problematic string and extract info
# fill the date info with either -15 or -99-99


cleanID <- function(){
  
  library(seqinr)
  library(stringr)
  
  file = read.fasta(file.choose())
  
  seq_name0 = attributes(file)$names
       seq0 = getSequence(file)
  
  
  # Deal with problematic string
  
  seq_name = gsub(" ", "_", seq_name0)
  
  seq_name = gsub("\\(", "-", seq_name)
  seq_name = gsub("\\)", "-", seq_name)
  seq_name = gsub("\\[", "-", seq_name)
  seq_name = gsub("\\]", "-", seq_name)
  
  seq_name = gsub("\\'", "", seq_name)
  seq_name = gsub(">", "", seq_name)  
  seq_name = gsub("\\.", "-", seq_name)  
  
  
  # Dissect the string 
  # for GISAID
  
  seq_name = gsub("_-Month_and_day_unknown-", "-99-99", seq_name)
  seq_name = gsub("_-Day_unknown-", "-15", seq_name)
  
  seq_name[which( endsWith(seq_name, "-") == "TRUE")] <- 
    paste0(seq_name[which( endsWith(seq_name, "-") == "TRUE")], "15")
  
  
  # Time
  
  d = "([0-9]{4})-([0-9]{2})-([0-9]{2})"
  
  time_raw <- str_match(seq_name, d)[,1]
  
  # Site of isolate
  
  s = "A/([A-Za-z0-9-_]+)/"
  
  site <- str_match(seq_name, s)[,2]
  
  
  # Deal with replicated
  
    duplicated_id = which(duplicated(seq_name) == "TRUE")
  duplicated_note = duplicated(seq_name)
  
  # loop for each replicated case
  
  if( length(duplicated_id) > 0 ){
    
    # for all the duplicated id
    
    for (i in 1: length(duplicated_id)){
      
      dup0 = which(match(seq_name, seq_name[duplicated_id[i]]) != "NA")
      
      # create a null vector to appendex  
      
         app = c("", paste0("_",letters), toupper(paste0("_",letters)))
      
       ap.id = seq_name[dup0]
      app.id = c()
      
      # loop to deal with multiple replicated
      # find the date info, insert labeling in the middle
      
      time_rep <- str_match(ap.id, d)[,1]
      
      for (k in 1: length(ap.id)){
        
        app.id[k] = sub(time_rep[k], paste0(app[k], time_rep[k]), ap.id[k])
      }
      
      # back to seq_name0  
      
      seq_name[dup0] = app.id
      
    }
  }     
  
  duplicated_id_ed = which(duplicated(seq_name) == "TRUE")
  
  
  if (length(duplicated_id_ed) > 0 ){
    
print("ERROR") 
    
  }else{
    
    # write fasta file
    
    write.fasta(seq0, 
                file.out = "~/Desktop/cleanID.fasta", 
                names = seq_name)
    
    # dataframe output
    
    identicalID = as.character(duplicated_note)
    
    fastaInfo <- data.frame(no = seq(1:length(seq_name)),seq_name, site, time_raw, identicalID)
    
    
return(fastaInfo)
print("DONE")
    
  }
  
}

# Subtree seq extraction #####


subtreseq<-function(){
  
  library(seqinr) 
  
  fasta0 = read.fasta(file.choose())
  seq.name0 = attributes(fasta0)$names
  seq0 = getSequence(fasta0)
  
  sub.tree = read.csv(file.choose(), header = FALSE, stringsAsFactors = FALSE)
  id.subtree = match(sub.tree[,1], seq.name0)
  
  if( any(NA %in% id.subtree) == "TRUE"){
    
    print("Mismatch")
    
  }else{
    
    seq.name0_subtree = seq.name0[id.subtree]
    seq0_subtree = seq0[id.subtree]
    
    write.fasta(seq0_subtree, file.out = "~/Desktop/subtree.fasta", names = seq.name0_subtree)
    
    print("Done")
    
  } }

subtreseq2<-function(){
  
  library(seqinr) 
  
  fasta0 = read.fasta(file.choose())
  seq.name0 = attributes(fasta0)$names
  seq0 = getSequence(fasta0)
  bug.id = grep("'", seq.name0)
  
  sub.tree=read.csv(file.choose())
  
  tipslabel.subtree0 = as.character(sub.tree[,1])
  startno = which(tipslabel.subtree0 == "\ttaxlabels") + 1
  endno = head(which(tipslabel.subtree0 == ";"), 1) -1
  
  tipslabel.subtree = tipslabel.subtree0[startno:endno]
  
  name.subtree = sapply(strsplit(tipslabel.subtree, split = "'", fixed = T), function(x)x[2])
  id.subtree= match(name.subtree, seq.name0)
  
  if( any(NA %in% id.subtree) == "TRUE"){
    
    bug = which(is.na(id.subtree))
    
    match.id = c()
    for(i in 1: length(bug)){
      match = grep(name.subtree[bug][i], seq.name0[bug.id])
      if ( length(match) == 1 ){ match.id[length(match.id) + 1] = match} }
    
    if (  ( length(match.id) >= 1) & ( length(match.id) == length(bug) ) ){
      
      id.subtree[bug] = bug.id[match.id]
      
      
      seq.name0.subtree = seq.name0[id.subtree]
      seq0.subtree = seq0[id.subtree]
      
      seq.name2.subtree = seq.name0[-id.subtree]
      seq2.subtree = seq0[-id.subtree]
      
      write.fasta(seq0.subtree, file.out = "subtree.fasta", names = seq.name0.subtree)
      write.fasta(seq2.subtree, file.out = "remain.fasta", names = seq.name2.subtree)
      
      print("Done")
      
    }else{
      problem = c(bug, name.subtree[bug])
      print( problem ) }  } else{
        
        seq.name0.subtree = seq.name0[id.subtree]
        seq0.subtree = seq0[id.subtree]
        
        seq.name2.subtree = seq.name0[-id.subtree]
        seq2.subtree = seq0[-id.subtree]
        
        write.fasta(seq0.subtree, file.out = "subtree.fasta", names = seq.name0.subtree)
        write.fasta(seq2.subtree, file.out = "remain.fasta", names = seq.name2.subtree)
        
        print("Done")
        
      } }


# Turn tip names into number ######

toNotip <- function(file){
  
  library(seqinr) 
  
  file = read.fasta(file.choose())
  
  seq.name0 = attributes(file)$names
  seq0 = getSequence(file)
  
  duplicated.id = which(duplicated(seq.name0) == "TRUE")
  
  ls.name0 = as.list(seq.name0[duplicated.id])
  ls.name = sapply(ls.name0, function(x){ paste0(x[1], "b")   })
  seq.name0[duplicated.id] = ls.name
  
  ls.name1 = seq(1:length(seq.name0))
  
  write.fasta(seq0, file.out = "noduplicated.fasta", names = ls.name1)
  return(seq.name0)
  
  
  print("Done")
  
}


# Curate the seq ####
# seqrep = 0 remove identical seq
# seqrep = 1 remove all seq with the same nt. info

curateSeq <- function(maxamb, minseq, seqrep){
  
  library(seqinr)
  library(stringr)
  
  file = read.fasta(file.choose())
  
  seq_name0 = attributes(file)$names
       seq0 = getSequence(file)
  
  # null vector for seq to-be-delect  
  
  tobedelect1 = c()
  tobedelect2 = c()
  
  for(i in 1: length(seq0)){
    
    ATCG = c("a", "t", "c", "g")
    
    # convert character to string
    
    seq_0 = c2s(seq0[[i]])
    seq_i = gsub("-", "", seq_0)
    seq_i = gsub("~", "", seq_i)
    
    # seq length and number of ambiguous nucleotide
    
    seqlth = length( s2c(seq_i) )
       amb = length( which(! s2c(seq_i) %in% ATCG ) )
    
    if( ( seqlth < minseq ) | ( amb > maxamb ) ){
      
      tobedelect1[length(tobedelect1) + 1 ] = i
      
    }
    
  }
  
print(length(tobedelect1))
  
  # Deal with seq duplication - tobedelect2
  
  dup = duplicated(sapply(seq0, 
                          function(x){
                            
                            y = c2s(x)
                            z = gsub("-", "", y)
                            return(z)
                          }
  ))
  
  tobedelect2 <- which(dup == "TRUE")
  
print(length(tobedelect2))
  
  # find out identical seq
  
  if (seqrep == 0){
    
    # Strain name
    
    n <- "A/([A-Za-z0-9-_]+)/([A-Za-z0-9-_]+)/([0-9]+)"
    strain <- str_match(seq_name0, n)[,1]
    
    strain_dup <- which( duplicated(strain) == "TRUE")
    
    # get interset of duplicated name of seq to find the identical data
    
    tobedelect2 <- intersect(tobedelect2, strain_dup)
    
  }else{
    
    # for each replicated seqs, leave the oldest strain
    
    seq_vec <- sapply(seq0, 
                      function(x){
                        y = c2s(x)
                        z = gsub("-", "", y)
                        return(z)
                      }
    )
    # extract year
    y = "_([0-9]{4})-"
    
    tobedelect = c()
    
    for(k in 1: length(tobedelect2)){
      
      # find all duplicated seq  
      dup0 <- which( match(seq_vec, seq_vec[tobedelect2[k]]) != "NA" )
      
      # old related seq      
      dup_oldest <- dup0[ which.min(str_match(seq_name0[dup0], y)[,2]) ]
      
      # add to a new vector
      tobedelect <- unique ( c(tobedelect, dup0[ which(!dup0 %in% dup_oldest == TRUE) ]) )
      
    }
    
    tobedelect2 <- tobedelect
  } 
  
  if ( ( length(tobedelect1) + length(tobedelect2) ) > 0){
    
    tobedelect <- sort( unique( c(tobedelect1, tobedelect2) ) )
    ramain <- seq(1:length(seq0))[-tobedelect]
    
    seq_name_out = seq_name0[ramain]
    seq_out = seq0[ramain]
    
    # write fasta file
    
write.fasta(seq_out, 
            file.out = "~/Desktop/curateSeq.fasta", 
            names = seq_name_out)           
    
    # extract fastaInfo
    # Time
    
    d = "([0-9]{4})-([0-9]{2})-([0-9]{2})"
    
    time_raw <- str_match(seq_name_out, d)[,1]
    
    # Site of isolate
    
    s = "A/([A-Za-z0-9-_]+)/"
    
    site <- str_match(seq_name_out, s)[,2]       
    
    fastaInfo <- data.frame(no = seq(1:length(seq_name_out)), seq_name_out, site, time_raw)
    
print("DONE")  
return(fastaInfo)
    
    
  }else{
    
print("DONE")
  }
  
}

# To no redundant .fasta ####

tonoredundant <- function(){
  
  file = read.fasta(file.choose())
  
  seq_name0 = attributes(file)$names
  seq0 = getSequence(file)
  
  # duplicated id
  duplicated_id = which(duplicated(seq_name0) == "TRUE")
  
  # duplicated seq
  dup = duplicated(sapply(seq0, 
                          function(x){
                            
                            y = c2s(x)
                            z = gsub("-", "", y)
                            return(z)
                          }
  ))
  
  tobedelect <- which(dup == "TRUE")
  
  # take intersection
  toberemoveid <- intersect(duplicated_id, tobedelect)
  ramain <- seq(1:length(seq0))[-toberemoveid]
  
  seq_name_out = seq_name0[ramain]
  seq_out = seq0[ramain] 
  
  write.fasta(seq_out, 
              file.out = "~/Desktop/nonredundant.fasta", 
              names = seq_name_out)         
  
  print(seq_name0[toberemoveid])    
  
}


# convert YYYY-MM-DD to YY.datee ####

phylo_date <- function(x){
  
  library(stringr)
  
  d = "([0-9]{4})-([0-9]{2})-([0-9]{2})"
  
  yr <- as.numeric(str_match(x, d)[,2])
  mo <- as.numeric(str_match(x, d)[,3])
  dy <- as.numeric(str_match(x, d)[,4])
  
  yr.0 <- paste0(yr, "-01-01")
  daydifference <- as.numeric(difftime(strptime(x, format = "%Y-%m-%d"),
                                       strptime(yr.0, format = "%Y-%m-%d"), 
                                       units = "days"))/365
  
  yr.daydifference = yr + daydifference
  return(yr.daydifference)
  
}


# convert YYYY-MM-DD to distributed table ####

epi_tomonth <- function(start, end, date_v){
  
  # start and end most be YYYYMM and numeric format 
  # date_V in "YYYY-MM-DD"
  
  library(stringr)
  
  # vector must be YYYY-MM-DD format
  d = "([0-9]{4})-([0-9]{2})-([0-9]{2})"
  mon = c("01", "02", "03", "04", "05", "06",
          "07", "08", "09", "10", "11", "12")
  
  # month distribution as table
  x = as.data.frame(
    table(
      as.numeric(
        paste0(
          str_match(date_v, d)[,2], str_match(date_v , d)[,3] 
        ))))
  
  x[,1] = as.integer(as.character(x$Var1))
  
  # parse input date
  start.y <- as.numeric(str_sub(as.character(start), 1, 4))
  start.m <- as.numeric(str_sub(as.character(start), 5, 6))
  
  end.y <- as.numeric(str_sub(as.character(end), 1, 4))
  end.m <- as.numeric(str_sub(as.character(end), 5, 6))
  
  diff.y <- end.y - start.y
  
  # y should be 
  
  if( diff.y == 0 ){
    
    y = as.integer( 
      paste0(start.y, mon[seq(start.m, end.m)] 
      ) )
  }
  
  if( diff.y == 1 ){
    
    y = as.integer(c(paste0(start.y, mon[seq(start.m, 12)]), 
                     paste0(end.y, mon[seq(1, end.m)]) 
    ))
  }
  
  if(diff.y > 1){
    
    l <- paste0(start.y, mon[seq(start.m, 12)])
    o <- paste0(end.y, mon[seq(1, end.m)])
    
    m <- as.list(seq( start.y + 1, end.y -1))
    n <- as.vector(sapply(m, 
                          function(x){ 
                            paste0(x, mon)
                          }))
    
    y <- c(l, n, o)
    
  }
  
  # find out what's difference  
  z <- setdiff(as.integer(y), x[,1])
  z_df <- data.frame(z, rep(0, length(z)))
  
  colnames(z_df) = colnames(x)
  
  z_df <- rbind(x, z_df)
  z_df <- z_df[order(z_df$Var1),]
  colnames(z_df)[1] = "Month"
  
  return(z_df)
}

# tip/branch label ####

findtaxa <- function(type, 
                     tree, 
                     targetid, 
                     target){
  
  # type 1 = branch coloring 
  # type 0 = tip shape
  # default branch color = gray
  
  library(ape)
  library(ggtree)
  
  # extract tree data
  tree.d <- fortify(tree)
  tree.d[, ncol(tree.d) + 1] <- gsub("'", "", tree.d$label)
  colnames(tree.d)[ncol(tree.d)] <- "taxaid"
  
  # for tip shape
  if (type == 0){
    
    shapetaxa <- data.frame(node = c(1:length(tree.d$isTip)), shapee = NA)
    
    for (i in 1: length(targetid)){
      
      shapetaxa$shapee[ grep(targetid[i], tree.d$taxaid) ] <- target[i]
      
    }
    
return(shapetaxa)
    
    # for branch colorring     
  }else {
    
    # new column
    
    tree.d[, ncol(tree.d) + 1] <- "gray"
    colnames(tree.d)[ncol(tree.d)] <- "colorr"
    
    # for branch extension
    
    edgematrix <- as.matrix(tree.d[,c(2,1)])
    
    # color grouping 
    
    group_color <- unique(target)
    
    for (i in 1: length(group_color)){
      
      # color as group to combine key word to targetno
      
      sub_color <- which(target == group_color[i])
      targetno <- c()
      
      for (t in 1: length(sub_color)){
        
        targetno <- unique( c(targetno, grep( targetid[ sub_color[t] ], tree.d$taxaid)) )
        
      }
      
      tobecolor = c()
      pre_targetno = length(targetno)
      post_targetno = 0
      
      # while loop 
      
      while( pre_targetno != post_targetno ){
        
        pre_targetno = length(targetno)
        
        for(k in 1:length(targetno)){
          
          # all sibiling 
          sibs <- edgematrix[
            which(edgematrix[,1] == 
                    edgematrix[which(edgematrix[,2] == targetno[k]),][1]),][,2]
          
          if (length(sibs) == 1){
            
            targetno = c(targetno, edgematrix[which(edgematrix[,2] == targetno[k]),][1])
            
          }else{
            
            if (length(which(sibs %in% targetno == "FALSE")) == 0){
              
              tobecolor = c(edgematrix[which(edgematrix[,2] == targetno[k]),][1], tobecolor)
              
              targetno = c(targetno, edgematrix[which(edgematrix[,2] == targetno[k]),][1])
            }
            
          }
          targetno = unique(targetno)
          tobecolor = unique(c(targetno, tobecolor))
        }
        
        post_targetno = length(targetno)
        
      }
      
      # coloring
      
      tree.d$colorr[tobecolor] <- group_color[i]
      
    }
return(tree.d)    
    
  }
  
}


# cleanID light version ####

cleanID_lt <- function(){
  
  library(seqinr)
  library(stringr)
  
  file = read.fasta(file.choose())
  
  seq_name0 = attributes(file)$names
  seq0 = getSequence(file)
  
  
  # Deal with problematic string
  
  seq_name = gsub(" ", "_", seq_name0)
  
  seq_name = gsub("\\(", "-", seq_name)
  seq_name = gsub("\\)", "-", seq_name)
  seq_name = gsub("\\[", "-", seq_name)
  seq_name = gsub("\\]", "-", seq_name)
  
  seq_name = gsub("\\'", "", seq_name)
  seq_name = gsub(">", "", seq_name)  
  seq_name = gsub("\\.", "-", seq_name)  
  
  
  # Dissect the string 
  # for GISAID
  
  seq_name = gsub("_-Month_and_day_unknown-", "-99-99", seq_name)
  seq_name = gsub("_-Day_unknown-", "-15", seq_name)
  
  seq_name[which( endsWith(seq_name, "-") == "TRUE")] <- 
    paste0(seq_name[which( endsWith(seq_name, "-") == "TRUE")], "15")
  
  
  # Time
  
  d = "([0-9]{4})-([0-9]{2})-([0-9]{2})"
  
  time_raw <- str_match(seq_name, d)[,1]
  
  # Site of isolate
  
  s = "A/([A-Za-z0-9-_]+)/"
  
  site <- str_match(seq_name, s)[,2]
  
  
  # Deal with replicated
  
  duplicated_id = which(duplicated(seq_name) == "TRUE")
  duplicated_note = duplicated(seq_name)
  
  # loop for each replicated case
  
  if( length(duplicated_id) > 0 ){
    
    # for all the duplicated id
    
    for (i in 1: length(duplicated_id)){
      
        dup0 = which(match(seq_name, seq_name[duplicated_id[i]]) != "NA")
      
      # create a null vector to appendex  
      
         app = c("", paste0("_",letters), toupper(paste0("_",letters)))
      
       ap.id = seq_name[dup0]
      app.id = c()
      
      # loop to deal with multiple replicated
      # find the date info, insert labeling in the middle
      
      for (k in 1: length(ap.id)){
        
        app.id[k] = paste0(ap.id[k], app[k])
      }
      
      # back to seq_name0  
      
      seq_name[dup0] = app.id
      
    }
  }     
  
  duplicated_id_ed = which(duplicated(seq_name) == "TRUE")
  
  
  if (length(duplicated_id_ed) > 0 ){
    
    print("ERROR") 
    
  }else{
    
    # write fasta file
    
    write.fasta(seq0, 
                file.out = "~/Desktop/cleanID.fasta", 
                names = seq_name)
    
    # dataframe output
    
    identicalID = as.character(duplicated_note)
    
    fastaInfo <- data.frame(no = seq(1:length(seq_name)),seq_name, site, time_raw, identicalID)
    
    
    return(fastaInfo)
    print("DONE")
    
  }
  
}


# curateSeq light version ####

curateSeq_lt <- function(maxamb, minseq, seqrep){
  
  library(seqinr)
  library(stringr)
  
  file = read.fasta(file.choose())
  
  seq_name0 = attributes(file)$names
  seq0 = getSequence(file)
  
  # null vector for seq to-be-delect  
  
  tobedelect1 = c()
  tobedelect2 = c()
  
  for(i in 1: length(seq0)){
    
    ATCG = c("a", "t", "c", "g")
    
    # convert character to string
    
    seq_0 = c2s(seq0[[i]])
    seq_i = gsub("-", "", seq_0)
    seq_i = gsub("~", "", seq_i)
    
    # seq length and number of ambiguous nucleotide
    
    seqlth = length( s2c(seq_i) )
    amb = length( which(! s2c(seq_i) %in% ATCG ) )
    
    if( ( seqlth < minseq ) | ( amb > maxamb ) ){
      
      tobedelect1[length(tobedelect1) + 1 ] = i
      
    }
    
  }
  
  print(length(tobedelect1))
  
  # Deal with seq duplication - tobedelect2
  
  dup = duplicated(sapply(seq0, 
                          function(x){
                            
                            y = c2s(x)
                            z = gsub("-", "", y)
                            return(z)
                          }
  ))
  
  tobedelect2 <- which(dup == "TRUE")
  
  print(length(tobedelect2))
  
  # find out identical seq
  
  if (seqrep == 0){
    
    # Strain name
    
    n <- "A/([A-Za-z0-9-_]+)/([A-Za-z0-9-_]+)/([A-Za-z0-9-_]+)"
    strain <- str_match(seq_name0, n)[,1]
    
    strain_dup <- which( duplicated(strain) == "TRUE")
    
    # get interset of duplicated name of seq to find the identical data
    
    tobedelect2 <- intersect(tobedelect2, strain_dup)
    
  }else{
    
    # find labeled & duplicated
    
    y <-  "_([A-Z0-9]{4})_([a-zA-Z]{1})"
    strain <- str_match(seq_name0, y)[,1]
    
    labeled <- which( strain != "NA" )
    
    tobedelect = unique( c(labeled, tobedelect2) )
    
    tobedelect2 <- tobedelect
  } 
  
  if ( ( length(tobedelect1) + length(tobedelect2) ) > 0){
    
    tobedelect <- sort( unique( c(tobedelect1, tobedelect2) ) )
    ramain <- seq(1:length(seq0))[-tobedelect]
    
    seq_name_out = seq_name0[ramain]
    seq_out = seq0[ramain]
    
    # write fasta file
    
    write.fasta(seq_out, 
                file.out = "~/Desktop/curateSeq.fasta", 
                names = seq_name_out)           
    
    
    print("DONE")  
    
  }else{
    
    print("DONE")
  }
  
}