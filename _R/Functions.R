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
      
       app = c("","a_", "b_", "c_", "d_", "e_", "f_", "g_", "h_", "i_", "j_")
      
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
  bug.id = grep("'", seq.name0)
  
  sub.tree=read.csv(file.choose())
  
  tipslabel.subtree0 = as.character(sub.tree[,1])
  startno = which( tipslabel.subtree0 == "\ttaxlabels" ) + 1
  endno = head( which(tipslabel.subtree0 == ";"), 1 ) - 1
  
  tipslabel.subtree = tipslabel.subtree0[ startno:endno ]
  
  name.subtree = sapply(strsplit(tipslabel.subtree, split = "'", fixed = T), function(x)x[2])
  id.subtree = match(name.subtree, seq.name0)
  
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
      
      write.fasta(seq0.subtree, file.out = "subtree.fasta", names = seq.name0.subtree)
      
      print("Done")
      
      
    }else{
      problem = c(bug, name.subtree[bug])
      print( problem ) }  } else{
        
        seq.name0.subtree = seq.name0[id.subtree]
        seq0.subtree = seq0[id.subtree]
        
        write.fasta(seq0.subtree, file.out = "subtree.fasta", names = seq.name0.subtree)
        
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

toNotip<-function(file){
  
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



