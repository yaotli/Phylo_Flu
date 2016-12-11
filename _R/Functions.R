##### getDescendants (from http://blog.phytools.org; 20120126)
#

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


##### 

########## Clean tree IDs ########## 

# require package seqinr and stringr
# deal with duplicated, problematic string and extract info
# fill the date info with either -15 or -99-99


IDcleaner <- function(){
  
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
      
      app = c("a_", "b_", "c_", "d_", "e_", "f_", "g_", "h_", "i_", "j_")
      
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
                file.out = "~/Desktop/IDcleaner.fasta", 
                names = seq_name)
    
    # dataframe output
    
    identicalID = as.character(duplicated_note)
    
    fastaInfo <- data.frame(no = seq(1:length(seq_name)),seq_name, site, time_raw, identicalID)
    
    
return(fastaInfo)
print("DONE")
    
  }
  
}





########## Subtree seq extraction ########## 


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


##########
#

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
