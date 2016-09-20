source("https://bioconductor.org/biocLite.R")
biocLite()

biocLite("phytools")

library(seqinr)


cb10112aligntrimul = read.fasta(file.choose())
attributes(cb10112aligntrimul)$names              # list the seq

seq.name = attributes(cb10112aligntrimul)$names              # list the seq
mx = getSequence(cb10112aligntrimul)

table.name = cbind(seq.name, Seq.name.no)

# 4952, 4953

Seq.name.no = c(1:4951, "GSGD", 4953:length(seq.name))    # assign the numbering
write.fasta(mx, file.out = "cb10112aligntrimul.fas", names = Seq.name.no)


#####
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

library(phytools)
cb10112aligntrimul.tree=read.tree(file.choose())
cb10112aligntrimul.tree$tip.label      # 1921
# 1264 for tip3668
edge.cb10112aligntrimul.tree = cb10112aligntrimul.tree$edge  

length(getDescendants(cb10112aligntrimul.tree, 11887))  # 11887 with 8195 descendant tips
length(getDescendants(cb10112aligntrimul.tree, 11276)) 

subtree.no = getDescendants(cb10112aligntrimul.tree, 11276)
subtree.label = cb10112aligntrimul.tree$tip.label[subtree.no]

subtree.ID = match(subtree.label, Seq.name.no)
submx = mx[subtree.ID]
sub.seq.name = seq.name[subtree.ID]


write.fasta(submx, file.out = "R8849sim.fas", names = subtree.label)
write.fasta(submx, file.out = "R8849.fas", names = sub.seq.name)



#######




tobedelect=c()
for(i in 1: length(submx)){
  ATCG = c("a","t","c","g")
  
  seq0 = c2s(submx[[i]])
  seq = gsub("-", replacement = "", seq0)
  
  seq1 = s2c(seq)
  amb = length(which( seq1 %in% ATCG == "FLASE"))
  
  lf = length(seq1)
  
  if ((lf< 1600)|(amb >5)){ 
    tobedelect[length(tobedelect) + 1] = i
    
  }
  
}


dup = duplicated(sapply(submx, function(x){
  
  seq0 = c2s(x)
  seq = gsub("-", replacement = "", seq0)
}))

tobedelect2=which(dup == "TRUE")


curated = sort(unique(c(tobedelect, tobedelect2)))
ramain = seq(1:length(submx))[-curated]

#####



submx4102 = submx[ramain]
submx4102.name = sub.seq.name[ramain]


write.fasta(submx4102, file.out = "R_sub4102.fas", names = submx4102.name)



###############

########## Arrange duplicated name ########## 



rmdup<-function(file){
  
  library(seqinr) 
  
  file = read.fasta(file.choose())
  
  seq.name0 = attributes(file)$names
  seq0 = getSequence(file)
  
  duplicated.id = which(duplicated(seq.name0) == "TRUE")
  
  ls.name0 = as.list(seq.name0[duplicated.id])
  ls.name = sapply(ls.name0, function(x){ paste0(x[1], "b")   })
  
  seq.name0[duplicated.id] = ls.name
  
  write.fasta(seq0, file.out = "noduplicated.fasta", names = seq.name0)
  
  print("Done")
  
}



# rmdup()


########## Subtree seq extraction ########## 
#
# Requirement: 
#     1. Original fasta for the tree
#     2. .csv file with well-aligned SUBtree tip labels
#     3. Packages: bioconductor::seqinr
#     4. Be sure that all tip names are in ''

subtreseq()

subtreseq<-function(){
  
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
      
      write.fasta(seq0.subtree, file.out = "subtree.fasta", names = seq.name0.subtree)
      
      print("Done")
      
      
    }else{
      problem=c(bug, name.subtree[bug])
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
      problem=c(bug, name.subtree[bug])
      print( problem ) }  } else{
        
        seq.name0.subtree = seq.name0[id.subtree]
        seq0.subtree = seq0[id.subtree]
        
        seq.name2.subtree = seq.name0[-id.subtree]
        seq2.subtree = seq0[-id.subtree]
        
        write.fasta(seq0.subtree, file.out = "subtree.fasta", names = seq.name0.subtree)
        write.fasta(seq2.subtree, file.out = "remain.fasta", names = seq.name2.subtree)
        
        print("Done")
        
      } }


subtreseq2()

