############### Strain name ############### 
#
#
# source("https://bioconductor.org/biocLite.R")

library(seqinr)
fasinput = read.fasta(file.choose())
attributes(fasinput)$names              # list the seq

  mx = getSequence(allfl)

Seq.name = c(1:3346, "GSGD", 3348:8507)    # assign the numbering

write.fasta(mx, file.out = "_name.fas", names = Seq.name)

library(phytools)

  test=read.tree(file.choose())
  test$edge  
  test$tip.label
  plotTree(test)
  nodelabels()
  
  x = c("KX014886", "KX523694")
  sTree = drop.tip(test, which(!test$tip.label %in% x))
  plotTree(sTree)
  
# library(geiger)


##########
  
  
  source("https://bioconductor.org/biocLite.R")
  biocLite()
  
  biocLite("phytools")
  
  library(seqinr)
  
  
  cb10112aligntrimul = read.fasta(file.choose())
  attributes(cb10112aligntrimul)$names                         # list the seq
  
  seq.name = attributes(cb10112aligntrimul)$names              # list the seq
  mx = getSequence(cb10112aligntrimul)
  
  table.name = cbind(seq.name, Seq.name.no)                    # 4952, 4953
  
  Seq.name.no = c(1:4951, "GSGD", 4953:length(seq.name))       # assign the numbering
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
  cb10112aligntrimul.tree$tip.label                              # 1921
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
      tobedelect[length(tobedelect) + 1] = i } }
  
  
  dup = duplicated(sapply(submx, function(x){
    
    seq0 = c2s(x)
    seq = gsub("-", replacement = "", seq0) }))
  
  tobedelect2=which(dup == "TRUE")
  
  
  curated = sort(unique(c(tobedelect, tobedelect2)))
  ramain = seq(1:length(submx))[-curated]
  
#####
  

  submx4102 = submx[ramain]
  submx4102.name = sub.seq.name[ramain]
  
  write.fasta(submx4102, file.out = "R_sub4102.fas", names = submx4102.name)
  
  



