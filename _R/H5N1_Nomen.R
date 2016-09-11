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



