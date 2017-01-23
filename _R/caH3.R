library(seqinr)
library(stringr)

# HA ####
# be care of ID from NCBI should remove space

# remove _A_/
# cleanID() and curateSeq() [n = 2944]

# align, trim

# remove: A/swine/Illinois/A00857131/2011_H3N2_20110924
#    1 - 1790 mammalian
# 1791 - 2939 avian
# 2940 - 2943 human

# check ID
     file = read.fasta(file.choose())
seq_name0 = attributes(file)$names
     seq0 = getSequence(file)

     id <- "A/([A-Za-z0-9-_]+)/([A-Za-z0-9-_]+)/([A-Za-z0-9-_]+)"
      a <- data.frame(str_match(seq_name0, id)[,1])
     seq_name0[which(is.na(a))] = c(
       
       "A/TBY/91_H3N8_1991-99-99",
       "A/Equine/Sweden/SVA111208SZ0077VIR167905/2011_H3N8_2011-12-07", 
       "A/Equine/Sweden/SVA111208SZ0077VIR167906/2011_H3N8_2011-12-07",
       "A/Equine/Sweden/SVA111212SZ0058VIR169816/2011_H3N8_2011-12-08", 
       "A/Equine/Sweden/SVA111209SZ0099VIR169118/2011_H3N8_2011-12-08",
       "A/Equine/Sweden/SVA111128SZ0073VIR160172/2011_H3N8_2011-11-27"
     )
     
     write.fasta(seq0, 
                 file.out = "~/Desktop/checkID.fasta", 
                 names = seq_name0)
     
# curate once 
# resulting n = 2844
# M:    1 - 1752
# A: 1753 - 2840
# H: 2841 - 2844
     
     curateSeq(5, 100, 1)
     
     file = read.fasta(file.choose())
     presavedname = attributes(file)$names
     
     
# FastTree
# Tree annotation ####

library(ggtree)
library(ape)
   
# sub_1273      
     
     sub_1273 <- read.tree(file.choose())
     ggtree(sub_1273)
     
          treedata_sub1273 <- fortify(sub_1273)
     treedata_sub1273[,10] <- gsub("'", "",treedata_sub1273$label)
     
     # color for canine + feline/ equine / swine [11]
     
      n_dog <- grep("anine", treedata_sub1273[,10])
      n_cat <- grep("feline", treedata_sub1273[,10]) 
      
      n_dog = unique(c(n_dog, n_cat))
     
    n_horse <- grep("equine", treedata_sub1273[,10])
    
    n_swine <- grep("swine", treedata_sub1273[,10])
    n_Swine <- grep("Swine", treedata_sub1273[,10])
    
    n_swine = unique(c(n_swine, n_Swine))
    
    treedata_sub1273[11] <- "gray"
    colnames(treedata_sub1273)[11] = "colorr"
    
    test = tiplabel(type = 1, tdata = treedata_sub1273, datacolumn = 11, n_dog, target = "red")
    test = tiplabel(type = 1, tdata = test, datacolumn = 11, n_horse, target = "green")
    test = tiplabel(type = 1, tdata = test, datacolumn = 11, n_swine, target = "orange")
    
    an1 = ggtree(sub_1273) %<+% test + aes(color = I(colorr), alpha = 0.5)
    
    
    # color for time [12]
     
    n_2016 <- grep("2016", treedata_sub1273[,10])
    n_2015 <- grep("2015", treedata_sub1273[,10])
    
    n_t1 <- unique(c(n_2016, n_2015))
    
    n_2014 <- grep("2014", treedata_sub1273[,10])
    n_2013 <- grep("2013", treedata_sub1273[,10])
    n_2012 <- grep("2012", treedata_sub1273[,10])
    
    n_t2 <- unique(c(n_2014, n_2013, n_2012))
    
    test[12] <- "gray"
    colnames(test)[12] = "colorrT"
    
    test = tiplabel(type = 1, tdata = test, datacolumn = 12, n_t1, target = "red")
    test = tiplabel(type = 1, tdata = test, datacolumn = 12, n_t2, target = "orange")
    
    ggtree(sub_1273) %<+% test + aes(color = I(colorrT), alpha = 0.5)
    
# sub_105
    
    sub_105 <- read.tree(file.choose())
    ggtree(sub_105)
     
    treedata_sub105 <- fortify(sub_105)
    treedata_sub105[,10] <- gsub("'", "",treedata_sub105$label)
    
    # color for canine + feline
     
    n_dog <- grep("canine", treedata_sub105[,10])
    n_cat <- grep("feline", treedata_sub105[,10]) 
    
    n_dog = unique(c(n_dog, n_cat))
    
    treedata_sub105[11] <- "gray"
    colnames(treedata_sub105)[11] = "colorr"
    
    test2 = tiplabel(type = 1, tdata = treedata_sub105, datacolumn = 11, n_dog, target = "red")
    
    an2 = ggtree(sub_105) %<+% test2 + aes(color = I(colorr), alpha = 0.5)
    
    # color time
    
    n_2016 <- grep("2016", treedata_sub105[,10])
    n_2015 <- grep("2015", treedata_sub105[,10])
    
    n_t1 <- unique(c(n_2016, n_2015))
    
    test2[12] <- "gray"
    colnames(test2)[12] = "colorrT"
    
    test2 = tiplabel(type = 1, tdata = test2, datacolumn = 12, n_t1, target = "red")
    
    an3 = ggtree(sub_105) %<+% test2 + aes(color = I(colorrT), alpha = 0.5)
    
multiplot(an1, an2, an3, ncol = 3)
     
     
# NA ####

N2_pool_7393 <- cleanID()

# remove _A_/

       file = read.fasta(file.choose())
  seq_name0 = attributes(file)$names
       seq0 = getSequence(file)
   seq_name = sub("_A_/", "",seq_name0)

   write.fasta(seq0, 
               file.out = "~/Desktop/cleanID2.fasta", 
               names = seq_name)

# curateseq()   
   
   N2_pool_7393 <- curateSeq(5, 1300, 1)
   
   
# read-in tree
   
   N2_sub3047 <- read.tree(file.choose())
   ggtree(N2_sub3047)
   
         treedata_N2_sub3047 <- fortify(N2_sub3047)
    treedata_N2_sub3047[,10] <- gsub("'", "",treedata_N2_sub3047$label)
   
  # colorring
  
  n_dog <- grep("canine", treedata_N2_sub3047[,10])
  n_cat <- grep("feline", treedata_N2_sub3047[,10]) 
  
  n_dog = unique(c(n_dog, n_cat))
  
  n_swine <- grep("swine", treedata_N2_sub3047[,10])
  n_Swine <- grep("Swine", treedata_N2_sub3047[,10])
  
  n_swine = unique(c(n_swine, n_Swine))
  
  treedata_N2_sub3047[11] <- "gray"
  colnames(treedata_N2_sub3047)[11] = "colorr"
  
  test3 = tiplabel(type = 1, tdata = treedata_N2_sub3047, datacolumn = 11, n_dog, target = "red")
  test3 = tiplabel(type = 1, tdata = test3, datacolumn = 11, n_swine, target = "orange")

  an3 = ggtree(N2_sub3047) %<+% test3 + aes(color = I(colorr), alpha = 0.5)
  
  anm = c(n_dog, n_swine)
  test3_canie <- data.frame(node = c(1:length(N2_sub3047$tip.label) ), shapee = NA) 
  test3_canie <- tiplabel(type = 0, tdata = test3_canie, datacolumn = 2, anm, target = 16)

  an4 = an3 %<+% test3_canie +  geom_tippoint(aes(shape = factor(shapee)), size = 2)
  

  # color time 
  
  n_2016 <- grep("2016", treedata_N2_sub3047[,10])
  n_2015 <- grep("2015", treedata_N2_sub3047[,10])
  
  n_t1 <- unique(c(n_2016, n_2015))

  test3[12] <- "gray"
  colnames(test3)[12] = "colorrT"
  
  test3 = tiplabel(type = 1, tdata = test3, datacolumn = 12, n_t1, target = "red")
  
  an5 = ggtree(N2_sub3047) %<+% test3 + aes(color = I(colorrT), alpha = 0.5)
  an6 = an5 %<+% test3_canie +  geom_tippoint(aes(shape = factor(shapee)), size = 2)
  
  

     
     
     
     
     