# Pre-clean the .fas file divided according to gene and host ####
  
  # source: an.fasta

library(seqinr)
library(stringr)

file <- read.fasta(file.choose())

  seq_name0 = attributes(file)$names
       seq0 = getSequence(file)
  
  # _{}          
  seq_name0[2] = gsub(pattern = "_\\{\\}", replacement = "", seq_name0[2])
  
  # not in a flu code
  fluform <- "A/([A-Za-z0-9_-]+)/([A-Za-z0-9_-]+)"
  
  tobeedit = which(is.na(str_match(seq_name0, fluform)[,1]) == TRUE)
  
  # 1-6: 
  # "|AX018718|_Equine_influenza_virus_H3N8_H3N8_--" 
  
  seq_name0[tobeedit][1:6] = 
  gsub(pattern = "Equine_influenza_virus_H3N8", 
       replacement = "A/equine/unknown", seq_name0[tobeedit][1:6]) 
  
  # 7  
  # "|AJ344023|_swine/Finistere/127/99_H3N2_1999--"
  
  seq_name0[tobeedit][7] = 
    gsub(pattern = "swine/", 
         replacement = "A/swine/", seq_name0[tobeedit][7]) 
  
  # 16 
  # "."
  
  seq_name0[tobeedit][16] = 
    gsub(pattern = "\\.", 
         replacement = "/", seq_name0[tobeedit][16]) 
  
  # 18, 19, 20, 22 
  
  seq_name0[tobeedit][c(18,19,20,22)] = 
    gsub(pattern = "_A_", 
         replacement = "_A/", seq_name0[tobeedit][c(18,19,20,22)]) 
  
  # 21
  
  seq_name0[tobeedit][21] = 
    gsub(pattern = "_#A_", 
         replacement = "_A/", seq_name0[tobeedit][21]) 
  
  # 23
  
  seq_name0[tobeedit][23] = 
    gsub(pattern = "_A_/", 
         replacement = "_A/", seq_name0[tobeedit][23]) 
  
  # Equine_Sweden_
  seq_name0 = 
    gsub(pattern = "Equine_Sweden_", replacement = "Equine/Sweden/", seq_name0)
  
  seq_name0 = 
    gsub(pattern = "_A//_", replacement = "_A_/_", seq_name0)
  
  
  write.fasta(seq0, 
              file.out = "~/Desktop/precleanID.fasta", 
              names = seq_name0)
  
# remove accession number ####
  
  # source: h3_human_6564_forsample.fasta
  
  file <- read.fasta(file.choose())
  
  seq_name0 = attributes(file)$names
       seq0 = getSequence(file)
   
  assessiono <- "\\|([A-Za-z0-9_-]+)\\|_"
  
  seq_name0 = 
  gsub(pattern = assessiono,
       replacement ="", 
       x = seq_name0)
  
  write.fasta(seq0, 
              file.out = "~/Desktop/precleanID.fasta", 
              names = seq_name0)
  
  
# add vx strain and random human seq ####
  
  # source: preC_h3_ha_an_comb.fasta
  # vx and rSeq combined: h3_pool_hu_111
  
  rSeq(n = 100, seed = 16)
  

  
# Clean ID and curate
  
  # source: addhu_preC_h3_ha_an
  # manually remove _{animal_name}
  
  cleanID()
  
  curateSeq(maxamb = 5, minseq = 1600, mode = 8)
  
  # delete: 11613
  # remain: 4298
  
# align and trim
  
  # MAFFT
  # source: curateSeq-8_h3_ha_an.fasta
  # manually remove i 
  
  curateSeq(maxamb = 5, minseq = 1000, mode = 2)
  
# FastTree
  
  # source: curateSeq-2_trim_h3_ha_an.fasta
  
  
  
  
  
  
  
   
  