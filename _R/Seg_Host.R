# divide the data by segment 
# by recognize _{}
#

library(seqinr)
library(stringr)

file <- read.fasta(file.choose())

  seq_name0 = attributes(file)$names
       seq0 = getSequence(file)
       
  for (k in 1: 8){
    
    endwithnumber <- paste0("_{", k, "}")
    
    assign( paste0("seg.", k), which( endsWith(seq_name0, endwithnumber) == "TRUE")  )
    
  }
  
  seg.9 <- setdiff(seq(1, length(seq0)), 
                   c(seg.1, seg.2, seg.3, seg.4, seg.5, seg.6, seg.7, seg.8) )
  
  
  str_sub(seq_name0, -4, -1) = ""
  
  for (k in 1: 9){
    
    segno <- get(paste0("seg.", k))
    
    if( length( segno ) > 0 ){
      
      seqname_out <- seq_name0[segno]
          seq_out <- seq0[segno]
          
      write.fasta(seq_out,
                  file.out = paste0("~/Desktop/Seq_", k, ".fasta"), 
                  name = seqname_out)    
      
      }
  }
  
print("DONE")  
  
# divide the data by species: human / animal
# by recognize _{}


file <- read.fasta(file.choose())

  seq_name0 = attributes(file)$names
       seq0 = getSequence(file)


    humannote <- "_\\{Human\\}"
  speciesnode <- "_\\{[A-Za-z0-9_-]+\\}"
  
  hu <- grep(humannote, seq_name0)
  
  if ( length(hu) > 0 ){
    
    an <- seq(1:length(seq0))[-hu]
    
  }else{
    
    an <- seq(1:length(seq0))
    
  }
  
  seq_name <- gsub(pattern = speciesnode, x = seq_name0, replacement = "")

  seq_an <- seq0[an]
  seq_hu <- seq0[hu]
  
  name_an <- seq_name[an]
  name_hu <- seq_name[hu]
  
  
  write.fasta(seq_an,
              file.out = "~/Desktop/an.fasta", 
              name = name_an)
  
  write.fasta(seq_hu,
              file.out = "~/Desktop/hu.fasta", 
              name = name_hu)
  
print("DONE")
  
  
  
  
  
  
  
  
  
  
  
