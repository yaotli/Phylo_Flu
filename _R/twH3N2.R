library(seqinr) 

# Read-in files

gisaidGB <- read.fasta("/Users/yaosmacbook/Desktop/twH3N2/gisaid_gb_4898.fasta")
gisaidTW <- read.fasta("/Users/yaosmacbook/Desktop/twH3N2/gisaid_tw_131.fasta")
  ncbiTW <- read.fasta("/Users/yaosmacbook/Desktop/twH3N2/ncbi_tw_31.fasta")


# ID cleaning for GISAID
  
file = gisaidGB

seq_name0 = attributes(file)$names
     seq0 = getSequence(file)

     
# Deal with replicated ID
     
duplicated_id = which(duplicated(seq_name0) == "TRUE")

if( length(duplicated_id) > 0 ){
  
  # for all the duplicated id
  
  for (i in 1: length(duplicated_id)){
    
    dup0 = which(match(seq_name0, seq_name0[duplicated_id[i]]) != "NA")
    
    # create a null vector to appendex  
    
     app = c("", "_b", "_c", "_d", "_e", "_f", "_g", "_h", "_i", "_j")
     
     ap.id = seq_name0[dup0]
    app.id = c()
    
    # loop to deal with multiple replicated
    
    for (k in 1: length(ap.id)){
      
      app.id[k] = paste0(ap.id[k], app[k])
      
    }
    
    # back to seq_name0  
    
    seq_name0[dup0] = app.id
    
  }
}     
     
     

# Deal with 
     



     