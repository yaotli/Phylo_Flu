library(seqinr) 

# ID cleaning ####

g_20170112_8135 <- cleanID()

# remove _A_/
         file = read.fasta(file.choose())
    seq_name0 = attributes(file)$names
         seq0 = getSequence(file)
         
    seq_name = sub("_A_/", "",seq_name0)
    
    write.fasta(seq0, 
                file.out = "~/Desktop/cleanID2.fasta", 
                names = seq_name)

n_20170112_3964 <- cleanID()


# seq curation 

# n = 6840
g_20170112_8135 <- curateSeq(5, 1500, 1)

# n = 3413
n_20170112_3964 <- curateSeq(5, 1500, 1)


pool_10653 <- curateSeq(5, 1500, 0)


# clade define file (from Dr. Gavin Smith)
# n = 238 

         file = read.fasta(file.choose())
    seq_name0 = attributes(file)$names
         seq0 = getSequence(file)

    seq_name = sub("\\{new\\}", "",seq_name0)
    
    seq_name = gsub(" ", "_", seq_name)
    seq_name = gsub("\\(", "-", seq_name)
    seq_name = gsub("\\)", "-", seq_name)
    seq_name = gsub("\\[", "-", seq_name)
    seq_name = gsub("\\]", "-", seq_name)
    seq_name = gsub("\\'", "", seq_name)
    seq_name = gsub(">", "", seq_name)  
    seq_name = gsub("\\.", "-", seq_name) 
    seq_name = gsub("\\|", "_", seq_name)  
    seq_name = gsub("\\{", "_", seq_name)  
    seq_name = gsub("\\}", "", seq_name)  
    
    write.fasta(seq0, 
                file.out = "~/Desktop/smallref_cleanID.fasta", 
                names = seq_name)
    
# align and trim file: pool_7982ref = pool_10653_curateSeq + smallref_cleanID
# another round of cleanID and curate
# resulting n = 7229    
    
    pool_7982ref_trim = cleanID()
    pool_7982ref_trim = curateSeq(99, 99, 0)
    

# read-in tree ####
    
    library(ape)
    library(ggtree)
    
    sub_7225 <- read.tree(file.choose())
    ggtree(sub_7225)
    
    
    
    
    
    
    
    
    
    

