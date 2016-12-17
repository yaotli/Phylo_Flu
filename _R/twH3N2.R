library(seqinr) 

# ID cleaning ####

   gisaidh3 <- cleanID()
pooledh3_tw <- cleanID()


# seq curation ####

# maxamb = 5, minseq = 1500, replicates are allowed in TW data
# 1647 / 15 reps

    gisaidh3_cr <- curateSeq(5, 1500, 1)
 pooledh3_tw_cr <- curateSeq(5, 1500, 0)
 
 
# align and trim
# remove apparent duplicated data (Taiwan)
 
library(seqinr)
library(stringr)
 
 # read in 
 file = read.fasta(file.choose())
 
 seq_name0 = attributes(file)$names
      seq0 = getSequence(file)
 
   toberemove <- pooledh3_tw_cr[-(8:11)]
 toberemoveid <-match(toberemove, seq_name0)
 
 ramain <- seq(1:length(seq0))[-toberemoveid]
 
 seq_name_out = seq_name0[ramain]
      seq_out = seq0[ramain] 
 
 write.fasta(seq_out, 
             file.out = "~/Desktop/out_seq.fasta", 
             names = seq_name_out)         
 
 
 # deal with replated between TW from gisaid pool
 # since ERROR from FastTree, remove additional 6 seqs
 # resulting file: allpooled_trim_mo2
 
 tonoredundant()
 
 
# 
 
 
 
 
 
 
 
 
 
 
 