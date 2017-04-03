library(seqinr)
library(stringr)

  # H3_avian (n = 1385)
  # H3_mammal(n = 2147)
  # H3_human (n = 11 )

# pool 3 .fasta (n = 3643)

# clean ID

cleanID()

# curateSeq

  # mode = 8
  # remain: 2880

curateSeq(maxamb = 5, minseq = 1600, mode = 8)

# pool: curateSeq-8_H3_3643.fasta and H3_sample

  # sample (n = 35)

# manually remove one "I" from poolsample_H3_2915.fasta
# MAFFT
# trim to the length of 1698 nt

# curateSeq2

  # remain: 2826 

curateSeq(maxamb = 10000, minseq = 0, mode = 2, vip = 35)

# FastTree
# edit tree

# random human isolates ####

  # H3_human (n = 6564)  

rSeq(n = 50, seed = 99)

# pool 3 .fasta (n = 3682)

# clean ID

cleanID()

# curateSeq

# mode = 8
# remain: 2918

curateSeq(maxamb = 5, minseq = 1600, mode = 8)

# pool: curateSeq-8_h3_3682.fasta and H3_sample

# sample (n = 35)

# manually remove one "I" from poolsample_H3_2915.fasta
# MAFFT
# trim to the length of 1698 nt

# curateSeq2

# remain: 2862 

curateSeq(maxamb = 10000, minseq = 0, mode = 2, vip = 35)

# FastTree
# edit tree










