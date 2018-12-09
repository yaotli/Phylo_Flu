# download criteria `query_20170926`
# replace blank with underscore for .fasta from NCBI 
# convert .xml to .csv for gisaid

require(seqinr)
require(stringr)
require(ape)
library(ggtree)

setwd("~/Phylo_Flu/H9N2/")
source("functions.R")

#gisaid
fas_ha_g <- "H9_G_168_20170926.fasta"
csv_g    <- "H9_G_168_20170926.csv"

#ncbi
fas_ha_n <- "H9_N_6607_20170926.fasta"


## examine ----------------

infolist.ha.g <- idInfo( fastaEx( fas_ha_g )$id, datasource = "g", g.csv = csv_g)

infolist.ha.n <- idInfo( fastaEx( fas_ha_n )$id, datasource = "n")

# [1] "ERROR in  year"

# [1] "AY738451_A/turkey/Neve_Ilan/90710/2000_H9N2_|Israel|_2000-5-30"                                                   
# [2] "AY738452_A/turkey/Givat_Haim/965/2002_H9N2_|Israel|_2002-3-17"                                                    
# [3] "AY738455_A/chicken/Talmei_Elazar/1304/03_H9N2_|Israel|_2003-5-05"                                                 
# [4] "KY785907_A/reassortant/H9N2:pH1N1_RGpassage7(quail/Hong_Kong/G1/1997_x_California/04/2009)_H9N2_|Hong_Kong|_NON--"

tem = fastaEx( fas_ha_n )$id
tem[ grep("AY738451", tem) ] = "AY738451_A/turkey/Neve_Ilan/90710/2000_H9N2_|Israel|_2000-05-30"
tem[ grep("AY738452", tem) ] = "AY738452_A/turkey/Givat_Haim/965/2002_H9N2_|Israel|_2002-03-17"
tem[ grep("AY738455", tem) ] = "AY738455_A/chicken/Talmei_Elazar/1304/03_H9N2_|Israel|_2003-05-05"
tem[ grep("KY785907", tem) ] = "KY785907_A/reassortant/H9N2:pH1N1_RGpassage7(quail/Hong_Kong/G1/1997_x_California/04/2009)_H9N2_|Hong_Kong|_--"

infolist.ha.n <- idInfo( tem, datasource = "n")


## combine ----------------

# n = 6775
seq_h9      <- c( fastaEx( fas_ha_n )$seq, fastaEx( fas_ha_g )$seq )

infolist_h9 <- list( c(infolist.ha.n[[1]], infolist.ha.g[[1]]),
                     c(infolist.ha.n[[2]], infolist.ha.g[[2]]),
                     c(infolist.ha.n[[3]], infolist.ha.g[[3]]), 
                     c(infolist.ha.n[[4]], infolist.ha.g[[4]]), 
                     c(infolist.ha.n[[5]], infolist.ha.g[[5]]),
                       seq_h9 )   


## remove dupicate strain and seq curation ----------------

# n = 6994
s_infolist_h9 <- strainSelect( infolist_h9 )


# format time 
s_infolist_h9[[4]] <- seqDate( s_infolist_h9[[4]] )


# save .fasta with clean id 

write.fasta( sequences = s_infolist_h9[[6]], 
             names     = paste0(s_infolist_h9[[1]], "_",
                                s_infolist_h9[[5]], "_|",
                                s_infolist_h9[[3]], "|_",
                                s_infolist_h9[[2]], "_",
                                s_infolist_h9[[4]]
                                ),
             file.out = "H9_6694.fasta"
             )

# NOTE 
grep("RG|reassort", s_infolist_h9[[5]], value = TRUE)
# [1] "reassortant_H9N2_pH1N1_RGpassage7_quail_Hong_Kong_G1_1997_x_California_04_2009_"
# [2] "chicken_Jiangsu_WJHRG_2012" 


# seq curation 

# n = 5034
c_infolist_h9 <- seqSelect( minlth = 1500, maxamb = 1, s_infolist_h9, rmdup = TRUE)

write.fasta( sequences = c_infolist_h9[[6]], 
             names     = paste0(c_infolist_h9[[1]], "_",
                                c_infolist_h9[[5]], "_|",
                                c_infolist_h9[[3]], "|_",
                                c_infolist_h9[[2]], "_",
                                c_infolist_h9[[4]] ),
             file.out  = "H9_5034.fasta")

## alignment, trim ----------------

system("mafft --reorder H9_5034.fasta > align_H9_5034.fasta")

# trim
# & manually check in BioEdit
trimtool( propblank = 0.9, filedir = "align_H9_5034.fasta" )

# further remove duplicated seq
# n = 4982
rmDup( fasfile = "trim_H9_5034_lth1680.fasta", year = c(1970, 2018), rmdup = TRUE)


system("~/./Fasttree -nt -nni 10 -spr 4 -gtr -cat 20 -gamma -notop H9_4982.fasta > H9_4982.tre")


## tree ----------------

# tree file - H9_4982.tre
# save as nwk format - H9_4982.nwk
# save as nex format - H9_4982_clade.tre


# inspect current samples 

list_H9_4982 <- taxaInfo( file = "H9_4982_clade.tre", useTree = TRUE)
# data.frame(table(list_H9_4982[[2]]))


CN <- c("China", "Hong_Kong")
sA <- c("Bangladesh", "Australia", "India", "Indonesia", "Malaysia", 
        "Myanmar", "Nepal", "New_Zealand", "Thailand", "Vietnam", 
        "Singapore")

nA <- c("Japan", "Russia", "South_Korea")
ME <- c("Afghanistan", "Egypt", "Iran", "Iraq", "Israel", 
        "Jordan", "Kuwait", "Lebanon", "Middle_East", "Pakistan", 
        "United_Arab_Emirates")

E  <- c("Austria", "Belgium", "Finland", "Georgia", "Germany", 
        "Ireland", "Italy", "Netherlands", "Norway", "Poland", 
        "Portugal","Saudi_Arabia", "Sweden", "Switzerland", 
        "United_Kingdom", "France", "Hungary")  
        
Na <- c("Argentina", "Canada", "Chile", "Mexico", "USA")
A  <- c("Libya", "Morocco", "Tunisia", "Zambia", "South_Africa")


# ggtree

rawtre <- read.tree("H9_4982.nwk")
gt     <- ggtree(rawtre, ladderize = FALSE)

geo.key = c( paste0( sA, collapse = "|"),
             paste0( CN, collapse = "|"),
             paste0( nA, collapse = "|"),
             paste0( ME, collapse = "|"),
             paste0( E, collapse = "|"),
             paste0( A, collapse = "|"),
             paste0( Na, collapse = "|"))

# http://colorbrewer2.org/?type=diverging&scheme=RdYlGn&n=7
col.key = c("#1a9850", "#91cf60", "#d9ef8b", "#fee08b", "#fc8d59", "#ffffbf", "#d73027")

geoann <- findtaxa(type     = 1, 
                   tree     = rawtre, 
                   targetid = geo.key,
                   target   = col.key)


gt %<+% geoann + aes(color = I(colorr))





# cladesampling 

cladeSampling( trefile = "./H9_4982_clade.tre", fasfile = "./H9_4982.fasta", seed = 999, 
               grid = 1.1, minBranchlth = TRUE, saveFasta = TRUE)

s.tre <- read.nexus("./H9_1612_s")
gt    <- ggtree(s.tre, ladderize = FALSE)



geo.key = c( paste0( sA, collapse = "|"),
             paste0( CN, collapse = "|"),
             paste0( nA, collapse = "|"),
             paste0( ME, collapse = "|"),
             paste0( E, collapse = "|"),
             paste0( A, collapse = "|"),
             paste0( Na, collapse = "|"))

col.key = c("#1a9850", "#91cf60", "#d9ef8b", "#fee08b", "#fc8d59", "#ffffbf", "#d73027")

geoann <- findtaxa(type     = 1, 
                   tree     = rawtre, 
                   targetid = geo.key,
                   target   = col.key)


gt %<+% geoann + aes(color = I(colorr))



