library(seqinr) 

# ID cleaning ####

  g_gb_pooled_22685 <- cleanID()
       g_gb_fg_4863 <- cleanID()
      tw_pooled_168 <- cleanID()


# seq curation ####

# maxamb = 5, minseq = 1500, replicates are allowed in TW data

    # 39; 1647; 3180
        g_gb_fg_cleanID <- curateSeq(5, 1500, 1)
    
    # 58; 7479; 15151
    g_gb_pooled_cleanID <- curateSeq(5, 1500, 1)
    
    # 0; 15; 157
      tw_pooled_cleanID <- curateSeq(5, 1500, 0)
    
    # 0; 129; 15186
            pooled_gbtw <- curateSeq(5, 1500, 0)
    
    # 0; 14; 3331      
         pooled_gbtw_fg <- curateSeq(5, 1500, 0)
      
      
# align and trim
 
# Distribution of TW data ####
 
library(ggplot2)
library(ggtree)
library(ape)
library(stringr)
 
 # read-in
 
 s3273_pooled_gbtw_fg_align <- read.tree(file.choose())
 ggtree(s3273_pooled_gbtw_fg_align)
 
 # tablized tree data
 
 treedata_s3273_pooled_gbtw_fg <- fortify(s3273_pooled_gbtw_fg_align)   
 
 # _ID
 
 treedata_s3273_id <- gsub("'", "",
                           treedata_s3273_pooled_gbtw_fg$label[which(treedata_s3273_pooled_gbtw_fg$label != "NA")])
 # _time
 
 d = "([0-9]{4})-([0-9]{2})-([0-9]{2})"
 treedata_s3273_date <- str_match(treedata_s3273_id, d)[,1]
 
 # _isolation place
 
 s = "A/([A-Za-z0-9-_]+)/"
 treedata_s3273_place <- str_match(treedata_s3273_id, s)[,2]
 
 # label TW isolates
 # number = 150
 
 treedata_s3273_TW <- rep(0, length(treedata_s3273_place))
 treedata_s3273_TW[grep("Tai",treedata_s3273_place)] <- 1
 treedata_s3273_id[which(treedata_s3273_TW == 1)]
 
 du = "/([0-9]{4})_([a-z])_([0-9{4}])"
 
 # eliminate seq with same info (ex: _a_)
 # resulting TW isolates: 128
 
 treedata_s3273_TWb = rep(0, length(treedata_s3273_place))
 treedata_s3273_TWb[ which(treedata_s3273_TW == 1)[-grep(du, treedata_s3273_id[which(treedata_s3273_TW == 1)])] ] <- 1
 
   s3273_TW_id <- treedata_s3273_id[which(treedata_s3273_TWb == 1)]
 s3273_TW_date <- str_match(s3273_TW_id, d)[,1]
 
 # s3273_TW_id[117] without month and date info
 # which(is.na(phylo_date(s3273_TW_date)) == TRUE)
 
 # continuous varibale of date
 s3273_TW_time <- phylo_date(s3273_TW_date)
 
 # deal with TW CDC chart
 # http://nidss.cdc.gov.tw/en/SingleDisease.aspx?dc=1&dt=4&disease=487a
 twcdcchart <- read.csv(file.choose())[,c(1,4)]
 colnames(twcdcchart) <- c("Month", "Case")
 
 # epi-curve
   
   llb <- rep(2009: 2016, each = 2 )
   llb[seq(1,15,2)] = paste0(llb[seq(1,15,2)], "/1")
   llb[seq(2,16,2)] = paste0(llb[seq(2,16,2)], "/7")
   llb <- gsub("20", "", llb)
   
   p<-ggplot(twcdcchart, aes(x = c(1:96), y = Case)) + 
   geom_point(size = 1.5) + 
   geom_line(size = 1 ) +
   theme_bw() + 
   scale_x_continuous(breaks = c(seq(1, 96, by = 6)), 
                      labels = llb) + 
   xlab("") +  
     
     theme(
     panel.border = element_rect(colour = "black", fill=NA, size=1),
     axis.title.y = element_text(size = 20),
     axis.text.y = element_text(size = 15))
   
       
     
     

 # make sample distribution 
 # only 115 seqs during 20090101 - 20161230 are available
   
   s3273_TW_time_dist <- epi_tomonth(200901, 201612, 
                                     s3273_TW_date[-117][which(s3273_TW_date[-117] > "2009")])
   
   q<-ggplot(s3273_TW_time_dist, aes(x = c(1:96), y = Freq)) + 
     geom_line(size = 1, color = "Blue") +
     theme_bw() + 
     scale_x_continuous(breaks = c(seq(1, 96, by = 6)), 
                        labels = c()) + 
     scale_y_continuous(breaks = c(seq(0,10, by = 2))) +
     xlab("") + 
     ylab("Seq") + 
     
     theme(
       panel.border = element_rect(colour = "black", fill=NA, size=1),
       axis.title.y = element_text(size = 20),
       axis.text.y = element_text(size = 15))
   
  # combine 2 fig
  multiplot(q, p, ncol =1) 
   
   
# locate TW strains on the tree ####
  
# treedata_s3273_TW  
# treedata_s3273_TWb  
  
  library(ggtree)
  library(stringr)
  
  # ggtree(s3273_pooled_gbtw_fg_align)
  # treedata_s3273_pooled_gbtw_fg = 6454 node; 3273 tip

  treedata_s3273_pooled_gbtw_fg[,10] <- gsub("'", "",
                                             treedata_s3273_pooled_gbtw_fg$label)
  # _time
  
  d = "([0-9]{4})-([0-9]{2})-([0-9]{2})"
  
  treedata_s3273_pooled_gbtw_fg[,11] <- str_match(treedata_s3273_pooled_gbtw_fg$label, d)[,2]
  treedata_s3273_pooled_gbtw_fg[,12] <- str_match(treedata_s3273_pooled_gbtw_fg$label, d)[,3]
  treedata_s3273_pooled_gbtw_fg[,13] <- str_match(treedata_s3273_pooled_gbtw_fg$label, d)[,4]
  
  # _isolation place
  
  s = "A/([A-Za-z0-9-_]+)/"
  treedata_s3273_pooled_gbtw_fg[,14] <- str_match(treedata_s3273_pooled_gbtw_fg$label, s)[,2]
  
  # deal with TW
  # 150 ID contain "Tai"
  
  treedata_s3273_pooled_gbtw_fg[,15] <- rep(0, 
                                            dim(treedata_s3273_pooled_gbtw_fg)[1])
  
  treedata_s3273_pooled_gbtw_fg[,15][ grep("Tai", treedata_s3273_pooled_gbtw_fg[,14])] <- 1
  
  
  # eliminate seq with same info
  # remaining = 128
  
  du = "/([0-9]{4})_([a-z])_([0-9{4}])"
  
  treedata_s3273_pooled_gbtw_fg[,16] <- rep(0, 
                                            dim(treedata_s3273_pooled_gbtw_fg)[1])
  
  
  treedata_s3273_pooled_gbtw_fg[, 16][ which(
    treedata_s3273_pooled_gbtw_fg[, 15] == 1)[-grep(du, 
                                                  treedata_s3273_pooled_gbtw_fg[, 10][which(
                                                    treedata_s3273_pooled_gbtw_fg[,15] == 1) ])] 
    ] <- 1
  
  treedata_s3273_pooled_gbtw_fg[,17] <- rep("grey", dim(treedata_s3273_pooled_gbtw_fg)[1])
  treedata_s3273_pooled_gbtw_fg[,17][which(treedata_s3273_pooled_gbtw_fg[,16] == 1)] = "red"
  
  colnames(treedata_s3273_pooled_gbtw_fg)[10:17] <- c("ID", "Year", "Month", "Day", "Place", "TW", "TWb", "TWcolor")
  
  
  # TW distribution 
  # note: %<+%, color = I, alpha, geom_tippoint
  
  shapee = c(rep(NA, 3273))
  shapee[s3273_TW_node[-117]] = 1
  treetipshape = data.frame(node=c(1:3273), shapee)
  
  c3 <- ggtree(s3273_pooled_gbtw_fg_align) %<+% treedata_s3273_pooled_gbtw_fg + aes(color = I(TWcolor), alpha = 0.5)
  t4 = c3 %<+% treetipshape + geom_tippoint(aes(shape = factor(shapee)))
  
# Temporal distribution of global H3N2 ####
  
  # create a continus temporal parameter: treedata_s3273_pooled_gbtw_fg$YYYYMM
  
  treedata_s3273_pooled_gbtw_fg[, 18] <- as.numeric(paste0(treedata_s3273_pooled_gbtw_fg$Year,
                                             treedata_s3273_pooled_gbtw_fg$Month))
  
  colnames(treedata_s3273_pooled_gbtw_fg)[18] <- "YYYYMM"
  
  # arbitrarily set epi season as (Nov - April)[0], inter-epi (May - Oct)[1] 
  episeason = c(seq(200610.1, 201610.1, by=100), seq(200704.1, 201604.1, by=100))
  episeason = sort(episeason)
  
  treedata_s3273_pooled_gbtw_fg[, 19] <- 0
  
  # designate each seq to each interval
  library(dplyr)
  
  for (i in 1: 21){
    
    if (i == 1){
      
      x = treedata_s3273_pooled_gbtw_fg %>%
        filter( YYYYMM < episeason[i] ) %>%
        select(node)
      
      x <- unlist(x, use.names = FALSE)
      treedata_s3273_pooled_gbtw_fg[, 19][x] = i
      
    }else{
      
      x = treedata_s3273_pooled_gbtw_fg %>%
        filter( episeason[i-1] < YYYYMM & YYYYMM < episeason[i]) %>%
        select(node)
      
      x <- unlist(x, use.names = FALSE)
      treedata_s3273_pooled_gbtw_fg[, 19][x] = i
      
    }
  }
    
  
  treedata_s3273_pooled_gbtw_fg[, 19][which(treedata_s3273_pooled_gbtw_fg$Month == 99) ] = 0
  
  # now we have "period"
  colnames(treedata_s3273_pooled_gbtw_fg)[19] <- "Period"
  table(treedata_s3273_pooled_gbtw_fg$Period[1:3273])

  # annual epidemic period on tree
  
  treedata_s3273_pooled_gbtw_fg[, 20] <- "#CCCCCC"
  treedata_s3273_pooled_gbtw_fg[, 20][which( treedata_s3273_pooled_gbtw_fg$Period != 0)] <- 
    rainbow(20)[ treedata_s3273_pooled_gbtw_fg$Period ]
    
  colnames(treedata_s3273_pooled_gbtw_fg)[20] <- "color_rb"
  
  
  c4 <-ggtree(s3273_pooled_gbtw_fg_align) %<+% 
    treedata_s3273_pooled_gbtw_fg + aes(color = I(color_rb), alpha = 0.5 )
  
  # based on season and inter-epidemic
  # orange = inter-epi; blue: epi
  
  treedata_s3273_pooled_gbtw_fg[, 21] <- "#CCCCCC"
  treedata_s3273_pooled_gbtw_fg[, 21][which( treedata_s3273_pooled_gbtw_fg$Period != 0)] <- 
    rep(c("orange", "blue"), 11)[ treedata_s3273_pooled_gbtw_fg$Period ]
  
  colnames(treedata_s3273_pooled_gbtw_fg)[21] <- "color_epi"
  
  c1<-ggtree(s3273_pooled_gbtw_fg_align) %<+% 
    treedata_s3273_pooled_gbtw_fg + aes(color = I(color_epi), alpha = 0.5 )  
  
  # TW temporal distribution
  
  treedata_s3273_pooled_gbtw_fg[, 22] <- "#CCCCCC"
  
  treedata_s3273_pooled_gbtw_fg[, 22][which( treedata_s3273_pooled_gbtw_fg$TWb != 0)] <- 
    treedata_s3273_pooled_gbtw_fg$color_epi[which( treedata_s3273_pooled_gbtw_fg$TWb != 0)]
    
  colnames(treedata_s3273_pooled_gbtw_fg)[22] <- "color_epi_TW"
  
  treetipshape = c(rep(NA, 3273))
  treetipshape[s3273_TW_node[-117]] = 16
  treetipshape = data.frame(node=c(1:3273), treetipshape)
  
  
  c2 <- ggtree(s3273_pooled_gbtw_fg_align) %<+% 
    treedata_s3273_pooled_gbtw_fg + aes(color = I(color_epi_TW), alpha = 0.5 ) 
  s2 = c2 %<+% treetipshape + geom_tippoint(aes(shape = factor(shapee)), size = 3)
  
  multiplot(c4,c1,s2, t4, ncol = 4)
  
  
# extract distance to root #####  
  
  # s3273_TW_id
  # s3273_TW_time
  
  
       s3273_TW_node <- treedata_s3273_pooled_gbtw_fg$node[
         which(treedata_s3273_pooled_gbtw_fg$TWb == 1) ]
  
       # to deal with month "99"  
       x = treedata_s3273_pooled_gbtw_fg$Period[
         which(treedata_s3273_pooled_gbtw_fg$TWb == 1) ] 
       x[117] = 1
       
        s3273_TW_epi <- rep(c(1,0), 11)[x]
  
  s3273_TW_distoroot <- dist.nodes(s3273_pooled_gbtw_fg_align)[3274, s3273_TW_node]
  
  s3273_TW = data.frame(s3273_TW_id, s3273_TW_time, s3273_TW_epi, s3273_TW_distoroot)
  colnames(s3273_TW) = c("ID", "Time", "Period", "DisToRoot")
  
  s3273_TW <- s3273_TW[-117,]
  
  
  # plot: time against distance to root
  ggplot(s3273_TW, aes(Time, DisToRoot, color = factor(Period)) ) + 
    geom_point() + 
    theme_bw() + 
    
    geom_point(color = "black", size = 3.2, alpha = 0.8, shape = 1) + 
    geom_point(size = 3, alpha = 0.8) + 
    
    theme(
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.title = element_text(face="bold"),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x = element_text(size = 15), 
      axis.text.y = element_text(size = 15), 
      
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 20), 
      legend.position="top") +
    
    xlab("") + ylab("Root-to-Tip Distance") +
    scale_x_continuous(breaks = seq(2007, 2016, by = 1) ) + 
    scale_color_discrete(name = "Period", labels=c("Epi", "Inter-Epi"))
  
  
# extract 127 seq from the original file ####
  
  library(seqinr) 
  
  # read-in: pooled_gbtw_fg_trim.fasta
      fasta_s3273_tree <- read.fasta(file.choose())
  fasta_seq_s3273_tree <- getSequence(fasta_s3273_tree)
   fasta_ID_s3273_tree <- attributes(fasta_s3273_tree)$names
    
   fasta_noTW_s3273_tree <- sort(match(s3273_TW_id[-117], fasta_ID_s3273_tree))
  
  # which(treedata_s3273_TWb == 1) n = 128
 
  write.fasta(fasta_seq_s3273_tree[fasta_noTW_s3273_tree], 
              file.out = "s3273_TW.fasta", 
              names = fasta_ID_s3273_tree[fasta_noTW_s3273_tree])
  
  
# BEAST outcome ####
  
  s3273_tw_mcc <- read.beast(file.choose())
  ggtree(s3273_tw_mcc, mrsd="2016-03-05") + theme_tree2()
  
  # s3273_TW
  
  treedata_s3273_TW_mcc = fortify(s3273_tw_mcc)
  
       t_color <- match(treedata_s3273_TW_mcc$label, s3273_TW$ID)[1:127]
  t_color_tree <- c("blue","orange")[s3273_TW$Period[t_color]+1]
  
  g = ggtree(s3273_tw_mcc, mrsd="2016-03-05") + theme_tree2() + 
    geom_range(range='length_0.95_HPD', color='red', alpha=.4, size=0.5) +
    geom_tippoint(color = t_color_tree) +
    scale_x_continuous(breaks = seq(2007, 2016, by=1), 
                       minor_breaks = seq(2007, 2016, by = 0.5)) +
    theme_tree2(panel.grid.major.x = element_line(color = "black", size = 0.2),
          panel.grid.minor.x = element_line(color = "grey", size = 0.1),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank() )
    
  # skyride
  
  skyridedata = read.csv(file.choose())
  
  g + geom_line(data = skyridedata, aes(x= Time, y=Median*10), 
                color = "gray", size = 1.5, alpha = 0.6) +
    
    geom_ribbon(data = skyridedata, aes(x = Time, ymin = Lower*10, ymax = Upper*10),  
                inherit.aes=FALSE, fill = "gray", alpha = 0.2) 
  

  
  
  
  