##### 1 [ML]: perfect test / Journal of Aquatic Animal Health 17:386-391, 2005 /
# http://www.webpages.uidaho.edu/~chrisw/research/prevalence/
#
# Function:
#


llprevr <- function(p, yes = c(0), no = c(0)){
  sumcheck <- sum(yes) + sum(no)
  if (sumcheck == 0) stop("Data must be entered for yes or no")
  if (sum(yes) == 0) { llmck <- sum(no) * log(1-p) }
  else { llmck <- sum(log(1-(1-p)^yes)) + sum(no) * log(1-p) }
  
  llmck
}
dprev <- function(yes = c(0),no = c(0),disp = "y",conf = .95){
  sumcheck <- sum(yes) + sum(no)
  if (sumcheck == 0) stop("Data must be entered for yes or no")
  if (sum(yes) == 0)
  {
    ucl <- 1 - exp(-qchisq(conf,1)/(2*sum(no)))
    result <- c(0., 0., ucl)
    if (disp == 'y') print("Lower 95% limit, MLE, Upper 95% limit = ")
    result
  }
  else if (sum(no) == 0)
  {
    tfct <- function(p)
    {
      sum(log(1-(1-p)^yes)) + qchisq(conf,1)/2
    }
    lcl <- uniroot(tfct, interval = c(0.0001, 1.))
    result <- c(lcl$root, 1., 1.)
    if (disp == "y") print("Lower 95% limit, MLE, Upper 95% limit = ")
    result
  }
  else
  {
    #print(yes) ; #print(no)
    llpmax <- optimize(llprevr, c(0., 1.), maximum = TRUE, yes = yes, no = no)
    #print(llpmax)
    
    mval <- llprevr(llpmax$maximum, yes = yes, no = no)
    tfct <- function(p)
    {
      llprevr(p, yes, no)-(mval-qchisq(conf,1)/2)
    }
    lcl <- uniroot(tfct, interval = c(0.00000001, llpmax$maximum))
    ucl <- uniroot(tfct, interval = c(llpmax$maximum, 0.99999999))
    result <- c(lcl$root, llpmax$maximum, ucl$root)
    if (disp == "y") print("Lower 95% limit, MLE, Upper 95% limit = ")
    result
  }
}


# Example:
#

# 4 test positive pools
# 16 test negative pools
# each pool has 5 samples

p = rep(5,4)
# [1] 5 5 5 5

n = rep(5,16)
# [1] 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5


dprev( y = c(p), n=c(n) )

# [1] "Lower 95% limit, MLE, Upper 95% limit = "
# [1] 0.01373742 0.04364423 0.09873630




##### 2 [Baysiean w MCMC]: imperfect test / Ecological Informatics 5:273-280, 2010
# http://www.webpages.uidaho.edu/~chrisw/research/prevalence/
#
# Function:
#


post1mat <- function(data,pia,pib,etaa,etab,lama,lamb,pi,eta,lambda){
  samples <- dim(data)[1]
  posprob <- numeric(samples)
  negprob <- numeric(samples)
  pos <- data[,3]
  neg <- data[,2]
  r   <- data[,1]
  llpart <- numeric(samples)
  for (i in 1:samples)
  {
    posprob[i] <- eta*(1-(1-pi)^r[i]) + (1-lambda)*(1-pi)^r[i]
    negprob[i] <- (1-eta)*(1-(1-pi)^r[i]) + lambda*(1-pi)^r[i]
    llpart[i] <- neg[i]*log(negprob[i]) + pos[i]*log(posprob[i]) 
  }
  lres <- sum(llpart) + (pia-1)*log(pi) + (pib-1)*log((1-pi)) +
    (etaa-1)*log(eta) + (etab-1)*log((1-eta)) + (lama-1)*log(lambda) + (lamb-1)*log((1-lambda))
  #  print("Density value =  ") ; print(exp(sum(llpart))) ;
  res  <- exp(lres)
  res
}
compenummat1a <- function(data,pia,pib,etaa,etab,lama,lamb){
  begin <- proc.time()
  postv <- array(0, dim=c(199,99,99))
  for (pi in seq(.005,.995,by=.005) )
  {  
    for (eta in seq(.01,.99,by=.01) )
    {
      for (lambda in seq(.01,.99,by=.01) )
      {
        i <- round((pi -.005)/.005 + 1)
        j <- round((eta -.01)/.01 + 1)
        k <- round((lambda -.01)/.01 + 1)
        postv[i,j,k] <- post1mat(data,pia,pib,etaa,etab,lama,lamb,pi,eta,lambda)
      } 
    }
  }
  pimarg <- apply(postv, c(1), sum)
  pimargnorm <- pimarg/sum(pimarg)
  print(dim(postv)) ; print(length(pimarg)) 
  print(cumsum(pimargnorm))
  plot(seq(.005,.995,by=.005),pimarg,type="l")
  # windows() ; 
  plot(seq(.005,.995,by=.005),pimargnorm,type="l",xlab="Prevalence",ylab="Posterior Density")
  # windows() ; 
  plot(1:80/200,pimargnorm[1:80],type="l")
  quant <- c(.005,.01,.025,.05,.10,.25,.5,.75,.90,.95,.975,.99,.995)
  qres <- approx(x=cumsum(pimargnorm),y=seq(.005,.995,by=.005),xout=quant)
  pimean <- t(seq(.005,.995,by=.005))%*%pimargnorm
  pivar <- t((seq(.005,.995,by=.005)-pimean)^2)%*%pimargnorm
  pistd <- sqrt(pivar)
  print(qres)
  print("Mean and SD for pi") ; print(c(pimean,pistd))
  end <- proc.time()
  print("total time in minutes is "); print( (end[3]-begin[3])/60) 
}

# Example:
#

# 3 test positive pools ( n = 5 )
# 1 test positive pools ( n = 6 )
# 15 test negative pools ( n = 5 )
#
# '''' design the alpha & beta of beta-dist.
#
# pia = prevalence alph (ex. 1.)
# pib = prevalence beta (1.)
# etaa = sensitivity alph (24)
# etab = sensitivity beta (6)
# lama = specificity alph (81)
# lama = specificity beta (9)
#
#
#   size  NO.ne No.po
#     5    15     3   
#     6     0     1
#     
# For Beta-dist
# prior mean = a/(a+b)
# observation = a + b -2 
# use the following function with m = mean ; n = observation 

beta.ab<-function(m, n){
  
  Valpha = m*(n + 2)
  Vbeta = n - Valpha + 2 
  
  r = c(Valpha, Vbeta)
  
  print("alpha, beta:")
  print(r)
  
}      

# beta.ab(0.9, 1000)

# a=beta.ab(0.9, 1000)[1]
# b=beta.ab(0.9, 1000)[2]

# x<- rbeta(100000, a, b)
# plot(density(x))          # hist(x)


mx = matrix(c(5,15,3,6,0,1), nrow = 2, ncol = 3, byrow = T)

pia = 1.
pib = 1.
etaa = 902
etab = 100
lama = 902
lamb = 100

compenummat1a(mx,pia,pib,etaa,etab,lama,lamb)


######################################## refine the functions to preML and preBay
#
#  Example:
#  
#   size  NO.Po No.Ne
#     5    3     15   
#     6    1      0
#
# 

test = matrix(c(5,3,15,6,1,0), nrow = 2, ncol = 3, byrow = T)

###### preML ###### 
#

preML<-function(mx){
  
  
  # preML modified from function llprevr and dprev (source: Journal of Aquatic Animal Health 17:386???391, 2005)
  
  
  llprevr <- function(p, yes = c(0), no = c(0)){
    sumcheck <- sum(yes) + sum(no)
    if (sumcheck == 0) stop("Data must be entered for yes or no")
    if (sum(yes) == 0) { llmck <- sum(no) * log(1-p) }
    else { llmck <- sum(log(1-(1-p)^yes)) + sum(no) * log(1-p) }
    
    llmck
  }
  
  disp = "y"
  conf = 0.95
  
  d = dim(mx)
  
  yes = c()
  no = c()
  
  
  for(i in 1: d[1] ){
    
    yy = rep(mx[i, 1], mx[i, 2])
    nn = rep(mx[i, 1], mx[i, 3])
    
    yes = c(yes, yy)
    no = c(no, nn)
    
  }
  
  
  sumcheck <- sum(yes) + sum(no)
  if (sumcheck == 0) stop("Data must be entered for yes or no")
  if (sum(yes) == 0)
  {
    ucl <- 1 - exp(-qchisq(conf,1)/(2*sum(no)))
    result <- c(0., 0., ucl)
    if (disp == 'y') print("Lower 95% limit, MLE, Upper 95% limit = ")
    result
  }
  else if (sum(no) == 0)
  {
    tfct <- function(p)
    {
      sum(log(1-(1-p)^yes)) + qchisq(conf,1)/2
    }
    lcl <- uniroot(tfct, interval = c(0.0001, 1.))
    result <- c(lcl$root, 1., 1.)
    if (disp == "y"){ 
      print("Lower 95% limit, MLE, Upper 95% limit = ")}
    return(result)
    
  }
  else
  {
    #print(yes) ; #print(no)
    llpmax <- optimize(llprevr, c(0., 1.), maximum = TRUE, yes = yes, no = no)
    #print(llpmax)
    
    mval <- llprevr(llpmax$maximum, yes = yes, no = no)
    tfct <- function(p)
    {
      llprevr(p, yes, no)-(mval-qchisq(conf,1)/2)
    }
    lcl <- uniroot(tfct, interval = c(0.00000001, llpmax$maximum))
    ucl <- uniroot(tfct, interval = c(llpmax$maximum, 0.99999999))
    result <- c(lcl$root, llpmax$maximum, ucl$root)
    if (disp == "y") {
      print ("Lower 95% limit, MLE, Upper 95% limit = ")}
    return(result)
    
  }
  
  
  
  
  
  
  
}


#
#         [,1] [,2] [,3]
#   [1,]    5    3   15
#   [2,]    6    1    0
#

preML(test)

#   [1] "Lower 95% limit, MLE, Upper 95% limit = "
#   [1] 0.01446416 0.04592731 0.10381140  



###### preBay ######
#
# no idea for prior (ex: prevalence), put 0.5 and 0 for Pm and Pn (noinformative prior)
# 
# Pm: mean for prevalence
# Pn: observations
# Sem: mean for sensitivity
# Sen: observations
# Spm: mean for specificity
# Spn: observations
#


preBay<-function(data, Pm, Pn, Sem, Sen, Spm, Spn){
  
  # preBay modified from function post1mat and compenummat1a (source: Ecological Informatics 5:273-280, 2010)
  
  beta.ab<-function(m, n){
    
    Valpha = m*(n + 2)
    Vbeta = n - Valpha + 2 
    
    r = c(Valpha, Vbeta)
    
    return(r)
    
  }      
  
  post1mat <- function(data,pia,pib,etaa,etab,lama,lamb,pi,eta,lambda){
    samples <- dim(data)[1]
    posprob <- numeric(samples)
    negprob <- numeric(samples)
    pos <- data[,2]
    neg <- data[,3]
    r   <- data[,1]
    llpart <- numeric(samples)
    for (i in 1:samples)
    {
      posprob[i] <- eta*(1-(1-pi)^r[i]) + (1-lambda)*(1-pi)^r[i]
      negprob[i] <- (1-eta)*(1-(1-pi)^r[i]) + lambda*(1-pi)^r[i]
      llpart[i] <- neg[i]*log(negprob[i]) + pos[i]*log(posprob[i]) 
    }
    lres <- sum(llpart) + (pia-1)*log(pi) + (pib-1)*log((1-pi)) +
      (etaa-1)*log(eta) + (etab-1)*log((1-eta)) + (lama-1)*log(lambda) + (lamb-1)*log((1-lambda))
    #  print("Density value =  ") ; print(exp(sum(llpart))) ;
    res  <- exp(lres)
    res
  }
  
  
  pia = beta.ab(Pm, Pn)[1]
  pib = beta.ab(Pm, Pn)[2]
  etaa = beta.ab(Sem, Sen)[1]
  etab = beta.ab(Sem, Sen)[2]
  lama = beta.ab(Spm, Spn)[1]
  lamb = beta.ab(Spm, Spn)[2]
  
  postv <- array(0, dim=c(199,99,99))
  
  for (pi in seq(.005,.995,by=.005) ){  
    for (eta in seq(.01,.99,by=.01) )
    {
      for (lambda in seq(.01,.99,by=.01) )
      {
        i <- round((pi -.005)/.005 + 1)
        j <- round((eta -.01)/.01 + 1)
        k <- round((lambda -.01)/.01 + 1)
        postv[i,j,k] <- post1mat(data,pia,pib,etaa,etab,lama,lamb,pi,eta,lambda)
      } 
    }
  }
  
  pimarg <- apply(postv, c(1), sum)
  pimargnorm <- pimarg/sum(pimarg)
  
  plot(seq(.005,.995,by=.005),pimargnorm,type="l",xlab="Prevalence",ylab="Posterior Density")
  
  quant <- c(.005,.01,.025,.05,.10,.25,.5,.75,.90,.95,.975,.99,.995)
  qres <- approx(x=cumsum(pimargnorm),y=seq(.005,.995,by=.005),xout=quant)
  pimean <- t(seq(.005,.995,by=.005))%*%pimargnorm
  
  pivar <- t((seq(.005,.995,by=.005)-pimean)^2)%*%pimargnorm
  pistd <- sqrt(pivar)
  
  
  ## 
  
  picd <- cumsum(pimargnorm)
  sp = seq(.005,.995,by=.005)
  
  # set 95% interval
  
  # 1 Baysiean Credibility Interval 
  
  if (picd[1] > 0.025){
    
    lo.bci = approx(c(0, picd[1]), c(0, sp[1]), xout = 0.025)$y
    
  }else{
    
    lo.bci = approx( c( picd[min(which(picd > 0.025))-1],  picd[min(which(picd > 0.025))]), 
                     c(   sp[min(which(picd > 0.025))-1] ,   sp[min(which(picd > 0.025))]), xout = 0.025)$y
    
    hi.bci = approx( c( picd[min(which(picd > 0.975))-1],  picd[min(which(picd > 0.975))]), 
                     c(   sp[min(which(picd > 0.975))-1] ,   sp[min(which(picd > 0.975))]), xout = 0.975)$y
    
  }
  
  
  # 2 High Density Interval 
  
  
  id.sp = which(picd < 0.05)
  
  if(length(id.sp) == 0 ){
    
    lo.hdi = 0
    hi.hdi = 0
    
  }else{
    
    id.sp.k = sapply(id.sp, function(x){
      y <- min(which(picd-picd[x] > 0.95))
      return(y)
    })
    
    lo.hdi = sp[id.sp[which.min(id.sp - id.sp.k)]]
    hi.hdi = sp[id.sp.k[which.min(id.sp - id.sp.k)]] }
  
  
  # 3 Median
  
  pimedian<-approx( c(  picd[min(which(picd > 0.5)) -1 ],   picd[min(which(picd > 0.5))] ), 
                    c(  sp[min(which(picd > 0.5)) -1 ], sp[min(which(picd > 0.5))] ), xout = 0.5 )$y
  
  
  result = c(pimean, pistd, pimedian, lo.bci, hi.bci, lo.hdi, hi.hdi)
  result = matrix(result, ncol=7, byrow = TRUE)
  colnames(result) = c("mean", "sd", "median", "Lower CI", "Higher CI", "Lower HPD", "Higher HPD")
  return(result)
}

#
# pm = 0.5
# pn = 0
# Sem = 0.90 
# Sen = 100
# Spm = 0.9
# Sen = 100
#

preBay(test, 0.5, 0, 0.9, 100, 0.9, 100)

#          mean         sd     median    Lower CI Higher CI Lower HPD Higher HPD
# [1,] 0.045849 0.02855483 0.03854346 0.003192145 0.1094946     0.005       0.13


###


test2 = matrix(c(5,15,3,6,0,1), nrow = 2, ncol = 3, byrow = T)


#   size  NO.Po No.Ne
#     5    15     3   
#     6     0     1



preBay(test2, 0.5, 0, 0.9, 100, 0.9, 100)


#           mean        sd    median  Lower CI Higher CI Lower HPD Higher HPD
# [1,] 0.5369014 0.2316967 0.4904218 0.1953796 0.9681858     0.005      0.945


preML(test2)

# [1] "Lower 95% limit, MLE, Upper 95% limit = "
# [1] 0.1546105 0.2621143 0.4026503



######################################## sensitivity test with preML and preBay

D = c()
for (i in 1: 10){
  
  assign(paste0("mmx",i), matrix( c(5,i*10,100-i*10), nrow = 1, ncol = 3, byrow = T ))
  
  mm = get(paste0("mmx",i))
  
  d = preML(mm)[2] - preBay(mm, 0.5, 0, 0.9, 100, 0.9, 100)[1] 
  
  D[length(D) + 1 ] = d
  
  print(i)
}


plot(D, xlab="positive sample (*10)", ylab = "E(preML) - E(preBay)")

######################################## input via .csv file 

csvf = read.csv(file.choose(),header = FALSE)
csvf = as.matrix(csvf)

# > csvf
# V1 V2 V3
# [1,]  5  5 35
# [2,] 10  1  5

preML(csvf)

# [1] "Lower 95% limit, MLE, Upper 95% limit = "
# [1] 0.009813239 0.024497547 0.049016848

preBay(csvf, 0.5, 0, 0.9, 100, 0.9, 100)

#            mean         sd     median    Lower CI Higher CI Lower HPD Higher HPD
# [1,] 0.01728512 0.01098294 0.01260368 0.000628662 0.1094946         0          0



