# Aim 2 functions

# Fit models with robust standard errors (practice-level cluster)
sglm <- function(d){
  # extract parameters and covariance matrix
  betas <- lapply(d, FUN=function(rr){ coef(rr) } )
  vars <- lapply(d, FUN=function(rr){ vcov(rr) } )
  
  
  summary( miceadds::pool_mi( qhat=betas, u=vars ) )
}

# Generate predicted probabilities
glmpp <- function(m, d, cp = F, sts = F, cvar = "cancer"){
  d1 <- d[[1]]
  nimp <- m %>% length()
  
  # Create a data matrix with potential combination of values we want to explore.
  # Also set up a tibble to accommodate those combinations where we can store
  # probabilities.
  sdat <- tibble(est.janmar2019 = numeric(), 
                 lci.janmar2019 = numeric(), 
                 uci.janmar2019 = numeric(),
                 
                 est.aprjul2019 = numeric(), 
                 lci.aprjul2019 = numeric(), 
                 uci.aprjul2019 = numeric(), 
                 
                 est.janmar2020 = numeric(), 
                 lci.janmar2020 = numeric(), 
                 uci.janmar2020 = numeric(),
                 
                 est.aprjul2020 = numeric(), 
                 lci.aprjul2020 = numeric(), 
                 uci.aprjul2020 = numeric(),
                 
                 est.d2019 = numeric(), 
                 lci.d2019 = numeric(), 
                 uci.d2019 = numeric(),
                 p.d2019 = numeric(),
                 
                 est.d2020 = numeric(), 
                 lci.d2020 = numeric(),
                 uci.d2020 = numeric(),
                 p.d2020 = numeric(),
                 
                 est.dd = numeric(), 
                 lci.dd = numeric(),
                 uci.dd = numeric(),
                 p.dd = numeric())
  
  if(!cp){
    # Set # possible combinations (marginal, 1)
    combos <- 1
  } else { # for conditional (cancer) probabilities use this one
    if(cvar == "cancer"){
      yp <- d1 %>% tidyr::expand(cnsc, cpro, cucc) %>%
        filter(cnsc + cpro + cucc < 2)
      
      cvar.t <- tibble(cnsc = numeric(),
                       cpro = numeric(), 
                       cucc = numeric())
      
      # race/ethnicity 
    } else if(cvar == "reth"){
      yp <- d1 %>% tidyr::expand(rblack, rhisp, rother) %>%
        filter(rblack + rhisp + rother < 2)
      
      cvar.t <- tibble(rblack = numeric(), rhisp = numeric(), rother = numeric())
      
      # denovo metastatic status
    } else if(cvar == "denov"){
      yp <- d1 %>% tidyr::expand(denovo_met)
      
      cvar.t <- tibble(denovo_met = factor())
      
      # insurance type
    } else if(cvar == "insur"){
      yp <- d1 %>% tidyr::expand(igov, iother) %>%
        filter(igov + iother < 2)
      
      cvar.t <- tibble(igov = numeric(), iother = numeric())
      
      # age
    } else if(cvar == "age"){
      yp <- d1 %>% tidyr::expand(agegt75)
      
      cvar.t <- tibble(agegt75 = numeric())
    }
    # Set # possible combinations
    combos <- nrow(yp)
  }
  
  
  for(k in 1:combos){ # loop through potential combinations
    sp.store <- tibble(
      est.janmar2019 = numeric(), 
      se.janmar2019 = numeric(), 
      
      est.aprjul2019 = numeric(), 
      se.aprjul2019 = numeric(), 
      
      est.janmar2020 = numeric(), 
      se.janmar2020 = numeric(), 
      
      est.aprjul2020 = numeric(), 
      se.aprjul2020 = numeric(), 
      
      est.d2019 = numeric(), 
      se.d2019 = numeric(), 
      
      est.d2020 = numeric(), 
      se.d2020 = numeric(),
      
      est.dd = numeric(), 
      se.dd = numeric() 
    )
    
    sp.pooled <- sp.store
    
    
    
    for(i in 1:length(m)){ # loop through each imputation
      # Grab all model values in imputed data i where our combination set is met.
      tt <- d[[i]]
      
      # Create TF vector for subset and get n in subset
      # Since no years/cancers/period obs are missing this will be the same
      # in every imputation
      if(!cp){
        subtf <- T
      } else {
        if(cvar == "cancer"){
          subtf <- ifelse(tt$cnsc == yp$cnsc[k] &
                            tt$cpro == yp$cpro[k] &
                            tt$cucc == yp$cucc[k], T, F)
        } else if(cvar == "reth"){
          subtf <- ifelse(tt$rblack == yp$rblack[k] &
                            tt$rhisp == yp$rhisp[k] &
                            tt$rother == yp$rother[k], T, F)
        } else if(cvar == "denov"){
          subtf <- ifelse(tt$denovo_met == yp$denovo_met[k], T, F)
        } else if(cvar == "insur"){
          subtf <- ifelse(tt$igov == yp$igov[k] &
                            tt$iother == yp$iother[k], T, F)
        } else if(cvar == "age"){
          subtf <- ifelse(tt$agegt75 == yp$agegt75[k], T, F)
        }
        
      }
      
      tt$sub <- subtf
      
      # Calculate predicted probs for imputation i across study periods.
      # We can directly use the model fit with cluster specified as 
      # stdGlm uses the same approach as glm.cluster to calculate robust
      # standard errors.
      sf <- stdGlm(m[[i]]$glm_res, data=tt, X="studypf", clusterid = "practiceid", 
                   subsetnew = sub)
      
      covm <- sf$vcov
      estv <- as.vector(sf$est)
      sev <- diag(covm)
      
      ## Get estimates and standard errors for diffs and diff in diffs
      # Diff 1: -janmar2019 + aprjul2020
      d1v <- c(-1, 1, 0, 0)
      est.d2019 <- estv[2] - estv[1]
      se.d2019 <- d1v %*% covm %*% d1v
      
      # Diff 2: -janmar2020 + aprjul2020
      d2v <- c(0, 0, -1, 1)
      est.d2020 <- estv[4] - estv[3]
      se.d2020 <- d2v %*% covm %*% d2v
      
      # Diff 3: -(-janmar2019 + aprjul2019) + (-janmar2020 + aprjul2020) 
      # = janmar2019 - aprjul2019 - janmar2020 + aprjul2020
      d3v <- c(1, -1, -1, 1)
      est.dd <- est.d2020 - est.d2019
      se.dd <- d3v %*% covm %*% d3v
      
      # Store estimates and variances 
      sp.store <- sp.store %>% 
        add_row(
          est.janmar2019 = estv[1], 
          se.janmar2019 = sev[1], 
          
          est.aprjul2019 = estv[2], 
          se.aprjul2019 = sev[2], 
          
          est.janmar2020 = estv[3], 
          se.janmar2020 = sev[3], 
          
          est.aprjul2020 = estv[4], 
          se.aprjul2020 = sev[4], 
          
          est.d2019 = est.d2019, 
          se.d2019 = se.d2019, 
          
          est.d2020 = est.d2020, 
          se.d2020 = se.d2020,
          
          est.dd = est.dd, 
          se.dd = se.dd 
        )
      print(paste("Completed iteration", (k-1)*nimp + i, "of", nimp*combos))
    } # end imputation loop
    
    # pool results
    p <- NULL
    pval <- NULL
    subn <- sum(tt$sub)
    sp.store <- sp.store %>% as.data.frame()
    for(i in 1:(ncol(sp.store)/2)){
      ind <- i*2
      pf <- pool.scalar(sp.store[,ind-1], sp.store[,ind])
      #, n = subn, 
      #k=length(m[[1]]$coefficients))
      
      # Invert probabilities if sts = T
      if(!sts){
        e <- pf$qbar
      } else if(i > 4){
        e <- -pf$qbar
      } else {
        e <- 1 - pf$qbar
      }
      
      lci <- e - 1.96*sqrt(pf$t)
      uci <- e + 1.96*sqrt(pf$t)
      
      # Calculate p-values for diffs and diff in diffs
      if(i>4){
        pval <- c(pval, 2*pt(-abs(e/sqrt(pf$t)), nrow(tt)))
      }
      
      # set lower limit of CI to 0
      if(lci < 0 & i <= 4){lci <- 0}
      
      # add data to vector
      p <- c(p, e, lci, uci)
    }
    
    sdat <- sdat %>% add_row( 
      est.janmar2019 = p[1], 
      lci.janmar2019 = p[2], 
      uci.janmar2019 = p[3],
      
      est.aprjul2019 = p[4], 
      lci.aprjul2019 = p[5], 
      uci.aprjul2019 = p[6], 
      
      est.janmar2020 = p[7], 
      lci.janmar2020 = p[8], 
      uci.janmar2020 = p[9],
      
      est.aprjul2020 = p[10], 
      lci.aprjul2020 = p[11], 
      uci.aprjul2020 = p[12],
      
      est.d2019 = p[13], 
      lci.d2019 = p[14], 
      uci.d2019 = p[15],
      p.d2019 = pval[1],
      
      est.d2020 = p[16], 
      lci.d2020 = p[17],
      uci.d2020 = p[18],
      p.d2020 = pval[2],
      
      est.dd = p[19], 
      lci.dd = p[20],
      uci.dd = p[21],
      p.dd = pval[3])
    
    # Store probabilities in tibble
    if(cp){
      if(cvar == "cancer"){
        cvar.t <- cvar.t %>% add_row(cnsc = yp$cnsc[k],
                                     cpro = yp$cpro[k], 
                                     cucc = yp$cucc[k])
      } else if(cvar == "reth"){
        cvar.t <- cvar.t %>% add_row(rblack = yp$rblack[k], rhisp = yp$rhisp[k], 
                                     rother = yp$rother[k])
      } else if(cvar == "denov"){
        cvar.t <- cvar.t %>% add_row(denovo_met = yp$denovo_met[k])
      } else if(cvar == "insur"){
        cvar.t <- cvar.t %>% add_row(igov = yp$igov[k], iother = yp$iother[k])
      } else if(cvar == "age"){
        cvar.t <- cvar.t %>% add_row(agegt75 = yp$agegt75[k])
      }
    } 
  }
  
  if(cp){
    sdat.f <- bind_cols(sdat, cvar.t)
  } else{
    sdat.f <- sdat
  }
  
  
  return(sdat.f)
}
