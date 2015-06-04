logMeanExpLogs = function(v)
{
	N = length(v)
	.Call("RLogMeanExpLogs", as.numeric(v), N, package="BayesSingleSub")
}



matrixDet = function(A, N, returnLog){
  
  out = .Call("matrixDetCall",          
              as.numeric(A),
              as.integer(N),
              as.integer(returnLog),
              package="BayesSingleSub")
  
  return(out)
  
}

expXphiAminusphiB = function(x, a, b, returnLog){
  
  out = .Call("expXphiAminusphiBCall",          
              as.numeric(x),
              as.numeric(a),
              as.numeric(b),
              as.integer(returnLog),
              package="BayesSingleSub")
  
  return(out)
  
}



log1mp = function(logp){
  
  out = .Call("log1mpCall",          
              as.numeric(logp),
              package="BayesSingleSub")
  
  return(out)
  
}

decimalplaces = function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}


findlogPriorp = function(iter, test, a, b, r1, r2){
  
  N = length(test)
  
  logPriorp = .Call("findlogPriorpCall",          
                    as.integer(iter),
                    as.integer(N),
                    as.integer(test),
                    as.numeric(a),
                    as.numeric(b),
                    as.numeric(r1),
                    as.numeric(r2),
                    package="BayesSingleSub")
  
  return(logPriorp)
  
}


findlogPriorpIT = function(iter, test, aB1, bB1, aB3, bB3, r1, r2, r3, r4){
  
  N = length(test)
  
  logPriorp = .Call("findlogPriorpITCall",          
                    as.integer(iter),
                    as.integer(N),
                    as.integer(test),
                    as.numeric(aB1),
                    as.numeric(bB1),
                    as.numeric(aB3),
                    as.numeric(bB3),
                    as.numeric(r1),
                    as.numeric(r2),
                    as.numeric(r3),
                    as.numeric(r4),
                    package="BayesSingleSub")
  
  return(logPriorp)
  
}


func_BF.nr.changed = function(N,logp0, logp1, a, b, r1, r2, iter){
  
  log.pchanged = rep(-100000,N+1)
  log.prior.pchanged = rep(-100000,N+1)
  prior.odds = vector(length=N)
  post.odds = vector(length=N)
  
  pb = txtProgressBar(min = 0, max = (N-1), style = 3)
  message("pb 2/2")
  
  for(i in 0:(2^N-1)){
    
    setTxtProgressBar(pb, i)
        
    for(j in 0:N){      
      
      if(sum(rev(as.numeric(intToBits(i)[1:N])))==j){
        
        log.pchanged[j+1] = log(2) + logMeanExpLogs(c(log.pchanged[j+1],
                                                      
                                                      logMeanExpLogs(rowSums(cbind(
                                                        logp1*matrix(rev(as.numeric(intToBits(i)[1:N])),ncol=N,nrow=dim(logp1)[1], byrow=TRUE), 
                                                        logp0*(1- matrix(rev(as.numeric(intToBits(i)[1:N])),ncol=N,nrow=dim(logp0)[1],byrow=TRUE))
                                                      )))
        ))
        
        log.prior.pchanged[j+1] = log(2) + logMeanExpLogs(c(log.prior.pchanged[j+1],
                                                            
                                                            findlogPriorp(iter, rev(as.numeric(intToBits(i)[1:N])), a, b, r1, r2)
                                                            
                                                            
        ))
      }
      
    }
  }
  
  close(pb)
  
  # Log odds
  log.post.odds = vector(length=N)
  log.prior.odds = vector(length=N)
  
  for(i in 2:(N+1)){
    log.post.odds[i-1] =  (log(length(i:(N+1))) + logMeanExpLogs(log.pchanged[i:(N+1)])) - 
      (log(length(1:(i-1))) + logMeanExpLogs(log.pchanged[1:(i-1)]))
    
    log.prior.odds[i-1] = (log(length(i:(N+1))) + logMeanExpLogs(log.prior.pchanged[i:(N+1)])) - 
      (log(length(1:(i-1))) + logMeanExpLogs(log.prior.pchanged[1:(i-1)]))
  }
  
  logBF = log.post.odds - log.prior.odds
  
  return(logBF)     
  
}



func_BF.nr.changedTrendJAGS = function(N, iterchains, muB1, sig2B1, muB3, sig2B3, mu_sig2_e, 
                                       aB1, bB1, aB3, bB3, r1, r2, r3, r4, iter){
  
  logp0Int = matrix(ncol=N, nrow = iterchains)
  logp1Int = matrix(ncol=N, nrow = iterchains)
  logp0Trend = matrix(ncol=N, nrow = iterchains)
  logp1Trend = matrix(ncol=N, nrow = iterchains)
  logp0IT = matrix(ncol=N, nrow = iterchains)
  logp1IT = matrix(ncol=N, nrow = iterchains)
  
  for(k in 1:iterchains){
    
    for(i in 1:N){
      
      left1 = (aB1 - muB1[k,i] / sqrt(exp(mu_sig2_e[k]))) /
        sqrt(sig2B1[k,i] /exp(mu_sig2_e[k]));
      
      right1 = (bB1 - muB1[k,i] / sqrt(exp(mu_sig2_e[k]))) /
        sqrt(sig2B1[k,i] /exp(mu_sig2_e[k]));
      
      left2 = (aB3 - muB3[k,i] / sqrt(exp(mu_sig2_e[k]))) /
        sqrt(sig2B3[k,i] /exp(mu_sig2_e[k]));
      
      right2 = (bB3 - muB3[k,i]/ sqrt(exp(mu_sig2_e[k]))) /
        sqrt(sig2B3[k,i] /exp(mu_sig2_e[k]));
      
      # Conditional posterior probabilities that subject i changed
      # Intercept
      logp0Int[k,i] = max(c(expXphiAminusphiB(0, right1, left1, 1), -1000000));
      logp1Int[k,i]= max(c(log1mp(logp0Int[k,i]), -1000000));
      
      #Trend
      logp0Trend[k,i] = max(c(expXphiAminusphiB(0, right2, left2, 1), -1000000));
      logp1Trend[k,i] = max(c(log1mp(logp0Trend[k,i]), -1000000));
      
      #Intercept + trend
      logp0IT[k,i] = logp0Int[k,i] + logp0Trend[k,i];
      logp1IT[k,i] = logp1Int[k,i] + logp1Trend[k,i];
      
    }
    
  }
  
  
  log.pchangedInt = rep(-100000,N+1)
  log.prior.pchangedInt = rep(-100000,N+1)
  log.pchangedTrend = rep(-100000,N+1)
  log.prior.pchangedTrend = rep(-100000,N+1)
  log.pchangedIT = rep(-100000,N+1)
  log.prior.pchangedIT = rep(-100000,N+1)
  
  prior.oddsInt = vector(length=N)
  post.oddsInt = vector(length=N)
  prior.oddsTrend = vector(length=N)
  post.oddsTrend = vector(length=N)
  prior.oddsIT = vector(length=N)
  post.oddsIT = vector(length=N)
  
  pb = txtProgressBar(min = 0, max = (N-1), style = 3)
  message("pb 2/2")
  
  for(i in 0:(2^N-1)){
    
    setTxtProgressBar(pb, i)
    
    for(j in 0:N){      
      
      if(sum(rev(as.numeric(intToBits(i)[1:N])))==j){
        
        log.pchangedInt[j+1] = log(2) + logMeanExpLogs(c(log.pchangedInt[j+1],
                                                         
                                                         logMeanExpLogs(rowSums(cbind(
                                                           logp1Int*matrix(rev(as.numeric(intToBits(i)[1:N])),ncol=N,nrow=dim(logp1Int)[1], byrow=TRUE), 
                                                           logp0Int*(1- matrix(rev(as.numeric(intToBits(i)[1:N])),ncol=N,nrow=dim(logp0Int)[1],byrow=TRUE))
                                                         )))
        ))
        
        log.pchangedTrend[j+1] = log(2) + logMeanExpLogs(c(log.pchangedTrend[j+1],
                                                           
                                                           logMeanExpLogs(rowSums(cbind(
                                                             logp1Trend*matrix(rev(as.numeric(intToBits(i)[1:N])),ncol=N,nrow=dim(logp1Trend)[1], byrow=TRUE), 
                                                             logp0Trend*(1- matrix(rev(as.numeric(intToBits(i)[1:N])),ncol=N,nrow=dim(logp0Trend)[1],byrow=TRUE))
                                                           )))
        ))
        
        log.pchangedIT[j+1] = log(2) + logMeanExpLogs(c(log.pchangedIT[j+1],
                                                        
                                                        logMeanExpLogs(rowSums(cbind(
                                                          logp1IT*matrix(rev(as.numeric(intToBits(i)[1:N])),ncol=N,nrow=dim(logp1IT)[1], byrow=TRUE), 
                                                          logp0IT*(1- matrix(rev(as.numeric(intToBits(i)[1:N])),ncol=N,nrow=dim(logp0IT)[1],byrow=TRUE))
                                                        )))
        ))
        
        log.prior.pchangedInt[j+1] = log(2) + logMeanExpLogs(c(log.prior.pchangedInt[j+1],
                                                               
                                                               findlogPriorp(iter, rev(as.numeric(intToBits(i)[1:N])), aB1, bB1, r1, r2)
                                                               
        ))
        
        log.prior.pchangedTrend[j+1] = log(2) + logMeanExpLogs(c(log.prior.pchangedTrend[j+1],
                                                                 
                                                                 findlogPriorp(iter, rev(as.numeric(intToBits(i)[1:N])), aB3, bB3, r3, r4)
                                                                 
        ))
        
        log.prior.pchangedIT[j+1] = log(2) + logMeanExpLogs(c(log.prior.pchangedIT[j+1],
                                                              
                                                              findlogPriorpIT(iter, rev(as.numeric(intToBits(i)[1:N])), aB1, bB1, aB3, bB3, r1, r2, r3, r4)
                                                              
        ))
        
      }
      
    }
  }
  
  close(pb)
  
  # Log odds
  log.post.oddsInt = vector(length=N)
  log.prior.oddsInt = vector(length=N)
  log.post.oddsTrend = vector(length=N)
  log.prior.oddsTrend = vector(length=N)
  log.post.oddsIT = vector(length=N)
  log.prior.oddsIT = vector(length=N)
  
  for(i in 2:(N+1)){
    log.post.oddsInt[i-1] =  (log(length(i:(N+1))) + logMeanExpLogs(log.pchangedInt[i:(N+1)])) - 
      (log(length(1:(i-1))) + logMeanExpLogs(log.pchangedInt[1:(i-1)]))
    
    log.post.oddsTrend[i-1] =  (log(length(i:(N+1))) + logMeanExpLogs(log.pchangedTrend[i:(N+1)])) - 
      (log(length(1:(i-1))) + logMeanExpLogs(log.pchangedTrend[1:(i-1)]))
    
    log.post.oddsIT[i-1] =  (log(length(i:(N+1))) + logMeanExpLogs(log.pchangedIT[i:(N+1)])) - 
      (log(length(1:(i-1))) + logMeanExpLogs(log.pchangedIT[1:(i-1)]))
    
    
    log.prior.oddsInt[i-1] = (log(length(i:(N+1))) + logMeanExpLogs(log.prior.pchangedInt[i:(N+1)])) - 
      (log(length(1:(i-1))) + logMeanExpLogs(log.prior.pchangedInt[1:(i-1)]))
    
    log.prior.oddsTrend[i-1] = (log(length(i:(N+1))) + logMeanExpLogs(log.prior.pchangedTrend[i:(N+1)])) - 
      (log(length(1:(i-1))) + logMeanExpLogs(log.prior.pchangedTrend[1:(i-1)]))
    
    log.prior.oddsIT[i-1] = (log(length(i:(N+1))) + logMeanExpLogs(log.prior.pchangedIT[i:(N+1)])) - 
      (log(length(1:(i-1))) + logMeanExpLogs(log.prior.pchangedIT[1:(i-1)]))
  }
  
  logBFInt = log.post.oddsInt - log.prior.oddsInt
  logBFTrend = log.post.oddsTrend - log.prior.oddsTrend
  logBFIT = log.post.oddsIT - log.prior.oddsIT
  
  BFs = cbind(logBFInt, logBFTrend, logBFIT)
  
  return(BFs)     
  
}



