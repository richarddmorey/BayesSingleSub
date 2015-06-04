
bayesMLint = function(y, subject, occ, int, testPatt=rep(0, length(unique(subject))), iter=10000, thin=1, r1=.5, r2=.5, 
                    beta_theta=5, trunc = 1,
                    a = -.0001, b=.0001, thres = 1e-12,
                    sdmetlogsig2z=1.5, sdmetmusig2z=.9, sdmettheta=.3,
                    mu_mu_sig2_z=0, sig2_mu_sig2_z=1, return.chains = TRUE, return.BFnrchanged = TRUE, raftery.check = FALSE){
  
  if(a >= b){
    stop("a should be smaller than b")
  }
  
  if(trunc <=0 | trunc > 1){
    stop("trunc should be larger than 0 and not exceed 1")
  }
  
  if(decimalplaces(iter/thin)>0){
    stop("iter should be a multiple of thin")
  }
  
  if(length(unique(int)) > 2){
    stop("variable int should consist of only two unique values")
  }
  
  if(length(unique(testPatt)) > 2){
    stop("variable testPatt should consist of only two unique values")
  }
  
  if(length(testPatt) != length(unique(subject))){
    stop("Length of testpatt should be equal to the number of subjects")
  }
  
  else{
    
    progress = TRUE
    
    subNumb = subject
    
    subject = as.integer(factor(subject))
    
    N = length(unique(subject))
    
    M1 = rep(0,N)
    M2 = rep(0,N)
    
    for(i in 1:length(y)){
      
      if(int[i]==min(int)){
        M1[subject[i]] = M1[subject[i]] + 1
      }
      if(int[i]==max(int)){
        M2[subject[i]] = M2[subject[i]] + 1
      }
      
    }
    
    M = M1+M2
    cumM = cumsum(M); cumM = c(0,cumM)
    cumMsqr = cumsum(M^2); cumMsqr = c(0,cumMsqr)
    distMat = vector(length = cumMsqr[N+1])
    
    for(i in 1:N){
      distMat[(cumMsqr[i]+1):cumMsqr[i+1]] = toeplitz(occ[(cumM[i]+1):cumM[i+1]] - min(occ[(cumM[i]+1):cumM[i+1]]))       
    }
    
    test = testPatt
    test[which(testPatt == max(testPatt))] = 1
    test[which(testPatt == min(testPatt))] = 0
    
    if(progress){
      progress = round(iter/100)
      pb = txtProgressBar(min = 0, max = as.integer(iter), style = 3)      
    }else{ 
      pb=NULL 
    }
    
    if(progress && return.BFnrchanged){ message("pb 1/2") }
    
    pbFun = function(samps){ if(progress) setTxtProgressBar(pb, samps)}
    
    out = .Call("gibbsSamplerCall",        
                as.numeric(y),
                as.integer(M1),
                as.integer(M2),
                as.integer(N),
                as.integer(test),
                as.integer(iter),
                as.integer(thin),
                as.numeric(r1),
                as.numeric(r2),
                as.numeric(beta_theta),
                as.numeric(trunc),
                as.numeric(a),
                as.numeric(b),
                as.numeric(thres),
                as.numeric(sdmetlogsig2z),
                as.numeric(sdmetmusig2z),
                as.numeric(sdmettheta),
                as.numeric(mu_mu_sig2_z),
                as.numeric(sig2_mu_sig2_z),
                as.numeric(distMat),
                progress, 
                pbFun,
                new.env(),
                package="BayesSingleSub")
    
    if(progress){ close(pb)}
    
    BFsubject = as.matrix(out[[1]])
    dim(BFsubject) = c(1,N)
    colnames(BFsubject) = c(paste("logbf subject",unique(subNumb)))
    
    BFpatt = as.matrix(out[[2]])
    colnames(BFpatt) = c("logbf pattern")
    
    chains = as.matrix(out[[3]])
    dim(chains) = c(iter/thin, 5*N+6)
    colnames(chains) = c(c(paste("mu",unique(subNumb))), c(paste("delta*sig2z",unique(subNumb))), c(paste("sig2z",unique(subNumb))), 
                         "mu_delta*sig2z", "musig2z", "sig2sig2z", "g1", "g2", "rho", c(paste("logp0",unique(subNumb))), c(paste("logp1",unique(subNumb))))
    
    
    # Check convergence
    
    if(raftery.check == FALSE){
      conv = "convergence not checked"
    }
    
    if(raftery.check == TRUE){
      
      raft = raftery.diag(chains)
      
      if(is.null(dim(raft[[2]]))){
        
        conv = 0
        
        message("Warning: More iterations may be needed for all parameter chains to converge")
        print(raft)
      }
      
      else{
        
        if(sum(raft[[2]][,2] > rep(iter/thin,5*N+6)) > 0){
          
          conv = 0
          
          message("Warning: More iterations may be needed for all parameter chains to converge")
          print(raft)
        }
        else{
          conv = 1
        }
      }
      
    }
    
    if(N > 39 & return.BFnrchanged == TRUE){
      return.BFnrchanged = FALSE
      message("BFnrchanged cannot be computed for N > 39")
    }
    
    if(return.BFnrchanged == TRUE){
      BFnrchanged = -1*as.matrix(func_BF.nr.changed(N,chains[,(3*N+7):(4*N+6)], chains[,(4*N+7):(5*N+6)], a,b, r1, r2, iter=10000))
      dim(BFnrchanged) = c(1,N)
      colnames(BFnrchanged) = c(paste("logbf at least",1:(N-1)),N)
    }
    
    
    if(thin > 1){
      accrates = "Acceptance rates only computed when thinning rate = 1" 
    }
    
    if(thin==1){
      accsig2 = vector(length=N)
      for (i in 1:N){
        accsig2[i] = mean(diff(chains[1:(iter/thin-1),2*N + i])!=0)
      }
      
      accmusig2 =  mean(diff(chains[1:(iter/thin-1),3*N+2])!=0)
      acctheta =  mean(diff(chains[1:(iter/thin-1),3*N+6])!=0)
      
      accrates = list(accmusig2,acctheta,accsig2)
      names(accrates) = c("accmusig2", "acctheta", "accsig2")
      
      print(accrates)
    }
    
    
    if(return.chains == TRUE & return.BFnrchanged == TRUE){
      
      ret = list(BFsubject, BFpatt, BFnrchanged, chains, conv, accrates)
      names(ret) = c("logbfSubject", "logbfPattern", "logbfAtleastX", "chains", "converged", "accRates")
      
      return(ret)
    }
    
    if(return.chains == TRUE & return.BFnrchanged == FALSE){
      
      ret = list(BFsubject, BFpatt, chains, conv, accrates)
      names(ret) = c("logbfSubject", "logbfPattern", "chains", "converged", "accRates")
      
      return(ret)
    }
    
    if(return.chains == FALSE & return.BFnrchanged == TRUE){
      
      ret = list(BFsubject, BFpatt, BFnrchanged, conv, accrates)
      names(ret) = c("logbfSubject", "logbfPattern", "logbfAtleastX", "converged", "accRates")
      
      return(ret)
    }
    
    if(return.chains == FALSE & return.BFnrchanged == FALSE){
      
      ret = list(BFsubject, BFpatt, conv, accrates)
      names(ret) = c("logbfSubject", "logbfPattern", "converged", "accRates")
      
      return(ret)
    }
    
  }
}



bayesMLtrend = function(y, subject, occ, int, intpoint = 0, testPatt=rep(0, length(unique(subject))), iter=80000, n.burnin = 100, thin = 20, 
                              r1=.5, r2=.5, r3=.1, r4=.1, beta_theta=5, trunc = 1,
                              aInt = -.0001, bInt=.0001, aTrend = -.0001, bTrend = .0001,
                              mu_mu_sig2_z=0, sig2_mu_sig2_z=1, return.chains = TRUE, return.BFnrchanged = TRUE, raftery.check = FALSE){
  
  if(aInt >= bInt | aTrend >= bTrend){
    stop("aInt should be smaller than bInt and aTrend should be smaller than bTrend")
  }
  
  if(length(unique(subject)) < 2){
    stop("Number of subjects should be at least 2; see the BayesSingleSub package for analysis of single subject data")
  }
    
  
  if(length(unique(int)) > 2){
    stop("variable int should consist of only two unique values")
  }
  
  if(length(unique(testPatt)) > 2){
    stop("variable testPatt should consist of only two unique values")
  }
  
  if(length(testPatt) != length(unique(subject))){
    stop("Length of testpatt should be equal to the number of subjects")
  }
  
  if(length(intpoint) != 1 & length(intpoint) != length(unique(subject)) ){
    stop("Length of intpoint should be equal to the number of subjects")
  }
  
  else{
    
    subNumb = subject
    
    subject = as.integer(factor(subject))
    
    N = length(unique(subject))
    
    M1 = rep(0,N)
    M2 = rep(0,N)
    
    for(i in 1:length(y)){
      
      if(int[i]==min(int)){
        M1[subject[i]] = M1[subject[i]] + 1
      }
      if(int[i]==max(int)){
        M2[subject[i]] = M2[subject[i]] + 1
      }
      
    }
    
    M = M1+M2
    Mmax = max(M)
    cumM = cumsum(M); cumM = c(0,cumM)
    cumMsqr = cumM^2
    
    t = vector(length = cumM[N+1])
    time = vector(length = cumM[N+1])
    txtime = vector(length = cumM[N+1])
    ones = vector(length = cumM[N+1])
    
     
    if(length(intpoint) == 1){
    
    for(i in 1:N){
      
      
      for(j in 1:M[i]){
        t[cumM[i]+j] = -.5
        time[cumM[i]+j] = occ[cumM[i]+j] - .5*(occ[cumM[i] + M1[i]] + occ[cumM[i] + M1[i] + 1]) 
        txtime[cumM[i]+j] = t[cumM[i]+j]*time[cumM[i]+j]
        ones[cumM[i]+j] = 1
      }
      for(j in  (M1[i]+1):M[i]){
        t[cumM[i]+j] = .5
        time[cumM[i]+j] = occ[cumM[i]+j] - .5*(occ[cumM[i] + M1[i]] + occ[cumM[i] + M1[i] + 1]) 
        txtime[cumM[i]+j] = t[cumM[i]+j]*time[cumM[i]+j]
        ones[cumM[i]+j] = 1
        }
      
      }
        
    }
    
    else{
      
        for(i in 1:N){
          
          for(j in 1:M[i]){
            t[cumM[i]+j] = -.5
            time[cumM[i]+j] = occ[cumM[i]+j] - intpoint[i] 
            txtime[cumM[i]+j] = t[cumM[i]+j]*time[cumM[i]+j]
            ones[cumM[i]+j] = 1
          }
          for(j in  (M1[i]+1):M[i]){
            t[cumM[i]+j] = .5
            time[cumM[i]+j] = occ[cumM[i]+j] - intpoint[i] 
            txtime[cumM[i]+j] = t[cumM[i]+j]*time[cumM[i]+j]
            ones[cumM[i]+j] = 1
          }
          
        }
        
     }
      
      
    
    
    # Obtain posterior distributions using JAGS
    
    
    # initialize data and prior parameters
    
    D = array(rep(NA,N*Mmax^2), c(Mmax,Mmax,N))
    
    for (i in 1:N){
      D[1:M[i],1:M[i],i] = toeplitz(occ[(cumM[i]+1):cumM[i+1]] - min(occ[(cumM[i]+1):cumM[i+1]]))
    }
    
    myData = list(Y=y, ones=ones, t=t, time = time, txtime = txtime, N=N, M=M, D=D, beta_theta = beta_theta, trunc=trunc, cumM = cumM, 
                  r1=r1, r2=r2, r3=r3, r4=r4, mu_mu_sig2_e = mu_mu_sig2_z, sig2_mu_sig2_e = sig2_mu_sig2_z)
    
    # list of all parameters to trace
    
    if(return.chains == TRUE){
      params = c("mu_B1","mu_B3","mu_sig2_e","sig2_sig2_e","theta","B1","B3","B0","B2","sig2_e","g1","g2","g3","g4",
                 "post.mu_B1","post.sig2_B1", "post.mu_B3","post.sig2_B3")
    }
    
    if(return.chains == FALSE){
      params = c("post.mu_B1","post.sig2_B1", "post.mu_B3","post.sig2_B3", "mu_sig2_e")
    }
    
    
    test = testPatt
    test[which(testPatt == max(testPatt))] = 1
    test[which(testPatt == min(testPatt))] = 0
    
    log_priorpnullIT= log(pmvt(lower=c(aInt, aTrend), upper=c(bInt, bTrend), sigma=diag(1,2)*c(r1+r2, r3+r4))[1])
    log_priorpaltIT = log1mp(log_priorpnullIT)
    log_prioroddsIT = log_priorpnullIT -log_priorpaltIT
    
    
    
    # fit the model
    
    fileLocation = system.file('jags', package = "BayesSingleSub") 
    
    if(return.BFnrchanged){ message("pb 1/2")}
    
    jagsfit <- jags(data=myData, inits=NULL, params, n.iter=iter, n.burnin=n.burnin, n.chains=1, n.thin=thin,
                    working.directory = fileLocation, model.file="JAGS_model3_trend.txt")
    
    
    jagsfit.mcmc = as.mcmc(jagsfit)
    jagsfit.matrix = as.matrix(jagsfit.mcmc)
    jagsfit.frame = as.data.frame(jagsfit.matrix)
     
    # Add chain names
    
    if(return.chains == TRUE){
      start.B0 = match(colnames(jagsfit.frame),"B0[1]")
      start.B1 = match(colnames(jagsfit.frame),"B1[1]")
      start.B2 = match(colnames(jagsfit.frame),"B2[1]")
      start.B3 = match(colnames(jagsfit.frame),"B3[1]")
      start.sig2_e = match(colnames(jagsfit.frame),"sig2_e[1]")
      start.post.mu_B1 = match(colnames(jagsfit.frame),"post.mu_B1[1]")
      start.post.sig2_B1 = match(colnames(jagsfit.frame),"post.sig2_B1[1]")
      start.post.mu_B3 = match(colnames(jagsfit.frame),"post.mu_B3[1]")
      start.post.sig2_B3 = match(colnames(jagsfit.frame),"post.sig2_B3[1]")
      
      
      B0 = as.matrix(jagsfit.frame[,which(start.B0==1):(which(start.B0==1)+N-1)])
      B1 = as.matrix(jagsfit.frame[,which(start.B1==1):(which(start.B1==1)+N-1)])
      B2 = as.matrix(jagsfit.frame[,which(start.B2==1):(which(start.B2==1)+N-1)])
      B3 = as.matrix(jagsfit.frame[,which(start.B3==1):(which(start.B3==1)+N-1)])
      sig2_e = as.matrix(jagsfit.frame[,which(start.sig2_e==1):(which(start.sig2_e==1)+N-1)])
      mu_B1 = jagsfit.frame$mu_B1
      mu_B3 = jagsfit.frame$mu_B3
      mu_sig2_e = jagsfit.frame$mu_sig2_e
      sig2_sig2_e = jagsfit.frame$sig2_sig2_e
      theta = jagsfit.frame$theta
      g1 = jagsfit.frame$g1
      g2 = jagsfit.frame$g2
      g3 = jagsfit.frame$g3
      g4 = jagsfit.frame$g4
      postmu_B1 = as.matrix(jagsfit.frame[,which(start.post.mu_B1==1):(which(start.post.mu_B1==1)+N-1)])
      postsig2_B1 = as.matrix(jagsfit.frame[,which(start.post.sig2_B1==1):(which(start.post.sig2_B1==1)+N-1)])
      postmu_B3 = as.matrix(jagsfit.frame[,which(start.post.mu_B3==1):(which(start.post.mu_B3==1)+N-1)])
      postsig2_B3 = as.matrix(jagsfit.frame[,which(start.post.sig2_B3==1):(which(start.post.sig2_B3==1)+N-1)])
    }
    
    if(return.chains == FALSE){
      start.post.mu_B1 = match(colnames(jagsfit.frame),"post.mu_B1[1]")
      start.post.sig2_B1 = match(colnames(jagsfit.frame),"post.sig2_B1[1]")
      start.post.mu_B3 = match(colnames(jagsfit.frame),"post.mu_B3[1]")
      start.post.sig2_B3 = match(colnames(jagsfit.frame),"post.sig2_B3[1]")
      
      postmu_B1 = as.matrix(jagsfit.frame[,which(start.post.mu_B1==1):(which(start.post.mu_B1==1)+N-1)])
      postsig2_B1 = as.matrix(jagsfit.frame[,which(start.post.sig2_B1==1):(which(start.post.sig2_B1==1)+N-1)])
      postmu_B3 = as.matrix(jagsfit.frame[,which(start.post.mu_B3==1):(which(start.post.mu_B3==1)+N-1)])
      postsig2_B3 = as.matrix(jagsfit.frame[,which(start.post.sig2_B3==1):(which(start.post.sig2_B3==1)+N-1)])
      
      mu_sig2_e = jagsfit.frame$mu_sig2_e
    }
    
  
    out = .Call("BayesFactorsTrendCall",
                as.numeric(aInt),
                as.numeric(bInt),
                as.numeric(aTrend),
                as.numeric(bTrend),
                as.numeric(as.vector(postmu_B1)),
                as.numeric(as.vector(postsig2_B1)),
                as.numeric(as.vector(postmu_B3)),
                as.numeric(as.vector(postsig2_B3)),
                as.numeric(as.vector(mu_sig2_e)),
                as.integer(test),
                as.integer(length(postmu_B1[,1])),
                as.integer(N),
                as.numeric(r1),
                as.numeric(r2),
                as.numeric(r3),
                as.numeric(r4),
                as.numeric(log_prioroddsIT),
                package="BayesSingleSub") 
    
 
    BFsubjectInt = as.matrix(out[[1]])
    dim(BFsubjectInt) = c(1,N)
    colnames(BFsubjectInt) = c(paste("logbfInt subject",unique(subNumb)))
    
    BFsubjectTrend = as.matrix(out[[2]])
    dim(BFsubjectTrend) = c(1,N)
    colnames(BFsubjectTrend) = c(paste("logbfTrend subject",unique(subNumb)))
    
    BFsubjectIT = as.matrix(out[[3]])
    dim(BFsubjectIT) = c(1,N)
    colnames(BFsubjectIT) = c(paste("logbfit subject",unique(subNumb)))
    
    BFpatt = as.matrix(out[[4]])
    dim(BFpatt) = c(1,3)
    colnames(BFpatt) = c("logbf pattern int", "logbf pattern trend", "logbf pattern it")
    
    if(return.chains == TRUE){
      chains = cbind(B0, B1, B2, B3, sig2_e, mu_B1, mu_B3, mu_sig2_e, sig2_sig2_e, g1, g2, g3, g4, theta)
      
      colnames(chains) = c(c(paste("mu",unique(subNumb))), c(paste("delta*sig2z",unique(subNumb))), c(paste("beta0",unique(subNumb))), c(paste("beta1*sig2z",unique(subNumb))), c(paste("sig2z",unique(subNumb))), 
                           "mu_delta*sig2z", "mu_beta1*sig2z", "musig2z", "sig2sig2z", "g1", "g2", "g3", "g4", "rho")
    }
    
    if(return.chains == FALSE){
      chains = cbind(postmu_B1, postsig2_B1, postmu_B3, postsig2_B3, mu_sig2_e)    
    }
    
    
    # Check convergence
    
    if(raftery.check == FALSE){
      conv = "convergence not checked"
    }
    
    if(raftery.check == TRUE){
      
      raft = raftery.diag(chains)
      
      if(is.null(dim(raft[[2]]))){
        
        conv = 0
        
        message("Warning: More iterations need to be drawn for all parameter chains to converge")
        print(raft)
      }
      
      else{
        
        if(sum(raft[[2]][,2] > rep(length(postmu_B1[,1]),dim(chains)[2])) > 0){
          
          conv = 0
          
          message("Warning: More iterations need to be drawn for all parameter chains to converge")
          print(raft)
        }
        else{
          conv = 1
        }
      }
      
    }
    
    
    if(N > 39 & return.BFnrchanged == TRUE){
      return.BFnrchanged = FALSE
      message("BFnrchanged cannot be computed for N > 39")
    }
    
    if(return.BFnrchanged == TRUE){
      
      BFnrchanged = -1*as.matrix(func_BF.nr.changedTrendJAGS(N, length(postmu_B1[,1]), postmu_B1, postsig2_B1, postmu_B3, postsig2_B3, mu_sig2_e, 
                                                          aInt, bInt, aTrend, bTrend, r1, r2, r3, r4, 10000))
      
      
      dim(BFnrchanged) = c(N,3)
      rownames(BFnrchanged) = c(paste("logbf at least",1:(N-1)),N)
      colnames(BFnrchanged) = c("Int", "Trend", "IT")
    }
    
    
    
    if(return.chains == TRUE & return.BFnrchanged == TRUE){
      
      ret = list(BFsubjectInt, BFsubjectTrend, BFsubjectIT, BFpatt, BFnrchanged, chains, conv)
      names(ret) = c("logbfSubjectInt", "logbfSubjectTrend", "logbfSubjectIT", "logbfPattern", "logbfAtleastX", "chains", "converged")
      
      return(ret)
    }
    
    if(return.chains == TRUE & return.BFnrchanged == FALSE){
      
      ret = list(BFsubjectInt, BFsubjectTrend, BFsubjectIT, BFpatt, chains, conv)
      names(ret) = c("logbfSubjectInt", "logbfSubjectTrend", "logbfSubjectIT", "logbfPattern", "chains", "converged")
      
      return(ret)
    }
    
    if(return.chains == FALSE & return.BFnrchanged == TRUE){
      
      ret = list(BFsubjectInt, BFsubjectTrend, BFsubjectIT, BFpatt, BFnrchanged, conv)
      names(ret) = c("logbfSubjectInt", "logbfSubjectTrend", "logbfSubjectIT", "logbfPattern", "logbfAtleastX", "converged")
      
      return(ret)
    }
    
    if(return.chains == FALSE & return.BFnrchanged == FALSE){
      
      ret = list(BFsubjectInt, BFsubjectTrend, BFsubjectIT, BFpatt, conv)
      names(ret) = c("logbfSubjectInt", "logbfSubjectTrend", "logbfSubjectIT", "logbfPattern", "converged")
            
      return(ret)
    }
    
  }
}




