OBFBF = function(Y,X,M,a.gamma=0.5,a.G=0.5){
  
  # INPUTS #
  # Y       = n by q matrix of data
  # X       = n by (p+1) matrix of covariates
  # M       = MCMC iterations
  # a.gamma = scalar in proposal prob 0->1 for gamma
  # a.G     = scalar in proposal prob 0->1 for G
  
  # OUTPUT #
  # chain.gamma = M by (p+1) chain for gamma
  # chain.G     = M by q(q-1)/2 chain for graph
  # t1          = computational time
  # accept      = M by 2 chain of accept/reject for 
  #               gamma (1st col) and G (2nd col)
      
      # preliminaries
      p = ncol(X)-1  # num of covariates (excluding covariates)
      q = ncol(Y)  # number of nodes
      n = nrow(Y)  # sample size
      
      ################################
      ## Install packages if needed ##
      ################################  
      
      library("igraph") # package for graphs              
      library("covreg") # package for matrix normal
      library("gRbase") # package for graphs
      
      ###########################################
      ## load function for marginal likelihood ##
      ###########################################
      #  ------------------------------------------------------------------------
      source("MarginalRegr.R")     
      source("CheckDecomp.R")

      #  ------------------------------------------------------------------------
      
      time.curr = proc.time()
      #################################
      ## memory allocations for MCMC ##
      #################################
      #  ------------------------------------------------------------------------
      chain.gamma = matrix(0,M,p+1)        # for gamma (variable selection)
      chain.G     = matrix(0,M,q*(q-1)/2)  # for G (covariance selection)
      accept      = matrix(FALSE,M,2)      # for acceptance/rejection
      #  ------------------------------------------------------------------------
      
      ########################
      ## MCMC preliminaries ##
      ########################
      #  ------------------------------------------------------------------------
      q.edge          = q*(q-1)/2                        # total number of edges
      vec             = function(x){x[lower.tri(x)]}     # define function for vec matrix
      r.gamma         = rbinom(M,1,a.gamma) 
      runif.gamma     = log(runif(M))                    # random variates for gamma steps
      r.G             = rbinom(M,1,a.G)                  
      runif.G         = log(runif(M))                    # random variates for G steps
      chain.gamma[1,] = 0; chain.gamma[1,1] = 1;
      Gamma           = chain.gamma[1,]                  # starting value for gamma
      chain.G[1,]     = 0 
      G               = chain.G[1,]                      # starting value for G
      qmax = 0; check = list(); check$qmax = 0
      #  ------------------------------------------------------------------------
      
      #####################
      ## MCMC iterations ##
      #####################
      for(i in 2:M){
        
        if(i%%100==0){
          print(chain.gamma[i-1,])
          print(chain.G[i-1,])
          print(paste("Iteration",i))
        } 
        ####################
        ## sampling gamma ##
        ####################
        #  ------------------------------------------------------------------------
        
        ## Propose a change in gamma
        p.gamma    = sum(Gamma)   # num of active regressors
        Gamma.star = Gamma; deltap = 0
        if(r.gamma[i]==1){# 0 --> 1
          if(p.gamma<(p+1)){
            chosen = sample.int(p,1,prob=1-Gamma[-1])
            Gamma.star[chosen+1]  = 1
            deltap                = 1  
          }} else if(r.gamma[i]==0){ # 1 --> 0
            if(p.gamma>1){
              chosen = sample.int(p,1,prob=Gamma[-1])
              Gamma.star[chosen+1] = 0
              deltap               = -1
            }}
        
        if(deltap!=0 & n>(p.gamma+deltap+qmax)){
          ## Compute priors of gamma and gamma*
          prior.gamma = -log(p+1)-lchoose(p,p.gamma-1)
          prior.star  = -log(p+1)-lchoose(p,p.gamma-1+deltap)
          
          ## Compute marginals (given G)
          m      = MarginalRegr(Y,X,Gamma,G)$M
          m.star = MarginalRegr(Y,X,Gamma.star,G)$M
          
          r = min(0,m.star - m + prior.star - prior.gamma + 
                    deltap*log((1-a.gamma)/a.gamma))           # acceptance probability
          if(runif.gamma[i]<=r){
            accept[i,1] = TRUE
            Gamma = Gamma.star             # accept/reject gamma
          } 
          p.gamma = p.gamma + deltap
        }
        chain.gamma[i,] = Gamma                              # store gamma
        
        #  ------------------------------------------------------------------------
        
        ################
        ## sampling G ##
        ################
        #  ------------------------------------------------------------------------
        ## Propose a change in G
        p.G    = sum(G)   # num of active links
        G.star = G; deltap = 0
        if(r.G[i]==1){ # 0 --> 1
          if(p.G<q.edge){
            chosen = sample.int(q.edge,1,prob=1-G)
            # print(chosen)
            # print(G)
            check = CheckDecomp(chosen,G)  # check you still decomposable
            if(check$out){
              G.star[chosen]  = 1
              deltap     = 1  
            }
          }} else if(r.G[i]==0){ # 1 --> 0
            if(p.G>0){
              chosen = sample.int(q.edge,1,prob=G)
              # print(chosen)
              # print(G)
              check = CheckDecomp(chosen,G)  # check you still decomposable
              if(check$out){
                G.star[chosen] = 0
                deltap         = -1
              }
            }}
        
        
        if(deltap!=0 & n>(p.gamma+check$qmax)){
          ## Compute priors of G and G*
          prior.G    = -log(q.edge+q)-lchoose(q.edge+q,p.G)
          prior.star = -log(q.edge+q)-lchoose(q.edge+q,p.G+deltap)
          
          ## Compute marginals (given Gamma)
          m      = MarginalRegr(Y,X,Gamma,G)$M
          m.star = MarginalRegr(Y,X,Gamma,G.star)$M    
          
          r = min(0,m.star - m + prior.star - prior.G + 
                    deltap*log((1-a.G)/a.G))             # acceptance probability
          if(runif.G[i]<=r){
            accept[i,2]    = TRUE
            G = G.star                   # accept/reject G
            qmax = check$qmax
          } 
        }
        chain.G[i,] = G                                # store G
        #  ------------------------------------------------------------------------
      }
      
      t1 = proc.time() - time.curr # comp time our method
      
      return(list(chain.gamma=chain.gamma,chain.G=chain.G,t1=t1,accept=accept))
    }
