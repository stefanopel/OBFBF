MarginalRegr = function(Y,X,Gamma,G){
  
  # INPUTS:
  # Y     = nxq matrix of observations
  # X     = mx(p+1) full matrix of predictors
  # Gamma = (p+1) 0/1 vector of chosen predictors
  # G     = qxq adj matrix of underlying DAG
  # OUTPUTS:
  
  ## preliminaries
  q     = dim(Y)[2]            # num of vertices
  n     = dim(Y)[1]            # sample size
  p     = sum(Gamma)-1           # num of active regressors
  
  # function for (log) multivariate gamma
  logMultiGamma = function(q,a){q*(q-1)/4*log(pi) + sum(lgamma((a+1-(1:q))/2))}

  # function for marginal of set of vertices
  MarginalSet = function(set){
    if(length(set)>0){
      set = as.numeric(set)
      YJ     = Y[,set]   # observations in the clique
      B.hat  = solve(t(X.gamma)%*%X.gamma,t(X.gamma)%*%t(t(YJ))) # OLS of regressors parameters
      E      = YJ-X.gamma%*%B.hat  # estimated residuals
      J.card = length(set) # clique cardinality
      
      ## compute multivariate gamma functions
      gamma.num = logMultiGamma(J.card,alpha+n-p-1-(q-J.card))
      gamma.den = logMultiGamma(J.card,alpha+n0-p-1-(q-J.card))
      
      ## compute m
      m = -(n-n0)*J.card/2*log(pi) + gamma.num - gamma.den + 
        J.card*(alpha+n0-(q-J.card))/2*log(n0/n) - (n-n0)/2*log(det(t(E)%*%E))      
    } else m=0
    return(m)
  }
  
  ## fix fractional par to suggested
  alpha = q - 1   # fix par to suggested
  n0    = p + 2   # fix par to suggested
  
  ## convert G in matrix form
  G.mat = matrix(0,q,q)
  G.mat[lower.tri(G.mat)]=G
  G.mat = G.mat + t(G.mat)
  #G.mat = G
  
  ## find cliques and separators
  g         = mpd(G.mat) #graph.adjacency(G.mat,mode="undirected") # convert G to graph
  clique    = g$cliques #maximal.cliques(g)
  n.clique  = length(clique)
  separator = g$separators #minimal.st.separators(g)
  n.sep     = length(separator)
  if(n.sep==0) separator[[1]] = integer(0)
  
  X.gamma = X[,Gamma==1]     # cancel from X regressors not selected

  ## find marginals of cliques and separators
  M.clique = sapply(clique,MarginalSet)
  M.sep    = sapply(separator,MarginalSet)
  
  # total marginal likelihood
  M = sum(M.clique) - sum(M.sep)
  
  return(list(M=M,M.clique=M.clique,M.sep=M.sep))
}
