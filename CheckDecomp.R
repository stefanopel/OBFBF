CheckDecomp = function(chosen,G){
  
  # OBJ: check if edge-perturbed graph is decomposable
  #      From conditions (c) and (d) section 2.1 Green Thomas
  # "Sampling decomposable graphs using a Markov chain on junction trees"
  
  # Inputs
  # chosen  = scalar index of chosen edge (out of q(q-1) choices) 
  # G       = q(q-1) vector representing the starting graph
  
  # output
  # out  = logical 1/0 for allowable or not step
  # case = character a - d cases in Figure 3 Green Thomas
  # C   = cliques of the current graph
  # S   = separators of the current graph
  
  ###################
  ## preliminaries ##
  ###################
  
  # output initialized to FALSE and NULL
  out = FALSE
  case = XS = YS = NULL
  
  q = (1+sqrt(1+8*length(G)))/2  # num of nodes (from q(q-1)/2=length(G))
  
  # define function to convert present G in matrix form
  convert = function(G){G.mat = matrix(0,q,q)
  G.mat[lower.tri(G.mat)] = G
  G.mat = G.mat + t(G.mat)}
  
  # matrix of present and proposed graphs
  connect = G[chosen]==0
  G.old   = convert(G)
  G.star  = G; G.star[chosen] = connect
  G.new   = convert(G.star)
  
  # find (x,y) of involved edge
  newedge = which(abs(G.new-G.old)==1,arr.ind=TRUE)
  newedge = newedge[newedge[,1]>newedge[,2],]
  
  # find cliques of current and proposed graphs
  g = gRbase::mpd(G.old) #graph.adjacency(G.old,mode="undirected") # convert G to graph
  C = g$clique #maximal.cliques(g)
  g.new = gRbase::mpd(G.new) #graph.adjacency(G.new,mode="undirected") # convert G to graph
  C.new = g.new$clique #maximal.cliques(g.new)
  
  # find separators 
  S     = g$separators #minimal.st.separators(g)
  S.new = g.new$separators #minimal.st.separators(g.new)
  
  # cliques indicators in current graph
  x.index   = sapply(C,function(x)newedge[1]%in%x) # cliques containing x
  y.index   = sapply(C,function(x)newedge[2]%in%x) # cliques containing y
  xy.index  = x.index & y.index                    # cliques containing both 
  xy.grid   = expand.grid(which(x.index),which(y.index))
  xy.intersect = apply(xy.grid,1,function(x)intersect(C[[x[1]]],C[[x[2]]]))
  if(!is.null(dim(xy.intersect))) xy.intersect = list(xy.intersect)
  adj.index = sapply(xy.intersect,function(x)length(x)>0)            # adjacent cliques
  if(length(adj.index)==0) adj.index = FALSE
  
  # cliques indicators in proposed graph
  x.new = sapply(C.new,function(x)newedge[1]%in%x) # cliques containing x
  y.new = sapply(C.new,function(x)newedge[2]%in%x) # cliques containing y
  xy.new  = x.new & y.new                    # cliques containing both
  xy.gridnew   = expand.grid(which(x.new),which(y.new))
  xy.internew = apply(xy.gridnew,1,function(x)intersect(C.new[[x[1]]],C.new[[x[2]]]))
  if(!is.null(dim(xy.internew))) xy.internew = list(xy.internew)
  adj.new = sapply(xy.internew,function(x)length(x)>0)            # adjacent cliques
  if(length(adj.new)==0) adj.new = FALSE
  
  if(!connect){  ## case of edge cancellation
    
    ## check condition D) Green Thomas
    if(sum(xy.index)==1){
      out = TRUE
      
      ## indentify case a) b) c) d) in Figure 3 Green Thomas
      if(length(xy.intersect)==1){case = "a"
      } else if(sum(sapply(xy.intersect,function(x)newedge[2]%in%x))<2){case = "b"
      } else if(sum(sapply(xy.intersect,function(x)newedge[1]%in%x))<2){case = "c"
      } else case = "d"
    }                   
  } else{        ## case of edge addition
    
    ## check condition C) Green Thomas
    if(sum(xy.new)==1){
      out = TRUE
      
      ## indentify case a) b) c) d) in Figure 3 Green Thomas
      if(length(xy.internew)==1){case = "a"
      } else if(sum(sapply(xy.internew,function(x)newedge[2]%in%x))<2){case = "b"
      } else if(sum(sapply(xy.internew,function(x)newedge[1]%in%x))<2){case = "c"
      } else case = "d"
    }

return(list(
    out     = out,
    case    = case,
    connect = connect,
    XS      = XS,
    YS      = YS,
    qmax    = max(sapply(C.new,length))
  ))  
}
