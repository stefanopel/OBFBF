"# OBFBF" 
    out = OBFBF(Y,X,M)
    
    ## posterior means
    (post.g = colMeans(out$chain.gamma[(M/2):M,]))
    (post.G = colMeans(out$chain.G[(M/2):M,]))
    
    # computational time
    (t1 = out$t1)
    
