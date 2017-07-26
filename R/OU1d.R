OU1d.loglik <- function(trait, phy, model=c("OUrandomRoot","OUfixedRoot"), parameters=NULL) {
  
  ## initialize  
  if (!inherits(phy, "phylo")) stop("object \"phy\" is not of class \"phylo\".")
  model = match.arg(model)	
  if (is.null(phy$edge.length)) stop("the tree has no branch lengths.")
  if (is.null(phy$tip.label)) stop("the tree has no tip labels.")	
  tol = 1e-10	
  phy = reorder(phy,"pruningwise")
  n <- length(phy$tip.label)
  N <- dim(phy$edge)[1]
  ROOT <- n + 1L
  anc <- phy$edge[, 1]
  des <- phy$edge[, 2]
  externalEdge <- des<=n
  
  ## Default parameters
  parameters.default = c(0,1,0,1,0)
  names(parameters.default) = c("ancestral.state", "sigma2", 
                                "optimal.value", "alpha", "sigma2_error")
  
  ## User defined parameters
  if (is.null(parameters)) { parameters = parameters.default } else { 
    if (class(parameters)!= "list") {
      stop("please specify parameters as a list().")
    } else {
      specified <- !c(is.null(parameters$ancestral.state), 
                      is.null(parameters$sigma2), 
                      is.null(parameters$optimal.value),
                      is.null(parameters$alpha),
                      is.null(parameters$sigma2_error))
      parameters.user <- c(parameters$ancestral.state,
                           parameters$sigma2,
                           parameters$optimal.value,
                           parameters$alpha,
                           parameters$sigma2_error)
      names(parameters.default) = c("ancestral.state", "sigma2",
                                    "optimal.value", "alpha", "sigma2_error")
      parameters <- parameters.default
      parameters[specified] <- parameters.user 
    }  			
  }
  p = list(ancestral.state = parameters[1],
           sigma2 = parameters[2],
           optimal.value = parameters[3],
           alpha = parameters[4],
           sigma2_error = parameters[5])
  
  ## preparing for general use of "parameter" for branch length transformation
  prm = list(alpha = p$alpha, sigma2_error = p$sigma2_error*2*p$alpha/p$sigma2)
  
  if (is.null(names(trait))) {
    warning("the data has no names, order assumed to be the same as tip labels in the tree.\n")
    data.names = phy$tip.label 
  }
  else data.names = names(trait)
  order = match(data.names, phy$tip.label)
  if (sum(is.na(order))>0) {
    warning("data names do not match with the tip labels.\n")
    names(trait) = data.names
  } else {
    tmp = trait
    names(trait) = phy$tip.label
    trait[order] = tmp[1:n]
  }
  
  dis = pruningwise.distFromRoot(phy)[1:n]
  Tmax = mean(dis)
  
  ### Check if alpha is too big or too small
  if ((prm$alpha < 1e-7/Tmax)||(prm$alpha > 50/Tmax)) return(-Inf)
  
  y = trait
  if (model=="OUrandomRoot") {
    X = as.matrix(rep(p$optimal.value, n))
  } else {
    X = as.matrix(p$ancestral.state*exp(-p$alpha*dis) + p$optimal.value*(1-exp(-p$alpha*dis)))
  } # X is the mean
  
  flag = 0 # flag and D are used for OU model if tree is not ultrametric:
  D = numeric(n)
  if (!is.ultrametric(phy)) {
    flag = 1
    dis = pruningwise.distFromRoot(phy)
    Tmax = max(dis[1:n])
    D = Tmax - dis[1:n]
    D = D - mean(D)
    phy$edge.length[externalEdge] <- phy$edge.length[externalEdge] + D[des[externalEdge]]
    ## phy is now ultrametric, with height Tmax:
    Tmax = Tmax + min(D)
  }
  
  d = ncol(X)
  ole= 4 + 2*d + d*d # output length  
  tree = transf.branch.lengths(phy,model,parameters=prm,
                               check.pruningwise=F,check.ultrametric=F,D=D,check.names=F)$tree
  if (flag) { # need diagonal terms in D
    y = exp(-p$alpha*D) * y # variables local to this function
    X = exp(-p$alpha*D) * X
  }
  tmp=.C("threepoint", as.integer(N),as.integer(n),as.integer(phy$Nnode),
         as.integer(1),as.integer(d),as.integer(ROOT),as.double(tree$root.edge),as.double(tree$edge.length),
         as.integer(des), as.integer(anc), as.double(as.vector(y)), as.double(as.vector(X)),
         result=double(ole))$result # tmp has, in this order:
  ## logdetV, 1'V^{-1}1, y'V^{-1}1, y'V^{-1}y, X'V^{-1}1, X'V^{-1}X, X'V^{-1}y
  comp = list(vec11=tmp[2], y1=tmp[3], yy=tmp[4], X1=tmp[5:(4+d)],
              XX=matrix(tmp[(5+d):(ole-d)], d,d),Xy=tmp[(ole-d+1):ole],logd=tmp[1])
  
  n2llh = as.numeric( n*log(2*pi) + comp$logd + n*log(p$sigma2/2/p$alpha) + 
                        2*p$alpha/p$sigma2*(comp$yy - 2*comp$Xy + comp$XX )) # -2 log-likelihood
  if (n2llh < 0 ) n2llh = Inf # -2 log-likelihood has to be positive
  
  if (flag)
    n2llh = n2llh + p$alpha * 2 * sum(D)
  
  return(unname(-n2llh/2))
}