pruningwise.branching.times <- function(phy) {   	
  ## calculates branching times = node ages, for ultrametric tree in pruningwise order
  ## !warning! no test that tree is in pruningwise order and ultrametric.
  ## branching time = age = time from the node to tips
  xx <- numeric(phy$Nnode) # times from root to nodes
  nt = length(phy$tip.label)
  interns <- which(phy$edge[, 2] > nt)
  for (i in rev(interns)) { # assumes preorder = 'intern'al nodes in reverse
    ## time of descendant = time of ancestor + branch length:
    xx[phy$edge[i,2] - nt] <- xx[phy$edge[i,1] - nt] + phy$edge.length[i]
  }
  depth <- xx[phy$edge[1, 1] - nt] + phy$edge.length[1] # total tree height
  xx <- depth - xx
  names(xx) <- if (is.null(phy$node.label)) (nt + 1):(nt + phy$Nnode) else phy$node.label
  return(xx)
}

pruningwise.distFromRoot <- function(phy) {
  ## distance from root to all nodes, for tree in pruningwise order
  ## !warning! no test that tree is in pruningwise order.
  nt = length(phy$tip.label)
  xx <- numeric(phy$Nnode + nt)
  for (i in length(phy$edge.length):1)
    xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i]
  names(xx) <- if (is.null(phy$node.label)) 1:(nt + phy$Nnode) else
                 c(phy$tip.label, phy$node.label)
  return(xx)
}



################################################

transf.branch.lengths <-
  function(phy,
           model = c("BM","OUrandomRoot","OUfixedRoot","lambda","kappa","delta","EB","trend","GC"),
           parameters = NULL, check.pruningwise = TRUE, check.ultrametric=TRUE, D=NULL, check.names = TRUE)
{
  if (!inherits(phy, "phylo")) stop("object \"phy\" is not of class \"phylo\".")
  model = match.arg(model)	
  if (model=="trend")
    if (is.ultrametric(phy))
      stop("the trend is unidentifiable for ultrametric trees.")	
  if (is.null(phy$edge.length)) stop("the tree has no branch lengths.")
  if (is.null(phy$tip.label)) stop("the tree has no tip labels.")	
  tol = 1e-10	
  if (check.pruningwise) phy = reorder(phy,"pruningwise")
  n <- length(phy$tip.label)
  N <- dim(phy$edge)[1]
  ROOT <- n + 1L
  anc <- phy$edge[, 1]
  des <- phy$edge[, 2]
  externalEdge = (des <= n)
  if (!is.null(phy$root.edge))
    if (phy$root.edge>0)
      stop("the tree is supposed to have no root edge (or of length 0).")

  ## Default parameters
  parameters.default = c(0,1,1,1,0,1,0)
  names(parameters.default) = c("alpha", "lambda", "kappa", "delta", "rate", "lambda_GC", "sigma2_error")

  ## User defined parameters
  if (is.null(parameters)) {
    parameters = parameters.default
  } else { 
    if (!inherits(parameters, "list")) {
      stop("please specify parameters as a list().")
    } else {
      specified <- !c(is.null(parameters$alpha),
                      is.null(parameters$lambda),
                      is.null(parameters$kappa),
                      is.null(parameters$delta),
                      is.null(parameters$rate),
                      is.null(parameters$lambda_GC),
                      is.null(parameters$sigma2_error))
      parameters.user <- c(parameters$alpha,
                           parameters$lambda,
                           parameters$kappa,
                           parameters$delta,
                           parameters$rate,
                           parameters$lambda_GC,
                           parameters$sigma2_error)
      names(parameters.default) = c("alpha", "lambda", "kappa", "delta", "rate", "lambda_GC", "sigma2_error")
      parameters <- parameters.default
      parameters[specified] <- parameters.user 
    }				
  }
  p = list(alpha = parameters[1],
    lambda = parameters[2],
    kappa = parameters[3],
    delta = parameters[4],
    rate = parameters[5],
    lambda_GC = parameters[6],
    sigma2_error = parameters[7])  # note that sigma2_error = true_sigma2_error/sigma2)

  root.edge = 0 # default, holds for most models. Assumes original tree has no root edge.
  diagWeight = rep(1,n)
  errEdge = rep(p$sigma2_error,n)

  ## BM model
  if (model %in% c("BM","trend")) {
    edge.length = phy$edge.length
  }	
  ## OU models
  OU = c("OUrandomRoot","OUfixedRoot")
  if (model %in% OU) {
    if (check.ultrametric){
      D = numeric(n) # adjustments to external branck lengths
      if (!is.ultrametric(phy)){
        flag = 1
        dis = pruningwise.distFromRoot(phy) # has all nodes
        D = max(dis[1:n]) - dis[1:n]
        D = D - mean(D)
        phy$edge.length[externalEdge] <- phy$edge.length[externalEdge] + D[des[externalEdge]]
      }
      ## phy is now ultrametric
    } else {
      if (is.null(D)) stop("Provide D if you choose check.ultrametric=F")
      if (length(D)!=n) stop("D should be a vector with one term for each tip in the tree")
      if (check.names) {
        if (is.null(names(D)))  stop("D is lacking names (tip labels)")
        ordr = match(phy$tip.label, names(D))
        if (sum(is.na(ordr))>0) stop("names of D do not match the tree tip labels.")
        D = D[ordr,drop=F]
      }
    }
    times <- pruningwise.branching.times(phy) # has internal nodes only
    Tmax <- max(times)
    alpha = p$alpha
    errEdge = errEdge*exp(-2*alpha*D[des[externalEdge]]) # adjust measurement errors for OU models
    ## OUrandomRoot model	
    if (model=="OUrandomRoot") {
      distFromRoot <-  exp(-2*alpha*times) # fixit: divide by 2 alpha??
      d1 = distFromRoot[anc-n] # distFromRoot has internal nodes only, not the n external nodes.
      d2 = numeric(N)
      d2[externalEdge]  = exp(-2*alpha*D[des[externalEdge]])
      d2[!externalEdge] = distFromRoot[des[!externalEdge]-n]
    }
    ## OUfixedRoot model
    if (model=="OUfixedRoot") {	
      distFromRoot <-  exp(-2*alpha*times)*(1 - exp(-2*alpha*(Tmax-times))) # fixit: divide by 2 alpha?
      d1 = distFromRoot[anc-n]
      d2 = numeric(N)
      d2[externalEdge] = exp(-2*alpha*D[des[externalEdge]]) * (1-exp(-2*alpha*(Tmax-D[des[externalEdge]])))
      d2[!externalEdge]= distFromRoot[des[!externalEdge]-n]
    }
    edge.length = d2 - d1
    root.edge = min(distFromRoot)
    diagWeight = exp(alpha*D)
  }
  ## lambda model
  if (model=="lambda") {
    lambda = p$lambda
    distFromRoot <- pruningwise.distFromRoot(phy)
    edge.length = phy$edge.length * lambda 
    edge.length[externalEdge] = edge.length[externalEdge] + (1-lambda)*distFromRoot[des[externalEdge]]
  }
  ## kappa model
  if (model=="kappa") {
    kappa = p$kappa
    edge.length = phy$edge.length^kappa
  }
  ## delta model
  if (model=="delta") {
    delta = p$delta
    distFromRoot <- pruningwise.distFromRoot(phy)
    depth = max(distFromRoot)
    edge.length = (distFromRoot[des]^delta - distFromRoot[anc]^delta)*depth^(1-delta)
  }
  ## early burst model
  if (model=="EB") {
    rate = p$rate
    if (rate==0) edge.length = phy$edge.length
    else {
      distFromRoot <- pruningwise.distFromRoot(phy)
      edge.length = (exp(rate*distFromRoot[des])-exp(rate*distFromRoot[anc]))/rate
    }			
  }
  ## GC model
  if (model=="GC") {
    lambda_GC <- p$lambda_GC
    coal_proba <- function(t) return(1 - exp(-t))
    distFromRoot <- pruningwise.distFromRoot(phy)
    edge.length <- phy$edge.length + (lambda_GC - 1) * (coal_proba(distFromRoot[des]) - coal_proba(distFromRoot[anc])) # note: this is zero for a zero length branch
    edge.length[externalEdge] <- edge.length[externalEdge] + 1 + (lambda_GC - 1) * (1 - coal_proba(distFromRoot[des[externalEdge]]))
  }
	
  edge.length[externalEdge] = edge.length[externalEdge] + errEdge # add measurement errors to the tree
  phy$edge.length = edge.length
  phy$root.edge = root.edge
  names(diagWeight) = phy$tip.label
  return(list(tree = phy, diagWeight = diagWeight))
}

################################################

three.point.compute <-
  function(phy, P, Q = NULL, diagWeight = NULL, 
           check.pruningwise = TRUE, check.names = TRUE, check.precision = TRUE) 
{
  ## For an extra tip Variance, like diagonal measurement error E:
  ## V = DVoD + E = D Vn D with Vn = Vo + D^{-1}ED^{-1} still from tree.
  ## D^{-1}ED^{-1} is diagonal --> add these terms to external branch length
  if (!inherits(phy, "phylo")) stop('object "phy" is not of class "phylo".')
  if (check.pruningwise)	phy = reorder(phy,"pruningwise")
  if ((!check.names)&(check.pruningwise))
    stop("check.names has to be TRUE unless check.pruningwise=FALSE")
  n <- length(phy$tip.label)
  N <- dim(phy$edge)[1]
  ROOT <- n + 1L
  root.edge = if (is.null(phy$root.edge)) 0 else phy$root.edge 
  anc <- phy$edge[, 1]
  des <- phy$edge[, 2]
  if (is.null(diagWeight)) {
    diagWeight = rep(1,n)
    names(diagWeight) = phy$tip.label
  } else {
    if (any(abs(diagWeight) <= .Machine$double.eps ^ 0.8) & check.precision)
      stop ("diagonal weights need to be non-zero.")
  }
  flag = 0
  if (is.null(Q)) {
    flag = 1
    Q = rep(1,n)
    names(Q) = phy$tip.label
  }
  P = as.matrix(P)
  Q = as.matrix(Q)
  if (check.names) {
    if (is.null(rownames(P))) stop("P needs to have row names.")
    ordr = match(phy$tip.label, rownames(P))
    if (sum(is.na(ordr))>0)
      stop("row names of P do not match the tree tip labels.")
    P = P[ordr,,drop=F]
    if (is.null(rownames(Q))) stop("Q needs to have row names.")
    ordr = match(phy$tip.label, rownames(Q))
    if (sum(is.na(ordr))>0)
      stop("row names of Q do not match the tree tip labels.")
    Q = Q[ordr,,drop=F]
    if (is.null(names(diagWeight))) stop("diagWeight needs to have names.")
    ordr = match(phy$tip.label, names(diagWeight))
    if (sum(is.na(ordr))>0)
      stop("names of diagWeight do not match the tree tip labels.")
    diagWeight = diagWeight[ordr,drop=F]
  }
  if (nrow(P)!=n)
    stop("the number of rows in P needs to be the same as the number of tips in the tree.")
  if (nrow(Q)!=n)
    stop("the number of rows in Q needs to be the same as the number of tips in the tree.")
  if (length(diagWeight)!=n)
    stop("the length of diagWeight needs to be the same as the number of tips in the tree.")
  ## now doing: Q'V^{-1}P = Q' (DVoD)^{-1} P = (D^{-1}Q)' Vo^{-1} (D^{-1}P)
  ##        and log|V|= log|Vo| + 2log|D|.
  P = cbind(rep(1,n),P)
  Q = cbind(rep(1,n),Q)
  P = P/diagWeight
  Q = Q/diagWeight
  colP = ncol(P)
  colQ = ncol(Q)
  nout = 2 + colP + colP^2 + colQ + colQ^2 + colP*colQ
  tmp=.C("threepoint", as.integer(N),as.integer(n),as.integer(phy$Nnode),
      as.integer(colP), as.integer(colQ), as.integer(ROOT), as.double(root.edge),
      as.double(phy$edge.length),as.integer(des), as.integer(anc),
      as.double(as.vector(P)), as.double(as.vector(Q)),
      result=double(nout))$result # P=y in threepoint.c, and Q=X
  logd=tmp[1] + 2*sum(log(diagWeight)) # vec11, P1 and Q1 not needed here
  PP=matrix(tmp[2+colP+ 1:(colP^2)], colP, colP)
  QQ=matrix(tmp[2+colP+colP^2+colQ+ 1:(colQ^2)], colQ, colQ)
  QP=matrix(tmp[2+colP+colP^2+colQ+colQ^2 + 1:(colP*colQ)], colQ, colP)

  if (flag==1) # Q absent
    return(list(vec11=PP[1,1], P1=PP[1,-1], PP=PP[-1,-1], logd=logd))
  return(list(vec11=PP[1,1], P1=PP[1,-1], PP=PP[-1,-1],
              Q1=QQ[1,-1], QQ=QQ[-1,-1], QP=QP[-1,-1], logd=logd)) 

}

add.replicates <- function(phy, data, species_id = "species", sample_id = "sample", eps = .Machine$double.eps^2) {
  if (!requireNamespace("phytools", quietly = TRUE)) {
    stop("Package 'phytools' is needed for function 'add.replicates'.", call. = FALSE)
  }
  if (is.null(data[[sample_id]])) {
    stop("The `data` data frame should contain a column named ", sample_id, " with sample ids. Plead adjust argument `sample_id` accordingly.")
  }
  if (length(data[[sample_id]]) != length(unique(data[[sample_id]]))){
    stop("The sample ids in column named ", sample_id, " should be unique identifiers of the samples.")
  }
  if (is.null(data[[species_id]])) {
    stop("The `data` data frame should contain a column named ", species_id, " with species names for each sample. Plead adjust argument `species_id` accordingly.")
  }
  data_tree_cor <- match(data[[species_id]], phy$tip.label)
  if (anyNA(data_tree_cor)) {
    # Species in data NOT in the tree
    stop("Species '", paste(unique(data[[species_id]][is.na(data_tree_cor)]), collapse = "', '"), "' are in the data but not in the tree. Please remove them from the data before proceding (or use a tree that include them)." )
  }
  tree_data_cor <- match(phy$tip.label, data[[species_id]])
  if (anyNA(tree_data_cor)) {
    # Species in tree NOT in the data
    warning("Species '", paste(unique(phy$tip.label[is.na(tree_data_cor)]), collapse = "', '"), "' are in the tree but not in the data. They will be droped from the final tree." )
  }
  ## Make tree
  tree_rep <- phy
  ## create unique labels for original tips to avoid confusion
  tmpname <- paste(sample(c(letters, 0:9), 10, replace = TRUE), collapse = "")
  tip_labels_original <- paste(phy$tip.label, tmpname, sep = "_")
  tree_rep$tip.label <- tip_labels_original
  # Add replicates
  for (tip_index in seq_along(phy$tip.label)) {
    all_rep_ids <- data[[sample_id]][data[[species_id]] == phy$tip.label[tip_index]]
    for (rep_id in rev(all_rep_ids)) {
      tree_rep <- phytools::bind.tip(tree_rep, tip.label = rep_id, edge.length = eps,
                                     where = which(tree_rep$tip.label == tip_labels_original[tip_index]))
    }
  }
  # Remove original tips
  tree_rep <- ape::drop.tip(tree_rep, tip_labels_original)
  # Function extracted from ape:::.compressTipLabel
  relabel <- function(y, ref) {
    n <- length(ref)
    label <- y$tip.label
    if (!identical(label, ref)) {
      if (length(label) != length(ref))
        stop("one tree has a different number of tips")
      ilab <- match(label, ref)
      if (any(is.na(ilab)))
        stop("one tree has different tip labels")
      ie <- match(1:n, y$edge[, 2])
      y$edge[ie, 2] <- ilab
    }
    y$tip.label <- ref
    y
  }
  # relabel to initial order
  tree_rep <- relabel(tree_rep, data[[sample_id]])
  # result
  return(tree_rep)
}

