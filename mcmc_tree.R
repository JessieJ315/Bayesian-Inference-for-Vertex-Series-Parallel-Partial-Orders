#
# tree (vsp) MCMC
#
# the updated mcmc with the binary-decomposition tree (BDT) representation of VSP.

# source('~/Desktop/POwTie/code/source-code/pofun.R')

# probability distributions

loglikP.vsp <- function(tree,Y,no_cores,r,model){
  
  # calculates the loglikelihood under the error-free ('P') observation model
  #
  # parameter
  # ---------
  # tree: tree (vsp). 
  #   The BDT of interest. 
  # Y: list. 
  #   The data list of varied lengths. 
  # no_cores: int.
  #   The number of cores used for parallel computing. NOT used in this model.
  # r: float. 
  #   The error rate. NOT used in this model. 
  # model: string.
  #   The type of queue-jumping model. NOT used in this model. 
  #
  # return
  # ------
  #   float. 
  #   The log-likelihood under the error-free ('P') model. 
  
  if (all(sapply(Y,noconflict.tree,tree=tree))){
    # if conflict: return -Inf
    loglkp <- function(y,tree){
      return(-log(nle.tree(sub.tree(tree,y))))
    }
    loglk=sum(sapply(Y,loglkp,tree))
  } else {loglk = -Inf}
  
  return(loglk)
}

loglikQ.vsp <- function(tree,Y,no_cores,r=0,phi=0,model='lkdup'){
  
  # calculates the loglikelihood under the Queue-Jumping ('Q') observation model, 
  # for the jumping-up (model='lkddown'), jumping-down (model='lkdup') and bi-direction (model='bi-directn') model. 
  #
  # parameter
  # ---------
  # tree: tree (vsp). 
  #   The BDT of interest. 
  # Y: list. 
  #   The data list of varied lengths. 
  # no_cores: int.
  #   The number of cores used for parallel computing.
  # r: float. 
  #   The error rate. The error probability p=1/(1+exp(-r)). 
  # phi: float.
  #   The hyperparameter for q (the probability for jumping-up). ONLY USED IN the bi-directn model.
  #   In Jumping-Up, q=1; Jumping-down, q=0. In the Bi-direction model, q=1/(1+exp(-phi)).
  # model: string.
  #   The type of queue-jumping model, takes value 'lkddown' (jumping-up), 'lkdup' (jumping-down) and 'bi-directn' (bi-directional).
  # 
  # return
  # ------
  #   float. 
  #   The QJ log-likelihood under model.
  
  Q.LEProb.vsp <- function(tree,le,p,q=0) {
    
    # calculates the likelihood on a single linear extension (le).
    # if tree is a vsp and le is one list then calculate the likeligood for the list
    # given the vsp-PO tree. Here p (err prob) and q (prob choose top down at a given le entry insertion)
    # adapts to the suborder so le and actors.in.tree have same content
    #
    # parameter
    # ---------
    # tree: tree (vsp). 
    #   The BDT of interest. 
    # le: vector. 
    #   The single data vector.
    # p: float in [0,1]. 
    #   The error probability p=1/(1+exp(-r)). 
    # q: float in [0,1]. 
    #   The probability of jumping up. ONLY used in bi-directional QJ. 
    #
    # return
    # ------
    #   flaot.
    #   The QJ likelihood of a siingle list le. 
    
    tree = sub.tree(tree,le)
    
    n=length(le)
    if (n==1) return(1)
    
    leaf.nodes=which(sapply(tree,is.leaf))
    actors.in.tree=sapply(leaf.nodes, function(x) tree[[x]]$actor)
    
    ntc=NA
    if (q>0) {
      top.i=leaf.nodes[which(actors.in.tree==le[1])]
      if (length(top.i)!=1) stop('err in QP.LEProb.vsp length(top.i)!=1')
      le.not=le[-1]
      tree.not=delete(tree,top.i)
      prob.top=p/n
      if (is.top.bot(tree,top.i,'top')) {
        ntc=nle.tree(tree); 
        prob.top=prob.top+(1-p)*nle.tree(tree.not)/ntc
      }
      top.fac=q*prob.top*QP.LEProb.vsp(tree.not,le.not,p,q)
    } else {top.fac=0}
    
    if (q<1) {
      bot.i=leaf.nodes[which(actors.in.tree==le[n])]
      if (length(bot.i)!=1) stop('err in QP.LEProb.vsp length(bot.i)!=1')
      le.nob=le[-n]
      tree.nob=delete(tree,bot.i)
      prob.bot=p/n
      if (is.top.bot(tree,bot.i,'bot')) {
        if (is.na(ntc)) {ntc=nle.tree(tree)}
        prob.bot=prob.bot+(1-p)*nle.tree(tree.nob)/ntc
      }
      bot.fac=(1-q)*prob.bot*QP.LEProb.vsp(tree.nob,le.nob,p,q)
    } else {bot.fac=0}
    
    return(top.fac+bot.fac)
  }
  
  p = 1/(1+exp(-r)) # define p
  
  if (model=='lkddown') {llkda = foreach(le=Y,.combine=c,.export=c('nle.tree','sub.tree','find.parents','find.actor','is.root','is.top.bot','is.zombie','QP.LEProb.vsp','delete','is.leaf')) %dopar% Q.LEProb.vsp(tree=tree,le=le,p=p,q=1)}
  if (model=='lkdup') {llkda = foreach(le=Y,.combine=c,.export=c('nle.tree','sub.tree','find.parents','find.actor','is.root','is.top.bot','is.zombie','QP.LEProb.vsp','delete','is.leaf')) %dopar% Q.LEProb.vsp(tree=tree,le=le,p=p,q=0)}
  if (model=='bi-directn') {
    q = 1/(1+exp(-phi))
    llkda = foreach(le=Y,.combine=c,.export=c('nle.tree','sub.tree','find.parents','find.actor','is.root','is.top.bot','is.zombie','QP.LEProb.vsp','delete','is.leaf')) %dopar% Q.LEProb.vsp(tree=tree,le=le,p=p,q=q)
  }
  
  return(sum(log(llkda)))
}

logp.vsp <- function(tree,psi){
  
  # calculates the prior probaility for a BDT. 
  #
  # parameter
  # ---------
  # tree: tree (vsp). 
  #   The BDT of interest. 
  # psi: float.
  #   The hyperparameter for q (the probability for S-node) - q=1/(1+exp(-psi)).
  # 
  # return
  # ------
  #   float. 
  #   The prior probaility for BDT (tree). 
  
  q = 1/(1+exp(-psi))
  Pcount=sum(sapply(tree,function(x) x$type=='P'))
  Scount=sum(sapply(tree,function(x) x$type=='S'))
  return(Scount*log(q/2)+Pcount*log(1-q))
}

# tree operations

local_counter <- function(tree,v,l){
  
  # counts the number of potential local moves and returns such list
  #
  # Parameter
  # ---------
  # tree: tree (vsp). 
  #   The BDT of interest. 
  # v: int. 
  #   Node index, the tail of the edge of interest. 
  # l: ind 
  #   Node index, the head of the edge of interest. 
  # 
  # Return
  # ------
  #   list. 
  # 
  #   List Attribute
  #   --------------
  #   el: list. 
  #     The edge list of the potential edges for local move. 
  #   counter: int. 
  #     The number of potential local moves. 
  
  v_up = tree[[v]]$parent
  v_down = setdiff(tree[[v]]$child,l)
  
  edge_ls = matrix(nrow=0,ncol=2)
  colnames(edge_ls) = c('from','to')
  
  # when v is not root
  if (!is.na(v_up)){
    edge_ls = rbind(edge_ls,c(v_up,setdiff(tree[[v_up]]$child,v)))
    if(!is.na(tree[[v_up]]$parent)){edge_ls = rbind(edge_ls,c(tree[[v_up]]$parent,v_up))}
  }
  # when v_down is not leaf
  if (length(tree[[v_down]]$child)>1){
    edge_ls = rbind(edge_ls,matrix(c(rep(v_down,2),tree[[v_down]]$child),ncol=2,byrow=FALSE))
  }
  
  counter = nrow(edge_ls)
  
  return(list(el = edge_ls, counter = counter))
}

reattach_edge <- function(tree,v,l,target_edge){
  
  # reattach edge v -> l to target edge
  # 
  # Parameter
  # ---------  
  # tree: tree (vsp). 
  #   The BDT of interest. 
  # v: int. 
  #   Node index, the tail of the edge of interest.
  # l: ind 
  #   Node index, the head of the edge of interest. 
  # target_edge: array, int. 
  #   The target to attach v->l to, includes 'from' and 'to'.
  # 
  # Return
  # ------
  #   tree (vsp). 
  #   The BDT with edge v -> l reattached to target edge. 
  
  v_up = tree[[v]]$parent
  v_down = setdiff(tree[[v]]$child,l)
  
  tree[[v_down]]$parent = v_up
  tree[[v_down]]$order=tree[[v]]$order
  if(!is.na(v_up)){
    tree[[v_up]]$child = c(v_down,setdiff(tree[[v_up]]$child,v))
    parents_v_old = find.parents(tree,v_up)
    for(parent in parents_v_old){
      tree[[parent]]$nc = tree[[parent]]$nc - tree[[l]]$nc
    }
  } 
  tree[[v]]$parent = target_edge[1]
  tree[[v]]$child = c(l,target_edge[2])
  tree[[target_edge[2]]]$parent = v
  if(!is.na(target_edge[1])){
    tree[[target_edge[1]]]$child = c(setdiff(tree[[target_edge[1]]]$child,target_edge[2]),v)
    parent_v_new = find.parents(tree,target_edge[1])
    for(parent in parent_v_new){
      tree[[parent]]$nc = tree[[parent]]$nc + tree[[l]]$nc
    }
  }
  tree[[v]]$nc = tree[[l]]$nc + tree[[target_edge[2]]]$nc
  
  tree[[v]]$order = tree[[target_edge[2]]]$order
  
  if(tree[[v]]$type=='S'){
    plus_node = sample(c(target_edge[2],l),1)
    tree[[plus_node]]$order = '+'
    tree[[setdiff(c(target_edge[2],l),plus_node)]]$order = '-'
  } else {
    tree[[target_edge[2]]]$order = tree[[l]]$order = ""
  }
  return(tree)
}

edge_ops <- function(tree,l,mode='global'){
  
  # performs edge operation on the vsp tree
  # 
  # parameter
  # ---------
  # tree: tree (vsp). 
  #   The BDT of interest. 
  # l: int.
  #   Node index. The head node to the edge to operate on. 
  # mode: str. 
  #   The type of edge operation, mode = 'local'/'global'.
  #
  # return
  # ------
  #   tree (vsp). 
  #   The updated BDT. 
  
  v = tree[[l]]$parent # id
  
  if(mode=='local'){
    lc = local_counter(tree,v,l)
    if(lc$counter==0){return(tree)}
    local_el = lc$el
    target_edge = local_el[sample(lc$counter,1),]
    tree = reattach_edge(tree,v,l,target_edge)
  }
  
  if(mode=='global'){
    choices = setdiff(1:length(tree),c(v,nodes.under(tree,l)))
    sample.vec <- function(x, ...) x[sample(length(x), ...)]
    target_l = sample.vec(choices,1)
    if(is.root(tree[[target_l]])|length(choices)==1) {
      target_edge = c(NA,target_l)
    } else {
      if(tree[[target_l]]$parent==v){target_edge=c(tree[[v]]$parent,target_l)} else {target_edge = c(tree[[target_l]]$parent,target_l)}
    }
    tree = reattach_edge(tree,v,l,target_edge)
  }
  return(tree)
}

# MCMC

mcmc.vsp <- function(tree,q,p,pq,Y,mode='lkddown',n.itr=1000,n.record=2*n,chain.index,model,write=TRUE){
  
  # performs MCMC for the vsp tree
  #
  # parameter
  # ---------
  # tree: tree (vsp). 
  #   The initial BDT. 
  # q: float, in [0,1].
  #   The initial probability of 'S' node in the vsp tree.
  # p: float, in [0,1].
  #   The inital error probability for the queue-jumping observation model. 
  # pq:float, in [0,1].
  #   The inital probability for jump-up over jump-down in the queue-jumping observation model (ONLY used in bi-directional QJ). 
  # Y: list. 
  #   The data list (can be of varied lengths). 
  # mode: string.
  #   The jumping type for the queue-jumping 'q' model. It takes value 'lkdup' - jumping down, 'lkddown' -jumping up (default), or 'bi-directn' - bi-directional queue-jumping.
  # n.itr: int. 
  #   The number of iteration for the MCMC. 
  # n.record: int. 
  #   The step size to record the results. We only record every n.record steps in the MCMC. 
  # chain.index: string. 
  #   The name of the run, used as the file name to store results. 
  # model: string. 
  #   The observation model considered, can be 'p' (the error-free observation model), 'q' (the queue-jumping observation model) and 'n' (the null observation model, used to check the vsp prior). 
  # write: bool.
  #   Whether to write the results in the current directory. Defult=TRUE.
  #
  # return
  # ------
  #   list. 
  #   The Result lists from MCMC. 
  #   attribute
  #   ---------
  #   tree: list. 
  #     The tree MCMC. 
  #   q: float, in [0,1].
  #     The MCMC of the probability of 'S' node in the vsp tree.
  #   p: float, in [0,1].
  #     The MCMC of the error probability for the queue-jumping observation model. 
  #   loglik: float.
  #     The MCMC of the log-likelihood. 
  
  if(model=='p'){loglik.vsp = loglikP.vsp}
  if(model=='q'){loglik.vsp = loglikQ.vsp}
  if(model=='n'){loglik.vsp = function(tree,Y,no_cores,r,phi,model){return(0)}} # for testing - returns prior
  
  ptm = proc.time()
  N = length(Y)
  n = (length(tree)+1)/2
  
  psi = log(q/(1-q))
  r = log(p/(1-p))
  phi = log(pq/(1-pq))
  
  no_cores <- 2 #detectCores()/4
  registerDoParallel(no_cores)
  
  tree_res = vector(mode= "list", length = n.itr/n.record)
  psi_res = vector(mode="numeric", length = n.itr/n.record)
  r_res = vector(mode="numeric", length = n.itr/n.record)
  phi_res = vector(mode="numeric", length = n.itr/n.record)
  loglk_res = vector(mode="numeric", length = (n.itr/n.record))
  
  tree_res[[1]]=tree
  psi_res[1]=psi
  r_res[1] = r
  loglk = loglik.vsp(tree,Y,no_cores,r,phi,model=mode); loglk_res[1] = loglk
  
  for(i in 2:n.itr){
    
    ## update on S&P
    ind = sample(which(sapply(tree,'[[','type')!=""),1)
    
    tree_temp = tree
    if (tree[[ind]]$type == "P"){
      tree_temp[[ind]]$type="S"
      pos_child = sample(2,1)
      tree_temp[[tree_temp[[ind]]$child[pos_child]]]$order = '+'
      tree_temp[[tree_temp[[ind]]$child[setdiff(1:2,pos_child)]]]$order = '-'
      loglk_temp = loglik.vsp(tree_temp,Y,no_cores,r,phi,model=mode)
      log_eta3 = loglk_temp+logp.vsp(tree_temp,psi)-loglk-logp.vsp(tree,psi)+log(2)
    } else {
      tree_temp[[ind]]$type="P"
      tree_temp[[tree_temp[[ind]]$child[1]]]$order = tree_temp[[tree_temp[[ind]]$child[2]]]$order = ''
      loglk_temp = loglik.vsp(tree_temp,Y,no_cores,r,phi,model=mode)
      log_eta3 = loglk_temp+logp.vsp(tree_temp,psi)-loglk-logp.vsp(tree,psi)-log(2)
    }
    
    if(log_eta3 > log(runif(1))){tree = tree_temp; loglk=loglk_temp}
    
    # update on tree
    ## update on the tree structure
    ### global move = 1
    
    l = sample(which(!sapply(tree,is.root)),1)
    tree_temp = edge_ops(tree,l,mode='global')
    
    loglk_temp = loglik.vsp(tree_temp,Y,no_cores,r,phi,model=mode)
    log_eta1 = loglk_temp+logp.vsp(tree_temp,psi)-loglk-logp.vsp(tree,psi)
    if(log_eta1 > log(runif(1))){tree = tree_temp; loglk=loglk_temp}
    
    ## local move = n
    if(n<=20){itr=n} else {itr=floor(n/10)}
    for (j in 1:itr){ # modified to floor(n/10) for large node number
      l = sample(which(!sapply(tree,is.root)),1)
      v = tree[[l]]$parent # id
      tree_temp = edge_ops(tree,l,mode='local')

      loglk_temp = loglik.vsp(tree_temp,Y,no_cores,r,phi,model=mode)
      log_eta5 = loglk_temp+logp.vsp(tree_temp,psi)-loglk-logp.vsp(tree,psi)+log(local_counter(tree,v,l)$counter)-log(local_counter(tree_temp,v,l)$counter)
      if(!is.na(log_eta5) & log_eta5 > log(runif(1))){tree = tree_temp; loglk=loglk_temp}
    }
    
    # update on q (prior psi ~ norm(0,1.5))
    psi_temp = rnorm(1,psi,1)
    log_eta2 = logp.vsp(tree,psi_temp)+dnorm(psi_temp,0,1.5,log=TRUE)-logp.vsp(tree,psi)-dnorm(psi,0,1.5,log=TRUE)
    if(log_eta2 > log(runif(1))){psi = psi_temp}
    
    # update on r if model='q' (prior r~norm(0,1.5))
    if (model=='q'){
      r_temp = rnorm(1,r,1)
      loglk_temp = loglik.vsp(tree,Y,no_cores,r_temp,phi,model=mode)
      log_eta4 = loglk_temp+dnorm(r_temp,0,1.5,log=TRUE)-loglk-dnorm(r,0,1.5,log=TRUE)
      if(log_eta4 > log(runif(1))){r = r_temp; loglk = loglk_temp}
    }
    
    # update on phi if mode='bi-directn' (prior phi~norm(0,1.5))
    if (mode == 'bi-directn'){
      phi_temp = rnorm(1,phi,1)
      loglk_temp = loglik.vsp(tree,Y,no_cores,r,phi_temp,model=mode)
      log_eta6 = loglk_temp+dnorm(phi_temp,0,1.5,log=TRUE)-loglk-dnorm(phi,0,1.5,log=TRUE)
      if(log_eta6 > log(runif(1))){phi = phi_temp; loglk = loglk_temp}
    }
    
    if (i %% n.record == 0){
      tree_res[[i/n.record]] = tree
      psi_res[i/n.record] = psi
      r_res[i/n.record] = r
      phi_res[i/n.record] = phi
      loglk_res[i/n.record] = loglk
    }
    if (i %% (n.itr/10) == 0){
      print(i)
      if (write==TRUE){
        Result = list(tree=tree_res,q=1/(1+exp(-psi_res)),p=1/(1+exp(-r_res)),pq=1/(1+exp(-phi_res)),loglk=loglk_res)
        save(Result, file=paste0(chain.index,".RData"))
      }
    }
  }
  return(list(tree=tree_res,q=1/(1+exp(-psi_res)),p=1/(1+exp(-r_res)),pq=1/(1+exp(-phi_res)),loglk=loglk_res))
}

