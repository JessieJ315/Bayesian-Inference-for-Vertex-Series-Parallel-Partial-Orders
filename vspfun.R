#
# VSP helper functions
#

## tree helper functions

node<-function(actor=NA,parent=NA,child=NA,type='',order='',nc=NA) {
  # creates a new node in a tree (vsp). default - creates an empty node to be filled in
  #
  # parameter
  # ---------
  #   actor: int. 
  #     The actor represented by this new node. default actor = NA. 
  #   parent: int. 
  #     The node index of its parent node. default parent = NA. 
  #   child: int. 
  #     The node index of its child node. default child = NA. 
  #   type: character. 
  #     The type of operation on internal node, type = "S" (Serial) or "P" (Parallel). default type = "". 
  #   order: character. 
  #     The order for serial operation, order = "+" or "-". Order "+" nodes beat the order "-" nodes. default order = "". 
  #   nc: int.
  #     The number of leaves covered by the node. default nc=NA. 
  # 
  # return
  # ------
  #   list.
  #   A new node in a tree (vsp).
  list(actor=actor,parent=parent,child=child,type=type,order=order,nc=nc)
}
is.leaf<-function(node) { 
  # test whether a node is a leaf node.
  # parameter
  # ---------
  # node: list.
  #   A node of a tree.
  # return
  # ------
  #   bool. 
  (!is.zombie(node) && is.na(node$child[1])) 
  } 
is.anst<-function(node) { 
  # test whether a node is an ansester node.
  # parameter
  # ---------
  # node: list.
  #   A node of a tree.
  # return
  # ------
  #   bool. 
  (!is.zombie(node) && !is.na(node$child[1])) 
  }
is.root<-function(node) { 
  # test whether a node is a root node.
  # parameter
  # ---------
  # node: list.
  #   A node of a tree.
  # return
  # ------
  #   bool. 
  (!is.zombie(node) && is.na(node$parent))
  }
is.zombie<-function(node) {
  # test whether a node is a zombie node.
  # parameter
  # ---------
  # node: list.
  #   A node of a tree.
  # return
  # ------
  #   bool. 
  identical(node,list())
  }
is.cherry<-function(h,i,j) {
  # test whether actors i & j are cherries of h (vsp).
  # parameter
  # ---------
  # h: matrix.
  #   The transitive closure of VSP h.
  # i,j: int in [n]. 
  #   The actors of interest. 
  # return
  # ------
  #   bool. 
  a=all(h[i,-c(i,j),drop=FALSE]==h[j,-c(i,j),drop=FALSE])
  b=all(h[-c(i,j),i,drop=FALSE]==h[-c(i,j),j,drop=FALSE])
  return(a&b)
}
nodes.under<-function(tree,v){
  # returns the indices of the nodes under v (including v)
  # 
  # parameter
  # ---------
  # tree: tree (vsp).
  #   The BDT of interest.
  # v: int. 
  #   The node (index v) of interest.
  #
  # return
  # ------
  #   vector. 
  #   The indices of all nodes under node v.
  if(all(is.na(tree[[v]]$child))) return(v)
  nodes = c(v,nodes.under(tree,tree[[v]]$child[1]),nodes.under(tree,tree[[v]]$child[2]))
  return(nodes)
}
find.parents <- function(tree,i){
  # returns the indices of all nodes above node i (including i)
  # 
  # parameter
  # ---------
  # tree: tree (vsp).
  #   The BDT of interest.
  # i: int. 
  #   The node (index i) of interest.
  #
  # return
  # ------
  #   vector. 
  #   The indices of all nodes above node i.
  if (is.na(i)) stop('The index i is NA.')
  parents.i = i
  while(!is.na(tree[[i]]$parent)){
    parents.i = c(parents.i,tree[[i]]$parent)
    i = tree[[i]]$parent
  }
  return(parents.i)
}
find.actor <- function(tree,actor){
  # returns the node index of an actor name 
  #
  # paremter
  # --------
  # tree: tree (vsp).
  #   The BDT of interest. 
  # actor: int.
  #   The actor index for the actor of interest.
  # 
  # return
  # ------
  #   int. 
  #   The node index in tree for the actor in intest.
  return(which(sapply(tree,function(node){if(is.null(node[[1]])){return(NA)} else {return(node$actor)}})==actor))
}
get.child.count<-function(tree,i=NA) {
  # get a 2 x length(tree) matrix (not counting zombie nodes)
  # - top row is the node index in tree[] 
  # - bot row is the number of leaves that node covers.
  # - a leaf covers one leaf
  # parameter
  # ---------
  # tree: tree (vsp). 
  #   The BDT of interest. 
  # i: int.
  #   Node index. Default i=NA (root).
  # return
  # ------
  #   int. 
  #   The total number of child in tree under node i.
  
  if (is.na(i)) i=which(sapply(tree,function(y) is.root(y)))
  
  if (is.anst(tree[[i]])) {
    gcc1<-get.child.count(tree,tree[[i]]$child[1])
    gcc2<-get.child.count(tree,tree[[i]]$child[2])
    gcci<-matrix(c(i,gcc1[2,1]+gcc2[2,1]),2,1)
    gcc<-cbind(gcci,gcc1,gcc2)
  } else {
    gcc=matrix(c(i,1),2,1)
  }
  return(gcc)
}
is.top.bot<-function(tree,i,extreme='top') {
  # test whether node i is an extreme node in tree
  # parameter
  # ---------
  # tree: tree (vsp).
  #   The BDT of interest.
  # i: int. 
  #   Node index. 
  # extreme: string.
  #   Type of extreme node. Takes value 'top' (top-node) or 'bottom' (bottom node).
  # return
  # ------
  #   Bool. 
  if (extreme=='top') {od='-'} else {od='+'} 
  #test if a leaf node of tree is a top or bottom node in the PO 
  if (!is.leaf(tree[[i]])) stop('is.top.bot input node i is not a leaf')
  if (is.root(tree[[i]])) return(TRUE) #tree has just one leaf which is also root
  j=i
  if (tree[[j]]$order==od) return(FALSE)
  while (!is.na(tree[[j]]$parent)) {
    j=tree[[j]]$parent
    if (tree[[j]]$order==od) return(FALSE)    
  }
  return(TRUE)
}
delete<-function(tree,i) {
  
  # deletes leaf node i from tree
  # reconnect tree - notice that the deleted leaf and parent node stay in tree-list but are unconnected zombies
  # 
  # parameter
  # ---------
  # tree: tree (vsp). 
  #   The vsp tree of interest. 
  # i: int. 
  #   The index of a leaf (!) node. 
  #
  # return
  # ------
  #   tree (vsp). 
  #   The vsp tree with leaf node i removed - i and its parent node stays as zombie nodes (list(NULL)). 
  
  pd=tree[[i]]$parent
  pp=tree[[pd]]$parent
  cd.i=which(tree[[pd]]$child!=i)
  cd=tree[[pd]]$child[cd.i]
  
  tree[[cd]]$parent=pp
  tree[[cd]]$order=tree[[pd]]$order
  if (!is.na(pp)) {
    pp.ci=which(tree[[pp]]$child==pd)
    tree[[pp]]$child[pp.ci]=cd
    parents.i = find.parents(tree,pp)
    for(parent in parents.i){
      tree[[parent]]$nc = tree[[parent]]$nc - 1
    }
  }
  tree[[pd]]=tree[[i]]=list(NULL) #make the deleted nodes zombies
  return(tree)
}

sub.tree <- function(tree,y){
  # creates a subtree given data list y
  #
  # parameter
  # ---------
  # tree: tree (vsp). 
  #   The vsp tree of interest. 
  # y: int, vector. 
  #   A data list representing the order relation between actors. 
  #
  # return
  # ------
  #   tree (vsp). 
  #   The subtree with unrelated nodes removed (zombie nodes). 
  n = (length(tree)+1)/2
  unused = setdiff(1:n,y)
  
  for (actor in unused){
    i = find.actor(tree,actor)
    tree = delete(tree,i)
  }
  return(tree)
}
noconflict.tree <- function(tree,y){
  # tests if a data list y violates the vsp tree
  #
  # parameter
  # ---------
  # tree: tree (vsp). 
  #   The vsp tree of interest. 
  # y: int, vector. 
  #   A data list representing the order relation between actors. 
  # 
  # return
  # ------
  #   bool. 
  #   If the data list y violates the vsp tree, TRUE for no violation and FALSE for violation. 
  tree = sub.tree(tree,y)
  i = find.actor(tree,y[1])
  test = is.top.bot(tree,i)
  
  if (!test) return(FALSE)
  
  while(length(y)>1){
    tree = delete(tree,i)
    y = y[-1]
    i = find.actor(tree,y[1])
    test = is.top.bot(tree,i)
    if (!test) return(FALSE)
  }
  return(TRUE)
}

rvsp.tree<-function(actors,prob.snode=0.5) {
  
  # creates a random vsp tree
  #
  # parameter
  # ---------
  # actors: vector, int. 
  #   A vector with integer values. The actors of interest, e.g. actors = 1:N. 
  # prob.snode: float. in [0,1]. 
  #   The probability of getting a type = "S" (Serial) internal node. 
  #
  # return
  # ------
  #   Tree (vsp), list of length 2N-1.
  #
  #   Attribute
  #   ---------
  #   actor: int. 
  #     The actor represented by node i (i in 1:N). If internal node, actor = NA. 
  #   parent: int. 
  #     The node index of parent node. If root, parent = NA. 
  #   child: int. 
  #     The node index of child node. If leaf node, child = NA. 
  #   type: character. 
  #     The type of operation on internal node, type = "S" (Serial) or "P" (Parallel). If leaf node, type = "". 
  #   order: character. 
  #     The order for serial operation, order = "+" or "-". Order "+" nodes beat the order "-" nodes. If root or parented by "P" nodes, order = "". 
  #   nc: int.
  #     The number of leaves covered by the node. A leaf node covers one leaf. 
  
  n=length(actors)
  tree=vector('list',2*n-1)
  tree[[1]]=node(actor=actors[1],parent=3)
  tree[[2]]=node(actor=actors[2],parent=3)
  tree[[3]]=node(child=c(1,2))
  v=3; root=3
  if (n>2) {
    for (k in 3:n) {
      a=sample(1:v,1)
      new.leaf=v+1;
      new.anst=v+2
      tree[[new.leaf]]=node(actor=actors[k],parent=new.anst)
      tree[[new.anst]]=node(child=c(a,new.leaf))
      if (a==root) {
        root=new.anst
      } else {
        ap=tree[[a]]$parent
        ca=which(tree[[ap]]$child==a)
        tree[[ap]]$child[ca]=new.anst
        tree[[new.anst]]$parent=ap
      }
      tree[[a]]$parent=new.anst
      v=v+2
    }
  }
  
  for (i in 1:(2*n-1)) {
    if (is.anst(tree[[i]])) {
      node.SP=sample(c('S','P'),1,replace=FALSE,prob=c(prob.snode,1-prob.snode))
      tree[[i]]$type=node.SP
      if (node.SP=='S') {
        o=sample(tree[[i]]$child);
        tree[[i]]$child=o
        a=o[1]; b=o[2]; tree[[a]]$order="+"; tree[[b]]$order="-"
      }
    }
  }
  # add children counts
  v=which(sapply(tree,function(x) is.root(x)))
  cc=get.child.count(tree,v)
  n=dim(cc)[2]
  for (i in 1:n) tree[[cc[1,i]]]$nc=cc[2,i]
  return(tree)
}

tree2bin<-function(tree,use.node.names=FALSE) {
  # converts a BDT to a binary matrix
  # parameter
  # ---------
  # tree: tree (vsp). 
  #   The BDT of interest. 
  # use.node.names: bool.
  #   Whether to use the supplied node names.
  # return
  # ------
  #   matrix. 
  nv=length(tree)
  h=matrix(0,nv,nv)
  nm=rep(NA,nv)
  for (i in 1:nv) {
    if (is.anst(tree[[i]])) {
      child=tree[[i]]$child
      h[i,child]<-1
      if (use.node.names) {
        nm[i]=i
      } else {
        nm[i]=paste(tree[[i]]$type, tree[[i]]$order, sep='') #substitute(list(a[b]),list(a=tree[[i]]$type,b=tree[[i]]$order))  #
      }
    } 
    if (is.leaf(tree[[i]])) {
      if (use.node.names) {
        nm[i]=i
      } else {
        nm[i]=paste(tree[[i]]$actor, tree[[i]]$order, sep='') #substitute(list(a[b]),list(a=tree[[i]]$actor,b=tree[[i]]$order))
      }
    }
  }
  
  row.names(h)<-colnames(h)<-nm
  
  #when deleting leaves we create zombie nodes - dont include these
  del=which(sapply(tree, function(x) is.zombie(x)))
  if (length(del)>0) h<-h[-del,-del]
  
  return(h)
}
showTREE<-function(tree=NULL,use.node.names=FALSE,...) {
  # plots a BDT
  # parameter
  # ---------
  # tree: tree (vsp).
  #   The BDT of interest' 
  # use.node.names: bool.
  #   Whether to use supplied node names. 
  m=tree2bin(tree,use.node.names)
  g<-graph_from_adjacency_matrix(m,mode ="directed")#as(m,"graphNEL")
  #h<-graph_from_adjacency_matrix(my.transitive.closure(m),mode ="directed") #ed 5-4-22 for mnem
  h<-g #seems like a bug in layout function puts vertices on top of one another
  plot(g,layout=layout_as_tree(h),...);     
}

nle.tree <- function(tree,v=NA){
  
  # counts the number of le's of the VSP (represented by tree)
  #
  # parameter
  # ---------
  # tree: tree (vsp). 
  #   The BDT of interset. 
  # v: int. 
  #   The vertex index. If v=NA, root is used.
  #
  # return
  # ------
  #   int.
  #   The number of linear extensiosn to vsp tree.
  
  if (is.na(v)) {
    v=which(sapply(tree,function(x) is.root(x)))
    # cc=get.child.count(tree,v)
    # n=dim(cc)[2]
    # for (i in 1:n) tree[[cc[1,i]]]$nc=cc[2,i]
  }
  
  if (is.leaf(tree[[v]])) {
    return(1)
  }
  
  c1=tree[[v]]$child[1]
  c2=tree[[v]]$child[2]
  nle1=nle.tree(tree,c1)
  nle2=nle.tree(tree,c2)
  tot=nle1*nle2
  
  if (tree[[v]]$type=='P') {    
    nc1=tree[[c1]]$nc
    nc2=tree[[c2]]$nc
    tot=tot*choose(nc1+nc2,nc1)
  }
  
  return(tot)
}

## conversion functions

po2tree<-function(h,tree=vector('list',0)) {
  
  # generates a BDT from a VSP.
  # 
  # parameter
  # ---------
  # h: matrix. 
  #   The transitive closure of the VSP. 
  # tree: tree (vsp), default tree = vector('list',0).
  #   The working BDT. 
  #
  # return
  # ------
  #   list. 
  #
  #   Attribute
  #   ---------
  #     h: matrix.  
  #       The unstandardised partial order (transitive closure).
  #     tree: tree (vsp). 
  #       The tree (vsp) corresponding to the input partial order. 
  # 
  # example
  # -------
  # N=6; actors=1:N; tree=rvsp(actors,0.5); h=tree2po(tree)
  # par(mfrow=c(1,2)); showDAG(transitive.reduction(h),edge.arrow.size=3/N); showTREE(tree,edge.arrow.size=3/N)

  n=dim(h)[1]
  if (length(tree)==0) {
    rn=rownames(h)
    tree=lapply(1:n,function(x) node(actor=as.numeric(rn[x])))
    rn=1:n
    rownames(h)<-colnames(h)<-rn 
  }
  rn=sapply(rownames(h),as.numeric)

  new.node.i=length(tree)+1
  if (n>2) {
    i=0; finished_1=FALSE
    while (!finished_1) {
      i=i+1
      if (i==n) stop('no cherries in po2tree()')
      j=i; finished_2=FALSE
      while (!finished_2 & j < n) {
        #ij can be a cherry if they have same relns to all other nodes
        j=j+1
        finished_1=finished_2=is.cherry(h,i,j)
      }
    }
  } else {i=1;j=2}

  i.i=rn[i]
  j.i=rn[j]
  tree[[i.i]]$parent=new.node.i
  tree[[j.i]]$parent=new.node.i
  if (h[i,j]==1) {
    tree[[new.node.i]]=node(child=c(i.i,j.i),type='S',order='');
    tree[[i.i]]$order='+'
    tree[[j.i]]$order='-'
  }
  if (h[j,i]==1) {
    tree[[new.node.i]]=node(child=c(j.i,i.i),type='S',order='');
    tree[[j.i]]$order='+'
    tree[[i.i]]$order='-'
  }
  if (h[i,j]==0 & h[j,i]==0 ) {
    tree[[new.node.i]]=node(child=c(j.i,i.i),type='P',order='');
  }

  h<-h[-j,-j,drop=FALSE]; 
  rn[i]<-new.node.i; rn<-rn[-j] 
  rownames(h)<-colnames(h)<-rn

  if (dim(h)[1]>1) {
    out=po2tree(h,tree)
    h=out$h
    tree=out$tree
  }
  
  # add children counts
  v=which(sapply(tree,function(x) is.root(x)))
  cc=get.child.count(tree,v)
  n=dim(cc)[2]
  for (i in 1:n) tree[[cc[1,i]]]$nc=cc[2,i]

  return(list(h=h,tree=tree))
}
tree2po<-function(tree,r=NA) {
  
  # generates the corresponding transitive closure from a vsp tree
  # 
  # parameter
  # ---------
  # tree: tree (vsp). 
  #   The BDT of interest with N actors. 
  # r: int., default r = NA. 
  #   The node index of the root. 
  # 
  # return
  # ------
  #   matrix (N x N). 
  #   The transitive closure of the corresponding partial order. 

  if (is.na(r)) r=which(sapply(tree,function(x) is.root(x)))

  if (is.leaf(tree[[r]])) {
    h=matrix(0,1,1); row.names(h)<-colnames(h)<-tree[[r]]$actor
    return(h)
  }

  c1=tree[[r]]$child[1]
  c2=tree[[r]]$child[2]
  if(tree[[c1]]$order=='+'|tree[[c1]]$order==''){h1=tree2po(tree,c1)} else {h1=tree2po(tree,c2)}
  if(tree[[c2]]$order=='-'|tree[[c2]]$order==''){h2=tree2po(tree,c2)} else {h2=tree2po(tree,c1)}
  
  n1=dim(h1)[1]
  n2=dim(h2)[1]
  is.S=(tree[[r]]$type=='S')+0
  ub=matrix(is.S,n1,n2)
  lb=t(0*ub)
  h.ub=cbind(h1,ub)
  h.lb=cbind(lb,h2)
  h=rbind(h.ub,h.lb); 
  row.names(h)<-colnames(h)<-c(row.names(h1),row.names(h2))

  return(standardise(h))
}

## helper functions to prior or observation model

QP.LEProb.vsp<-function(tree,le,p,q) {
  
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

## result analysis

dagdepth<-function(tr) {
  # compute dag depth from transitive reduction
  # parameter
  # ---------
  # tr: matrix in [n]x[n].
  #   transitive reduction of a po.
  # return
  # ------
  #   int. 
  #   The depth of po.
  if (length(tr)==1) {return(1)}
  n<-dim(tr)[1]
  cs<-apply(tr,2,sum)
  if (sum(cs)==0) {return(1)}
  csi<-(cs==0)
  bs<-apply(tr,1,sum)
  bsi<-(bs==0)
  free<-which(bsi&csi)
  k<-length(free)
  if (k==n) {return(1)}
  if (k>0) { 
    tr<-tr[-free,-free]
    cs<-apply(tr,2,sum)
    csi<-(cs==0)
    bs<-apply(tr,1,sum)
    bsi<-(bs==0)
  }
  tops<-which(csi)
  bots<-which(bsi)
  if (length(bots)>length(tops)) {tops<-bots}
  return(1+dagdepth(tr[-tops,-tops]))
}

plot.vsp= function(Result, burn_in = 0){
  # plots the result of MCMC
  #
  # parameter
  # ---------
  # Result: list
  #   The MCMC result, including the result tree, p, q, pq and loglk MCMC.
  # burn_in: int.
  #   The burn-in period. 
  #
  # return
  # ------
  #   plots.
  #   result related plots, including traceplots to each parameter, the most probable partial order and the concensus order.
  N = sum(Result$loglk<0); N
  # traceplots  
  # pdf('result_plot.pdf')
  par(mfrow=c(4,1),mar = c(2,5,2,2))
  plot(Result$loglk[1:N],type='l',main='Log-likelihood',ylab='log-likelihood',xlab='',cex.lab=1.5,cex.main=1.5,family='Times')
  abline(v = burn_in, col='red')
  plot(Result$q[1:N],type='l',main='Probability to Series Operation P(S)',ylab='P(S)',cex.lab=1.5,cex.main=1.5,xlab='',family='Times')
  abline(v = burn_in, col='red')
  plot(Result$p[1:N],type='l',main='Error Probability p',ylab='p',xlab='',cex.lab=1.5,cex.main=1.5,family='Times')
  abline(v = burn_in, col='red')
  plot(Result$pq[1:N],type='l',main=expression(paste('Jumping-Up Probability ',phi)),ylab=expression(phi),xlab='',cex.lab=1.5,cex.main=1.5,family='Times')
  abline(v = burn_in, col='red')
  # posteriors
  q_dis = ggplot(melt(data.frame(prior=1/(1+exp(-rnorm(N-burn_in,sd=1.5))),posterior = Result$q[(burn_in+1):N])),aes(x=value,y=..density..,group=variable,color=variable))+
    geom_density(alpha=0.4,adjust = 3)+
    labs(title="Prior/Posterior Distributiosn for P(S)",y="", x = "P(S)")+
    theme_classic()+
    theme(legend.position="none",plot.title = element_text(hjust = 0.5),legend.title=element_blank(),text=element_text(family='Times'))
  p_dis = ggplot(melt(data.frame(prior=1/(1+exp(-rnorm(N-burn_in,sd=1.5))),posterior = Result$p[(burn_in+1):N])),aes(x=value,y=..density..,group=variable,color=variable))+
    geom_density(alpha=0.4,adjust = 3)+
    labs(title="Prior/Posterior distributiosn for Error Probability p",y="", x = "p")+
    theme_classic()+
    theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),legend.title=element_blank(),text=element_text(family='Times'))
  pq_dis = ggplot(melt(data.frame(prior=1/(1+exp(-rnorm(N-burn_in,sd=1.5))),posterior = Result$pq[(burn_in+1):N])),aes(x=value,y=..density..,group=variable,color=variable))+
    geom_density(alpha=0.4,adjust = 3)+
    labs(title="Prior/Posterior distributiosn for Jump-Up Probability pq",y="", x = "pq")+
    theme_classic()+
    theme(legend.position="bottom",plot.title = element_text(hjust = 0.5),legend.title=element_blank(),text=element_text(family='Times'))
  grid.arrange(q_dis,p_dis,pq_dis)
  po_posterior = lapply(Result$tree[(burn_in):N],tree2po)
  n = nrow(po_posterior[[1]])
  po_posterior_tr = lapply(po_posterior,transitive.reduction)
  # depth distribution
  depth_posterior = sapply(po_posterior_tr,dagdepth)
  ggplot(melt(data.frame(posterior=depth_posterior)),aes(x=value,y=..density..,group=variable,color=variable))+
    geom_density(alpha=0.4,adjust = 2.5)+
    labs(title="Posterior Depth Distribution",y="", x = "Depth")+
    theme_classic()+
    theme(legend.position="none",plot.title = element_text(hjust = 0.5),legend.title=element_blank(),text=element_text(family='Times'))
  # most probable partial orders
  str2po=function(x,n){
    matrix(as.numeric(strsplit(x," ")[[1]]),nrow=n)
  }
  PO_unlist=table(sapply(po_posterior,paste,collapse = " "))
  PO_mp=str2po(names(which(PO_unlist==max(PO_unlist))),n)
  par(mfrow=c(1,3))
  showDAG(transitive.reduction(PO_mp))
  showDAG(transitive.reduction(str2po(names(sort(PO_unlist,decreasing = TRUE)[2]),n)))  
  showDAG(transitive.reduction(str2po(names(sort(PO_unlist,decreasing = TRUE)[3]),n)))
  mtext("Top 3 VSPs",side=3,line=-2,outer=TRUE,cex=1.5)
  # concensus po
  concensus_po = Reduce('+',po_posterior)/length(po_posterior) # po_posterior transitive closure
  showDAGcon <- function(con_po,threshold=.5,label_threshold=.9,size = 10){
    con_po1 = (con_po>threshold)+0.0
    con_po2 = (con_po>label_threshold)+0.0
    con_po1 = transitive.reduction(con_po1)
    con_po2 = transitive.reduction(con_po2)
    g = graph_from_adjacency_matrix(con_po1,mode ="directed")#as(m,"graphNEL")
    g_label = graph_from_adjacency_matrix(con_po2)
    el1 <- apply(get.edgelist(g), 1, paste, collapse="-")
    el2 <- apply(get.edgelist(g_label), 1, paste, collapse="-")
    E(g)$color <- ifelse(el1 %in% el2, "red", "darkgrey")
    g = set_edge_attr(g,name='weight',value=-1)
    dist = (-distances(g,v=(V(g)),mode='out'))
    layers = apply(dist,1,max)
    # max per row
    layers = max(layers) - layers
    plot(g,layout=layout_with_sugiyama(g)$layout,vertex.size=size,edge.arrow.size=.2,main='concensus order')
  }
  par(mfrow=c(1,1))
  showDAGcon(concensus_po,threshold = .5)
}

## others

standardise<-function(PO) {
  
  # standardise a partial order adjacency matrix
  #
  # parameter
  # ---------
  # PO: acjacency matrix. 
  #   The partial order adjacency matrix. 
  # 
  # return
  # ------
  #   matrix. 
  #   The standardised adjacency matrix with row and columns sorted.
  
  o=order(as.numeric(row.names(PO)))  
  return(PO[o,o])
}

showDAG<-function(m=NULL,...) {
  # plots a partial order from an adjacency matrix
  # parameter
  # ---------
  # m: matrix.
  #   an adjacency matrix
  g<-graph_from_adjacency_matrix(m,mode ="directed")#as(m,"graphNEL")
  #h<-graph_from_adjacency_matrix(transitive.closure(m,mat=TRUE,loops=FALSE),mode ="directed")
  h<-g #seems like a bug in layout function puts vertices on top of one another
  plot(g,layout=layout_with_sugiyama(h)$layout,...);   	 
}

transitive.reduction = function (g){
  # (from 'nem' package) calculates the transitive reduction of a dag. 
  # some version of R has trouble installing the 'nem' package. we provide the key functions here. 
  # parameter
  # ---------
  # g: matrix [n]x[n]. 
  #   an adjacency matrix.
  # return
  # ------
  #   matrix. 
  #   The transitive reduction of g.
  if (!(class(g)[1] %in% c("matrix", "graphNEL"))) 
    stop("Input must be an adjacency matrix or graphNEL object")
  if (class(g)[1] == "graphNEL") {
    g = as(g, "matrix")
  }
  g = transitive.closure(g, mat = TRUE)
  g = g - diag(diag(g))
  type = (g > 1) * 1 - (g < 0) * 1
  for (y in 1:nrow(g)) {
    for (x in 1:nrow(g)) {
      if (g[x, y] != 0) {
        for (j in 1:nrow(g)) {
          if ((g[y, j] != 0) & sign(type[x, j]) * sign(type[x, 
                                                            y]) * sign(type[y, j]) != -1) {
            g[x, j] = 0
          }
        }
      }
    }
  }
  g
}

transitive.closure = function (g, mat = FALSE, loops = TRUE) {  
  # (from 'nem' package) calculates the transitive closrue of a dag. 
  # some version of R has trouble installing the 'nem' package. we provide the key functions here. 
  # parameter
  # ---------
  # g: matrix [n]x[n]. 
  #   an adjacency matrix.
  # return
  # ------
  #   matrix. 
  #   The transitive closure of g.
  if (!(class(g)[1] %in% c("graphNEL", "matrix")))         
    stop("Input must be either graphNEL object or adjacency matrix")    
  g <- as(g, "matrix")    
  n <- ncol(g)    
  matExpIterativ <- function(x, pow, y = x, z = x, i = 1) {        
    while (i < pow) {            
      z <- z %*% x            
      y <- y + z            
      i <- i + 1}        
    return(y)}    
  h <- matExpIterativ(g, n)    
  h <- (h > 0) * 1    
  dimnames(h) <- dimnames(g)    
  if (!loops)         
    diag(h) <- rep(0, n)    
  else diag(h) <- rep(1, n)    
  if (!mat)         
    h <- as(h, "graphNEL")    
  return(h)
  
  if (!(class(g)[1] %in% c("matrix", "graphNEL")))         
    stop("Input must be an adjacency matrix or graphNEL object")    
  if (class(g)[1] == "graphNEL") {g = as(g, "matrix")}    
  g = transitive.closure(g, mat = TRUE)    
  g = g - diag(diag(g))    
  type = (g > 1) * 1 - (g < 0) * 1    
  for (y in 1:nrow(g)) {        
    for (x in 1:nrow(g)) {            
      if (g[x, y] != 0) {                
        for (j in 1:nrow(g)) {                  
          if ((g[y, j] != 0) & sign(type[x, j]) * sign(type[x,y]) * sign(type[y, j]) != -1) {g[x, j] = 0}}}}}    
  g}






