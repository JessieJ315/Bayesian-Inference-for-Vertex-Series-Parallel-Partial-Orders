#
# mcmc example 
#
# load sample data Y
load('sample_data.RData')
# load functions & packages
source('vspfun.R')
source('mcmc_tree.R')
library(igraph)
library(Zseq)
library(foreach)
library(doParallel)
library(ggplot2)
library(gridExtra)
library(reshape2)

# 
n = length(unique(unlist(Y))); n # number of actors
N = length(Y); N # number of lists in Y

# initial values
actors = 1:n
tree = rvsp.tree(actors,prob.snode = 0) # an empty vsp

K = 5000 # number of iterations
Result = mcmc.vsp(tree=tree,q=.5,p=.5,pq=.5,Y,mode='lkddown',n.itr=K,n.record=10,chain.index='test',model='q',write=TRUE)

# result plots
plot.vsp(Result) 
