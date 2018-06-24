# Load the package and dependencies (ape, phytools, corpcor, subplex, spam)
library(mvMORPH)

# Use a specified random number seed for reproducibility
set.seed(14)

tree<-pbtree(n=1000)

# Simulate two selective regimes
state<-as.vector(c(rep("Forest",600),rep("Savannah",400))); names(state)<-tree$tip.label

# Make the tree with mapped states using SIMMAP
tree<-make.simmap(tree, state, model="ER", nsim=1)

# Plot the phylogeny with the mapped discrete trait
col<-c("blue","orange"); names(col)<-c("Forest","Savannah")
plotSimmap(tree,col, fsize=0.6, node.numbers=FALSE, lwd=3, pts=FALSE)


# 2 Random traits evolving along the phylogeny as a two-optimum OU
set.seed(101)
alpha<-matrix(c(1.1,-0.9,-0.9,1),2)
sigma<-matrix(c(0.35,0.06,0.06,0.35),2)
theta<-c(5.5,5.1,1.2,1.4)
data<-mvSIM(tree, param=list(sigma=sigma, alpha=alpha, ntraits=2, mu=theta, names_traits=c("limb.length","limb.width")), model="OUM", nsim=1)
# We now fit a single optimum and a multiple optimum Ornstein-Uhlenbeck process to each trait separately (i.e., univariate analysis), and then a single and a multiple OU process to both traits analyzed simultaneously (i.e., multivariate analysis):
# Fitting the Ornstein Uhlenbeck on the whole tree
trait1_OU1<- mvOU(tree, data[,1], model="OU1", diagnostic=FALSE, echo=FALSE)
trait2_OU1<- mvOU(tree, data[,2], model="OU1", diagnostic=FALSE, echo=FALSE)

# Fitting the Ornstein Uhlenbeck with multiple optimums
trait1_OUM<- mvOU(tree, data[,1], model="OUM", diagnostic=FALSE, echo=FALSE)
trait2_OUM<- mvOU(tree, data[,2], model="OUM", diagnostic=FALSE, echo=FALSE)
# Compare the AIC values between models fit
AIC(trait1_OUM); AIC(trait1_OU1)
AIC(trait2_OUM); AIC(trait2_OU1)

OUM<- mvOU(tree, data, model="OUM")

OU1<- mvOU(tree, data, model="OU1")


AIC(OUM); AIC(OU1)
