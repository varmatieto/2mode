# Load tnet
library(tnet)

# Load Davis' (1940) Southern Women Dataset
data(tnet)
net <- Davis.Southern.women.2mode
net

str(net)
head(net)

net <- as.tnet(net, type="binary two-mode tnet")
str(net)


# Calculate the reinforcement coefficient (Robins and Alexander, 2004)
reinforcement_tm(net)

# Calculate the global coefficient (Opsahl, 2012)
clustering_tm(net)

# Calculate the local coefficient (Opsahl, 2012)
clustering_local_tm(net)

############################
data(tnet)
net <- Davis.Southern.women.2mode
netx<-net
netx

netx$V2<-netx$V2+100
netx

netx <- as.tnet(netx, type="binary two-mode tnet")
str(netx)


# Calculate the reinforcement coefficient (Robins and Alexander, 2004)
reinforcement_tm(netx)

# Calculate the global coefficient (Opsahl, 2012)
clustering_tm(netx)

# Calculate the local coefficient (Opsahl, 2012)
clu_loc<-clustering_local_tm(netx)
clu_loc

clu_loc[order(clu_loc$lc),  ]

library (igraph)

netg<- graph.data.frame(netx)
netg

##
tkplot(netg)

# set prj type <- TRUE
V(netg)$type <- V(netg)$name %in% net[,1]
netgg<-bipartite.projection(netg)


neta<-netgg$proj2
netb<-netgg$proj1

tkplot(neta)
neta

#############################################################

###############################################################

biotech<-read.table("~/GitHub/2mode/data/biotech.txt",sep=";", 
                    header=T) #, stringsAsFactors = FALSE
biotech

str(biotech)

biotechx <- as.tnet(biotech, type="binary two-mode tnet")
str(biotechx)

pippo<-biotechx$i
biotechx$i<-as.integer(biotechx$i)
biotechx$p<-as.integer(biotechx$p)+100
str(biotechx)


# Calculate the reinforcement coefficient (Robins and Alexander, 2004)
reinforcement_tm(biotechx)

# Calculate the global coefficient (Opsahl, 2012)
clustering_tm(biotech)

# Calculate the local coefficient (Opsahl, 2012)
clu_loc<-clustering_local_tm(biotechx)
clu_loc

clu_loc[order(clu_loc$lc),  ]

zz<-clu_loc$node[!is.na(clu_loc$lc)]


levels(pippo)[zz]
