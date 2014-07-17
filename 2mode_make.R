
biotech<-read.table("biotech.txt",sep=";", 
                  header=T, stringsAsFactors = FALSE)
head(biotech)
str(biotech)
biotech[order(biotech$prj),]
biotech[order(biotech$sbj),]

myprj<-unique(biotech$prj)

# prjsbj represents sbj for each prj

prjsbj<-sapply(myprj, function(x) biotech$sbj[biotech$prj==x])
length(prjsbj)
names(prjsbj)

sbj2<-table(biotech$sbj)[table(biotech$sbj)>1]
names(sbj2)

library (igraph)

gbiotech <- graph.data.frame(biotech)
gbiotech

##
tkplot(gbiotech)
mylay<-tkplot.getcoords(1)


##

# set prj type <- TRUE
V(gbiotech)$type <- V(gbiotech)$name %in% biotech[,1]



################
v_col<-rep("gold",87)

v_col[!V(gbiotech)$type]<-c("red")
v_col[V(gbiotech)$name %in% names(sbj2)]<-c("brown")

plot.igraph(gbiotech, layout=layout.kamada.kawai, edge.arrow.size=0.1,
            vertex.label.cex=.7,vertex.color=v_col,)

######################

gbiotech

nnbio<-bipartite.projection(gbiotech)


snbio<-nnbio$proj2
pnbio<-nnbio$proj1

tkplot(snbio)
snbio

mylay2<-tkplot.getcoords(4)


################
v_col<-rep("gold",78)

v_col[V(snbio)$name %in% names(sbj2)]<-c("brown")

plot.igraph(snbio, layout=mylay2, edge.arrow.size=0.1,
            vertex.label.cex=.7,vertex.color=v_col,)

######################

snbio

##
df_snbio<-get.data.frame(snbio)
df_pnbio<-get.data.frame(pnbio)
str(df_snbio)

write.table(df_snbio,file = "sbj_3bio_el.txt",  sep = " ",
            eol = "\n", na = "NA", dec = ",", row.names = TRUE)
write.table(df_pnbio,file = "prj_3bio_el.txt",  sep = " ",
            eol = "\n", na = "NA", dec = ",", row.names = TRUE)

str(df_snbio)
head(df_snbio)

###


sbj2to<-df_snbio[df_snbio$weight>1,]

unique(append(sbj2to$from,sbj2to$to))

# ----------------------------------------------
snbio 

V(snbio)$name

# clustering coefficient
tt<-transitivity (snbio,type=c("local" )) 


tt[tt<1]
V(snbio)$name[tt<1]
V(snbio)$name[tt<1]%in%names(sbj2)


# betweenness

bb<-betweenness(snbio, normalized = T, weights=NA  )#freeman

V(snbio)$name[bb>0]

cc<-constraint(snbio)
summary(cc)

plot(tt,bb)

# clique.community FUNCTION ##############

clique.community <- function(graph, k) {
  clq <- maximal.cliques(graph) # my correction with maximal
  edges <- c()
  for (i in seq_along(clq)) {
    for (j in seq_along(clq)) {
      if ( length(unique(c(clq[[i]], clq[[j]]))) == k+1 ) {
        edges <- c(edges, c(i,j)-1)
      }
    }
  }
  clq.graph <- simplify(graph(edges))
  V(clq.graph)$name <- seq_len(vcount(clq.graph))
  comps <- decompose.graph(clq.graph)
  
  lapply(comps, function(x) {
    unique(unlist(clq[ V(x)$name ]))
  })
}


cliques(snbio,min=4,max=18)
length(largest.cliques(snbio)[[1]])
maximal.cliques.count(snbio)


# maximal clique  
myclique<-maximal.cliques(snbio)
myclique_list<-sapply(myclique, function(x) V(snbio)$name[x])
length(myclique_list)

##########PLOT ####################
myblocks=length(myclique)
mycolors<-rainbow(myblocks)
v_col<-rep("white",78)


j<-0
for (i in 1:myblocks){
  j<-j+1
  v_col[myclique[[i]] ]<- mycolors[j]
}
v_col[V(snbio)$name %in% names(sbj2)]<-c("white")

plot.igraph(snbio, layout=mylay2, edge.arrow.size=0.1,
            vertex.label.cex=.7,vertex.color=v_col,)






# clique  percolation is a community detection method developed by Gergely Palla
cli_com<-clique.community (snbio,18)

myclique_perco<-sapply(cli_com, function(x) V(snbio)$name[x])
myclique_perco

############################################
c<-matrix(0,10,2)
for (i in names(prjsbj)){
 
  
  for (j in 1:length(myclique_list)) {
    
    if (mean(myclique_list[[j]]%in%prjsbj[[i]])==1){
     c[j,1]<-j 
     c[j,2]<-i 
    }
    
  }
        
}
 c 
  
##################################################
d<-matrix(0,length (prjsbj),length(myclique_perco))

for (i in 1:length(myclique_perco) ) {
  
  
  for (j in 1:length (prjsbj) ) {
    
    if (mean(prjsbj[[j]]  %in% myclique_perco[[i]]) == 1) {
            d[j,i]<-names(prjsbj)[j]
    }
    
  }
  
}

d # da capire??

###############################################
f<-list ()

for (i in names (prjsbj) ) {

  for (j in names (prjsbj) ) {
    
    if ( i != j) {
      
      f[[i]][[j]]<-prjsbj[[i]]  %in% prjsbj[[j]]
    }
      
    
 
    
  }
  
}
f


g<-list ()
for (i in names (prjsbj)) {
  
  g[[i]]<-sapply( f[[i]], function(x) prjsbj[[i]][x] )
}

g

unique(unlist (g))

names(sbj2) %in% unique(unlist (g))
