
biotech<-read.table("~/GitHub/2mode/data/biotech.txt",sep=";", 
                  header=T, stringsAsFactors = FALSE)
head(biotech)
str(biotech)
# biotech[order(biotech$prj),]
# biotech[order(biotech$sbj),]

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
# tkplot(gbiotech)
# mylay<-tkplot.getcoords(1)


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

# tkplot(snbio)
# snbio

# mylay2<-tkplot.getcoords(4)


################
v_col<-rep("gold",78)

v_col[V(snbio)$name %in% names(sbj2)]<-c("brown")
co<-layout.fruchterman.reingold(snbio)

plot.igraph(snbio, layout=co, edge.arrow.size=0.1,
            vertex.label.cex=.7,vertex.color=v_col,)

######################

snbio

##
df_snbio<-get.data.frame(snbio)
df_pnbio<-get.data.frame(pnbio)

str(df_snbio)

#write.table(df_snbio,file = "sbj_3bio_el.txt",  sep = " ",
#            eol = "\n", na = "NA", dec = ",", row.names = TRUE)
#write.table(df_pnbio,file = "prj_3bio_el.txt",  sep = " ",
#            eol = "\n", na = "NA", dec = ",", row.names = TRUE)

### subnetwork of doubles


sbj2to<-df_snbio[df_snbio$weight>1,]

sbj2tot<-append(sbj2to$from,sbj2to$to)

sbj2pu<-unique(sbj2tot)
table(sbj2tot)



g2 <- induced.subgraph(graph=snbio,vids=names(sbj2))

v_col<-rep("white",12)
v_col[V(g2)$name%in%sbj2pu]<-c("red")


plot(g2,vertex.color=v_col )
# ----------------------------------------------
snbio 

E(snbio)$weight


# clustering coefficient
tt<-transitivity (snbio,type=c("local" ))
str(tt)
dd<-degree(snbio)

tt1<-cbind(V(snbio)$name,round(tt,2),dd)
tt1<-as.data.frame(tt1, stringsAsFactors = FALSE)
str(tt1)

nprj1<-cbind(names(table(biotech$sbj)), as.numeric(table(biotech$sbj)))
nprj1<-as.data.frame(nprj1, stringsAsFactors = FALSE)
str(nprj1)

ttnprj<-merge (tt1,nprj1, by="V1", stringsAsFactors = FALSE)
str(ttnprj)
ttnprj<-cbind(ttnprj, c("no2prj"), stringsAsFactors = FALSE)
str(ttnprj)
colnames(ttnprj)<-c("subjects","cluindx","degree","nprj","s2ub")
colnames(ttnprj)

zz<-which(ttnprj$nprj>1)
ttnprj$s2ub[zz]  <-ttnprj$subjects[zz]

plot(ttnprj$nprj ,ttnprj$degree )

library (ggplot2)

ggplot(ttnprj, aes(cluindx, nprj) ) + 
    geom_point(size = 10, color="gold",show_guide = FALSE) + 
    scale_x_discrete(breaks=seq(0, 1, 0.10)) +
    ggtitle("clustering index VS n.projects") +
    ylab("n.prj") + xlab("index") +
    geom_text(aes(label=s2ub), size=3, , hjust=0.5, vjust=3, color="red")



# betweenness

bb<-betweenness(snbio, normalized = T, weights=NA  )#freeman

V(snbio)$name[bb>0]

cc<-constraint(snbio)
summary(cc)

plot(tt,bb)


#cliques(snbio,min=4,max=18)
length(largest.cliques(snbio)[[1]])
maximal.cliques.count(snbio)


# maximal clique  
myclique<-maximal.cliques(snbio)
myclique_list<-sapply(myclique, function(x) V(snbio)$name[x])
length(myclique_list)

##########PLOT ####################
blockmax=length(myclique)
mycolors<-rainbow(blockmax)
v_col<-rep("white",78)


j<-0
for (i in 1:blockmax){
  j<-j+1
  v_col[myclique[[i]] ]<- mycolors[j]
}
#v_col[V(snbio)$name %in% names(sbj2)]<-c("white")

plot.igraph(snbio, layout=layout.fruchterman.reingold, edge.arrow.size=0.1,
            vertex.label.cex=.7,vertex.color=v_col)





library (igraph)


####################################################################
# clique.community FUNCTION ##############

clique.community <- function(graph, k) {
    clq <- cliques(graph, min=k, max=k) # my correction with maximal
    edges <- c()
    for (i in seq_along(clq)) {
        for (j in seq_along(clq)) {
            if ( length(unique(c(clq[[i]], clq[[j]]))) == k+1 ) {
                edges <- c(edges, c(i,j))
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

# clique  percolation is a community detection method developed by Gergely Palla

cli_com3<-clique.community (snbio,3)
cli_com4<-clique.community (snbio,4)
cli_com5<-clique.community (snbio,5)


length(cli_com3)
mean(unlist(cli_com3)%in%c(1:78))

sapply(cli_com4, function (x) V(snbio)$name[x])

blockcpm=length(cli_com4)
mycolors<-rainbow(blockcpm)
v_col<-rep("white",78)


j<-0
for (i in 1:blockcpm){
    j<-j+1
    v_col[cli_com4[[i]] ]<- mycolors[j]
}
#v_col[V(snbio)$name %in% names(sbj2)]<-c("white")

plot.igraph(snbio, layout=layout.fruchterman.reingold, edge.arrow.size=0.1,
            vertex.label.cex=.7,vertex.color=v_col)






######################################################################



myclique_perco<-sapply(cli_com, function(x) V(snbio)$name[x])
myclique_perco

blockper=length(myclique_perco)
mycolors<-rainbow(blockper)
v_col<-rep("white",78)


j<-0
for (i in 1:blockper){
  j<-j+1
  v_col[cli_com[[i]] ]<- mycolors[j]
}
#v_col[V(snbio)$name %in% names(sbj2)]<-c("white")

plot.igraph(snbio, layout=layout.fruchterman.reingold, edge.arrow.size=0.1,
            vertex.label.cex=.7,vertex.color=v_col)


#####################################################
# COHESIVE BLOCKS
#
#
################################################


cohblk_snbio<-cohesive.blocks(snbio)

print(cohblk_snbio) #it says everything you need



length(cohblk_snbio)
block_snbio<-blocks(cohblk_snbio)

plot(hierarchy(cohblk_snbio))
plotHierarchy(cohblk_snbio)

ssblock_snbio<-cohesion(cohblk_snbio)

maxcohesion(cohblk_snbio)


plot(cohblk_snbio,snbio)


table(ssblock_snbio)

blockGraphs(cohblk_snbio)

V(snbio)$name[block_snbio[[7]]]


# plotting level 4 
block4<-block_snbio[4:11]

blockcpm=length(block4)
mycolors<-rainbow(blockcpm)
v_col<-rep("white",78)


j<-0
for (i in 1:blockcpm){
  j<-j+1
  v_col[block4[[i]] ]<- mycolors[j]
}
co<-layout.fruchterman.reingold(snbio)
plot.igraph(snbio, layout=co, edge.arrow.size=0.1,
            vertex.label.cex=.7,vertex.color=v_col)



########### check clique and projects #################################
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
myclique_list[[7]]




########### check clique and projects TWO ################################
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
myclique_list[[7]]




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
      
        if (mean(prjsbj[[i]]  %in% prjsbj[[j]]) >0) {
            
            f[[i]][[j]]<-prjsbj[[i]]  %in% prjsbj[[j]]
        }
    }
    
  }
  
}

f

ff<-list ()
f[[1]]

for (i in names (prjsbj)) {
    
    ff[[i]]<-sapply( f[[i]], function(x) sum(x) )
}

ff
str(ff)
table(unlist(ff))



g<-list ()
for (i in names (prjsbj)) {
  
  g[[i]]<-sapply( f[[i]], function(x) prjsbj[[i]][x] )
}

g[[1]][[1]]

unique(unlist (g))

choose(8,2)

names(sbj2) %in% unique(unlist (g))

