bixtech<-read.table("~/GitHub/2mode/data/minibix.txt",sep=";", 
                    header=T, stringsAsFactors = FALSE)
head(bixtech)
str(bixtech)
# bixtech[order(bixtech$prj),]
# bixtech[order(bixtech$sbj),]

myprj<-unique(bixtech$prj)

# prjsbj represents sbj for each prj

prjsbj<-sapply(myprj, function(x) bixtech$sbj[bixtech$prj==x])
length(prjsbj)
names(prjsbj)

sbj2<-table(bixtech$sbj)[table(bixtech$sbj)>1]
names(sbj2)

library (igraph)

gbixtech <- graph.data.frame(bixtech)
gbixtech

# set prj type <- TRUE
V(gbixtech)$type <- V(gbixtech)$name %in% bixtech[,1]




nnbix<-bipartite.projection(gbixtech)


snbix<-nnbix$proj2
pnbix<-nnbix$proj1

snbix


V(snbix)$name

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


mypcm<-clique.community(snbix,3)


sapply(mypcm, function (x) V(snbix)$name[x])
prjsbj

########################################################



k=3
clq <- cliques(snbix, min=k, max=k)

class(clq)
length(clq) # 374

clq[1:10]

edges <- c()
for (i in 1:3) {
    for (j in 1:3) {
        if ( length(unique(c(clq[[i]], clq[[j]]))) == k+1 ) {
            edges <- c(edges, c(i,j))
        }
    }
}


edges
clq [1:3]

edges [1:110]

length(edges)

clq.graph <- simplify(graph(edges))

plot(clq.graph )

V(clq.graph)$name <- seq_len(vcount(clq.graph))

comps <- decompose.graph(clq.graph)
comps[1]
unlist(clq[ V(comps[[1]])$name ])
clq[1:3]


lapply(comps, function(x) {
    unique(unlist(clq[ V(x)$name ]))})
