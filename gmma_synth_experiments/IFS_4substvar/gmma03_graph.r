# Copyright (C) 2017-2019 Kristoffer Enoe Johansson <kristoffer.johansson@bio.ku.dk>
# This file (gmma03_graph.r) is part of the GMMA project

options(width=200, digits=4, stringsAsFactors=F)

args = commandArgs(trailingOnly=TRUE)
if (interactive()) {
    # If data is excluded in the initial fit, also exclude them here
    file="gmma_fit_init.rda"
    # ...or use all data
    # file="gmma_structured.rda"
} else if (length(args) < 1) {
    print("")
    print("usage: Rscript gmma03_graph.r <gmma_structured.rda>")
    quit(save="no")
} else {
    file = args[1]
}
print(sprintf("Read %s",file))
load(file)

nmut = length(mutant[,1])
nsubst = length(subst[,1])
nres = length(residue[,1])

################################################################################
# Exclude substitutions from global analysis
################################################################################
subst$gmma = "use"
mutant$gmma = "use"

# Function to set mutant$gmma and subst$gmma
exclude_data = function(si, tag, overwrite=FALSE) {
    subst_loc = subst
    mutant_loc = mutant
    if (! tag %in% levels(subst_loc$gmma)) {
        levels(subst_loc$gmma) = c(levels(subst_loc$gmma), tag)
        stopifnot(! tag %in% levels(mutant_loc$gmma))
        levels(mutant_loc$gmma) = c(levels(mutant_loc$gmma), tag)	
    }
    mutants_removed = 0
    if (overwrite) {
        subst_loc[si,"gmma"] = tag
        mi = unique(unlist(subst_mut_indices[si]))
	mutant_loc[mi,'gmma'] = tag
	mutants_removed = length(mi)
    } else {
        sii = which(subst_loc[si,"gmma"] == "use")
        subst_loc[si[sii],"gmma"] = tag
        mi = unique(unlist(subst_mut_indices[si[sii]]))
	mii = which(mutant_loc[mi,"gmma"] == "use")
	mutant_loc[mi[mii],"gmma"] = tag
	mutants_removed = length(mi)
    }
    assign("subst", subst_loc, envir=.GlobalEnv)
    assign("mutant", mutant_loc, envir=.GlobalEnv)
    return(mutants_removed)
}

if ('init_m' %in% colnames(subst)) {
    si = which(subst$init_m == 'nonsense')
    mutants_removed = exclude_data(si, "nonsense")
    print(sprintf("Excluded %d nonsense substitutions and %d dependent variants", length(si), mutants_removed))

    si = which(subst$init_m == 'ifs')
    mutants_removed = exclude_data(si, "ifs")
    print(sprintf("Excluded %d ifs substitutions and %d dependent variants", length(si), mutants_removed))

    print("subst$gmma:")
    print(table(subst$gmma))
    print("mutant$gmma:")
    print(table(mutant$gmma))
} else {
    print("Do not exclude substitutions from graph")
}

################################################################################
# Functions
################################################################################
require(Matrix)
require(igraph)

# disable_subst_and_adj_mut = function(graph, subst_list, gmma_tag=0) {
#     s_mask = (vertex_attr(graph, "gmma", subst_list) > 0)
disable_subst_and_adj_mut = function(graph, subst_list, gmma_tag) {
    s_mask = (vertex_attr(graph, "gmma", subst_list) == "use")
    graph = set_vertex_attr(graph, "gmma", index=subst_list[s_mask], value=gmma_tag)
    adj_mut = c()
    for (sn in subst_list) {
        adj_mut = append(adj_mut, neighbors(graph, sn)$name)
    }
    # am_mask = (vertex_attr(graph, "gmma", adj_mut) > 0)
    am_mask = (vertex_attr(graph, "gmma", adj_mut) == "use")
    print(sprintf("Disabled %d substitutions (%d already disabled), %d adjecent mutants (%d already disabled)",length(subst_list),sum(! s_mask),length(adj_mut),sum(! am_mask)))
    graph = set_vertex_attr(graph, "gmma", index=adj_mut[am_mask], value=gmma_tag)
    return(graph)
}

# Consider to use gmma_tag as the "use" tag, e.g. gmma1, gmma2 and then use "hanging" and "disconnected" for discarded nodes
# disable_disconnected_and_hanging = function(graph, gmma_tag=0) {
disable_disconnected_and_hanging = function(graph, gmma_tag) {
    # Assign gmma_tag to subst nodes that are disconnected from the largest connected graph or hanging
    # Only the gmma attribute is changed, not the graph edges
    # Hanging substitution nodes (substitution only in one variant) and the adjecent variant does
    #     not contribute to the global fit because the unrestrained ddG parameter may adapt to any
    #     change in the other parameteers. These parameters are completely overfitted and may be
    #     calculated after the global fit of intangled parameters.
    pass = 1
    hanging_subst = c(0,0)
    while (length(hanging_subst) > 0) {
        print(sprintf("==== Pass %d detetion of disconnected and hanging nodes",pass))

        # Subgraph of non-disabled nodes to identify disconnected and hanging
        # subgraph = induced.subgraph(graph, (vertex_attr(graph,"gmma") > 0))
        subgraph = induced.subgraph(graph, (vertex_attr(graph,"gmma") == "use"))
    
        # Disable disconnected nodes
        dg = decompose.graph(subgraph)
	gmma_indices = c()
        if (length(dg) > 1 ) {
            disconnected = c()
            # for (gi in seq(2,length(dg))) {
            for (gi in seq(length(dg))) {
                nodes = V(dg[[gi]])$name
		if (length(nodes) > 2) {
                    gmma_indices = append(gmma_indices, gi)
		} else {
                    # mask = vertex_attr(graph, "gmma", nodes) > 0
	            mask = vertex_attr(graph, "gmma", nodes) == "use"
                    disconnected = append(disconnected, nodes[mask])
		}
	    }
	    stopifnot(length(gmma_indices)==1)
            # Disable subgraph disconnected in graph 
            # graph = set_vertex_attr(graph, "gmma", index=disconnected, value=gmma_tag)
            graph = set_vertex_attr(graph, "gmma", index=disconnected, value="disconnected")
            print(sprintf("Disabled %d disconnected nodes:",length(disconnected)))
            print(disconnected)
        } else {
            print("Graph is connected")
        }

        # Find hanging substitution nodes (hanging variants are fine)
        i_hang_subst = which(degree(subgraph, V(subgraph))==1 & V(subgraph)$type)
        hanging_subst = vertex_attr(subgraph, "name", i_hang_subst)
        print(sprintf("Found %d hanging substitution node(s)",length(hanging_subst)))
        # graph = disable_subst_and_adj_mut(graph, hanging_subst, gmma_tag)
        graph = disable_subst_and_adj_mut(graph, hanging_subst, "hanging")
	
        pass = pass+1
    }
    return(graph)
}

################################################################################
# Network construction
################################################################################
print("Construct graph")
# Construct sparce incidence matrix, with mut_list vs subst indices (row vs column)
i = unlist(mapply(function(i,n) {rep(i,each=n)}, seq(length(mut_list)), sapply(mut_list,length)))
incidence = sparseMatrix(i=i, j=unlist(mut_subst_indices), x=1) # the 'x=1' results in a dgCMatrix class being created, otherwise igraph will crash on the boolean ngCMatrix type (missing col 3 with values)
colnames(incidence) = rownames(subst)
rownames(incidence) = rownames(mutant)

graph = graph_from_incidence_matrix(incidence)

subst_names = vertex_attr(graph, "name", which(V(graph)$type))
graph = set_vertex_attr(graph, "gmma", index=subst_names, value=subst[subst_names,'gmma'])

mutant_names = vertex_attr(graph, "name", which(! V(graph)$type))
graph = set_vertex_attr(graph, "gmma", index=mutant_names, value=mutant[mutant_names,'gmma'])
graph = set_vertex_attr(graph, "signal", index=mutant_names, value=mutant[mutant_names,'signal'])

# # Check if mutants or substitutions have been assigned with zero/negative values to gmma column
# table(vertex_attr(graph,"gmma"))
# # subst_names[which(vertex_attr(graph,"gmma",subst_names) < 1)]
# subst_names[which(vertex_attr(graph,"gmma",subst_names) != "use")]

################################################################################
# Consider to first detect isolated subgraphs of some minimum size and replace "use" with "gmma1", "gmma2", etc
################################################################################

################################################################################
# Network cleaning
################################################################################
print("Clean disconnected and hanging nodes")
# graph = disable_disconnected_and_hanging(graph, gmma_tag="exclude")
graph = disable_disconnected_and_hanging(graph)

# This returns the largest graph and disregards the case where the disconnected nodes constitutes a connected subgraph.
# In the future, these should have gmma=="use" to signal that a global fit should be carried out for each of these subgraphs

# Take the gmma info back to the data frames
subst$gmma = vertex_attr(graph, "gmma", rownames(subst))
mutant$gmma = vertex_attr(graph, "gmma", rownames(mutant))

print("subst$gmma:")
print(table(subst$gmma))
print("mutant$gmma:")
print(table(mutant$gmma))
print("")

################################################################################
# Re-parameterization
################################################################################
# If substitutions only occur in groups they can only be estimated as one parameter for the group
# Consider the global fit and individual estiamtions as two separate groups
subst$gmma_reparam = NA

find_co_occuring = function(graph, obs=NA) {
    subst_mask=vertex_attr(graph,"type")

    # Matrix which is TRUE if two substitutions co-occur just once ()
    print("generating neighbor mask")
    neighbor_mask = distances(graph, v=which(subst_mask), to=which(subst_mask), weights=NA) < 3
    # This should be symmetric so that taking the upper triangle is ok
    stopifnot(all(neighbor_mask == t(neighbor_mask)))

    # Data per substitution
    print("generating obs vector")
    obs = degree(graph, which(subst_mask))
    ns = length(obs)

    # Matrix with obs in rows
    m = matrix(obs, nrow=ns, ncol=ns, byrow=T)

    # Upper triangle elements are TRUE if a pair of substitutions have the same number of observations and co-occurs at least once
    mmask = m == t(m) & upper.tri(m, diag=F) & neighbor_mask
    # In 1D, which elements are TRUE
    mi = which(mmask)

    # If no pairs are found, return empty list
    if (length(mi) < 1) return(list())

    # A new matrix with the indices of the upper triangle elements that are TRUE
    pairs = matrix(c(row(m)[mi],col(m)[mi]), ncol=2, nrow=length(mi) )

    oi=unique(as.vector(pairs))
    oi=oi[order(oi)]
    
    print("generating neighbor cache")
    n_cache = lapply(names(obs)[oi], neighbors, graph=graph)
    names(n_cache) = names(obs)[oi]

    print("Search for pairs that always co-occur")
    co_pi = c()
    for (pi in seq(dim(pairs)[1])) {
        n1 = n_cache[names(obs)[pairs[pi,1]]]
        n2 = n_cache[names(obs)[pairs[pi,2]]]
        if (all(unname(n1)[[1]] == unname(n2)[[1]])) {
	    co_pi = append(co_pi,pi)
        }
    }

    # Further clustering
    ret = list()
    ngroup = 0
    all_s = c()
    for (pi in co_pi) {
        s1 = names(obs)[pairs[pi,1]]
	s2 = names(obs)[pairs[pi,2]]
        if (s1 %in% all_s | s2 %in% all_s) {
	    # The above can only be true if ngroup>0 so seq(ngroup=0) will never happen
	    for (li in seq(ngroup)) {
	        if (s1 %in% ret[[li]] | s2 %in% ret[[li]]) {
		    print(paste("add",s1,s2))
		    if (! s1 %in% ret[[li]]) { ret[[li]] = append(ret[[li]], s1) }
		    if (! s2 %in% ret[[li]]) { ret[[li]] = append(ret[[li]], s2) }
		    break
		}
	    }
	} else {
	    ngroup = ngroup+1
	    ret[[ngroup]] = c(names(obs)[pairs[pi,]])
	}
        all_s = append(all_s, s1)
        all_s = append(all_s, s2)	
    }
    return(ret)
}

# # Version using projections which are slow for many edges
# # Consider all pairs of substitutions that co-occur, i.e. edges in the subst projected network
# # The graph of individually estimated substitutions
# # graph0 = induced.subgraph(graph, (vertex_attr(graph,"gmma") == 0))
# graph0 = induced.subgraph(graph, (vertex_attr(graph,"gmma") == "exclude"))
# graph0_proj = bipartite_projection(graph0, multiplicity=T)
# graph0_subst = graph0_proj$proj2
# edge_mat = ends(graph0_subst, E(graph0_subst))
# for (ei in seq(dim(edge_mat)[1])) {
#     if (all(neighbors(graph0, edge_mat[ei,1]) == neighbors(graph0, edge_mat[ei,2]))) {
#         print(edge_mat[ei,])
#     }
# }

group_id = 0
# for (gi in seq(0,max(subst$gmma))) { gmma==0 is no longer a fitting group
# for (gi in seq(max(subst$gmma))) { # gi also used in nested loop!
#     print(sprintf("Check for co-occuring subst in fitting group %d",gi))
    # Graph of each GMMA
    # graph_gmma = induced.subgraph(graph, (vertex_attr(graph,"gmma") == gi))
    graph_gmma = induced.subgraph(graph, (vertex_attr(graph,"gmma") == "use"))
    # Groups of co-occuring substitutions in this GMMA
    co_groups = find_co_occuring(graph_gmma)
    print(sprintf("Found %d groups of fully co-occuring substitutions",length(co_groups)))
    # Assign an id which is unique across all fitting sets
    for (gi in seq_along(co_groups)) {
        group_id = group_id+1
        subst[co_groups[[gi]],'gmma_reparam'] = group_id
    }
# }

################################################################################
# Checks and tests
################################################################################
print("")
# print("Check that all subst with gmma > 0 have data with gmma > 0")
print("Check that all subst with gmma==use have data with gmma==use")
# obs_gmma1 = sapply(subst_mut_indices, function(l) { sum(mutant[l,'gmma'] > 0) })
obs_gmma1 = sapply(subst_mut_indices, function(l) { sum(mutant[l,'gmma'] == "use") })
# # Any subst without such data and with gmma > 0?
# si = which(obs_gmma1 == 0 & subst$gmma > 0)
# Any subst without such data and with gmma==use?
si = which(obs_gmma1 == 0 & subst$gmma == "use")
print(si)
stopifnot(length(si) == 0)

# this is checked before the global fit
# # print("Check that all mutants with gmma > 0 have estimated substitution ddG's")
# # all_subst_estimated = sapply(mut_subst_indices, function(l) { all(subst[l,"gmma"]>0) })
# # stopifnot(all(all_subst_estimated[which(mutant$gmma > 0)]))
# print("Check that all mutants with gmma==use have estimated substitution ddG's")
# all_subst_estimated = sapply(mut_subst_indices, function(l) { all(subst[l,"gmma"]=="use") })
# stopifnot(all(all_subst_estimated[which(mutant$gmma == "use")]))

# # This paragrapf may be skiped
# # gmma==0 is no longer a fitting group since the initial fits are made before this graph analysis
# print("")
# print("Check that all subst with gmma >= 0 have data with gmma >= 0")
# obs_gmma0 = sapply(subst_mut_indices, function(l) { sum(mutant[l,'gmma'] >= 0) })
# # Any subst without such data and with gmma >= 0?
# si = which(obs_gmma0 == 0 & subst$gmma >= 0)
# print(si)
# stopifnot(length(si) == 0)

################################################################################
# Dump
################################################################################
# REMEMBER that 'graph' still contains the disconnected and hanging nodes
graph_clean = induced.subgraph(graph, V(graph)[V(graph)$gmma=="use"] )

# consider to save more, e.g. the un-cleaned graph or original data frames like residue
# only subst$gmma, subst$gmma_reparam and mutant$gmma are changed in this script
# version=2 is necessary for R versions <3.5 to be able to read it (binf version 3.1)
save(mutant, subst, graph_clean, file="gmma_graph.rda", version=2)
print("Saved gmma_graph.rda")

################################################################################
# Plots
################################################################################
quit()

# Make a plot of all subst nodes with weighted edges and colored according to the gmma attrbute, can I reduce the number of edges somehow..?
# gp = bipartite_projection(graph)
gcp = bipartite_projection(graph_clean)

breaks = seq(0,1000,10)
h  = hist(degree(graph_clean)[which(V(graph_clean)$type)], breaks=breaks, plot=F)
hp = hist(degree(gcp$proj2), breaks=breaks, plot=F)

# Degree distribution of subst nodes and projected subst nodes
quartz(height=6, width=7)
plot(0,0,col="white", xlab="Substitution node degree", ylab="log(counts)", xlim=c(1,1000), ylim=c(0,6))
lines(h$mids, log(h$counts), type="o", pch=20, col=2)
lines(hp$mids, log(hp$counts), type="o", pch=20, col=3)
legend("topright", c("Substitution nodes","Projected substitution graph"), pch=20, lty=1, col=c(2,3))
quartz.save(file="subst_degree.png", type="png")

g = gcp$proj2
table(E(g)$weight)


subst$gmma = vertex_attr(graph, "gmma", rownames(subst))
# remove disconnected 
g = delete_vertices(g, which(degree(g)==0))
# color via attribute according to subst$gmma
V(g)$color = "lightgray"
V(g)[V(g)$gmma=="disconnected"]$color = "red"
V(g)[V(g)$gmma=="hanging"]$color = "blue"
V(g)[V(g)$gmma=="nonsense"]$color = "yellow"
V(g)[V(g)$gmma=="ifs"]$color = "magenta"

# Large Graph Layout (http://lgl.sourceforge.net)
l = layout_with_lgl(g)
ll = l
ll[,1] = (ll[,1]-mean(ll[,1]))/sd(ll[,1])
ll[,2] = (ll[,2]-mean(ll[,2]))/sd(ll[,2])

V(g)$obs = subst[V(g)$name,"obs"]
V(g)$ddG = subst[V(g)$name,"init_ddG"]
V(g)$ddG[is.na(V(g)$ddG)] = 10

plot(g, rescale=F, layout=ll*.5, vertex.size=2*sqrt(V(g)$obs)/6, vertex.label = NA, edge.width=(E(g)$weight)/4*(E(g)$weight > 1), egde.color=gray.colors(max(E(g)$weight),start=0.7,end=0.4)[E(g)$weight])
# plot(g, rescale=F, layout=ll*.5, vertex.size=20/(V(g)$ddG+4),    vertex.label = NA, edge.width=(E(g)$weight)/4*(E(g)$weight > 1), egde.color=gray.colors(max(E(g)$weight),start=0.7,end=0.4)[E(g)$weight])

# Show where nodes are in the network by making a lightgray contour of the network and then only plot selected nodes with labels
# vertex.frame.color="lightgray", edge.color="white"


#=====================================
#===  Plot the disabled subgraphs  ===
vi = which(vertex_attr(graph,"gmma") != "use" )
disabled_vertices = vertex_attr(graph, "name", vi)
disabled_graph = induced_subgraph(graph, vi)

par(ask=T)
plotted_vertices = c()
for (v in disabled_vertices) {
    # if (v %in% plotted_vertices) {print(sprintf("%s already plotted",v)); next} else {print(sprintf("%s plot",v))}
    if (v %in% plotted_vertices) next
    # find names of all disabled vertices that are connected to v
    cv = vertex_attr(disabled_graph, "name", subcomponent(disabled_graph, v))
    print(sprintf("Found %d connected disabled nodes from %s",length(cv),v))
    # print(paste(cv, collapse=" "))
    plotted_vertices = append(plotted_vertices, cv)
    if (length(cv) <= 2) { print(paste(V(graph)[cv]$gmma,collapse=" ")); next }
    # make a graph of these and the first eighbors
    edge_list = E(graph)[from(cv)]  # see ?"[.igraph.es"
    g = subgraph.edges(graph, edge_list)
    print(g)
    
    V(g)$color = "lightgray"
    V(g)[V(g)$gmma=="disconnected"]$color = "red"
    V(g)[V(g)$gmma=="hanging"]$color = "blue"
    V(g)[V(g)$gmma=="nonsense"]$color = "yellow"
    V(g)[V(g)$gmma=="ifs"]$color = "magenta"
    V(g)$shape=c("circle","square")[V(g)$type+1]

    if (length(V(g)) < 1000 ) {
        print("Plotting...")
        l = layout_with_fr(g)
        n = V(g)$name
        n[which(!V(g)$type)] = NA  # no label on mutant nodes
        plot(g, layout=l, vertex.label=n, vertex.size=c(7,14)[V(g)$type+1], vertex.label.cex=.6, vertex.label.font=2)
        V(g)[setdiff(V(g)$name,cv)]$shape = "none"
        plot(g, layout=l, vertex.label=n, vertex.size=c(7,14)[V(g)$type+1], vertex.label.cex=.6, vertex.label.font=2, add=T)
        print("...plot done")
    } else {
        # gp = bipartite_projection(g)
	# l = layout_with_fr(gp$proj2)
	# plot(gp$proj2, layout=l, vertex.label=NA, vertex.size=14, vertex.label.cex=.6, vertex.label.font=2)
	# # edge.width=sqrt(E(gp$proj2)$weight)

        # g2 = induced.subgraph(graph, cv)
        # V(g2)$color = "lightgray"
        # V(g2)[V(g2)$gmma=="exclude"]$color = "red"
        # V(g2)[V(g2)$gmma=="nonsense"]$color = "yellow"
        # V(g2)[V(g2)$gmma=="ifs"]$color = "magenta"
        # V(g2)$shape=c("circle","square")[V(g2)$type+1]
	# print(g2)
	# # l = layout_with_fr(g2)
        # # plot(g2, layout=l, vertex.label=n, vertex.size=c(7,14)[V(g2)$type+1], vertex.label.cex=.6, vertex.label.font=2)

	# g2p = bipartite_projection(g2)
	# print(g2p$proj2)
	# l = layout_with_fr(g2p$proj2)
	# # plot(g2p$proj2, layout=l, vertex.size=8, vertex.label.cex=.5, vertex.label.font=2)
	# plot(g2p$proj2, rescale=F, layout=l*0.2, vertex.size=5*sqrt(degree(g2p$proj2)), vertex.label.cex=.5, vertex.label.font=2, edge.width=(E(g2p$proj2)$weight)*(E(g2p$proj2)$weight > 1))
    }
}
