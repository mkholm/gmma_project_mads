options(width=200)
pcol=c(seq(4,7),seq(9,13),17,19,20,21,27,28,29)

args = commandArgs(trailingOnly=TRUE)
if (interactive()) {
    scan_var = NA
} else if (length(args) > 0) {
    scan_var = as.numeric(args[1])
    print(paste("Using scanning varianble =",scan_var))
} else {
    scan_var = NA
    print("No arguments to make_synthetic.r")
}

# Reference directory to previously saved synthetic data from, not gmma_result and Rosetta stuff
ref_dir = "."

aa_one = "ACDEFGHIKLMNPQRSTVWY"

# # Get reference GMMA results
# ref_gmma_file="../gmma_result.rda"
# print(sprintf("Loading %s",ref_gmma_file))
# load(ref_gmma_file)

# Generate all single mutants from a amino acid sequence string
all_single = function(aa_seq) {
    nres = nchar(aa_seq)
    resi = rep(seq(nres), each=20)
    resn = rep(strsplit(aa_seq,"")[[1]], each=20)
    rest = rep(strsplit(aa_one,"")[[1]], nres)
    mask = resn != rest
    return( paste0(resn[mask], resi[mask], rest[mask]) )
}

generate_n_mut = function(n, n_mut, pos_l, oversample=5) {
    # This is uniform in positions but uniform in substitutions if a position occurs the same number of times in pos_l
    #   as there are different substitutions at that position
    df = data.frame(mut_01 = sample(pos_l, oversample*n, replace=T))
    # Make n_mut columns of positions
    for (nm in seq(2,n_mut)) {
        df[,sprintf("mut_%02d",nm)] = sample(pos_l, oversample*n, replace=T)
    }
    # Only keep rows with substitutions at different positions
    unq_pos = apply(df, MARGIN=1, function(v) {length(unique(v))} )
    df = df[which(unq_pos==n_mut),]
    if (nrow(df) >= n) { 
        return(df[1:n,])
    } else {
	return(NA)
    }
}

build_random_lib = function(lib_compo, subst_list) {
    # lib_compo should have two columns that gives the number of variants "nvar" with a given number of substitutions "nmut"
    stopifnot("nmut" %in% colnames(lib_compo) & "nvar" %in% colnames(lib_compo))
    # each nmut sould only be given once
    stopifnot(length(unique(lib_compo$nmut)) == length(lib_compo$nmut))
    # substitutions should all be unique
    stopifnot(length(unique(subst_list)) == length(subst_list))

    # list of positions for substitutions
    pos_l = as.integer(substr(subst_list, 2, nchar(subst_list)-1))
    nres = max(pos_l)
    print(sprintf("Protein of %d residues with %d og %d possible single substitutions given", nres, length(subst_list), nres*19))
    
    # subst_list indices to lookup substitutions at a given position
    resi_index = lapply(seq(nres), function(i) {which(pos_l == i)} )

    # double mut calculation based on resi_index
    pos_subst_n = sapply(resi_index, length)
    m = outer(pos_subst_n, pos_subst_n)
    n_double = sum(m[upper.tri(m,diag=F)])
    print(sprintf("Substitution list can make %d of %d possible double mutants", n_double, (nres^2-nres)/2*19^2))
    
    var_list = c()
    if (1 %in% lib_compo$nmut) {
        i = which(lib_compo$nmut == 1)
	stopifnot(lib_compo[i,"nvar"] <= length(subst_list))
	var_list = sample(subst_list, lib_compo[i,"nvar"], replace=F)
	print(sprintf("Added %d of %d single mutants", lib_compo[i,"nvar"], length(subst_list)))
	# remove row with single mutants
	lib_compo = lib_compo[-i,]
    }
    for (mi in seq(nrow(lib_compo))) {
        n_mut = lib_compo[mi,"nmut"]
	n_var = lib_compo[mi,"nvar"]
	print(paste(n_mut, n_var))
	
	# data frame containing position to substitute in each variant
        var_df = generate_n_mut(n_var, n_mut, pos_l)

	# sort substitutions by position
	col_names = colnames(var_df)
	var_df = as.data.frame(t(apply(var_df, MARGIN=1, sort)))
	colnames(var_df) = col_names
	# order by first two substitutions beause I don't know how to order by all (these are always double or higher substituted variants)
	var_df[order(var_df$mut_01, var_df$mut_02),]

	# add random substitution index for each position
        for (cn in colnames(var_df)) {
            var_df[,paste0(cn,"_si")] = sapply(var_df[,cn], function(p) {sample(seq_along(resi_index[[p]]),1)})
        }

        # fill in substitutions at each position
	f = function(l,n_mut) {
	    paste( sapply(seq(n_mut), function(i){subst_list[resi_index[[l[i]]][l[i+n_mut]]]}), collapse=":" )
	}
        var = apply(var_df, MARGIN=1, f, n_mut=n_mut)
	var = unname(var)

        var_list = c(var_list, var)
	print(sprintf("Added %d %d-mutants, total variants %d", length(var), n_mut, length(var_list)))
    }
    return(var_list)
}

make_pairs = function(subst_list) {
    # Make all pair combinations excluding combinations of substitutions at same position
    pos_list = as.integer(substr(subst_list, 2, nchar(subst_list)-1))
    nres = max(pos_list)
    # subst_list indices to lookup substitutions at a given position
    resi_index = lapply(seq(nres), function(i) {which(pos_list == i)} )
    # lower triangle coordinates of position-pairs matrix
    lower_indices = data.frame(i=rep(seq(2,nres),seq(1,nres-1)), j=sequence(1:(nres-1)))
    # pairs of subst_list indices for all possible combinations of substitutions (not substitutions at same position)
    pairs = unlist(apply(lower_indices, MARGIN=1, function(l) { c(outer(resi_index[[l[2]]], resi_index[[l[1]]], paste, sep=":")) }))
    pairs_list = strsplit(pairs, ":")
    pairs_df = data.frame(i=as.numeric(sapply(pairs_list, "[", 1)), j=as.numeric(sapply(pairs_list, "[", 2)))
    return( paste(subst_list[pairs_df$i], subst_list[pairs_df$j], sep=":") )
}

# poor solution to make complete n subst var sets, but works: make all subst pairs, triples etc and concatenate (function "make_complete()")
# (will work for Laub data (only up to 4 subst/var), will not work for Hollfelder data ! (up to 6 subst per var, and different sets of subst for different positions!))
make_triples = function(subst_list){
  
  
  # Make all pair combinations excluding combinations of substitutions at same position
  pos_list = as.integer(substr(subst_list, 2, nchar(subst_list)-1))
  nres = max(pos_list)
  # subst_list indices to lookup substitutions at a given position
  resi_index = lapply(seq(nres), function(i) {which(pos_list == i)} )
  
  # modify below here for higher order subst var
  # lower triangle coordinates of position-pairs matrix
  lower_indices = t(data.frame(combn(unique(pos_list), 3))) # replace 2 with n
  # pairs of subst_list indices for all possible combinations of substitutions (not substitutions at same position)
  triples = unlist(apply(lower_indices, MARGIN=1, function(l) { apply((expand.grid(resi_index[[l[3]]], resi_index[[l[2]]], resi_index[[l[1]]])), 1, paste, collapse = ":")  }))
  triples_list = strsplit(triples, ":")
  triples_df = data.frame(i=as.numeric(sapply(triples_list, "[", 1)), j=as.numeric(sapply(triples_list, "[", 2)), k=as.numeric(sapply(triples_list, "[", 3)))
  return(paste((paste(subst_list[triples_df$i], subst_list[triples_df$j], sep=":")), subst_list[triples_df$k], sep = ":")  )
}

make_quads = function(subst_list){
  
  # Make all pair combinations excluding combinations of substitutions at same position
  pos_list = as.integer(substr(subst_list, 2, nchar(subst_list)-1))
  nres = max(pos_list)
  # subst_list indices to lookup substitutions at a given position
  resi_index = lapply(seq(nres), function(i) {which(pos_list == i)} )
  
  # modify below here for higher order subst var
  # lower triangle coordinates of position-pairs matrix
  lower_indices = t(data.frame(combn(unique(pos_list), 4))) # replace 2 with n
  # pairs of subst_list indices for all possible combinations of substitutions (not substitutions at same position)
  quads = unlist(apply(lower_indices, MARGIN=1, function(l) { apply((expand.grid(resi_index[[l[4]]], resi_index[[l[3]]], resi_index[[l[2]]], resi_index[[l[1]]])), 1, paste, collapse = ":")  }))
  quads_list = strsplit(quads, ":")
  quads_df = data.frame(i=as.numeric(sapply(quads_list, "[", 1)), j=as.numeric(sapply(quads_list, "[", 2)), k=as.numeric(sapply(quads_list, "[", 3)), l=as.numeric(sapply(quads_list, "[", 4)))
  return(paste(paste((paste(subst_list[quads_df$i], subst_list[quads_df$j], sep=":")), subst_list[quads_df$k], sep = ":"), subst_list[quads_df$l], sep = ":"))
}



make_complete <- function(subst_list, nmut){
  vars <- c()
  if (1 %in% nmut){
    vars <- append(vars, subst_list)
  }
  if (2 %in% nmut){
    vars <- append(vars, make_pairs(subst_list))
  }
  if (3 %in% nmut){
    vars <- append(vars, make_triples(subst_list))
  }
  if (4 %in% nmut){
    vars <- append(vars, make_quads(subst_list))
  }
  
  return(vars)
}




order_variants = function(variant_list) {
    # function to order substitutions in a variant and filter out variants with substitutions at same position
    order_var = function(v) {
        vl = strsplit(v,":")[[1]]
	p = as.numeric( substr(vl,2,nchar(vl)-1) )
	if (length(p)==length(unique(p))) {
	    return(paste(vl[order(p,decreasing=F)], collapse=":"))
	} else {
	    return("")
	}
    }
    new_list = sapply(variant_list, order_var)
    return(new_list[nchar(new_list) > 0])
    
    # unname also ?
    # could also order list here so single mut are first according to resi-taa, then double mut etc
}

# # Don't add terminal MS and K residues to keep Sarkisyan numbering
# gfp_ref="KGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELY"

# pyr1_wt="MPSELTPEERSELKNSIAEFHTYQLDPGSCSSLHAQRIHAPPELVWSIVRRFDKPQTYKHFIKSCSVEQNFEMRVGCTRDVIVISGLPANTSTERLDILDDERRVTGFSIIGGEHRLTNYKSVTTVHRFEKENRIWTVVLESYVVDMPEGNSEDDTRMFADTVVKLNLQKLATVAEAMARNSGDGSGSQVT"

laub_wt <- "AVST"

# #######################################################################################################################################
# #######################################################################################################################################
# ##
# ##  SUBSTITUTION GENERATION
# ##
# #######################################################################################################################################
# #######################################################################################################################################
# subst_param = list()
# extra = list()
# 
# 
# # 
# # 
# # # gfp_all_single = all_single(gfp_ref)
# # 
# # # List of substitutions
# # # subst_list = c("A1d", "A1e", "B2d", "B2e", "B2f", "C4c") # test list
# # # subst_list = gfp_all_single
# # # subst_list = sample(gfp_all_single, 1000, replace=F)
# # 
# # # # Sarkisyan substitutions
# # # subst_param$subst_mask = subst$gmma=="use"
# # # subst_list = rownames(subst[which(subst_param$subst_mask),])
# # 
# # # # Library with destabilizing pair-combinatorial backgrounds, aka s77
# # # stab_subst = rownames(subst_sig[1:50,])
# # # taa = substr(gfp_all_single,nchar(gfp_all_single),nchar(gfp_all_single)) 
# # # wtaa = substr(gfp_all_single,1,1) 
# # # iI2V = which(wtaa=="I" & taa=="V")
# # # iV2A = which(wtaa=="V" & taa=="A")
# # # iV2A = iV2A[iV2A != which(gfp_all_single == "V161A")]
# # # destab_subst = gfp_all_single[c(iI2V,iV2A)]
# # # destab_subst = destab_subst[! destab_subst %in% stab_subst]
# # # subst_list = c(stab_subst, destab_subst)
# # # # store for later usage
# # # extra[["stab_subst"]] = stab_subst
# # # extra[["destab_subst"]] = destab_subst
# # 
# # 
# # 
# # # Trimodal synthetic substitution effects, aka s75
# # subst_param$wt = laub_wt
# # subst_param$n_stab = 0
# # subst_param$n_neu = 19
# # subst_param$n_destab = 57
# # subst_param$n_ifs = 0
# # subst_list = sample(all_single(laub_wt), subst_param$n_stab+subst_param$n_neu+subst_param$n_destab+subst_param$n_ifs, replace=FALSE)
# # print(sprintf("Generated %d substotutions of Laub with %d stabilizing, %d neutral and %d destabilizing, and %d IFS",
# #               length(subst_list),subst_param$n_stab,subst_param$n_neu,subst_param$n_destab, subst_param$n_ifs))
# # 
# # 
# # 
# # 
# # # # scan trimodal distribution of n_stab,n_neu,n_destab
# # # extra$scan_var = round(scan_var)
# # # n = extra$scan_var
# # # print(sprintf(" RESULT_SCAN : %d IFSs", n))
# # # 
# # # 
# # # # Trimodal synthetic substitution effects, aka s75
# # # subst_param$wt = laub_wt
# # # subst_param$n_stab = 0
# # # subst_param$n_neu = 19 - n 
# # # subst_param$n_destab = 57
# # # subst_param$n_ifs = n
# # # subst_list = sample(all_single(laub_wt), subst_param$n_stab+subst_param$n_neu+subst_param$n_destab+subst_param$n_ifs, replace=FALSE)
# # # print(sprintf("Generated %d substotutions of Laub with %d stabilizing, %d neutral and %d destabilizing, and %d IFS",
# # #               length(subst_list),subst_param$n_stab,subst_param$n_neu,subst_param$n_destab, subst_param$n_ifs))
# # 
# # 
# # # Sort substitutions according to position and new amino acid
# # stopifnot(length(subst_list) == length(unique(subst_list)))
# # subst_list = subst_list[order( as.numeric(substr(subst_list,2,nchar(subst_list)-1)), substr(subst_list,nchar(subst_list),nchar(subst_list)))]
# # 
# # save(subst_list, subst_param, extra, file="synth_subst.rda")
# # # load(paste0(ref_dir,"/synth_subst.rda"))
# # 
# # #############################################
# # # Generate or calculate substitution effects (ddG's)
# # #############################################
# # ddG_param = list()
# # ddG = rep(NA, length(subst_list))
# # names(ddG) = subst_list
# # 
# # # # Sarkisyan-like data
# # # ddG[] = subst[which(subst_param$subst_mask),"ddG_glob"]
# # # # replace uncertain ddG's with random values drawn from a normal approximation of GFP ddG's
# # # subst_uncert = rownames(subst[which(subst_param$subst_mask & is.na(subst$rank)),])
# # # ddG[subst_uncert] = rnorm(length(subst_uncert), mean=mean(subst_sig$ddG_glob), sd=sd(subst_sig$ddG_glob))
# # 
# # 
# # # # scan ddG sd
# # # extra$scan_var = as.numeric(scan_var)
# # # n = extra$scan_var
# # # print(sprintf(" RESULT_SCAN : ddG sd  : %f  ", n))
# # 
# # # Trimodal synthetic substitution effects, aka s75
# # ddG_param$mean_stab = -1.0
# # ddG_param$mean_neu = 0.0
# # ddG_param$mean_destab = 2.0
# # ddG_param$mean_ifs = 10.0
# # ddG_param$sd = 0.4
# # ddG_param$cat = sample(c(rep(1,subst_param$n_stab), rep(2,subst_param$n_neu), rep(3,subst_param$n_destab), rep(4,subst_param$n_ifs)))
# # ddG[which(ddG_param$cat==1)] = rnorm(subst_param$n_stab, ddG_param$mean_stab, ddG_param$sd)
# # ddG[which(ddG_param$cat==2)] = rnorm(subst_param$n_neu, ddG_param$mean_neu, ddG_param$sd)
# # ddG[which(ddG_param$cat==3)] = rnorm(subst_param$n_destab, ddG_param$mean_destab, ddG_param$sd)
# # ddG[which(ddG_param$cat==4)] = rnorm(subst_param$n_ifs, ddG_param$mean_ifs, ddG_param$sd)
# # print(sprintf("Generated %d ddGs with observed mean %.3f and std %.3f",length(ddG),mean(ddG),sd(ddG)))
# # 
# # 
# # 
# # 
# # 
# # save(subst_list, subst_param, extra, ddG, ddG_param, file="synth_ddG.rda")
# # # load(paste0(ref_dir,"/synth_ddG.rda"))
# # 
# 
# #############################################
# print("Generate all possible pair couplings on ddG level")
# #############################################
# coupling_param = list()
# pairs = make_pairs(subst_list)
# couplings = rep(NA, length(pairs))
# names(couplings) = pairs
# 
# # Turn couplings off
# coupling_param$type = "off"
# couplings[] = rep(0.0, length(pairs))
# 
# # # Normal distributed couplings
# # coupling_param$type = "normal"
# # coupling_param$mean = 0.2
# # coupling_param$sd = 0.2
# # print(sprintf("Using normal-distributed pair-couplings with mean %.2f and sd %.2f",coupling_param$mean,coupling_param$sd))
# # couplings[] = rnorm(length(pairs), mean=coupling_param$mean, sd=coupling_param$sd)
# 
# print(sprintf("Summary of %d pair couplings",length(couplings)))
# print(summary(couplings))
# 
# save(subst_list, subst_param, extra, ddG, ddG_param, couplings, coupling_param, file="synth_subst.rda")
load(paste0(ref_dir,"/synth_subst.rda"))


#######################################################################################################################################
#######################################################################################################################################
##
##  VARIANT GENERATION
##
#######################################################################################################################################
#######################################################################################################################################
print("Generate variants")


# extra$scan_var = round(scan_var)
# n = extra$scan_var
# print(sprintf(" RESULT_SCAN : Added multiply of variants: %d ",n))


n_subst = length(subst_list)

# lib_compo = data.frame(nmut=c(1,2,3), nvar=c(n_subst, n_subst*n*0.5, n_subst*n*0.5))
# lib_compo = data.frame(nmut=c(1,2), nvar=c(n_subst, n_subst*n*0.5))
# lib_compo = data.frame(nmut=c(1,2), nvar=c(76, 2166))
# lib_compo = data.frame(nmut=c(1,2,3), nvar=c(n_subst, 2166, 27436))
lib_compo = data.frame(nmut=c(1,2,3,4), nvar=c(76, 2166, 27436, 130321))
# lib_compo = data.frame(nmut=c(1,2,3,4), nvar=c(n_subst, n_subst*2, n_subst*5, n_subst*10))

# random libarary
# var = build_random_lib(lib_compo, subst_list)






# # Sarkisyan variants, should match subst list with gmma=="use" and be connected etc.
# im = which(mutant$gmma=="use" & mutant$subst!="")
# var = mutant[im,"subst"]
# # var = order_variants(var)
# lib_compo = as.data.frame( table(mutant[im,"N_sub"]) )
# colnames(lib_compo) = c("nmut","nvar")

# # Combinatorial library
# stopifnot("stab_subst" %in% names(extra) & "destab_subst" %in% names(extra))
# destab_bg = c(extra[["destab_subst"]], make_pairs(extra[["destab_subst"]]))
# var = outer(extra[["stab_subst"]], destab_bg, paste, sep=":")
# var = order_variants(var)
# nmut = sapply(strsplit(var,":"), length)
# lib_compo = data.frame(table(nmut))
# colnames(lib_compo) = c("nmut","nvar")
# print(sprintf("Generated %d variants from %d substitutions (%d stabilizing) and %d destabilized backgrounds",
#               length(var), length(subst_list), length(extra[["stab_subst"]]), length(destab_bg)))





# # Combinatorial library
# stopifnot("cat" %in% names(ddG_param))
# # make pairs of all destabilizing substitutions
# # destab_bg =  make_pairs(subst_list[which(ddG_param$cat==3)])
# destab_bg =  c(subst_list[which(ddG_param$cat==3)], make_pairs(subst_list[which(ddG_param$cat==3)]))
# # down-sample destabilizing backgrounds
# si = sample(seq_along(destab_bg), length(destab_bg)/2, replace=F)
# destab_bg = destab_bg[si[order(si,decreasing=F)]]
#
# # make library of two destabilizing and one stabilizing or neutral
# var = outer(subst_list[which(ddG_param$cat!=3)], destab_bg, paste, sep=":")
# #var = unname(order_variants(var))
# var = c(destab_bg, unname(order_variants(var)))
#
# # determine n-mut distribution
# nmut = sapply(strsplit(var,":"), length)
# lib_compo = data.frame(table(nmut))
# colnames(lib_compo) = c("nmut","nvar")
# print(sprintf("Generated %d variants from %d substitutions (%d stabilizing+neutral) and %d destabilized backgrounds",
#               length(var), length(subst_list), sum(ddG_param$cat!=3), length(destab_bg)))
#







# complete combinatroial library of all possible combinations of subst/var (disregard nvar in lib_compo)
var <- make_complete(subst_list, lib_compo$nmut)




save(var, lib_compo, subst_list, file="synth_var.rda")

# # check if loaded variant are from same list of substitutions
# subst_list_preload = subst_list
# load(paste0(ref_dir,"/synth_var.rda"))
# stopifnot(all(subst_list_preload == subst_list))


#############################################
print("Calculate pair couplings per variant")
#############################################
ddG_pairs = rep(0.0, length(var))

all_pairs = function(s) {
    l=strsplit(s,":")[[1]]
    n=length(l)
    if (n==2) {
        return(s)
    } else if (n>2) {
        i = rep(seq(2,n),seq(1,n-1))
	j = sequence(1:(n-1))
	paste(l[j], l[i], sep=":")
    }
}

if (coupling_param$type != "off") {
    print(sprintf("Generating couplings of type: %s",coupling_param$type))
    
    # Generate substitution pairs for each variant
    var_subst_pairs = lapply(var, all_pairs)

    # Couplings on all pairs - there may be many 
    print(sprintf("Lookup and sum %d couplings",sum(sapply(var_subst_pairs, length))))
    ddG_pairs = sapply(var_subst_pairs, function(l) { sum(couplings[l]) } )
    print("done summing")
    vp = unlist(var_subst_pairs)
    uvp = unique(vp)
    print(sprintf("Variant library probes %d pairs with %d of %d possible unique pairs",length(vp),length(uvp),length(couplings)))
    print(sprintf("Couplings effect per variant: mean %.3f std %.3f",mean(ddG_pairs),sd(ddG_pairs)))
    print(sprintf("Unique observed couplings: mean %.3f std %.3f (input %.3f and %.3f)",mean(couplings[uvp]),sd(couplings[uvp]),coupling_param$mean,coupling_param$sd))
} else {
    print("Couplings are off")
}

# The coupling lookup list assumes that variant names are ordered according to position
stopifnot(all(! is.na(ddG_pairs)))


#############################################
print("Calculate variant stabilities")
#############################################

# dG_wt = cp$dG_wt # use estimated WT stability from Sarkisyan data
# dG_wt = 0.2756

# # scan dGwt
# extra$scan_var = scan_var
# n = extra$scan_var * 0.1
# print(sprintf(" RESULT_SCAN : dG_wt : %f  ", n))
# 
# dG_wt = n
# 

# Binary signal noise
extra$scan_var = round(scan_var)
n = extra$scan_var # FNR = n * 0.01

ddG_add = sapply(strsplit(var,":"), function(subst_list){ sum(ddG[subst_list]) })

# adjust dGwt to keep f_active constant with changing scan var
# quantile(ddG_add, f_active) + dG_wt = 0
# f_active = 0.0112 # 1659 active + n % FN of 160k variants
f_active = (1659/(1-n*0.01))/160000 # 1659 active + n % FN of 160k variants
dG_wt <- -(quantile(ddG_add, f_active)) 


dG = dG_wt + ddG_add + ddG_pairs



# Plot stability distributions
com_break = seq(floor(min(c(ddG_add,ddG_pairs,dG))), ceiling(max(c(ddG_add,ddG_pairs,dG))), length.out=101)
# quartz(width=8, height=10)

# rqeq for proper plot naming
extra$scan_var = round(scan_var)
n = extra$scan_var

# save plot for every scan var value run, collect in separate directory
if(!(dir.exists("stability_distributions"))){ 
  dir.create("stability_distributions")
}
# quartz(width=7, height=7.5)
png(sprintf("stability_distributions/stability_distributions_%f.png", extra$scan_var), width = 8, height = 10, units = "in", res = 100)

# png("stability_distributions.png", width=8, height=10, units = "in", res = 100)
par(mfcol=c(3,1))
s = sprintf("Additive effect per variant mean %.2f sd %.2f (per subst %.2f and %.2f)", mean(ddG_add), sd(ddG_add), mean(ddG), sd(ddG))
hist(ddG_add, breaks=com_break, main=s)
s = sprintf("Coupling per variant mean %.2f sd %.2f (per pair %.2f and %.2f)", mean(ddG_pairs), sd(ddG_pairs), coupling_param$mean, coupling_param$sd)
hist(ddG_pairs, breaks=com_break, main=s)
s = sprintf("Total per variant mean %.2f sd %.2f, dG_wt %.1f",mean(dG),sd(dG),dG_wt)
hist(dG, breaks=com_break, main=s)
# quartz.save("stability_distributions.png", type="png")
dev.off()

save(subst_list, subst_param, extra, ddG, ddG_param, couplings, coupling_param, lib_compo, var, dG, dG_wt, ddG_pairs, file="synth_bib.rda")
# load(paste0(ref_dir,"/synth_bib.rda"))


#######################################################################################################################################
#######################################################################################################################################
##
##  ASSAY SIGNAL GENERATION
##
#######################################################################################################################################
#######################################################################################################################################
print("Calculate variant brightness from stabilities")

signal_param = list()

# # Sarkisayn parameters
# signal_param$type = "Sigmoid Sarkisyan"
# signal_param$B_max = cp$B_max
# signal_param$B_D = cp$B_D
# signal_param$RT = settings$RT

# # Sigmoid, same as Sarkisyan fit
# B_pred_glob = function(dG) { e = exp(dG/signal_param$RT); (signal_param$B_max + signal_param$B_D*e) /(1.0+e) }
# B_synth = B_pred_glob

# Log-linear, Looks a lot like Sarkisyan GMMA
# signal_param$type = "Log-linear"
# B_synth = function(dG) { i = which(dG < 1.5); B = rep(5.0,length(dG)); B[i] = -10*dG[i]+20; log(B) } 
# B_synth = function(dG) { i = which(dG < 0.0); B = rep(5.0,length(dG)); B[i] = -10*dG[i]+20; log(B) } # this

# Log-linear, more bimodal
# signal_param$type = "Log-linear bimodal"
# Dis-continuert i ovregangen ser ud til at ddG'erne taber skalaen dvs. dG_wt er dårligt bestemt hvilket giver mening
# B_synth = function(dG) { i = which(dG < 0); B = rep(5.0,length(dG)); B[i] = -2*dG[i]+30; log(B) }

# Binary
signal_param$type = "Binary"
B_synth = function(dG) { i = which(dG > 0.0); B = rep(1.0,length(dG)); B[i] = 0.0; return(B)}

B_wt = B_synth(dG_wt)
B = B_synth(dG)


#############################################
print("Genearte Cauchy and misclassification noise on brightness level")
#############################################
# signal_param$noise_scale = 0.05
# noise = rcauchy(length(var), location=0, scale=signal_param$noise_scale)
# print("Noise summary (before capping)")
# print(summary(noise))

# # Sarkisyan inspired mis-classification noise, move variants between active and inactive
# nvar = length(var)
# misclass_signal = cp$B_mid-cp$B_D

# # High residual, active variants predicted to be inactive
# signal_param$misclass_rate_active = sum(mutant$residual < -misclass_signal, na.rm=T)/nrow(mutant)
# i_act = which(B > cp$B_mid)
# ii_act = sample(seq_along(i_act), signal_param$misclass_rate_active*nvar, replace=F)
# # B[i_act[ii_inact]] = cp$B_D
# noise[i_act[ii_act]] = noise[i_act[ii_act]] + cp$B_D - B[i_act[ii_act]]

# # High negative residual, inactive variants predicted to be active
# signal_param$misclass_rate_inact = sum(mutant$residual > misclass_signal, na.rm=T)/nrow(mutant)
# i_inact = which(B < cp$B_mid)
# ii_inact = sample(seq_along(i_inact), signal_param$misclass_rate_inact*nvar, replace=F)
# # B[i_inact[ii_inact]] = wt$signal
# noise[i_inact[ii_inact]] = noise[i_inact[ii_inact]] + wt$signal - B[i_inact[ii_inact]]

# print(sprintf("Mis-classification noise: Active -> inactive %.3f%% or %d variants. Inactive -> active %.3f%% or %d variants",
#               signal_param$misclass_rate_active*100, floor(signal_param$misclass_rate_active*nvar),
# 	      signal_param$misclass_rate_inact*100,  floor(signal_param$misclass_rate_inact*nvar)))


# Binary signal noise
extra$scan_var = round(scan_var)
n = extra$scan_var
print(sprintf("FN rate: %d", n))

signal_param$misclass_rate_active = n * 0.01

# signal_param$misclass_rate_active = 0.07
i_active = which(B > 0.5)
n_active = round(length(i_active)*signal_param$misclass_rate_active)
i_misclass_active = sample(i_active, n_active)

signal_param$misclass_rate_inact = 0.0001
i_inact = which(B < 0.5)
n_inact = round(length(i_inact)*signal_param$misclass_rate_inact)
i_misclass_inact = sample(i_inact, n_inact)

noise = rep(0.0, length(var))
noise[i_misclass_active] = -1
noise[i_misclass_inact] = +1
print(sprintf("Misclassified %d active and %d inactive variants",n_active,n_inact))

# # Bounds on noise
# signal_param$noise_low_lim = 0.0
# # Limit noise in lower bound, if brightness plus noise is lower than limit, set to limit
# i = which(B + noise < signal_param$noise_low_lim)
# noise[i] = signal_param$noise_low_lim - B[i]
# print(sprintf("Returned %d signal values lower than %.3f to %.3f",length(i),signal_param$noise_low_lim,signal_param$noise_low_lim))

# # Limit noise in high bound, if brightness plus noise is higher than limit, remove noise
# signal_param$noise_high_lim = 1.0
# i = which(B + noise > signal_param$noise_high_lim)
# noise[i] = 0
# print(sprintf("Returned %d signal values higher than %.3f to model value (zero noise)",length(i),signal_param$noise_high_lim))

B = B + noise

# Perhaps replace brightness from sarkisyan - unlikely to be many higher order variants

# Plot signal distributions
com_break = seq(floor(min(c(B,noise))), ceiling(max(c(B,noise))), length.out=101)
# quartz(widt=8, height=7)

# save plot for every scan var value run, collect in separate directory
if(!(dir.exists("signal_distributions"))){ 
  dir.create("signal_distributions")
}
# quartz(width=7, height=7.5)
png(sprintf("signal_distributions/signal_distributions_%f.png", extra$scan_var), width = 8, height = 7, units = "in", res = 100)
# png("signal_distributions.png", width = 8, height = 7, units = "in", res = 100)
par(mfcol=c(2,1))
s = sprintf("Brightness mean %.2f sd %.2f", mean(B), sd(B))
hist(B, breaks=com_break, main=s)
s = sprintf("Noise mean %.2f sd %.2f", mean(noise), sd(noise))
hist(noise, breaks=com_break, main=s)
# quartz.save("signal_distributions.png", type="png")
dev.off()

# plot synthetic data with model fittet to GFP data
# quartz(width=8, height=5)
png("synthetic_data.png", width=8, height=5, units = "in", res = 100)
if (signal_param$type == "Binary") { jitter=rnorm(length(B), 0.0, 0.02) } else { jitter=0} 
plot(dG, B+jitter, pch=".", xlab="stability", ylab="signal")
x = seq(-10,20,.1)
# B_pred_glob = function(dG) { e = exp(dG/settings$RT); (cp$B_max + cp$B_D*e) /(1.0+e) }
# lines(x, B_pred_glob(x), col=2, lwd=3)
lines(x, B_synth(x), col=2, lwd=2)
points(dG_wt, B_wt, pch=16, col=2)
# quartz.save("synthetic_data.png", type="png")
dev.off()


#############################################
# Dump everything
#############################################

save(subst_list, subst_param, ddG, ddG_param, couplings, coupling_param, lib_compo, var, dG, dG_wt, ddG_pairs, extra,
     B_synth, B_wt, B, noise, signal_param, file="make_synthetic.rda")

# Dump synthetic variant data in prism format
f = file("gmma_synthetic.txt", "wt")
write("# Synthetic data for GMMA", f)
write("# ", f)
# Set WT std to noise_cauchy_scale because this is used for ddG_gmma error estimate
write.table(data.frame(variant=c("WT",var), n_syn=1, bright=c(B_wt,B), std=c(0.0,noise)), file=f, row.names=F, quote=F)
close(f)


# Dump substitution effects for comparison
f = file("gmma_synthetic_true.txt", "wt")
write("# True substitution effect used for making synthetic GMMA data", f)
write("# ", f)
write.table(data.frame(variant=subst_list, ddG=ddG, row.names=NULL), file=f, row.names=F, quote=F)
close(f)


# Plot library composition and frection active
# First, build multi-mutant data frame - could use lib_compo here
n_mut = sapply(strsplit(var, ":"), length)
max_mut = max(n_mut)
nmut_df = data.frame(n_mut=seq(max_mut), n=0, n_act=0, f_act=NA)
t = table(n_mut)
# t_act = table(n_mut[which(B >= cp$B_mid)])
t_act = table(n_mut[which(B >= 0.5)])
for (i in nmut_df$n_mut) {
    if (! is.na(t[as.character(i)]))     nmut_df[i,"n"] = t[as.character(i)]
    if (! is.na(t_act[as.character(i)])) nmut_df[i,"n_act"] = t_act[as.character(i)]
    nmut_df[i,"f_act"] = nmut_df[i,"n_act"] / nmut_df[i,"n"]
}
nlib = sum(nmut_df[,"n"])

# plot distribution of effects, library composition, and number of variants per substitution
# quartz(height=10, width=7)

# save plot for every scan var value run, collect in separate directory
if(!(dir.exists("lib_compo"))){
    dir.create("lib_compo")
}
png(sprintf("lib_compo/%f_lib_compo.png", extra$scan_var), height = 10, width = 7, units = "in", res = 100)


# png("lib_compo.png", height = 10, width = 7, units = "in", res = 100)
par(mfcol=c(3,1), mar=c(5,3,4,2)+.1)
hist(ddG, breaks=50)
s = sprintf("Fraction active %.2f of %d",sum(nmut_df$n_act)/nlib,nlib)
plot(0, 0, col="white", xlim=c(1,max_mut), ylim=c(0,1.3), xaxp=c(1,max_mut,max_mut-1), xlab="N-mutant", ylab="", main=s)
points(nmut_df$n_mut, nmut_df$n/nlib, type="b", pch=20, lwd=2, col=2)
points(nmut_df$n_mut, nmut_df$f_act, type="b", pch=20, lwd=2, col=4)
legend("topleft", c("Composition","Fraction active"), pch=20, col=c(2,4), ncol=2)

par(mar=c(4,5,2,2)+.1)
t = table(unlist(strsplit(var, ":")))
t = t[order(t,decreasing=T)]
barplot(t, las=2, cex.names=.6, ylab="Observations")

# quartz.save("lib_compo.png", type="png")
dev.off()
