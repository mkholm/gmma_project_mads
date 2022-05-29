# Copyright (C) 2017-2019 Kristoffer Enoe Johansson <kristoffer.johansson@bio.ku.dk>
# This file (gmma05_analysis.r) is part of the GMMA project

options(width=300, digits=4, stringsAsFactors=F)
pcol=c(seq(6),9,10,11,12,17,19,20,21,27,28,29)

args = commandArgs(trailingOnly=TRUE)
if (interactive()) {
    file = "gmma_fit_global.rda"
    # file = "gmma_fit_global_assigned.rda"
} else if (length(args) < 1) {
    print("")
    print("usage: Rscript gmma06_analysis.r <gmma_fit_global.rda>")
    quit(save="no")
} else {
    file = args[1]
}
print(sprintf("Read %s",file))
load(file)

# Wild-type stability
cp$dG_wt_init = cp$dG_wt
if ('dG_wt' %in% names(global_fit$par)) cp$dG_wt = global_fit$par['dG_wt']

# Maximum fluorescence base line
cp$B_max_init = cp$B_max
if ('B_max' %in% names(global_fit$par)) cp$B_max = global_fit$par['B_max']

# Minimum fluorescence base line
cp$B_D_init = cp$B_D
if ('B_D' %in% names(global_fit$par)) cp$B_D = global_fit$par['B_D']

print(sprintf("Base line initial fit: dG_wt %6.2f, B_max %6.2f, B_D %6.2f",cp$dG_wt_init,cp$B_max_init,cp$B_D_init))
print(sprintf("Base line global fit : dG_wt %6.2f, B_max %6.2f, B_D %6.2f",cp$dG_wt,cp$B_max,cp$B_D))

nmut = length(mutant[,1])
nsubst = length(subst[,1])
nres = length(residue[,1])

print("=== GMMA Settings ===")
print(settings)


################################################################################
# Global fitted ddG in subst and mutant
################################################################################
print(sprintf("Global fit result (%d): %s",global_fit$info,global_fit$message))
global_dof = length(global_fit$fvec)-length(global_fit$par)
print(sprintf("Estimated %d parameters from %d data points (%d degrees of freedom)",length(global_fit$par),length(global_fit$fvec),global_dof))
print(sprintf("Sum-of-squares %.2f in %d iterations",global_fit$deviance,global_fit$niter))

# chi_sq_red = 1/(length(global_fit$fvec)-length(global_fit$par)) * global_fit$deviance/wt$sd^2
chi_sq = global_fit$deviance/wt$sd^2
print(sprintf("Chi-squared: %.3f (using WT uncertainty %.2f)",chi_sq,wt$sd))
chi_sq_red = chi_sq/global_dof
print(sprintf("Reduced chi-squared: %.3f",chi_sq_red))

subst$ddG_glob = NA
mutant$dG_glob = NA
mutant$residual = NA
si = which(subst$gmma == "use")
subst[si,'ddG_glob'] = global_fit$par[rownames(subst[si,])]
mi = which(mutant$gmma == "use")
mutant[mi,'dG_glob'] = sapply(mut_subst_indices[mi], function(l) {sum(subst[l,'ddG_glob'], na.rm=F)}) + cp$dG_wt
res = global_fit$fvec
names(res) = names(mut_ddG_indices)
mutant[mi,'residual'] = res[rownames(mutant[mi,])]
wt$residual = res["WT"]

# I may change this!
B_pred_init = function(dG) { e = exp(dG/settings$RT); (cp$B_max_init + cp$B_D_init*e) /(1.0+e) }
B_pred_glob = function(dG) { e = exp(dG/settings$RT); (cp$B_max + cp$B_D*e) /(1.0+e) }

# # Compare initial and global fit
# quartz(width=10, height=4)
# par(mfcol=c(1,2))
# x=seq(-20,20,.1)
# # if (! is.null(mutant$dG_init))  {
#     plot(mutant$dG_init, mutant$signal, pch=".", xlim=c(-8,12),
#          main=sprintf("Initial fit, dG_wt %.2f",cp$dG_wt_init), xlab="Mutant stability [kcal/mol]", ylab="Brightness")
# # } else if (! is.null(mutant$ddG_init))  {
# #     plot(cp$dG_wt_init+mutant$ddG_init, mutant$signal, pch=".", xlim=c(-8,12),
# #          main=sprintf("Initial fit, dG_wt %.2f",cp$dG_wt_init), xlab="Mutant stability [kcal/mol]", ylab="Brightness")
# # } else {
# #     plot(0,0,main="Missing mutant$ddG_init or mutant$dG_init")
# # }
# lines(x, B_pred_init(x), col="red", lwd=2)

quartz(width=6, height=4)
x=seq(-20,30,.1)
# plot(mutant$dG_glob, mutant$signal, pch=".", xlim=c(-8,12), main=sprintf("Global fit, dG_wt %.2f",cp$dG_wt), xlab="Mutant stability [kcal/mol]", ylab="Brightness")
# plot(mutant$dG_glob, mutant$signal, pch=".", main=sprintf("Global fit, dG_wt %.2f",cp$dG_wt), xlab="Mutant stability [kcal/mol]", ylab="Brightness")
plot(mutant$dG_glob, jitter(mutant$signal, factor = 0.5), pch=".", main=sprintf("Global fit, dG_wt %.2f",cp$dG_wt), xlab="Variant stability [kcal/mol]", ylab="Signal", yaxt = "n")
axis(side=2, at=c(0,1)) # Signal only takes value 0 or 1
lines(x, B_pred_glob(x), col="red", lwd=2)

quartz.save("global_fit.png",type="png")

# Insert with dG_N vs N to show extrapolation to N=0 and slope <ddG>

# get_si = function(tag, col="i", verbose=F) {
#     # put col="rank" to get subst_sig indices
#     ri = which(residue[,tag] != "-")
#     if (residue[117,tag] == "L") ri = ri[ri!=117]
#     if (residue[118,tag] == "V") ri = ri[ri!=118]
#     # sn = paste(residue[ri,'wt'], rownames(residue)[ri], residue[ri,tag], sep="")
#     sn = unname(unlist(mapply(function(wt,i,s) { paste(wt,i,unlist(strsplit(s,"")),sep="") }, residue[ri,'wt'], rownames(residue)[ri], residue[ri,tag])))
#     si = which(substr(rownames(subst),1,nchar(rownames(subst))-settings$taa_letters+1) %in% sn)
#     # si = subst[sn,col]
#     # if (verbose) {
#     #     print(sprintf("Found %d subst for tag %s with %d missing:",length(sn),tag,sum(is.na(si))))
#     #     print(paste(sn, collapse=", "))
#     # }
#     # si[!is.na(si)]
#     return(si)
# }

# # Plot ddG global vs init correlation w marked known mutations
# quartz(width=6, height=6)
# plot(subst$ddG_glob, subst$init_ddG, pch='.')
# points(subst[get_si('pdb'),'ddG_glob'], subst[get_si('pdb'),'init_ddG'], pch='x', col=5)
# points(subst[get_si('do'),'ddG_glob'], subst[get_si('do'),'init_ddG'], pch='x', col=3)
# points(subst[get_si('now'),'ddG_glob'], subst[get_si('now'),'init_ddG'], pch='x', col=4)
# points(subst[get_si('ts'),'ddG_glob'], subst[get_si('ts'),'init_ddG'], pch='x', col=6)
# points(subst[get_si('sf'),'ddG_glob'], subst[get_si('sf'),'init_ddG'], pch='x', col=2)
# legend("topleft", c("SF","TS","DoBox","Now","PDB"), pch='x', col=c(2,6,3,4,5))
# quartz.save("ddG_global_vs_init.png", type="png")

################################################################################
# Uncertainties from global fit - NEW VERSION
################################################################################
# See minpack.lm:::summary.nls.lm for standard error calculation (methods(class=class(global_fit)) or methods("summary"))
# Condition number of a matrix is the ratio of the largest singular value to the smallest singular value. A matrix is ill-conditioned if the condition number is very high,
#   usually indicating that (i) the lowest singular value is orders of magnitude smaller than the highest one, and (ii) columns/rows of the matrix are heavily correlated
#   with each other leading to redundancies and a matrix that is pretending to be of a higher rank than it truly is.
# Hessian encodes the second derivatives of a function with respect to all pairs of variables. So, if there are n inputs to a function, the gradient is n-dimensional and
#   the Hessian is nxn-dimensional. In machine learning, the inputs are usually the features and the function is usually a loss function we are trying to minimize. When
#   the Hessian is ill-conditioned, it means the basins of the loss function have contours that are very long “ellipsoids” rather than being close to “circular”.
h = global_fit$hessian
print(sprintf("Reciprocal condition number of Hessian: %.2g",rcond(h)))
print("If this is very low, consider reparamerizing some variables to avoid covarying parameters")

# print("Calculating correlations in Hessian matrix")
# hcor = cor(h,t(h))
# # How does parameters correlate with the global parameter
# df = cbind(subst$obs, hcor["dG_wt",][rownames(subst)])
# rownames(df) = rownames(subst) # fill in gaps
# df = df[order(df[,2]),]
# # the reference stability correlates more with substitutions with many observations ... obviously ...

# The covariance matrix may be approximated (I hope!) by the inverse of the Hessian
ih = tryCatch(chol2inv(chol(h)), error = function(e) {print("Cannot calculate inverse of full hessian"); print(e); NA})
# if (all(is.na(ih))) {
#     # h =	...do something clever...
#     h[h<1e-12] = h[h<1e-12] + 1e-12
#     ih = tryCatch(chol2inv(chol(h)), error = function(e) {print("Cannot calculate inverse of tinkered hessian"); print(e); NA})
# }

# Parameter derivatives, i.e. d-param/d-residual, in a list that match the Hessian 
pdh = sqrt(diag(ih))
names(pdh) = names(global_fit$par)
subst$se = pdh[rownames(subst)]
# Matrix inversion may result in very large numbers
subst[which(subst$se > 100),'se'] = Inf

# print("Calculating standard errors of WT fit")
# stderr_wt = sqrt(abs(diag(solve(fit_wt$hessian))))
# n_glob_par = 0
# se_glob = list()
# for (tag in c('dG_wt','B_max','B_D')) {
#     if (tag %in% names(pdh)) {
#         se_glob[tag] = pdh[tag]
#         n_glob_par = n_glob_par+1
# 	print(sprintf("Std error of %s estimated from global fit",tag))
#     } else if (tag %in% names(stderr_wt)) {
#         se_glob[tag] = stderr_wt[tag]
# 	print(sprintf("Std error of %s estimated from WT fit",tag))
#     } else {
#         print(sprintf("Cannot estimate error for %s",tag))
#     }
# }
# print("Std. error of global parameters")
# print(se_glob)

n_glob_par = 0
se_glob = list()
for (tag in c('dG_wt','B_max','B_D')) {
    if (tag %in% names(pdh)) {
        se_glob[tag] = pdh[tag]
        n_glob_par = n_glob_par+1
    }
}
print("Std. error of global parameters")
print(se_glob)

# Remove global parameters and check
if (n_glob_par > 0) pdh = pdh[-seq(length(pdh)+1-n_glob_par,length(pdh))]
i_pdh = which(sapply(names(pdh), function(n) {! n %in% rownames(subst)}))
if (length(i_pdh) > 0) { print("WARNING: Globally fitted parameter(s) not in subst nor recognised as global parameters:"); print(names(pdh)[i_pdh]) }

si_pds  = which(sapply(rownames(subst), function(n) {! n %in% names(pdh)}))
print(sprintf("Substitutions not in global fit: %d",length(si_pds)))
# print(rownames(subst)[si_pds])

si_use = which(subst$gmma=="use")
mi_use = which(mutant$gmma=="use")
obs_use = sapply(si_use, function(si) { length(intersect(subst_mut_indices[[si]], mi_use)) })
obserr = 1/sqrt(obs_use)

# make sure that hessian rows match subst rows before assigning new columns in subst
stopifnot(all( rownames(subst)[si_use] == names(obserr) ))

# The uncertainty based to the variantion of the measurement of the reference
# Includes a factor 3 to get the 99.7% percentile (2 would be 95%)
# The effect categories (subst$eff) depends on this with neutral being zero +- stderrmeas
subst$stderr_meas = subst$se * wt$sd *3 

# Uncertainty of estimation
subst$stderr_meas_est = NA
subst[si_use,"stderr_meas_est"] = subst$stderr_meas[si_use] * obserr

# The uncertainty based on the fit of the model
subst$stderr_fit = subst$se * sqrt(global_fit$deviance/global_dof)
# Uncertainty of estimation
subst$stderr_fit_est = NA
subst[si_use,"stderr_fit_est"] = subst$stderr_fit[si_use] * obserr

res_sq = global_fit$fvec^2
names(res_sq) = names(mut_ddG_indices)

dpf = 1 - length(si_use)/length(mi_use)
print(sprintf("data-to-parameters factor: %.3f (should be close to one but smaller)",dpf))

# Instead of scaling the curvature with the global average residual, scale with the average residual of the points that are
#   directly coupled to the parameter: 'subset error'
calc_suberr = function(si) {
    mi = intersect(subst_mut_indices[[si]], mi_use)
    mn = rownames(mutant)[mi]
    # The 2 is taken out of thin air, one is the fitted ddG and dpf accounts for the parameters partially estimated from this data
    #   incl. the reference stability. Perhaps one is more correct but it favors estimation from only 2 data points which I find to be an overfit
    # dof = length(mn)-2
    if (length(mn) <= 2) return(NA)
    dof = length(mn)*dpf-2
    stopifnot(subst[si,"obs"] > dof)
    # return( sqrt(sum(res_sq[mn]) /dof /chi_sq_red) )
    # Relative mean residual: if 1 the subset of data has the same residuals as the entire model
    return( sqrt( sum(res_sq[mn])/dof / (global_fit$deviance/global_dof) ) )
}

print("Calculate data subset uncertainties")
subst$suberr = NA
subst$suberr[si_use] = sapply(si_use, calc_suberr)

# Uncertaintyfrom number of observations. This should compensate for the case where a suberr is low by chance simply because there are few points to fit.
# E.g. a suberr may be ~0.5 whereas obserr is 0.5 for n=4 and 0.25 for n=16, i.e. suberr*obserr is 0.25 for a lucky fit of 4 points and 0.25 for a normal fit of 16 points
# 4 is arbitrary because I believe that obs=16 is where this factor should be one (could also be 5...). This factor is not in the '_est' unceratinties
subst$obserr = NA
subst$obserr[si_use] = 4.0/sqrt(obs_use)

# An uncertainty to indicate if the ddG_glob value is well estimated
subst$stderr_subfit_est = subst$stderr_fit_est * subst$suberr

################################################################################
# Calculate ddG for hanging nodes 
################################################################################
# NOTE: substitutions without error estimate are discarded later!!

print(sprintf("Estimating stability and std errors for %d substitutions with gmma=='hanging'",sum(subst$gmma=='hanging')))
for (mi in which(mutant$gmma == "hanging")) {
    si = mut_subst_indices[[mi]]
    # which does not have ddG assigned
    sii = which(is.na(subst[si,'ddG_glob']))
    nest = length(sii)
    if (nest < 1) {
        # Subst in hanging variants may in rare cases be in more than one
	print(sprintf("WARNING: No substitutions in hanging variant %s needs to estimate hanging ddG",rownames(mutant)[mi]))
        next
    }
    # ddG sum of those not assigned
    ddG_other = sum(subst[si,'ddG_glob'], na.rm=T)
    # If there's a NA in the error of the ddG's the error of the calculated ddG is also NA (remove hanging subst error which is always NA)
    # max_se = max(subst[si,'stderr_meas'], na.rm=F)
    # max_se = max(subst[si[which(!is.na(subst[si,'ddG_glob']))],'stderr_meas'], na.rm=F)
    # sum_var = sqrt(sum(subst[si[which(!is.na(subst[si,'ddG_glob']))],'se']^2)) * wt$sd * 3
    sum_var = sqrt(sum(subst[si[which(!is.na(subst[si,'ddG_glob']))],'stderr_meas']^2))
    # Only assign non-zero ddG if activity and the sign of dG does not match
    if (sign(0.5-mutant[mi,'active']) == sign(cp$dG_wt+ddG_other)) {
        ddG_missing = 0.0
    # Only calculate ddG in switching region
    } else if (mutant[mi,'signal'] <= cp$B_D+3*wt$sd) {
        ddG_deactivating = -2*cp$dG_wt  # Say it takes twice the stability of the protein to kill it properly 
        ddG_missing = max(ddG_deactivating - ddG_other, 0.0) # if others destabilize enough for deactivation assume no effect
    # } else if (mutant[mi,'signal'] > cp$B_max) {
    } else if (mutant[mi,'signal'] > cp$B_max-3*wt$sd) {
        ddG_missing = min(-1.0*ddG_other, 0.0)               # if others substitutions are stabilizing assume no effect
    } else {
        ddG_missing = settings$RT*log((cp$B_max-mutant[mi,'signal'])/(mutant[mi,'signal']-cp$B_D)) - cp$dG_wt - ddG_other
    }
    
    # If more substitutions are unknown, simply split missing ddG among them. In principle the two are reparametrized to one effect that cannot be split.
    subst[si[sii],'ddG_glob'] = ddG_missing/nest
    # trim to allowed ddG region
    if (ddG_missing/nest > settings$glob_ddG_max) {
        subst[si[sii],'ddG_glob'] = settings$glob_ddG_max
    } else if (ddG_missing/nest < settings$glob_ddG_min) {
        subst[si[sii],'ddG_glob'] = settings$glob_ddG_min
    }
    
    # subst[si[sii],'stderr_meas'] = max_se
    subst[si[sii],'stderr_meas'] = sum_var
    mutant[mi,"dG_glob"] = cp$dG_wt + sum(subst[si,'ddG_glob']) 
    mutant[mi,"residual"] = mutant[mi,"signal"] - B_pred_glob(mutant[mi,"dG_glob"])
    print(sprintf("Added ddG = %6.2f pm %.2f (original %.2f) to subst %4d %5s based on mutant %5d %s of brightness %.2f, dG_wt+ddG_other = %.2f, dG = %.2f and residual %.2f",
        subst[si[sii],'ddG_glob'], sum_var, ddG_missing/nest, si[sii], rownames(subst)[si[sii]], mi, paste(mut_list[[mi]], collapse=","),
	mutant[mi,'signal'], cp$dG_wt+ddG_other, mutant[mi,"dG_glob"], mutant[mi,"residual"]))
	
    # Report if more subst for this mutant was not determined in global fit (independent of how many was estimated here)
    if (sum(subst[si,'gmma']=='hanging') > 1) {
        i = which(subst[si,'gmma']=='hanging')
        print(sprintf("WARNING: Cluster of %d coupled substitutions outside global fit: %s",length(i),paste(rownames(subst[si[i],]),collapse=" ")))
    }
    # stopifnot(abs(subst[si[sii],'ddG_glob']) < 1e-9 | sign(mutant[mi,"dG_glob"]) == sign(0.5-mutant[mi,'active']) )
}
print(sprintf("Total estimated ddG's %d out of %d. Total missing %d", sum(! is.na(subst$ddG_glob)), nsubst, sum(is.na(subst$ddG_glob))))
# print(sprintf("Total estimated std. errors %d out of %d. Total missing %d", sum(! is.na(subst$se)), nsubst, sum(is.na(subst$se))))
print(sprintf("Total estimated std. errors %d out of %d. Total missing %d", sum(! is.na(subst$stderr_meas)), nsubst, sum(is.na(subst$stderr_meas))))


################################################################################
# Substitutions with reliable ddG fit
################################################################################
# Substitutions for which we believe the fit
# se cut between 1 and 2 makes a difference
min_obs = 20
max_se = 2.0
# # original selection
# i_sig = which(!is.na(subst$se) & subst$se < max_se & subst$obs > min_obs)

# Brings discards from 962 to 886 with max_see=0.2
max_see = 0.2
# i_sig = which(!is.na(subst$se) & subst$se/sqrt(subst$obs) < max_see)
# i_sig = which(!is.na(subst$stderr) & subst$stderr*subst$suberr*subst$obserr < max_see)

# Value 0.05 based on manual inspection, e.g. max_see=0.1 gives L5M and F221C as rank 1 and 3 which I do not believe.
# max_sfe = 0.035 reproduces Emma & Philips set1 with 990 subst_sig, 0.04 the GMMA top10 with 1034 (not fluo. top10), and 0.05 1107 subst_sig
max_sfe = 0.05
i_sig = which(!is.na(subst$stderr_meas_est) & subst$stderr_subfit_est < max_sfe)

subst_sig = subst[i_sig,]
nsubst_sig = length(subst_sig$i)
subst_sig = subst_sig[order(subst_sig$ddG_glob),]
# subst_sig = subst_sig[order(subst_sig$ddG_glob+subst_sig$stderr_meas),]
subst$rank = rep(NA,nsubst)
if (nsubst_sig > 0) {
    subst_sig$rank = seq(nsubst_sig)
    subst[rownames(subst_sig),'rank'] = subst_sig$rank
}

quartz(height=6, width=6)
plot(subst$se, 1/sqrt(subst$obs), col=(subst$ddG_glob<0)+1, xlim=c(0,4), ylim=c(0,.6), pch=20, main="Significant substitutions")
abline(v=max_se, h=1/sqrt(min_obs), col=8)
x=seq(0,4,.1)
lines(x,max_see/x, col=8)
points(subst_sig$se, 1/sqrt(subst_sig$obs), col=3, pch=1)
legend("topright", c("ddG_glob < 0","ddG_glob > 0","subst_sig","decision bounds"), col=c(2,1,3,8), pch=c(20,20,1,NA), lty=c(NA,NA,NA,1), bg="white")
quartz.save("err_obs.png",type="png")

print(sprintf("Discarding %d substitutions (out of %d) with high uncertainties",nsubst-length(subst_sig$i),nsubst))
print("Top25 substitutions")
print(head(subst_sig, n=25))

# si = union(get_si('sf'),get_si('do'))
# print(sprintf("Substitutions from SuperFolder and Do & Boxer GFP (%d in Top25, %d in Top50)",sum(subst[si,'rank']<=25, na.rm=T),sum(subst[si,'rank']<=50, na.rm=T)))
# print(subst[si,])


################################################################################
# Substitution effect
################################################################################
# Stability effect of each substitution: unknown, neutral, stab, destab, (very_destab, ifs)
# To be used in the position categories
stability_effect = rep("neu",nsubst)
stability_effect[which(is.na(subst$rank))] = "unknown"  # also covers ddG_glob==NA

# Having fixed boundaries on the neutral category is similar to defining a flat-rate uncertainty that we believe more than other uncertainties
# stability_effect[which(subst$ddG_glob < -.2 & ! is.na(subst$rank))] = "stab"
# stability_effect[which(subst$ddG_glob > .2 & ! is.na(subst$rank))] = "destab"
stability_effect[which(subst$ddG_glob+subst$stderr_meas < 0 & ! is.na(subst$rank))] = "stab"
stability_effect[which(subst$ddG_glob-subst$stderr_meas > 0 & ! is.na(subst$rank))] = "destab"

# stability_effect[which(subst$ddG_glob-subst$stderr_meas > abs(cp$dG_wt) & ! is.na(subst$rank))] = "very_destab"
# Overwrite some with rank==NA
# stability_effect[which(subst$gmma=="hanging" & subst$ddG_glob < -1.0)] = "stab"
# stability_effect[which(subst$gmma=="hanging" & subst$ddG_glob > 1.0)] = "destab"
stability_effect[which(subst$gmma=="hanging" & subst$ddG_glob+subst$stderr_meas < 0.0)] = "stab"
stability_effect[which(subst$gmma=="hanging" & subst$ddG_glob-subst$stderr_meas > 0.0)] = "destab"
stability_effect[which(subst$ddG_glob > abs(cp$dG_wt))] = "destab"  # Also covers subst with Inf or NA errors - consider limit 1.5 or abs(cp$dG_wt)
# stability_effect[which(sapply(subst$note, grep, pattern="ifs")==1)] = "ifs"

subst$eff = factor(stability_effect, levels=c("stab","destab","neu","unknown"))
subst_sig$eff = subst[rownames(subst_sig),"eff"]

# Stability effect categories
ssi_stab = which(subst_sig$eff == "stab")
ssi_neut = which(subst_sig$eff == "neu")
ssi_dest = which(subst_sig$eff == "destab")
# check that all subst_sig should be in the above 3 categories
ssi = c(ssi_stab,ssi_neut,ssi_dest)
stopifnot(ssi[order(ssi)] == seq(nsubst_sig))
# stopifnot(ssi == seq(nsubst_sig))
print(sprintf("Of %d sig. subst., %d (%.1f%%)  %d (%.1f%%)  %d (%.1f%%)  are stabilizing, neutral and destabilizing respectively.",
    nsubst_sig, length(ssi_stab),length(ssi_stab)*100.0/nsubst_sig, length(ssi_neut),length(ssi_neut)*100.0/nsubst_sig, length(ssi_dest),length(ssi_dest)*100.0/nsubst_sig))


################################################################################
# Plot fit
################################################################################
print("Method of ddG_glob vs. method of ddG_init")
print(xtabs(~ subst$gmma + subst$init_m))
print("Substantial fit vs. method of ddG_init")
sig = c("sig","un-sig")[is.na(subst$rank)+1]
print(xtabs(~ sig + subst$init_m))



















# from make_synth_binary
###############################################################

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

#######################################################################################################################################
#######################################################################################################################################
##
##  SUBSTITUTION GENERATION
##
#######################################################################################################################################
#######################################################################################################################################
subst_param = list()
extra = list()


# mock scan_var to save individual run outptus plots
# extra$scan_var = round(scan_var)
# n = extra$scan_var
# print(sprintf(" RESULT_SCAN : %d mock scan ", n))




# gfp_all_single = all_single(gfp_ref)

# List of substitutions
# subst_list = c("A1d", "A1e", "B2d", "B2e", "B2f", "C4c") # test list
# subst_list = gfp_all_single
# subst_list = sample(gfp_all_single, 1000, replace=F)

# # Sarkisyan substitutions
# subst_param$subst_mask = subst$gmma=="use"
# subst_list = rownames(subst[which(subst_param$subst_mask),])

# # Library with destabilizing pair-combinatorial backgrounds, aka s77
# stab_subst = rownames(subst_sig[1:50,])
# taa = substr(gfp_all_single,nchar(gfp_all_single),nchar(gfp_all_single)) 
# wtaa = substr(gfp_all_single,1,1) 
# iI2V = which(wtaa=="I" & taa=="V")
# iV2A = which(wtaa=="V" & taa=="A")
# iV2A = iV2A[iV2A != which(gfp_all_single == "V161A")]
# destab_subst = gfp_all_single[c(iI2V,iV2A)]
# destab_subst = destab_subst[! destab_subst %in% stab_subst]
# subst_list = c(stab_subst, destab_subst)
# # store for later usage
# extra[["stab_subst"]] = stab_subst
# extra[["destab_subst"]] = destab_subst



# Trimodal synthetic substitution effects, aka s75
subst_param$wt = laub_wt
subst_param$n_stab = 0
subst_param$n_neu = 19
subst_param$n_destab = 57
subst_param$n_ifs = 0
subst_list = sample(all_single(laub_wt), subst_param$n_stab+subst_param$n_neu+subst_param$n_destab+subst_param$n_ifs, replace=FALSE)
print(sprintf("Generated %d substotutions of Laub with %d stabilizing, %d neutral and %d destabilizing, and %d IFS",
              length(subst_list),subst_param$n_stab,subst_param$n_neu,subst_param$n_destab, subst_param$n_ifs))



# 
# # scan trimodal distribution of n_stab,n_neu,n_destab
# extra$scan_var = round(scan_var)
# n = extra$scan_var
# print(sprintf(" RESULT_SCAN : %d IFSs", n))
# 
# 
# # Trimodal synthetic substitution effects, aka s75
# subst_param$wt = laub_wt
# subst_param$n_stab = 0
# subst_param$n_neu = 19 
# subst_param$n_destab = 57 - n
# subst_param$n_ifs = n
# subst_list = sample(all_single(laub_wt), subst_param$n_stab+subst_param$n_neu+subst_param$n_destab+subst_param$n_ifs, replace=FALSE)
# print(sprintf("Generated %d substotutions of Laub with %d stabilizing, %d neutral and %d destabilizing, and %d IFS",
#               length(subst_list),subst_param$n_stab,subst_param$n_neu,subst_param$n_destab, subst_param$n_ifs))


# Sort substitutions according to position and new amino acid
stopifnot(length(subst_list) == length(unique(subst_list)))
subst_list = subst_list[order( as.numeric(substr(subst_list,2,nchar(subst_list)-1)), substr(subst_list,nchar(subst_list),nchar(subst_list)))]

save(subst_list, subst_param, extra, file="synth_subst.rda")
# load(paste0(ref_dir,"/synth_subst.rda"))

#############################################
# Generate or calculate substitution effects (ddG's)
#############################################
ddG_param = list()
ddG = rep(NA, length(subst_list))
names(ddG) = subst_list

# # Sarkisyan-like data
# ddG[] = subst[which(subst_param$subst_mask),"ddG_glob"]
# # replace uncertain ddG's with random values drawn from a normal approximation of GFP ddG's
# subst_uncert = rownames(subst[which(subst_param$subst_mask & is.na(subst$rank)),])
# ddG[subst_uncert] = rnorm(length(subst_uncert), mean=mean(subst_sig$ddG_glob), sd=sd(subst_sig$ddG_glob))


# # scan ddG sd
# extra$scan_var = as.numeric(scan_var)
# n = extra$scan_var
# print(sprintf(" RESULT_SCAN : ddG sd  : %f  ", n))

# Trimodal synthetic substitution effects, aka s75
ddG_param$mean_stab = -1.0
ddG_param$mean_neu = 0.0
ddG_param$mean_destab = 2.0
ddG_param$mean_ifs = 10.0
ddG_param$sd = 0.4
ddG_param$cat = sample(c(rep(1,subst_param$n_stab), rep(2,subst_param$n_neu), rep(3,subst_param$n_destab), rep(4,subst_param$n_ifs)))
ddG[which(ddG_param$cat==1)] = rnorm(subst_param$n_stab, ddG_param$mean_stab, ddG_param$sd)
ddG[which(ddG_param$cat==2)] = rnorm(subst_param$n_neu, ddG_param$mean_neu, ddG_param$sd)
ddG[which(ddG_param$cat==3)] = rnorm(subst_param$n_destab, ddG_param$mean_destab, ddG_param$sd)
ddG[which(ddG_param$cat==4)] = rnorm(subst_param$n_ifs, ddG_param$mean_ifs, ddG_param$sd)
print(sprintf("Generated %d ddGs with observed mean %.3f and std %.3f",length(ddG),mean(ddG),sd(ddG)))




#############################################################################################################























# Histogram of ddG_glob
# quartz(height=5, width=6)

hist(subst$ddG_glob, breaks=seq(settings$glob_ddG_min,settings$glob_ddG_max,.2), xlim=c(settings$glob_ddG_min,settings$glob_ddG_max), col="gray90", main="Global ddG fit", xlab = "ddGs [kcal/mol]")
h=hist(subst_sig$ddG_glob, breaks=seq(settings$glob_ddG_min,settings$glob_ddG_max,.2), plot=F)
plot(h, add=T, col="gray70")
legend("topright", c("Global fit","Low uncertainty"), col=c("gray90","gray70"), pch=15)
lines(density(ddG))

# xaxis()
# quartz.save("ddG_hist.png",type="png")

quartz(height=5, width=6)
plot(density(ddG), xlim=(c(-2,10)), col = "red", main = "ddGs values Laub real and synthetic data", xlab = "ddGs [kcal/mol]")
lines(density(subst$ddG_glob), xlim=(c(-2,10)), col = "blue")
legend("topright", c("Synthetic data","Laub data"), col = c("red", "blue"), pch = 1)

quartz.save("laub_real_synth_ddGs_dist.png",type="png")





# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ################################################################################
# # Position categories
# ################################################################################
# residue$fit = rep("",nres)    # Which observed substitutions (residue$subst) could be fitted reliably
# residue$ifs = rep("",nres)    # Which observed substitutions (residue$subst) are IFS
# residue$opt = rep('-',nres)   # Optimal identity from GMMA
# # residue$opt_sig = rep('-',nres)   # Optimal identity from GMMA
# residue$opt0 = rep('-', nres) # Optimal identity from intial fit
# residue$top = rep(NA,nres)    # Min rank among substitutions
# residue$n_max = rep(0,nres)   # Max number of substitutions among active variants
# # residue$cat_old = rep(NA,nres)    # Category of position
# residue$cat = rep(NA,nres)    # Category of position
# 
# # pos_eff = t(sapply(seq(nres), function(resi) {table(subst[which(subst$resi==resi),'eff'])}))
# 
# # Notes to implement
# # If a position has only inactive and subst with >100 obs (needed to detect ifs by statistical means) are ifs the position may be strictly conserved
# 
# # Assign residue cathegory
# si_ifs = which(sapply(subst$note, function(s) {"ifs" %in% strsplit(s," ")[[1]]} ))
# for (resi in seq(nres)) {
#     # Assignments
#     si = which(subst$resi==resi)
#     residue[resi,'fit'] = paste(subst[intersect(si,which(!is.na(subst$rank))),'taa'], collapse="")
#     residue[resi,'ifs'] = paste(subst[intersect(si,si_ifs),'taa'], collapse="")
#     if (any(!is.na(subst[si,'rank']))) residue[resi,'top'] = min(subst[si,'rank'], na.rm=T)
# 
#     mi = unique(unlist(subst_mut_indices[si]))
#     # i_active = which(mutant[mi,'signal'] > 2.5)
#     mii_active = which(mutant[mi,'active'] > 0.5)
#     residue[resi,'n_max'] = max(c(0, mutant[mi[mii_active],'N_sub']), na.rm=T)
#     si = which(subst$resi==resi)
#     if (length(si) > 0) {
#         if (min(subst[si,'init_ddG'], na.rm=T) < 0.0) residue[resi,'opt0'] = subst[ si[which.min(subst[si,'init_ddG'])], 'taa' ]
#         if (any(subst[si,'eff'] == "stab")) {
# 	    si_stab = si[which(subst[si,'eff']=="stab")]
# 	    residue[resi,'opt'] = subst[ si_stab[which.min(subst[si_stab,'ddG_glob'])], 'taa' ]
# 	}
#     }
# 
#     t = table(subst[si,"eff"])
#     t_known_sum = sum(t[c("stab","destab","neu")])
#     if      (t["stab"]   > 1)                              residue[resi,'cat'] = 9  # more than one stabilizing
#     else if (t["stab"]   > 0)                              residue[resi,'cat'] = 8  # any stabilizing
#     else if (t["neu"]    > 0 & t["neu"]    >  t["destab"]) residue[resi,'cat'] = 6  # tolerant position with most neutral
#     else if (t["destab"] > 0 & t["destab"] >= t["neu"])    residue[resi,'cat'] = 0  # stable position with most destabilizing
#     else if (t["unknown"] == sum(t))                       residue[resi,'cat'] = 5  # this covers sum(t)==0
#     else                                                   residue[resi,'cat'] = 4  # what didn't I think about?
#     
# }
# 
# #     wt           subst na_mut nd_mut av sf do now      pdb opt opt0 ifs     fit n_max cat top
# # 1    K         *EMNQRT    734    595  -  -  Q   -        A   -    -       EMQRT     9   2 144
# f = file("residue.txt", "wt")
# write('"# Global Multi-Mutant Analysis per position"',f)
# write('"# wt: Wild-type (eGFP) amino acid at position"',f)
# write('"# subst: Attempted substitutions at position"',f)
# write('"# active: Number of active variants with this position mutated"',f)
# write('"# inactive: Number of inactive variants with this position mutated"',f)
# write('"# av: Amino acid in Aequorea victoria (even-wilder type)"',f)
# write('"# sf,do,now: Amino acid in SuperFolderGFP, Do & Boxer split-GFP and Sarkisyan & Kondrashov NowGFP"',f)
# write('"# pdb: Amino acid observed in homologs from the PDB"',f)
# write('"# opt: Optimal amino acid according to GMMA"',f)
# write('"# opt0: Optimal amino acid according to initial fit"',f)
# write('"# ifs: Amino acids (among attempted) that inactivates irreversibly at this position "',f)
# write('"# fit: Amino acids (among attempted) that resulted in an reliable GMMA fit"',f)
# write('"# top: Min rank among substitutions"',f)
# write('"# n_max: Max number of substitutions among active variants"',f)
# write('"# cat: Mutational category of position:"',f)
# write('"#     9   Hot-spot: All stabilizing, no IFS"',f)
# write('"#     8   Half or more stabilizing, no IFS"',f)
# write('"#     7   Less than half stabilizing and possible some IFS but most variant are active"',f)
# write('"#     6   At lest one stabilizing substitution"',f)
# write('"#     5   Nothing is known"',f)
# write('"#     4   No fitted substitutions but some IFS"',f)
# write('"#     3   Nothing is known but few active"',f)
# write('"#     2   Less than half are very destabilizing and some variants are active"',f)
# write('"#     1   Not all very destabilizing"',f)
# write('"#     0   Conserved: All fitted substitutions are very detabilizing"',f)
# residue_named = cbind(resi=seq(nres), residue)
# 
# # # clean fields without data
# # ri = which(residue_named$wt=="x")
# # residue_named = residue_named[-ri,]
# # print(sprintf("Removed %d rows from residue data frame without substitutions: %s",length(ri),paste(ri,collapse=" ")))
# 
# write.table(residue_named,f)
# close(f)
# print("Dumped residue.txt")
# 
# for (ri in seq(nres)) {
#     si = which(subst$resi == ri)
#     ssi = which(subst_sig$resi == ri)
#     if (residue[ri,"cat"] == 0) {
#         subst[si,"note"] = paste(subst[si,"note"], "posS", sep=" ")
#         subst_sig[ssi,"note"] = paste(subst_sig[ssi,"note"], "posS", sep=" ")
#     } else if (residue[ri,"cat"] %in% c(7,8,9)) {
#         subst[si,"note"] = paste(subst[si,"note"], "posU", sep=" ")
#         subst_sig[ssi,"note"] = paste(subst_sig[ssi,"note"], "posU", sep=" ")
#     } else if (residue[ri,"cat"] %in% c(6)) {
#         subst[si,"note"] = paste(subst[si,"note"], "posT", sep=" ")
#         subst_sig[ssi,"note"] = paste(subst_sig[ssi,"note"], "posT", sep=" ")
#     }
# }
# 
# # Dump subst and mut_list in excel raedable format
# subst_named = cbind(s=rownames(subst), subst)
# write.table(subst_named, "subst.csv", sep=";", row.names=F)
# write.table(mutant, "mutant.csv", sep=";", row.names=F)
# print("Dumped subst.csv and mutant.csv")
# 
# # print("Table of residue$cat and correlation with burial2")
# # print(table(residue$cat))
# # print(table(residue[,c("burial2","cat")]))
# 
# print("")
# ri = which(residue$cat == 0)
# print("Stabil positions (resi+2)")
# print(paste(residue[ri,"wt"], ri+2, sep="", collapse=", "))
# print("PyMol notation (resi+2):")
# print(paste(ri+2, collapse="+"))
# 
# print("")
# ri = which(residue$cat == 6)
# print("Tolerant positions (resi+2)")
# print(paste(residue[ri,"wt"], ri+2, sep="", collapse=", "))
# print("PyMol notation (resi+2):")
# print(paste(ri+2, collapse="+"))
# 
# print("")
# ri = which(residue$cat %in% c(7,8,9))
# print("Positions with engineering potential (resi+2)")
# print(paste(residue[ri,"wt"], ri+2, sep="", collapse=", "))
# print("PyMol notation (resi+2):")
# print(paste(ri+2, collapse="+"))
# 
# print("")
# print("stabilizing substitution at unstable positions")
# si = which(sapply(subst$note, function(s){length(grep(s,pattern="posU"))>.5}) & subst$eff=="stab")
# si = si[order(subst[si,"ddG_glob"])]
# print(paste(strsplit(wt$seq, "")[[1]][subst[si,"resi"]], subst[si,"resi"]+2, subst[si,"taa"], sep="", collapse=", "))
# print("stabilizing substitution at possible unstable positions")
# si = which(sapply(subst$note, function(s){length(grep(s,pattern="posPU"))>.5}) & subst$eff=="stab")
# si = si[order(subst[si,"ddG_glob"])]
# print(paste(strsplit(wt$seq, "")[[1]][subst[si,"resi"]], subst[si,"resi"]+2, subst[si,"taa"], sep="", collapse=", "))
# 
# print("Effect of sig. substitutions from Gly:")
# si = which(strsplit(wt$seq, "")[[1]][subst_sig$resi] == "G")
# print(table(subst_sig[si,"eff"]))
# 
# print("Efect of sig. substitutions from His:")
# si = which(strsplit(wt$seq, "")[[1]][subst_sig$resi] == "H")
# print(table(subst_sig[si,"eff"]))
# 
# ################################################################################
# # Dump results of post processing
# ################################################################################
# # version=2 is necessary for R versions <3.5 to be able to read it (binf version 3.1)
# save(global_fit, mut_ddG_indices, cp, settings, fit_wt, wt, mut_list,
#      mutant, subst, subst_sig, residue, mut_subst_indices, subst_mut_indices, res_mut_indices,
#      file="gmma_result.rda", version=2)
# print("Post-processing done, dumped gmma_result.rda.")
# 
# 
# # Dump substitution effects
# f = file("prism_gmma_000_GFP_synthetic.txt", "wt")
# write("# --------------------", f)
# write("# version: 1", f)
# write("# protein:", f)
# write("#     name: GFP", f)
# write("#     organism: Aequorea victoria (Jellyfish)", f)
# write(paste("#     sequence:",wt$seq), f)
# write("#     uniprot: P42212", f)
# write("# gmma:", f)
# write(sprintf("#     dG_wt: %.3f",cp$dG_wt), f)
# write("# variants:", f)
# write(paste("#     number:",nrow(subst)), f)
# write("# columns:", f)
# write("#     ddG: Calculated or random stability effects", f)
# write("# --------------------", f)
# write("#", f)
# write.table(data.frame(variant=rownames(subst), ddG=subst$ddG_glob, std=subst$stderr_meas, std_mean=subst$stderr_subfit_est, ddG_init=subst$init_ddG,
#                        obs_act=subst$active, obs_inact=subst$inactive, rank=subst$rank, row.names=NULL), file=f, row.names=F, quote=F)
# close(f)
