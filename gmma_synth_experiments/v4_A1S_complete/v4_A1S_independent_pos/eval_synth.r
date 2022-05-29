options(width=300, digits=4, stringsAsFactors=F)

load("gmma_result.rda")
load("make_synthetic.rda")

gmma = read.table("prism_gmma_000_GFP_synthetic.txt", header=T)
synth = read.table("gmma_synthetic_true.txt", header=T)

common_subst = intersect(gmma$variant, synth$variant)
print(sprintf("GMMA has %d subst., synth. data %d subst. with %d subst. in common",length(gmma$variant), length(synth$variant), length(common_subst)))
stopifnot(length(gmma$variant) == length(synth$variant))

# rownames(synth) = synth$variant
gmma$synth = synth[match(gmma$variant,synth$variant),"ddG"]

rs_all = cor(gmma[,"ddG"], gmma[,"synth"], method="spearman", use="complete.obs")
rp_all = cor(gmma[,"ddG"], gmma[,"synth"], method="pearson", use="complete.obs") 
std_all = sd(gmma[,"ddG"] - gmma[,"synth"], na.rm=T)/sqrt(2)
mae_all = mean(abs(gmma[,"ddG"] - gmma[,"synth"]), na.rm=T)

std_mean_cut = 1.0
i = which(gmma$std_mean < std_mean_cut)

ni = setdiff(seq(nrow(gmma)), i)
rs = cor(gmma[i,"ddG"], gmma[i,"synth"], method="spearman")
rp = cor(gmma[i,"ddG"], gmma[i,"synth"], method="pearson") 
std = sd(gmma[i,"ddG"] - gmma[i,"synth"])/sqrt(2)
mae = mean(abs(gmma[i,"ddG"] - gmma[i,"synth"]))

# The first may be sensitive to high-ddG outliers which are not always expected to be well determined
fit1 = lm(ddG ~ synth, data=gmma)
# fit2 = glm(ddG ~ synth, family=gaussian, data=gmma)
# fit3 = lm(ddG ~ synth, data=gmma[which(gmma$ddG < -cp$dG_wt),])


# save plot for every scan var value run, collect in separate directory
if(!(dir.exists("eval_synth"))){ 
  dir.create("eval_synth")
}
# quartz(width=7, height=7.5)
png(sprintf("eval_synth/eval_synth_%f.png", extra$scan_var), width = 7, height = 7.5, units = "in", res = 100)
# png("eval_synth.png", width = 7, height = 7.5, units = "in", res = 100)
s = sprintf("dG_wt = %.2f (true %.2f)", cp$dG_wt, dG_wt)
plot(gmma[,"synth"], gmma[,"ddG"], ylim=c(min(gmma[,"ddG"], na.rm=T)-1.0, max(gmma[,"ddG"], na.rm=T)+1.0), main=s, xlab="True ddG", ylab="Estimated ddG")
i = which(gmma[,"std"] >0.001)
arrows(x0=gmma[i,"synth"], y0=gmma[i,"ddG"]-gmma[i,"std"], y1=gmma[i,"ddG"]+gmma[i,"std"], length=0.01, angle=90, code=3, lwd=2)
i = which(gmma[,"std_mean"] >0.001)
arrows(x0=gmma[i,"synth"], y0=gmma[i,"ddG"]-gmma[i,"std_mean"], y1=gmma[i,"ddG"]+gmma[i,"std_mean"], length=0.02, angle=90, code=3, lwd=2)
i = which(gmma[,"std"] >0.001)
arrows(x0=gmma[i,"synth"], y0=gmma[i,"ddG"]-gmma[i,"std"], y1=gmma[i,"ddG"]+gmma[i,"std"], length=0.01, angle=90, code=3, lwd=2)
abline(0,1)
# mark high error ddG's
points(gmma[ni,"synth"], gmma[ni,"ddG"], pch=20, col=2)
x = c(min(gmma$synth),max(gmma$synth))
lines(x, x*fit1$coefficients["synth"] + fit1$coefficients["(Intercept)"], col=2, lwd=3)
legend("topleft", c(sprintf("Points std>%.2f: %d (%d, NA's %d)", std_mean_cut, length(i), nrow(gmma), sum(is.na(gmma$ddG))), sprintf("Pearson %.3f (%.3f)",rp,rp_all),
                    sprintf("Spearman %.3f (%.3f)",rs,rs_all), sprintf("MAE %.3f (%.3f)",mae,mae_all),
                    sprintf("Std %.3f (%.3f)",std,std_all)))
legend("bottomright", sprintf("lm all ddG: %.3f*ddG + %.3f",fit1$coefficients["synth"], fit1$coefficients["(Intercept)"]), lwd=c(3,2,1), col=c(2,3,4))

# Regional correlation - mostly for tri-modal synthetic data
region_cor_plot = function(x0, x1, y_offset) {
  i = which(gmma[,"synth"] < x1 & gmma[,"synth"] > x0)
  if(length(i) > 1){
    mx = mean(gmma[i,"synth"])
    my = mean(gmma[i,"ddG"], na.rm=T)
    arrows(x0=x0, x1=x1, y0=my+y_offset, angle=90, length=0.1, code=3)
    correl = cor(gmma[i,"ddG"], gmma[i,"synth"], method="pearson", use="complete.obs") # leads to error for if no values exist in interval.
    text(mx, my+y_offset*1.2, sprintf("Pearson %.2f",correl))
    print(sprintf("Printed correlation af x,y = %.2f, %.2f",mx, my+y_offset*1.2))
    return(correl)    
  }
  else
    return(NA) # return NA if no values exist in interval given
}

# automatic text position y-offset from linear fit1
ddG_mode_sep = 1.0
offest = 0.8 * ddG_mode_sep*fit1$coefficients["synth"] + fit1$coefficients["(Intercept)"]


rp_neustab = region_cor_plot(-10, 0.5, offest) # correlation excluding high ddG outliers

##################
# regional spearman correlation
region_cor_plot_s = function(x0, x1, y_offset) {
  i = which(gmma[,"synth"] < x1 & gmma[,"synth"] > x0)
  if(length(i) > 1){
    mx = mean(gmma[i,"synth"])
    my = mean(gmma[i,"ddG"], na.rm=T)
    arrows(x0=x0, x1=x1, y0=my+y_offset, angle=90, length=0.1, code=3)
    correl = cor(gmma[i,"ddG"], gmma[i,"synth"], method="spearman", use="complete.obs") # leads to error for if no values exist in interval.
    text(mx, my+y_offset*1.2, sprintf("Pearson %.2f",correl))
    print(sprintf("Printed correlation af x,y = %.2f, %.2f",mx, my+y_offset*1.2))
    return(correl)    
  }
  else
    return(NA) # return NA if no values exist in interval given
}


##################




# quartz.save(sprintf("eval_synth/eval_synth%f.png", extra$scan_var), type="png")
dev.off()


# save values but do not print in eval_synth plot
# mock plot to reuse corr function without plotting in main eval_synth plots
plot(cars)
rp_stab = region_cor_plot(-10, -0.5, offest) 
rp_neu = region_cor_plot(-0.5, 0.5, -offest)
rp_destab = region_cor_plot(0.5, 10, -offest)


rs_stab = region_cor_plot_s(-10, -0.5, offest) 
rs_neu = region_cor_plot_s(-0.5, 0.5, -offest)
rs_destab = region_cor_plot_s(0.5, 10, -offest)

rs_neustab = region_cor_plot_s(-10, 0.5, offest) # correlation excluding high ddG outliers


# Plot library composition and frection active
# First, build multi-mutant data frame - could use lib_compo here
n_mut = sapply(strsplit(var, ":"), length)
max_mut = max(n_mut)
nmut_df = data.frame(n_mut=seq(max_mut), n=0, n_act=0, f_act=NA)
t = table(n_mut)
t_act = table(n_mut[which(B >= cp$B_mid)])
for (i in nmut_df$n_mut) {
  if (! is.na(t[as.character(i)]))     nmut_df[i,"n"] = t[as.character(i)]
  if (! is.na(t_act[as.character(i)])) nmut_df[i,"n_act"] = t_act[as.character(i)]
  nmut_df[i,"f_act"] = nmut_df[i,"n_act"] / nmut_df[i,"n"]
}
nlib = sum(nmut_df[,"n"])
frac_act = sum(nmut_df$n_act)/nlib

# # # save plot for every scan var value run, collect in separate directory
# # if(!(dir.exists("library"))){ 
# #     dir.create("library")
# # }
# # # quartz(width=7, height=7.5)
# # png(sprintf("library/library_%f.png", extra$scan_var), width = 5, height = 7, units = "in", res = 100)
# # quartz(height=5, width=7)
# png("library.png", height=5, width=7, units = "in", res = 100)
# s = sprintf("Fraction active %.2f of %d",frac_act,nlib)
# plot(0, 0, col="white", xlim=c(1,max_mut), ylim=c(0,1.3), xaxp=c(1,max_mut,max_mut-1), xlab="N-mutant", ylab="", main=s)
# points(nmut_df$n_mut, nmut_df$n/nlib, type="b", pch=20, lwd=2, col=2)
# points(nmut_df$n_mut, nmut_df$f_act, type="b", pch=20, lwd=2, col=4)
# legend("topleft", c("Composition","Fraction active"), pch=20, col=c(2,4), ncol=2)
# # quartz.save("library.png", type="png")
# dev.off()

# print(sprintf("Spearman, all points: %.3f ",rs_all))
# print(sprintf("Pearson, all points: %.3f ",rp_all))
# print(sprintf("Pearson, low uncertainty points: %.3f ",rp))
# print(sprintf("Pearson, stabilizing: %.3f ",rp_stab))
# print(sprintf("Pearson, neutral: %.3f ",rp_neu))
# print(sprintf("Pearson, destabilizing: %.3f ",rp_destab))
# print(sprintf("Line fit slope %.3f intersect %.3f ",fit1$coefficients["synth"],fit1$coefficients["(Intercept)"]))
# print(sprintf("Fraction active: %.3f ",frac_act))



 
# scan_var output
print(" RESULT_LAB   dG_wt_fit dG_wt_true      B_D    B_max  cor_s  cor_all  cor_low  cor_p_stab  cor_p_neu  cor_p_dest cor_s_stab  cor_s_neu  cor_s_dest cor_p_neustab cor_s_neustab    slope    inter  frac_act  rss scan p1_-1_ddG p1_-1_rank p1_0_ddG p1_0_rank p1_1_ddG p1_1_rank p2_-1_ddG p2_-1_rank p2_0_ddG p2_0_rank p2_1_ddG p2_1_rank p3_-1_ddG p3_-1_rank p3_0_ddG p3_0_rank p3_1_ddG p3_1_rank p4_-1_ddG p4_-1_rank p4_0_ddG p4_0_rank p4_1_ddG p4_1_rank ")
print(sprintf(" RESULT_VAL %7.3f %7.3f  %7.3f  %7.3f  %5.3f    %5.3f    %5.3f     %5.3f    %5.3f     %5.3f  %5.3f    %5.3f     %5.3f  %7.3f  %7.3f  %7.3f  %7.3f   %5.3f  %5.3f %5.3f  %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f",
              cp$dG_wt, dG_wt, cp$B_D, cp$B_max,
              rs_all, rp_all,
              rp, rp_stab, rp_neu, rp_destab,
              rs_stab, rs_neu, rs_destab,
              rp_neustab, rs_neustab,
              fit1$coefficients["synth"], fit1$coefficients["(Intercept)"], frac_act, global_fit$deviance, extra$scan_var,
              
              # fitted ddG value and rank of subst init as -1 (gaussian sampled ddG value very unlikely to be exactly -1) on pos1 (the first because gmma sorted by subst pos)
              gmma[which(gmma$synth == -1), "ddG" ][1],gmma[which(gmma$synth == -1), "rank" ][1],               
              gmma[which(gmma$synth == 0), "ddG" ][1],gmma[which(gmma$synth == 0), "rank" ][1],
              gmma[which(gmma$synth == 1), "ddG" ][1],gmma[which(gmma$synth == 1), "rank" ][1],
              
              # pos 2
              gmma[which(gmma$synth == -1), "ddG" ][2],gmma[which(gmma$synth == -1), "rank" ][2],               
              gmma[which(gmma$synth == 0), "ddG" ][2],gmma[which(gmma$synth == 0), "rank" ][2],
              gmma[which(gmma$synth == 1), "ddG" ][2],gmma[which(gmma$synth == 1), "rank" ][2],
              
              # pos 3 and 4 "free" replicates of posiiton 2
              gmma[which(gmma$synth == -1), "ddG" ][3],gmma[which(gmma$synth == -1), "rank" ][3],               
              gmma[which(gmma$synth == 0), "ddG" ][3],gmma[which(gmma$synth == 0), "rank" ][3],
              gmma[which(gmma$synth == 1), "ddG" ][3],gmma[which(gmma$synth == 1), "rank" ][3],
              
              gmma[which(gmma$synth == -1), "ddG" ][4],gmma[which(gmma$synth == -1), "rank" ][4],               
              gmma[which(gmma$synth == 0), "ddG" ][4],gmma[which(gmma$synth == 0), "rank" ][4],
              gmma[which(gmma$synth == 1), "ddG" ][4],gmma[which(gmma$synth == 1), "rank" ][4]
              
              
              )) 


