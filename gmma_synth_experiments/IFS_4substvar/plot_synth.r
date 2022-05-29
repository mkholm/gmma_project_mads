# load gmma scan processed output, make plots





# plot scan against all other parameters at sime time, split grid? --> quick overview
# 9 params + scan
plot_summary <- function(d){
  {
    par(mfrow = c(4,4))
    
    plot(d$scan, d$dG_wt_fit)
    plot(d$scan, d$dG_wt_true)
    plot(d$scan, d$cor_s)
    plot(d$scan, d$cor_all)
    # plot(d$scan, d$cor_low)
  
    plot(d$scan, d$cor_p_stab)
    plot(d$scan, d$cor_p_neu)
    plot(d$scan, d$cor_p_dest)
    plot(d$scan, d$cor_p_neustab)
    
    
    plot(d$scan, d$cor_s_stab)
    plot(d$scan, d$cor_s_neu)
    plot(d$scan, d$cor_s_dest)
    plot(d$scan, d$cor_s_neustab)
    
    
    plot(d$scan, d$slope)
    plot(d$scan, d$inter)
    plot(d$scan, d$frac_act)
    plot(d$scan, d$RSS)
  
    # follow ddGs value of specific subst
    # plot(d$scan, d$rank_ddGs_fit_stab)
    # plot(d$scan, d$ddGs_fit_stab)
    # 
        
    mtext("synth scan summary",side = 3,
          line = - 2,
          outer = TRUE)
  }
}

d <- read.table("scan_dat_pruned.txt", header = TRUE)

png("scan_summary.png", width = 10, height = 14, units = "in", res = 200)
plot_summary(d)
dev.off()

# make nice plots of interestning correlations
{
png(sprintf("correlations_summary_%s.png", rev(strsplit(getwd(), "/")[[1]])[1]))

par(mfrow = c(2,2))
# correlations
# global correlations
plot(d$scan, d$cor_all,
     xlab = "Number of IFSs", ylab = "r", main = "r, ddGs global")
plot(d$scan, d$cor_s,
     xlab = "Number of IFSs", ylab = "ρ", main = "ρ, ddGs global")

# type = "l" 

# correlations of interest: true stabilising and neutral subst correlations (ddG true [-10;0.5])
plot(d$scan, d$cor_p_neustab, 
     xlab = "Number of IFSs", ylab = "r", main = "r, true ddGs [-10;0.5]")
plot(d$scan, d$cor_s_neustab, 
     xlab = "Number of IFSs", ylab = "ρ", main = "ρ, true ddGs [-10;0.5]")

# mtext("synth scan summary",side = 3,
#       line = - 2,
#       outer = TRUE)

dev.off()
}

png(sprintf("GOF_%s.png", rev(strsplit(getwd(), "/")[[1]])[1]))
# model goodness of fit (as residual sum of squares)
plot(d$scan, d$RSS, type = "l"
     ,xlab = "Number of IFSs", ylab = "Residual sum of squares", main = "Model fit")
dev.off()
