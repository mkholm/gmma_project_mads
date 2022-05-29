# load gmma scan processed output, make plots





# plot scan against all other parameters at sime time, split grid? --> quick overview
# 9 params + scan
# plot_summary <- function(d){
#   {
#     par(mfrow = c(5,6))
#     # par(mfrow = c(3,3))
#     plot(d$scan, d$dG_wt_fit)
#     plot(d$scan, d$dG_wt_true)
#     plot(d$scan, d$cor_s)
#     plot(d$scan, d$cor_all)
#     # plot(d$scan, d$cor_low)
#   
#     plot(d$scan, d$cor_p_stab)
#     plot(d$scan, d$cor_p_neu)
#     plot(d$scan, d$cor_p_dest)
#     
#     plot(d$scan, d$cor_s_stab)
#     plot(d$scan, d$cor_s_neu)
#     plot(d$scan, d$cor_s_dest)
#     
#     
#     plot(d$scan, d$slope)
#     plot(d$scan, d$inter)
#     plot(d$scan, d$frac_act)
#   
#     # follow ddGs value of specific subst
#     plot(d$scan, d$rank_ddGs_fit_stab)
#     plot(d$scan, d$ddGs_fit_stab)
# 
#         
#     mtext("synth scan summary",side = 3,
#           line = - 2,
#           outer = TRUE)
#   }
# }
plot_summary_all <- function(d){
  {
    par(mfrow = c(5,6))
    # vapply(colnames(d), FUN = function(par){plot(d$scan, d$par)})
    for(i in seq(1, length(colnames(d))-20, 1)){
      
      plot(d$scan, d[,colnames(d)[i]], xlab = "scan", ylab = (colnames(d)[i]))
    }
    
    mtext("synth scan summary",side = 3,
          line = - 2,
          outer = TRUE)
  }
}



d <- read.table("scan_dat_pruned.txt", header = TRUE)

png("scan_summary.png", width = 16, height = 14, units = "in", res = 200)
plot_summary_all(d)
dev.off()



# plot more, compare probed subst ddGs at position scanned (pos1) and position static (pos2)

png("A1S_true_stab.png", width = 16, height = 4, units = "in", res = 200)
par(mfrow = c(1,6))
# vapply(colnames(d), FUN = function(par){plot(d$scan, d$par)})
for(i in seq(19, 19+5, 1)){
  print(i)
  plot(d$scan, d[,colnames(d)[i]], xlab = "scan", ylab = (colnames(d)[i]), col = "red")
  points(d$scan, d[,colnames(d)[i + 6]], xlab = "scan", ylab = (colnames(d)[i + 6]), col = "blue")
}
mtext("plot title",side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

# 
# png("A1S_true_stab2.png", width = 16, height = 4, units = "in", res = 200)
par(mfrow = c(1,2))

plot(d$scan, d$p1_.1_ddG, xlab = "scan", ylab = "p1_-1_ddG", col = "red", )
points(d$scan, d$p1_0_ddG, xlab = "scan", ylab = "p1_0_ddG", col = "blue")
points(d$scan, d$p1_1_ddG, xlab = "scan", ylab = "p1_1_ddG", col = "orange")

points(d$scan, d$p2_.1_ddG, xlab = "scan", ylab = "p2_-1_ddG", col = "red", pch = 4)
points(d$scan, d$p2_0_ddG, xlab = "scan", ylab = "p2_0_ddG", col = "blue", pch = 4)
points(d$scan, d$p2_1_ddG, xlab = "scan", ylab = "p2_1_ddG", col = "orange", pch = 4)

# 
# 
# mtext("plot title",side = 3,
#       line = - 2,
#       outer = TRUE)
# dev.off()



