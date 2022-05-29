library(ggplot2)
library(tidyverse)
library(ggrepel)
library(reshape2)
library(cowplot)


load("gmma_result.rda")
subst_all <- subst

subst_sig <- subst_all[subst_all$active > 0, ]

# scatter posiion dependtent ddGs, f_active
# exclude non_sig subst: S3K, A1C, A1K, V2P, V2W
ggplot(data = subst_sig, aes(x = log10(active+1e-10) , y = ddG_glob, label = rownames(subst_sig), color = resi))+
  scale_colour_gradientn(colours=rainbow(4))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)+
  geom_text_repel(size = 2, max.overlaps = 40)+
  labs(y = "ddGs_gmma", x = "log10(active count)")+
  theme_cowplot()
ggsave("gmma_laub_all_scatter_factive_ddGs.png", device = "png", height = 4, width = 6)



# # include non_sig subst: S3K, A1C, A1K, V2P, V2W
# ggplot(data = subst_all, aes(x = log10(active+1e-10) , y = ddG_glob, label = rownames(subst_sig), color = resi))+
#   scale_colour_gradientn(colours=rainbow(4))+
#   geom_point()+
#   geom_smooth(method = "lm", se = FALSE)+
#   geom_text_repel(size = 2, max.overlaps = 40)+
#   labs(y = "ddGs_gmma", x = "log10(active count)")+
#   theme_bw()




# position wise ddGs distribution boxplot
#########################################
png("gmma_laub_all_ddGs_boxplot_pos.png")
boxplot(subst$ddG_glob ~ subst$resi, xlab = "Residue position", ylab = "ddGs fitted gmma", main = "Position specific ddGs distribution" )
dev.off()



# suppressor plot top 5 gmma ddGs scores
########################################
# fraction binding across all variants
# "supressor plot"
binding <- mutant[mutant$signal == 1, ]

f_bind <- table(binding$N_sub) / table(mutant$N_sub)

f_bind <- data.frame(f_bind)
colnames(f_bind) <- c("n_subst", "f_bind")
f_bind$n_subst <- as.numeric(f_bind$n_subst) - 1
f_bind <-cbind(f_bind, subst = "WT")


# only variants with given gmma subst 
subst_f_bind <- function(s){
  sub <- mutant[str_detect(mutant$sub, s), ]
  sub_bind <- sub[sub$signal == 1, ] 
  
  f_sub <- table(factor(sub_bind$N_sub, levels = 1:4)) / table(sub$N_sub) # remember to count 0 if not binding, but if fraction very low consider counting raw numbers
  
  f_sub <- data.frame(f_sub)
  colnames(f_sub) <- c("n_subst", "f_bind")
  f_sub$n_subst <- as.numeric(f_sub$n_subst) - 1
  f_sub <- cbind(f_sub, subst = s)
  
  return(f_sub)
}


# do plot for chosen subset of subst

f_bind_plot <- function(subst_subset){
  # prepare for plotting
  
  var_fs <- lapply(subst_subset, subst_f_bind)
  var_fs <- do.call(rbind, var_fs)
  fs_plot <- rbind(f_bind, var_fs)
  
  # plot
  ggplot(fs_plot, aes(x = n_subst, y = f_bind, col = subst))+
    geom_point()+
    geom_line()+
    # scale_y_continuous(trans = "log10")+ # scale y axis log10 for easier separation of lines
    labs(x ="Number of additional substitutions", y = "Fraction of variants binding")+
    scale_color_manual("Data subset", values = c("green", "orange", "blue", "red", "magenta", "black"))+
    scale_linetype_manual(values = c("solid","solid","solid","solid","solid","dashed"))+
    
    theme_cowplot()
  ggsave("suppressorplot_gmma_laub_top5_2.png", height = 5, width = 8)
  
  # smooth curve as gmma figure 1f
  
  
}


# # give subst to be included in plot 
# ###################################
# 
# # all subst 
# all_subst <- rownames(subst)
# 
# # neutral (most stabilising subst)
# top_subst <- rownames(subst[subst$eff == "neu", ])
# 
# # top5 subst
top_ddG_subst <- subst[order(subst$ddG_glob, decreasing = FALSE), ]
top5_subst <- rownames(top_ddG_subst[1:5, ])

# # top10
# top10_subst <- rownames(top_ddG_subst[1:10, ])
# 
# 
# ###################################

f_bind_plot(top5_subst)

