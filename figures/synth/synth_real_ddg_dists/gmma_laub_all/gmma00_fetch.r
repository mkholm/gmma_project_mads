
# All combinations of 4 amino acids
all4 = read.table("all4.seq")
podgornaia15 = read.csv2("podgornaia15_sup_table1.csv")
colnames(podgornaia15) = c("binding")

# All varints found to bind in Podgornaia15
d = data.frame(var = all4$V2, bind = 0)
rownames(d) = d$var
d[podgornaia15$binding,"bind"] = 1

# Variants are given as a string of 4 amino acids that are varied in the study
wt = c("A","V","S","T")
var2subst = function(v) {
    v_list = strsplit(v,"")[[1]]
    all_subst = paste0(wt, c(1,2,3,4), v_list)
    paste0( all_subst[which(v_list != wt)], collapse=":")
}

# Translate the string variants to a colon-separated string of substitutions
d$subst = sapply(d$var, var2subst)

# Mark variants without substitutions as wild-types
d[which(d$subst == ""),"subst"] = "WT"

# Column with number of substitutions
d$nsubst = sapply(d$subst, function(s) { length(strsplit(s,":")[[1]]) } )

# Sort variants based on number of substitutions so single mutants come first
d = d[order(d$nsubst),]

# Other data have a two additional columns, number of observations (e.g. dna variants) and uncertainty from averaging these
d$dummy1 = 1
d$dummy2 = 1.0

# Write a data file
write.table(d[,c("subst","dummy1","bind","dummy2")], file="podgornaia15.txt", row.names=F, col.names=F, quote=F)
