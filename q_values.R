
# q values (08102015)
data = read.table("/Users/chris/Desktop/p_vals.txt")
no_nas = data[is.finite(data[,3]) == TRUE,]
hist(as.matrix(no_nas[,3]), breaks = 1000)
library(fdrtool)
pvals = as.vector(no_nas[,3])
fdr = fdrtool( pvals, statistic = 'pvalue') #, color.figure=FALSE)
qvals = as.matrix(fdr$qval) # estimated Fdr values 
hist(qvals)
final_df = cbind(no_nas, qvals)
names(final_df) = cbind('gene_name', 'num_isoforms', 'pval', 'iso1_version', 'qval')
filtered = (final_df[final_df$qval < 0.05, ])
write.table(filtered, 'Desktop/sigQvals_pt05_1fpkm_ILR_1014015.txt', row.names = F, quote = F)