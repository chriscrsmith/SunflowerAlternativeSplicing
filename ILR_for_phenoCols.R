
library(robCompositions)
options(digits=22)
my_data = read.table("/Users/chris/Dropbox/Alt_splicing_clean_files/RQTL_mapping09102015/FPKM_cols09102015.txt", header = T)
dim(my_data) # looks good, 124 rows of RILs, and 284 cols of genes with isoforms 1 and 2




# now, for every pair of isoforms, want to do an ILR transform of the data
new_matrix = c()
new_gene = T
for (col in 1:dim(my_data)[2])
{
	
	if (new_gene == T)
	{
		new_gene = F
		iso1s = my_data[col]
	}
	
	else
	{
		new_gene = T
		iso2s = my_data[col]
		x = isomLR( cbind(iso1s, iso2s) )	
			
		iso_name_split = strsplit( names(my_data)[col], '_')
		first = iso_name_split[[1]][1]
		second = iso_name_split[[1]][2]
		gene_name = as.matrix(paste(first, second, sep = '_'))
		new_col = rbind(gene_name, x)
		new_matrix = cbind(new_matrix, new_col)
	}
}
dim(new_matrix)


#write.table(new_matrix, 'Desktop/temp.txt', sep = '\t', row.names = F, col.names = F, quote = F)




