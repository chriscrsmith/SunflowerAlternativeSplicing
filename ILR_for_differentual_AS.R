

########## ANALYSIS #############
library(robCompositions)
options(digits=22)
my_file = read.table('Desktop/retention_deletion_1tpm_2iso_proportions_10292015_withVersions.txt')
pop = c('w', 'd', 'w', 'w', 'd', 'd')
transformed_df = c() # simply for observing the distribution 


tempwildcount = 0
tempdomcount = 0
df = c()
xs = c() # this vector is simply for collecting transformed values for visualization
i1s = c() # vector of isoform1 proportions to view distribution 
new_gene = TRUE # true means currently open to new gene
gene_num = 0
for (line in 1:dim(my_file)[1])
{
	print(c('line', line, 'out of', dim(my_file)[1]))
	newline = my_file[line,]
	# this bit is saving proportion isoform1 to view distribution
	newline[1]
	
	# if gene line, then simply record the gene name and number of isoforms
	if (new_gene == TRUE)
	{		
		gene_name = newline[2]
		gene_num = gene_num +1
		print(c('gene_num', gene_num))
		num_isos = newline[3]
		iso_versions = c(newline[11], newline[12])
		individual_counter = 0 
		new_gene = FALSE
		fpkm_matrix = c()
	}
	
	# if data line, then add data to growing matrix
	else
	{
		individual_counter = individual_counter + 1		
		new_row = c()
		for (iso in 1:as.integer(num_isos) )
		{
			fpkm = newline[iso]
			new_row = c(new_row, as.numeric(as.matrix(fpkm)) )
		}
		fpkm_matrix = rbind(fpkm_matrix, new_row)
		#this bit is just for visualizing the distribution
		i1s = c(i1s, as.numeric(as.matrix(newline[1]))/as.numeric(as.matrix(newline[2])) )
	}
	
	# if end of data lines, then do ILR transform and stats test
	# then switch 'new_gene' back to true to get ready for the next gene
	colinear = FALSE
	if (individual_counter == 6) # six individuals
	{
		new_gene = TRUE
		x = isomLR(fpkm_matrix)
		xs = c(xs, x)
		transformed_df = c(transformed_df, x[1,1]) # simply for observing the distribution
		transformed_df = c(transformed_df, x[2,1]) 
		transformed_df = c(transformed_df, x[3,1]) 
		transformed_df = c(transformed_df, x[4,1]) 
		transformed_df = c(transformed_df, x[5,1]) 
		transformed_df = c(transformed_df, x[6,1]) 

		# stats
		if (dim(x)[2] == 1)
		{
			test = t.test(x ~ pop)
			p = as.numeric( as.matrix( (as.character(test)[3]) ) )
			print( c(mean(fpkm_matrix[c(1,3,4), 1]), mean(fpkm_matrix[c(2,5,6), 1] ) ) )
			if ( mean(fpkm_matrix[c(1,3,4), 1]) > mean(fpkm_matrix[c(2,5,6), 1]) ) # getting wild/dom version
			{
				wild = iso_versions[1]
				dom = iso_versions[2]	
				tempwildcount = tempwildcount +1
			}
			else
			{
				wild = 	iso_versions[2]
				dom = iso_versions[1]
				tempdomcount = tempdomcount +1
			}			
			names(wild) = 'wild'
			names(dom) = 'dom'
			datum = data.frame(as.character(as.matrix(gene_name)), as.integer(num_isos), p, data.frame(wild), data.frame(dom) )#, stringsAsFactors=FALSE)
			print(datum)
			df = rbind(df, datum)
		}
		else
		{
			# checking for linear independence of columns in x
			c = cor(x)
			for (i in 1:length(c) )
			{
				current = c[i]
				if (current != 1)
				{
					if (abs(current) > .5)
					{
						colinear = TRUE					
					}
				}
			}
			if (colinear == FALSE)
			{	
#				print(cor(x))
#				print(summary(manova(x ~ pop)) )
				sum = summary(manova(x ~ pop))
				split = strsplit(as.character(sum)[4], ',')
				p = as.numeric(as.matrix(split[[1]][11]))
				datum = ( c(as.character(as.matrix(gene_name)), as.integer(num_isos), p) )
				df = rbind(df, datum)
			}
		}
	}
}
#write.table(df, 'Desktop/p_vals_v2.txt', quote = F, row.names = F, col.names = F)








