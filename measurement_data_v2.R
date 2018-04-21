




rqtl = read.table('/Users/chris/Dropbox/Kane_lab/Alt_splicing/rqtl_table_04062016_justSplicingPhenos.txt', header = T, sep = '\t')[3:102,]  # two blank header lines
rils = read.table('/Users/chris/Documents/Kane_lab/Alt_splicing_project/Morph_analysis/rils_just_ids.txt', sep = '\t')
rnaPCA = prcomp(rqtl)
rnaPCA = prcomp(rqtl, center = T, scale = T)
summary(rnaPCA) # 6.4%, 4.9%, 4.9%, 4.5%, 3.2%... 7PCs > 3%, 14PC >2%, 32PCs > 1%
rnaData = cbind(rils, rnaPCA$x)
names(rnaData) = c("RIL", seq(1,100))




	######## germination data
germData = read.table('/Users/chris/Dropbox/Kane_lab/Alt_splicing/Morphology/germinationData_CS.txt', sep = '\t', header = T)[1:128,] # cutting last four rows, which aren't rils
df <- merge(rnaData, germData, by = "RIL", sort = F)
	# proportion germinated
X = (df[,103]/df[,102])
PCnum = 1; Y = as.matrix(df[, 1+PCnum]); summary(lm( Y ~ X )) # ****
plot(X ~ Y)






