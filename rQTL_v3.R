
options(digits=22)
library(qtl)

# current (03292016) table, ILR transformed, 200+ genes or something
my_data = read.cross(format = "csv",  dir = "Dropbox/Kane_lab/Alt_splicing", file = "rqtl_table_04062016.txt", genotypes = c("1238", "hetero", "HA89") , sep= "\t") 
my_data = read.cross(format = "csv",  dir = "Dropbox/Kane_lab/Alt_splicing", file = "rqtl_NonSignificant_12292017.txt", genotypes = c("1238", "hetero", "HA89") , sep= "\t") 



plot.map(my_data)
summary(my_data)

#calculating genotype probabilities with step size 1cM and assuming given rate of genotyping errors=0.001
my_data_jitter = jittermap(my_data)
summary(my_data_jitter)
genotype_probs = calc.genoprob(my_data_jitter, step=0, error.prob=0.001)
genes = names(genotype_probs$pheno)

big_ch_list = c()
for (i in 1:length(genes))
{
  	  print( c("Phenotype", i, "out of", length(genes), ":", genes[i]) )
	  ch_list <- c()	
	  pos_list = c()
	  other_counter = 0
	  
	  # getting the significance threshold using permutations
	  int_map = scanone(genotype_probs, pheno.col=i, n.perm = 1000)
   	  # Threshold from the permutation test
	  sum = summary(int_map, alpha=c(0.05, 0.10))
	  	  	  
	  #interval mapping function:
	  int_map = scanone(genotype_probs, pheno.col=i)
	  for (j in 1:length(int_map$lod))
	  {
	    #if (int_map$lod[j] == Inf)
	    #{
	    # 	int_map$lod[j] = 100
	    #}
	   
	    #making list of putative map locations
	    if (int_map$lod[j] > sum[1])
	    {	
	    	other_counter = other_counter +1
	    	#print( c(int_map$chr[j], int_map$pos[j]) )
			ch = int_map$chr[j]
			pos = int_map$pos[j]
	    	ch_list = c(ch_list, ch)
	    	pos_list = c(pos_list, pos)
	    	score = int_map$lod[j]
	    	big_ch_list = rbind(big_ch_list, cbind(genes[i], ch, pos, score) )
	    }
	  }

	  # plotting	  
	  #par(mfrow=c(2,1), mar=c(5,5,1,5))    
	  #plot(int_map, ylab="LOD score", xlab=genes[i])
	  #abline(sum[1],0, col = 'red')

	  # analyzing variance explained by putative snps
#	  if ( length(pos_list) > 0)
#	  {
#	  	print("imputing for chromosomes:")
#	  	print(ch_list)
#	  	impute <- sim.geno(genotype_probs, step=2, n.draws=128, err=0.001) #imputations
#	  	qtl <- makeqtl(impute, chr=ch_list, pos=pos_list)
#	  	plot(qtl)
#	  }

	  #	user input between plots
	 #cat ("Press [enter] to continue")
     #line <- readline() #makes you press enter before continuing

}	
write.table(big_ch_list, 'Desktop/qtl_NonSignificant_1000perm_12292017.txt', sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)










####################################################
# qtl effect size analysis and cis vs trans 7.17.16
####################################################
options(digits=22) # indicating decimal precision
library(qtl) 
#genomeLocations = read.table('sigQvals_pt05_1fpkm_ILR_04062016_cMpositions.txt', header = F)
genomeLocations = read.table('/Users/chris/Dropbox/Kane_lab/Alt_splicing/sigQvals_pt05_1fpkm_ILR_04062016_cMpositions.txt', header = F)
genomeLocations = read.table('/Users/chris/Dropbox/Kane_lab/Alt_splicing/nonSig_firstHit_01022018_cMpositions.txt', header = F)
geneN = paste(t(data.frame(strsplit(as.matrix(genomeLocations[,8]), "_")))[,1], t(data.frame(strsplit(as.matrix(genomeLocations[,8]), "_")))[,2], sep = '_')
geneC = genomeLocations[,1]
geneP = genomeLocations[,3]
locationMat = cbind(geneN, geneC, geneP)
# current (03292016) table, ILR transformed, 200+ genes or something
#my_data = read.cross(format = "csv",  dir = ".", file = "rqtl_table_04062016.txt", genotypes = c("1238", "hetero", "HA89") , sep= "\t") # reads in data
my_data = read.cross(format = "csv",  dir = "/Users/chris/Dropbox/Kane_lab/Alt_splicing/", file = "rqtl_table_04062016.txt", genotypes = c("1238", "hetero", "HA89") , sep= "\t") # reads in data
par(bg = 'grey'); plot.map(my_data) # map of marker positions
summary(my_data) # doesn't tell us much
#calculating genotype probabilities with step size 1cM and assuming given rate of genotyping errors=0.001
	# no longer jittering because I want to keep the positions
#my_data_jitter = jittermap(my_data) # jitters the positions- I think because the decimals aren't precise enough to distinguish marker positions
genotype_probs = calc.genoprob(my_data, step=0, error.prob=0.001) # don't really know, but calculates "true" genotype probabilities given the multipoint marker data
genes = names(genotype_probs$pheno)

big_list = c()
imputed <- sim.geno(genotype_probs, step=2, n.draws=128, err=0.001) # imputed genpotypes (I think) anyways important
all_Qs = c()
for (i in 1:length(genes))
{
  	  print( c("Phenotype", i, "out of", length(genes), ":", genes[i]) )
	  ch_list <- c()	
	  pos_list = c()
	  other_counter = 0
	  
	  # getting the significance treshold using permutations
	  int_map = scanone(genotype_probs, pheno.col=i, n.perm = 10) # single-qtl scan with permutations for getting signifcance threshold
	  out = scanone(genotype_probs, pheno.col=i) # single-qtl scan of actual data
   	  # Threshold from the permutation test
	  threshold = summary(int_map, alpha=c(0.05, 0.10))[1]
	  # Results above the 0.05 threshold
	  lodPeaks = summary(out, perms=int_map, alpha=0.05) # this good command gives you the MAX of each significant LOD peak. Just what I need.

	  # for-loop for checking for multiple peaks per chromosome (god I hate R... already had to do this in python)
	  # we are identifying a "peak" as the maximum-lod markers separated by at least 10cM with lod below the threshold
	  # I think I want to redo with the PEAKS being at least 10cM apart (and some dead space)
	  Qs = c()
	  for (chrm in 1:17)
	  {
		  maxLoc = c(0,0,0)
		  deadZone = F
		  thisChrom = out[out$chr == chrm,]
		  for (pos in 1:dim(thisChrom)[1])
		  {
	  		LODscore = thisChrom[pos,]$lod
			if (LODscore > threshold) # if we're inside an LOD peak
			{
				deadZone = F
				if (LODscore > maxLoc[3])
				{
					maxLoc = thisChrom[pos,]
				}
			}
			else # if we're outside an LOD peak
			{
				if (deadZone == F) # if this is the first location outside of peak, begin dead zone
				{
					deadZone = T
					deadStart = thisChrom[pos,]$pos
				}
				else # if we're already inside a dead zone, keep track of how long it is
				{
					currentPos = thisChrom[pos,]$pos
					if ( (currentPos - deadStart) > 10) # (alt version) hits within 10cM of DEAD space are the same qtl

					{
						if (maxLoc[3] > 0) # this if is keeping us from adding non-lod hits 
						{
							Qs = rbind(Qs, maxLoc)
						}
						maxLoc = c(0,0,0)
					}
				}
			}
	  	  }
	  }
	  	
	  if ( is.null(dim(Qs)) )
	  {
	  	print('no peaks')
	  }
	  else
	  {
		  # find qtl regions: within 15% of max	
		  	# so, need a table of start and end locations with LOD scores within 15% of each peak
	  		# I'll use while loops, inching away from each peak until less than 15%. Yeah
		  qtlLocationTable = c()
		  for (peak in 1:dim(Qs)[1])
		  {
		  	relavantHits = out[out$chr == Qs[peak, ]$chr, ]
			sameIndex = which(relavantHits$pos == Qs[peak, ]$pos)
			startIndex = sameIndex  
			endIndex = sameIndex
			# first while loop going backward
			backwardsExit = F
			while (backwardsExit == F)
			{
				newIndex = startIndex - 1
				newHit = relavantHits[newIndex, ]
				if (newIndex < 1) # real indices only
				{
					backwardsExit = T 
				}
				else if (newHit$lod >= (.85*Qs[peak, ]$lod) ) # if within 15% of max, save it
				{
					startIndex = newIndex
				}
				else # else, you're done. Found the range (in the backwards direction)
				{
					backwardsExit = T
				}
			}
			# second while loop going forward
			forwardsExit = F
			while (forwardsExit == F)
			{
				newIndex = endIndex + 1
				newHit = relavantHits[newIndex, ]
				if (newIndex > dim(relavantHits)[1] ) # real indices only
				{
					backwardsExit = T 
				}
				else if (newHit$lod >= (.85*Qs[peak, ]$lod) )
				{
					endIndex = newIndex
				}
				else
				{
					forwardsExit = T
				}				
			}
			startPos = relavantHits[startIndex, ]$pos
			endPos = relavantHits[endIndex, ]$pos
			qtlLocationTable = rbind(qtlLocationTable, c(startPos, endPos))
			colnames(qtlLocationTable) = c('start', 'end')
		  }
		  Qs = cbind(as.matrix(Qs), qtlLocationTable )
		  
		  # interval mapping function for estimating variance explained
		  # qtl <- makeqtl(imputed, chr=lodPeaks[,1], pos=lodPeaks[,2] ) # using one peak per chromosome
		  qtl <- makeqtl(imputed, chr=as.numeric(Qs[,1]), pos=as.numeric(Qs[,2]) ) # multiple peaks per chromosome, separated by at least 10cM
		  out.fq <- fitqtl(imputed, qtl=qtl, pheno.col = i) # lod in output is different than scanone() output because its a multiple qtl model
		  
		  #mqmplot.circle(genotype_probs,  out[c(2000:2010), ])
	  		# (1) takes the imputed genotype probabilities (?)
		  	# (2) the qtl of interest (?)
		  	# (3) the phenotype column we're interested
		  if (dim(Qs)[1] == 1)
		  {
		 	varExplained = out.fq$result.full[1,5]			  	
		  }
		  else 
		  {
		 	varExplained = out.fq$result.drop[,4]	
		  }
		  
		  # lastly, for each peak, is it cis or trans (qtl region within 1 cM?)
		  # need 1. the gene location. 
		  # need 2. the qtl locations. Have these in "Qs" variable
		  cisTrans = c()
		  for (peak in 1:dim(Qs)[1])
		  {
		  	#print(Qs[peak,])
		  	gN = genes[i]
		  	gL = locationMat[locationMat[,1] == gN, 3]	# start
		  	gC = locationMat[locationMat[,1] == gN, 2]	# start
			if (length(gL) == 1) # would be false if this gene has no genomic position (we couldn't identify it)
			{
				trans = T
				if ( gC == Qs[peak,1] )
				{
					if ((abs( as.numeric(gL) - as.numeric(Qs[peak,4]) ) < 1) | (abs( as.numeric(gL) - as.numeric(Qs[peak,5]) ) < 1) )
					{
						if (gC == Qs[peak,1])
							{
								trans = F
								print('cis')					
							}
					}
				}
				if (trans == F)
				{
					cisTrans = c(cisTrans, 'cis')
				}
				else
				{
					cisTrans = c(cisTrans, 'trans')
				}
			}
			else
			{
				cisTrans = rep('NA', dim(Qs)[1])					
			}
		  }
		  
		  geneName = rep(genes[i], dim(Qs)[1])			
		  datum = cbind( geneName, data.frame(Qs), varExplained, cisTrans)
		  big_list = rbind(big_list, datum)
	  }
}	
write.table(big_list, 'noSig_varExplained_1000perm_01022018.txt', sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)







############### interaction tests #############
addint(imputed, qtl=qtl, pheno.col = i, formula=y~Q1+Q2+Q3+Q4) # lod in output is different than scanone() output because its a multiple qtl model







########### plotting big genetic map with gene and qtl locations ##############
options(digits=22)
library(qtl)
myRegs = read.table('Dropbox/Kane_lab/Alt_splicing/varExplained_1000perm_09132016.txt')
myRegs = read.table('Desktop/tempPleio18.txt')
myQs = na.omit(myRegs)[, 1:6]
myQtls <- makeqtl(imputed, chr=as.numeric(myQs[,2]), pos=as.numeric(myQs[,3]) ) # multiple peaks per chromosome, separated by at least 10cM

genomeLocations = read.table('/Users/chris/Dropbox/Kane_lab/Alt_splicing/sigQvals_pt05_1fpkm_ILR_04062016_cMpositions.txt', header = F)
genomeLocations = read.table('/Users/chris/Desktop/tempGenes.txt', header = F)
geneN = paste(t(data.frame(strsplit(as.matrix(genomeLocations[,4]), "_")))[,1], t(data.frame(strsplit(as.matrix(genomeLocations[,4]), "_")))[,2], sep = '_')
geneC = genomeLocations[,1]
geneP = genomeLocations[,3]
locationMat = cbind(geneN, geneC, geneP)
myLs = c()
for (gene in 1:dim(locationMat)[1])
{
			print(gene)
		  	gN = locationMat[gene, 1]
		  	qtls = myQs[myQs[,1] == gN, 3]	# start

			if (length(qtls) != 0) # would be false if this gene has no genomic position (we couldn't identify it)
			{
		  	gL = locationMat[locationMat[,1] == gN, 2:3]	# start

				myLs = rbind(myLs, gL)
			}
}
myGeneLocs <- makeqtl(imputed, chr=as.numeric(myLs[,1]), pos=as.numeric(myLs[,2]) ) # multiple peaks per chromosome, separated by at least 10cM

#		  plot(myGeneLocs, justdots = T, col = 'red')
#		      par(new=T)
#		  plot(myQtls, justdots = T, col = 'blue')

# circos !    circular plots

# need to do this by hand, unfortunately. But not that bad!
genoTable = read.table('Dropbox/Kane_lab/Alt_splicing/precise_genoColsFiltered_11092015.txt')
par(mar = c(5,5,1,1))
plot(1, type="n", xlab="chromosome", ylab="Location (cM)", axes = F, xlim = c(0.5,17.5), ylim = c(100,0), cex.lab = 1.5)
axis(side = 1, at = seq(1,17), cex.axis = 1.5)
axis(side = 2, at = c(100, 80, 60,40,20,0), cex.axis = 1.5)
for (chrom in 1:17)
{
	print(chrom)
	cur = genoTable[, genoTable[2,] == chrom]
	curGenes = myLs[myLs[,1] == chrom,]
	curQtl = myQs[myQs[,2] == chrom,]
	minCm = min(as.numeric(as.matrix(cur[3,])))
	maxCm = max(as.numeric(as.matrix(cur[3,])))
	segments(chrom, 0, chrom, maxCm-minCm, lwd = 1.5) # max minus min to start each at 0
	# add every genetic marker
	for (marker in 1:dim(cur)[2])
	{
		segments(chrom-0.15, as.numeric(as.matrix(cur[3,marker]))-minCm, chrom+0.15, as.numeric(as.matrix(cur[3,marker]))-minCm, lwd = 1) 
	}	
	# add every gene location
	if (dim(curGenes)[1] > 0){
	for (L in 1:dim(curGenes)[1])
	{
		points(chrom-0.15, as.numeric(as.matrix(curGenes[L,2]))-minCm, col = 'red', pch = 16, cex = 2)
	}}
	# add every qtl
	if (dim(curQtl)[1] > 0){
	for (q in 1:dim(curQtl)[1])
	{
		print(curQtl[q,3])
		points(chrom+0.15, as.numeric(as.matrix(curQtl[q,3]))-minCm, col = 'blue', pch = 16, cex = 2)
		for (L2 in 1:dim(myLs)[1])
		{
			toChrom = as.numeric(myLs[L2,1])
			arrows(chrom+0.15, as.numeric(as.matrix(curQtl[q,3]))-minCm, toChrom+0.25, as.numeric(as.matrix(myLs[L2,2]))-minCm, lwd = 1, col = 'green') 
		}	
	}}
}
legend(12.5,80, c("gene", "qtl"), col = c('red', 'blue'), lty = c(1,1,1), pch = c(16, 16), cex = 2)











