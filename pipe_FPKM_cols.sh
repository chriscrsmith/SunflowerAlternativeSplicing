

SAMPLE_LIST=$1

#python Scripts/FPKM_columns.py Clean_gene_lists/Significant/sigQvals_pt05_1fpkm_ILR_04062016.txt Trinity_output/TPMorFPKMs/RIL_fpkms/INX299_D08EFACXX_2_RIL008_CAAAAG_fpkm.txt Filtered_TPM_tables/filtered_1TPM_03092016.txt Trinity_output/Trinity_singleline.fasta TPM_columns/bamSizes.txt 1 # significant gene list
python Scripts/FPKM_columns.py Clean_gene_lists/BlastVerification/goodSplicing_cleanBlast_04062016_minusSignificant.txt Trinity_output/TPMorFPKMs/RIL_fpkms/INX299_D08EFACXX_2_RIL008_CAAAAG_fpkm.txt Filtered_TPM_tables/filtered_1TPM_03092016.txt Trinity_output/Trinity_singleline.fasta TPM_columns/bamSizes.txt 1 # non-significant gene list

for sample in $(cat $SAMPLE_LIST)
do
    #echo $sample
    #python Scripts/FPKM_columns.py Clean_gene_lists/Significant/sigQvals_pt05_1fpkm_ILR_04062016.txt Trinity_output/TPMorFPKMs/RIL_fpkms/$sample"_"fpkm.txt Filtered_TPM_tables/filtered_1TPM_03092016.txt Trinity_output/Trinity_singleline.fasta TPM_columns/bamSizes.txt 2 # significant genes
    python Scripts/FPKM_columns.py Clean_gene_lists/BlastVerification/goodSplicing_cleanBlast_04062016_minusSignificant.txt Trinity_output/TPMorFPKMs/RIL_fpkms/$sample"_"fpkm.txt Filtered_TPM_tables/filtered_1TPM_03092016.txt Trinity_output/Trinity_singleline.fasta TPM_columns/bamSizes.txt 2 # non significant genes
    
done
