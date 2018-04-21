

 #just once
#done   bwa index ref.fa
#done   samtools faidx ref.fa  


SAMPLE_LIST=$1  # line separated list of sample names (must correspond to fasta names somehow)
DATA_PATH=$2  # path to the directory containing the sample directories
REF=$3

echo "using $SAMPLE_LIST..."
echo "using $DATA_PATH..."
echo "using $REF"
echo

for sample in $(cat $SAMPLE_LIST)
do
     #making sample directories
    echo $sample
    echo "mkdir $sample"_"files"
    mkdir $sample"_"files

     #change into sample dir
    echo "cd $sample"_"files"
    cd $sample"_"files

     #link .fq's into the sample dir's
    echo "ln -s $DATA_PATH/Parent_$sample/$sample"_"forward_paired.fq"
    ln -s $DATA_PATH/Parent_$sample/$sample"_"forward_paired.fq
    echo "ln -s $DATA_PATH/Parent_$sample/$sample"_"reverse_paired.fq"
    ln -s $DATA_PATH/Parent_$sample/$sample"_"reverse_paired.fq 

     # trimmomatic
#    java -jar /data4/AltSplicing09182015/Scripts/Software/Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 $sample"_"forward_paired.fq $sample"_"reverse_paired.fq Parent_HA89A/HA89A_forward_paired.fq Parent_HA89A/HA89A_forward_unpaired.fq Parent_HA89A/HA89A_reverse_paired.fq Parent_HA89A/HA89A_reverse_unpaired.fq ILLUMINACLIP:contaminants.fa:2:30:10 LEADING:10 TRAILING:10 MINLEN:90

     #alignment and SNP calling
    echo "bwa mem $REF $sample"_"forward_paired.fq $sample"_"reverse_paired.fq 1> $sample.sam 2> $sample.bwa.stderr" 
    bwa mem $REF $sample"_"forward_paired.fq $sample"_"reverse_paired.fq 1> $sample.sam 2> $sample.bwa.stderr

    echo "samtools view -b -o $sample.bam -S $sample.sam 2> $sample.sam_view.stderr"
    samtools view -b -o $sample.bam -S $sample.sam 2> $sample.sam_view.stderr

    echo "samtools sort $sample.bam $sample.bam.sorted 2> $sample.sam_sort.stderr"
    samtools sort $sample.bam $sample.bam.sorted 2> $sample.sam_sort.stderr

    echo "samtools index $sample.sorted.bam 2> $sample.sam_index.stderr"
    samtools index $sample.bam.sorted.bam 2> $sample.sam_index.stderr

    ### new bcftools version (bcftools call replaces bcftools view) ###
#    echo "samtools mpileup -P ILLUMINA -u -g -t DP -f ../bwa/masked_high_coveragesoap_assembly1_k47_min1000bp_scafs.fasta 189../bwa/bam_bai/aln*sorted.bam | bcftools call -c -v -o out1.vcf"

     #switch out of sample dir
    echo "cd ../"
    cd ..
    echo

done


 #once overall? OR just the parents
#samtools mpileup -gf HA412_noAltSplice_400bpmin_chriss_genes_replaced.fa */*sorted.bam > parents.bcf
#bcftools view -vcg parents.bcf  > parents.vcf 
