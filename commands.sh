#Basic shell commands that were run during the initial manual analysis. These have not been automated as of yet
 
# Add creating required environments here

# ncbi environment for downloading the reads using fasterq-dump requires sra-tools (ideally >= v3)

# List of software in ebola environment - this can be trimmed to a smaller list as needed
STAR
samtools bcftools
bamtools
gatk4
r-base
r-curl
tophat
hisat2


# Activate ncbi environment and download test files from two samples

#Step 1 - download files

fasterq-dump -3 -e 100 -p SRR4888788; pigz *fastq
fasterq-dump -3 -e 100 -p SRR4888657; pigz *fastq


# Step 2 - create indexes
#Create the star mapping index 
STAR --runThreadN 100 --runMode genomeGenerate --genomeDir genomeTest/ --genomeFastaFiles genome/GCA_000001405.15_GRCh38_full_analysis_set.fna --sjdbGTFfile genome/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf --sjdbOverhang 100

# Step 3 - 2 step mapping of the samples using STAR

#This one failed
time STAR --readFilesIn SRR4888788_1.fastq.gz SRR4888788_2.fastq.gz --runMode alignReads --runThreadN 100 --genomeDir genomeTest/ --genomeFastaFiles genome/GCA_000001405.15_GRCh38_full_analysis_set.fna --readFilesCommand zcat --outFileNamePrefix pass1 --outSAMtype BAM Unsorted --outBAMcompression -1

#Trying
time STAR --readFilesIn SRR4888788_1.fastq.gz SRR4888788_2.fastq.gz --runMode alignReads --runThreadN 100 --genomeDir genomeTest/ --readFilesCommand zcat --outFileNamePrefix pass1 --outSAMtype BAM Unsorted --outBAMcompression -1

# This took ~13 minutes to run on the first pass. Don't think I actually need the resulting bam file so it can be deleted to save space

# Add filtering of splice junctions here - recommended otherwise the second pass will take a long time due to the number of splice sites added
#check the manual for what each column in out.tab corresponds to

cat pass1SJ.out.tab | awk '{if ($5 > 0 && $7 > 2 && ($8 < $7)) print}' | grep -v "^chrM" > pass1SJ.out.filtered.tab

# Once all samples have had a first pass reindex the genome with the new, filtered splice sites
# Create new genome index with splice junctions
mkdir splice
mv *tab splice/.
time STAR --runThreadN 100 --runMode genomeGenerate --genomeDir genomeTestSplice/ --genomeFastaFiles genome/GCA_000001405.15_GRCh38_full_analysis_set.fna --sjdbGTFfile genome/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf --sjdbOverhang 100 --sjdbFileChrStartEnd splice/*tab

# Reindexing took 18 minutes

# Remove the initial bam files
rm pass1*bam

# With a single, unfiltered sj tab file this took 18 minutes to run

# Run second pass with the splice junctions
time STAR --readFilesIn SRR4888788_1.fastq.gz SRR4888788_2.fastq.gz --runMode alignReads --runThreadN 100 --genomeDir genomeTestSplice/ --readFilesCommand zcat --outFileNamePrefix pass2 --outSAMtype BAM Unsorted --outBAMcompression -1

# Second pass took 12 minutes so not significantly different from first

# At 24 minutes per sample with 2 passes it will take 40 hours to run 100 samples, or 1.6 days. This does not include any time for downloading the raw fastq reads. At 30 minutes per sample for 2 passess it will need 2 days to create 100 bam files and that's before I event get to SNP calling.

# Once I have all of the bam files need to snp call with gatk
samtools sort pass2Aligned.out.bam -o pass2.sorted.bam -@ 100

gatk MarkDuplicates -I pass2.sorted.bam -M pass2.mkdup.log -O pass2.mkdup.bam # This took 49 minutesd!

time gatk --java-options "-Xmx500G" MarkDuplicates -I pass2.sorted.bam -M pass2.mkdup.log -O pass2.mkdup.500G.bam

#Switching to the MarkDuplicatesSpark algorithm cuts this down to 3 minutes (as it failed but at least gave me a usable error! Missing rea group tags)
time gatk --java-options "-Xmx500G" MarkDuplicatesSpark -I pass2Aligned.out.bam -M pass2.mkdup.log -O pass2.mkdupSPARK.500G.bam


time gatk --java-options "-Xmx500G" AddOrReplaceReadGroups --I pass2Aligned.out.bam -O pass2.RG.bam -PL Illumina -PU ABC -SM p1 -LB a
# This took 18 minutes to run

time gatk --java-options "-Xmx500G" MarkDuplicatesSpark -I pass2.RG.bam -M pass2.mkdup.log -O pass2.mkdupSPARK.500G.bam


library(dyplr)
byID=str_split_i(rownames(sampleTable), "P", 1)
byID=str_split_i(byID, "-", -1)
sampleTable$condition=byID

plot(pcaData[,"PC1"],pcaData[,"PC2"], pch=19, cex=0.1)
k=1
j=1
for (i in names(table(byID))){
	points(pcaData[pcaData$byID == i,"PC1"],pcaData[pcaData$byID == i,"PC2"], pch=k, col=j, cex=2)
	k=k+1
	j=j+1
	if (k > 25){k=1}
}



docker run --rm --gpus all --volume /testIn:/workdir --volume /testOut:/outputdir \
    --workdir /work \
    nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1 \
    pbrun rna_fq2bam \
    --in-fq /workdir/*R1* /workdir/*R2* \
    --genome-lib-dir /workdir/genome/ \
    --output-dir /outputdir/ \
    --ref /workdir/GCA_000001405.15_GRCh38_full_analysis_set.fna \
    --out-bam /outputdir/test1.bam \
    --read-files-command zcat



for i in *_R1.fastq.gz; do
	j=${i/_R1/_R2}
	k=${i/_R1.fastq.gz/}	
time docker run --rm --gpus all --volume $(pwd):/workdir --volume $(pwd):/outputdir \
    nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1 \
    pbrun rna_fq2bam \
    --in-fq /workdir/$i /workdir/$j \
    --output-dir /outputdir/ \
	--genome-lib-dir /workdir/genome2/star/ \
    --ref /workdir/genome2/GCA_000001405.15_GRCh38_full_analysis_set.fna \
    --out-bam /outputdir/$k.bam \
    --read-files-command zcat \
	--num-threads 56 \
	--two-pass-mode Basic
done



gshawli:ifp

tmux
while read line; do echo $line; sudo pigz -p 25 $line; done <RBtest1




tmux
mamba activate ncbi

while read line; do echo $line; fasterq-dump -3 -m 1000 -e 10 $line; pigz -p *fastq; echo $line >> done.ae; done <../split_PRJNA577693ae