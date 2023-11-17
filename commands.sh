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






mamba activate clara
time STAR --runThreadN 150 --runMode genomeGenerate --genomeDir . --genomeFastaFiles Homo_sapiens_assembly19_1000genomes_decoy.fasta --sjdbGTFfile star.gencode.v19.transcripts.patched_contigs.gtf --sjdbOverhang 100

#Running from ~clara_parabricks
for i in SRR10307483_1.fastq.gz; do
	j=${i/_1.fastq.gz/_2.fastq.gz}
	k=${i/_1.fastq.gz/}	
time docker run --rm --gpus all --volume $(pwd):/workdir --volume $(pwd):/outputdir \
    nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1 \
    pbrun rna_fq2bam \
    --in-fq /workdir/$i /workdir/$j \
    --output-dir /outputdir/ \
	--genome-lib-dir /workdir/testIn/GATK_resources/ \
    --ref /workdir/testIn/GATK_resources/Homo_sapiens_assembly19_1000genomes_decoy.fasta \
    --out-bam /outputdir/$k.bam \
    --read-files-command zcat \
	--num-threads 100 \
	--two-pass-mode Basic
done

# The version of clara parabricks with splitncigar is no longer available, using gatk for next step
mamba activate gatk
time gatk --java-options "-Xmx100G" SplitNCigarReads -I SRR10307483.bam -O SRR10307483.sncr.bam -R testIn/GATK_resources/Homo_sapiens_assembly19_1000genomes_decoy.fasta


time docker run --rm --gpus all --volume $(pwd):/workdir --volume $(pwd):/outputdir \
    nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1 \
	pbrun bqsr \
	--ref /workdir/testIn/GATK_resources/Homo_sapiens_assembly19_1000genomes_decoy.fasta \
	--in-bam /workdir/SRR10307483.sncr.bam \
	--knownSites /workdir/testIn/GATK_resources/Homo_sapiens_assembly19_1000genomes_decoy.known_indels.vcf \
	--knownSites /workdir/testIn/GATK_resources/Mills_and_1000G_gold_standard.indels.b37.sites.vcf \
	--out-recal-file /outputdir/SRR10307483.recal

# Don't need this step as can be run during haplotype caller
#time docker run --rm --gpus all --volume $(pwd):/workdir --volume $(pwd):/outputdir \
#    nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1 \
#	pbrun applybqsr \
#	--ref /workdir/testIn/GATK_resources/Homo_sapiens_assembly19_1000genomes_decoy.fasta \
#	--in-bam /workdir/SRR10307483.sncr.bam \
#	--in-recal-file /workdir/SRR10307483.recal \
#	--out-bam /outputdir/SRR10307483.recal.bam \
#	--num-threads 100

time docker run --rm --gpus all --volume $(pwd):/workdir --volume $(pwd):/outputdir \
    nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1 \
	pbrun haplotypecaller \
	--gvcf \
	--ref /workdir/testIn/GATK_resources/Homo_sapiens_assembly19_1000genomes_decoy.fasta \
	--in-bam /workdir/SRR10307483.recal.bam \
	--in-recal-file /workdir/SRR10307483.recal \
	--out-variants /outputdir/SRR10307483.g.vcf \
	--num-htvc-threads 100 \
	--rna





for i in SRR10307489_1.fastq.gz; do
	j=${i/_1.fastq.gz/_2.fastq.gz}
	k=${i/_1.fastq.gz/}	
time docker run --rm --gpus all --volume $(pwd):/workdir --volume $(pwd):/outputdir \
    nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1 \
    pbrun rna_fq2bam \
    --in-fq /workdir/$i /workdir/$j \
    --output-dir /outputdir/ \
	--genome-lib-dir /workdir/testIn/GATK_resources/ \
    --ref /workdir/testIn/GATK_resources/Homo_sapiens_assembly19_1000genomes_decoy.fasta \
    --out-bam /outputdir/$k.bam \
    --read-files-command zcat \
	--num-threads 100 \
	--two-pass-mode Basic
done

# The version of clara parabricks with splitncigar is no longer available, using gatk for next step
mamba activate gatk
time gatk --java-options "-Xmx100G" SplitNCigarReads -I SRR10307489.bam -O SRR10307489.sncr.bam -R testIn/GATK_resources/Homo_sapiens_assembly19_1000genomes_decoy.fasta


time docker run --rm --gpus all --volume $(pwd):/workdir --volume $(pwd):/outputdir \
    nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1 \
	pbrun bqsr \
	--ref /workdir/testIn/GATK_resources/Homo_sapiens_assembly19_1000genomes_decoy.fasta \
	--in-bam /workdir/SRR10307489.sncr.bam \
	--knownSites /workdir/testIn/GATK_resources/Homo_sapiens_assembly19_1000genomes_decoy.known_indels.vcf \
	--knownSites /workdir/testIn/GATK_resources/Mills_and_1000G_gold_standard.indels.b37.sites.vcf \
	--out-recal-file /outputdir/SRR10307489.recal

time docker run --rm --gpus all --volume $(pwd):/workdir --volume $(pwd):/outputdir \
    nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1 \
	pbrun haplotypecaller \
	--gvcf \
	--ref /workdir/testIn/GATK_resources/Homo_sapiens_assembly19_1000genomes_decoy.fasta \
	--in-bam /workdir/SRR10307489.sncr.bam \
	--in-recal-file /workdir/SRR10307489.recal \
	--out-variants /outputdir/SRR10307489.g.vcf \
	--num-htvc-threads 100 \
	--rna



















mamba activate ncbi
while read line; do echo $line; fasterq-dump -3 -m 1000 -e 10 $line; pigz -p 10 *fastq; echo $line >> done.ae; done <../split_PRJNA577693ae



# Notes from the gatk best references github page which list the location of the various reference, known sites and gtf files. May be out of date

"##_COMMENT2": "REFERENCE FILES",
 "RNAseq.refFasta": https://storage.googleapis.com/gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.fasta
 "RNAseq.refFastaIndex": https://storage.googleapis.com/gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.fasta.fai
 "RNAseq.refDict": https://storage.googleapis.com/gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.dict

 "##_COMMENT4": "RESOURCE FILES",
 "RNAseq.dbSnpVcf": https://storage.googleapis.com/gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.dbsnp138.vcf
 "RNAseq.dbSnpVcfIndex" https://storage.googleapis.com/gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.dbsnp138.vcf.idx
 
 "RNAseq.knownVcfs": https://storage.googleapis.com/gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Mills_and_1000G_gold_standard.indels.b37.sites.vcf https://storage.googleapis.com/gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.known_indels.vcf

 "RNAseq.knownVcfsIndices":    https://storage.googleapis.com/gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Mills_and_1000G_gold_standard.indels.b37.sites.vcf.idx https://storage.googleapis.com/gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.known_indels.vcf.idx

 "RNAseq.annotationsGTF": https://storage.googleapis.com/gatk-test-data/intervals/star.gencode.v19.transcripts.patched_contigs.gtf