!/bin/bash -i

#Start of script to automate SNP calling using clara parabricks gpu accelerated version of gatk - must be run on a server with active gpu for which docker and clara parabricks have already been set up. For current purposes use borgcube.

# Set the script running from the directory with the fastq files. For now the location of all of the other files are assumed to be in a subdirectory in that location called testIn. Will update later.

# This also assumes that the reference for star alignment has already been generated.

rm all.samples
for s in *_1.fastq.gz; do
echo $s >> all.samples
done

while read i; do
	echo "Starting sample" $i "- rna_fq2bam"
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

# The version of clara parabricks with splitncigar is no longer available, using gatk for next step
mamba activate gatk
echo "SplitNCigarReads"
time gatk --java-options "-Xmx100G" SplitNCigarReads -I $k.bam -O $k.sncr.bam -R testIn/GATK_resources/Homo_sapiens_assembly19_1000genomes_decoy.fasta

echo "BQSR"
time docker run --rm --gpus all --volume $(pwd):/workdir --volume $(pwd):/outputdir \
    nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1 \
	pbrun bqsr \
	--ref /workdir/testIn/GATK_resources/Homo_sapiens_assembly19_1000genomes_decoy.fasta \
	--in-bam /workdir/$k.sncr.bam \
	--knownSites /workdir/testIn/GATK_resources/Homo_sapiens_assembly19_1000genomes_decoy.known_indels.vcf \
	--knownSites /workdir/testIn/GATK_resources/Mills_and_1000G_gold_standard.indels.b37.sites.vcf \
	--out-recal-file /outputdir/$k.recal

echo "Haplotypecaller"
time docker run --rm --gpus all --volume $(pwd):/workdir --volume $(pwd):/outputdir \
    nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1 \
	pbrun haplotypecaller \
	--gvcf \
	--ref /workdir/testIn/GATK_resources/Homo_sapiens_assembly19_1000genomes_decoy.fasta \
	--in-bam /workdir/$k.sncr.bam \
	--in-recal-file /workdir/$k.recal \
	--out-variants /outputdir/$k.g.vcf \
	--num-htvc-threads 100 \
	--rna

rm -f *bam*
rm -f *bai*
rm -f *recal

done <all.samples