## clone repository
cd /workspaces/class-variantcalling
mkdir -p analysis
cd analysis

## sym link so we do not change the repository itself
mkdir -p raw_data
cd raw_data
ln -s /workspaces/class-variantcalling/datasets-class-variantcalling/reads/*.gz .
cd ..
mkdir -p alignment
cd alignment

## now we can perform the alignment with BWA

bwa mem \
-t 2 \
-R "@RG\tID:sim\tSM:normal\tPL:illumina\tLB:sim" \
/workspaces/class-variantcalling/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
/workspaces/class-variantcalling/analysis/raw_data/normal_1.000+disease_0.000_1.fq.gz \
/workspaces/class-variantcalling/analysis/raw_data/normal_1.000+disease_0.000_2.fq.gz \
| samtools view -@ 8 -bhS -o normal.bam -

## Real time: 176.099 sec; CPU: 256.669 sec

bwa mem \
-t 2 \
-R "@RG\tID:sim\tSM:disease\tPL:illumina\tLB:sim" \
/workspaces/class-variantcalling/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
/workspaces/class-variantcalling/analysis/raw_data/normal_0.000+disease_1.000_1.fq.gz \
/workspaces/class-variantcalling/analysis/raw_data/normal_0.000+disease_1.000_2.fq.gz \
| samtools view -@ 8 -bhS -o disease.bam -

## Real time: 173.232 sec; CPU: 256.204 sec


# sort the bam file
samtools sort -o normal_sorted.bam normal.bam
samtools sort -o disease_sorted.bam disease.bam

# index the bam file
samtools index normal_sorted.bam
samtools index disease_sorted.bam


# Marking duplicates

gatk MarkDuplicates \
-I normal_sorted.bam \
-M normal_metrics.txt \
-O normal_md.bam

gatk MarkDuplicates \
-I disease_sorted.bam \
-M disease_metrics.txt \
-O disease_md.bam


### recalibrating

gatk BaseRecalibrator \
   -I normal_md.bam \
   -R /workspaces/class-variantcalling/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   --known-sites /workspaces/class-variantcalling/datasets_reference_only/gatkbundle/dbsnp_144.hg38_chr21.vcf.gz \
   --known-sites /workspaces/class-variantcalling/datasets_reference_only/gatkbundle/Mills_and_1000G_gold_standard.indels.hg38_chr21.vcf.gz \
   -O normal_recal_data.table

gatk BaseRecalibrator \
   -I disease_md.bam \
   -R /workspaces/class-variantcalling/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   --known-sites /workspaces/class-variantcalling/datasets_reference_only/gatkbundle/dbsnp_144.hg38_chr21.vcf.gz \
   --known-sites /workspaces/class-variantcalling/datasets_reference_only/gatkbundle/Mills_and_1000G_gold_standard.indels.hg38_chr21.vcf.gz \
   -O disease_recal_data.table


#### Apply recalibration

gatk ApplyBQSR \
   -R /workspaces/class-variantcalling/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -I normal_md.bam \
   --bqsr-recal-file normal_recal_data.table \
   -O normal_recal.bam

gatk ApplyBQSR \
   -R /workspaces/class-variantcalling/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -I disease_md.bam \
   --bqsr-recal-file disease_recal_data.table \
   -O disease_recal.bam


tar -zvcf alignments.tar.gz *_recal.b*


### variant calling

cd /workspaces/class-variantcalling
mkdir -p analysis/variants
cd analysis/variants


## first single sample discovery

gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R /workspaces/class-variantcalling/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -I /workspaces/class-variantcalling/analysis/alignment/normal_recal.bam \
   -O normal.g.vcf.gz \
   -ERC GVCF

gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R /workspaces/class-variantcalling/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -I /workspaces/class-variantcalling/analysis/alignment/disease_recal.bam \
   -O disease.g.vcf.gz \
   -ERC GVCF

## then consolidate the 2 files

mkdir -p tmp

### on AMD64 this code ######
## combine the files into one
gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
      -V normal.g.vcf.gz \
      -V disease.g.vcf.gz \
      --genomicsdb-workspace-path compared_db \
      --tmp-dir /workspaces/class-variantcalling/analysis/variants/tmp \
      -L chr21

### on ARM64 (Mac M1 chip) this code
## combine the files into one
 gatk CombineGVCFs \
   -R /workspaces/class-variantcalling/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -V normal.g.vcf.gz \
   -V disease.g.vcf.gz \
   -O cohort.g.vcf.gz

### on AMD64 this code ######
### finally we can call the genotypes jointly
gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R /workspaces/class-variantcalling/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -V gendb://compared_db \
   --dbsnp /workspaces/class-variantcalling/datasets_reference_only/gatkbundle/dbsnp_146.hg38_chr21.vcf.gz \
   -O results.vcf.gz



### on ARM64 (Mac M1 chip) this code
### finally we can call the genotypes jointly
gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R /workspaces/class-variantcalling/datasets_reference_only/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -V cohort.g.vcf.gz \
   --dbsnp /workspaces/class-variantcalling/datasets_reference_only/gatkbundle/dbsnp_146.hg38_chr21.vcf.gz \
   -O results.vcf.gz

   
#### ANNOTATE THE SAMPLE

mkdir -p /workspaces/class-variantcalling/analysis/variants/cache
cd /workspaces/class-variantcalling/analysis/variants

## download hg38 (UCSC) version of database
snpEff download -v hg38 -dataDir /workspaces/class-variantcalling/analysis/variants/cache

### to execute snpeff we need to contain the memory
snpEff -Xmx4g ann -dataDir /workspaces/class-variantcalling/analysis/variants/cache -v hg38 results.vcf.gz >results_ann.vcf


### filter variants

cat results_ann.vcf | grep "#CHROM" | cut -f 10-

grep "#" results_ann.vcf >filtered_variants.vcf
cat results_ann.vcf | grep HIGH | perl -nae 'if($F[10]=~/0\/0/ && $F[9]=~/1\/1/){print $_;}' >>filtered_variants.vcf
cat results_ann.vcf | grep HIGH | perl -nae 'if($F[10]=~/0\/0/ && $F[9]=~/0\/1/){print $_;}' >>filtered_variants.vcf

sudo conda install bioconda::snpsift

SnpSift extractFields \
-s "," -e "." \
filtered_variants.vcf \
"CHROM" "POS" "ID" "GEN[disease].GT" "GEN[normal].GT" ANN[*].GENE ANN[*].EFFECT

SnpSift extractFields \
-s "," -e "." \
filtered_variants.vcf \
"CHROM" "POS" "ID" "REF" "ALT" "GEN[*].GT" ANN[0].GENE ANN[0].EFFECT
