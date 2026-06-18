##### Please modify to assign right path to the files accordingly #####

### 1. RNA-seq data preprocessing commands (HSC/CD4/CD8) ###
# Adapter trimming and quality filtering (fastp)

PREFIX=Sample1
THREADS=8

### Input files
fastq_R1=/path/to/raw/${PREFIX}_R1.fastq.gz
fastq_R2=/path/to/raw/${PREFIX}_R2.fastq.gz

### Output files
trimmed_R1=/path/to/output/trimmed/${PREFIX}_R1.trimmed.fastq.gz
trimmed_R2=/path/to/output/trimmed/${PREFIX}_R2.trimmed.fastq.gz

qc_html=/path/to/output/qc/${PREFIX}_fastp.html
qc_json=/path/to/output/qc/${PREFIX}_fastp.json

fastp
--in1 $fastq_R1
--in2 $fastq_R2
--out1 $trimmed_R1
--out2 $trimmed_R2
--detect_adapter_for_pe
--cut_front
--cut_tail
--length_required 40
--correction
--thread $THREADS
--html $qc_html
--json $qc_json

# Alignment to reference genome

genome_dir=/path/to/STAR_index/
gtf=/path/to/Mus_musculus.GRCm38.102.gtf
star_prefix=/path/to/output/bam/${PREFIX}

STAR
--runThreadN $THREADS
--genomeDir $genome_dir
--readFilesIn $trimmed_R1 $trimmed_R2
--readFilesCommand zcat
--sjdbGTFfile $gtf
--outFileNamePrefix $star_prefix
--outSAMtype BAM SortedByCoordinate
--quantMode GeneCounts

# Post-alignment processing (sorting + indexing)

$bam_sorted=${star_prefix}Aligned.sortedByCoord.out.bam
$bam_final=/path/to/output/final/${PREFIX}.sorted.bam

sambamba sort
-t $THREADS
-o $bam_final
$bam_sorted

sambamba index
-t $THREADS
$bam_final


### 2. ATACseq data preprocessing commands ###
## Prepare data
PREFIX=Sample1

fastq_R1=${PREFIX}_R1.fastq.gz
fastq_R2=${PREFIX}_R2.fastq.gz

# Cutadapt
cutadapt -j 8 \
  -a CTGTCTCTTATACACATCT \
  -A CTGTCTCTTATACACATCT \
  -O 5 -q 20,20 -m 20 --pair-filter=any \
  -o ${PREFIX}_R1.trimmed.fastq.gz -p ${PREFIX}_R2.trimmed.fastq.gz \
  $fastq_R1 $fastq_R2

# Bowtie2 and Samtools sort
bowtie2 -X 2000 \
  -x /path/to/bowtie2_index/mm10 \
  -1 ${PREFIX}_R1.trimmed.fastq.gz \
  -2 ${PREFIX}_R2.trimmed.fastq.gz \
  --threads $THREADS 2> ${PREFIX}.align.log \
| samtools sort -@ <M> -O BAM -l 0 -m 1G \
  -o ${PREFIX}.align.sorted.bam

# Samtools filtering, sort, mark duplicates, index
samtools view -b -@ $THREADS \
  ${PREFIX}.align.sorted.bam \
  chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
  chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 \
  chrX chrY \
  > ${PREFIX}.keepChr.bam

sambamba sort -t $THREADS -o ${PREFIX}.sorted.bam ${PREFIX}.keepChr.bam

sambamba markdup -t <THREADS_MARKDUP> \
  ${PREFIX}.sorted.bam \
  ${PREFIX}.markdup.bam

sambamba index -t <THREADS_INDEX> ${PREFIX}.markdup.bam

# macs2 callpeak
macs2 callpeak -t ${PREFIX}.markdup.bam -f BAMPE -g mm -n ${PREFIX} --nomodel --keep-dup all \
    --outdir outdir 2> ${PREFIX}.macs2.log


# HOMER motif analysis
CLOSED_DAR_BED= #convert closed DARs to BED format
OPEN_DAR_BED= #convert open DARs to BED format
ALL_PEAKS_BED= #convert all peaks to BED format
OUTDIR= #output directory

findMotifsGenome.pl ${CLOSED_DAR_BED} mm10 ${OUTDIR}/HOMER_closed_DAR \
  -size 100 \
  -bg ${ALL_PEAKS_BED} \
  -allowTargetOverlaps \
  -keepOverlappingBg \
  -len 8,10,12 \
  -p ${THREADS} \
  > ${OUTDIR}/HOMER_closed_DAR.log 2>&1

findMotifsGenome.pl ${OPEN_DAR_BED} mm10 ${OUTDIR}/HOMER_open_DAR \
  -size 100 \
  -bg ${ALL_PEAKS_BED} \
  -allowTargetOverlaps \
  -keepOverlappingBg \
  -len 8,10,12 \
  -p ${THREADS} \
  > ${OUTDIR}/HOMER_open_DAR.log 2>&1


### 3. EMseq data preprocessing commands ###
PREFIX=Sample1
fastq_R1=/path/to/raw/${PREFIX}_R1.fastq.gz
fastq_R2=/path/to/raw/${PREFIX}_R2.fastq.gz

trimmed_dir=/path/to/output/trimmed

## Trim Galore paired-end adapter trimming and quality filtering

~/tools/TrimGalore-0.6.10/trim_galore \
  --paired \
  --clip_R1 10 --clip_R2 10 \
  --three_prime_clip_R1 10 --three_prime_clip_R2 10 \
  -q 20 \
  -O 1 \
  --illumina \
  -j $THREADS \
  -o $trimmed_dir \
  --fastqc \
  $fastq_R1 $fastq_R2

## Bismark alignment to bisulfite genome
$genome_folder=/path/to/BismarkGenome
$aligned_dir=/path/to/output/aligned
$tmp_dir=/path/to/tmp

./tools/Bismark-0.24.2/bismark
  --temp_dir $tmp_dir \
  --bowtie2 \
  --bam \
  --dovetail \
  -p $THREADS \
  $genome_folder \
  -1 $trimmed_dir/$(basename $fastq_R1 .fastq.gz)_val_1.fq.gz \
  -2 $trimmed_dir/$(basename $fastq_R2 .fastq.gz)_val_2.fq.gz \
  -o $aligned_dir

## Bismark methylation extraction
$methylation_dir=/path/to/output/methylation

./tools/Bismark-0.24.2/bismark_methylation_extractor \
  -p --parallel $THREADS \
  --cytosine_report \
  --bedGraph \
  -o $methylation_dir \
  --gzip \
  --genome_folder $genome_folder \
  $aligned_dir/${PREFIX}_val_1_bismark_bt2_pe.bam
