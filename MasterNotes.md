#This is the master notes
#Sample Data: Palsa Frozen Replicate A: SAMN08784143
#UNDER mjd 356
$module load anaconda3
$conda create -n sra_env -c bioconda sra-tools
$conda activate sra_env
$prefetch SAMN08784143
$fasterq-dump *.sra
$gzip *.fastq
$srun --pty bash
$module load fastqc
$fastqc -j
$mkdir -p fastqc_out
$fastqc -o fastqc_out SRR6996005.sra_1.fastq.gz
$fastqc -o fastqc_out SRR6996005.sra_2.fastq.gz
$ls fastqc_out
#Trimmomatic Slurm:
  
#!/bin/bash
#SBATCH --job-name=”trim_MS”
#SBATCH --output="%x.o%j"
#SBATCH --mail-type=END,FAIL --mail-user=mjd356@georgetown.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=8G

shopt -s expand_aliases
module load trimmomatic

BASE=/home/mjd356/virome_project

R1=$BASE/fastq_files/raw/SRR6996005.sra_1.fastq.gz
R2=$BASE/fastq_files/raw/SRR6996005.sra_2.fastq.gz
ADAPTERS=$BASE/adapters/TruSeq3-PE.fa
OUTDIR=$BASE/fastq_files/trimmed

trimmomatic PE -threads $SLURM_CPUS_PER_TASK \
    $R1 $R2 \
    $OUTDIR/SRR6996005_forward_paired.fastq.gz $OUTDIR/SRR6996005_forward_unpaired.fastq.gz \
    $OUTDIR/SRR6996005_reverse_paired.fastq.gz $OUTDIR/SRR6996005_reverse_unpaired.fastq.gz \
ILLUMINACLIP:$ADAPTERS:2:30:10 \ SLIDINGWINDOW:4:20 MINLEN:50

$fastqc -o fastqc_out SRR6996005_forward_paired.fastq.gz SRR6996005_forward_unpaired.fastq.gz SRR6996005_reverse_paired.fastq.gz SRR6996005_reverse_unpaired.fastq.gz

