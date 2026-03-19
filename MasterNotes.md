#This is the master notes
#Sample Data: Palsa Frozen Replicate A: SAMN08784143
#UNDER mjd356
$mkdir virome_project
$cd virome_project
$mkdir adapters
$mkdir reads
$mkdir reads/raw
$mkdir reads/trimmed
$mkdir fastqc_out
$mkdir fastqc_out/raw fastqc_out/trimmed
$mkdir log+scripts
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
$nano Trim_MS_1.txt
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

$module load mamba/
#environment manager
$mamba create -y -n megahit-env -c conda-forge -c bioconda megahit
#create an environment that contains megahit, and name the environment "megahit-env."
$mkdir assembled
$nano assembled/megahit_MS.txt

#SLURM SCRIPT

#!/bin/bash
#SBATCH --job-name=megahit_MS2   	# how job appears in the queue
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8  
#SBATCH --mail-type=END, FAIL
#SBATCH --mail-user=mjd356@georgetown.edu
#SBATCH --mem=32G                         
#SBATCH --time=03:00:00                   
#SBATCH --output=/home/mjd356/virome_project/logs+scripts/megahit_test1.%j.out      
#SBATCH --error=/home/mjd356/virome_project/logs+scripts/megahit_test1.%j.err       

#note %j = job ID

# ==== Load mamba/conda module (students: no need to change) ====
module load mamba
source $(mamba info --base)/etc/profile.d/conda.sh

# Activate the environment where you had MEGAHIT installed
conda activate megahit-env

# ==== Set paths and filenames (students: edit this block!) ====

# Directory where the cleaned reads live
READDIR=/home/mjd356/virome_project/reads/trimmed

# Input read files (paired-end)
READ1=${READDIR}/SRR6996005_forward_paired.fastq.gz
READ2=${READDIR}/SRR6996005_reverse_paired.fastq.gz

# Output directory (give it a name, it will be created by MEGAHIT)
OUTDIR=/home/mjd356/virome_project/assembled/megahit_out

# ==== Run MEGAHIT ====

megahit \
  -1 ${READ1} \
  -2 ${READ2} \
  -t ${SLURM_CPUS_PER_TASK} \
  -o ${OUTDIR}

echo "Done. Contigs should be in: ${OUTDIR}/final.contigs.fa"

$-q -u mjd356

$module load mamba/
$mamba activate megahit-env
$mamba install -c bioconda seqkit
$grep ">" final.contigs.fa | wc
#59,754 contigs
$less final.contigs.fa
#Looks like DNA
$seqkit stats -a final.contigs.fa
#longest contig: 185,956
#minimum: 200
#sum_len: 30,261,886
#avg_len: 506.4
#Q1: 340, Q2: 397, Q3:  495 sum_gap: 0 N50: 463, N50_snum: 1559  GC%: 55.1 all otehr values are zero


$ module load mamba
$ mamba create -y -n vs2-env -c conda-forge -c bioconda virsorter
$ mamba activate vs2-env
$ rm -rf db 					# just in case there is a failed attempt before, 
$ virsorter setup -d db -j 4		# run setup for the database

$cd logs+scripts
$nano virsorter_MS.txt

SLURM SCRIPT

#!/bin/bash
#SBATCH --job-name=virsorter_MS  
#SBATCH --nodes=1
#SBATCH --mail-type=END,FAIL 
#SBATCH --mail-user=mjd356@georgetown.edu
#SBATCH --cpus-per-task=8                 
#SBATCH --mem=20G                         
#SBATCH --time=03:00:00                  
#SBATCH --output=/home/mjd356/virome_project/logs+scripts/virsorter.out      
#SBATCH --error=/home/mjd356/virome_project/logs+scripts/virsorter.err    

# ==== Load mamba (students: no need to change) ====
module load mamba
source $(mamba info --base)/etc/profile.d/conda.sh

# Activate the environment where you had VirSorter2 installed
mamba activate vs2-env

# ==== Set paths and filenames (students: edit this block!) ====
#set up directories
INDIR=/home/mjd356/virome_project/assembled/	 #directory where input will come from
OUTROOT=/home/mjd356/virome_project/virsorter/	 #directory output will go
mkdir -p "${OUTROOT}"					 #new directory to be created for output files

SAMPLE_ID= SAMN08784143                #just the basic sample name (sample2 ?)
DBDIR=/home/mjd356/virome_project/db           #just the basic sample name (sample2 ?)
INPUT="${INDIR}/final.contigs.fa"			 #contig file name/location
OUTDIR="${OUTROOT}/vs2-${SAMPLE_ID}" 			 #where you’ll find the outputs
mkdir -p "${OUTDIR}"


# ==== Run virsorter2 with >5kb cutoff and DNA virus categories first
echo "Running VirSorter2 on ${INPUT}"
virsorter run \
  -w "${OUTDIR}" \
  -i "${INPUT}" \
  --db-dir “${DBDIR}”\
  --keep-original-seq \
  --include-groups dsDNAphage,NCLDV,ssDNA \
  --min-length 5000

echo "Done."

$cd ..
$sbatch virsorter_MS.txt








