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

$module load checkv						#its available as a module on the HPC
$checkv download_database ./				#make sure you’re in your checkv folder!

$nano checkvslurm

SLURM script:

#!/bin/bash
#SBATCH --job-name=checkv_MS
#SBATCH --output=/home/mjd356/virome_project/logs+scripts/checkv-%j.out
#SBATCH --error=/home/mjd356/virome_project/logs+scripts/checkv-%j.err
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mjd356@georgetown.edu


# ==== Load checkv program module (students: no need to change) ====

module load checkv


# ==== Set variables, paths, and filenames (students: edit this block!) ====

CHECKVDB="/home/mjd356/virome_project/checkv/checkv-db-v1.5"

SAMPLE_ID="vOTUs"
INPUT="/home/mjd356/virome_project/votus/votus_final.fna"
OUTDIR="/home/mjd356/virome_project/checkv/${SAMPLE_ID}"

mkdir -p "${OUTDIR}"


# ==== run checkv (students: no need to change. The second line is the command) ====
echo "Running CheckV on ${INPUT}"
checkv end_to_end "${INPUT}" "${OUTDIR}" -d "${CHECKVDB}" -t ${SLURM_CPUS_PER_TASK}
echo "Done."

$ls/checkv


$less quality_summary_votus.tsv.  

Complete 1 k141_59834
Not determined 2
Low Quality 18 
Medium Quality 1

$ gcloud storage cp gs://gu-biology-dept-class/ClassProject/votus_10kb_6samples.fna [location]
#run vOTU against class contigs 

$mkdir bowtie2
$cp votus_10kb_6samples.fna bowtie2
$ srun --pty bash  #connect to compute node
$ module load bowtie2
$ bowtie2-build votus_10kb_6samples.fna votu_index
$ exit

nano bowtieslurm

SLURM SCRIPT:


#!/bin/bash
#SBATCH --jobname=bowtie2_vOTUs                 
#SBATCH --output=/home/mjd356/virome_project/logs+scripts/bowtie-%j.out
#SBATCH --error=/home/mjd356/virome_project/logs+scripts/bowtie-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mjd356@georgetown.edu
#SBATCH --time=8:00:00                                
#SBATCH --mem=16G                        

# ---------SET UP----------
SAMPLE= ”sample2_mjd356”  	#whatever your sample # is!
INDEX="/home/mjd356/virome_project/bowtie2/votu_index"
OUTPUTDIR="/home/mjd356/virome_project/bowtie2/${SAMPLE}"

# --------- LOAD MODULES ----------
module purge
module load bowtie2/2.5.4

# --------- RUN BOWTIE2 AND PIPE TO SAMTOOLS ----------

# First make output and log directories; move into OUTPUTDIR
mkdir -p "${OUTPUTDIR}"
cd "${OUTPUTDIR}"
mkdir -p logs

echo "Running bowtie2 on sample ${SAMPLE}"

bowtie2 -p 8 -x "${INDEX}" -1 "/home/mjd356/virome_project/reads/trimmed/SRR6996005_forward_paired.fastq.gz" -2 "/home/mjd356/virome_project/reads/trimmed/SRR6996005_reverse_paired.fastq.gz" \
| samtools view -bS - > "${SAMPLE}.bam"

echo "Finished running bowtie2 and performing compression"

#---------sort and index files
echo "Sorting"
samtools sort "${SAMPLE}.bam" > "${SAMPLE}_sorted.bam"

echo "Indexing"
samtools index "${SAMPLE}_sorted.bam"

echo "Finished ${SAMPLE}"


$ gcloud storage cp [file] gs://gu-biology-dept-class/ClassProject/bam


# -------- DATA VISUALIZATION ---------

Downloaded https://console.cloud.google.com/storage/browser/_details/gu-biology-dept-class/ClassProject/votus_6samples_coverm.tsv?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&project=gcp-gu-hpc-medusa from class bucket

install.packages("tidyverse")
install.packages("readr")
read.csv("ClassProject_votus_6samples_coverm.tsv")

# ==== EDIT THESE THREE LINES ====
filename <- "ClassProject_votus_6samples_coverm.tsv"  # Excel file name
tpm_threshold <- 10                                  # keep vOTUs with max TPM > this
heatmap_colors <- c("#440154", "#31688e", "#35b779", "#fde725")

# =================================
library(readr)
library(readxl)
library(pheatmap)

# 1. Read the Excel file
cov <- read_tsv(filename)

# 2. Keep Contig and TPM columns only
tpm_cols <- grepl("TPM$", names(cov))
cov_tpm <- cov[ , c("Contig", names(cov)[tpm_cols])]

# Remove S1k141_26921_full
cov_tpm <- subset(cov_tpm, Contig != "S1k141_26921||full")

# 3. Optional filter: drop very low-abundance vOTUs
cov_tpm$max_tpm <- apply(cov_tpm[ , -1], 1, max, na.rm = TRUE)
cov_tpm <- subset(cov_tpm, max_tpm > tpm_threshold)
cov_tpm$max_tpm <- NULL

# If everything got filtered out (threshold too high), warn and stop
if (nrow(cov_tpm) == 0) {
  stop("No vOTUs passed the TPM threshold. Try lowering tpm_threshold.")
}

# 4. Make matrix for heatmap (rows = vOTUs, cols = samples)
mat <- as.matrix(cov_tpm[ , -1])
rownames(mat) <- cov_tpm$Contig

# 5. Log-transform for nicer color scaling
mat_log <- log10(mat + 1)

# 6. Draw heat map
pheatmap(mat_log,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "none",
         color = colorRampPalette(heatmap_colors)(100),
         fontsize_row = 8,
         fontsize_col = 10,
         main = "vOTU relative abundance (log10 TPM + 1)")
