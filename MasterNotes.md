# Metagenomic Viral Analysis Pipeline


## Dataset
- Source: NCBI SRA
- BioSample: SAMN08784143
- Run ID: SRR6996005

## Workflow Overview
1. Download raw sequencing data
2. Assess read quality with FastQC
3. Trim adapters and low-quality bases with Trimmomatic
4. Assemble reads into contigs with MEGAHIT
5. Identify viral contigs with VirSorter2
6. Filter and cluster viral contigs into vOTUs
7. Assess viral genome quality with CheckV
8. Map reads back to reference vOTUs with Bowtie2
9. Analyze abundance and diversity in R

## Step 1: Set up project directories and download sequencing reads

The first step is to create a reproducible project structure so that raw reads, trimmed reads, logs, and outputs are separated. This makes the workflow easier to follow and reduces path errors in downstream scripts. We then download the sequencing data from the SRA database and convert it into FASTQ format, which stores the raw reads and their quality scores.

-------------------
```bash
# Create a main project directory to store all files for this workflow
mkdir virome_project

# Move into the project directory so all subsequent files are organized here
cd virome_project

# Create a folder to store adapter sequences used during trimming
mkdir adapters

# Create a parent directory for sequencing reads
mkdir reads

# Create a subdirectory for raw (unprocessed) sequencing reads
mkdir reads/raw

# Create a subdirectory for trimmed (processed) reads after quality control
mkdir reads/trimmed

# Create a directory to store FastQC output reports
mkdir fastqc_out

# Create a subdirectory for FastQC results from raw reads
mkdir fastqc_out/raw

# Create a subdirectory for FastQC results from trimmed reads
mkdir fastqc_out/trimmed

# Create a directory to store log files and SLURM scripts for reproducibility
mkdir logs_scripts


# Load the Anaconda module, which provides access to conda environments and bioinformatics tools
module load anaconda3

# Create a new conda environment named "sra_env" and install SRA tools from the bioconda channel
conda create -n sra_env -c bioconda sra-tools

# Activate the newly created environment so SRA tools can be used
conda activate sra_env


# Download sequencing data from the NCBI SRA using the BioSample accession
prefetch SAMN08784143

# Convert the downloaded .sra file(s) into FASTQ format (raw sequencing reads with quality scores)
fasterq-dump *.sra

# Compress FASTQ files to save disk space (these files are typically very large)
gzip *.fastq


# Move the compressed FASTQ files into the "raw" reads directory for organization
mv *.fastq.gz reads/raw/
```
--------------------

## Step 2: Run FastQC on raw reads

FastQC is used to evaluate the quality of the raw sequencing reads before any downstream analysis. This is important because poor-quality bases, adapter contamination, and abnormal sequence duplication can negatively affect trimming, assembly, and viral detection. By inspecting the raw reads first, we can justify why trimming is needed.

--------------------
```bash
# Start an interactive session on a compute node (instead of running on the login node)
# This is important because computational tools like FastQC should run on compute nodes in an HPC environment
srun --pty bash

# Load the FastQC module, which contains the software for quality control analysis of sequencing reads
module load fastqc

# Run FastQC and specify the output directory using -o
# The backslash "\" allows the command to continue onto the next line for readability
fastqc -o fastqc_out/raw \
# Input file 1: forward reads (paired-end read 1)
reads/raw/SRR6996005.sra_1.fastq.gz \
# Input file 2: reverse reads (paired-end read 2)
reads/raw/SRR6996005.sra_2.fastq.gz
```
--------------------

## Step 3: Trim adapters and low-quality bases with Trimmomatic

Trimming removes adapter contamination and low-quality sequence from the reads. This improves downstream assembly because low-quality bases can create assembly errors and reduce confidence in viral contigs. We use paired-end trimming so that the relationship between the two reads from the same DNA fragment is preserved whenever possible.

--------------------
```bash
nano trim_reads.slurm
```
--------------------

--------------------
###Slurm Script: Parameter choices (ILLUMINACLIP, SLIDINGWINDOW, MINLEN) were selected to balance read quality and data retention. Adapter clipping removes non-biological sequences, SLIDINGWINDOW:4:20 trims low-quality regions based on a common Q20 threshold (≈1% error rate), and MINLEN:50 ensures reads are long enough to be useful for assembly while discarding overly short, unreliable fragments.

```bash
#!/bin/bash
#SBATCH --job-name=trim_MS
#SBATCH --output=/home/mjd356/virome_project/logs_scripts/trim_%j.out
#SBATCH --error=/home/mjd356/virome_project/logs_scripts/trim_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mjd356@georgetown.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=8G

module load trimmomatic

BASE=/home/mjd356/virome_project
R1=$BASE/reads/raw/SRR6996005.sra_1.fastq.gz
R2=$BASE/reads/raw/SRR6996005.sra_2.fastq.gz
ADAPTERS=$BASE/adapters/TruSeq3-PE.fa
OUTDIR=$BASE/reads/trimmed

trimmomatic PE -threads ${SLURM_CPUS_PER_TASK} \
  $R1 $R2 \
  $OUTDIR/SRR6996005_forward_paired.fastq.gz $OUTDIR/SRR6996005_forward_unpaired.fastq.gz \
  $OUTDIR/SRR6996005_reverse_paired.fastq.gz $OUTDIR/SRR6996005_reverse_unpaired.fastq.gz \
  ILLUMINACLIP:$ADAPTERS:2:30:10 \
  SLIDINGWINDOW:4:20 \
  MINLEN:50
```
--------------------


--------------------
# Submit the SLURM job script to the HPC scheduler to run the trimming step on a compute node
```bash
sbatch trim_reads.slurm
```
--------------------

## Step 4: Evaluate trimmed reads

After trimming, FastQC is run again to confirm that read quality improved and that adapter contamination was reduced. This provides evidence that the trimming step worked as intended and that the cleaned reads are ready for assembly.

--------------------
```bash
# Load the FastQC module so the quality control tool is available
module load fastqc

# Run FastQC on all trimmed read files and save reports to the trimmed output directory
fastqc -o fastqc_out/trimmed \
reads/trimmed/SRR6996005_forward_paired.fastq.gz \
reads/trimmed/SRR6996005_reverse_paired.fastq.gz \
reads/trimmed/SRR6996005_forward_unpaired.fastq.gz \
reads/trimmed/SRR6996005_reverse_unpaired.fastq.gz
```
--------------------

## Step 5: Assemble reads into contigs with MEGAHIT

Metagenomic reads are too short to identify viral genomes directly, so they must first be assembled into longer contiguous sequences called contigs. Assembly uses overlapping k-mers from the reads to reconstruct longer DNA fragments. This step is essential because tools like VirSorter2 work much better on assembled contigs than on raw reads.
--------------------
```bash
# Load the mamba module, which is a fast package manager for creating bioinformatics environments
module load mamba

# Create a new environment called "megahit-env" and install MEGAHIT from conda-forge and bioconda
mamba create -y -n megahit-env -c conda-forge -c bioconda megahit

# Open a new file to write the SLURM script that will run MEGAHIT on the cluster
nano megahit_run.slurm
```
--------------------
--------------------
SLURM SCRIPT
```bash
#!/bin/bash
#SBATCH --job-name=megahit_MS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mjd356@georgetown.edu
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/home/mjd356/virome_project/logs_scripts/megahit_%j.out
#SBATCH --error=/home/mjd356/virome_project/logs_scripts/megahit_%j.err

module load mamba
source $(mamba info --base)/etc/profile.d/conda.sh
conda activate megahit-env

READDIR=/home/mjd356/virome_project/reads/trimmed
READ1=${READDIR}/SRR6996005_forward_paired.fastq.gz
READ2=${READDIR}/SRR6996005_reverse_paired.fastq.gz
OUTDIR=/home/mjd356/virome_project/assembled/megahit_out

mkdir -p /home/mjd356/virome_project/assembled

megahit \
  -1 ${READ1} \
  -2 ${READ2} \
  -t ${SLURM_CPUS_PER_TASK} \
  -o ${OUTDIR}

echo "Done. Contigs are in ${OUTDIR}/final.contigs.fa"
```
--------------------
--------------------
```bash
sbatch megahit_run.slurm
```
--------------------
## Step 6: Inspect assembly statistics

**After assembly, it is useful to summarize the contigs to understand whether the assembly is fragmented or contains long sequences suitable for viral detection. Metrics such as total number of contigs, N50, and maximum contig length help evaluate assembly quality.
**--------------------
```bash
# Load mamba module to manage conda environments and install tools
module load mamba

# Initialize conda in the current shell so environments can be activated
source $(mamba info --base)/etc/profile.d/conda.sh

# Activate the environment where MEGAHIT is installed
conda activate megahit-env

# Install SeqKit (a tool for working with FASTA/FASTQ files) into the current environment
mamba install -c bioconda seqkit


# Move into the directory containing assembled contigs
cd /home/mjd356/virome_project/assembled/megahit_out

# Count the number of contigs by counting FASTA headers (lines starting with ">")
grep ">" final.contigs.fa | wc -l

# Generate detailed assembly statistics (e.g., total length, N50, min/max contig length, GC content)
seqkit stats -a final.contigs.fa
```
--------------------
## Step 7: Identify viral contigs with VirSorter2

VirSorter2 is used to distinguish viral contigs from non-viral metagenomic contigs. This is necessary because the assembly contains sequences from many organisms, not just viruses. VirSorter2 uses viral hallmark genes, gene-content patterns, and similarity to known viral sequences to predict which contigs are likely viral. We also apply a minimum length cutoff because longer contigs provide more genomic context and reduce false positives.
--------------------
```bash
# Load mamba to manage environments and install bioinformatics tools
module load mamba

# Create a new environment called "vs2-env" and install VirSorter2 from conda-forge and bioconda
mamba create -y -n vs2-env -c conda-forge -c bioconda virsorter

# Activate the VirSorter2 environment so the tool can be used
mamba activate vs2-env

# Remove any existing database directory to avoid conflicts or corrupted installs
rm -rf db

# Download and set up the VirSorter2 database using 4 CPU threads
virsorter setup -d db -j 4
```
--------------------

```bash
nano virsorter_run.slurm
```
--------------------

--------------------
SLURM SCRIPT
```bash
#!/bin/bash
#SBATCH --job-name=virsorter_MS
#SBATCH --nodes=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mjd356@georgetown.edu
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=03:00:00
#SBATCH --output=/home/mjd356/virome_project/logs_scripts/virsorter_%j.out
#SBATCH --error=/home/mjd356/virome_project/logs_scripts/virsorter_%j.err

module load mamba
source $(mamba info --base)/etc/profile.d/conda.sh
mamba activate vs2-env

INPUT=/home/mjd356/virome_project/assembled/megahit_out/final.contigs.fa
OUTROOT=/home/mjd356/virome_project/virsorter
SAMPLE_ID=SAMN08784143
DBDIR=/home/mjd356/virome_project/db
OUTDIR=${OUTROOT}/vs2-${SAMPLE_ID}

mkdir -p "${OUTDIR}"

virsorter run \
  -w "${OUTDIR}" \
  -i "${INPUT}" \
  --db-dir "${DBDIR}" \
  --keep-original-seq \
  --include-groups dsDNAphage,NCLDV,ssDNA \
  --min-length 5000

echo "Done."
```
--------------------
--------------------
```bash
sbatch virsorter_run.slurm
```
--------------------
## Step 8: Assess viral genome quality with CheckV

VirSorter2 predicts which contigs are viral, but it does not fully evaluate how complete those genomes are. CheckV is used after viral detection to estimate genome completeness and classify viral contigs into categories such as complete, medium quality, or low quality. This helps determine how reliable the final viral dataset is.
--------------------
--------------------
```bash
# Load the CheckV module, which is used to assess viral genome quality and completeness
module load checkv

# Create a directory for CheckV files (the -p ensures no error if it already exists)
mkdir -p /home/mjd356/virome_project/checkv

# Move into the CheckV directory so all outputs and databases are stored here
cd /home/mjd356/virome_project/checkv

# Download the CheckV reference database into the current directory
checkv download_database ./
```
--------------------
```bash
nano checkv_run.slurm
```

SLURM SCRIPT
--------------------
```bash
#!/bin/bash
#SBATCH --job-name=checkv_MS
#SBATCH --output=/home/mjd356/virome_project/logs_scripts/checkv_%j.out
#SBATCH --error=/home/mjd356/virome_project/logs_scripts/checkv_%j.err
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mjd356@georgetown.edu

module load checkv

CHECKVDB="/home/mjd356/virome_project/checkv/checkv-db-v1.5"
INPUT="/home/mjd356/virome_project/votus/votus_final.fna"
OUTDIR="/home/mjd356/virome_project/checkv/vOTUs"

mkdir -p "${OUTDIR}"

checkv end_to_end "${INPUT}" "${OUTDIR}" -d "${CHECKVDB}" -t ${SLURM_CPUS_PER_TASK}
```
--------------------
```bash
sbatch checkv_run.slurm
```
--------------------
## Step 9: Map reads back to reference vOTUs with Bowtie2

After identifying and filtering viral sequences, reads are mapped back to a reference set of vOTUs to estimate abundance. This step is important because viral detection alone only tells us which viral genomes are present, while mapping lets us quantify how abundant each viral population is in the sample.
--------------------
```bash
# Create a directory for Bowtie2 files (reference genomes, index files, outputs)
mkdir -p /home/mjd356/virome_project/bowtie2

# Move into the Bowtie2 directory to organize all mapping-related files here
cd /home/mjd356/virome_project/bowtie2


# Download the class pooled vOTU reference file from Google Cloud storage
gcloud storage cp gs://gu-biology-dept-class/ClassProject/votus_10kb_6samples.fna .


# Start an interactive session on a compute node (required for heavier computations)
srun --pty bash

# Load Bowtie2, the aligner used to map reads to reference sequences
module load bowtie2

# Build a Bowtie2 index from the vOTU reference FASTA file
# This creates indexed files that allow fast alignment of reads
bowtie2-build votus_10kb_6samples.fna votu_index

# Exit the interactive compute session and return to the login node
exit
```
--------------------
```bash
nano bowtie2_run.slurm
```
--------------------
```bash
#!/bin/bash
#SBATCH --job-name=bowtie2_vOTUs
#SBATCH --output=/home/mjd356/virome_project/logs_scripts/bowtie_%j.out
#SBATCH --error=/home/mjd356/virome_project/logs_scripts/bowtie_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mjd356@georgetown.edu
#SBATCH --time=08:00:00
#SBATCH --mem=16G

module purge
module load bowtie2/2.5.4
module load samtools

SAMPLE="sample2_mjd356"
INDEX="/home/mjd356/virome_project/bowtie2/votu_index"
OUTPUTDIR="/home/mjd356/virome_project/bowtie2/${SAMPLE}"

mkdir -p "${OUTPUTDIR}"
cd "${OUTPUTDIR}"

bowtie2 -p 8 \
  -x "${INDEX}" \
  -1 "/home/mjd356/virome_project/reads/trimmed/SRR6996005_forward_paired.fastq.gz" \
  -2 "/home/mjd356/virome_project/reads/trimmed/SRR6996005_reverse_paired.fastq.gz" \
| samtools view -bS - > "${SAMPLE}.bam"

samtools sort "${SAMPLE}.bam" > "${SAMPLE}_sorted.bam"
samtools index "${SAMPLE}_sorted.bam"
```
--------------------
```bash
sbatch bowtie2_run.slurm
```
--------------------

## Step 10: Viral abundance and diversity analysis in R

The final abundance table is used to calculate ecological diversity metrics such as richness and Shannon diversity. Richness measures how many vOTUs are detected in a sample, while Shannon diversity also accounts for how evenly abundance is distributed across vOTUs. These analyses help connect the sequencing workflow to ecological interpretation.

Link to workflow

https://rpubs.com/skar/1414271












