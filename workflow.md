# Two technical replicates were pooled
cat *_R1_*.fastq.gz > ../pooled_reads/TPG_MSP_01_R1.fastq.gz
cat *_R2_*.fastq.gz > ../pooled_reads/TPG_MSP_01_R2.fastq.gz

# Create a sample sheet
conda activate bactopia

bactopia prepare --path pooled_reads --fastq-ext '.fastq.gz' > bactopia_sample_sheet.txt

# Run Bactopia
bactopia \
    --samples bactopia_sample_sheet.txt \
    -profile docker \
    --shovill_assembler spades \
    --ask_merlin

bactopia \
    --wf bracken \
    --kraken2_db /home/maxcu/RDH/databases/k2_plus_pf \
    -profile docker \
    --bactopia bactopia

bactopia \
    --wf checkm \
    -profile docker \
    --bactopia bactopia

# Make a directory for tree_manuscript
mkdir -p tree_manuscript

# Copy the reads to the target dir
cp bactopia/TPG_MSP_01/main/qc/*.gz tree_manuscript/

# Enter the directory
cd tree_manuscript

# Make a directory for our assemblies
mkdir -p assemblies

# Copy our assembly
cp ../bactopia/TPG_MSP_01/main/assembler/TPG_MSP_01.fna.gz assemblies

# Copy genomes from our list to the target dir
# Note the list contains representative reseq genomes for clusters as per BacSort N50 designation
for genome in $(cat genomes_for_inclusion.txt | cut -f2 | sort -u); do
    cp /mnt/d/sequence_runs/2025/March/NTM_25012024/bacsort/Bacsort/bacsort_results/assemblies/Mycolicibacterium/${genome} assemblies/
done

# Activate mashtree
conda activate mashtree

# Copy our own genomes
cp snp_analyses/fastani/fastas/D5313916.fna tree_manuscript/assemblies/

# Copy our reference genomes of interest
cp snp_analyses/fastani/*.fasta tree_manuscript/assemblies/

# Rename them
mv assemblies/AP023287.fasta assemblies/AP023287.fna
mv assemblies/NZ_AP022586.fasta assemblies/NZ_AP022586.fna

# Make a mashtree of the assemblies
mashtree assemblies/*fna --mindepth 0 --outmatrix mashdistance.tsv> tree.dnd

# Write the filenames to a file
ls assemblies/*fna > assemblies.txt

# Deactivate mashtree
conda deactivate

# Activate fastani
conda activate fastani

# Make a directory for our fastani results
mkdir -p fastani

# Run fastANI
fastANI -t 16 --ql assemblies.txt --rl assemblies.txt -o fastANI_output.txt

## Prepare reads for upload

# Activate hostile for dehosting
conda activate hostile

# Run hostile to dehost the reads
hostile clean --fastq1 TPG_MSP_01_R1.fastq.gz --fastq2 TPG_MSP_01_R2.fastq.gz --index ~/RDH/databases/hostile/human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401 --threads 16 --out-dir TPG_MSP_01_hostile --force

```
11:54:54 INFO: Hostile version 1.1.0. Mode: paired short read (Bowtie2)
11:54:54 INFO: Found custom index /home/maxcu/RDH/databases/hostile/human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401
11:54:54 INFO: Cleaning…
11:55:00 INFO: Cleaning complete
[
    {
        "version": "1.1.0",
        "aligner": "bowtie2",
        "index": "/home/maxcu/RDH/databases/hostile/human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401",
        "options": [],
        "fastq1_in_name": "TPG_MSP_01_R1.fastq.gz",
        "fastq1_in_path": "/mnt/d/sequence_runs/2025/March/NTM_25012024/tree_manuscript/TPG_MSP_01_R1.fastq.gz",
        "fastq1_out_name": "TPG_MSP_01_R1.clean_1.fastq.gz",
        "fastq1_out_path": "TPG_MSP_01_hostile/TPG_MSP_01_R1.clean_1.fastq.gz",
        "reads_in": 1780494,
        "reads_out": 1779550,
        "reads_removed": 944,
        "reads_removed_proportion": 0.00053,
        "fastq2_in_name": "TPG_MSP_01_R2.fastq.gz",
        "fastq2_in_path": "/mnt/d/sequence_runs/2025/March/NTM_25012024/tree_manuscript/TPG_MSP_01_R2.fastq.gz",
        "fastq2_out_name": "TPG_MSP_01_R2.clean_2.fastq.gz",
        "fastq2_out_path": "TPG_MSP_01_hostile/TPG_MSP_01_R2.clean_2.fastq.gz"
    }
]
```