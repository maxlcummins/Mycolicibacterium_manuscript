# Mycolicibacterium Manuscript Resources

This repository contains scripts and data for analyzing Mycolicibacterium genomic data, focusing on phylogenetic relationships and genome similarity comparisons.

## 📂 Repository Contents

| File | Description |
|------|-------------|
| `workflow.md` | Step-by-step bioinformatics pipeline documentation |
| `assemblies.txt` | List of genome assemblies used in analyses |
| `tree.dnd` | Phylogenetic tree file generated by Mashtree |
| `mashdistance.tsv` | Distance matrix from Mashtree analysis |
| `fastANI_output.txt` | Average Nucleotide Identity (ANI) comparison results |
| `genomes_for_inclusion.txt` | List of reference genomes included in analyses |
| `species_ids.txt` | Mapping of genome IDs to species names |
| `Treeplot.R` | R script to generate phylogenetic tree visualization with ANI heatmap |
| `Supplementary_table_1` | Table containing ANI values and species information |

## 🔧 Requirements

- **Conda** for environment management
- The following **Conda environments**:
    - bactopia
    - mashtree
    - fastani
    - hostile
- **R** with the following packages:
    - tidyverse
    - ggtree
    - ggtreeExtra

## 🧬 Workflow Overview

### 1. Data Preparation
- Pool technical replicates
- Create sample sheets for Bactopia

### 2. Genome Assembly and Quality Control
- Run Bactopia with SPAdes assembler
- Execute Bracken for taxonomic classification
- Perform CheckM for assembly quality assessment

### 3. Phylogenetic Analysis
- Collect assemblies in `assemblies/` directory
- Generate Mashtree (`tree.dnd`) and distance matrix (`mashdistance.tsv`)

### 4. Genome Similarity Analysis
- Calculate Average Nucleotide Identity with FastANI

### 5. Read Clean-up
- Use Hostile for dehosting sequencing reads before genome upload

### 6. Visualization and Analysis
- Use `Treeplot.R` to generate an integrated figure showing:
    - Phylogenetic tree with species labels
    - Color-coded ANI heatmap displaying genome similarity
- Create supplementary table with ANI values between genomes

## 📋 Usage

1. Clone this repository:
```bash
git clone https://github.com/maxlcummins/Mycolicibacterium_manuscript.git
cd Mycolicibacterium_manuscript
```

2. Follow the steps in `workflow.md` to replicate analyses or adapt for your own data.

### Key commands:

**Generate phylogenetic tree:**
```bash
conda activate mashtree
mashtree assemblies/*fna --mindepth 0 --outmatrix mashdistance.tsv > tree.dnd
```

**Calculate ANI values:**
```bash
conda activate fastani
fastANI -t 16 --ql assemblies.txt --rl assemblies.txt -o fastANI_output.txt
```

**Generate visualization:**
```bash
Rscript Treeplot.R
```