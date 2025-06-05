# 🧬 RNA-seq Analysis Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20WSL-blue.svg)](https://github.com/microsoft/WSL)
[![RNA-seq](https://img.shields.io/badge/Analysis-RNA--seq-green.svg)](https://en.wikipedia.org/wiki/RNA-Seq)

A comprehensive, production-ready RNA-seq analysis pipeline featuring **STAR alignment**, **RSEM quantification**, and **automated result merging**. Designed for high-throughput processing with enhanced quality control and beautiful output formatting.

## ✨ Features

- 🚀 **High-Performance Alignment**: STAR aligner with optimized parameters
- 🧬 **Dual Quantification**: Both RSEM and featureCounts for comprehensive analysis
- 📊 **Quality Control**: Extensive QC metrics and reporting
- 🔄 **Batch Processing**: Automated processing of multiple samples
- 📈 **Result Merging**: Beautiful CSV outputs ready for downstream analysis
- 🎯 **Deduplication Support**: Optional PCR duplicate removal
- 🛡️ **Error Handling**: Robust error checking and recovery
- 📱 **Progress Tracking**: Real-time progress bars and status updates

## 🎯 Pipeline Overview

```mermaid
graph LR
    A[📥 FASTQ Files] --> B[⚙️ Install Dependencies]
    B --> C[📚 Index Genome]
    C --> D[🧬 STAR Alignment]
    D --> E[📊 RSEM Quantification]
    E --> F[📈 Result Merging]
    F --> G[📋 Final CSV Reports]
    
    style A fill:#e1f5fe
    style B fill:#fff3e0
    style C fill:#f3e5f5
    style D fill:#e8f5e8
    style E fill:#fff9c4
    style F fill:#fce4ec
    style G fill:#e0f2f1
```

## 📋 Requirements

### System Requirements
- **OS**: Linux or Windows Subsystem for Linux (WSL)
- **RAM**: Minimum 16GB (32GB+ recommended for large genomes)
- **Storage**: ~50GB for mouse genome indices, ~100GB for human
- **CPU**: Multi-core processor (8+ cores recommended)

### Software Dependencies
- **Conda/Miniconda**: Package manager
- **STAR**: RNA-seq aligner
- **RSEM**: RNA-seq quantification
- **Samtools**: BAM file processing
- **Python 3.9+**: For result merging scripts

## 🚀 Quick Start

### Step 1: Clone Repository
```bash
git clone https://github.com/yourusername/rnaseq-pipeline.git
cd rnaseq-pipeline
chmod +x *.sh
```

### Step 2: Install Dependencies
```bash
./install.sh
```

**What this does:**
- Creates conda environment `rnaseq`
- Installs STAR, RSEM, samtools, and QC tools
- Verifies all installations

### Step 3: Download Reference Files
```bash
# Example for mouse genome (GRCm39)
mkdir -p /path/to/genome
cd /path/to/genome

# Download genome FASTA
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M37/GRCm39.genome.fa.gz
gunzip GRCm39.genome.fa.gz

# Download GTF annotation
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M37/gencode.vM37.chr_patch_hapl_scaff.annotation.gtf.gz
gunzip gencode.vM37.chr_patch_hapl_scaff.annotation.gtf.gz
```

### Step 4: Create Genome Indices
```bash
# Edit paths in index_genome.sh first!
nano index_genome.sh  # Update GENOME_DIR, GENOME_FASTA, GTF_FILE paths

# Run indexing (takes 30-60 minutes)
./index_genome.sh
```

### Step 5: Prepare Your Data
Organize your FASTQ files:
```
input_directory/
├── sample1_R1_001.fastq.gz
├── sample1_R2_001.fastq.gz
├── sample2_R1_001.fastq.gz
├── sample2_R2_001.fastq.gz
└── ...
```

### Step 6: Run Analysis Pipeline
```bash
# Edit paths in run.sh first!
nano run.sh  # Update input/output directories and index paths

# Activate conda environment
conda activate rnaseq

# Run the pipeline
./run.sh
```

### Step 7: Merge Results
```bash
# Edit paths in results.py first!
nano results.py  # Update results directory and GTF file paths

# Run result merger
python results.py
```

## 📁 Directory Structure

```
rnaseq-pipeline/
├── 📜 README.md                    # This file
├── ⚙️ install.sh                   # Install dependencies
├── 📚 index_genome.sh             # Create genome indices
├── 🧬 run.sh                      # Main analysis pipeline
└── 📊 results.py                  # Merge results into CSV files

```

## 🔧 Configuration

### Key Parameters to Modify

**index_genome.sh:**
```bash
GENOME_DIR="/path/to/your/genome"
GENOME_FASTA="$GENOME_DIR/genome.fa"
GTF_FILE="$GENOME_DIR/annotation.gtf"
```

**run.sh:**
```bash
INPUT_DIR="/path/to/fastq/files"
OUTPUT_DIR="/path/to/results"
STAR_INDEX_DIR="/path/to/STAR_index"
RSEM_INDEX_PREFIX="/path/to/RSEM_index/prefix"
```

**results.py:**
```python
results_dir = "/path/to/results"
gtf_file_path = "/path/to/annotation.gtf"
output_dir = "/path/to/merged_results"
```

## 📊 Output Files

### Pipeline Outputs (per sample)
- `sample_Aligned.sortedByCoord.out.bam` - Aligned reads
- `sample_rsem.genes.results` - Gene-level quantification
- `sample_rsem.isoforms.results` - Isoform-level quantification
- `sample_featureCounts.txt` - Alternative gene counts
- `sample_Log.final.out` - Alignment statistics

### Merged Results
- `1_merged_star_gene_counts.csv` - 🌟 STAR gene counts matrix
- `2_merged_rsem_gene_expected_counts.csv` - 🧬 RSEM gene counts
- `3_merged_rsem_isoform_expected_counts.csv` - 🧬 RSEM isoform counts
- `4_merged_rsem_isoform_percentages.csv` - 📊 Isoform percentage usage
- `5_merged_rsem_gene_fpkm.csv` - 📈 Gene FPKM values
- `6_merged_rsem_gene_tpm.csv` - 📈 Gene TPM values
- `diagnostic_sample_completeness.csv` - 🔍 Quality diagnostic report

## 🎛️ Advanced Options

### Enable Deduplication
The pipeline supports PCR duplicate removal:
```bash
# During run.sh execution, choose option 1 when prompted:
# "Do you want to perform deduplication? (1/2): 1"
```

### Performance Tuning
Adjust thread counts in `run.sh`:
```bash
STAR_THREADS=24          # STAR alignment threads
RSEM_THREADS=24          # RSEM quantification threads
BAM_SORT_RAM=64          # RAM for BAM sorting (GB)
```

### Strand Specificity
Configure strand-specific protocols:
```bash
STRAND_SPECIFICITY="unstranded"  # Options: unstranded, forward, reverse
```

## 🔍 Quality Control

### Built-in QC Features
- **Alignment Statistics**: Mapping rates, multi-mappers, unmapped reads
- **Quantification Metrics**: Gene detection rates, isoform diversity
- **Sample Completeness**: Cross-sample data quality assessment
- **Deduplication Reports**: PCR duplicate rates (if enabled)

### QC Output Locations
```
qc_reports/
├── sample1/
│   ├── sample1_star.log
│   ├── sample1_rsem.log
│   ├── sample1_flagstat.txt
│   └── sample1_bamstats.txt
└── ...
```

## 🚨 Troubleshooting

### Common Issues

**"No samples found" Error:**
```bash
# Ensure FASTQ files follow naming convention:
*_R1_001.fastq.gz and *_R2_001.fastq.gz
```

**Memory Issues:**
```bash
# Reduce memory usage in index_genome.sh:
MAX_RAM=8000000000  # 8GB instead of 15GB
```

**Permission Denied:**
```bash
# Make scripts executable:
chmod +x *.sh
```

**Conda Environment Issues:**
```bash
# Reset environment:
conda remove -n rnaseq --all
./install.sh
```

## 📈 Performance Benchmarks

| Genome | Samples | Time | RAM Usage | Storage |
|--------|---------|------|-----------|---------|
| Mouse (GRCm39) | 10 | ~4 hours | 32GB | ~200GB |
| Human (GRCh38) | 10 | ~8 hours | 64GB | ~400GB |
| Arabidopsis | 10 | ~1 hour | 16GB | ~50GB |

*Benchmarks on 24-core system with NVMe storage*

## 🙏 Acknowledgments

- **STAR**: Alexander Dobin et al. - RNA-STAR aligner
- **RSEM**: Bo Li and Colin Dewey - RSEM quantification
- **GENCODE**: Genome annotation consortium
- **Community**: RNA-seq analysis community for best practices

## 🔗 Related Projects

- [nf-core/rnaseq](https://github.com/nf-core/rnaseq) - Nextflow RNA-seq pipeline
- [STAR](https://github.com/alexdobin/STAR) - Original STAR aligner
- [RSEM](https://github.com/deweylab/RSEM) - Original RSEM quantifier

---

<div align="center">

**⭐ If this pipeline helped your research, please give it a star! ⭐**

Made with ❤️ for the community

</div>
