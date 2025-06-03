#!/bin/bash

echo "Installing STAR and RSEM for RNA-seq analysis"
echo "=============================================="

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda not found. Please install Miniconda or Anaconda first."
    exit 1
fi

# Create or activate rnaseq environment
echo "Setting up conda environment..."
if conda env list | grep -q "rnaseq"; then
    echo "rnaseq environment already exists, activating..."
    source activate rnaseq
else
    echo "Creating new rnaseq environment..."
    conda create -n rnaseq python=3.9 -y
    source activate rnaseq
fi

echo "Installing required software packages..."

# Install STAR
echo "Installing STAR aligner..."
conda install -c bioconda star -y

# Install RSEM
echo "Installing RSEM..."
conda install -c bioconda rsem -y

# Install samtools (required for BAM processing)
echo "Installing samtools..."
conda install -c bioconda samtools -y

# Install additional useful tools
echo "Installing additional tools..."
conda install -c bioconda fastqc multiqc -y

echo ""
echo "Verifying installations..."

# Check STAR
if command -v STAR &> /dev/null; then
    star_version=$(STAR --version 2>/dev/null)
    echo "✓ STAR installed: $star_version"
else
    echo "✗ STAR installation failed"
    exit 1
fi

# Check RSEM
if command -v rsem-calculate-expression &> /dev/null; then
    rsem_version=$(rsem-calculate-expression --version 2>&1 | head -1)
    echo "✓ RSEM installed: $rsem_version"
else
    echo "✗ RSEM installation failed"
    exit 1
fi

# Check samtools
if command -v samtools &> /dev/null; then
    samtools_version=$(samtools --version | head -1)
    echo "✓ Samtools installed: $samtools_version"
else
    echo "✗ Samtools installation failed"
    exit 1
fi

echo ""
echo "Installation complete!"
echo "To activate this environment in the future, run:"
echo "conda activate rnaseq"
echo ""
echo "Next step: Run the genome indexing script"
