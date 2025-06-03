#!/bin/bash

echo "Complete STAR and RSEM Genome Indexing"
echo "======================================"

# Configuration - MODIFY THESE PATHS AS NEEDED
GENOME_DIR="/mnt/d/genome"
GENOME_FASTA="$GENOME_DIR/GRCm39.genome.fa"
GTF_FILE="$GENOME_DIR/gencode.vM37.chr_patch_hapl_scaff.annotation.gtf"

# Output directories
STAR_INDEX_DIR="$GENOME_DIR/STAR_index"
RSEM_INDEX_DIR="$GENOME_DIR/RSEM_index"
RSEM_INDEX_PREFIX="$RSEM_INDEX_DIR/mouse_rsem"

# WSL-optimized settings
THREADS=8
MAX_RAM=15000000000  # 15GB for WSL stability

echo "Configuration:"
echo "Genome FASTA: $GENOME_FASTA"
echo "GTF file: $GTF_FILE"
echo "STAR index output: $STAR_INDEX_DIR"
echo "RSEM index output: $RSEM_INDEX_PREFIX"
echo "Threads: $THREADS"
echo "Max RAM: $(($MAX_RAM / 1000000000))GB"
echo ""

# WSL environment setup
export TMPDIR=/tmp
ulimit -n 65536
ulimit -v unlimited

# Validate input files
echo "=== INPUT VALIDATION ==="
echo "Validating input files..."
if [[ ! -f "$GENOME_FASTA" ]]; then
    echo "ERROR: Genome FASTA file not found: $GENOME_FASTA"
    exit 1
fi

if [[ ! -f "$GTF_FILE" ]]; then
    echo "ERROR: GTF file not found: $GTF_FILE"
    exit 1
fi

if [[ ! -r "$GENOME_FASTA" ]] || [[ ! -r "$GTF_FILE" ]]; then
    echo "ERROR: Cannot read input files"
    exit 1
fi

echo "✓ Input files validated"

# Validate GTF format
echo "Validating GTF format..."
gtf_lines=$(head -1000 "$GTF_FILE" | grep -v "^#" | wc -l)
if [ $gtf_lines -eq 0 ]; then
    echo "ERROR: GTF file appears to be empty or only contains headers"
    exit 1
fi

# Check for required GTF attributes
if ! head -1000 "$GTF_FILE" | grep -q 'gene_id\|transcript_id'; then
    echo "ERROR: GTF file missing required gene_id or transcript_id attributes"
    exit 1
fi

echo "✓ GTF format validated"

# Display file information
fasta_size=$(du -h "$GENOME_FASTA" | cut -f1)
gtf_size=$(du -h "$GTF_FILE" | cut -f1)
available_mem=$(free -g | awk 'NR==2{printf "%.0f", $7}')

echo "File sizes: Genome=${fasta_size}, GTF=${gtf_size}"
echo "Available memory: ${available_mem}GB"
echo ""

# Install required software
echo "=== SOFTWARE INSTALLATION ==="
echo "Installing required software..."
conda install -c bioconda star=2.7.10a rsem bowtie samtools -y

# Verify installations
echo "Verifying software installations..."
if ! command -v STAR &> /dev/null; then
    echo "ERROR: STAR installation failed"
    exit 1
fi

if ! command -v rsem-prepare-reference &> /dev/null; then
    echo "ERROR: RSEM installation failed"
    exit 1
fi

if ! command -v bowtie &> /dev/null; then
    echo "ERROR: Bowtie installation failed"
    exit 1
fi

STAR_VERSION=$(STAR --version 2>/dev/null)
RSEM_VERSION=$(rsem-calculate-expression --version 2>&1 | head -1)
BOWTIE_VERSION=$(bowtie --version 2>&1 | head -1)

echo "✓ STAR version: $STAR_VERSION"
echo "✓ RSEM version: $RSEM_VERSION"
echo "✓ Bowtie version: $BOWTIE_VERSION"
echo ""

# Create output directories
echo "Preparing output directories..."
mkdir -p "$STAR_INDEX_DIR"
mkdir -p "$RSEM_INDEX_DIR"

# Step 1: Create STAR index
echo ""
echo "=== STEP 1: STAR INDEX CREATION ==="
if [[ -f "$STAR_INDEX_DIR/genomeParameters.txt" ]]; then
    echo "STAR index already exists, skipping creation..."
    echo "✓ STAR index found - Size: $(du -sh "$STAR_INDEX_DIR" | cut -f1)"
else
    echo "Creating STAR genome index..."
    echo "This will take 30-60 minutes..."
    echo "Started at: $(date)"

    start_time=$(date +%s)

    STAR --runThreadN $THREADS \
         --runMode genomeGenerate \
         --genomeDir "$STAR_INDEX_DIR" \
         --genomeFastaFiles "$GENOME_FASTA" \
         --sjdbGTFfile "$GTF_FILE" \
         --sjdbOverhang 149 \
         --limitGenomeGenerateRAM $MAX_RAM \
         --genomeSAindexNbases 14 \
         --genomeChrBinNbits 18 \
         --limitIObufferSize 50000000 50000000 \
         --outTmpDir "/tmp/STARtmp_$$"

    star_exit_code=$?
    star_time=$(($(date +%s) - start_time))

    # Clean up temp directory
    rm -rf "/tmp/STARtmp_$$" 2>/dev/null

    if [ $star_exit_code -ne 0 ] || [[ ! -f "$STAR_INDEX_DIR/genomeParameters.txt" ]]; then
        echo "ERROR: STAR index creation failed"
        exit 1
    fi

    echo "✓ STAR index created successfully in ${star_time} seconds"
    echo "Index size: $(du -sh "$STAR_INDEX_DIR" | cut -f1)"
fi

# Step 2: Create RSEM index
echo ""
echo "=== STEP 2: RSEM INDEX CREATION ==="

# Check if RSEM index already exists
if [[ -f "${RSEM_INDEX_PREFIX}.grp" && -f "${RSEM_INDEX_PREFIX}.ti" && -f "${RSEM_INDEX_PREFIX}.transcripts.fa" ]]; then
    echo "RSEM index already exists!"
    echo "Existing files:"
    ls -la "${RSEM_INDEX_PREFIX}".*
    echo ""
    read -p "Do you want to overwrite RSEM index? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Keeping existing RSEM index..."
    else
        echo "Removing existing RSEM index..."
        rm -f "${RSEM_INDEX_PREFIX}".*
    fi
fi

# Create RSEM index if not exists or user chose to overwrite
if [[ ! -f "${RSEM_INDEX_PREFIX}.grp" ]]; then
    echo "Creating RSEM reference index..."
    echo "This will take 20-45 minutes..."
    echo "Started at: $(date)"

    start_time=$(date +%s)

    # Create RSEM index with proper error handling
    echo "Running rsem-prepare-reference..."
    
    # Method 1: With Bowtie indices (full functionality)
    rsem-prepare-reference \
        --gtf "$GTF_FILE" \
        --bowtie \
        --num-threads $THREADS \
        "$GENOME_FASTA" \
        "$RSEM_INDEX_PREFIX"

    rsem_exit_code=$?
    rsem_time=$(($(date +%s) - start_time))

    echo ""
    echo "RSEM finished with exit code: $rsem_exit_code"
    echo "Duration: ${rsem_time} seconds"

    # Check if RSEM succeeded by looking for essential files
    if [ $rsem_exit_code -eq 0 ] && [[ -f "${RSEM_INDEX_PREFIX}.grp" ]]; then
        echo "✓ RSEM index created successfully in ${rsem_time} seconds"
        echo "Index size: $(du -sh "$RSEM_INDEX_DIR" | cut -f1)"
    else
        echo "! RSEM indexing encountered issues (exit code: $rsem_exit_code)"
        echo "Checking for essential RSEM files..."
        
        # Check which files were created
        echo "Files created in RSEM directory:"
        ls -la "$RSEM_INDEX_DIR"/ 2>/dev/null || echo "No files found"
        
        # Check for minimum required files
        if [[ -f "${RSEM_INDEX_PREFIX}.grp" ]]; then
            echo "✓ RSEM group file exists - basic functionality available"
        else
            echo "✗ RSEM group file missing - trying fallback method..."
            
            # Alternative: Try without bowtie if the above fails
            echo "Trying RSEM without bowtie indices (reference files only)..."
            
            rsem-prepare-reference \
                --gtf "$GTF_FILE" \
                --num-threads $THREADS \
                "$GENOME_FASTA" \
                "$RSEM_INDEX_PREFIX"
            
            if [[ -f "${RSEM_INDEX_PREFIX}.grp" ]]; then
                echo "✓ RSEM reference files created without bowtie indices"
                echo "Note: You can still use STAR+RSEM workflow"
            else
                echo "✗ RSEM reference creation failed completely"
                exit 1
            fi
        fi
    fi
else
    echo "✓ RSEM index found - Size: $(du -sh "$RSEM_INDEX_DIR" | cut -f1)"
fi

# Verification
echo ""
echo "============================================"
echo "INDEX CREATION COMPLETE"
echo "============================================"

echo "STAR index: $STAR_INDEX_DIR"
if [[ -f "$STAR_INDEX_DIR/genomeParameters.txt" ]]; then
    echo "✓ STAR index verified - $(du -sh "$STAR_INDEX_DIR" | cut -f1)"
else
    echo "✗ STAR index missing"
    exit 1
fi

echo ""
echo "RSEM index: $RSEM_INDEX_PREFIX"
if [[ -f "${RSEM_INDEX_PREFIX}.grp" ]]; then
    echo "✓ RSEM index verified - $(du -sh "$RSEM_INDEX_DIR" | cut -f1)"
    
    # Check what files are available
    echo ""
    echo "Available RSEM files:"
    [[ -f "${RSEM_INDEX_PREFIX}.grp" ]] && echo "  ✓ Group file (.grp)"
    [[ -f "${RSEM_INDEX_PREFIX}.ti" ]] && echo "  ✓ Transcript info (.ti)"
    [[ -f "${RSEM_INDEX_PREFIX}.transcripts.fa" ]] && echo "  ✓ Transcripts FASTA"
    [[ -f "${RSEM_INDEX_PREFIX}.seq" ]] && echo "  ✓ Sequence file (.seq)"
    [[ -f "${RSEM_INDEX_PREFIX}.chrlist" ]] && echo "  ✓ Chromosome list"
    [[ -f "${RSEM_INDEX_PREFIX}.idx.fa" ]] && echo "  ✓ Index FASTA"
    
    # Check for bowtie indices
    if [[ -f "${RSEM_INDEX_PREFIX}.1.ebwt" && -s "${RSEM_INDEX_PREFIX}.1.ebwt" ]]; then
        echo "  ✓ Bowtie indices available (full RSEM functionality)"
        RSEM_FUNCTIONALITY="FULL"
    else
        echo "  ! Bowtie indices not available (STAR+RSEM workflow still works)"
        RSEM_FUNCTIONALITY="BASIC"
    fi
    
else
    echo "! RSEM index incomplete - will use STAR-only quantification"
    RSEM_FUNCTIONALITY="NONE"
fi

echo ""
echo "=== KEY INDEX FILES ==="
echo "STAR files:"
ls -la "$STAR_INDEX_DIR"/genomeParameters.txt "$STAR_INDEX_DIR"/SA "$STAR_INDEX_DIR"/Genome 2>/dev/null

echo ""
echo "RSEM files:"
ls -la "${RSEM_INDEX_PREFIX}".grp "${RSEM_INDEX_PREFIX}".ti "${RSEM_INDEX_PREFIX}".transcripts.fa 2>/dev/null | head -3

echo ""
echo "=== SUMMARY ==="
total_time=$(date)
echo "Completed at: $total_time"

# Calculate total sizes
if [[ -d "$STAR_INDEX_DIR" ]]; then
    star_size=$(du -sh "$STAR_INDEX_DIR" | cut -f1)
else
    star_size="N/A"
fi

if [[ -d "$RSEM_INDEX_DIR" ]]; then
    rsem_size=$(du -sh "$RSEM_INDEX_DIR" | cut -f1)
else
    rsem_size="N/A"
fi

echo "Index sizes: STAR=${star_size}, RSEM=${rsem_size}"

# Status summary
echo ""
echo "🎯 INDEX STATUS:"
echo "✓ STAR: Ready for alignment"

case $RSEM_FUNCTIONALITY in
    "FULL")
        echo "✓ RSEM: Full functionality (with Bowtie indices)"
        echo "✓ Recommended workflow: STAR+RSEM"
        ;;
    "BASIC")
        echo "✓ RSEM: Basic functionality (without Bowtie indices)"
        echo "✓ Recommended workflow: STAR+RSEM (transcriptome alignment)"
        ;;
    "NONE")
        echo "! RSEM: Not available"
        echo "✓ Alternative: Use STAR GeneCounts for quantification"
        ;;
esac

echo ""
echo "🚀 Ready for batch alignment!"

# Next steps
echo ""
echo "=== NEXT STEPS ==="
echo "1. Use your batch alignment script with these indices"
echo "2. STAR index: $STAR_INDEX_DIR"
echo "3. RSEM index: $RSEM_INDEX_PREFIX"
echo ""
echo "Example alignment command:"
echo "STAR --genomeDir $STAR_INDEX_DIR --quantMode TranscriptomeSAM ..."
echo "rsem-calculate-expression --alignments ... $RSEM_INDEX_PREFIX ..."
echo ""
echo "✅ All indices are ready for RNA-seq analysis!"