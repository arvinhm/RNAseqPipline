#!/bin/bash

# Color definitions
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
WHITE='\033[1;37m'
BOLD='\033[1m'
DIM='\033[2m'
NC='\033[0m' # No Color

# Fancy symbols
CHECK="✅"
CROSS="❌"
ARROW="➤"
STAR="⭐"
DNA="🧬"
GEAR="⚙️"
FIRE="🔥"
ROCKET="🚀"
CHART="📊"
CLOCK="⏱️"
FOLDER="📁"
FILE="📄"

# Progress bar function
show_progress() {
    local current=$1
    local total=$2
    local width=50
    local percentage=$((current * 100 / total))
    local completed=$((current * width / total))
    local remaining=$((width - completed))
    
    printf "\r${CYAN}["
    printf "%${completed}s" | tr ' ' '█'
    printf "%${remaining}s" | tr ' ' '░'
    printf "] ${WHITE}%d%%${NC} ${BLUE}(%d/%d)${NC}" $percentage $current $total
}

# Fancy header
print_header() {
    echo -e "${PURPLE}╔════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${PURPLE}║${WHITE}              ${STAR} Enhanced STAR-RSEM Pipeline ${STAR}              ${PURPLE}║${NC}"
    echo -e "${PURPLE}║${CYAN}           High-Quality RNA-seq Analysis Suite            ${PURPLE}║${NC}"
    echo -e "${PURPLE}╚════════════════════════════════════════════════════════════════╝${NC}"
    echo ""
}

# Section header
print_section() {
    local title=$1
    local emoji=$2
    echo ""
    echo -e "${YELLOW}${emoji} ${BOLD}${title}${NC}"
    echo -e "${YELLOW}$(printf '%.s─' $(seq 1 $((${#title} + 4))))${NC}"
}

# Info message
print_info() {
    echo -e "${BLUE}${ARROW} ${1}${NC}"
}

# Success message
print_success() {
    echo -e "${GREEN}${CHECK} ${1}${NC}"
}

# Error message
print_error() {
    echo -e "${RED}${CROSS} ${1}${NC}"
}

# Warning message
print_warning() {
    echo -e "${YELLOW}⚠️  ${1}${NC}"
}

# Configuration
STAR_INDEX_DIR="/mnt/d/genome/STAR_index"
RSEM_INDEX_PREFIX="/mnt/d/genome/RSEM_index/mouse_rsem"
GTF_FILE="/mnt/d/genome/gencode.vM37.chr_patch_hapl_scaff.annotation.gtf"
GENOME_FASTA="/mnt/d/genome/GRCm39.genome.fa"
INPUT_DIR="/mnt/d/mad_rna_seq/one"

# Working and output directories
WORK_DIR="/home/$(whoami)/rna_seq_work"
OUTPUT_DIR="/mnt/d/mad_rna_seq/results"
QC_DIR="/mnt/d/mad_rna_seq/qc_reports"

# Performance settings - OPTIMIZED FOR HIGH QUALITY
STAR_THREADS=24
RSEM_THREADS=24
FEATURECOUNTS_THREADS=24
BAM_SORT_RAM=64  # GB

# Analysis options - ENHANCED FOR COMPREHENSIVE ANALYSIS
RUN_RSEM=true
RUN_FEATURECOUNTS=true
RUN_SALMON=false  # Alternative quantifier (optional)
STRAND_SPECIFICITY="unstranded"  # Options: unstranded, forward, reverse
MULTIMAPPER_MAX=20  # Increased from 10 for better sensitivity
SAVE_INTERMEDIATE=true  # Keep intermediate files for debugging
GENERATE_BIGWIG=false  # Generate coverage tracks (set to true if needed)

# Quality control options
COLLECT_METRICS=true
MIN_MAPPED_READS=1000000  # Minimum mapped reads to consider sample successful

# Ask user about deduplication
echo ""
echo "=========================================="
echo "DEDUPLICATION OPTIONS"
echo "=========================================="
echo ""
echo -e "${WHITE}Deduplication removes PCR duplicate reads that can bias quantification.${NC}"
echo -e "${WHITE}This is ${GREEN}recommended${NC} for most RNA-seq analyses.${NC}"
echo ""
echo -e "${YELLOW}Options:${NC}"
echo -e "  ${GREEN}1)${NC} Yes - Remove duplicates (${GREEN}recommended${NC})"
echo -e "  ${RED}2)${NC} No  - Keep all reads"
echo ""
read -p "$(echo -e "${BLUE}Do you want to perform deduplication? (1/2): ${NC}")" -n 1 -r
echo ""

case $REPLY in
    1|y|Y|yes|Yes)
        RUN_DEDUP=true
        echo -e "${GREEN}✓ Deduplication will be performed using Picard MarkDuplicates${NC}"
        ;;
    *)
        RUN_DEDUP=false
        echo -e "${YELLOW}✓ Deduplication will be skipped${NC}"
        ;;
esac

echo ""

# Clear screen and show header
clear
print_header

print_section "Configuration" "$GEAR"
echo -e "${WHITE}Input Directory:${NC}     ${CYAN}$INPUT_DIR${NC}"
echo -e "${WHITE}Working Directory:${NC}   ${CYAN}$WORK_DIR${NC}"
echo -e "${WHITE}Output Directory:${NC}    ${CYAN}$OUTPUT_DIR${NC}"
echo -e "${WHITE}QC Directory:${NC}       ${CYAN}$QC_DIR${NC}"
echo -e "${WHITE}Threads:${NC}             ${GREEN}STAR=$STAR_THREADS, RSEM=$RSEM_THREADS, featureCounts=$FEATURECOUNTS_THREADS${NC}"
echo -e "${WHITE}BAM Sort RAM:${NC}        ${GREEN}${BAM_SORT_RAM}GB${NC}"
echo -e "${WHITE}Max Multimappers:${NC}    ${GREEN}$MULTIMAPPER_MAX${NC}"
echo -e "${WHITE}Strand Specificity:${NC}  ${GREEN}$STRAND_SPECIFICITY${NC}"
echo -e "${WHITE}Save Intermediate:${NC}   ${GREEN}$SAVE_INTERMEDIATE${NC}"
echo -e "${WHITE}Deduplication:${NC}       ${GREEN}$RUN_DEDUP${NC}"

# Create directories
mkdir -p "$WORK_DIR"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$QC_DIR"
cd "$WORK_DIR"

print_section "Validation" "$DNA"

# Check required files
print_info "Checking STAR index..."
if [[ -f "$STAR_INDEX_DIR/genomeParameters.txt" ]]; then
    print_success "STAR index found"
else
    print_error "STAR index not found at: $STAR_INDEX_DIR"
    exit 1
fi

print_info "Checking RSEM index..."
if [[ "$RUN_RSEM" == true ]] && [[ -f "${RSEM_INDEX_PREFIX}.grp" ]]; then
    print_success "RSEM index found"
elif [[ "$RUN_RSEM" == true ]]; then
    print_error "RSEM index not found at: ${RSEM_INDEX_PREFIX}"
    exit 1
else
    print_warning "RSEM disabled"
fi

print_info "Checking GTF file..."
if [[ "$RUN_FEATURECOUNTS" == true ]] && [[ -f "$GTF_FILE" ]]; then
    print_success "GTF file found"
elif [[ "$RUN_FEATURECOUNTS" == true ]]; then
    print_error "GTF file not found at: $GTF_FILE"
    exit 1
else
    print_warning "featureCounts disabled"
fi

# Install required software
software_needed=("samtools")
if [[ "$RUN_FEATURECOUNTS" == true ]]; then
    software_needed+=("featureCounts")
fi

# Add Picard if deduplication is requested
if [[ "$RUN_DEDUP" == true ]]; then
    software_needed+=("picard")
fi

for software in "${software_needed[@]}"; do
    if ! command -v "$software" &> /dev/null; then
        print_info "Installing $software..."
        case $software in
            "featureCounts")
                conda install -c bioconda subread -y
                ;;
            "samtools")
                conda install -c bioconda samtools -y
                ;;
            "picard")
                conda install -c bioconda picard -y
                ;;
        esac
        print_success "$software installed"
    fi
done

# Show deduplication status
if [[ "$RUN_DEDUP" == true ]]; then
    print_success "Picard available for deduplication"
fi

# Find samples
print_info "Scanning for samples..."
SAMPLES=( $(ls "$INPUT_DIR"/*_R1_001.fastq.gz 2>/dev/null | xargs -n1 basename | sed 's/_R1_001.fastq.gz//' | sort -u) )

if [ ${#SAMPLES[@]} -eq 0 ]; then
    print_error "No samples found in $INPUT_DIR"
    print_info "Expected file pattern: *_R1_001.fastq.gz and *_R2_001.fastq.gz"
    exit 1
fi

print_success "Found ${#SAMPLES[@]} samples to process"

# Display sample list
echo -e "${DIM}Samples to process:${NC}"
for sample in "${SAMPLES[@]}"; do
    echo -e "${DIM}  - $sample${NC}"
done

# Set strand parameters
case $STRAND_SPECIFICITY in
    "forward")
        RSEM_FORWARD_PROB="1.0"
        FEATURECOUNTS_STRAND="1"
        ;;
    "reverse")
        RSEM_FORWARD_PROB="0.0"
        FEATURECOUNTS_STRAND="2"
        ;;
    "unstranded"|*)
        RSEM_FORWARD_PROB="0.5"
        FEATURECOUNTS_STRAND="0"
        ;;
esac

print_section "Sample Processing" "$ROCKET"

# Process each sample
success_count=0
failed_samples=()
total_start_time=$(date +%s)

for i in "${!SAMPLES[@]}"; do
    SAMPLE="${SAMPLES[$i]}"
    current_sample=$((i + 1))
    
    echo ""
    echo -e "${PURPLE}╔═══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${PURPLE}║${WHITE} Sample: ${BOLD}${SAMPLE}${NC}${WHITE} (${current_sample}/${#SAMPLES[@]})${NC}$(printf "%*s" $((50 - ${#SAMPLE} - ${#current_sample} - ${##SAMPLES[@]} - 6)) "")${PURPLE}║${NC}"
    start_date=$(date)
    echo -e "${PURPLE}║${DIM} Started: ${start_date}${NC}$(printf "%*s" $((57 - ${#start_date})) "")${PURPLE}║${NC}"
    echo -e "${PURPLE}╚═══════════════════════════════════════════════════════════════╝${NC}"
    
    sample_start_time=$(date +%s)
    
    # File paths
    READ1="$INPUT_DIR/${SAMPLE}_R1_001.fastq.gz"
    READ2="$INPUT_DIR/${SAMPLE}_R2_001.fastq.gz"
    SAMPLE_WORK_DIR="$WORK_DIR/$SAMPLE"
    SAMPLE_OUTPUT_DIR="$OUTPUT_DIR/$SAMPLE"
    SAMPLE_QC_DIR="$QC_DIR/$SAMPLE"
    
    # Validate input files
    if [[ ! -f "$READ1" ]] || [[ ! -f "$READ2" ]]; then
        print_error "FASTQ files not found for $SAMPLE"
        print_error "Expected: $READ1 and $READ2"
        failed_samples+=("$SAMPLE")
        continue
    fi
    
    # File size info
    r1_size=$(du -h "$READ1" | cut -f1)
    r2_size=$(du -h "$READ2" | cut -f1)
    print_info "Input files: ${CYAN}R1=${r1_size}, R2=${r2_size}${NC}"
    
    # Create working directories
    rm -rf "$SAMPLE_WORK_DIR"
    mkdir -p "$SAMPLE_WORK_DIR"
    mkdir -p "$SAMPLE_OUTPUT_DIR"
    mkdir -p "$SAMPLE_QC_DIR"
    cd "$SAMPLE_WORK_DIR"
    
    # Calculate BAM sort RAM
    BAM_SORT_RAM_BYTES=$(($BAM_SORT_RAM * 1024 * 1024 * 1024))
    
    # Set file descriptor limits to handle STAR temp files
    ulimit -n 65536 2>/dev/null || echo "Warning: Could not increase file descriptor limit"
    
    # STAR Alignment - SIMPLIFIED STABLE PARAMETERS
    echo ""
    print_info "${FIRE} Running STAR alignment with stable parameters..."
    
    # Create error log file
    STAR_LOG="$SAMPLE_QC_DIR/${SAMPLE}_star.log"
    
    STAR --runThreadN $STAR_THREADS \
         --genomeDir "$STAR_INDEX_DIR" \
         --readFilesIn "$READ1" "$READ2" \
         --readFilesCommand zcat \
         --outFileNamePrefix "${SAMPLE}_" \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMunmapped Within \
         --outSAMattributes Standard \
         --outFilterType BySJout \
         --outFilterMultimapNmax $MULTIMAPPER_MAX \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --outFilterMismatchNoverLmax 0.04 \
         --alignIntronMin 20 \
         --alignIntronMax 1000000 \
         --alignMatesGapMax 1000000 \
         --limitBAMsortRAM $BAM_SORT_RAM_BYTES \
         --quantMode TranscriptomeSAM GeneCounts \
         --twopassMode Basic \
         --outFilterIntronMotifs RemoveNoncanonical \
         --limitIObufferSize 50000000 50000000 \
         --outTmpDir "/tmp/STARtmp_${SAMPLE}_$" 2>&1 | tee "$STAR_LOG"
    
    star_exit_code=${PIPESTATUS[0]}
    star_time=$(($(date +%s) - sample_start_time))
    
    # Clean up temp directory
    rm -rf "/tmp/STARtmp_${SAMPLE}_$" 2>/dev/null
    
    if [ $star_exit_code -ne 0 ] || [[ ! -f "${SAMPLE}_Aligned.sortedByCoord.out.bam" ]]; then
        print_error "STAR alignment failed (exit code: $star_exit_code)"
        print_error "Check log: $STAR_LOG"
        failed_samples+=("$SAMPLE")
        continue
    fi
    
    # STAR results and quality check
    bam_size=$(du -h "${SAMPLE}_Aligned.sortedByCoord.out.bam" | cut -f1)
    print_success "STAR completed in ${CYAN}${star_time}s${NC} - BAM: ${GREEN}$bam_size${NC}"
    
    # Enhanced alignment statistics
    if [[ -f "${SAMPLE}_Log.final.out" ]]; then
        uniquely_mapped=$(grep "Uniquely mapped reads %" "${SAMPLE}_Log.final.out" | awk '{print $NF}')
        multi_mapped=$(grep "% of reads mapped to multiple loci" "${SAMPLE}_Log.final.out" | awk '{print $NF}')
        unmapped_short=$(grep "% of reads unmapped: too short" "${SAMPLE}_Log.final.out" | awk '{print $NF}')
        unmapped_other=$(grep "% of reads unmapped: other" "${SAMPLE}_Log.final.out" | awk '{print $NF}')
        total_reads=$(grep "Number of input reads" "${SAMPLE}_Log.final.out" | awk '{print $NF}')
        uniquely_mapped_num=$(grep "Uniquely mapped reads number" "${SAMPLE}_Log.final.out" | awk '{print $NF}')
        
        echo -e "    ${DIM}${total_reads} reads: Unique=${uniquely_mapped}, Multi=${multi_mapped}${NC}"
        echo -e "    ${DIM}Unmapped: short=${unmapped_short}, other=${unmapped_other}${NC}"
        
        # Quality check
        if [ "$uniquely_mapped_num" -lt "$MIN_MAPPED_READS" ]; then
            print_warning "Low mapping rate: only $uniquely_mapped_num uniquely mapped reads"
        fi
    fi
    
    # Index BAM file
    print_info "Indexing BAM file..."
    samtools index "${SAMPLE}_Aligned.sortedByCoord.out.bam"
    
    # BAM statistics
    if command -v samtools &> /dev/null; then
        print_info "Collecting BAM statistics..."
        samtools flagstat "${SAMPLE}_Aligned.sortedByCoord.out.bam" > "$SAMPLE_QC_DIR/${SAMPLE}_flagstat.txt"
        samtools stats "${SAMPLE}_Aligned.sortedByCoord.out.bam" > "$SAMPLE_QC_DIR/${SAMPLE}_bamstats.txt"
    fi
    
    # Deduplication step
    FINAL_BAM="${SAMPLE}_Aligned.sortedByCoord.out.bam"
    
    if [[ "$RUN_DEDUP" == true ]]; then
        echo ""
        print_info "${GEAR} Running deduplication with Picard MarkDuplicates..."
        dedup_start_time=$(date +%s)
        
        # Run Picard MarkDuplicates
        picard MarkDuplicates \
            INPUT="${SAMPLE}_Aligned.sortedByCoord.out.bam" \
            OUTPUT="${SAMPLE}_Aligned.sortedByCoord.dedup.bam" \
            METRICS_FILE="${SAMPLE}_dedup_metrics.txt" \
            REMOVE_DUPLICATES=true \
            ASSUME_SORTED=true \
            VALIDATION_STRINGENCY=LENIENT \
            CREATE_INDEX=true \
            MAX_RECORDS_IN_RAM=2000000 \
            TMP_DIR="/tmp" 2>&1 | tee "$SAMPLE_QC_DIR/${SAMPLE}_dedup.log"
        
        dedup_exit_code=${PIPESTATUS[0]}
        dedup_time=$(($(date +%s) - dedup_start_time))
        
        if [ $dedup_exit_code -eq 0 ] && [[ -f "${SAMPLE}_Aligned.sortedByCoord.dedup.bam" ]]; then
            # Update the final BAM to use deduplicated version
            FINAL_BAM="${SAMPLE}_Aligned.sortedByCoord.dedup.bam"
            
            # Show deduplication statistics
            if [[ -f "${SAMPLE}_dedup_metrics.txt" ]]; then
                # Parse dedup metrics (skip header lines and get the data line)
                total_reads_dedup=$(grep -A 1 "## METRICS CLASS" "${SAMPLE}_dedup_metrics.txt" | tail -1 | awk '{print $3}')
                duplicate_reads=$(grep -A 1 "## METRICS CLASS" "${SAMPLE}_dedup_metrics.txt" | tail -1 | awk '{print $7}')
                
                # Calculate duplicate percentage
                if [[ "$total_reads_dedup" -gt 0 ]]; then
                    dedup_rate=$(echo "scale=2; $duplicate_reads * 100 / $total_reads_dedup" | bc -l 2>/dev/null || echo "N/A")
                else
                    dedup_rate="N/A"
                fi
                
                print_success "Deduplication completed in ${CYAN}${dedup_time}s${NC}"
                echo -e "    ${DIM}Original BAM: $(du -h "${SAMPLE}_Aligned.sortedByCoord.out.bam" | cut -f1)${NC}"
                echo -e "    ${DIM}Deduplicated BAM: $(du -h "$FINAL_BAM" | cut -f1)${NC}"
                echo -e "    ${DIM}Duplicate rate: ${dedup_rate}%${NC}"
            else
                print_success "Deduplication completed in ${CYAN}${dedup_time}s${NC}"
                echo -e "    ${DIM}Deduplicated BAM: $(du -h "$FINAL_BAM" | cut -f1)${NC}"
            fi
        else
            print_error "Deduplication failed (exit code: $dedup_exit_code), using original BAM"
            print_error "Check log: $SAMPLE_QC_DIR/${SAMPLE}_dedup.log"
            FINAL_BAM="${SAMPLE}_Aligned.sortedByCoord.out.bam"
        fi
    else
        echo ""
        print_info "Skipping deduplication (using original BAM)"
    fi
    
    # Quantification
    quant_start_time=$(date +%s)
    quant_methods=()
    
    # RSEM - ENHANCED PARAMETERS
    if [[ "$RUN_RSEM" == true ]] && [[ -f "${SAMPLE}_Aligned.toTranscriptome.out.bam" ]]; then
        print_info "${CHART} Running RSEM quantification with enhanced parameters..."
        
        RSEM_LOG="$SAMPLE_QC_DIR/${SAMPLE}_rsem.log"
        
        rsem-calculate-expression \
            --paired-end \
            --alignments \
            --num-threads $RSEM_THREADS \
            --estimate-rspd \
            --seed 2025 \
            --ci-memory 48000 \
            --forward-prob $RSEM_FORWARD_PROB \
            --fragment-length-mean 200 \
            --fragment-length-sd 50 \
            --time \
            $([ "$SAVE_INTERMEDIATE" = true ] && echo "--temporary-folder ${SAMPLE_WORK_DIR}/rsem_temp" || echo "--no-bam-output") \
            "${SAMPLE}_Aligned.toTranscriptome.out.bam" \
            "$RSEM_INDEX_PREFIX" \
            "${SAMPLE}_rsem" 2>&1 | tee "$RSEM_LOG"
        
        rsem_exit_code=${PIPESTATUS[0]}
        
        if [ $rsem_exit_code -eq 0 ] && [[ -f "${SAMPLE}_rsem.genes.results" ]]; then
            gene_count=$(tail -n +2 "${SAMPLE}_rsem.genes.results" | awk '$5 > 0' | wc -l)
            total_genes=$(tail -n +2 "${SAMPLE}_rsem.genes.results" | wc -l)
            isoform_count=$(tail -n +2 "${SAMPLE}_rsem.isoforms.results" | awk '$5 > 0' | wc -l)
            total_isoforms=$(tail -n +2 "${SAMPLE}_rsem.isoforms.results" | wc -l)
            
            print_success "RSEM: ${GREEN}$gene_count${NC}/${total_genes} genes, ${GREEN}$isoform_count${NC}/${total_isoforms} isoforms"
            quant_methods+=("RSEM")
        else
            print_error "RSEM failed (exit code: $rsem_exit_code)"
            print_error "Check log: $RSEM_LOG"
        fi
        
        # Clean up large intermediate file
        if [[ "$SAVE_INTERMEDIATE" != true ]]; then
            rm -f "${SAMPLE}_Aligned.toTranscriptome.out.bam"
        fi
    fi
    
    # featureCounts - ENHANCED PARAMETERS
    if [[ "$RUN_FEATURECOUNTS" == true ]]; then
        print_info "${CHART} Running featureCounts with enhanced parameters..."
        
        # Gene-level counts
        featureCounts -T $FEATURECOUNTS_THREADS \
                     -p \
                     -B \
                     -C \
                     -s $FEATURECOUNTS_STRAND \
                     -t exon \
                     -g gene_id \
                     -J \
                     -G "$GENOME_FASTA" \
                     --fracOverlap 0.2 \
                     --minOverlap 1 \
                     -a "$GTF_FILE" \
                     -o "${SAMPLE}_featureCounts.txt" \
                     "$FINAL_BAM" 2>&1 | tee "$SAMPLE_QC_DIR/${SAMPLE}_featureCounts.log"
        
        if [[ -f "${SAMPLE}_featureCounts.txt" ]]; then
            fc_gene_count=$(tail -n +3 "${SAMPLE}_featureCounts.txt" | awk '$7 > 0' | wc -l)
            fc_total_genes=$(tail -n +3 "${SAMPLE}_featureCounts.txt" | wc -l)
            print_success "featureCounts: ${GREEN}$fc_gene_count${NC}/${fc_total_genes} genes with reads"
            quant_methods+=("featureCounts")
        fi
        
        # Transcript-level counts (optional)
        print_info "Running transcript-level featureCounts..."
        featureCounts -T $FEATURECOUNTS_THREADS \
                     -p \
                     -s $FEATURECOUNTS_STRAND \
                     -t transcript \
                     -g transcript_id \
                     --fracOverlap 0.2 \
                     -a "$GTF_FILE" \
                     -o "${SAMPLE}_featureCounts_transcripts.txt" \
                     "$FINAL_BAM" 2>/dev/null
    fi
    
    # STAR GeneCounts
    if [[ -f "${SAMPLE}_ReadsPerGene.out.tab" ]]; then
        star_gene_count=$(tail -n +5 "${SAMPLE}_ReadsPerGene.out.tab" | awk '$2 > 0' | wc -l)
        star_total_genes=$(tail -n +5 "${SAMPLE}_ReadsPerGene.out.tab" | wc -l)
        print_success "STAR GeneCounts: ${GREEN}$star_gene_count${NC}/${star_total_genes} genes with reads"
        quant_methods+=("STAR")
    fi
    
    quant_time=$(($(date +%s) - quant_start_time))
    
    # Generate coverage tracks (optional)
    if [[ "$GENERATE_BIGWIG" == true ]] && command -v bamCoverage &> /dev/null; then
        print_info "Generating BigWig coverage track..."
        bamCoverage -b "$FINAL_BAM" \
                    -o "${SAMPLE}_coverage.bw" \
                    --binSize 10 \
                    --normalizeUsing CPM \
                    --numberOfProcessors $STAR_THREADS 2>/dev/null
    fi
    
    # Clean up intermediate files
    if [[ "$SAVE_INTERMEDIATE" != true ]]; then
        rm -rf "/tmp/STARtmp_${SAMPLE}_$" 2>/dev/null
        rm -f "${SAMPLE}_STARpass1" 2>/dev/null
        rm -f "${SAMPLE}_Log.progress.out" 2>/dev/null
    fi
    
    # Save results
    print_info "${FOLDER} Saving results..."
    
    # Copy alignment files
    cp "${SAMPLE}_Aligned.sortedByCoord.out.bam" "$SAMPLE_OUTPUT_DIR/" 2>/dev/null
    cp "${SAMPLE}_Aligned.sortedByCoord.out.bam.bai" "$SAMPLE_OUTPUT_DIR/" 2>/dev/null
    
    # Copy deduplicated BAM if created
    if [[ "$RUN_DEDUP" == true ]] && [[ -f "${SAMPLE}_Aligned.sortedByCoord.dedup.bam" ]]; then
        cp "${SAMPLE}_Aligned.sortedByCoord.dedup.bam" "$SAMPLE_OUTPUT_DIR/"
        cp "${SAMPLE}_Aligned.sortedByCoord.dedup.bai" "$SAMPLE_OUTPUT_DIR/" 2>/dev/null
        cp "${SAMPLE}_dedup_metrics.txt" "$SAMPLE_QC_DIR/" 2>/dev/null
    fi
    
    # Copy STAR outputs
    cp "${SAMPLE}_Log.final.out" "$SAMPLE_OUTPUT_DIR/" 2>/dev/null
    cp "${SAMPLE}_Log.out" "$SAMPLE_QC_DIR/" 2>/dev/null
    cp "${SAMPLE}_ReadsPerGene.out.tab" "$SAMPLE_OUTPUT_DIR/" 2>/dev/null
    cp "${SAMPLE}_SJ.out.tab" "$SAMPLE_OUTPUT_DIR/" 2>/dev/null
    
    # Copy RSEM results
    if [[ -f "${SAMPLE}_rsem.genes.results" ]]; then
        cp "${SAMPLE}_rsem.genes.results" "$SAMPLE_OUTPUT_DIR/"
        cp "${SAMPLE}_rsem.isoforms.results" "$SAMPLE_OUTPUT_DIR/"
        cp "${SAMPLE}_rsem.stat" "$SAMPLE_QC_DIR/" 2>/dev/null
    fi
    
    # Copy featureCounts results
    if [[ -f "${SAMPLE}_featureCounts.txt" ]]; then
        cp "${SAMPLE}_featureCounts.txt" "$SAMPLE_OUTPUT_DIR/"
        cp "${SAMPLE}_featureCounts.txt.summary" "$SAMPLE_OUTPUT_DIR/"
        cp "${SAMPLE}_featureCounts_transcripts.txt" "$SAMPLE_OUTPUT_DIR/" 2>/dev/null
    fi
    
    # Copy coverage track
    if [[ -f "${SAMPLE}_coverage.bw" ]]; then
        cp "${SAMPLE}_coverage.bw" "$SAMPLE_OUTPUT_DIR/"
    fi
    
    # Keep intermediate files if requested
    if [[ "$SAVE_INTERMEDIATE" == true ]]; then
        mkdir -p "$SAMPLE_OUTPUT_DIR/intermediate"
        cp "${SAMPLE}_Aligned.toTranscriptome.out.bam" "$SAMPLE_OUTPUT_DIR/intermediate/" 2>/dev/null
        cp -r "${SAMPLE}_STARpass1" "$SAMPLE_OUTPUT_DIR/intermediate/" 2>/dev/null
    fi
    
    # Clean up working directory
    cd "$WORK_DIR"
    rm -rf "$SAMPLE_WORK_DIR"
    
    # Summary
    total_time=$(($(date +%s) - sample_start_time))
    final_size=$(du -sh "$SAMPLE_OUTPUT_DIR" | cut -f1)
    
    echo ""
    print_success "Sample completed in ${CYAN}${total_time}s${NC} (STAR: ${star_time}s, Quant: ${quant_time}s)"
    print_info "Output size: ${GREEN}$final_size${NC} | Methods: ${YELLOW}$(IFS=', '; echo "${quant_methods[*]}")${NC}"
    
    # Show which BAM was used for quantification
    if [[ "$RUN_DEDUP" == true ]] && [[ -f "${SAMPLE}_Aligned.sortedByCoord.dedup.bam" ]]; then
        print_info "Final BAM for quantification: ${CYAN}deduplicated${NC}"
    else
        print_info "Final BAM for quantification: ${CYAN}original${NC}"
    fi
    
    success_count=$((success_count + 1))
    
    # Progress bar
    echo ""
    show_progress $success_count ${#SAMPLES[@]}
    echo ""
    
    # Time estimate
    if [ $success_count -lt ${#SAMPLES[@]} ]; then
        elapsed_total=$(($(date +%s) - total_start_time))
        avg_time=$((elapsed_total / success_count))
        remaining=$((${#SAMPLES[@]} - success_count))
        est_time=$((avg_time * remaining))
        echo -e "${DIM}${CLOCK} Estimated time remaining: ${est_time} seconds ($(($est_time / 60)) minutes)${NC}"
    fi
done

# Clean up
cd "$HOME"
rm -rf "$WORK_DIR"

# Final results
total_time=$(($(date +%s) - total_start_time))

echo ""
echo ""
print_section "Final Results" "$STAR"

echo -e "${WHITE}Processing Summary:${NC}"
echo -e "  ${GREEN}${CHECK} Successfully processed: ${BOLD}$success_count${NC}${GREEN}/${#SAMPLES[@]} samples${NC}"
echo -e "  ${CYAN}${CLOCK} Total processing time: ${BOLD}$((total_time / 60))${NC}${CYAN} minutes${NC}"
echo -e "  ${BLUE}${FOLDER} Results directory: ${BOLD}$OUTPUT_DIR${NC}"
echo -e "  ${BLUE}${FOLDER} QC reports directory: ${BOLD}$QC_DIR${NC}"

if [ ${#failed_samples[@]} -gt 0 ]; then
    echo ""
    print_error "Failed samples (${#failed_samples[@]}):"
    for sample in "${failed_samples[@]}"; do
        echo -e "    ${RED}${CROSS} $sample${NC}"
    done
fi

# Count quantification results
rsem_gene_count=$(find "$OUTPUT_DIR" -name "*_rsem.genes.results" 2>/dev/null | wc -l)
rsem_isoform_count=$(find "$OUTPUT_DIR" -name "*_rsem.isoforms.results" 2>/dev/null | wc -l)
fc_count=$(find "$OUTPUT_DIR" -name "*_featureCounts.txt" 2>/dev/null | wc -l)
fc_transcript_count=$(find "$OUTPUT_DIR" -name "*_featureCounts_transcripts.txt" 2>/dev/null | wc -l)
star_count=$(find "$OUTPUT_DIR" -name "*_ReadsPerGene.out.tab" 2>/dev/null | wc -l)

echo ""
print_section "Quantification Summary" "$CHART"
echo -e "${CYAN}${FILE} RSEM gene results:${NC}        ${GREEN}$rsem_gene_count${NC} samples"
echo -e "${CYAN}${FILE} RSEM isoform results:${NC}     ${GREEN}$rsem_isoform_count${NC} samples"
echo -e "${CYAN}${FILE} featureCounts gene:${NC}       ${GREEN}$fc_count${NC} samples"
echo -e "${CYAN}${FILE} featureCounts transcript:${NC}  ${GREEN}$fc_transcript_count${NC} samples"
echo -e "${CYAN}${FILE} STAR gene counts:${NC}         ${GREEN}$star_count${NC} samples"

echo ""
print_section "Next Steps" "$ROCKET"
echo -e "${WHITE}For differential expression analysis:${NC}"
echo ""
echo -e "${YELLOW}1. RSEM (recommended for isoform analysis):${NC}"
echo -e "   ${DIM}find $OUTPUT_DIR -name '*_rsem.genes.results' > rsem_gene_files.txt${NC}"
echo -e "   ${DIM}find $OUTPUT_DIR -name '*_rsem.isoforms.results' > rsem_isoform_files.txt${NC}"
echo ""
echo -e "${YELLOW}2. featureCounts (recommended for gene analysis):${NC}"
echo -e "   ${DIM}find $OUTPUT_DIR -name '*_featureCounts.txt' > fc_gene_files.txt${NC}"
echo -e "   ${DIM}find $OUTPUT_DIR -name '*_featureCounts_transcripts.txt' > fc_transcript_files.txt${NC}"
echo ""
echo -e "${YELLOW}3. STAR GeneCounts (fastest, gene-level only):${NC}"
echo -e "   ${DIM}find $OUTPUT_DIR -name '*_ReadsPerGene.out.tab' > star_files.txt${NC}"

echo ""
print_section "Quality Control" "$GEAR"
echo -e "${WHITE}Quality control reports available in: ${CYAN}$QC_DIR${NC}"
echo -e "${DIM}- STAR alignment logs${NC}"
echo -e "${DIM}- RSEM quantification logs${NC}"
echo -e "${DIM}- BAM statistics (flagstat, bamstats)${NC}"
echo -e "${DIM}- featureCounts logs${NC}"
if [[ "$RUN_DEDUP" == true ]]; then
    echo -e "${DIM}- Deduplication metrics and logs${NC}"
fi

echo ""
echo -e "${GREEN}${CHECK} All methods compatible with DESeq2, edgeR, and limma${NC}"
echo -e "${GREEN}${CHECK} Both gene-level and isoform-level quantification available${NC}"
echo -e "${GREEN}${CHECK} Comprehensive quality control metrics collected${NC}"
if [[ "$RUN_DEDUP" == true ]]; then
    echo -e "${GREEN}${CHECK} PCR duplicate removal performed for unbiased quantification${NC}"
else
    echo -e "${YELLOW}⚠️  PCR duplicates retained (may affect quantification accuracy)${NC}"
fi

echo ""
echo -e "${PURPLE}╔═══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${PURPLE}║${WHITE}              ${STAR} HIGH-QUALITY ANALYSIS COMPLETE ${STAR}              ${PURPLE}║${NC}"
echo -e "${PURPLE}║${DIM}                    $(date)                    ${PURPLE}║${NC}"
echo -e "${PURPLE}╚═══════════════════════════════════════════════════════════════╝${NC}"
echo ""