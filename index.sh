#!/bin/bash

# 🎨 Beautiful Color Palette
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
WHITE='\033[1;37m'
BOLD='\033[1m'
DIM='\033[2m'
BLINK='\033[5m'
UNDERLINE='\033[4m'
BG_BLUE='\033[44m'
BG_GREEN='\033[42m'
BG_RED='\033[41m'
BG_YELLOW='\033[43m'
NC='\033[0m' # No Color

# 🎭 Beautiful Unicode Symbols
DNA="🧬"
STAR="⭐"
ROCKET="🚀"
GEAR="⚙️"
CHECK="✅"
CROSS="❌"
WARNING="⚠️"
FIRE="🔥"
DIAMOND="💎"
CROWN="👑"
MAGIC="✨"
HOURGLASS="⏳"
CLOCK="🕐"
FOLDER="📁"
FILE="📄"
MAGNIFY="🔍"
TOOLS="🛠️"
CHART="📊"
PACKAGE="📦"
SHIELD="🛡️"
LIGHTNING="⚡"
TELESCOPE="🔭"

# 🎨 Beautiful Progress Bar Function
show_progress() {
    local current=$1
    local total=$2
    local prefix=$3
    local emoji=$4
    local width=60
    local percentage=$((current * 100 / total))
    local completed=$((current * width / total))
    local remaining=$((width - completed))
    
    # Color gradient based on progress
    local color
    if [ $percentage -lt 25 ]; then
        color=$RED
    elif [ $percentage -lt 50 ]; then
        color=$YELLOW
    elif [ $percentage -lt 75 ]; then
        color=$BLUE
    else
        color=$GREEN
    fi
    
    printf "\r${emoji} ${WHITE}${prefix}${NC} ${color}["
    printf "%${completed}s" | tr ' ' '█'
    printf "${DIM}%${remaining}s${NC}" | tr ' ' '░'
    printf "${color}] ${WHITE}%3d%%${NC} ${DIM}(%d/%d)${NC}" $percentage $current $total
}

# 🎪 Animated Loading Function
show_loading() {
    local message=$1
    local duration=$2
    local emoji=$3
    local chars="⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"
    local delay=0.1
    local temp_file="/tmp/loading_$$"
    
    # Start background process
    (sleep $duration; echo "done" > $temp_file) &
    local bg_pid=$!
    
    local i=0
    while [[ ! -f $temp_file ]]; do
        local char=${chars:$((i % ${#chars})):1}
        printf "\r${emoji} ${CYAN}${char}${NC} ${WHITE}${message}${NC}${BLINK}...${NC}"
        sleep $delay
        ((i++))
    done
    
    # Clean up
    kill $bg_pid 2>/dev/null
    rm -f $temp_file
    printf "\r${CHECK} ${GREEN}${message} completed!${NC}\n"
}

# 🎨 Fancy Header with ASCII Art
print_banner() {
    clear
    echo -e "${PURPLE}╔══════════════════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${PURPLE}║${WHITE}                                                                              ${PURPLE}║${NC}"
    echo -e "${PURPLE}║${CYAN}    ${DNA}${MAGIC}  ${BOLD}STAR & RSEM Genome Indexing Pipeline${NC}  ${MAGIC}${DNA}${CYAN}    ${PURPLE}║${NC}"
    echo -e "${PURPLE}║${WHITE}                                                                              ${PURPLE}║${NC}"
    echo -e "${PURPLE}║${YELLOW}           ${CROWN} Production-Ready RNA-seq Index Builder ${CROWN}           ${PURPLE}║${NC}"
    echo -e "${PURPLE}║${WHITE}                                                                              ${PURPLE}║${NC}"
    echo -e "${PURPLE}║${DIM}              Optimized for High-Performance Computing               ${PURPLE}║${NC}"
    echo -e "${PURPLE}║${DIM}                   Enhanced with Beautiful Progress UI                ${PURPLE}║${NC}"
    echo -e "${PURPLE}║${WHITE}                                                                              ${PURPLE}║${NC}"
    echo -e "${PURPLE}╚══════════════════════════════════════════════════════════════════════════════╝${NC}"
    echo ""
}

# 🎯 Beautiful Section Headers
print_section() {
    local title=$1
    local emoji=$2
    local color=$3
    
    echo ""
    echo -e "${color}${emoji} ═══════════════════════════════════════════════════════${NC}"
    echo -e "${color}${emoji} ${BOLD}${title}${NC}"
    echo -e "${color}${emoji} ═══════════════════════════════════════════════════════${NC}"
    echo ""
}

# 🎨 Status Messages with Beautiful Formatting
print_status() {
    local type=$1
    local message=$2
    local emoji=$3
    
    case $type in
        "info")
            echo -e "${BLUE}${emoji} ${WHITE}${message}${NC}"
            ;;
        "success")
            echo -e "${GREEN}${emoji} ${BOLD}${message}${NC}"
            ;;
        "warning")
            echo -e "${YELLOW}${emoji} ${BOLD}WARNING:${NC} ${YELLOW}${message}${NC}"
            ;;
        "error")
            echo -e "${RED}${emoji} ${BOLD}ERROR:${NC} ${RED}${message}${NC}"
            ;;
        "highlight")
            echo -e "${PURPLE}${emoji} ${BOLD}${message}${NC}"
            ;;
    esac
}

# 🎭 Beautiful Configuration Display
show_config() {
    echo -e "${BG_BLUE}${WHITE} 📋 CONFIGURATION OVERVIEW ${NC}"
    echo ""
    echo -e "${CYAN}${FOLDER} Genome Directory:${NC}     ${WHITE}$GENOME_DIR${NC}"
    echo -e "${CYAN}${FILE} Genome FASTA:${NC}         ${WHITE}$GENOME_FASTA_ORIGINAL${NC}"
    echo -e "${CYAN}${FILE} GTF Annotation:${NC}       ${WHITE}$GTF_FILE${NC}"
    echo -e "${CYAN}${STAR} STAR Index:${NC}           ${WHITE}$STAR_INDEX_DIR${NC}"
    echo -e "${CYAN}${DNA} RSEM Index:${NC}            ${WHITE}$RSEM_INDEX_PREFIX${NC}"
    echo -e "${CYAN}${LIGHTNING} Threads:${NC}              ${GREEN}$THREADS${NC}"
    echo -e "${CYAN}${CHART} Max RAM:${NC}              ${GREEN}$(($MAX_RAM / 1000000000))GB${NC}"
    echo ""
}

# 🎪 File Size with Beautiful Formatting
format_size() {
    local file=$1
    if [[ -f "$file" ]]; then
        local size=$(du -h "$file" | cut -f1)
        echo -e "${GREEN}${size}${NC}"
    else
        echo -e "${RED}N/A${NC}"
    fi
}

# 🎨 Time Formatter
format_time() {
    local seconds=$1
    local hours=$((seconds / 3600))
    local minutes=$(((seconds % 3600) / 60))
    local secs=$((seconds % 60))
    
    if [ $hours -gt 0 ]; then
        echo -e "${GREEN}${hours}h ${minutes}m ${secs}s${NC}"
    elif [ $minutes -gt 0 ]; then
        echo -e "${GREEN}${minutes}m ${secs}s${NC}"
    else
        echo -e "${GREEN}${secs}s${NC}"
    fi
}

# Configuration - MODIFY THESE PATHS AS NEEDED
GENOME_DIR="/mnt/d/genome"
GENOME_FASTA_ORIGINAL="$GENOME_DIR/GRCh38.p14.genome.fa.gz"
GTF_FILE="$GENOME_DIR/gencode.v47.chr_patch_hapl_scaff.annotation.gtf"

# Output directories
STAR_INDEX_DIR="$GENOME_DIR/STAR_index_human"
RSEM_INDEX_DIR="$GENOME_DIR/RSEM_index_human"
RSEM_INDEX_PREFIX="$RSEM_INDEX_DIR/human_rsem"

# WSL-optimized settings
THREADS=24
MAX_RAM=65000000000  # 65GB for WSL stability

# 🎨 Beautiful Banner
print_banner

# 🎯 Configuration Display
show_config

# WSL environment setup
export TMPDIR=/tmp
ulimit -n 65536
ulimit -v unlimited

# 🔍 INPUT VALIDATION
print_section "Input Validation & Preparation" "$MAGNIFY" "$BLUE"

print_status "info" "Scanning for input files..." "$TELESCOPE"

# Validate input files with beautiful progress
files_to_check=("$GENOME_FASTA_ORIGINAL" "$GTF_FILE")
total_files=${#files_to_check[@]}

for i in "${!files_to_check[@]}"; do
    file="${files_to_check[$i]}"
    show_progress $((i+1)) $total_files "Validating files" "$SHIELD"
    sleep 0.5
    
    if [[ ! -f "$file" ]]; then
        echo ""
        print_status "error" "File not found: $file" "$CROSS"
        exit 1
    fi
done

echo ""
print_status "success" "All input files validated!" "$CHECK"

# 🎭 Handle compressed genome file with beautiful UI
if [[ "$GENOME_FASTA_ORIGINAL" == *.gz ]]; then
    print_status "highlight" "Detected compressed genome file" "$PACKAGE"
    GENOME_FASTA="${GENOME_FASTA_ORIGINAL%.gz}"
    
    if [[ -f "$GENOME_FASTA" ]]; then
        print_status "success" "Uncompressed version already exists" "$CHECK"
        
        if [[ "$GENOME_FASTA_ORIGINAL" -nt "$GENOME_FASTA" ]]; then
            print_status "warning" "Compressed file is newer than uncompressed version" "$WARNING"
            echo -e "${YELLOW}${HOURGLASS} Would you like to re-extract? ${WHITE}(y/N):${NC} \c"
            read -n 1 -r
            echo
            if [[ $REPLY =~ ^[Yy]$ ]]; then
                print_status "info" "Re-extracting compressed genome file..." "$HOURGLASS"
                
                # Animated extraction with progress
                echo -e "${CYAN}${HOURGLASS} Extracting genome file...${NC}"
                start_time=$(date +%s)
                
                # Show animated progress during extraction
                (gunzip -c "$GENOME_FASTA_ORIGINAL" > "$GENOME_FASTA") &
                extract_pid=$!
                
                # Beautiful loading animation
                chars="⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"
                i=0
                while kill -0 $extract_pid 2>/dev/null; do
                    char=${chars:$((i % ${#chars})):1}
                    printf "\r${LIGHTNING} ${CYAN}${char}${NC} ${WHITE}Extracting genome sequences${NC}${BLINK}...${NC}"
                    sleep 0.1
                    ((i++))
                done
                
                wait $extract_pid
                extract_time=$(($(date +%s) - start_time))
                
                printf "\r${CHECK} ${GREEN}Genome extraction completed in $(format_time $extract_time)!${NC}\n"
                
                if [[ -f "$GENOME_FASTA" ]]; then
                    extracted_size=$(format_size "$GENOME_FASTA")
                    print_status "success" "Extracted size: $extracted_size" "$DIAMOND"
                else
                    print_status "error" "Failed to extract genome file" "$CROSS"
                    exit 1
                fi
            fi
        fi
    else
        print_status "info" "Extracting compressed genome file..." "$HOURGLASS"
        echo -e "${DIM}   Source: $GENOME_FASTA_ORIGINAL${NC}"
        echo -e "${DIM}   Target: $GENOME_FASTA${NC}"
        
        start_time=$(date +%s)
        
        # Beautiful extraction with progress
        (gunzip -c "$GENOME_FASTA_ORIGINAL" > "$GENOME_FASTA") &
        extract_pid=$!
        
        # Loading animation
        chars="⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"
        i=0
        while kill -0 $extract_pid 2>/dev/null; do
            char=${chars:$((i % ${#chars})):1}
            printf "\r${LIGHTNING} ${CYAN}${char}${NC} ${WHITE}Extracting genome sequences${NC}${BLINK}...${NC}"
            sleep 0.1
            ((i++))
        done
        
        wait $extract_pid
        extract_time=$(($(date +%s) - start_time))
        
        printf "\r${CHECK} ${GREEN}Genome extraction completed in $(format_time $extract_time)!${NC}\n"
        
        if [[ -f "$GENOME_FASTA" ]]; then
            extracted_size=$(format_size "$GENOME_FASTA")
            print_status "success" "Extracted size: $extracted_size" "$DIAMOND"
        else
            print_status "error" "Failed to extract genome file" "$CROSS"
            exit 1
        fi
    fi
else
    GENOME_FASTA="$GENOME_FASTA_ORIGINAL"
    print_status "success" "Genome file is already uncompressed" "$CHECK"
fi

print_status "info" "Using genome file: $GENOME_FASTA" "$FILE"

# 🧬 GTF Validation with Progress
print_status "info" "Validating GTF format..." "$DNA"

# Progress bar for GTF validation
for i in {1..20}; do
    show_progress $i 20 "Checking GTF format" "$MAGNIFY"
    sleep 0.05
done

echo ""

gtf_lines=$(head -1000 "$GTF_FILE" | grep -v "^#" | wc -l)
if [ $gtf_lines -eq 0 ]; then
    print_status "error" "GTF file appears to be empty or only contains headers" "$CROSS"
    exit 1
fi

if ! head -1000 "$GTF_FILE" | grep -q 'gene_id\|transcript_id'; then
    print_status "error" "GTF file missing required gene_id or transcript_id attributes" "$CROSS"
    exit 1
fi

print_status "success" "GTF format validated successfully!" "$CHECK"

# 📊 File Information Display
echo ""
echo -e "${BG_GREEN}${WHITE} 📊 FILE INFORMATION ${NC}"
echo ""

fasta_size=$(format_size "$GENOME_FASTA")
gtf_size=$(format_size "$GTF_FILE")
available_mem=$(free -g | awk 'NR==2{printf "%.0f", $7}')

echo -e "${CYAN}${FILE} Genome FASTA:${NC}         $fasta_size"
echo -e "${CYAN}${FILE} GTF Annotation:${NC}       $gtf_size"
echo -e "${CYAN}${CHART} Available Memory:${NC}      ${GREEN}${available_mem}GB${NC}"

# 📦 SOFTWARE INSTALLATION
print_section "Software Installation & Verification" "$PACKAGE" "$YELLOW"

print_status "info" "Installing required bioinformatics tools..." "$TOOLS"

# Progress bar for software installation
echo -e "${CYAN}${PACKAGE} Installing STAR, RSEM, Bowtie, and Samtools...${NC}"
conda install -c bioconda star=2.7.10a rsem bowtie samtools -y > /dev/null 2>&1

# Verification with progress
software_list=("STAR" "rsem-prepare-reference" "bowtie" "samtools")
total_software=${#software_list[@]}

echo ""
print_status "info" "Verifying software installations..." "$SHIELD"

for i in "${!software_list[@]}"; do
    software="${software_list[$i]}"
    show_progress $((i+1)) $total_software "Verifying installations" "$SHIELD"
    sleep 0.3
    
    if ! command -v "$software" &> /dev/null; then
        echo ""
        print_status "error" "$software installation failed" "$CROSS"
        exit 1
    fi
done

echo ""

# Get versions with beautiful display
STAR_VERSION=$(STAR --version 2>/dev/null)
RSEM_VERSION=$(rsem-calculate-expression --version 2>&1 | head -1)
BOWTIE_VERSION=$(bowtie --version 2>&1 | head -1)

echo -e "${BG_GREEN}${WHITE} 🛠️  INSTALLED SOFTWARE VERSIONS ${NC}"
echo ""
echo -e "${GREEN}${CHECK} STAR:${NC}        ${WHITE}$STAR_VERSION${NC}"
echo -e "${GREEN}${CHECK} RSEM:${NC}        ${WHITE}$RSEM_VERSION${NC}"
echo -e "${GREEN}${CHECK} Bowtie:${NC}      ${WHITE}$BOWTIE_VERSION${NC}"

# 📁 Directory Preparation
print_status "info" "Preparing output directories..." "$FOLDER"
mkdir -p "$STAR_INDEX_DIR"
mkdir -p "$RSEM_INDEX_DIR"
print_status "success" "Output directories created!" "$CHECK"

# ⭐ STAR INDEX CREATION
print_section "STAR Genome Index Creation" "$STAR" "$GREEN"

if [[ -f "$STAR_INDEX_DIR/genomeParameters.txt" ]]; then
    existing_size=$(du -sh "$STAR_INDEX_DIR" | cut -f1)
    print_status "success" "STAR index already exists - Size: $existing_size" "$CHECK"
else
    print_status "highlight" "Starting STAR genome indexing..." "$STAR"
    echo -e "${DIM}   This process will take 30-60 minutes...${NC}"
    echo -e "${DIM}   Started at: $(date)${NC}"
    echo ""
    
    start_time=$(date +%s)
    
    # Beautiful progress indicator for STAR
    echo -e "${STAR} ${WHITE}Running STAR genome indexing...${NC}"
    
    # Run STAR with output redirected
    (
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
             --outTmpDir "/tmp/STARtmp_$$" > /dev/null 2>&1
    ) &
    
    star_pid=$!
    
    # Beautiful animated progress
    chars="⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"
    i=0
    start_display=$(date +%s)
    
    while kill -0 $star_pid 2>/dev/null; do
        char=${chars:$((i % ${#chars})):1}
        elapsed=$(($(date +%s) - start_display))
        printf "\r${STAR} ${CYAN}${char}${NC} ${WHITE}Building STAR index${NC} ${DIM}($(format_time $elapsed) elapsed)${NC}${BLINK}...${NC}"
        sleep 0.2
        ((i++))
    done
    
    wait $star_pid
    star_exit_code=$?
    star_time=$(($(date +%s) - start_time))
    
    # Clean up temp directory
    rm -rf "/tmp/STARtmp_$$" 2>/dev/null
    
    printf "\r"
    
    if [ $star_exit_code -ne 0 ] || [[ ! -f "$STAR_INDEX_DIR/genomeParameters.txt" ]]; then
        print_status "error" "STAR index creation failed (exit code: $star_exit_code)" "$CROSS"
        exit 1
    fi
    
    # Success with beautiful formatting
    index_size=$(du -sh "$STAR_INDEX_DIR" | cut -f1)
    print_status "success" "STAR index created successfully!" "$CHECK"
    echo -e "${GREEN}   ${CLOCK} Time taken: $(format_time $star_time)${NC}"
    echo -e "${GREEN}   ${CHART} Index size: $index_size${NC}"
fi

# 🧬 RSEM INDEX CREATION
print_section "RSEM Reference Index Creation" "$DNA" "$PURPLE"

# Check if RSEM index already exists
if [[ -f "${RSEM_INDEX_PREFIX}.grp" && -f "${RSEM_INDEX_PREFIX}.ti" && -f "${RSEM_INDEX_PREFIX}.transcripts.fa" ]]; then
    echo -e "${YELLOW}${WARNING} RSEM index already exists!${NC}"
    echo ""
    echo -e "${DIM}Existing files:${NC}"
    ls -la "${RSEM_INDEX_PREFIX}".* | head -3
    echo ""
    echo -e "${YELLOW}${HOURGLASS} Do you want to overwrite RSEM index? ${WHITE}(y/N):${NC} \c"
    read -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        print_status "info" "Keeping existing RSEM index..." "$CHECK"
    else
        print_status "info" "Removing existing RSEM index..." "$CROSS"
        rm -f "${RSEM_INDEX_PREFIX}".*
    fi
fi

# Create RSEM index if not exists or user chose to overwrite
if [[ ! -f "${RSEM_INDEX_PREFIX}.grp" ]]; then
    print_status "highlight" "Starting RSEM reference indexing..." "$DNA"
    echo -e "${DIM}   This process will take 20-45 minutes...${NC}"
    echo -e "${DIM}   Started at: $(date)${NC}"
    echo ""
    
    start_time=$(date +%s)
    
    # Run RSEM with beautiful progress
    echo -e "${DNA} ${WHITE}Running RSEM reference preparation...${NC}"
    
    (
        rsem-prepare-reference \
            --gtf "$GTF_FILE" \
            --bowtie \
            --num-threads $THREADS \
            "$GENOME_FASTA" \
            "$RSEM_INDEX_PREFIX" > /dev/null 2>&1
    ) &
    
    rsem_pid=$!
    
    # Beautiful animated progress
    chars="⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"
    i=0
    start_display=$(date +%s)
    
    while kill -0 $rsem_pid 2>/dev/null; do
        char=${chars:$((i % ${#chars})):1}
        elapsed=$(($(date +%s) - start_display))
        printf "\r${DNA} ${PURPLE}${char}${NC} ${WHITE}Building RSEM reference${NC} ${DIM}($(format_time $elapsed) elapsed)${NC}${BLINK}...${NC}"
        sleep 0.2
        ((i++))
    done
    
    wait $rsem_pid
    rsem_exit_code=$?
    rsem_time=$(($(date +%s) - start_time))
    
    printf "\r"
    
    if [ $rsem_exit_code -eq 0 ] && [[ -f "${RSEM_INDEX_PREFIX}.grp" ]]; then
        # Success with beautiful formatting
        index_size=$(du -sh "$RSEM_INDEX_DIR" | cut -f1)
        print_status "success" "RSEM index created successfully!" "$CHECK"
        echo -e "${GREEN}   ${CLOCK} Time taken: $(format_time $rsem_time)${NC}"
        echo -e "${GREEN}   ${CHART} Index size: $index_size${NC}"
    else
        print_status "warning" "RSEM indexing encountered issues (exit code: $rsem_exit_code)" "$WARNING"
        print_status "info" "Checking for essential RSEM files..." "$MAGNIFY"
        
        if [[ -f "${RSEM_INDEX_PREFIX}.grp" ]]; then
            print_status "success" "RSEM group file exists - basic functionality available" "$CHECK"
        else
            print_status "warning" "Trying fallback method without bowtie indices..." "$TOOLS"
            
            # Fallback method
            rsem-prepare-reference \
                --gtf "$GTF_FILE" \
                --num-threads $THREADS \
                "$GENOME_FASTA" \
                "$RSEM_INDEX_PREFIX" > /dev/null 2>&1
            
            if [[ -f "${RSEM_INDEX_PREFIX}.grp" ]]; then
                print_status "success" "RSEM reference files created without bowtie indices" "$CHECK"
                print_status "info" "You can still use STAR+RSEM workflow" "$DNA"
            else
                print_status "error" "RSEM reference creation failed completely" "$CROSS"
                exit 1
            fi
        fi
    fi
else
    existing_size=$(du -sh "$RSEM_INDEX_DIR" | cut -f1)
    print_status "success" "RSEM index found - Size: $existing_size" "$CHECK"
fi

# 🎊 FINAL VERIFICATION & SUMMARY
print_section "Final Verification & Summary" "$ROCKET" "$CYAN"

# Verify STAR index
print_status "info" "Verifying STAR index..." "$STAR"
if [[ -f "$STAR_INDEX_DIR/genomeParameters.txt" ]]; then
    star_size=$(du -sh "$STAR_INDEX_DIR" | cut -f1)
    print_status "success" "STAR index verified - $star_size" "$CHECK"
else
    print_status "error" "STAR index missing" "$CROSS"
    exit 1
fi

# Verify RSEM index
print_status "info" "Verifying RSEM index..." "$DNA"
if [[ -f "${RSEM_INDEX_PREFIX}.grp" ]]; then
    rsem_size=$(du -sh "$RSEM_INDEX_DIR" | cut -f1)
    print_status "success" "RSEM index verified - $rsem_size" "$CHECK"
    
    # Check available files with beautiful display
    echo ""
    echo -e "${BG_PURPLE}${WHITE} 🧬 AVAILABLE RSEM FILES ${NC}"
    echo ""
    
    [[ -f "${RSEM_INDEX_PREFIX}.grp" ]] && echo -e "${GREEN}${CHECK} Group file (.grp)${NC}"
    [[ -f "${RSEM_INDEX_PREFIX}.ti" ]] && echo -e "${GREEN}${CHECK} Transcript info (.ti)${NC}"
    [[ -f "${RSEM_INDEX_PREFIX}.transcripts.fa" ]] && echo -e "${GREEN}${CHECK} Transcripts FASTA${NC}"
    [[ -f "${RSEM_INDEX_PREFIX}.seq" ]] && echo -e "${GREEN}${CHECK} Sequence file (.seq)${NC}"
    [[ -f "${RSEM_INDEX_PREFIX}.chrlist" ]] && echo -e "${GREEN}${CHECK} Chromosome list${NC}"
    [[ -f "${RSEM_INDEX_PREFIX}.idx.fa" ]] && echo -e "${GREEN}${CHECK} Index FASTA${NC}"
    
    # Check for bowtie indices
    if [[ -f "${RSEM_INDEX_PREFIX}.1.ebwt" && -s "${RSEM_INDEX_PREFIX}.1.ebwt" ]]; then
        echo -e "${GREEN}${CHECK} Bowtie indices available (full RSEM functionality)${NC}"
        RSEM_FUNCTIONALITY="FULL"
    else
        echo -e "${YELLOW}${WARNING} Bowtie indices not available (STAR+RSEM workflow still works)${NC}"
        RSEM_FUNCTIONALITY="BASIC"
    fi
    
else
    print_status "warning" "RSEM index incomplete - will use STAR-only quantification" "$WARNING"
    RSEM_FUNCTIONALITY="NONE"
fi

# 🎉 BEAUTIFUL SUCCESS BANNER
echo ""
echo -e "${GREEN}╔══════════════════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${GREEN}║${WHITE}                                                                              ${GREEN}║${NC}"
echo -e "${GREEN}║${MAGIC}${CROWN}  ${BOLD}INDEX CREATION COMPLETED SUCCESSFULLY!${NC}  ${CROWN}${MAGIC}${GREEN}                ${GREEN}║${NC}"
echo -e "${GREEN}║${WHITE}                                                                              ${GREEN}║${NC}"
echo -e "${GREEN}╚══════════════════════════════════════════════════════════════════════════════╝${NC}"

# 📊 Final Statistics
echo ""
echo -e "${BG_CYAN}${WHITE} 📊 FINAL STATISTICS ${NC}"
echo ""

total_time=$(date)
echo -e "${CYAN}${CLOCK} Completed at:${NC}        ${WHITE}$total_time${NC}"

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

echo -e "${CYAN}${STAR} STAR Index Size:${NC}      ${GREEN}$star_size${NC}"
echo -e "${CYAN}${DNA} RSEM Index Size:${NC}       ${GREEN}$rsem_size${NC}"

# 🎯 Status Summary with Beautiful Icons
echo ""
echo -e "${BG_GREEN}${WHITE} 🎯 INDEX STATUS SUMMARY ${NC}"
echo ""
echo -e "${GREEN}${CHECK} STAR: Ready for alignment${NC}"

case $RSEM_FUNCTIONALITY in
    "FULL")
        echo -e "${GREEN}${CHECK} RSEM: Full functionality (with Bowtie indices)${NC}"
        echo -e "${GREEN}${ROCKET} Recommended workflow: STAR+RSEM${NC}"
        ;;
    "BASIC")
        echo -e "${YELLOW}${CHECK} RSEM: Basic functionality (without Bowtie indices)${NC}"
        echo -e "${YELLOW}${ROCKET} Recommended workflow: STAR+RSEM (transcriptome alignment)${NC}"
        ;;
    "NONE")
        echo -e "${RED}${WARNING} RSEM: Not available${NC}"
        echo -e "${BLUE}${ROCKET} Alternative: Use STAR GeneCounts for quantification${NC}"
        ;;
esac

# 🚀 Next Steps with Beautiful Formatting
echo ""
echo -e "${BG_YELLOW}${WHITE} 🚀 NEXT STEPS ${NC}"
echo ""
echo -e "${YELLOW}1.${NC} ${WHITE}Use your batch alignment script with these indices${NC}"
echo -e "${YELLOW}2.${NC} ${CYAN}STAR index:${NC} ${WHITE}$STAR_INDEX_DIR${NC}"
echo -e "${YELLOW}3.${NC} ${CYAN}RSEM index:${NC} ${WHITE}$RSEM_INDEX_PREFIX${NC}"

echo ""
echo -e "${DIM}Example alignment commands:${NC}"
echo -e "${DIM}STAR --genomeDir $STAR_INDEX_DIR --quantMode TranscriptomeSAM ...${NC}"
echo -e "${DIM}rsem-calculate-expression --alignments ... $RSEM_INDEX_PREFIX ...${NC}"

# 💡 Cleanup Tip
echo ""
if [[ "$GENOME_FASTA" != "$GENOME_FASTA_ORIGINAL" ]]; then
    echo -e "${BG_BLUE}${WHITE} 💡 CLEANUP TIP ${NC}"
    echo ""
    echo -e "${BLUE}${LIGHTNING} After successful indexing, you can optionally remove the large${NC}"
    echo -e "${BLUE}${LIGHTNING} uncompressed genome file to save space:${NC}"
    echo -e "${DIM}   rm '$GENOME_FASTA'${NC}"
    cleanup_size=$(format_size "$GENOME_FASTA")
    echo -e "${BLUE}${LIGHTNING} This will save $cleanup_size of disk space${NC}"
fi

# 🎊 Final Success Message
echo ""
echo -e "${MAGIC}${GREEN}✨ All indices are ready for RNA-seq analysis! ✨${NC}${MAGIC}"
echo ""
