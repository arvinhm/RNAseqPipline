#!/usr/bin/env python3
"""
🧬 RNA-seq Results Merger 🧬
Merges RSEM and STAR (featureCounts) results from multiple samples
Creates 6 beautiful CSV output files with gene names and transcript names as annotations
"""

import os
import pandas as pd
import re
from pathlib import Path
import glob
import sys
from datetime import datetime

# Color codes for beautiful terminal output
class Colors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    PURPLE = '\033[35m'
    YELLOW = '\033[33m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    BLUE = '\033[34m'

def print_banner():
    """Print a beautiful banner"""
    banner = f"""
{Colors.HEADER}{'='*80}{Colors.ENDC}
{Colors.BOLD}{Colors.CYAN}
    🧬 RNA-seq Results Merger Pro 🧬
    
    ┌─────────────────────────────────────────────────────────┐
    │  🔬 Merging RSEM & STAR Results into Beautiful CSVs 📊  │
    │  🚀 Processing Multiple Samples Simultaneously        │
    │  ✨ Creating 6 Output Files with Gene Annotations     │
    └─────────────────────────────────────────────────────────┘
{Colors.ENDC}
{Colors.HEADER}{'='*80}{Colors.ENDC}
"""
    print(banner)

def print_step(step_num, title, emoji="🔄"):
    """Print a colorful step header"""
    print(f"\n{Colors.BOLD}{Colors.BLUE}{'='*20} Step {step_num}: {emoji} {title} {'='*20}{Colors.ENDC}")

def print_success(message, emoji="✅"):
    """Print success message in green"""
    print(f"{Colors.OKGREEN}{Colors.BOLD}{emoji} {message}{Colors.ENDC}")

def print_info(message, emoji="ℹ️"):
    """Print info message in cyan"""
    print(f"{Colors.OKCYAN}{emoji} {message}{Colors.ENDC}")

def print_warning(message, emoji="⚠️"):
    """Print warning message in yellow"""
    print(f"{Colors.WARNING}{emoji} WARNING: {message}{Colors.ENDC}")

def print_error(message, emoji="❌"):
    """Print error message in red"""
    print(f"{Colors.FAIL}{emoji} ERROR: {message}{Colors.ENDC}")

def print_processing(sample, file_type, emoji="⚙️"):
    """Print processing message"""
    print(f"  {Colors.CYAN}{emoji} Processing {Colors.BOLD}{sample}{Colors.ENDC}{Colors.CYAN}: {file_type}{Colors.ENDC}")

def print_progress_bar(current, total, prefix="Progress", emoji="🔄"):
    """Print a colorful progress bar"""
    percent = int(100 * (current / float(total)))
    bar_length = 50
    filled_length = int(bar_length * current // total)
    bar = '█' * filled_length + '-' * (bar_length - filled_length)
    
    color = Colors.GREEN if percent == 100 else Colors.CYAN
    print(f'\r{color}{emoji} {prefix}: |{bar}| {percent}% ({current}/{total}){Colors.ENDC}', end='')
    if current == total:
        print()  # New line when complete

def parse_gtf_attributes(attribute_string):
    """
    🧬 Parse the attributes field of a GTF line.
    Returns a dictionary of key-value pairs.
    """
    attributes = {}
    for attr in attribute_string.split(";"):
        attr = attr.strip()
        if not attr:
            continue
        parts = attr.split(" ", 1)
        if len(parts) == 2:
            key = parts[0]
            value = parts[1].strip().strip('"')
            attributes[key] = value
    return attributes

def create_gtf_mappings(gtf_file_path):
    """
    📚 Create mappings from GTF file:
    - gene_id to gene_name
    - transcript_id to transcript_name  
    - transcript_id to gene_id
    """
    print_step(1, "Parsing GTF Annotations", "📚")
    print_info(f"Reading GTF file: {gtf_file_path}")
    
    gene_id_to_name = {}
    transcript_id_to_name = {}
    transcript_id_to_gene_id = {}
    
    total_lines = 0
    processed_lines = 0
    
    # Count total lines first
    with open(gtf_file_path, 'r') as f:
        total_lines = sum(1 for line in f if not line.startswith("#"))
    
    with open(gtf_file_path, 'r') as gtf_file:
        for line_num, line in enumerate(gtf_file, 1):
            if line.startswith("#"):
                continue
            
            processed_lines += 1
            if processed_lines % 10000 == 0:
                print_progress_bar(processed_lines, total_lines, "Parsing GTF", "📖")
            
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
                
            feature_type = fields[2]
            attribute_str = fields[8]
            attributes = parse_gtf_attributes(attribute_str)
            
            # Extract gene mappings
            if feature_type == "gene":
                gene_id = attributes.get('gene_id')
                gene_name = attributes.get('gene_name')
                if gene_id and gene_name:
                    gene_id_to_name[gene_id] = gene_name
            
            # Extract transcript mappings
            elif feature_type == "transcript":
                gene_id = attributes.get('gene_id')
                transcript_id = attributes.get('transcript_id')
                transcript_name = attributes.get('transcript_name')
                
                if transcript_id and transcript_name:
                    transcript_id_to_name[transcript_id] = transcript_name
                if transcript_id and gene_id:
                    transcript_id_to_gene_id[transcript_id] = gene_id
    
    print_progress_bar(total_lines, total_lines, "Parsing GTF", "📖")
    print_success(f"Found {len(gene_id_to_name):,} gene mappings", "🧬")
    print_success(f"Found {len(transcript_id_to_name):,} transcript mappings", "🧬")
    
    return gene_id_to_name, transcript_id_to_name, transcript_id_to_gene_id

def get_sample_folders(results_dir):
    """
    📁 Get list of sample folders (directories that contain analysis results)
    """
    print_step(2, "Discovering Sample Folders", "📁")
    
    sample_folders = []
    all_folders = [item for item in os.listdir(results_dir) if os.path.isdir(os.path.join(results_dir, item))]
    
    for i, item in enumerate(all_folders, 1):
        print_progress_bar(i, len(all_folders), "Scanning folders", "🔍")
        
        item_path = os.path.join(results_dir, item)
        # Check if this folder contains expected files
        featurecounts_files = glob.glob(os.path.join(item_path, "*_featureCounts.txt"))
        rsem_gene_files = glob.glob(os.path.join(item_path, "*_rsem.genes.results"))
        rsem_isoform_files = glob.glob(os.path.join(item_path, "*_rsem.isoforms.results"))
        
        if featurecounts_files and rsem_gene_files and rsem_isoform_files:
            sample_folders.append(item)
            print_success(f"Valid sample folder: {item}", "📂")
        else:
            print_warning(f"Incomplete data in folder: {item}")
    
    if sample_folders:
        print_success(f"Found {len(sample_folders)} valid sample folders! 🎉")
        for folder in sorted(sample_folders):
            print_info(f"  📊 {folder}")
    else:
        print_error("No valid sample folders found!")
    
    return sorted(sample_folders)

def process_featurecounts_files(results_dir, sample_folders, gene_id_to_name):
    """
    ⭐ Process featureCounts files and create merged gene count matrix
    """
    print_step(3, "Processing STAR featureCounts Files", "⭐")
    
    # Collect all unique genes across all samples first
    all_genes = set()
    sample_data = {}
    
    # First pass: collect all genes and sample data
    for i, sample in enumerate(sample_folders, 1):
        print_progress_bar(i, len(sample_folders), "Reading STAR files", "⭐")
        
        sample_dir = os.path.join(results_dir, sample)
        featurecounts_files = glob.glob(os.path.join(sample_dir, "*_featureCounts.txt"))
        
        if not featurecounts_files:
            print_warning(f"No featureCounts file found for {sample}")
            sample_data[sample] = {}
            continue
            
        featurecounts_file = featurecounts_files[0]
        print_processing(sample, f"featureCounts: {os.path.basename(featurecounts_file)}")
        
        # Read featureCounts file (skip comment lines)
        df = pd.read_csv(featurecounts_file, sep='\t', comment='#')
        
        # Extract gene counts (last column contains the counts)
        count_column = df.columns[-1]  # Last column contains counts
        
        sample_counts = {}
        for _, row in df.iterrows():
            gene_id = row['Geneid']
            count = row[count_column]
            
            # Convert gene_id to gene_name
            gene_name = gene_id_to_name.get(gene_id, gene_id)
            all_genes.add(gene_name)
            sample_counts[gene_name] = count
        
        sample_data[sample] = sample_counts
    
    print_info(f"🧬 Found {len(all_genes):,} unique genes across all samples")
    
    # Second pass: create complete matrix with zeros for missing genes
    print_info("🔧 Creating complete matrix (filling missing genes with 0)...")
    merged_data = {}
    
    for gene_name in all_genes:
        merged_data[gene_name] = {}
        for sample in sample_folders:
            # Use 0 if gene is missing from this sample
            merged_data[gene_name][sample] = sample_data.get(sample, {}).get(gene_name, 0)
    
    # Convert to DataFrame
    result_df = pd.DataFrame(merged_data).T
    result_df = result_df.reindex(columns=sample_folders, fill_value=0)
    
    # Report missing gene statistics
    missing_stats = {}
    for sample in sample_folders:
        missing_count = sum(1 for gene in all_genes if sample_data.get(sample, {}).get(gene, 0) == 0)
        missing_stats[sample] = missing_count
    
    max_missing = max(missing_stats.values())
    min_missing = min(missing_stats.values())
    avg_missing = sum(missing_stats.values()) / len(missing_stats)
    
    print_info(f"📊 Missing gene statistics:")
    print_info(f"   • Max missing per sample: {max_missing:,} genes")
    print_info(f"   • Min missing per sample: {min_missing:,} genes") 
    print_info(f"   • Average missing per sample: {avg_missing:.1f} genes")
    
    if max_missing > 0:
        print_warning(f"Some samples are missing up to {max_missing:,} genes - filled with zeros")
    
    print_success(f"Created STAR matrix: {len(result_df):,} genes × {len(sample_folders)} samples", "📊")
    return result_df

def process_rsem_gene_files(results_dir, sample_folders, gene_id_to_name, value_column, file_type):
    """
    🧬 Process RSEM gene files and create merged matrix for specified column
    """
    # Collect all unique genes across all samples first
    all_genes = set()
    sample_data = {}
    
    # First pass: collect all genes and sample data
    for i, sample in enumerate(sample_folders, 1):
        print_progress_bar(i, len(sample_folders), f"Reading {file_type}", "🧬")
        
        sample_dir = os.path.join(results_dir, sample)
        rsem_files = glob.glob(os.path.join(sample_dir, "*_rsem.genes.results"))
        
        if not rsem_files:
            print_warning(f"No RSEM genes file found for {sample}")
            sample_data[sample] = {}
            continue
            
        rsem_file = rsem_files[0]
        print_processing(sample, f"RSEM genes: {os.path.basename(rsem_file)}")
        
        df = pd.read_csv(rsem_file, sep='\t')
        
        sample_values = {}
        for _, row in df.iterrows():
            gene_id = row['gene_id']
            value = row[value_column]
            
            # Convert gene_id to gene_name
            gene_name = gene_id_to_name.get(gene_id, gene_id)
            all_genes.add(gene_name)
            sample_values[gene_name] = value
        
        sample_data[sample] = sample_values
    
    print_info(f"🧬 Found {len(all_genes):,} unique genes across all samples")
    
    # Second pass: create complete matrix with zeros for missing genes
    print_info("🔧 Creating complete matrix (filling missing genes with 0)...")
    merged_data = {}
    
    for gene_name in all_genes:
        merged_data[gene_name] = {}
        for sample in sample_folders:
            # Use 0 if gene is missing from this sample
            merged_data[gene_name][sample] = sample_data.get(sample, {}).get(gene_name, 0)
    
    # Convert to DataFrame
    result_df = pd.DataFrame(merged_data).T
    result_df = result_df.reindex(columns=sample_folders, fill_value=0)
    
    # Report missing gene statistics
    missing_stats = {}
    for sample in sample_folders:
        missing_count = sum(1 for gene in all_genes if sample_data.get(sample, {}).get(gene, 0) == 0)
        missing_stats[sample] = missing_count
    
    if missing_stats:
        max_missing = max(missing_stats.values())
        min_missing = min(missing_stats.values())
        avg_missing = sum(missing_stats.values()) / len(missing_stats)
        
        print_info(f"📊 Missing gene statistics for {file_type}:")
        print_info(f"   • Max missing per sample: {max_missing:,} genes")
        print_info(f"   • Min missing per sample: {min_missing:,} genes") 
        print_info(f"   • Average missing per sample: {avg_missing:.1f} genes")
        
        if max_missing > 0:
            print_warning(f"Some samples are missing up to {max_missing:,} genes - filled with zeros")
    
    print_success(f"Created {file_type} matrix: {len(result_df):,} genes × {len(sample_folders)} samples", "📊")
    return result_df

def process_rsem_isoform_files(results_dir, sample_folders, transcript_id_to_name, 
                              transcript_id_to_gene_id, value_column, file_type):
    """
    🧬 Process RSEM isoform files and create merged matrix for specified column
    """
    # Collect all unique transcripts across all samples first
    all_transcripts = set()
    sample_data = {}
    
    # First pass: collect all transcripts and sample data
    for i, sample in enumerate(sample_folders, 1):
        print_progress_bar(i, len(sample_folders), f"Reading {file_type}", "🧬")
        
        sample_dir = os.path.join(results_dir, sample)
        rsem_files = glob.glob(os.path.join(sample_dir, "*_rsem.isoforms.results"))
        
        if not rsem_files:
            print_warning(f"No RSEM isoforms file found for {sample}")
            sample_data[sample] = {}
            continue
            
        rsem_file = rsem_files[0]
        print_processing(sample, f"RSEM isoforms: {os.path.basename(rsem_file)}")
        
        df = pd.read_csv(rsem_file, sep='\t')
        
        sample_values = {}
        for _, row in df.iterrows():
            transcript_id = row['transcript_id']
            value = row[value_column]
            
            # Convert transcript_id to transcript_name
            transcript_name = transcript_id_to_name.get(transcript_id, transcript_id)
            all_transcripts.add(transcript_name)
            sample_values[transcript_name] = value
        
        sample_data[sample] = sample_values
    
    print_info(f"🧬 Found {len(all_transcripts):,} unique transcripts across all samples")
    
    # Second pass: create complete matrix with zeros for missing transcripts
    print_info("🔧 Creating complete matrix (filling missing transcripts with 0)...")
    merged_data = {}
    
    for transcript_name in all_transcripts:
        merged_data[transcript_name] = {}
        for sample in sample_folders:
            # Use 0 if transcript is missing from this sample
            merged_data[transcript_name][sample] = sample_data.get(sample, {}).get(transcript_name, 0)
    
    # Convert to DataFrame
    result_df = pd.DataFrame(merged_data).T
    result_df = result_df.reindex(columns=sample_folders, fill_value=0)
    
    # Report missing transcript statistics
    missing_stats = {}
    for sample in sample_folders:
        missing_count = sum(1 for transcript in all_transcripts if sample_data.get(sample, {}).get(transcript, 0) == 0)
        missing_stats[sample] = missing_count
    
    if missing_stats:
        max_missing = max(missing_stats.values())
        min_missing = min(missing_stats.values())
        avg_missing = sum(missing_stats.values()) / len(missing_stats)
        
        print_info(f"📊 Missing transcript statistics for {file_type}:")
        print_info(f"   • Max missing per sample: {max_missing:,} transcripts")
        print_info(f"   • Min missing per sample: {min_missing:,} transcripts") 
        print_info(f"   • Average missing per sample: {avg_missing:.1f} transcripts")
        
        if max_missing > 0:
            print_warning(f"Some samples are missing up to {max_missing:,} transcripts - filled with zeros")
    
    print_success(f"Created {file_type} matrix: {len(result_df):,} transcripts × {len(sample_folders)} samples", "📊")
    return result_df

def calculate_isoform_percentages(results_dir, sample_folders, transcript_id_to_name, 
                                transcript_id_to_gene_id):
    """
    📊 Calculate isoform percentages relative to gene expression
    """
    # Collect all unique transcripts across all samples first
    all_transcripts = set()
    sample_data = {}
    
    # First pass: collect all transcripts and calculate percentages per sample
    for i, sample in enumerate(sample_folders, 1):
        print_progress_bar(i, len(sample_folders), "Reading percentage data", "📊")
        
        sample_dir = os.path.join(results_dir, sample)
        rsem_files = glob.glob(os.path.join(sample_dir, "*_rsem.isoforms.results"))
        
        if not rsem_files:
            print_warning(f"No RSEM isoforms file found for {sample}")
            sample_data[sample] = {}
            continue
            
        rsem_file = rsem_files[0]
        print_processing(sample, f"Calculating percentages: {os.path.basename(rsem_file)}")
        
        df = pd.read_csv(rsem_file, sep='\t')
        
        # Calculate gene-level sums for this sample
        gene_sums = {}
        for _, row in df.iterrows():
            transcript_id = row['transcript_id']
            expected_count = row['expected_count']
            gene_id = transcript_id_to_gene_id.get(transcript_id)
            
            if gene_id:
                if gene_id not in gene_sums:
                    gene_sums[gene_id] = 0
                gene_sums[gene_id] += expected_count
        
        # Calculate percentages for this sample
        sample_percentages = {}
        for _, row in df.iterrows():
            transcript_id = row['transcript_id']
            expected_count = row['expected_count']
            gene_id = transcript_id_to_gene_id.get(transcript_id)
            
            transcript_name = transcript_id_to_name.get(transcript_id, transcript_id)
            all_transcripts.add(transcript_name)
            
            # Calculate percentage
            if gene_id and gene_sums[gene_id] > 0:
                percentage = (expected_count / gene_sums[gene_id]) * 100
            else:
                percentage = 0
            
            sample_percentages[transcript_name] = percentage
        
        sample_data[sample] = sample_percentages
    
    print_info(f"🧬 Found {len(all_transcripts):,} unique transcripts across all samples")
    
    # Second pass: create complete matrix with zeros for missing transcripts
    print_info("🔧 Creating complete percentage matrix (filling missing transcripts with 0)...")
    merged_data = {}
    
    for transcript_name in all_transcripts:
        merged_data[transcript_name] = {}
        for sample in sample_folders:
            # Use 0 if transcript is missing from this sample
            merged_data[transcript_name][sample] = sample_data.get(sample, {}).get(transcript_name, 0)
    
    # Convert to DataFrame
    result_df = pd.DataFrame(merged_data).T
    result_df = result_df.reindex(columns=sample_folders, fill_value=0)
    
    # Report missing transcript statistics
    missing_stats = {}
    for sample in sample_folders:
        missing_count = sum(1 for transcript in all_transcripts if sample_data.get(sample, {}).get(transcript, 0) == 0)
        missing_stats[sample] = missing_count
    
    if missing_stats:
        max_missing = max(missing_stats.values())
        min_missing = min(missing_stats.values())
        avg_missing = sum(missing_stats.values()) / len(missing_stats)
        
        print_info(f"📊 Missing transcript statistics for percentages:")
        print_info(f"   • Max missing per sample: {max_missing:,} transcripts")
        print_info(f"   • Min missing per sample: {min_missing:,} transcripts") 
        print_info(f"   • Average missing per sample: {avg_missing:.1f} transcripts")
        
        if max_missing > 0:
            print_warning(f"Some samples are missing up to {max_missing:,} transcripts - filled with zeros")
    
    print_success(f"Created percentage matrix: {len(result_df):,} transcripts × {len(sample_folders)} samples", "📊")
    return result_df

def diagnose_missing_data(results_dir, sample_folders, output_dir):
    """
    🔍 Diagnose and report missing data patterns across samples
    """
    print_step(9, "Diagnosing Missing Data Patterns", "🔍")
    
    # Analyze featureCounts files
    print_info("📊 Analyzing STAR featureCounts completeness...")
    star_gene_counts = {}
    
    for sample in sample_folders:
        sample_dir = os.path.join(results_dir, sample)
        featurecounts_files = glob.glob(os.path.join(sample_dir, "*_featureCounts.txt"))
        
        if featurecounts_files:
            df = pd.read_csv(featurecounts_files[0], sep='\t', comment='#')
            star_gene_counts[sample] = len(df)
        else:
            star_gene_counts[sample] = 0
    
    # Analyze RSEM gene files
    print_info("🧬 Analyzing RSEM gene completeness...")
    rsem_gene_counts = {}
    
    for sample in sample_folders:
        sample_dir = os.path.join(results_dir, sample)
        rsem_files = glob.glob(os.path.join(sample_dir, "*_rsem.genes.results"))
        
        if rsem_files:
            df = pd.read_csv(rsem_files[0], sep='\t')
            rsem_gene_counts[sample] = len(df)
        else:
            rsem_gene_counts[sample] = 0
    
    # Analyze RSEM isoform files
    print_info("🧬 Analyzing RSEM isoform completeness...")
    rsem_isoform_counts = {}
    
    for sample in sample_folders:
        sample_dir = os.path.join(results_dir, sample)
        rsem_files = glob.glob(os.path.join(sample_dir, "*_rsem.isoforms.results"))
        
        if rsem_files:
            df = pd.read_csv(rsem_files[0], sep='\t')
            rsem_isoform_counts[sample] = len(df)
        else:
            rsem_isoform_counts[sample] = 0
    
    # Create diagnostic report
    diagnostic_data = {
        'Sample': sample_folders,
        'STAR_Genes': [star_gene_counts.get(s, 0) for s in sample_folders],
        'RSEM_Genes': [rsem_gene_counts.get(s, 0) for s in sample_folders],
        'RSEM_Isoforms': [rsem_isoform_counts.get(s, 0) for s in sample_folders]
    }
    
    diagnostic_df = pd.DataFrame(diagnostic_data)
    
    # Calculate statistics
    print_info("📈 Sample completeness statistics:")
    for col in ['STAR_Genes', 'RSEM_Genes', 'RSEM_Isoforms']:
        values = diagnostic_df[col]
        print_info(f"   • {col}:")
        print_info(f"     - Range: {values.min():,} - {values.max():,}")
        print_info(f"     - Mean: {values.mean():.0f}")
        print_info(f"     - Std: {values.std():.0f}")
    
    # Save diagnostic report
    diagnostic_file = os.path.join(output_dir, "diagnostic_sample_completeness.csv")
    diagnostic_df.to_csv(diagnostic_file, index=False)
    print_success(f"Saved diagnostic report: {os.path.basename(diagnostic_file)}", "📋")
    
    # Identify problematic samples
    star_threshold = diagnostic_df['STAR_Genes'].quantile(0.25)  # Lower quartile
    rsem_gene_threshold = diagnostic_df['RSEM_Genes'].quantile(0.25)
    rsem_isoform_threshold = diagnostic_df['RSEM_Isoforms'].quantile(0.25)
    
    problematic_samples = []
    for _, row in diagnostic_df.iterrows():
        issues = []
        if row['STAR_Genes'] < star_threshold:
            issues.append("Low STAR genes")
        if row['RSEM_Genes'] < rsem_gene_threshold:
            issues.append("Low RSEM genes")
        if row['RSEM_Isoforms'] < rsem_isoform_threshold:
            issues.append("Low RSEM isoforms")
        
        if issues:
            problematic_samples.append((row['Sample'], issues))
    
    if problematic_samples:
        print_warning("🚨 Samples with potential data completeness issues:")
        for sample, issues in problematic_samples:
            print_warning(f"   • {sample}: {', '.join(issues)}")
    else:
        print_success("All samples have consistent data completeness! ✨")
    
    return diagnostic_df

def save_csv_with_style(df, filepath, file_description):
    """
    💾 Save DataFrame as CSV with style and progress indication
    """
    print_info(f"💾 Saving {file_description}...")
    df.to_csv(filepath, index=True)
    file_size = os.path.getsize(filepath) / (1024 * 1024)  # MB
    print_success(f"Saved: {os.path.basename(filepath)} ({file_size:.2f} MB)", "💾")

def main():
    """
    🚀 Main function to orchestrate the RNA-seq merging process
    """
    # Print beautiful banner
    print_banner()
    
    # Configuration with colorful display
    print_step(0, "Configuration Setup", "⚙️")
    results_dir = "/mnt/d/mad_rna_seq/results"
    gtf_file_path = "/mnt/d/genome/gencode.vM37.chr_patch_hapl_scaff.annotation.gtf"
    output_dir = "/mnt/d/mad_rna_seq/merged_results"
    
    print_info(f"📁 Results directory: {results_dir}")
    print_info(f"📚 GTF file: {gtf_file_path}")
    print_info(f"💾 Output directory: {output_dir}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    print_success(f"Output directory ready: {output_dir}")
    
    try:
        # Step 1: Parse GTF file for annotations
        gene_id_to_name, transcript_id_to_name, transcript_id_to_gene_id = create_gtf_mappings(gtf_file_path)
        
        # Step 2: Find sample folders
        sample_folders = get_sample_folders(results_dir)
        if not sample_folders:
            print_error("No valid sample folders found! Exiting...")
            return
        
        # Create output files with progress tracking
        output_files = [
            ("1_merged_star_gene_counts.csv", "🌟 STAR Gene Counts"),
            ("2_merged_rsem_gene_expected_counts.csv", "🧬 RSEM Gene Expected Counts"),
            ("3_merged_rsem_isoform_expected_counts.csv", "🧬 RSEM Isoform Expected Counts"),
            ("4_merged_rsem_isoform_percentages.csv", "📊 RSEM Isoform Percentages"),
            ("5_merged_rsem_gene_fpkm.csv", "📈 RSEM Gene FPKM"),
            ("6_merged_rsem_gene_tpm.csv", "📈 RSEM Gene TPM")
        ]
        
        # Step 3: Process files and create merged matrices
        
        # File 1: STAR gene counts (featureCounts)
        print_step(3, output_files[0][1], "🌟")
        star_counts = process_featurecounts_files(results_dir, sample_folders, gene_id_to_name)
        output_file = os.path.join(output_dir, output_files[0][0])
        save_csv_with_style(star_counts, output_file, output_files[0][1])
        
        # File 2: RSEM expected gene counts
        print_step(4, output_files[1][1], "🧬")
        rsem_gene_counts = process_rsem_gene_files(results_dir, sample_folders, gene_id_to_name, 
                                                 'expected_count', 'RSEM Gene Counts')
        output_file = os.path.join(output_dir, output_files[1][0])
        save_csv_with_style(rsem_gene_counts, output_file, output_files[1][1])
        
        # File 3: RSEM isoform expected counts
        print_step(5, output_files[2][1], "🧬")
        rsem_isoform_counts = process_rsem_isoform_files(results_dir, sample_folders, 
                                                       transcript_id_to_name, transcript_id_to_gene_id, 
                                                       'expected_count', 'RSEM Isoform Counts')
        output_file = os.path.join(output_dir, output_files[2][0])
        save_csv_with_style(rsem_isoform_counts, output_file, output_files[2][1])
        
        # File 4: RSEM isoform percentages
        print_step(6, output_files[3][1], "📊")
        rsem_isoform_percentages = calculate_isoform_percentages(results_dir, sample_folders,
                                                               transcript_id_to_name, transcript_id_to_gene_id)
        output_file = os.path.join(output_dir, output_files[3][0])
        save_csv_with_style(rsem_isoform_percentages, output_file, output_files[3][1])
        
        # File 5: RSEM gene FPKM
        print_step(7, output_files[4][1], "📈")
        rsem_gene_fpkm = process_rsem_gene_files(results_dir, sample_folders, gene_id_to_name, 
                                               'FPKM', 'RSEM Gene FPKM')
        output_file = os.path.join(output_dir, output_files[4][0])
        save_csv_with_style(rsem_gene_fpkm, output_file, output_files[4][1])
        
        # File 6: RSEM gene TPM
        print_step(8, output_files[5][1], "📈")
        rsem_gene_tpm = process_rsem_gene_files(results_dir, sample_folders, gene_id_to_name, 
                                              'TPM', 'RSEM Gene TPM')
        output_file = os.path.join(output_dir, output_files[5][0])
        save_csv_with_style(rsem_gene_tpm, output_file, output_files[5][1])
        
        # Step 9: Diagnose missing data patterns
        diagnostic_df = diagnose_missing_data(results_dir, sample_folders, output_dir)
        
        # Beautiful summary
        print(f"\n{Colors.BOLD}{Colors.GREEN}{'='*80}{Colors.ENDC}")
        print(f"{Colors.BOLD}{Colors.GREEN}🎉 SUCCESS! RNA-seq Merger Completed Successfully! 🎉{Colors.ENDC}")
        print(f"{Colors.BOLD}{Colors.GREEN}{'='*80}{Colors.ENDC}")
        
        print(f"\n{Colors.CYAN}📊 Processing Summary:{Colors.ENDC}")
        print_success(f"Processed {len(sample_folders)} samples")
        print_success(f"Created 6 beautiful CSV files + 1 diagnostic report")
        print_success(f"Output location: {output_dir}")
        
        print(f"\n{Colors.YELLOW}📁 Files Created:{Colors.ENDC}")
        for filename, description in output_files:
            filepath = os.path.join(output_dir, filename)
            file_size = os.path.getsize(filepath) / (1024 * 1024)  # MB
            print_info(f"  {description}")
            print(f"    {Colors.CYAN}📄 {filename} ({file_size:.2f} MB){Colors.ENDC}")
        
        # Show diagnostic file
        diagnostic_file = "diagnostic_sample_completeness.csv"
        diagnostic_path = os.path.join(output_dir, diagnostic_file)
        if os.path.exists(diagnostic_path):
            file_size = os.path.getsize(diagnostic_path) / (1024)  # KB
            print_info(f"  🔍 Sample Completeness Diagnostic")
            print(f"    {Colors.CYAN}📄 {diagnostic_file} ({file_size:.1f} KB){Colors.ENDC}")
        
        print(f"\n{Colors.MAGENTA}✨ Ready for downstream analysis! ✨{Colors.ENDC}")
        print(f"\n{Colors.YELLOW}💡 Pro Tips:{Colors.ENDC}")
        print_info("📋 Check diagnostic_sample_completeness.csv for data quality insights")
        print_info("🔍 Missing values (now filled with 0) indicate genes not detected in those samples")
        print_info("📊 This is normal RNA-seq behavior - some genes aren't expressed in all conditions")
        
        # Show completion time
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"\n{Colors.BLUE}🕒 Completed at: {current_time}{Colors.ENDC}")
        
    except Exception as e:
        print_error(f"An error occurred: {str(e)}")
        import traceback
        print(f"{Colors.FAIL}{traceback.format_exc()}{Colors.ENDC}")

if __name__ == "__main__":
    main()