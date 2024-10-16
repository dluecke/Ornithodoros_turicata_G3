import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
from PyPDF2 import PdfMerger

# Function to run a command and check its output
def run_command(command):
    result = subprocess.run(command, shell=True, check=True, text=True, capture_output=True)
    if result.returncode != 0:
        print(f"Command failed with exit code {result.returncode}")
        print(result.stdout)
        print(result.stderr)
    return result

# Parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description='Align reads to reference genome and plot depth of coverage.')
    parser.add_argument('-r', '--reference', required=True, help='Reference genome file')
    parser.add_argument('-s', '--samples', required=True, help='File with sample names and paths to filtered fastq files')
    parser.add_argument('-c', '--chromosomes', nargs='+', required=True, help='List of chromosomes to include in analysis')
    return parser.parse_args()

# Read sample data from file
def read_sample_data(samples_file):
    samples = []
    with open(samples_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            sample_name = parts[0]
            read_files = parts[1:]
            samples.append((sample_name, read_files))
    return samples

# Calculate average coverage in 1Mb windows
def calculate_window_coverage(df, window_size=1000000):
    df = df.copy()  # Ensure we are working with a copy of the DataFrame
    df['window'] = (df['pos'] // window_size) * window_size
    window_coverage = df.groupby(['chr', 'window'])['depth'].mean().reset_index()
    return window_coverage

def main():
    args = parse_arguments()
    reference_genome = args.reference
    sample_file = args.samples
    chromosomes = set(args.chromosomes)

    samples = read_sample_data(sample_file)
    results_dir = "results"
    os.makedirs(results_dir, exist_ok=True)

    # Step 1: Align reads to the reference genome using minimap2
    for sample_name, filtered_files in samples:
        sam_output = f"{results_dir}/{sample_name}.sam"
        run_command(f"minimap2 -ax sr -t 20 {reference_genome} {' '.join(filtered_files)} > {sam_output}")

    # Step 2: Convert SAM to BAM, sort, and index
    for sample_name, _ in samples:
        bam_output = f"{results_dir}/{sample_name}.bam"
        sorted_bam_output = f"{results_dir}/{sample_name}_sorted.bam"
        run_command(f"samtools view -b -F 256 {results_dir}/{sample_name}.sam > {bam_output}")  # Filter for primary mappings
        run_command(f"samtools sort {bam_output} -o {sorted_bam_output}")
        run_command(f"samtools index {sorted_bam_output}")
        os.remove(f"{results_dir}/{sample_name}.sam")  # Remove the .sam file

    # Step 3: Generate mapping statistics using samtools flagstat
    for sample_name, _ in samples:
        bam_file = f"{results_dir}/{sample_name}_sorted.bam"
        flagstat_output = f"{results_dir}/{sample_name}_flagstat.txt"
        run_command(f"samtools flagstat {bam_file} > {flagstat_output}")

    # Step 3.1: Summarize flagstat results using multiqc
    run_command(f"multiqc {results_dir} -o {results_dir}")

    # Step 4: Generate depth of coverage data
    coverage_data = {}
    for sample_name, _ in samples:
        coverage_output = f"{results_dir}/{sample_name}_coverage.txt"
        run_command(f"samtools depth -a {sample_name}_sorted.bam > {coverage_output}")
        coverage_data[sample_name] = pd.read_csv(coverage_output, sep='\t', header=None, names=['chr', 'pos', 'depth'])

    # Step 5: Normalize the depth of coverage by average genome-wide depth
    normalized_coverage = {}
    average_coverage = []
    for sample_name in coverage_data:
        coverage = coverage_data[sample_name]
        genome_avg_depth = coverage[coverage['chr'].isin(chromosomes)]['depth'].mean()

        # Handle short contigs separately
        window_coverage = calculate_window_coverage(coverage[coverage['chr'].isin(chromosomes)])
        short_contigs = coverage.groupby('chr').filter(lambda x: x['pos'].max() < 1000000)
        short_contigs_avg = short_contigs.groupby('chr')['depth'].mean().reset_index()
        short_contigs_avg['window'] = 0
        short_contigs_avg = short_contigs_avg.rename(columns={'depth': 'normalized_depth'})
        short_contigs_avg['normalized_depth'] = short_contigs_avg['normalized_depth'] / genome_avg_depth

        window_coverage = window_coverage[~window_coverage['chr'].isin(short_contigs_avg['chr'])]
        window_coverage['normalized_depth'] = window_coverage['depth'] / genome_avg_depth

        normalized_coverage[sample_name] = pd.concat([window_coverage, short_contigs_avg])

        avg_depth = normalized_coverage[sample_name].groupby('chr')['normalized_depth'].mean().reset_index()
        avg_depth = avg_depth.rename(columns={'normalized_depth': sample_name})
        average_coverage.append(avg_depth)

    # Merge all average coverage dataframes on 'chr' column
    avg_coverage_df = average_coverage[0]
    for df in average_coverage[1:]:
        avg_coverage_df = pd.merge(avg_coverage_df, df, on='chr', how='outer')

    avg_coverage_df.to_csv(f"{results_dir}/normalized_coverage.csv", sep='\t', index=False)

    # Determine maximum position for each chromosome for setting x-axis limits
    max_positions = {}
    for chr in chromosomes:
        max_pos = max([coverage_data[sample_name][coverage_data[sample_name]['chr'] == chr]['pos'].max() for sample_name in coverage_data])
        max_positions[chr] = max_pos

    pdf_files = []
    for chromosome in chromosomes:
        plt.figure(figsize=(12, 6))
        for sample_name in normalized_coverage:
            chr_data = normalized_coverage[sample_name][normalized_coverage[sample_name]['chr'] == chromosome]
            if not chr_data.empty:
                plt.plot(chr_data['window'], chr_data['normalized_depth'], label=f'{sample_name}')
        
        plt.axhline(y=1, color='black', linestyle='dashed')
        plt.xlabel('Position (Mb)')
        plt.ylabel('Normalized Depth')
        plt.title(f'Depth of Coverage for Chromosome {chromosome}')
        plt.legend()
        plt.xlim(0, max_positions[chromosome])
        
        # Determine x-ticks at 10 Mb intervals and label every 50 Mb
        min_pos = 0
        max_pos = max_positions[chromosome]
        xticks = np.arange(min_pos, max_pos + 1, 10 * 1e6)
        xlabels = [f'{tick/1e6}' if tick % (50 * 1e6) == 0 else '' for tick in xticks] # Convert to Mb
        
        plt.xticks(xticks, labels=xlabels)  
        
        plt.ylim(0, 3)
        
        # Save the plot as SVG and PDF
        svg_path = f"{results_dir}/depth_of_coverage_{chromosome}.svg"
        pdf_path = f"{results_dir}/depth_of_coverage_{chromosome}.pdf"
        plt.savefig(svg_path)
        plt.savefig(pdf_path)
        pdf_files.append(pdf_path)
        
        plt.close()
    
    # Merge all individual PDF files into one
    merger = PdfMerger()
    for pdf in pdf_files:
        merger.append(pdf)
    
    merged_pdf_path = f"{results_dir}/merged_depth_of_coverage.pdf"
    merger.write(merged_pdf_path)
    merger.close()

if __name__ == "__main__":
    main()
