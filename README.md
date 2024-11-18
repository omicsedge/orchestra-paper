# Orchestra Workflow: Simulation, Training and Inference

[![Paper](https://img.shields.io/badge/bioRxiv-10.1101%2F2023.09.11.557177v1-blue)](https://www.biorxiv.org/content/10.1101/2023.09.11.557177v1)

Orchestra is a pipeline for genetic ancestry inference using deep learning. This README provides setup and usage instructions.

## Setup Instructions

### 1. Download and Prepare Reference Files

Download and process reference FASTA files (~600MB) from NCBI:

```bash
for chr in {1..22}; do
    # Download chromosome files
    wget --timestamping https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr$chr.fna.gz
    
    # Process files
    bgzip -d chr$chr.fna.gz
    head -n 1 chr$chr.fna > header
    sed -e 's/[NRWYMKBSV]/A/g' -e '1d' -e 'y/actg/ACTG/' chr$chr.fna > chr$chr.fa
    cat header chr$chr.fa > reference_files/fasta/$chr.fa
    
    # Cleanup
    rm header chr$chr.fa chr$chr.fna
done
```

### 2. Reference Files Structure

The following files are required:

| File/Directory | Description |
|---------------|-------------|
| `reference_files/fasta/*.fa` | Chromosome FASTA files (chr1-22) |
| `reference_files/QC_ancestry_associated_SNP_set.hg38.keep` | QC ancestry-associated SNP set |
| `reference_files/Ancestry_regions.hg38.txt` | Ancestry regions definition file |
| `reference_files/example_data/source_panel.vcf.gz` | Source panel for simulation |
| `reference_files/example_data/sample_map.tsv` | Population structure definitions |
| `reference_files/example_data/inference_panel.gz` | Test data for inference |

#### File Format Specifications

<details>
<summary><strong>FASTA Files</strong> (*.fa)</summary>

- One file per chromosome (chr1.fa to chr22.fa)
- Contains reference genome sequence
- Processed from NCBI reference files
</details>

<details>
<summary><strong>Ancestry SNP Set</strong> (QC_ancestry_associated_SNP_set.hg38.keep)</summary>

- Space-separated text file
- Contains 1,202,443 ancestry-informative variants
- Format: `CHROM POS REF ALT`
- Based on hg38/GRCh38 assembly
</details>

<details>
<summary><strong>Ancestry Regions</strong> (Ancestry_regions.hg38.txt)</summary>

- Tab-separated text file
- Defines genomic regions for ancestry analysis
- Format: `CHROM START_POS END_POS`
- Based on hg38/GRCh38 assembly
</details>

<details>
<summary><strong>Source Panel</strong> (source_panel.vcf.gz)</summary>

- Compressed VCF format
- Contains genetic variants for simulation
- Required fields: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
</details>

<details>
<summary><strong>Sample Map</strong> (sample_map.tsv)</summary>

- Tab-separated values
- Defines population structure
- Required columns: Sample ID, Population, Super Population
</details>

### 3. Build Docker Images

```bash
make build
```

### 4. Run Simulation Pipeline

```bash
mkdir -p output_simulation
docker run --rm \
    -v $(pwd)/reference_files:/reference_files \
    -v $(pwd)/output_simulation:/output_simulation \
    simulation \
    -sc 1 -ec 22 \
    -sp /reference_files/example_data/source_panel.vcf.gz \
    -sm /reference_files/example_data/sample_map.tsv \
    -v "example-0.01" \
    -t "random" \
    -nt 2 \
    -o /output_simulation
```

<details>
<summary><strong>Simulation Parameters</strong></summary>

| Parameter | Description |
|-----------|-------------|
| `-sc`, `--start-chromosome` | Start chromosome number |
| `-ec`, `--end-chromosome` | End chromosome number |
| `-sp`, `--source-panel` | Source panel VCF path |
| `-sm`, `--sample-map` | Sample map TSV path |
| `-v`, `--version` | Version identifier |
| `-t`, `--type` | Simulation type |
| `-nt`, `--num-threads` | Number of threads |
| `-o`, `--output` | Output directory |
</details>

### 5. Train Models

Process chromosomes in pairs (smallest 19-22 chromosomes grouped together):

```bash
mkdir -p output_training
for chr in "1 2" "3 4" "5 6" "7 8" "9 10" "11 12" "13 14" "15 16" "17 18" "19 22"; do
    start_chr=$(echo $chr | cut -d' ' -f1)
    end_chr=$(echo $chr | cut -d' ' -f2)

    docker run --rm \
        -v $(pwd)/reference_files:/reference_files \
        -v $(pwd)/output_simulation:/output_simulation \
        -v $(pwd)/output_training:/output_training \
        training \
        -sd /output_simulation \
        -sc $start_chr \
        -ec $end_chr \
        -ws 600 \
        -l 3 \
        -o /output_training \
        -v "example-0.01" \
        -e 100

    echo "✓ Completed chromosomes $start_chr-$end_chr"
done
```

<details>
<summary><strong>Training Parameters</strong></summary>

| Parameter | Description |
|-----------|-------------|
| `-sd`, `--simulation-dir` | Simulation data directory |
| `-sc`, `--start-chromosome` | Start chromosome |
| `-ec`, `--end-chromosome` | End chromosome |
| `-ws`, `--window-size` | Processing window size |
| `-l`, `--level` | Model complexity level |
| `-o`, `--output` | Output directory |
| `-v`, `--version` | Version identifier |
| `-e`, `--epochs` | Training epochs |
</details>

### 6. Run Inference Pipeline

```bash
mkdir -p output_inference
docker run --rm \
    -v $(pwd)/reference_files:/reference_files \
    -v $(pwd)/output_training:/output_training \
    -v $(pwd)/output_inference:/output_inference \
    inference \
    -p /reference_files/example_data/inference_panel.gz \
    -o /output_inference \
    -m /output_training/example-0.01
```

<details>
<summary><strong>Inference Parameters</strong></summary>

| Parameter | Description |
|-----------|-------------|
| `-p`, `--panel` | Inference panel path |
| `-o`, `--output` | Output directory |
| `-m`, `--model` | Trained model directory |
</details>
