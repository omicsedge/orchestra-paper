
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
    cat header chr$chr.fa > data/fasta/$chr.fa
    
    # Cleanup
    rm header chr$chr.fa chr$chr.fna
done
```

### 2. Reference Files Structure

The following files are required:

| File/Directory | Description |
|---------------|-------------|
| `data/fasta/*.fa` | Chromosome FASTA files (chr1-22) |
| `reference_files/QC_ancestry_associated_SNP_set.hg38.keep` | QC ancestry-associated SNP set |
| `reference_files/Ancestry_regions.hg38.txt` | Ancestry regions defined for our custom panel |
| `data/toy_example/Source_panel.vcf.gz` | Source panel for simulation |
| `data/toy_example/SampleTable.forTraining.txt` | Population structure definitions |
| `data/toy_example/Admixed_Mexicans.target_panel.vcf.gz` | Test data for inference |

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
<summary><strong>Source Panel</strong> (Source_panel.vcf.gz)</summary>

- Compressed VCF format
- Contains genetic variants for simulation
- Required fields: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
</details>

<details>
<summary><strong>Sample Map</strong> (SampleTable.forTraining.txt)</summary>

- Tab-separated values
- Defines population structure
- Required columns: Sample ID, Population, Super Population
</details>

### 3. Set Environment Variable

```bash
export EXPERIMENT_NAME="example-0.01"
```

### 4. Build Docker Images

```bash
make build
```

### 5. Run Simulation Pipeline

```bash
docker run --rm \
    -v $(pwd)/data:/data \
    -v $(pwd)/results:/results \
    orchestra simulation\
    -sc 1 -ec 22 \
    -sp /data/toy_example/Source_panel.vcf.gz \
    -sm /data/toy_example/SampleTable.forTraining.txt \
    -v $EXPERIMENT_NAME \
    -t "random" \
    -nt 2 \
    -o /results/simulation
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

### 6. Train Models

Process chromosomes in pairs (smallest 19-22 chromosomes grouped together):

```bash
for chr in "1 2" "3 4" "5 6" "7 8" "9 10" "11 12" "13 14" "15 16" "17 18" "19 22"; do
    start_chr=$(echo $chr | cut -d' ' -f1)
    end_chr==$(echo $chr | cut -d' ' -f2)

    docker run --rm \
        -v $(pwd)/data:/data \
        -v $(pwd)/results:/results \
        orchestra training \
        -sd /results/simulation \
        -sc $start_chr \
        -ec $end_chr \
        -ws 600 \
        -l 3 \
        -v $EXPERIMENT_NAME \
        -e 100 \
        -o /results/training 

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

### 7. Run Inference Pipeline

```bash
docker run --rm \
    -v $(pwd)/data:/data \
    -v $(pwd)/results:/results \
    orchestra inference \
    -p /data/toy_example/Admixed_Mexicans.target_panel.vcf.gz \
    -m /results/training/$EXPERIMENT_NAME \
    -o /results/inference
```

<details>
<summary><strong>Inference Parameters</strong></summary>

| Parameter | Description |
|-----------|-------------|
| `-p`, `--panel` | Inference panel path |
| `-o`, `--output` | Output directory |
| `-m`, `--model` | Trained model directory |
</details>

## References
* Lerga-Jaso, J., Novković, B., Unnikrishnan, D., Bamunusinghe, V., Hatorangan, M.R., Manson, C., Pedersen, H., Osama, A., Terpolovsky, A., Bohn, S., De Marino, A., Mahmoud, A.A., Bircan, K.O., Khan, U., Grabherr, M.G., Yazdi, P.G. Retracing Human Genetic Histories and Natural Selection Using Precise Local Ancestry Inference. bioRxiv 2023.09.11.557177; doi: https://doi.org/10.1101/2023.09.11.557177
* Cuadros-Espinoza, S., Laval, G., Quintana-Murci, L., Patin, E. The genomic signatures of natural selection in admixed human populations. Am. J. Hum. Genet. 109, 710-726 (2022). doi: 10.1016/j.ajhg.2022.02.011; pmid: 35259336

# Non-Commercial Use License
### Version 1.0

## NOTICE
This software is provided free of charge for **academic research use only**. Any use by **commercial entities, for-profit organizations, or consultants** is strictly prohibited without prior authorization. For inquiries about commercial licensing, contact **pyazdi@gmail.com**.
