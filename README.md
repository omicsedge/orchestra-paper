# Orchestra Workflow: simulation, training and inference
## Paper: https://www.biorxiv.org/content/10.1101/2023.09.11.557177v1

## Setup Instructions

### 1. Download and Prepare Reference Files

Download and process reference FASTA files (~600MB) from NCBI, converting them to the format required for SLiM simulation:

```bash
for chr in {1..22}; do
    # Download chromosome files from NCBI
    wget --timestamping https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr$chr.fna.gz
    
    # Decompress and process the files
    bgzip -d chr$chr.fna.gz
    head -n 1 chr$chr.fna > header
    sed -e 's/[NRWYMKBSV]/A/g' -e '1d' -e 'y/actg/ACTG/' chr$chr.fna > chr$chr.fa
    cat header chr$chr.fa > reference_files/fasta/$chr.fa
    
    # Clean up temporary files
    rm header chr$chr.fa chr$chr.fna
done
```

### 2. Build Docker Images

Build the required Docker images for simulation, training, and inference:

```bash
make build
```

### 3. Run Simulation Pipeline

Execute the simulation using Docker with the following parameters:

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

### 4. Train Models

Process chromosomes in pairs (with chromosomes 19-22 grouped together) using the training pipeline:

```bash
mkdir -p output_training
for chr in "1 2" "3 4" "5 6" "7 8" "9 10" "11 12" "13 14" "15 16" "17 18" "19 22"; do
    # Extract chromosome range
    start_chr=$(echo $chr | cut -d' ' -f1)
    end_chr=$(echo $chr | cut -d' ' -f2)

    # Run training container
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

    echo "Done with chromosomes $start_chr-$end_chr, model saved in output_training/"
done
```

Each training run processes a pair of chromosomes with the following parameters:
- Window size: 600
- Training level: 3
- Validation split: 1%
- Output directory: ./output_training

### 5. Run Inference Pipeline

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
