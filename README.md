# ancestry-paper

## Setup Instructions

### 1. Download and Prepare Reference Files

# Download reference fasta files (~600MB)
# Adjust them for SLiM simulation
```
for chr in {1..22}; do
    wget --timestamping https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr$chr.fna.gz
    bgzip -d chr$chr.fna.gz
    head -n 1 chr$chr.fna > header
    sed -e 's/[NRWYMK]/A/g' -e '1d' -e 'y/actg/ACTG/' chr$chr.fna > chr$chr.fa
    cat header chr$chr.fa > reference_files/fasta/$chr.fa
    rm header chr$chr.fa chr$chr.fna
done
```

# Build docker images for simulation, training, and inference
# Docker images are built in the Makefile with tags 'simulation', 'training', and 'inference'
```
make build
```
# Run simulation pipeline

```
docker run --rm -v $(pwd)/reference_files:/reference_files -v $(pwd)/output:/output simulation \
        -sc 1 -ec 21 \
        -sp /reference_files/example_data/source_panel.vcf.gz \
        -sm /reference_files/example_data/sample_map.tsv \
        -d "Peruvians 0.01" \
        -t "random" \
        -nt 2 \
        -o /output
```

# Loop through chromosome pairs (1-2, 3-4, etc. up to 19-20-21-22)
# Small chromosomes are processed together
```
for chr in "1 2" "3 4" "5 6" "7 8" "9 10" "11 12" "13 14" "15 16" "17 18" "19 20 21 22"; do
    # Extract the start and end chromosome numbers from the pair
    start_chr=$(echo $chr | cut -d' ' -f1)
    end_chr=$(echo $chr | cut -d' ' -f2)

    # Docker parameters:
    #   --rm: Remove container after execution
    #   -v: Mount local sample-data directory to /training in container
    #
    # Training script parameters:
    #   -sd: Training dataset directory
    #   -sc: Start chromosome
    #   -ec: End chromosome
    #   -ws: Window size (600)
    #   -l: Level (3)
    #   -o: Output directory
    #   -v: Validation split (0.01)
    docker run --rm -v $(pwd)/sample-data:/training toy-dataset-training \
        python train.py -sd training/toy-training-dataset -sc $start_chr -ec $end_chr -ws 600 -l 3 -o output -v 0.01

    # Print completion message for current chromosome pair
    echo "Done with $start_chr $end_chr, model saved in output/"
done

```