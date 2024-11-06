# Download reference fasta files (~1.5GB)
wget --timestamping http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
mkdir -p reference-files/fasta && tar -xzvf hg38.chromFa.tar.gz
for i in {1..22}; do
    mv chroms/chr${i}.fa reference-files/fasta/hg38_chr${i}.fa
done
rm -rf chroms

# Build docker image for simulation
docker build -t simulation simulation

# Build Docker image named 'toy-dataset-training' from the 'training' directory
docker build -t training training

# Loop through chromosome pairs (1-2, 3-4, etc. up to 19-20-21-22)
# Small chromosomes are processed together
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

