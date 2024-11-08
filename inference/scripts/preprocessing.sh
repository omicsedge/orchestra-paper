#!/bin/bash
set -e  # Exit immediately if any command fails

# Accept command line arguments
inference_panel=$1  # First argument: path to the inference panel VCF file
chr=$2              # Second argument: chromosome number (appears unused in script)

# Index the VCF file
bcftools index -f $inference_panel

# Print the number of samples in the inference panel
echo "Samples for inference: $(bcftools query -l $inference_panel | wc -l)"

# Create and move to a temporary directory
temp_dir=$(mktemp -d)
cd $temp_dir
echo "Working directory: $temp_dir"

# Create necessary subdirectories
mkdir -p chromosomes samples commands

# Split VCF by chromosome (1-22) and process in parallel
for chr in {1..22}; do
    nohup bcftools view -Oz -r $chr $inference_panel -o chromosomes/"chr$chr.vcf.gz" && bcftools index chromosomes/"chr$chr.vcf.gz" &
done

wait  # Wait for all chromosome splitting jobs to complete

# Extract list of all samples
bcftools query -l $inference_panel > all.txt

# Split samples into 5 roughly equal groups
split -l $(( ($(wc -l < all.txt) + 4) / 5 )) all.txt samples/samples_

# Rename split files with numeric suffixes (01-05)
a=1
for file in samples/samples_*; do
  mv "$file" "samples/samples_$(printf "%02d" $a).csv"
  let a++
done

# Print sample counts for first and last batch as sanity check
echo "Number of sampels in first batch: " $(cat samples/samples_01.csv | wc -l)
echo "Number of sampels in last batch: " $(cat samples/samples_05.csv | wc -l)

# Generate commands to split VCF by sample groups
[ -f "commands.sh" ] && rm "commands.sh"
for i in {1..5}; do
  path="dataset/sample_list_$i"
  mkdir -p $path
  for chr in {1..22}; do
      # Create bcftools command to extract samples for each chromosome
      echo "nohup bcftools view -S samples/samples_0${i}.csv -Oz chromosomes/"chr$chr.vcf.gz" > $path/chr${chr}.vcf.gz &" >> commands.sh
  done
done

# Split commands into 5 separate shell scripts for parallel execution
rm command_*  # Remove any existing command files
split -l $(( ($(wc -l < commands.sh) + 4) / 5 )) commands.sh command_

# Process each command file
for i in command_*; do
   # Add shebang and make executable
   echo '#!/bin/bash' > aux
   cat $i >> aux
   mv aux $i
   chmod +x $i

   # Execute the command file
   echo "Start shell script ./$i"
   ./$i
done

wait  # Wait for all processes to complete
