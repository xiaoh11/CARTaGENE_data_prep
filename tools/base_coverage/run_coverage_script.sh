#!/bin/bash
#SBATCH --job-name=coverage_full
#SBATCH --output=coverage.log

#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=25G
#SBATCH --time=24:00:00

source ~/ENV/bin/activate
module load tabix


# Set the input and output directories
input_dir=~/projects/rrg-vmooser/xiaoh11/bravo/data/basis/vcfs
output_dir=~/projects/rrg-vmooser/xiaoh11/bravo/data/runtime/coverage/full
prune025_dir=~/projects/rrg-vmooser/xiaoh11/bravo/data/runtime/coverage/bin_0.25
prune050_dir=~/projects/rrg-vmooser/xiaoh11/bravo/data/runtime/coverage/bin_0.50
prune075_dir=~/projects/rrg-vmooser/xiaoh11/bravo/data/runtime/coverage/bin_0.75
prune100_dir=~/projects/rrg-vmooser/xiaoh11/bravo/data/runtime/coverage/bin_1.00

# Loop over the chromosomes
for chr in {1..22} X; do
  input_file="$input_dir/chrom${chr}_output.vep.vcf.gz"
  output_file="$output_dir/${chr}.full.tsv.gz"

  pruned025_output_file="$prune025_dir/${chr}.bin_0.25.tsv.gz"
  pruned050_output_file="$prune050_dir/${chr}.bin_0.50.tsv.gz"
  pruned075_output_file="$prune075_dir/${chr}.bin_0.75.tsv.gz"
  pruned100_output_file="$prune100_dir/${chr}.bin_1.00.tsv.gz"

  # Run the command for each chromosome
  python create_coverage_vcf.py -i "$input_file" | bgzip -c > "$output_file"
  # tabix index
  tabix -s 1 -b 2 -e 3 "$output_file"

  python ../py_tools/prune.py -i "$output_file" -l 0.25 -o "$pruned025_output_file"
  tabix -s 1 -b 2 -e 3 "$pruned025_output_file"

  python ../py_tools/prune.py -i "$output_file" -l 0.50 -o "$pruned050_output_file"
  tabix -s 1 -b 2 -e 3 "$pruned050_output_file"

  python ../py_tools/prune.py -i "$output_file" -l 0.75 -o "$pruned075_output_file"
  tabix -s 1 -b 2 -e 3 "$pruned075_output_file"
  
  python ../py_tools/prune.py -i "$output_file" -l 1.00 -o "$pruned100_output_file"
  tabix -s 1 -b 2 -e 3 "$pruned100_output_file"
done

