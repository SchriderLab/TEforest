#!/bin/bash
# Define a function to display usage and exit with an error code
usage() {
    echo "Usage: $0 -o <output_directory> -c <fasta_file> -r <ref_te_locations> -@ <threads> -g <reference genome fasta> -1 <fastq1> -2 <fastq2> -n genome_name"
    exit 1
}

# Check if the script was run with at least one argument
if [ $# -eq 0 ]; then
    usage
fi

while getopts o:c:@:d:g:1:2:r:n: flag
do
    case "${flag}" in
        o) output_directory=${OPTARG};;
        c) fasta_file=${OPTARG};;
        @) threads=${OPTARG};;
        g) ref_genome=${OPTARG};;
        1) fq1=${OPTARG};;
        2) fq2=${OPTARG};;
        r) ref_te_locations=${OPTARG};;
        n) genome_name=${OPTARG};;
    esac
done

# Check if threads is a number
if ! [[ "$threads" =~ ^[0-9]+$ ]]; then
    echo "Error: The '-@' option (threads) must be a valid number."
    usage
fi

# Function to check if a path is a valid directory
is_valid_directory() {
    if [ ! -d "$1" ]; then
        echo "Error: '$1' is not a valid directory."
        usage
    fi
}

# Function to check if a path is a valid file
is_valid_file() {
    if [ ! -f "$1" ]; then
        echo "Error: '$1' is not a valid file."
        usage
    fi
}

# Function to check if a file is a valid FASTA file
is_valid_fasta() {
    if ! head -n 1 "$1" | grep -q '^>'; then
        echo "Error: '$1' is not a valid FASTA file."
        usage
    fi
}

# Function to check if a file is a valid FASTQ file
is_valid_fastq() {
    if ! head -n 1 "$1" | grep -q '^@'; then
        echo "Error: '$1' is not a valid FASTQ file."
        usage
    fi
}

# Function to check if a path is a valid directory
create_valid_directory() {
    if [ ! -d "$1" ]; then
        echo "Creating directory: $1"
        mkdir -p "$1" || {
            echo "Error: Unable to create the directory '$1'."
            exit 1
        }
    else
        echo "Output directory '$1' already exists. Files will be overwritten."
    fi
}

# Check if all other options are real directories or files
create_valid_directory "$output_directory"
is_valid_file "$fasta_file"
is_valid_fasta "$fasta_file"
is_valid_file "$ref_genome"
is_valid_file "$fq1"
#is_valid_fastq "$fq1"
is_valid_file "$fq2"
#is_valid_fastq "$fq2"

echo "-o: Path to output directory: $output_directory";
echo "-c: Path to consensus fasta: $fasta_file";
echo "-@: Number of threads: $threads";
echo "-1: Path to fastq1: $fq1";
echo "-2: Path to fastq2: $fq2";

# Extract the genome name after the first underscore
#genome_name=$(echo "$fq1" | awk -F'/' '{split($NF,a,"_"); print a[1]}')
#if [[ "$fq1" == */* ]]; then
#    genome_name=$(echo "$fq1" | awk -F'/' '{split($NF,a,"_"); print a[1]}')
#else
#    genome_name=$(echo "$fq1" | awk -F'_' '{print $1}')
#fi
# Print the extracted genome name
echo "Genome name: $genome_name"



#kindof hard coded rn, just doing this to be concise when testing.
echo "-d: Path to reference genome: $ref_genome";

# Check if the file exists
if [ ! -f "$fasta_file" ]; then
  echo "ERROR: $fasta_file not found."
  exit 1
fi

# Initialize an array to store the sequence names
names=()

# Read the FASTA file line by line
while IFS= read -r line; do
  # Check if the line starts with '>'
  if [[ $line == ">"* ]]; then
    # Extract the sequence name by removing the '>' character
    name="${line:1}"
    names+=("$name")
  fi
done < "$fasta_file"

# Print the list of names
#echo "List of names:"
#printf '%s\n' "${names[@]}"

#add reference TEs to TE fasta file
awk '{OFS="\t"; print $1, $2, $3, $4, $5, $6}' ${ref_te_locations} | bedtools getfasta -fi ${ref_genome} -bed - -fo ${output_directory}/extracted_sequences.fasta -nameOnly
cat ${output_directory}/extracted_sequences.fasta ${fasta_file} > ${output_directory}/tefasta.fasta
sed 's/-/_/g' ${output_directory}/tefasta.fasta > ${output_directory}/tefasta2.fasta

echo "Aligning reads to TE seqs"
# Map reads to te seqs, filter unmapped reads (4), sort
bwa-mem2 index ${output_directory}/tefasta2.fasta
bwa-mem2 mem -t $threads ${output_directory}/tefasta2.fasta ${fq1} ${fq2} | samtools view -@ $threads -F 4 -b - | samtools sort -@ $threads -o ${output_directory}/${genome_name}_to_alltes_alignment.bam -
samtools index ${output_directory}/${genome_name}_to_alltes_alignment.bam

# Iterate over each name provided as arguments
for name in ${names[@]}; do
  samtools view -@ ${threads} -H ${output_directory}/${genome_name}_to_alltes_alignment.bam | grep '^@SQ' | cut -f 2 | cut -d ':' -f 2 | sed 's/-/_/g' | grep "${name}$" > ${output_directory}/current_tes.txt
  # Retrieve read names
  samtools view -@ ${threads} -b ${output_directory}/${genome_name}_to_alltes_alignment.bam $(cat "${output_directory}/current_tes.txt" | tr '\n' ' ') > ${output_directory}/${genome_name}_to_${name}_alignment.bam

  samtools view -@ ${threads} -h ${output_directory}/${genome_name}_to_${name}_alignment.bam | awk '{if ($1 ~ /^@/) next; print $1}' > ${output_directory}/readnames_${name}.txt

  # Use read names to retrieve reads from original files
  # These commands ensure each read in the pair is retrieved
  # Change depending on the format of the file
  #sed 's/\(\SRR[0-9]*\.[0-9]*\)/@\1 \n@\1\//' ${output_directory}/readnames_${name}.txt | grep -A 3 -F -f - ${fq1}| sed '/--/d' > ${output_directory}/${name}_${genome_name}_1.fq 
  #sed 's/\(\SRR[0-9]*\.[0-9]*\)/@\1 \n@\1\//' ${output_directory}/readnames_${name}.txt | grep -A 3 -F -f - ${fq2}| sed '/--/d' > ${output_directory}/${name}_${genome_name}_2.fq 
  #sed 's/^/@/' ${output_directory}/readnames_${name}.txt | grep -A 3 -x -F -f - ${fq1} > ${output_directory}/${name}_${genome_name}_1.fq 
  #sed 's/^/@/' ${output_directory}/readnames_${name}.txt | grep -A 3 -x -F -f - ${fq2} > ${output_directory}/${name}_${genome_name}_2.fq 
  seqkit grep -f ${output_directory}/readnames_${name}.txt ${fq1} -o ${output_directory}/${name}_${genome_name}_1.fq
  seqkit grep -f ${output_directory}/readnames_${name}.txt ${fq2} -o ${output_directory}/${name}_${genome_name}_2.fq

  # Remap te mapping reads to the reference genome
  bwa-mem2 mem -t $threads ${ref_genome} ${output_directory}/${name}_${genome_name}_1.fq ${output_directory}/${name}_${genome_name}_2.fq | samtools view -@ $threads -b - | samtools sort -@ $threads -o ${output_directory}/${name}_to_ISO1.bam -
  
  #index te bam file for easy viewing later
  samtools index ${output_directory}/${name}_to_ISO1.bam
done

touch ${output_directory}/completed.txt
echo "Candidate region mapping complete!"