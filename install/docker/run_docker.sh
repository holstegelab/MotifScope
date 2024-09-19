#check that all the required arguments are provided
if [ "$#" -ne 3 ]; then
	echo "Usage: $0 <input.fa> <population.txt> <output_prefix>"
	exit 1
fi

# Set up variables
#set output dir to current working directory
#OUTPUT_DIR=$(pwd)
#convert $1 to absolute path
FASTA_FILE=$(realpath $1)
#convert $2 to absolute path
POPULATION_FILE=$(realpath $2)
OUTPUT_DIR=$(dirname "$3")       # Extract directory part from $3
OUTPUT_PREFIX=$(basename "$3")   # Extract file prefix from $3

if [[ ! "$OUTPUT_DIR" = /* ]]; then
    OUTPUT_DIR=$(realpath "$OUTPUT_DIR")
fi

# Check if the output directory exists
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Error: Output directory '$OUTPUT_DIR' does not exist."
    exit 1
fi

echo "Input file: ${FASTA_FILE}"
echo "Population file: ${POPULATION_FILE}"
echo "Results will be stored in: ${OUTPUT_DIR}"

# Motifscope case-control analysis command
docker run -it --rm \
        -v ${FASTA_FILE}:/run_input.fa \
        -v ${POPULATION_FILE}:/run_population.txt \
        -v ${OUTPUT_DIR}:/output_dir \
        motifscope \
        -i /run_input.fa \
        -mink 2 \
        -maxk 10 \
        -p /run_population.txt \
        -o /output_dir/${OUTPUT_PREFIX} \
		-msa POAMotif

