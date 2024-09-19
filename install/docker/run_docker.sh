# Set up variables
OUTPUT_DIR=cwd

#check that all the required arguments are provided
if [ "$#" -ne 3 ]; then
	echo "Usage: $0 <input.fa> <population.txt> <output_prefix>"
	exit 1
fi

# Motifscope case-control analysis command
docker run -it --rm \
        -v $1:/run_input.fa \
        -v ${OUTPUT_DIR}:/output_dir \
        -v $2:/run_population.txt \
        motifscope \
        -i /run_input.fa \
        -mink 2 \
        -maxk 10 \
        -p /run_population.txt \
        -o /output_dir/$3 \
	-msa POAMotif

