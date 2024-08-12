# Set up variables
INPUT_FA=/root/example/example.fa
OUTPUT_DIR=/root/example/
POPULATION=/root/example/population.txt

# TREAT case-control analysis command
docker run -it --rm \
        -v ${INPUT_FA}:/run_input.fa \
        -v ${OUTPUT_DIR}:/output_dir \
        -v ${POPULATION}:/run_population.txt \
        motifscope \
        --sequence-type reads \
        -i /run_input.fa \
        -mink 2 \
        -maxk 10 \
        -p /run_population.txt \
        -o /output_dir/example \
	-msa POAMotif

