#local assemlies were done with otter: https://github.com/holstegelab/otter

BED=rfc1.bed
BAM=HG002.bam #aligned to GRCh38
OUTPUT=HG002.rfc1.assembly.fa

#local assembly
otter assemble -b $BED -o 31 $BAM > $OUTPUT

#header adjustment
PREFIX=$(echo "$OUTPUT" | cut -d '.' -f 1)
awk '/^>/{sub(/>/, "&" ++i "#")} 1' "$OUTPUT" > temp.fa
awk -v prefix="$PREFIX" '/^>/ {sub(/^>/, ">" prefix "#"); print; next} 1' "temp.fa" > $OUTPUT
rm temp.fa
