
# MotifScope
A tool for motif annotation and visualization in tandem repeats.

## Installation
- To install with conda
  ```
  cd install/conda
  sh INSTALL.sh
  ```
- To install with docker
  ```
  cd install/docker
  docker build --network host -t motifscope .
  ```
  The example command is provided in run_docker.sh
  <br>
### Usage
- For running MotifScope on a set of sequences (reads or assemblies) without population data:
  ```bash
  motifscope --sequence-type reads [-i input.fa] [-mink 2] [-maxk 10] [-o output.prefix]
  ```
    - The header of the sequences in ```input.fa``` should start with ```>sample_name#read_number#```, for example, ```>HG002#1#```. <br>
    <br>
- For running MotifScope on local assemblies with population data:
  ```bash
  motifscope --sequence-type assembly [-i input.fa] [-mink 2] [-maxk 10] [-p population.txt] [-o output.prefix]
  ```
    - The header of the sequences in ```input.fa``` should start with ```>sample_name#hap_number#```, for example, for HG002, it should either start with ```>HG002#1#``` or ```>HG002#2#```. <br>
    - The population file ```population.txt``` should be a tab separated file with the first column being the sample names and the second column being the population. <br>
    <br>
    
- For running MotifScope on a single sequence:
  ```bash
  motifscope --sequence-type single [-i input.fa] [-mink 2] [-maxk 10] [-o output.prefix]
  ```
<br>

- To run multiple sequence alignment on the compressed representation of the sequence, set ```-msa``` to ```POAMotif``` (recommended, aligns based on motifs) or ```POANucleotide``` (aligns based on nucleotides). <br>

- To run the algorithm with a set of known motifs, set ```-g``` to ```True``` and provide the motifs with ```-motifs motifs.txt```. The motif file ```motifs.txt``` should contain the motifs separated with tab. <br>

- To use random categorical colors for motifs, set ```-e``` to ```random```. To project motifs onto a color scale, set ```-e``` to ```UMAP``` or ```MDS``` for dimension reduction based on motif similarities. <br>

- To characterize motif composition without generating a figure, set ```-figure``` to ```False```. <br>

- To use the reverse complement of input fasta, set ```-reverse``` to ```True```. <br>

## Output
- The repeat compositions are output in a fasta file. For example,
```bash
>HG002#2#JAHKSD010000034.1:9910981-9913041/rc
G1 A1 G1 C1 A2 G1 A1 C1 T1 C1 T1 G1 T3 C1 A2 AAAAG12 A1 AAAAG1 C1 A1 T1 G1 T2 C1 T1 A3 G1 A1 G1
```
The motifs are separated by spaces. Each string represents a motif, and the following number indicates how many consecutive copies of that motif occur. <br>

- The motif summary per sequence is output in a tab-separated file. The first column is the sequence header, the second column is the motif, the third column is the amount of sequence covered by the motif, and the fourth column is the count of the motif. 
