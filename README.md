
# MotifScope
A tool for motif annotation and visualization in tandem repeats.

Motifscope is also available online at [https://motifscope.holstegelab.eu](https://motifscope.holstegelab.eu).

## Installation
- To install with conda
  ```
  cd install/conda
  sh INSTALL.sh
  ```
  Conda will install an environment called motifscope, in which the necessary dependencies are installed. <br>
  The conda environment is activated by executing ```conda activate motifscope``` in the shell. <br>

  See the usage section on how to run MotifScope once the conda environment is activated. <br>

- To install with docker
  ```
  cd install/docker
  sh build.sh
  ```
  Docker will create an image called motifscope, in which the necessary dependencies are installed. <br>

  An example command for running motifscope within this docker image is available in run_docker.sh

  To run it:
  ```
     sh run_docker.sh ../../example/example_sequence.fa ../../example/example_population.txt
  ```
  <br>
  
### Usage
- For running MotifScope on a set of sequences (reads or assemblies):
  ```bash
  motifscope  [-i input.fa] [-mink 2] [-maxk 10] [-o output.prefix]
  ```
    <br>
- To annotate sequences with class labels, one can use the -p option to provide an annotation file. 
  ```bash
  motifscope [-i input.fa] [-mink 2] [-maxk 10] [-p classes.txt] [-o output.prefix]
  ```
  The class information will be shown as a separate color-coded column in the figure.

    - The header of the sequences in ```input.fa``` should start with ```>sample#hap_number#```, for example, for HG002, it could start with ```>HG002#1#``` . <br>
    - The class annotation file ```classes.txt``` should be a tab separated file with the first column being the sample ids and the second column being the sample class. E.g. HG002	EUR <br>
    - The class annotation file can contain a header, which should read ```sample	<class_name>```. The label of the second column <class_name> can be adapted, and will be shown in the figure. When there is no header, the default is 'population'. <br>

    <br>
    
- To disable sequence clustering and the dendrogram (e.g. in case of a single sequence), use the -c option:
  ```bash
  motifscope -c False [-i input.fa] [-mink 2] [-maxk 10] [-o output.prefix]
  ```
<br>

- To run multiple sequence alignment on the compressed representation of the sequence, set ```-msa``` to ```POAMotif``` (aligns complete motifs) or ```POANucleotide``` (aligns nucleotides). <br>

- To guide the algorithm with a set of known motifs, provide the motifs with ```-motifs motifs.txt```. The motif file ```motifs.txt``` should contain the motifs separated with a tab. <br>

- To use random categorical colors for motifs, set ```-e``` to ```random```. To project motifs onto a color scale, set ```-e``` to ```UMAP``` or ```MDS``` for dimension reduction based on motif similarities. <br>

- To characterize motif composition without generating a figure, set ```-figure``` to ```False```. <br>

- To use the reverse complement of the input fasta, set ```-reverse``` to ```True```. <br>

## Output
- The repeat compositions are output in a fasta file. For example,
```bash
>HG002#2#JAHKSD010000034.1:9910981-9913041/rc
G1 A1 G1 C1 A2 G1 A1 C1 T1 C1 T1 G1 T3 C1 A2 AAAAG12 A1 AAAAG1 C1 A1 T1 G1 T2 C1 T1 A3 G1 A1 G1
```
The motifs are separated by spaces. Each string represents a motif, and the following number indicates how many consecutive copies of that motif occur. <br>

- The motif summary per sequence is output in a tab-separated file. The first column is the sequence header, the second column is the motif, the third column is the amount of sequence covered by the motif, and the fourth column is the count of the motif. 
