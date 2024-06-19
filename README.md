# MotifScope
A tool for motif annotation and visualization in tandem repeats.

## Installation
- Clone the repository with the following command:
  ```bash
  git clone https://github.com/holstegelab/MotifScope.git
  ```
- Install the conda environment:
  ```bash
  conda env create -f environment.yml
  ```
- Activate the envrionment:
  ```bash
  conda activate motifscope
  ```
- Install pylibsais (https://github.com/holstegelab/pylibsais)
  ```bash
  sh INSTALL.sh
  ```
- Install MotifScope
  ```bash
  python setup.py install
  ```
### Usage
- For running MotifScope on a set of sequences (reads or assemblies) without population data:
```bash
motifscope --sequence-type reads [-i input.fa] [-mink 2] [-maxk 10] [-t title] 
```
- For running MotifScope on sequences with population data:
```bash
motifscope --sequence-type assembly [-i input.fa] [-mink 2] [-maxk 10] [-t title] [-p population.txt]
```
- For running MotifScope on a single sequence:
```bash
motifscope --sequence-type single [-i input.fa] [-mink 2] [-maxk 10] [-t title] 
```
The description of the sequences in ```input.fa``` should start with ```>sample_name#```. <br>

The population file ```population.txt``` should be a tab separated file with the first column being the sample names and the second column being the population. 

To run multiple sequence alignment on the compressed representation of the sequence, set ```-msa``` to ```True```. <br>

To run the algorithm with a set of known motifs, set ```-m``` to ```True``` and provide the motifs with ```-motifs motifs.txt```. The motif file ```motifs.txt``` should contain the motifs separated with tab. <br>

To use random categorical colors for motifs, set ```-e``` to ```random```. To project motifs onto a continuous color scale, set ```-e``` to ```UMAP``` or ```MDS``` for dimension reduction based on motif similarities. <br>

The repeat compositions are output in a fasta file. For example,
```bash
>HG002#2#JAHKSD010000034.1:9910981-9913041/rc
G1 A1 G1 C1 A2 G1 A1 C1 T1 C1 T1 G1 T3 C1 A2 AAAAG12 A1 AAAAG1 C1 A1 T1 G1 T2 C1 T1 A3 G1 A1 G1
```
The motifs are separated by spaces. Each string represents a motif, and the following number indicates how many consecutive copies of that motif occur.
