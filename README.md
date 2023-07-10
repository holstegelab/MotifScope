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

### Usage
- For running MotifScope on a set of sequences (reads or assemblies) without population data:
```bash
python motifscope.py --sequence-type reads [-i input.fa] [-mink 2] [-maxk 10] [-t title] 
```
- For running MotifScope on sequences with population data:
```bash
python motifscope.py --sequence-type assembly [-i input.fa] [-mink 2] [-maxk 10] [-t title] [-p population.txt]
```
The description of the sequences in ```input.fa``` should start with ```>sample_name#```. <br>
The population file ```population.txt``` should be a tab separated file with the first column being the sample names and the second column being the population. 

To run multiple sequence alignment on the compressed representation of the sequence, set ```-msa``` to ```True```. <br>
To run the algorithm with a set of known motifs, set ```-m``` to ```True``` and provide the motifs with ```-motifs motifs.txt```. The motif file ```motifs.txt``` should contain the motifs separated with tab. 
