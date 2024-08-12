
echo "installing conda environment for MotifSope"

conda env create -f environment.yml
echo "environment successfully installed"

echo "activating environment"
conda init bash
source ~/.bash_profile
conda activate motifscope
echo "environment activated"

echo "installing pylibsais"
cd ../../
git clone https://github.com/holstegelab/pylibsais.git
cd pylibsais
./setup.py build
./setup.py install
echo "pylibsais successfully installed"

cd ../install/conda
python setup.py install

cd ../../
