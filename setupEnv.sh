#SOURCE THIS BEFORE RUNNING PHENO CODE
export PYTHONPATH=/data/atlas/atlasdata2/DiHiggsSharedSoftware/YODA/lib/python2.7/site-packages:$PYTHONPATH
export PATH=/afs/cern.ch/sw/lcg/external/Boost/1.48.0_python2.7/x86_64-slc6-gcc46-opt//include/boost-1_48:$PATH

export PATH=/data/atlas/atlasdata2/DiHiggsSharedSoftware/LHAPDF616/local/bin:$PATH
export PATH=/data/atlas/atlasdata2/DiHiggsSharedSoftware/SHERPA/bin:$PATH
export PATH=/data/atlas/atlasdata2/DiHiggsSharedSoftware/Pythia/pythia8201/bin:$PATH
export PATH=/data/atlas/atlasdata2/DiHiggsSharedSoftware/YODA/YODA-1.5.0/bin/:$PATH
export PATH=/data/atlas/atlasdata2/DiHiggsSharedSoftware/FastJet/fastjet-install/bin:$PATH
export LD_LIBRARY_PATH=/data/atlas/atlasdata2/DiHiggsSharedSoftware/YODA/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/data/atlas/atlasdata2/DiHiggsSharedSoftware/Pythia/pythia8201/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/afs/cern.ch/sw/lcg/external/HepMC/2.06.08/x86_64-slc6-gcc47-opt/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/data/atlas/atlasdata2/DiHiggsSharedSoftware/LHAPDF616/local/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/data/atlas/atlasdata2/DiHiggsSharedSoftware/SHERPA/lib:$LD_LIBRARY_PATH

module load gcc/4.9.3
module load python/2.7.13__gcc49
source /data/atlas/atlasdata2/DiHiggsSharedSoftware/HH4B-Pheno-venv/bin/activate

