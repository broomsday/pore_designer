PROTEINMPNNDIR="/ProteinMPNN"

cd ${PROTEINMPNNDIR}
wget -q -P . https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash ./Mambaforge-Linux-x86_64.sh -b -p ${PROTEINMPNNDIR}/conda
rm Mambaforge-Linux-x86_64.sh
. "${PROTEINMPNNDIR}/conda/etc/profile.d/conda.sh"
export PATH="${PROTEINMPNNDIR}/conda/condabin:${PATH}"
conda create -p $PROTEINMPNNDIR/proteinmpnn-conda python=3.10 -y
conda activate $PROTEINMPNNDIR/proteinmpnn-conda
conda update -n base conda -y
conda install -y pytorch==1.12.1 torchvision==0.13.1 torchaudio==0.12.1 cudatoolkit=11.3 -c pytorch
conda deactivate
