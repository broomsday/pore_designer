# Install base utilities
FROM nvidia/cuda:11.7.1-runtime-ubuntu22.04

RUN apt-get update \
    && apt-get install -y build-essential \
    && apt-get install -y wget \
    && apt-get install -y curl \
    && apt-get install -y git \
    && apt-get install -y unzip \
    && apt-get install -y openssh-server \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

# setup ProteinMPNN
RUN git clone https://github.com/dauparas/ProteinMPNN.git
WORKDIR /ProteinMPNN
SHELL ["/bin/bash", "--login", "-c"]
RUN conda create --name mlfold
RUN source activate mlfold
RUN conda install -y pytorch torchvision torchaudio cudatoolkit=11.3 -c pytorch
RUN source deactivate

# setup AlphaFold colab (above needs to be working though)
WORKDIR /AlphaFold
RUN wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh
RUN bash install_colabbatch_linux.sh
RUN source activate /AlphaFold/localcolabfold/colabfold-conda
RUN conda install -y cuda -c nvidia/label/cuda-11.6.0
RUN source deactivate

# setup pore-designer
WORKDIR /
RUN pip install poetry

RUN git clone https://github.com/broomsday/pore_designer.git
WORKDIR /pore_designer

SHELL ["/bin/bash", "--login", "-c"]
RUN poetry config virtualenvs.create false \
    && poetry install 

ENV PROTEINMPNN="/ProteinMPNN"
ENV ALPHAFOLD_ENV="/AlphaFold/localcolabfold/colabfold-conda"
ENV CONDA_SOURCE="/opt/conda/etc/profile.d/conda.sh"

# ready for design
WORKDIR /workspace
ENTRYPOINT ["tail", "-f", "/dev/null"]