# Install base utilities
FROM nvidia/cuda:11.8.0-runtime-ubuntu22.04

RUN apt-get update \
    && apt-get install -y build-essential \
    && apt-get install -y wget \
    && apt-get install -y curl \
    && apt-get install -y git \
    && apt-get install -y unzip \
    && apt-get install -y openssh-server \
    && apt-get install -y python3.10 \
    && apt-get install -y python3.10-dev \ 
    && apt-get install -y python3.10-distutils \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH

# setup AlphaFold
WORKDIR /AlphaFold
RUN wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh
SHELL ["/bin/bash", "--login", "-c"]
RUN bash install_colabbatch_linux.sh
RUN conda install -y cuda -c nvidia/label/cuda-11.8.0

# setup pore-designer
WORKDIR /
RUN git clone https://github.com/broomsday/pore_designer.git
WORKDIR /pore_designer
SHELL ["/bin/bash", "--login", "-c"]
RUN pip install poetry
#RUN poetry config virtualenvs.create false \   # this was causing issues with python 3.12 being used to install
RUN poetry config virtualenvs.create true \
    && poetry install 
ENV PROTEINMPNN="/ProteinMPNN"
ENV PROTEINMPNN_ENV="/ProteinMPNN/proteinmpnn-conda"
ENV ALPHAFOLD_ENV="/AlphaFold/localcolabfold/colabfold-conda"
ENV CONDA_SOURCE="/opt/conda/etc/profile.d/conda.sh"

# setup ProteinMPNN
WORKDIR /
RUN git clone https://github.com/dauparas/ProteinMPNN.git
RUN cp /pore_designer/scripts/install_proteinmpnn.sh /ProteinMPNN/install_proteinmpnn.sh
WORKDIR /ProteinMPNN
SHELL ["/bin/bash", "--login", "-c"]
RUN bash install_proteinmpnn.sh

# ready for design
WORKDIR /workspace
ENTRYPOINT ["tail", "-f", "/dev/null"]