# Install base image and update OS utilities
FROM nvidia/cuda:11.8.0-devel-ubuntu22.04

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
    && apt-get install -y python3-pip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# pull in the installer utility scripts
WORKDIR /
RUN git clone https://github.com/broomsday/external_tools.git

# setup AlphaFold
WORKDIR /AlphaFold
#RUN git clone https://github.com/YoshitakaMo/localcolabfold.git ./
RUN cp /external_tools/install/install_alphafold* /AlphaFold/
SHELL ["/bin/bash", "--login", "-c"]
#RUN bash install_colabbatch_linux.sh
#RUN cp /external_tools/install/fix_alphafold.sh /AlphaFold/fix_alphafold.sh
#RUN bash fix_alphafold.sh
RUN bash install_alphafold_step1.sh
RUN bash install_alphafold_step2.sh
RUN bash install_alphafold_step3.sh
RUN bash install_alphafold_step4.sh

# setup ProteinMPNN
WORKDIR /ProteinMPNN
RUN git clone https://github.com/dauparas/ProteinMPNN.git .
RUN cp /external_tools/install/install_proteinmpnn* /ProteinMPNN/
SHELL ["/bin/bash", "--login", "-c"]
RUN bash install_proteinmpnn_step1.sh
RUN bash install_proteinmpnn_step2.sh

# force redoing the pore-designer portion each build
ARG CACHEBUST
RUN echo "$CACHEBUST"

# setup environment variables
ENV PROTEINMPNN="/ProteinMPNN"
ENV PROTEINMPNN_ENV="/ProteinMPNN/proteinmpnn-conda"
ENV ALPHAFOLD_ENV="/AlphaFold/localcolabfold/colabfold-conda"
ENV CONDA_SOURCE="/AlphaFold/localcolabfold/conda/etc/profile.d/conda.sh"
# TODO: should just source the proteinmpnn and alphafold condas when needed
# TODO: above requires change to the code in proteinmpnn.py and alphafold.py when building the launch scripts

# setup pore-designer
WORKDIR /pore_designer
RUN git clone https://github.com/broomsday/pore_designer.git .
SHELL ["/bin/bash", "--login", "-c"]
RUN pip install poetry
RUN poetry config virtualenvs.create true \
    && poetry install 

# get ready for design
WORKDIR /workspace
ENTRYPOINT ["tail", "-f", "/dev/null"]