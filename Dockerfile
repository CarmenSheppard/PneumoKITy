# Use Miniconda as the base image
FROM continuumio/miniconda3

# Set up Conda channels for BioConda and Conda-Forge and create environment with specific versions
RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict && \
    conda create -n myenv python=3.9 numpy=1.22.0 pandas=1.3.3 mash sqlalchemy=1.4.22 -y && \
    conda clean -a



# Activate the Conda environment for subsequent commands
SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]

# Set working directory inside the container
WORKDIR /app

# Make sure input/output directories exist
RUN mkdir -p /app/input /app/output

# Add pneumokity files to app
COPY . /app/

# Set entry point
ENTRYPOINT ["conda", "run", "-n", "myenv", "python", "pneumokity.py"]
