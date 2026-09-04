# load: basic Docker image
FROM continuumio/miniconda3
#python:3.11.4

# Create a "work folder"
WORKDIR /app

# Copy project in this work folder
COPY . /app

# RUN: Create environment with dependencies
RUN conda env create -f /app/requirements/environment.yml

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "VibToolsCondaENV", "/bin/bash", "-c"]

# RUNS:
# Install conda build for package installation
# Install pytest
# Install local package (editable) via conda develop
# Activate corresponding Conda environment in container

RUN conda install -c conda-forge conda-build \
    && conda install -c conda-forge pytest=7.2.0 \
    && conda develop src/ \
    && echo "source activate VibToolsCondaENV" > ~/.bashrc
