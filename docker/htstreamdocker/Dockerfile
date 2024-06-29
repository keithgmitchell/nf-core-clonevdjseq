# Use an official Python runtime as a parent image
FROM continuumio/miniconda3

# Install wget
RUN apt-get update && apt-get install -y wget

# Download and extract HTStream
RUN wget https://github.com/s4hts/HTStream/releases/download/v1.3.3/HTStream_v1.3.3.tar.gz && \
    tar xvf HTStream_v1.3.3.tar.gz && \
    rm HTStream_v1.3.3.tar.gz

# # Install Bioconda
# RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
#     chmod +x Miniconda3-latest-Linux-x86_64.sh && \
#     ./Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
#     rm Miniconda3-latest-Linux-x86_64.sh

# Add Bioconda channels
RUN /opt/conda/bin/conda config --add channels defaults && \
    /opt/conda/bin/conda config --add channels bioconda && \
    /opt/conda/bin/conda config --add channels conda-forge

# Install HTStream
RUN /opt/conda/bin/conda install htstream

# Define the command to run HTStream
# just enter bash
CMD ["/bin/bash"]