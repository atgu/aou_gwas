FROM --platform=linux/amd64 us-central1-docker.pkg.dev/broad-mpg-gnomad/images/vep105:latest

USER root

# Ensure the VEP executable is executable and create a symlink at /vep
RUN chmod +x /opt/vep/src/ensembl-vep/vep && \
    ln -s /opt/vep/src/ensembl-vep/vep /vep

# Install Python 3 and pip if they're not already installed
RUN apt-get update && apt-get install -y python3 python3-pip && \
    rm -rf /var/lib/apt/lists/*

# Copy your vep.py script into the container
COPY vep.py /hail-vep/vep.py