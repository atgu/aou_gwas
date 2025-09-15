# 1. Force the amd64 SAIGE 1.4.4 base image
FROM --platform=linux/amd64 wzhou88/saige:1.4.4

USER root
RUN apt-get update && \
    apt-get install -y --no-install-recommends wget tar libgomp1 && \
    rm -rf /var/lib/apt/lists/*

# 2. Install KING
RUN wget -O /tmp/king.tar.gz https://www.kingrelatedness.com/Linux-king.tar.gz && \
    tar -xzvf /tmp/king.tar.gz -C /usr/local/bin && \
    chmod +x /usr/local/bin/king && \
    rm /tmp/king.tar.gz

# 3. Verify without breaking the build
#    king --help will exit non-zero unless you give it data,
#    so we append "|| true" to swallow that exit code.
RUN king --help || true && \
    Rscript -e "library(SAIGE); message('SAIGE version: ', packageVersion('SAIGE'))"

# Use the exact Artifact Registry path you specified
# IMAGE="us-central1-docker.pkg.dev/aou-neale-gwas/saige/saige-king:0.0.1"

# Authenticate Docker to Artifact Registry
# gcloud auth configure-docker us-central1-docker.pkg.dev --quiet

# Build the image (assumes Dockerfile is in the current directory)
# docker build -t "${IMAGE}" .


# Push the image
# docker push "${IMAGE}"