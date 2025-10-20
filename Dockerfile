# =====================================
# Stage 1 — Build environment
# =====================================
FROM python:3.11-slim AS builder

ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    wget \
    bcftools \
    libgsl-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b && \
    rm /tmp/miniconda.sh

# Configure Bioconda channels
RUN conda config --system --add channels defaults && \
    conda config --system --add channels bioconda && \
    conda config --system --add channels conda-forge

# Accept tos
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Install rust-bio-tools and varlociraptor from Bioconda
RUN conda install -y rust-bio-tools varlociraptor

# Set working directory
WORKDIR /app

# Copy Poetry configuration first (for caching)
COPY pyproject.toml poetry.lock* ./

# Install Poetry
RUN pip install --no-cache-dir poetry

# Install only dependencies (no packaging step)
RUN poetry install --no-root --no-interaction --no-ansi --without dev

# Copy your project files
COPY preprocess.py main.sh ./

# Make the main script executable
RUN chmod +x main.sh

# Environment variables
ENV PYTHONUNBUFFERED=1

# Default command
ENTRYPOINT ["./main.sh"]
