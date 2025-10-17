# =====================================
# Stage 1 — Build environment
# =====================================
FROM python:3.11-slim AS builder

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    wget \
    bcftools \
    tabix \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install Rust (needed for rbt + varlociraptor)
RUN curl https://sh.rustup.rs -sSf | bash -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Install rbt (Rust BioTools) and varlociraptor from crates.io
RUN cargo install rust-bio-tools
RUN cargo install varlociraptor

# Set working directory
WORKDIR /app

# Copy Poetry configuration first (for caching)
COPY pyproject.toml poetry.lock* ./

# Install Poetry
RUN pip install --no-cache-dir poetry

# Install only dependencies (no packaging step)
RUN poetry install --no-root --no-interaction --no-ansi

# Copy your project files
COPY preprocess.py main.sh ./

# Make the main script executable
RUN chmod +x main.sh

# =====================================
# Stage 2 — Runtime image
# =====================================
FROM python:3.11-slim

# Copy system binaries from builder
COPY --from=builder /usr/bin/bcftools /usr/bin/bcftools
COPY --from=builder /usr/bin/tabix /usr/bin/tabix
COPY --from=builder /root/.cargo/bin /usr/local/bin

# Copy environment from builder
COPY --from=builder /usr/local/lib/python3.11 /usr/local/lib/python3.11
COPY --from=builder /usr/local/bin /usr/local/bin

# Copy application
WORKDIR /app
COPY --from=builder /app /app

# Environment variables
ENV PATH="/usr/local/bin:${PATH}"
ENV PYTHONUNBUFFERED=1

# Default command
ENTRYPOINT ["./main.sh"]
