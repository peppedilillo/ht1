FROM ubuntu:18.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Update and install Python 3.6 and pip
RUN apt-get update && apt-get install -y \
    python3.6 \
    python3-pip \
    python3.6-dev \
    build-essential \
    libssl-dev \
    libffi-dev \
    && rm -rf /var/lib/apt/lists/*

# Set Python 3.6 as the default python3
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.6 1

# Upgrade pip to match requirements
RUN python3 -m pip install --upgrade pip==9.0.1

# Set working directory
WORKDIR /workspace

# Copy requirements file
COPY . /workspace/

# Install Python dependencies
RUN pip3 install ".[dev]"

# Default command
CMD ["/bin/bash"]
