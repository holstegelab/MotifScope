# Use the official Python image from the Docker Hub
FROM python:3.10-slim

# Set the working directory
WORKDIR /app

# Copy the requirements file into the container
COPY requirements.txt .

# Install the dependencies
RUN apt-get update && \
        apt-get install -y --no-install-recommends \
        curl \
        gcc \
        build-essential \
        libffi-dev \
        libssl-dev \
        git \
        zlib1g-dev \
        libpng-dev \
        && \
        pip install --no-cache-dir -r requirements.txt && \
        rm -rf /var/lib/apt/lists/*

# Install TRF
# !! You can remove these lines
RUN git clone https://github.com/holstegelab/pylibsais.git /opt/pylibsais
# Set up build directory
WORKDIR /opt/pylibsais
# Configure and build TRF
RUN ./setup.py build
# Make it system-wide available
RUN ./setup.py install

# Clone the Motifscope GitHub repository
# !! You can remove this line, and change it to MotifScope
RUN git clone https://github.com/holstegelab/MotifScope.git /app/MotifScope

WORKDIR /app/MotifScope

RUN python install/docker/setup.py install

# Change working directory to the cloned repository
# Modify this line -- changing treat to motifscope should work
WORKDIR /app/MotifScope

# Clean image by removing unnecessary programs
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        libssl-dev \
        libffi-dev \
        && apt-get remove -y build-essential \
        && apt-get autoremove -y \
        && rm -rf /var/lib/apt/lists/*

# Copy the rest of your application files

# Define the default command to run your application
# Modify the mofflowind depending on how you run motifscope
ENTRYPOINT ["python", "motifscope"]
