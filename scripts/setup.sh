#!/bin/bash

# Create necessary directories for Single - cell RNA - seq project
BASE_DIR = "$(pwd)"

# Create main directories
mkdir -p "${DATA_DIR}/raw"
mkdir -p "${DATA_DIR}/metadata"

mkdir -p "${RESULTS_DIR}/
mkdir -p 

mkdir -p "${RESULT_DIR}/
mkdir -p "${BASE_DIR}/logs"

# Set appropriate permissions
chmod -R 755 "${BASE_DIR}"

# Check and downloads required tools
TOOLS_DIR = "${BASE_DIR}/tools"
mkdir -p "${TOOLS_DIR}"

