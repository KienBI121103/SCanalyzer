#!/bin/bash

# Create necessary directories for Single - cell RNA - seq project
BASE_DIR = "$(pwd)"

# Create main directories
mkdir -p "${DATA_DIR}/raw"
mkdir -p "${DATA_DIR}/metadata"
mkdir -p "${RESULTS_DIR}/