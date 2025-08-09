#!/bin/bash
set -e

IMAGE_NAME="dna-chaos:latest"
CONTAINER_NAME="dna-chaos-server"

# Build image only if it doesn't exist
if [[ "$(docker images -q $IMAGE_NAME 2> /dev/null)" == "" ]]; then
  echo "[+] Building Docker image..."
  docker build -t $IMAGE_NAME .
else
  echo "[+] Docker image $IMAGE_NAME already exists. Skipping build."
fi

echo "[+] Removing old container (if exists)..."
docker rm -f $CONTAINER_NAME || true

echo "[+] Starting container..."
docker run -d \
  --name $CONTAINER_NAME \
  -p 8000:8000 \
  $IMAGE_NAME

echo "[+] Server is running at: http://localhost:8000"