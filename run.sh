
---

### **`run.sh`**
```bash
#!/bin/bash
set -e

IMAGE_NAME="dna-chaos:latest"
CONTAINER_NAME="dna-chaos-server"

echo "[+] Building Docker image..."
docker build -t $IMAGE_NAME .

echo "[+] Removing old container (if exists)..."
docker rm -f $CONTAINER_NAME || true

echo "[+] Starting container..."
docker run -d \
  --name $CONTAINER_NAME \
  -p 8000:8000 \
  $IMAGE_NAME

echo "[+] Server is running at: http://localhost:8000"
