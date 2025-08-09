# DNA Crypto

A DNA sequence–inspired encryption/decryption framework that:
- Maps binary image data to synthetic DNA sequences.
- Applies chaotic splicing.
- Encodes DNA to protein sequences.
- Compresses the result.
- Encrypts with XChaCha20-Poly1305 (AEAD) for confidentiality + integrity.
- Adds HMAC for early tampering detection.

## Features
- **CLI** for local encryption/decryption
- **FastAPI microservice** for cloud deployment
- Docker-ready

## Installation

```bash
pip install -r requirements.txt
