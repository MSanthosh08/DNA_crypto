import os
import pytest
from core.aead_cipher import encrypt_aead, decrypt_aead

def test_aead_round_trip():
    key = os.urandom(32)  # XChaCha20-Poly1305 requires 256-bit key
    plaintext = b"Hello AEAD encryption!"
    aad = b"metadata-header"

    ciphertext, nonce = encrypt_aead(key, plaintext, aad)
    decrypted = decrypt_aead(key, ciphertext, nonce, aad)

    assert decrypted == plaintext, "AEAD decryption failed"
