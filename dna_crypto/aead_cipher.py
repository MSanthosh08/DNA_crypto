# dna_crypto/aead_cipher.py
import nacl.utils
from nacl.bindings import (
    crypto_aead_xchacha20poly1305_ietf_encrypt,
    crypto_aead_xchacha20poly1305_ietf_decrypt
)
import hashlib

def _derive_key(key_str: str) -> bytes:
    # 32-byte key from user key string
    return hashlib.sha256(key_str.encode()).digest()

def aead_encrypt(plaintext: bytes, key_str: str, aad: bytes=b'') -> bytes:
    key = _derive_key(key_str)  # 32 bytes
    nonce = nacl.utils.random(24)
    ct = crypto_aead_xchacha20poly1305_ietf_encrypt(plaintext, aad, nonce, key)
    return nonce + ct  # prepend nonce

def aead_decrypt(blob: bytes, key_str: str, aad: bytes=b'') -> bytes:
    key = _derive_key(key_str)
    nonce = blob[:24]
    ct = blob[24:]
    return crypto_aead_xchacha20poly1305_ietf_decrypt(ct, aad, nonce, key)
