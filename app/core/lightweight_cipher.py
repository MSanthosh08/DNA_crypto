import hmac, hashlib

def _keystream(key: str, img_hash: bytes, nbytes: int):
    keyb = key.encode()
    out = b''
    counter = 0
    while len(out) < nbytes:
        data = img_hash + counter.to_bytes(8, 'big')
        out += hmac.new(keyb, data, hashlib.sha256).digest()
        counter += 1
    return out[:nbytes]

def xor_encrypt_bytes(plaintext: bytes, key: str, img_hash: bytes) -> bytes:
    ks = _keystream(key, img_hash, len(plaintext))
    return bytes([p ^ k for p,k in zip(plaintext, ks)])

def xor_decrypt_bytes(ciphertext: bytes, key: str, img_hash: bytes) -> bytes:
    return xor_encrypt_bytes(ciphertext, key, img_hash)
