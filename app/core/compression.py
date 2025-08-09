import zlib

def compress_bytes(b: bytes) -> bytes:
    return zlib.compress(b, level=6)

def decompress_bytes(b: bytes) -> bytes:
    return zlib.decompress(b)
