# app/core/hmac_utils.py
import hmac
import hashlib

def compute_hmac(data: bytes, key: bytes) -> bytes:
    """Compute HMAC-SHA256 for the given data."""
    return hmac.new(key, data, hashlib.sha256).digest()

def verify_hmac(data: bytes, key: bytes, expected_hmac: bytes) -> bool:
    """Verify the HMAC-SHA256 of the given data."""
    computed = compute_hmac(data, key)
    return hmac.compare_digest(computed, expected_hmac)
