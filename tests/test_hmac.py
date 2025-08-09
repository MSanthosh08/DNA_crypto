import os
import pytest
from core.hmac_utils import generate_hmac, verify_hmac

def test_hmac_validation():
    key = os.urandom(32)
    data = b"Important data payload"

    hmac_value = generate_hmac(key, data)
    assert verify_hmac(key, data, hmac_value) is True

def test_hmac_tampering():
    key = os.urandom(32)
    data = b"Important data payload"

    hmac_value = generate_hmac(key, data)
    tampered_data = b"Important data payloAd"  # changed one letter

    assert verify_hmac(key, tampered_data, hmac_value) is False
