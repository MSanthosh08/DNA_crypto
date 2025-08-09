import os
from app.main import decrypt_image

def test_decrypt_with_wrong_key(tmp_path):
    # Dummy test for integrity check
    encrypted_file = "data/kodak/temp.enc"
    decrypted_img = tmp_path / "decrypted.png"
    wrong_key = "wrongkey"

    try:
        decrypt_image(encrypted_file, wrong_key, str(decrypted_img))
    except Exception:
        assert True  # Expected to fail with wrong key or HMAC check
    else:
        assert False, "Decryption should have failed with wrong key"
