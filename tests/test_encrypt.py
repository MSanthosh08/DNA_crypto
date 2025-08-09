import os
import filecmp
from app.core import utils
from app.main import encrypt_image, decrypt_image

def test_encrypt_decrypt_cycle():
    input_img = "data/kodak/kodim01.png"
    encrypted_file = "data/kodak/temp.enc"
    decrypted_img = "data/kodak/decrypted.png"
    key = "testkey123"

    # Encrypt
    encrypt_image(input_img, key, encrypted_file)

    assert os.path.exists(encrypted_file), "Encrypted file was not created"

    # Decrypt
    decrypt_image(encrypted_file, key, decrypted_img)

    assert os.path.exists(decrypted_img), "Decrypted image was not created"

    # Compare byte-by-byte
    assert filecmp.cmp(input_img, decrypted_img, shallow=False), \
        "Decrypted image does not match original"

    # Cleanup
    os.remove(encrypted_file)
    os.remove(decrypted_img)
