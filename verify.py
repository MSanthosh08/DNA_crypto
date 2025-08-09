import hashlib, numpy as np
from PIL import Image
orig = np.array(Image.open("data/kodak/kodim01.png").convert("RGB"), dtype=np.uint8).tobytes()
dec  = np.array(Image.open("decrypted.png").convert("RGB"), dtype=np.uint8).tobytes()
print("orig:", hashlib.sha256(orig).hexdigest())
print("dec :", hashlib.sha256(dec).hexdigest())
print("equal?", orig==dec)