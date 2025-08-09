import numpy as np
from PIL import Image
import hashlib

def image_to_bytes_and_hash(path: str):
    """Return (raw_bytes, shape, sha256_hash) for RGB image."""
    img = Image.open(path).convert('RGB')
    arr = np.array(img, dtype=np.uint8)
    h, w, c = arr.shape
    raw = arr.tobytes()
    hsh = hashlib.sha256(raw).digest()
    return raw, (h, w, c), hsh

def bytes_to_image_and_save(raw_bytes: bytes, shape: tuple, out_path: str):
    """Save raw bytes (RGB) back to image file (shape = (h,w,c))."""
    import numpy as np
    from PIL import Image
    h, w, c = shape
    arr = np.frombuffer(raw_bytes, dtype=np.uint8)
    arr = arr.reshape((h, w, c))
    img = Image.fromarray(arr, 'RGB')
    img.save(out_path)
