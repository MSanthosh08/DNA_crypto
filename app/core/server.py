# server.py
from fastapi import FastAPI, UploadFile, File, Form, HTTPException
import shutil
from core import utils
import tempfile
import os
from main import encrypt_image, decrypt_image

app = FastAPI(title="DNA-Chaos Encryptor")

@app.post("/encrypt")
async def api_encrypt(file: UploadFile = File(...), key: str = Form(...), no_compress: bool = Form(False), no_splice: bool = Form(False)):
    # save upload to temp file
    tf = tempfile.NamedTemporaryFile(delete=False, suffix=".png")
    with tf as out:
        shutil.copyfileobj(file.file, out)
    out_enc = tf.name + ".enc"
    try:
        encrypt_image(tf.name, key, out_enc, no_compress=no_compress, no_splice=no_splice)
        return {"encrypted_path": out_enc}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    finally:
        # caller can collect out_enc; we keep uploaded tmp for debugging
        pass

@app.post("/decrypt")
async def api_decrypt(file: UploadFile = File(...), key: str = Form(...)):
    tf = tempfile.NamedTemporaryFile(delete=False)
    with tf as out:
        shutil.copyfileobj(file.file, out)
    out_dec = tf.name + ".png"
    try:
        decrypt_image(tf.name, key, out_dec)
        return {"decrypted_path": out_dec}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
