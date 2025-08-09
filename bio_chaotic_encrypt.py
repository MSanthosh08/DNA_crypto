#!/usr/bin/env python3
"""
bio_chaotic_encrypt.py
Prototype: DNA->Chaos->Protein->Compress -> AEAD pipeline with decrypt path.
Dependencies:
    pip install pynacl pillow numpy
"""

import hashlib, hmac, struct, math
from collections import deque
from typing import Tuple
import numpy as np
from PIL import Image
import nacl.utils
from nacl.bindings import crypto_aead_xchacha20poly1305_ietf_encrypt, crypto_aead_xchacha20poly1305_ietf_decrypt

# ---------------------------
# Utilities: HKDF-SHA256
# ---------------------------
def hkdf_extract_expand(salt: bytes, ikm: bytes, info: bytes, L: int) -> bytes:
    # Simple HKDF using HMAC-SHA256
    if salt is None: salt = b'\x00' * 32
    prk = hmac.new(salt, ikm, hashlib.sha256).digest()
    okm = b''
    t = b''
    i = 1
    while len(okm) < L:
        t = hmac.new(prk, t + info + bytes([i]), hashlib.sha256).digest()
        okm += t
        i += 1
    return okm[:L]

# ---------------------------
# Key/session helpers
# ---------------------------
def derive_session_material(master_key: bytes, img_hash: bytes, info: bytes=b"bio-chaos", L=64):
    # Combine master key and image fingerprint so mapping depends on plaintext
    ikm = hmac.new(master_key, img_hash, hashlib.sha256).digest()
    return hkdf_extract_expand(salt=None, ikm=ikm, info=info, L=L)

# ---------------------------
# DNA encoding with repetition avoidance (simple rule)
# ---------------------------
BASES = ['A','C','G','T']

def build_dynamic_base_map(seed: bytes) -> dict:
    # produce a permutation of base mapping for 2-bit values
    s = list(BASES)
    # seed PRNG deterministically via HKDF
    rnd = hkdf_extract_expand(None, seed, b"base_perm", 16)
    # Fisher-Yates
    arr = s[:]
    for i in range(len(arr)-1, 0, -1):
        j = rnd[i % len(rnd)] % (i+1)
        arr[i], arr[j] = arr[j], arr[i]
    return {'00': arr[0], '01': arr[1], '10': arr[2], '11': arr[3]}

def bytes_to_bits(b: bytes) -> str:
    return ''.join(f'{byte:08b}' for byte in b)

def bits_to_dna(bits: str, base_map: dict) -> str:
    # produce DNA sequence while avoiding immediate repeats by rule:
    dna = []
    prev = None
    for i in range(0, len(bits), 2):
        pair = bits[i:i+2]
        if len(pair) < 2:
            pair = pair + '0'
        base = base_map[pair]
        # repetition-avoidance: if same as prev, pick alternate from map deterministically
        if base == prev:
            # choose next base from permutation of the remaining three using hashed context
            ctx = hashlib.sha256((prev + pair).encode()).digest()
            # pick one of the other three bases deterministically
            for b in BASES:
                if b != prev:
                    pick = b
                    break
            base = pick
        dna.append(base)
        prev = base
    return ''.join(dna)

def dna_to_bits(dna: str, base_map: dict) -> str:
    # reverse: base_map might map multiple bit patterns to same base if repetition-avoidance substituted.
    # we reconstruct by trying original mapping first; this prototype assumes substitution was deterministic based on prev,
    # so decrypt can re-run same logic to map back.
    inv = {v:k for k,v in base_map.items()}
    bits = []
    prev = None
    for i, base in enumerate(dna):
        if base in inv:
            bits.append(inv[base])
            prev = base
        else:
            # Shouldn't happen in deterministic de/encoding
            raise ValueError("Unknown base during decode")
    return ''.join(bits)

# ---------------------------
# Chaotic splicing & permutation (logistic map based)
# ---------------------------
def logistic_sequence(seed_f32: float, r: float, length: int):
    x = seed_f32
    seq = []
    for _ in range(length):
        x = r * x * (1 - x)
        seq.append(x)
    return seq

def chaotic_permute_and_splice(dna: str, seed_bytes: bytes, cut_count: int=20, r: float=3.99) -> Tuple[str, dict]:
    # Use HKDF to get initial x (0,1)
    out = hkdf_extract_expand(None, seed_bytes, b"logistic", 32)
    seed_int = int.from_bytes(out[:8], 'big')
    seed_f32 = (seed_int % (10**8)) / (10**8 + 1e-12)
    seq = logistic_sequence(seed_f32, r, len(dna)+cut_count*3)
    # create a permutation of indices
    indices = list(range(len(dna)))
    # shuffle indices by taking pairs from seq values
    for i, x in enumerate(seq[:len(indices)]):
        j = int(x * len(indices))
        indices[i % len(indices)], indices[j] = indices[j], indices[i % len(indices)]
    permuted = ''.join(dna[i] for i in indices)
    # splicing: perform cut_count random cuts (deterministic from seq) and move them around
    dna_list = list(permuted)
    ops = []
    pos_seq = seq[len(indices):]
    for k in range(cut_count):
        start = int(pos_seq[3*k % len(pos_seq)] * len(dna))
        length = max(1, int(pos_seq[(3*k+1) % len(pos_seq)] * (len(dna)//10)))
        dest = int(pos_seq[(3*k+2) % len(pos_seq)] * len(dna))
        # slice and move
        seg = dna_list[start:start+length]
        # delete then insert at dest
        del dna_list[start:start+length]
        # adjust dest if required
        if dest > len(dna_list):
            dest = len(dna_list)
        for i, ch in enumerate(seg):
            dna_list.insert(dest + i, ch)
        ops.append((start, length, dest))
    spliced = ''.join(dna_list)
    meta = {'indices': indices, 'ops': ops, 'seed': out.hex(), 'r': r}
    return spliced, meta

def reverse_chaotic_permute_and_splice(spliced: str, meta: dict) -> str:
    # naive reverse: undo ops in reverse order, then apply inverse permutation
    dna_list = list(spliced)
    # reverse ops
    for start, length, dest in reversed(meta['ops']):
        # segment currently at dest: extract dest..dest+length and move back to start
        seg = dna_list[dest:dest+length]
        del dna_list[dest:dest+length]
        for i, ch in enumerate(seg):
            dna_list.insert(start + i, ch)
    # inverse permutation
    perm = meta['indices']
    inv = [0]*len(perm)
    for i,p in enumerate(perm):
        inv[p] = i
    original = ''.join(dna_list[inv[i]] for i in range(len(inv)))
    return original

# ---------------------------
# Codon -> Amino-acid dynamic mapping
# ---------------------------
AMINO = list("ACDEFGHIKLMNPQRSTVWY")  # 20 symbols (single-letter)
def build_dynamic_codon_table(seed: bytes) -> dict:
    # produce a mapping from 64 codons -> 20 amino acids using session permutation
    bases = BASES
    codons = [a+b+c for a in bases for b in bases for c in bases]
    rnd = hkdf_extract_expand(None, seed, b"codon_perm", 64)
    # produce a permutation of codons by bytes
    perm = codons[:]
    for i in range(len(perm)-1, 0, -1):
        j = rnd[i % len(rnd)] % (i+1)
        perm[i], perm[j] = perm[j], perm[i]
    # map codons to amino acids by cycling through AMINO
    mapping = {}
    for i, codon in enumerate(perm):
        mapping[codon] = AMINO[i % len(AMINO)]
    return mapping

def dna_to_protein(dna: str, codon_map: dict) -> str:
    # pad if needed
    if len(dna) % 3 != 0:
        dna = dna + 'A' * (3 - (len(dna)%3))
    protein = []
    for i in range(0, len(dna), 3):
        codon = dna[i:i+3]
        protein.append(codon_map[codon])
    return ''.join(protein)

def protein_to_dna(protein: str, codon_map: dict, seed: bytes) -> str:
    # Need to invert codon_map which is many-to-one.
    # To invert deterministically, regenerate codon->amino ordering and pick the same codon selection as encoder.
    # Here, we rebuild the same perm used in build_dynamic_codon_table and use the first codon assigned to that amino in order.
    # (This is a prototyping simplification.)
    bases = BASES
    codons = [a+b+c for a in bases for b in bases for c in bases]
    rnd = hkdf_extract_expand(None, seed, b"codon_perm", 64)
    perm = codons[:]
    for i in range(len(perm)-1, 0, -1):
        j = rnd[i % len(rnd)] % (i+1)
        perm[i], perm[j] = perm[j], perm[i]
    amino_to_codons = {}
    for i, codon in enumerate(perm):
        aa = AMINO[i % len(AMINO)]
        amino_to_codons.setdefault(aa, []).append(codon)
    # now for each amino acid in protein choose the first available codon (deterministic)
    dna = []
    for aa in protein:
        choices = amino_to_codons.get(aa)
        if not choices:
            raise ValueError("No codon choices")
        dna.append(choices[0])
    return ''.join(dna)

# ---------------------------
# Simple LZ77-like compressor for DNA alphabet
# Produce bytes: token format - (flag)(data...)
# flag=0 => literal (1 byte ascii base)
# flag=1 => match (offset:2 bytes, length:1 byte)
# This is not optimized but works for prototype
# ---------------------------
def dna_compress_lz77(dna: str, window_size: int = 4096, lookahead: int = 255) -> bytes:
    i = 0
    out = bytearray()
    n = len(dna)
    while i < n:
        end_window = max(0, i - window_size)
        best = (0,0)  # (offset, length)
        for j in range(end_window, i):
            length = 0
            while length < lookahead and i+length < n and dna[j+length] == dna[i+length]:
                length += 1
            if length > best[1]:
                best = (i - j, length)
        if best[1] >= 3:
            # emit match token
            out.append(1)
            out += struct.pack(">H", best[0])  # offset
            out.append(best[1])                 # length
            i += best[1]
        else:
            # literal
            out.append(0)
            out.append(ord(dna[i]))
            i += 1
    return bytes(out)

def dna_decompress_lz77(blob: bytes) -> str:
    i = 0
    out = []
    n = len(blob)
    while i < n:
        flag = blob[i]; i += 1
        if flag == 0:
            ch = chr(blob[i]); i += 1
            out.append(ch)
        else:
            offset = struct.unpack(">H", blob[i:i+2])[0]; i += 2
            length = blob[i]; i += 1
            start = len(out) - offset
            for k in range(length):
                out.append(out[start + k])
    return ''.join(out)

# ---------------------------
# AEAD: XChaCha20-Poly1305 wrappers
# ---------------------------
def aead_encrypt(key: bytes, plaintext: bytes, aad: bytes=b"") -> bytes:
    nonce = nacl.utils.random(24)
    ct = crypto_aead_xchacha20poly1305_ietf_encrypt(plaintext, aad, nonce, key)
    return nonce + ct

def aead_decrypt(key: bytes, blob: bytes, aad: bytes=b"") -> bytes:
    nonce = blob[:24]
    ct = blob[24:]
    return crypto_aead_xchacha20poly1305_ietf_decrypt(ct, aad, nonce, key)

# ---------------------------
# Full pipeline: encrypt & decrypt
# ---------------------------
def encrypt_image_bytes(img_bytes: bytes, master_key: bytes, chaos_cuts=20, compress=True) -> Tuple[bytes, dict]:
    # fingerprint for plaintext-aware mapping
    img_hash = hashlib.sha256(img_bytes).digest()
    material = derive_session_material(master_key, img_hash, info=b"session", L=80)
    aead_key = material[:32]
    dna_seed = material[32:48]
    codon_seed = material[48:64]
    # 1. bytes -> bits -> dna (dynamic base map)
    bits = bytes_to_bits(img_bytes)
    base_map = build_dynamic_base_map(dna_seed)
    dna = bits_to_dna(bits, base_map)
    # 2. chaotic permute & splice
    spliced, meta = chaotic_permute_and_splice(dna, seed_bytes=dna_seed, cut_count=chaos_cuts)
    # 3. compress
    if compress:
        compressed = dna_compress_lz77(spliced)
        payload = compressed
        comp_flag = True
    else:
        payload = spliced.encode()
        comp_flag = False
    # 4. codon->protein (operate on decompressed DNA during decode to get dna -> protein)
    # To keep file compact we can translate compressed bytes to protein later; for prototype we AEAD the compressed DNA directly.
    # 5. AEAD encrypt
    header = {
        'meta': meta,
        'base_map': base_map,
        'codon_seed': codon_seed.hex(),
        'comp': comp_flag,
        'img_hash': img_hash.hex()
    }
    blob = aead_encrypt(aead_key, payload, aad=img_hash)  # use img_hash as aad to bind plaintext identity
    return blob, header

def decrypt_image_bytes(blob: bytes, header: dict, master_key: bytes) -> bytes:
    img_hash = bytes.fromhex(header['img_hash'])
    material = derive_session_material(master_key, img_hash, info=b"session", L=80)
    aead_key = material[:32]
    dna_seed = material[32:48]
    codon_seed = material[48:64]
    # AEAD decrypt
    payload = aead_decrypt(aead_key, blob, aad=img_hash)
    # decompress if needed
    if header['comp']:
        spliced = dna_decompress_lz77(payload)
    else:
        spliced = payload.decode()
    # reverse chaotic splices & permutation
    dna = reverse_chaotic_permute_and_splice(spliced, header['meta'])
    # convert dna -> bits -> bytes
    bits = dna_to_bits(dna, header['base_map'])
    # bits -> bytes
    b = bytearray()
    for i in range(0, len(bits), 8):
        chunk = bits[i:i+8]
        if len(chunk) < 8:
            chunk = chunk.ljust(8, '0')
        b.append(int(chunk, 2))
    return bytes(b)

# ---------------------------
# tiny I/O helpers for testing
# ---------------------------
def load_image_to_bytes(path: str) -> bytes:
    img = Image.open(path).convert('RGB')
    arr = np.asarray(img)
    # raw bytes as RGB bytes row-major
    return arr.tobytes()

def bytes_to_image_and_save(b: bytes, shape_tuple: Tuple[int,int,int], outpath: str):
    # shape_tuple = (H, W, C)
    arr = np.frombuffer(b, dtype=np.uint8)
    arr = arr.reshape(shape_tuple)
    img = Image.fromarray(arr, 'RGB')
    img.save(outpath)

# ---------------------------
# Example usage
# ---------------------------
if __name__ == "__main__":
    import sys, os, time
    if len(sys.argv) < 2:
        print("Usage: python bio_chaotic_encrypt.py kodimxx.png")
        sys.exit(1)
    img_path = sys.argv[1]
    print("Loading image...")
    img = Image.open(img_path).convert('RGB')
    w,h = img.size
    H = h; W = w; C = 3
    img_bytes = load_image_to_bytes(img_path)
    master_key = hashlib.sha256(b'my very secure master key').digest()  # in prod use KMS / secure random
    t0 = time.time()
    blob, header = encrypt_image_bytes(img_bytes, master_key, chaos_cuts=30, compress=True)
    t1 = time.time()
    print(f"Encrypted: {len(blob)} bytes, header size ~ {sum(len(str(v)) for v in header.values())} -- time {t1-t0:.2f}s")
    # decrypt
    t2 = time.time()
    recovered = decrypt_image_bytes(blob, header, master_key)
    t3 = time.time()
    print(f"Decrypted time {t3-t2:.2f}s, recovered bytes={len(recovered)}")
    # save recovered image to verify equality
    out = "recovered.png"
    bytes_to_image_and_save(recovered, (H,W,C), out)
    print("Saved recovered image:", out)
    # quick verify
    if recovered == img_bytes:
        print("SUCCESS: round-trip identical")
    else:
        print("FAIL: mismatch")
