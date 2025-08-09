#!/usr/bin/env python3
import argparse
import struct
from dna_crypto import utils, dna_mapping, chaotic_splicing, protein_encoding, compression, lightweight_cipher
import hashlib

MAGIC = b'DNAC'   # 4 bytes magic
VERSION = 1       # 4 bytes

def encrypt_image(input_path: str, key: str, output_path: str):
    raw_bytes, shape, img_hash = utils.image_to_bytes_and_hash(input_path)

    # 1. bytes -> DNA (bases string)
    dna = dna_mapping.binary_to_dna(raw_bytes, key)

    # 2. chaotic splice (permutation)
    spliced = chaotic_splicing.splice_sequence(dna, key, img_hash)

    # 3. codon (3 bases) -> 64-symbol bytes (session-permuted)
    codon_bytes, codon_meta = protein_encoding.dna_to_codon_bytes(spliced, key, img_hash)

    # 4. compress (lossless)
    compressed = compression.compress_bytes(codon_bytes)

    # 5. lightweight cipher (XOR keystream)
    cipher = lightweight_cipher.xor_encrypt_bytes(compressed, key, img_hash)

    # 6. write file: header + payload
    # header: MAGIC(4) | VERSION(4) | H(4) | W(4) | C(4) | IMG_HASH(32)
    header = MAGIC + struct.pack(">I", VERSION) + struct.pack(">III", shape[0], shape[1], shape[2]) + img_hash
    with open(output_path, "wb") as f:
        f.write(header)
        f.write(cipher)

    print(f"[+] Encrypted saved to {output_path}")

def decrypt_image(input_path: str, key: str, output_path: str):
    with open(input_path, "rb") as f:
        header = f.read(4 + 4 + 4*3 + 32)
        if len(header) < (4+4+12+32):
            raise ValueError("File too short or corrupt")
        magic = header[:4]
        if magic != MAGIC:
            raise ValueError("Not a DNAC file")
        offset = 4
        version = struct.unpack(">I", header[offset:offset+4])[0]; offset += 4
        h, w, c = struct.unpack(">III", header[offset:offset+12]); offset += 12
        img_hash = header[offset:offset+32]
        payload = f.read()

    # 1. decrypt
    compressed = lightweight_cipher.xor_decrypt_bytes(payload, key, img_hash)

    # 2. decompress
    codon_bytes = compression.decompress_bytes(compressed)

    # 3. codon bytes -> dna (uses same key/img_hash to rebuild mapping)
    spliced = protein_encoding.codon_bytes_to_dna(codon_bytes, key, img_hash)

    # 4. reverse chaotic splice
    dna = chaotic_splicing.unsplice_sequence(spliced, key, img_hash)

    # 5. dna -> bytes
    recovered_bytes = dna_mapping.dna_to_binary(dna, key)

    # verify hash
    if hashlib.sha256(recovered_bytes).digest() != img_hash:
        raise ValueError("Integrity check failed: image hash mismatch")

    utils.bytes_to_image_and_save(recovered_bytes, (h, w, c), output_path)
    print(f"[+] Decrypted saved to {output_path}")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--encrypt", help="Path to image to encrypt")
    p.add_argument("--decrypt", help="Path to encrypted file to decrypt")
    p.add_argument("--key", required=True, help="Session key (string)")
    p.add_argument("--output", required=True, help="Output path")
    args = p.parse_args()

    if args.encrypt:
        encrypt_image(args.encrypt, args.key, args.output)
    elif args.decrypt:
        decrypt_image(args.decrypt, args.key, args.output)
    else:
        p.print_help()
