#!/usr/bin/env python3
import argparse
import struct
import hashlib
import hmac
from dna_crypto import utils, dna_mapping, chaotic_splicing, protein_encoding, compression
from dna_crypto import aead_cipher

MAGIC = b'DNAC'
VERSION = 1

# HMAC key derivation
def _hmac_key(key: str) -> bytes:
    return hashlib.sha256(key.encode() + b':hmac').digest()

def _compute_hmac(hmac_key: bytes, data: bytes) -> bytes:
    return hmac.new(hmac_key, data, hashlib.sha256).digest()

def encrypt_image(input_path: str, key: str, output_path: str, no_compress=False, no_splice=False):
    raw_bytes, shape, img_hash = utils.image_to_bytes_and_hash(input_path)

    # bytes -> DNA
    dna = dna_mapping.binary_to_dna(raw_bytes, key)

    # chaotic splice (skip if requested)
    if no_splice:
        spliced = dna
    else:
        spliced = chaotic_splicing.splice_sequence(dna, key, img_hash)

    # codon -> bytes
    codon_bytes, _meta = protein_encoding.dna_to_codon_bytes(spliced, key, img_hash)

    # compress (skip if requested)
    if no_compress:
        compressed = codon_bytes
    else:
        compressed = compression.compress_bytes(codon_bytes)

    # AEAD encrypt compressed payload (include header-less aad = img_hash)
    aead_blob = aead_cipher.aead_encrypt(compressed, key, aad=img_hash)

    # header: MAGIC | VERSION | h | w | c | img_hash(32)
    header = MAGIC + struct.pack(">I", VERSION) + struct.pack(">III", shape[0], shape[1], shape[2]) + img_hash

    # compute HMAC over header + aead_blob
    hmac_key = _hmac_key(key)
    tag = _compute_hmac(hmac_key, header + aead_blob)

    with open(output_path, "wb") as f:
        f.write(header)
        f.write(aead_blob)
        f.write(tag)

    print(f"[+] Encrypted saved to {output_path} (no_compress={no_compress} no_splice={no_splice})")


def decrypt_image(input_path: str, key: str, output_path: str):
    with open(input_path, "rb") as f:
        header = f.read(4 + 4 + 12 + 32)
        if len(header) < (4+4+12+32):
            raise ValueError("file too short")
        magic = header[:4]
        if magic != MAGIC:
            raise ValueError("Bad file magic")
        offset = 4
        version = struct.unpack(">I", header[offset:offset+4])[0]; offset += 4
        h, w, c = struct.unpack(">III", header[offset:offset+12]); offset += 12
        img_hash = header[offset:offset+32]
        rest = f.read()
    # last 32 bytes are HMAC tag
    if len(rest) < 32:
        raise ValueError("missing hmac")
    aead_blob = rest[:-32]
    tag = rest[-32:]

    # verify HMAC
    hmac_key = _hmac_key(key)
    expected = _compute_hmac(hmac_key, header + aead_blob)
    if not hmac.compare_digest(expected, tag):
        raise ValueError("HMAC verification failed: file tampered or wrong key")

    # AEAD decrypt
    compressed = aead_cipher.aead_decrypt(aead_blob, key, aad=img_hash)

    # try decompress; if zlib fails we assume it wasn't compressed:
    try:
        codon_bytes = compression.decompress_bytes(compressed)
        used_compression = True
    except Exception:
        codon_bytes = compressed
        used_compression = False

    # codon bytes -> dna
    spliced = protein_encoding.codon_bytes_to_dna(codon_bytes, key, img_hash)

    # unsplice
    dna = chaotic_splicing.unsplice_sequence(spliced, key, img_hash)

    # dna -> raw bytes
    recovered = dna_mapping.dna_to_binary(dna, key)

    # verify img hash
    if hashlib.sha256(recovered).digest() != img_hash:
        raise ValueError("Image integrity failed after decryption")

    utils.bytes_to_image_and_save(recovered, (h, w, c), output_path)
    print(f"[+] Decrypted saved to {output_path} (compression used: {used_compression})")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--encrypt", help="image path to encrypt")
    p.add_argument("--decrypt", help="encrypted file to decrypt")
    p.add_argument("--key", required=True, help="session key")
    p.add_argument("--output", required=True, help="output path")
    p.add_argument("--no-compress", action="store_true", help="skip compression")
    p.add_argument("--no-splice", action="store_true", help="skip chaotic splicing")
    args = p.parse_args()

    if args.encrypt:
        encrypt_image(args.encrypt, args.key, args.output, no_compress=args.no_compress, no_splice=args.no_splice)
    elif args.decrypt:
        decrypt_image(args.decrypt, args.key, args.output)
    elif args.encrypt and args.decrypt:
        p.print_help()
        raise ValueError("You must specify either --encrypt or --decrypt")
    else:
        p.print_help()
        raise ValueError("You cannot specify both --encrypt and --decrypt at the same time")
