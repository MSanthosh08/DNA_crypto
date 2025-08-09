# Map each codon (3 bases) to a unique byte 0..63 using a session permutation.
import itertools, hmac, hashlib

BASES = ['A','C','G','T']
CODONS = [''.join(p) for p in itertools.product(BASES, repeat=3)]

def _codon_permutation(key: str, img_hash: bytes):
    seed = hmac.new(key.encode(), img_hash, hashlib.sha256).digest()
    # derive permutation by Fisher-Yates over CODONS using seed-based PRNG
    perm = list(range(64))
    # PRNG similar to chaotic_splicing
    counter = 0
    def rnd_int(max_inclusive):
        nonlocal counter
        out = hmac.new(seed, counter.to_bytes(8,'big'), hashlib.sha256).digest()
        counter += 1
        return int.from_bytes(out, 'big') % (max_inclusive+1)
    i = 64
    while i > 1:
        i -= 1
        r = rnd_int(i)
        perm[i], perm[r] = perm[r], perm[i]
    return perm

def dna_to_codon_bytes(dna_seq: str, key: str, img_hash: bytes) -> (bytes, dict):
    # pad dna to multiple of 3 with 'A' (lossless because we store shape/hash)
    if len(dna_seq) % 3 != 0:
        dna_seq = dna_seq + 'A' * (3 - (len(dna_seq) % 3))
    perm = _codon_permutation(key, img_hash)
    # codon_index_map: codon -> permuted byte
    # base ordering of CODONS is 0..63
    codon_to_byte = {}
    for i, idx in enumerate(perm):
        codon_to_byte[CODONS[idx]] = i  # map actual codon to byte value i
    out = bytearray()
    for i in range(0, len(dna_seq), 3):
        codon = dna_seq[i:i+3]
        out.append(codon_to_byte[codon])
    meta = {'perm': perm}  # not strictly needed (rederived from key+img_hash)
    return bytes(out), meta

def codon_bytes_to_dna(codon_bytes: bytes, key: str, img_hash: bytes) -> str:
    perm = _codon_permutation(key, img_hash)
    # invert mapping: byte i -> CODONS[ perm[i] ]
    inv = [None]*64
    for i, idx in enumerate(perm):
        inv[i] = CODONS[idx]
    dna = ''.join(inv[b] for b in codon_bytes)
    return dna
