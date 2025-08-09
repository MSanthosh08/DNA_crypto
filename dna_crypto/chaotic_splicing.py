# deterministic permutation via seeded Fisher-Yates
import hmac, hashlib

def _seeded_prng(key_bytes: bytes):
    # returns a generator function that produces integer bytes via HMAC-SHA256 counter
    counter = 0
    def rnd_bytes(n):
        nonlocal counter
        out = b''
        while len(out) < n:
            ctr_b = counter.to_bytes(8, 'big')
            out += hmac.new(key_bytes, ctr_b, hashlib.sha256).digest()
            counter += 1
        return out[:n]
    return rnd_bytes

def _permutation(length: int, key: str, img_hash: bytes):
    seed = hmac.new(key.encode(), img_hash, hashlib.sha256).digest()
    rnd = _seeded_prng(seed)
    arr = list(range(length))
    # Fisher-Yates using PRNG bytes
    i = length
    while i > 1:
        i -= 1
        # get random int in [0, i]
        rbytes = rnd(8)
        r = int.from_bytes(rbytes, 'big') % (i+1)
        arr[i], arr[r] = arr[r], arr[i]
    return arr

def splice_sequence(dna_seq: str, key: str, img_hash: bytes) -> str:
    perm = _permutation(len(dna_seq), key, img_hash)
    return ''.join(dna_seq[i] for i in perm)

def unsplice_sequence(spliced: str, key: str, img_hash: bytes) -> str:
    perm = _permutation(len(spliced), key, img_hash)
    inv = [0]*len(perm)
    for i,p in enumerate(perm):
        inv[p] = i
    return ''.join(spliced[inv[i]] for i in range(len(inv)))
