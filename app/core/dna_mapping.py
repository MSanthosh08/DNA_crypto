def binary_to_dna(raw_bytes: bytes, key: str) -> str:
    bits = ''.join(f'{b:08b}' for b in raw_bytes)
    mapping = {'00':'A','01':'C','10':'G','11':'T'}
    dna = ''.join(mapping[bits[i:i+2]] for i in range(0, len(bits), 2))
    return dna

def dna_to_binary(dna_seq: str, key: str) -> bytes:
    rev = {'A':'00','C':'01','G':'10','T':'11'}
    bits = ''.join(rev[b] for b in dna_seq)
    if len(bits) % 8 != 0:
        bits = bits.ljust(((len(bits)//8)+1)*8, '0')
    out = bytearray()
    for i in range(0, len(bits), 8):
        out.append(int(bits[i:i+8], 2))
    return bytes(out)
