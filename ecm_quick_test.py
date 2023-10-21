"""
A naive implementation of cgbn_stage1
"""

import math

def to_mont(a, modulus):
    pass

def from_mont(a, modulus):
    pass

def double_add_v2(q, u, w, v, d, modulus):
    t = (v + w) % modulus
    v = (v - w) % modulus

    w = (u + q) % modulus
    u = (u - q) % modulus

    CB = (t * u) % modulus
    DA = (v * w) % modulus

    AA = (w * w) % modulus
    BB = (u * u) % modulus

    q = (AA * BB) % modulus

    K = (AA - BB) % modulus

    # d = sigma / 2^32
    dK = (K * d) % modulus

    u = (BB + dK) % modulus
    u = (K * u) % modulus
    w = (DA + CB) % modulus

    v = (DA - CB) % modulus
    w = (w * w) % modulus
    v = (v * v) % modulus
    v = (2*v) % modulus

    return (q, u), (w, v)

def kernel(bits, bRange, sigma, modulus, aX, aZ, bX, bZ):
    inv32 = pow(2 ** 32, -1, modulus)
    assert (2 ** 32 * inv32) % modulus == 1
    d = (sigma * inv32) % modulus

    swapped = 0

    print((aX, aZ), (bX, bZ))
    for b in bRange:
        nth = len(bits) - 1 - b
        bit = bits[nth]

        if bit != swapped:
            swapped = 1 - swapped
            aX, bX = bX, aX
            aZ, bZ = bZ, aZ

        (aX, aZ), (bX, bZ) = double_add_v2(aX, aZ, bX, bZ, d, modulus)
        print("\t", b, bit, "\t", (aX, aZ), (bX, bZ))

    if swapped:
        aX, bX = bX, aX
        aZ, bZ = bZ, aZ


    return aX, aZ, bX, bZ


if 0:
    B = math.prod([8,9,5,7])
    bits = list(map(int, reversed(bin(B)[2:])))
    sigma = 136
    N = 1009*1013

    print(B, bits)
    a = (2, 1)

    # Old way 64*d + 8
    b = (9, (8 + 64 * sigma * pow(2 ** 32, -1, N)) % N)
    t1 = kernel(bits, range(1, len(bits)), sigma, N, *a, *b)
    print()

    # Test
    b = (2, 1)
    a = (0, 1)
    t2 = kernel(bits, range(len(bits)), sigma, N, *a, *b)

    print()
    for x_final, z_final in [t1[0:2]]: #, t2[0:2]]:
        print(x_final, z_final)
        try:
            inverted = pow(z_final, -1, N)
        except:
            print(x_final, z_final, "\t", math.gcd(z_final, N))


N = 1009 * 1013
d = (136 * pow(2**32, -1, N)) % N
if 0:
    for a in range(N):
        for b in range(N):
            t = double_add_v2(2, 1,   a, b, d, N)
            if t[1] == (2, 1):
                print("HI", a, b, "->", t)
                exit()

print(double_add_v2(2, 1, 0, 1, d, N))
