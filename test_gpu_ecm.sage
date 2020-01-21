import random
import re
import subprocess
import sys

# Start all sigma scans from this value
SIGMA_0 = 1000

# Currently GPU can only do param=3
GPU_PARAM = 3

def GroupOrderParam0(p, sigma):
    K = GF(p)
    v = (4*sigma) * K.gen()
    u = (sigma^2 - 5) * K.gen()
    x = u^3
    b = 4*x*v
    a = (v-u)^3*(3*u+v)
    A = a/b-2
    x = x/v^3
    b = x^3 + A*x^2 + x
    return EllipticCurve([0,b*A,0,b^2,0])

# From parametrizations.c
def GroupOrderParam2(p,sigma):
   K = GF(p)
   E = EllipticCurve(K,[0,36])
   P = sigma*E(-3,3)
   x,y = P.xy()
   x3 = (3*x+y+6)/(2*(y-3))
   A = -(3*x3^4+6*x3^2-1)/(4*x3^3)
   d = K((A+2)/4)

   a = K(4 *d - 2)
   b = K(16*d + 2)
   return EllipticCurve(K,[0,b*a,0,b^2,0])

def GroupOrderParam1(p, sigma):
    K = GF(p)
    inv2_64 = inverse_mod(2 ^ 64, p)
    s = K(sigma^2 * inv2_64)
    A = K(4 *s - 2);
    b = K(16*s + 2);
    return EllipticCurve([0,b*A,0,b^2,0])

def GroupOrderParam3(p, sigma):
    K = GF(p)
    inv2_32 = inverse_mod(2 ^ 32, p)
    s = K(sigma * inv2_32)
    A = K(4 *s - 2);
    b = K(16*s + 2);
    return EllipticCurve([0,b*A,0,b^2,0])

def GroupOrder(param, prime, sigma):
    if param == 0:
        ec = GroupOrderParam0(prime, sigma)
    elif param == 1:
        ec = GroupOrderParam1(prime, sigma)
    elif param == 2:
        ec = GroupOrderParam2(prime, sigma)
    elif param == 3:
        ec = GroupOrderParam3(prime, sigma)
    else:
        raise ValueError("Unknown param: " + str(param))

    return ec.order()

def testInternal():
    # Verify code is working.

    # echo "78257675131877111603" | ecm -sigma "0:1057" 5000 15000-16000
    assert GroupOrder(0, 78257675131877111603, 1057) == \
        2^3 * 3 * 19 * 61 * 149 * 499 * 563 * 4337 * 15497

    # echo "78257675131877111603" | ecm -sigma "1:1511" 17000 25000-26000
    assert GroupOrder(1, 78257675131877111603, 1511) == \
        2^6 * 3 * 5 * 7 * 193 * 293 * 491 * 16651 * 25189

    # echo "78257675131877111603" | ecm -sigma "2:1272" 2000 88000-89000
    assert GroupOrder(2, 78257675131877111603, 1272) == \
        2^2 * 3^6 * 13 * 29 * 491 * 1117 * 1471 * 88237

    # echo "78257675131877111603" | ecm -sigma "3:1203" 6000 73000-74000
    assert GroupOrder(3, 78257675131877111603, 1203) == \
        2^3 * 3^3 * 13 * 139 * 587 * 821 * 5693 * 73079

    # From Zimmermann, https://homepages.cwi.nl/~herman/Zimmermann.pdf
    assert GroupOrder(0, 322410908070969630339041359359164154612901586904078700184707, 20041348) == \
        2^4 * 3^2 * 391063 * 1197631 * 82011967 * 126033329 * 1926338723 * 4654300159 * 51585518429

    # From David Broadhurst, https://lists.gforge.inria.fr/pipermail/ecm-discuss/2005-September/003790.html
    assert GroupOrder(0, 73372650975767950626979890709193208431269141871367229612025497, 175923) == \
        2^2 * 3^2 * 13 * 41 * 3389 * 3989 * 1662013 * 2782993 * 5013037 * 94921033 * 1144363489 * 112303943380877

    # From David Broadhurst, https://lists.gforge.inria.fr/pipermail/ecm-discuss/2005-September/003792.html
    assert GroupOrder(0, 2580118483169716809210552261225054520765090617558895237, 161957884569085) == \
        2^2 * 3 * 1483 * 91381 * 103231 * 239587 * 1151317 * 1186033 * 1611697 * 4199071 * 6941601157

    # From David Broadhurst, https://lists.gforge.inria.fr/pipermail/ecm-discuss/2005-September/003802.html
    assert GroupOrder(0, 6314722182591714308391592266483806595696758378370807102207443753223500809, 2481305347) == \
        2^3 * 3^6 * 11 * 13^2 * 17^4 * 31^2 * 53^2 * 163 * 449 * 853^2 * 3923^2 * 7489 * 11113 * \
        23459^2 * 116531 * 1016891 * 580801721


    print "Implementation successfully tested"


def smallGroupOrders(prime, param, B1, sigma_count):
    ''' Find list of sigmas that will find prime at with B1 <= B1'''
    found = []

    for sigma in range(SIGMA_0, SIGMA_0 + sigma_count):
        order = GroupOrder(param, prime, sigma)
        assert 1 <= order < 2 * prime, (prime, param, sigma, order)

        f = factor(order)
        assert len(f) >= 1, (order, f)

        for p, k in f[:-1]:
            if p ** k > B1:
                break
        else:
            if f[-1][0] ** f[-1][1] <= B1:
                found.append(sigma)

    return found


def stage1Tests(N_size, prime_size, B1, param, sigma_count, seed=None):
    '''
    Generate N such that many sigmas (SIGMA_0:SIGMA_0+sigma_count) have factors
    Verify sigmas found factor in stage 1.
    '''
    assert param == GPU_PARAM, ("GPU only supports param=%d" % GPU_PARAM)
    assert prime_size < N_size
    assert N_size <= 1020
    assert prime_size > 20

    prime_count = N_size // prime_size
    print
    print "Testing GPU stage1: N = %d x %d bits primes @ B1=%d " % (
        prime_count, prime_size, B1)

    import random
    if seed is None:
        seed = random.randrange(2 ^ 32)
    print "\tusing seed:", seed
    random.seed(seed)

    N = []
    factor_by_sigma = {}
    for pi in range(prime_count):
        for test in range(100):
            r = ZZ(random.randint(2 ^ (prime_size-1), 2 ^ prime_size))
            prime = Primes().next(r)
            if prime not in N:
                N.append(prime)
                break
        else:
            raise ValueError("Can't find enought primes at prime_size=%d" % prime_size)

        sigmas = smallGroupOrders(prime, param, B1, sigma_count)
        print "\t%2d: %20d found @ B1=%d by %s" % (pi, prime, B1, sigmas)
        for sigma in sigmas:
            if sigma not in factor_by_sigma:
                factor_by_sigma[sigma] = 1
            factor_by_sigma[sigma] *= prime

    if len(factor_by_sigma) == 0:
        raise ValueError("No primes would be found in step1, lower prime_size or increase B1")

    N_log2 = log(prod(N), 2).n()
    assert N_log2 < N_size, (N_size, N_log2, prime_size, pprime_count)

    N_str = "*".join(map(str, N))
    print
    print "N=" + N_str
    print
    print "Should find:"
    for sigma, factor in sorted(factor_by_sigma.items()):
        print "\tsigma %d => %d" % (sigma, factor)

    # TODO pass in $ECM
    cmd = "echo %s | ./ecm -gpu -gpucurves %d -sigma %d:%d %d" % (
        N_str, sigma_count, param, SIGMA_0, B1)

    print
    print (cmd)

    try:
        output = subprocess.check_output(cmd, shell=True, universal_newlines=True)
        assert False, "Should have factors and had non-zero return"
    except subprocess.CalledProcessError as e:
        assert e.returncode in (2, 6, 8, 10, 14), e.returncode
        lines = e.output.split("\n")

    found_factors = {}
    for line in lines:
        match = re.search("factor ([0-9]*) found.*-sigma [0-4]:([0-9]*)\)", line)
        if match:
            factor, sigma = map(int, match.groups())
            assert sigma not in found_factors
            found_factors[sigma] = factor

    success = found_factors == factor_by_sigma

    all_sigmas = set(factor_by_sigma.keys()) | set(found_factors.keys())
    for sigma in sorted(all_sigmas):
        theory = factor_by_sigma.get(sigma, 0)
        practice = found_factors.get(sigma, 0)
        if theory != practice:
            print "sigma=%d Expected to find %d, found %d" % (sigma, theory, practice)
            assert success == False

    print
    print "Result Matched" if success else "Bad Result"
    if not success:
        sys.exit(1)


if __name__ == "__main__":
    #testInternal()

    # TODO pass -gpucurves, -reps, -seed

    # For Nvidia gtx 970
    sigma_count = 32 #832

    # Test N = 1000 bits, compossed of ~20 50-bit primes, found at B1=1e4
    stage1Tests(1000, 30, 10^4, GPU_PARAM, sigma_count)

    # Test larger primes at B1=1e5
    stage1Tests(1000, 60, 10^5, GPU_PARAM, sigma_count)

