# script to spot check ecm -gpu stage 1 implementation.
#
# Computes order of a elliptic curve mod a prime to determine if that prime 
# will be found at a specified B1 level, then runs ecm -gpu to verify that
# those primes are indeed found. 
#
# to perform a quick test: 5 iterations of 32 curves at B1=10000 and B1=100000
#  $ sage check_gpuecm.sage ./ecm
# for even faster check
#  $ sage check_gpuecm.sage --iterations 1 ./ecm
#
# to test with ECM_GPU_NB_DIGITS=16 and extra verbose output
#  $ sage check_gpuecm.sage --nbits 512 -v -v ./ecm
# 
# to test more iterations (~12 minutes) or more curves per iteration
#  $ sage check_gpuecm.sage --iterations 20 ./ecm
#  $ sage check_gpuecm.sage --gpucurves 128 ./ecm
#
# to see all options  
#  $ sage check_gpuecm.sage -h
#
# On error, will exit with return code 1 and stdout
#    'Wrong result for seed=<X>'
# after checking with a clean build please email Seth Troisi or ecm-discuss
# with seed (and any non standard options you may have used).


# To suppress '***   Warning: increasing stack size to 2000000.' from pari.
pari.allocatemem(4000000, silent=True)

import argparse
import random
import re
import subprocess
import sys

parser = argparse.ArgumentParser(description='Spot check for ecm -gpu stage1')

parser.add_argument('ecm_cmd', type=str,
    help='which ecm (e.g. ecm, ./ecm) to run')

parser.add_argument('--iterations', type=int, default=5,
    help='Number of tests (small + large) to perform')

parser.add_argument('-c', '--gpucurves', type=int, default=32,
    help='number of curves to test in a batch [default: 32]')

parser.add_argument('--B1', type=int, default=10000,
    help='B1 to test at [default: 10,000]')

parser.add_argument('--seed', type=int, default=None,
    help='Random seed')

parser.add_argument('--nbits', '-n', type=int, default=1024,
    help='Only needed if ECM_GPU_NB_DIGITS was adjusted')

parser.add_argument('--timing', action='store_true',
    help='Producing timing information')

parser.add_argument('--verbose', '-v', action='count', default=1,
    help='Print more output (pass -v -v for even more)')
parser.add_argument('--quiet', '-q',
    action='store_const', const=0, dest='verbose',
    help='Suppress most output')



# Currently GPU can only do param=3
GPU_PARAM = 3

FACTOR_FOUND_RE = re.compile(
    'factor ([0-9]*) found in Step 1.*-sigma [0-4]:([0-9]*)\)')


def CurveParam0(p, sigma):
    K = GF(p)
    v = K(4*sigma)
    u = K(sigma^2 - 5)
    x = u^3
    b = 4*x*v
    a = (v-u)^3*(3*u+v)
    A = a/b-2
    x = x/v^3
    b = x^3 + A*x^2 + x
    return EllipticCurve(K,[0,A/b,0,1/b^2,0])

def PointOrderS(p, s):
    K = GF(p)
    A = K(4 *s - 2)
    b = K(16*s + 2)
    E = EllipticCurve(K,[0,A/b,0,1/b^2,0])
    return E(2/b,1/b) # x0=2, y0=1

# From parametrizations.c
def CurveParam2(p, sigma):
    K = GF(p)
    E = EllipticCurve(K,[0,36])
    P = sigma*E(-3,3)
    x,y = P.xy()
    x3 = (3*x+y+6)/(2*(y-3))
    A = -(3*x3^4+6*x3^2-1)/(4*x3^3)
    d = K((A+2)/4)
    return PointOrderS(p, d)

def CurveParam1(p, sigma):
    K = GF(p)
    return PointOrderS(p, K(sigma^2 / 2^64))

def CurveParam3(p, sigma):
    K = GF(p)
    return PointOrderS(p, K(sigma / 2^32))

def GroupOrder(param, prime, sigma):
    if param == 0:
        ec = CurveParam0(prime, sigma)
    elif param == 1:
        ec = CurveParam1(prime, sigma)
    elif param == 2:
        ec = CurveParam2(prime, sigma)
    elif param == 3:
        ec = CurveParam3(prime, sigma)
    else:
        raise ValueError('Unknown param: ' + str(param))

    return ec.order()


def testInternal():
    # Verify code is working.

    # echo '78257675131877111603' | ecm -sigma '0:4528' 3100 8000-9000
    assert GroupOrder(0, 78257675131877111603, 4528) == \
        2^3 * 3 * 71 * 563 * 1531 * 2153 * 3011 * 8219
    # echo '78257675131877111603' | ecm -sigma '1:3396' 4800 8000-9000
    assert GroupOrder(1, 78257675131877111603, 3396) == \
        2 * 3 * 223 * 271 * 811 * 821 * 4799 * 8443
    # echo '78257675131877111603' | ecm -sigma '2:1801' 2100 9000-9300
    assert GroupOrder(2, 78257675131877111603, 1801) == \
        2^8 * 3^3 * 23 * 41 * 47 * 67 * 101 * 2039 * 9257
    # echo '78257675131877111603' | ecm -sigma '3:2012' 6000 8000-9000
    assert GroupOrder(3, 78257675131877111603, 2012) == \
        2^2 * 5^3 * 43^2 * 71 * 83 * 139 * 5779 * 8941

    # echo '1082500099132634560519' | ecm -sigma '0:6677' 1200 5000-6000
    assert GroupOrder(0, 1082500099132634560519, 6677) == \
        2^2 * 3 * 5^2 * 7 * 139 * 677 * 887 * 947 * 1123 * 5807
    # echo '1082500099132634560519' | ecm -sigma '1:1800' 4000 6000-7000
    assert GroupOrder(1, 1082500099132634560519, 1800) == \
        2^5 * 13 * 17 * 79 * 701 * 2647 * 3347 * 6367
    # echo '1082500099132634560519' | ecm -sigma '2:2966' 2000 7000-8000
    assert GroupOrder(2, 1082500099132634560519, 2966) == \
        2^3 * 3 * 29 * 31^2 * 61 * 109 * 487 * 1709 * 7499
    # echo '1082500099132634560519' | ecm -sigma '3:1600' 2000 3000-3100
    assert GroupOrder(3, 1082500099132634560519, 1600) == \
        3^3 * 7 * 23^2 * 37 * 67 * 71 * 1297 * 1933 * 3067

    # From Zimmermann, https://homepages.cwi.nl/~herman/Zimmermann.pdf
    assert GroupOrder(0, 322410908070969630339041359359164154612901586904078700184707, 20041348) == \
        2^4 * 3^2 * 391063 * 1197631 * 82011967 * 126033329 * 1926338723 * 4654300159 * 51585518429

    # From David Broadhurst, https://sympa.inria.fr/sympa/arc/ecm-discuss/2005-09/msg00022.html
    assert GroupOrder(0, 2580118483169716809210552261225054520765090617558895237, 161957884569085) == \
        2^2 * 3 * 1483 * 91381 * 103231 * 239587 * 1151317 * 1186033 * 1611697 * 4199071 * 6941601157

    # From David Broadhurst, https://sympa.inria.fr/sympa/arc/ecm-discuss/2005-09/msg00020.html
    assert GroupOrder(0, 73372650975767950626979890709193208431269141871367229612025497, 175923) == \
        2^2 * 3^2 * 13 * 41 * 3389 * 3989 * 1662013 * 2782993 * 5013037 * 94921033 * 1144363489 * 112303943380877

    # From David Broadhurst, https://sympa.inria.fr/sympa/arc/ecm-discuss/2005-09/msg00032.html
    assert GroupOrder(0, 6314722182591714308391592266483806595696758378370807102207443753223500809, 2481305347) == \
        2^3 * 3^6 * 11 * 13^2 * 17^4 * 31^2 * 53^2 * 163 * 449 * 853^2 * 3923^2 * 7489 * 11113 * \
        23459^2 * 116531 * 1016891 * 580801721

    if args.verbose:
        print('Implementation successfully tested\n')


def smallestGroupOrder(prime, param, sigma_0, num_curves):
    '''
    Print smallest B1 (and B2) needed to find a curve from sigma_0 ... sigma_0 + curves
    At the end return the smallest sigma
    '''
    # This is not used directly in check_gpuecm.sage, but is very useful for building tests
    found_at = None
    smallest = 10 ^ 100
    for sigma in range(sigma_0, sigma_0 + num_curves):
        order = GroupOrder(param, prime, sigma)
        assert 1 <= order < 2 * prime, (prime, param, sigma, order)

        f = factor(order)
        # Largest prime
        exp = sorted(p ** k for p, k in f)
        # B1 to find this factor is 2nd smallest if B2 = 500 * B1 finds large factor
        min_b1 = exp[-2] if len(exp) >= 2 and exp[-1] < exp[-2] * 500 else exp[-1]
        if min_b1 <= smallest:
            print ("\tFound by sigma: %d with B1=%d, B2=%d" % (sigma, min_b1, exp[-1]))
            smallest = min_b1
            found_at = sigma
    return found_at


def smallGroupOrders(prime, param, sigma_0, test_B1, num_curves):
    '''Find list of sigmas that will find prime with B1 >= test_B1'''
    for sigma in range(sigma_0, sigma_0 + num_curves):
        order = GroupOrder(param, prime, sigma)
        assert 1 <= order < 2 * prime, (prime, param, sigma, order)

        f = factor(order)
        assert len(f) >= 1, (order, f)

        for p, k in f:
            if p ^ k > test_B1:
                break
        else:
            yield sigma


def findPrimesOfSize(count, prime_size):
    '''Find count primes with prime_size bits'''
    primes = set()
    for pi in range(count):
        for test in range(100):
            r = ZZ(random.randint(2 ^ (prime_size-1), 2 ^ prime_size))
            prime = Primes().next(r)
            if prime not in primes:
                primes.add(prime)
                break
        else:
            raise ValueError("Can't find enought primes at prime_size=%d" %
                prime_size)
    return sorted(primes)


def expectedFactorsBySigma(args, primes, param, sigma_0, B1):
    '''Calculate which primes will be found by which curves'''
    factor_by_sigma = {}
    for i, prime in enumerate(primes):
        # This can be slower than the actual ECM!
        sigmas = smallGroupOrders(prime, param, sigma_0, B1, args.gpucurves)
        if args.verbose > 2:
            print('\t%2d: %20d found @ B1=%d by %s' % (i, prime, B1, sigmas))
        for sigma in sigmas:
            if sigma not in factor_by_sigma:
                factor_by_sigma[sigma] = 1
            factor_by_sigma[sigma] *= prime
    return factor_by_sigma


def verifyFactorsFoundBySigma(
        args, primes, param, sigma_0, B1, factor_by_sigma):
    '''Verify the expected factors where found by ECM'''

    if not any(found for found in factor_by_sigma.values() if found > 1):
        raise ValueError(
            'No primes would be found in stage 1, '
            'lower prime_size or increase B1(%d)' % B1)

    N_log2 = log(prod(primes), 2).n()
    assert N_log2 < args.nbits, (args.nbits, N_log2)

    N_str = '*'.join(map(str, primes))
    if args.verbose > 1:
        sigma_str = ', '.join(map(str, sorted(factor_by_sigma)))
        print('\tSigmas with factors: %s' % sigma_str)

    echo_cmd = 'echo %s | ' % N_str
    ecm_cmd = '%s -gpu -gpucurves %d -sigma %d:%d %d 0' % (
        args.ecm_cmd, args.gpucurves, param, sigma_0, B1)


    if args.verbose > 2:
        print('\t' + echo_cmd + ecm_cmd)
    elif args.verbose:
        print('\t' + ecm_cmd)

    try:
        output = subprocess.check_output(
            echo_cmd + ecm_cmd, shell=True, universal_newlines=True)
        assert False, 'Should have factors and had non-zero return'
    except subprocess.CalledProcessError as e:
        assert e.returncode in (2, 6, 8, 10, 14), e.returncode
        lines = e.output.split('\n')

    found_factors = {}
    for line in lines:
        match = FACTOR_FOUND_RE.search(line)
        if match:
            f, sigma = map(int, match.groups())
            assert sigma not in found_factors
            found_factors[sigma] = f

    all_sigmas = set(factor_by_sigma.keys()) | set(found_factors.keys())
    for sigma in sorted(all_sigmas):
        theory = factor_by_sigma.get(sigma, 1)
        practice = found_factors.get(sigma, 1)
        if theory != practice:
            if theory % practice == 0:
                print('sigma=%d Expected to find %d, found %d' %
                    (sigma, theory, practice))
            elif practice % theory == 0:
                extra = practice / theory
                f = factor(GroupOrder(param, extra, sigma))
                print('\tExtra factor (%d) found by sigma=%d '
                      'expected order=%s' % (extra, sigma, f))
            else:
                print('MAJOR MISMATCH: %d vs %d' % (
                    factor(theory), factor(practice)))

    expected_curves = len(factor_by_sigma)
    perfect_match = factor_by_sigma == found_factors
    if perfect_match:
        if args.verbose:
            print('Results matched exactly (%d curves found factors)' %
                expected_curves)
    else:
        print('Wrong results for seed=%d' % seed)
        print('\t' + echo_cmd + ecm_cmd)
        sys.exit(1)

    if args.verbose:
        print('')

    return len(found_factors)


def stage1Tests(args, prime_size, B1, param, seed):
    '''
    Generate N such that many sigmas (sigma_0:sigma_0+args.gpucurves) have factors
    Verify sigmas found factor in stage 1.
    '''
    assert param == GPU_PARAM, ('GPU only supports param=%d' % GPU_PARAM)
    assert prime_size < args.nbits
    assert args.nbits <= 1020
    assert prime_size > 20

    prime_count = args.nbits // prime_size
    if args.verbose:
        print('Testing GPU stage1: N = %d x %d bits primes @ B1=%d ' % (
            prime_count, prime_size, B1))

    if args.verbose > 1:
        print('\tusing seed: %s' % seed)

    random.seed(seed)
    sigma_0 = random.randrange(1000, 2^31)

    primes = findPrimesOfSize(prime_count, prime_size)
    factor_by_sigma = expectedFactorsBySigma(args, primes, param, sigma_0, B1)

    return verifyFactorsFoundBySigma(
        args, primes, param, sigma_0, B1, factor_by_sigma)


def overflowTest(args, param, seed):
    '''
    Generate N such that N is VERY close to nbits
    Verify small factors found
    '''

    random.seed(seed)
    sigma_0 = random.randrange(1000, 2^31)

    # Multiply a handful of small primes that "should" be found by each sigma
    # Then pad out N with a giant prime
    primes = findPrimesOfSize(count=256//12, prime_size=12)
    expected = prod(primes)

    pad_limit = 2 ** args.nbits// expected
    large_prime = Primes().next(ZZ(random.randint(pad_limit // 2, pad_limit)))

    N_log2 = log(expected * large_prime, 2).n()
    # within 1 bits of the limit
    assert 0 < args.nbits - N_log2 < 1, (args.nbits, N_log2)
    if args.verbose:
        print ("Checking overflow with log2(N) = %.2f" % N_log2)

    # Expect all the small primes to be found by all sigmas
    factor_by_sigma = {}
    for sigma in range(sigma_0, sigma_0 + args.gpucurves):
        factor_by_sigma[sigma] = expected

    return verifyFactorsFoundBySigma(
        args, primes + [large_prime], param, sigma_0, args.B1, factor_by_sigma)


def timingTest(args, N_sizes):
    '''Produce some timing information on ecm'''
    primes = [Primes().next(ZZ(int(2 ** (n - 0.5)))) for n in N_sizes]

    # Doesn't matter
    sigma = "%d:%d" % (GPU_PARAM, 12)

    for size, prime in zip(N_sizes, primes):
        cmd = 'echo %s | %s -gpu -gpucurves %d -sigma %s %d 0' % (
            prime, args.ecm_cmd, args.gpucurves, sigma, args.B1)
        try:
            output = subprocess.check_output(cmd, shell=True, universal_newlines=True)
            timing = min(line for line in output.split('\n') if line.startswith("Computing "))
            timing = timing[timing.index("took") + 4:].strip()
        except subprocess.CalledProcessError as e:
            timing = "error"
        debug = 'N=%d bits, B1=%d, curves=%d, ' % (size, args.B1, args.gpucurves)
        print (debug, timing)


if __name__ == '__main__':
    args = parser.parse_args()

    testInternal()

    seed = args.seed
    if seed is None:
        seed = random.randrange(2 ^ 32)

    if args.timing:
        timingTest(args, [250, 500, 1000, 1500, 2000])
        exit()

    # GPU needs 6 bits for carry / temp results
    args.nbits -= 6
    overflowTest(args, GPU_PARAM, seed)

    if args.iterations:
        found = 0
        for i in range(args.iterations):
            # Test smallish primes (40 bits = 12 digit) at B1 (default: 10^4)
            found += stage1Tests(args, 40, args.B1, GPU_PARAM, seed)
            seed += 1

            # Test larger primes at 10xB1
            found += stage1Tests(args, 60, 10*args.B1, GPU_PARAM, seed)
            seed += 1

        print('Results matched in %d tests (%d curves found factors)' %
            (2*args.iterations, found))
