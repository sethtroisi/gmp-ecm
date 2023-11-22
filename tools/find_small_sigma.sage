# script to find a curve with small B1/B2 given a prime factor
#
# Copyright 2023
# Seth Troisi
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; see the file COPYING.  If not, see
# http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
# 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
#
# Computes order of a elliptic curve mod a prime and determine neede B1/B2
#
# to find a good sigma for param 3
#  $ sage find_small_sigma.sage 59649589127497217 --sigma 3
#
# With optional starting sigma and param 0
#  $ sage find_small_sigma.sage 59649589127497217 --sigma 0:1000
#
# Only search 1000 curves
#  $ sage find_small_sigma.sage 59649589127497217 --curves 1000
#
# With max B2
#  $ sage find_small_sigma.sage 59649589127497217 --sigma 0:1000 -B2lim 1e6
#
# to see all options
#  $ sage find_small_sigma.sage


# To suppress '***   Warning: increasing stack size to 2000000.' from pari.
pari.allocatemem(4000000, silent=True)

import argparse
from operator import itemgetter

parser = argparse.ArgumentParser(description='Spot check for ecm -gpu stage1')

parser.add_argument('prime', type=int,
    help='prime factor to find')

parser.add_argument('-c', '--curves', type=int, default=10000,
    help='Number of curves to test')

parser.add_argument('--B1lim', type=int, default=1000000,
    help='B1 Limit [default: 1,000,000]')

parser.add_argument('--B2lim', type=int, default=1000000000,
    help='B1 Limit [default: 1,000,000,000]')

parser.add_argument('--sigma', type=str,
    help='starting sigma or param:sigma')

parser.add_argument('--verbose', '-v', action='count', default=1,
    help='Print more output (pass -v -v for even more)')
parser.add_argument('--quiet', '-q',
    action='store_const', const=0, dest='verbose',
    help='Suppress most output')


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


def findSmallGroupOrder(prime, param, sigma_0, num_curves, B1lim, B2lim):
    def best_B1(r):
        return min(r, key=operator.itemgetter(0))[2]
    def best_B2(r):
        return min(r, key=operator.itemgetter(1))[2]
    def best_conforming(r):
        best =min(((b2, s) for b1, b2, s, _ in r if b1 < B1lim and b2 < B2lim), default=None)
        return best[1] if best else None
    def best_conforming_B1(r):
        best = min(((b2, s) for b1, b2, s, _ in r if b1 < B1lim), default=None)
        return best[1] if best else None
    def best_conforming_B2(r):
        best = min(((b1, s) for b1, b2, s, _ in r if b2 < B2lim), default=None)
        return best[1] if best else None

    def get_b1_b2_order(sigma):
        order = GroupOrder(param, prime, sigma)
        assert 1 <= order < 2 * prime, (prime, param, sigma, order)

        f = factor(order)
        min_b2 = max(p ** e for p, e in f)
        min_b1 = max((p ** e for p, e in f if p != min_b2), default=min_b2)
        return min_b1, min_b2, sigma, list(f)

    bests = [
        best_B2, best_conforming, best_conforming_B1, best_conforming_B2
    ]

    results = [get_b1_b2_order(sigma_0)]
    b1, b2, s, _ = results[0]
    print(f"First results sigma {param}:{s} => B1={b1:,}, B2={b2:,}")

    for sigma in range(sigma_0 + 1, sigma_0 + num_curves):
        r = get_b1_b2_order(sigma)
        b1, b2, s, _ = r
        results.append(r)

        for best_fn in bests:
           best = best_fn(results)
           if best == s:
                name = best_fn.__name__.replace('_', ' ')
                print(f"With sigma={param}:{s} | New {name}: B1={b1:,}  B2={b2:,}")
                break

    print("\n")
    for best_fn in [best_B1] + bests:
       best = best_fn(results)
       if best:
         b1, b2, s, _ = min(r for r in results if r[2] == best)
         name = best_fn.__name__.replace('_', ' ')
         print(f"With sigma={param}:{s} | {name}: B1={b1:,}  B2={b2:,}")

    best = best_conforming(results)
    if best:
      return best

    best = best_B2(results)
    return best


if __name__ == '__main__':
    args = parser.parse_args()

    testInternal()

    param = 0
    sigma_0 = 0

    if not args.sigma:
        sigma_0 = 100
    elif args.sigma.isnumeric():
        sigma_0 = int(args.sigma)
    elif ":" in args.sigma:
        param, sigma_0 = map(int, args.sigma.split(":"))
    else:
        raise ValueError("invalid --sigma {}".format(args.sigma))

    print("Searching up to {} curves over {}:{} that finds".format(
        args.curves, param, sigma_0))
    print("{} with B1 <= {:,}, B2 <= {:,}".format(
        args.prime, args.B1lim, args.B2lim))


    findSmallGroupOrder(args.prime, param, sigma_0, args.curves, args.B1lim, args.B2lim)

