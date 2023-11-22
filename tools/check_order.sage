# script to determine B1/B2 given a prime and sigma
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
#  $ sage check_order.sage <param:sigma> <prime>


# To suppress '***   Warning: increasing stack size to 2000000.' from pari.
pari.allocatemem(4000000, silent=True)

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


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("\t{sys.argv[0]} takes two arguments: <param:sigma> <prime>")

    curve, prime = sys.argv[1:]
    assert ":" in curve, "param:sigma should contain a :"
    assert prime.isnumeric(), "prime should be a number"

    param, sigma = map(int, curve.split(":"))
    prime = int(prime)

    order = GroupOrder(param, prime, sigma)
    factors = order.factor()
    print(factors)
    factors = dict(factors).items()
    B2 = max(f for f, c in factors if c == 1)
    B1 = max(f ** c for f, c in factors if f != B2)
    print(f"Requires {B1=} {B2=}")
