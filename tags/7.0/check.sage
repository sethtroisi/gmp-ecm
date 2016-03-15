def FindGroupOrder(p,s):
   K = GF(p)
   v = K(4*s)
   u = K(s^2-5)
   x = u^3
   b = 4*x*v
   a = (v-u)^3*(3*u+v)
   A = a/b-2
   x = x/v^3
   b = x^3 + A*x^2 + x
   E = EllipticCurve(K,[0,b*A,0,b^2,0])
   return factor(E.cardinality())

def FindGroupOrderA(p,A):
   K = GF(p)
   d = K((A+2)/4)
   a = K(4*d-2)
   b = K(16*d+2)
   E = EllipticCurve(K,[0,a/b,0,1/b^2,0])
   return factor(E.cardinality())

# for parameter sigma = 1:s
def FindGroupOrderParam1(p,s):
   return FindGroupOrderA (p, 4*s^2/2^64-2)

# for parameter sigma = 2:s
def FindGroupOrderParam2(p,s):
   K = GF(p)
   E = EllipticCurve(K,[0,36])
   P = s*E(-3,3)
   x,y = P.xy()
   x3 = (3*x+y+6)/(2*(y-3))
   A = -(3*x3^4+6*x3^2-1)/(4*x3^3)
   d = K((A+2)/4)
   a = K(4*d-2)
   b = K(16*d+2)
   E = EllipticCurve(K,[0,a/b,0,1/b^2,0])
   return factor(E.cardinality())

# for parameter sigma = 3:s
def FindGroupOrderParam3(p,s):
   K = GF(p)
   d = K(s/2^32)
   A = 4*d-2
   B = 4*A+10
   E = EllipticCurve(K,[0,A/B,0,1/B^2,0])
   return factor(E.cardinality())

def FindGroupOrderParam (p, sigma, param):
   if param == 0:
      return FindGroupOrder (p, sigma)
   elif param == 1:
      return FindGroupOrderParam1 (p, sigma)
   elif param == 2:
      return FindGroupOrderParam2 (p, sigma)
   elif param == 3:
      return FindGroupOrderParam3 (p, sigma)
   else:
      print "Invalid parametrization: ", param
      raise ValueError

# check if the prime p is found with B1,B2,param,sigma, or raises an error
# check_found_aux ("./ecm", 31622776601683800097, 11000, 1873422, 1, 800667805)
# check_found_aux ("./ecm", 31622776601683800097, 11000, 1873422, 1, 800667806)
def check_found_aux (ecm, p, B1, B2, param, sigma):
   f = open("/tmp/inxyz", "w")
   f.write(str(p) + "\n")
   f.close()
   f = open("/tmp/doitxyz", "w")
   f.write(ecm + " -param " + str(param) + " -sigma " + str(sigma) + " " + str(B1) + " " + str(B2) + " < /tmp/inxyz > /tmp/outxyz\n")
   f.close()
   os.system("chmod +x /tmp/doitxyz")
   os.system("/tmp/doitxyz")
   f = open("/tmp/outxyz", "r")
   l = f.readlines()
   f.close()
   n = len(l)
   if l[n-1] <> 'Found input number N\n':
      print "prime p=", p, "not found with B1=", B1, "B2=", B2, "param=", param, "sigma=", sigma
      raise ValueError

def is_found(l, B1, B2):
   n = len(l)
   if l[n-1][0] > B2:
      return False
   for i in range(n-2,-1,-1):
      if l[i][0]^l[i][1] > B1:
         return False
   return True
  
# check if a prime p is found with bounds B1 and B2,
# for parameter 'param' and sigma in [sigma_min,sigma_max-1]
# check_found ("./ecm", 31622776601683800097, 11000, 1873422, 0, 1000)
# check_found ("./ecm", 31622776601683800097, 11000, 1873422, 1, 1000)
# check_found ("./ecm", 31622776601683800097, 11000, 1873422, 2, 1000)
# check_found ("./ecm", 31622776601683800097, 11000, 1873422, 3, 1000)
def check_found (ecm, p, B1, B2, param, sigma_max):
   assert (is_prime (p))
   e2 = 0
   e3 = 0
   tries = 0
   for sigma in range(sigma_max):
      try:
         l = FindGroupOrderParam (p, sigma, param)
      except ArithmeticError:
         continue
      tries += 1
      assert (l[0][0] == 2)
      e2 += l[0][1]
      if l[1][0] == 3:
         e3 += l[1][1]
      if is_found (l, B1, B2):
         # check the factor is really found
         check_found_aux (ecm, p, B1, B2, param, sigma)
   print tries, 1.0*e2/tries, 1.0*e3/tries, 2.0^(e2/tries)*3.0^(e3/tries)

# check all parametrizations 0, 1, 2, 3
# check_found_all ("./ecm", 31622776601683800097, 11000, 1873422, 1000)
def check_found_all (ecm, p, B1, B2, sigma_max):
   for param in range(4):
      check_found (ecm, p, B1, B2, param, sigma_max)
