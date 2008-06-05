/* Simple PARI script to make test numbers for P-1 and P+1. For each prime q 
   in [B1, B2] it prints the n*p, where p is the smallest prime of the 
   form k*q+c. I.e. run with
   echo "makesmooth(1000000, 1001000, 1, 1)" | gp -p 1001000 -q makesmooth | ecm -pm1 100 1001000
   and check that P-1 stage 2 finds all the input numbers as factors.

   To test P+1 properly, a suitable x0 must be generated. The makesmooth_x0()
   function prints it along with the prime. The output can be used like
   echo "makesmooth_x0(1000000, 1001000, -1, 1000000000000273)" | gp -p 1001000 -q makesmooth | while read N X0; do echo $N | ./ecm -pp1 -x0 $X0 1000 1001000; done
   Alternatively, the function makesmooth_fixed_x0() can be used which outputs
   only primes so that x0^2-4 is a quadratic non-residue, i.e. by
   echo "makesmooth_fixed_x0(1000000, 1001000, -1, 3, 1000000000000273)" | gp -p 1001000 -q makesmooth | ecm -pm1 -x0 3 100 1001000
   will make sure that for all produced primes, 3^2-4 = 5 is a QNR and that
   GMP-ECM with -x0 3 parameter will work properly.

   The parameter n can be used to test that the code finds actual factors,
   not just the input number. 
   1000000000000273 +- 1 has relatively large prime factors, making it useful 
   as a value for n. */

/* Find the smallest prime of the form p = q*i+c */
makesmooth_one (q, c) = {local(i); i = 1; while (!isprime(q*i+c), i++); return(q*i+c)}

/* Find the smallest prime of the form p = q*i+c so that D is a QNR (mod p) */
makesmooth_one_D (q, c, D) = {local(p); p=q+c; while (!isprime(p) || kronecker(D, p) != -1, p+=q); return(p)}

find_x0 (p) = {local(i); i = 3; while (kronecker(i^2-4, p) == 1, i++); return(i);}
makesmooth(B1, B2, c, n) = {local(i, q); q = nextprime(B1); while (q <= B2, p = makesmooth_one (q, c); print(p * max(n,1)); q = nextprime(q + 1))}
makesmooth_x0(B1, B2, c, n) = {local(i, q); q = nextprime(B1); while (q <= B2, p = makesmooth_one (q, c); print(p * max(n,1), " ", find_x0(p)); q = nextprime(q + 1))}
makesmooth_fixed_x0(B1, B2, c, x0, n) = {local(i, q); q = nextprime(B1); while (q <= B2, p = makesmooth_one_D (q, c, x0^2-4); print(p * max(n,1)); q = nextprime(q + 1))}
