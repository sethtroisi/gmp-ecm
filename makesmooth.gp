/* Simple PARI script to make test numbers for P-1 and P+1. For each prime q 
   in [B1, B2] it prints the n*p, where p is the smallest prime of the 
   form k*q+c. I.e. run with
   echo "makesmooth(1000000, 1001000, 1, 1)" | gp -p 1001000 -q makesmooth | ecm -pm1 100 1001000
   and check that P-1 stage 2 finds all the input numbers as factors.
   To test P+1 properly, a suitable x0 must be generated. The makesmooth_x0()
   function prints it along with the prime. The output can be used like
   echo "makesmooth_x0(1000000, 1001000, -1, 1000000000000273)" | gp -p 1001000 -q makesmooth | while read N X0; do echo $N | ./ecm -pp1 -x0 $X0 1000 1001000; done
   The parameter n can be used to test that the code finds actual factors,
   not just the input number. 
   1000000000000273 +- 1 has relatively large prime factors, making it useful 
   as a value for n. */

makesmooth_one (q, c) = {local(i); i = 1; while (!isprime(q*i+c), i++); return(q*i+c)}
find_x0 (p) = {local(i); i = 3; while (kronecker(i^2-4, p) == 1, i++); return(i);}
makesmooth(B1, B2, c, n) = {local(i); forprime (q = B1, B2, p = makesmooth_one (q, c); print(p * max(n,1)))}
makesmooth_x0(B1, B2, c, n) = {local(i); forprime (q = B1, B2, p = makesmooth_one (q, c); print(p * max(n,1), " ", find_x0(p)))}
