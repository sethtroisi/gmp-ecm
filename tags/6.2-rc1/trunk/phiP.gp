
largest_primefactor(n) = vecmax(factorint(n)~[1,])

/* We examine P values with eulerphi(P) in [mini, maxi], 
   maxi = mini + phimult*(len - 1) */
/* We assume that both phiP and mini are divisible my phimult, so the array 
   entry phi[i] contains a P so that phiP = eulerphi(P) and 
   (phiP - mini) / phimult + 1 = i. (phi[i] contains zero if no such P was 
   found). Conversely, phiP = (i - 1) * phimult + mini */

/* mini = 150000000; */
/* len = 10000000; */
/* multP = 3*5*7*11; */

make_phiP (mini, len, multP, oldbest) = 
{
  local (phimult, maxi, phi, P, phiP, best);
  phimult = eulerphi(multP);
  best = oldbest;
  if (mini % phimult != 0, 
    error("mini = ", mini, " is not a multiple of phimult = ", phimult));
  maxi = mini + phimult * (len - 1);
  phi = vector(len);

  /* For each candidate odd P value <= ~4*maxi and a multiple of multP, 
     if eulerphi(P) is in the [mini, maxi] range, store the P value at 
     phi[i] with i = (eulerphi(P) - mini) / phimult + 1 */

  for (j = 1, floor(maxi*2 / multP), 
    P = multP*(j+j+1);  
    phiP = eulerphi(P); 
    if (phiP % phimult != 0, 
      error("phiP = ", phiP, " is not a multiple of phimult = ", phimult)
    ); 
    if (mini <= phiP && phiP <= maxi, phi[(phiP - mini) / phimult + 1] = P));

  /* Go through the array and report large P, P/phi(P) combinations.
     best contains the maximal value of P * phiP seen so far */

  for (i = 1, len, 
    P = phi[i]; 
    phiP = (i - 1) * phimult + mini; 
    if (P * phiP > best * 1.05 && largest_primefactor(phiP) < 13, 
      print(P, " = ", factorint(P), ", ", phiP, " = ", factorint(phiP), " P/phi(P) = ", precision(1. * P / phiP, 10), " ",precision(1. * P * phiP / best, 10)); 
      best = P * phiP 
    );
  );

  return (best);
}