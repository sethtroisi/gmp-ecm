
largest_primefactor(n) = vecmax(factorint(n)~[1,])

/* We examine P values with eulerphi(P) in [mini, maxi], 
   maxi = mini + multphi*(len - 1) */
/* We assume that P is divisible by multP and that 
   phiP and mini are divisible my multphi, so the array 
   entry phi[i] contains a P such that phiP = eulerphi(P) and 
   (phiP - mini) / multphi + 1 = i. (phi[i] contains zero if no such 
   P was found). Conversely, phiP = (i - 1) * multphi + mini */

/* mini = 150000000; */
/* len = 10000000; */
/* multP = 3*5*7*11; */

make_phiP (mini, len, multP, oldbest) = 
{
  local (multphi, maxi, phi, P, phiP, best, multP2);
  multphi = eulerphi(multP);
  best = oldbest;
  if (mini % multphi != 0, 
    error("mini = ", mini, " is not a multiple of multphi = ", multphi));
  maxi = mini + multphi * (len - 1);
  phi = vector(len);

  /* For each candidate odd P value < 4*maxi and a multiple of multP, 
     if eulerphi(P) is in the [mini, maxi] range, store the P value at 
     phi[i] with i = (eulerphi(P) - mini) / multphi + 1 */
  multP2 = 2*multP;
  P = 3*multP;

  while (P < maxi*4, 
    phiP = eulerphi(P); 
    if (phiP % multphi != 0, 
      error("phiP = ", phiP, " is not a multiple of multphi = ", multphi)
    ); 
    if (mini <= phiP && phiP <= maxi, phi[(phiP - mini) / multphi + 1] = P);
    P += multP2;
  );

  /* Go through the array and report large P, P/phi(P) combinations.
     best contains the maximal value of P * phiP seen so far */

  for (i = 1, len, 
    P = phi[i]; 
    phiP = (i - 1) * multphi + mini; 
    if (P * phiP > best * 1.05 && largest_primefactor(phiP) <= 19,
      print(P, " = ", factorint(P), ", ", phiP, " = ",  factorint(phiP), " P/phi(P) = ", 1. * P / phiP, " ", 1. * P * phiP / best); 
      best = P * phiP 
    );
  );

  return (best);
}
