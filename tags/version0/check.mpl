# stage2
stage2 := proc(p, d, a0, k) local x, F, inva, i, v, t, u, a, dF, G, j, H;
   a := a0;
   F := 1;
   dF := numtheory[phi](d) / 2;
   inva := 1/a mod p;
   for i from 1 to d by 6 do
      if igcd(i,d)=1 then
         v := Power(a,i) mod p + Power(inva, i) mod p;
         F := F * (x + v)
      fi
   od;
   lprint(F);
   F := Expand(F) mod p;
   lprint(F);

   a := Power(a, d);
   inva := Power(inva, d);
   t := 1;
   u := 1;

   for i to k do
      G := 1;
      for j from 1 to dF - 1 do
         t := t * a mod p;
         u := u * inva mod p;
         v := t + u mod p;
         G := G * (x + v);
      od;
      lprint("roots of G:", G);
      G := Expand(G) mod p;
      lprint("G:", G);
      if i = 1 then H := G
      else
	 lprint ("H:", H);
         G := Expand(G * H) mod p;
         lprint("G*H: ", G);
         H := Rem(G, F, x) mod p;
         lprint("G*H mod F:", H);
      fi
   od;
end:

# auxiliary memory needed for karatsuba
M := proc(K) option remember; local l;
   if K=1 then 0
   else
      l := iquo(K+1, 2);
      max(2*l-1+l,2*l-2+M(l))
   fi;
end: