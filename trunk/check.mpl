# memory used by toomcook3
M := proc(len) local l;
option remember;
   l := iquo(len + 2, 3);
   4 * l + max (M(l), 1)
end:
M(0):=0:
M(1):=0:
M(2):=1:
M(4):=5:

pm1_stage1 := proc(n, a0, B1)
local p, a, q;
   p := 2;
   a := a0;
   while p <= B1 do
      q := 1;
      while q*p <= B1 do
         q := q*p;
         a := Power(a, p) mod n;
      od;
      p := nextprime(p);
   od;
   a
end:

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

# cf Williams, Math. of Comp. 39, 1982
# pp1_stage1(328006342451, 5, 7043); # 2^235+1
# pp1_stage1(6215074747201, 5, 199729); # 2^297+1
# pp1_stage1(8857714771093, 3, 49207); # 2^418+1
# pp1_stage1(236344687097, 3, 55001); # 2^602+1
# pp1_stage1(87251820842149, 5, 170249); # 2^642+1 
# pp1_stage1(719571227339189, 4, 57679); # 3^134+1 
# pp1_stage1(5468575720021, 6, 175759); # 3^161+1
# pp1_stage1(49804972211, 5, 268757); # 3^185-1
# pp1_stage1(329573417220613, 3, 101573); # 5^94+1
# pp1_stage1(4866979762781, 4, 97609); # 6^59-1
# pp1_stage1(187333846633, 3, 9851); # 6^111-1
# pp1_stage1(332526664667473, 3, 111919); # 6^132+1
# pp1_stage1(265043186297, 3, 152791); # 7^133+1
# pp1_stage1(207734163253, 3, 4211); # 7^231+1
# pp1_stage1(225974065503889, 5, 8243); # 10^102+1
# pp1_stage1(660198074631409, 5, 115679); # 12^81-1
# pp1_stage1(563215815517, 3, 109849); # 12^183+1
# pp1_stage1(409100738617, 3, 70957); # fib(247)
# pp1_stage1(7901346123803597, 3, 18307); # fib(313)
# pp1_stage1(5648966761, 15, 100519); # fib(367)
# pp1_stage1(14279673833, 3, 823); # fib(387)
# pp1_stage1(1795220677069, 6, 159931); # fib(483)
# pp1_stage1(1250839826281, 5, 4673); # fib(495)
# pp1_stage1(2192843129417, 3, 38803); # fib(531)
# pp1_stage1(10424083697, 3, 131); # fib(549)
pp1_stage1 := proc(n, a0, B1) local p, a, q;
   if igcd(a0^2-4, n)<>1 then ERROR("igcd(a0^2-4, n)<>1") fi;
   if isprime(n) and numtheory[jacobi](a0^2-4, n)=1 then
      lprint("Warning: jacobi(a0^2-4, n) = 1")
   fi;
   p := 2;
   a := a0;
   while p <= B1 do
      q := p;
      while q*p <= B1 do q:=q*p od;
      a := PowerPP1(a, q, n);
      p := nextprime(p);
   od;
   igcd(a-2, n)
end:

PowerPP1 := proc(P0, p, n) local l, i, P, Q, R;
   i := p;
   l := NULL;
   while i > 1 do l:=i,l; i:=iquo(i+1, 2) od;
   P := P0;
   Q := 2;
   for i in [l] do
      # go from ceil(i/2) to i
      if i mod 2 = 0 then # (i,i-1) to (2i,2i-1)
         Q := P*Q-P0 mod n;
         P := P^2-2 mod n;
      else # (i,i-1) to (2i-1,2i-2)
         P := P*Q-P0 mod n;
         Q := Q^2-2 mod n;
      fi
   od;
   P
end:
