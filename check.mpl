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

Powering["P-1"] := proc(a,i,p) Power(a,i) mod p end:
Select["P-1"] := proc(i) evalb(i mod 6=1) end:
Select["P+1"] := proc(i) member(i mod 6,{1,5}) end:

# stage2, method="P-1" or "P+1"
stage2 := proc(p, d, a0, k, method)
local x, F, inva, i, v, t, u, a, dF, G, j, H, ij;
   a := a0;
   F := 1;
   dF := numtheory[phi](d) / 2;
   if method = "P-1" then inva := 1/a mod p fi;
   for i from 1 to d by 2 do
      if Select[method](i) and igcd(i,d)=1 then
         v := Powering[method](a,i,p);
         if method="P-1" then v := v + Powering[method](inva,i,p) fi;
         F := F * (x + v)
      fi
   od;
   F := Expand(F) mod p;

   a := Powering[method](a, d, p);
   t := 1;
   if method="P-1" then
      inva := Powering[method](inva, d, p);
      u := 1;
   fi;

   lprint("B2=", k*(dF-1)*d);

   ij := 0;
   for i to k do
      G := 1;
      for j from 1 to dF - 1 do
         ij := ij + 1;
         t := Powering[method](a, ij, p);
         v := t;
         if method="P-1" then
            u := Powering[method](inva, ij, p);
            v := v + u;
         fi;
         G := G * (x + v);
      od;
      G := Expand(G) mod p;
      if i = 1 then H := G
      else
         G := Expand(G * H) mod p;
         H := Rem(G, F, x) mod p;
      fi;
   od;
   G := Gcd(F, H) mod p;
   if degree(G)<>0 then lprint("****** Found factor in stage 2: ", p) fi;
end:

# number of scalar multiplies from karatsuba
karatsuba := proc(n) option remember;
   if n<=1 then n
   else 2*procname(ceil(n/2))+procname(floor(n/2))
   fi
end:

# number of scalar multiplies from toomcook3
toomcook3 := proc(n) option remember; local l, k;
   if n <= 2 or n=4 then karatsuba(n)
   else
      l := iquo(n + 2, 3);
      k := n - 2*l;
      4*procname(l) + procname(k)
   fi
end:

toomcook3_2 := proc(n) option remember; local l0, l1, l2;
   if n<=2 or n=4 then karatsuba(n)
   else
      l2 := iquo(n + 2, 3);
      l1 := iquo(n + 1, 3);
      l0 := n - l2 - l1;
      3*procname(l2) + procname(l0) + procname(l1)
   fi
end:

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

# number of scalar multiplies from toomcook4
toomcook4 := proc(n) option remember; local l, k;
   if n<=3 or member(n,{5,6,9,17,18,26,27,77,78,79,80,81}) then toomcook3(n)
   else
      l := iquo(n + 3, 4);
      k := n - 3 * l;
      6*procname(l) + procname(k)
   fi
end:

toomcook4_2 := proc(n) option remember; local l2, l1, l0;
   if n<=3 or member(n,{5,6,9,17,18,25,26,27,77,78,79,80,81}) then toomcook3_2(n)
   else
      l2 := iquo(n + 3, 4);
      l1 := iquo(n - 2*l2 + 1, 2);
      if member(n,{22,30,41,42,45,54,57,70,73,82,85,86,94}) then l1:=l2 fi;
      l0 := n - 2*l2 - l1;
      5*procname(l2) + procname(l1) + procname(l0)
   fi
end:

short_mul := proc(n) option remember; local k;
   if member(n,{1,4,14,15,16,$56..64,$221..256}) then toomcook4(n)
   else
      k := 1; while 2*k < n do k:=2*k od;
      toomcook4(k) + 2*procname(n-k)
   fi
end:

# assumes k>=l
list_mul := proc(k, l)
   if k=l then toomcook4(l)
   else toomcook4(l) + list_mul(max(k-l, l), min(k-l, l))
   fi
end:

# number of multiplies of PolyFromRoots
PolyFromRoots := proc(k) option remember; local l, m;
   if k<=1 then 0
   elif k=2 then 1
   else
      m := iquo(k, 2);
      l := k - m;
      procname(l) + procname (m) + list_mul(l, m)
   fi
end:

# PolyFromRoots with optimal cutoff point
# (depends on the multiplication used, and the way list_mul deals with
#  operands of different degree)
PolyFromRoots_opt := proc(k) option remember; local l, m, s, c, cmin;
   if k<=1 then 0
   elif k=2 then 1
   else
      cmin := infinity;
      for m from 1 to iquo(k,2) do
         l := k - m;
         c := procname(l) + procname (m) + list_mul(l, m);
         if c<cmin then
            cmin := c;
            s := {l};
         elif c=cmin then s:=s union {l}
         fi
      od;
      lprint(k, s);
      cmin
   fi
end:

# number of multiplies in RecursiveDivision
RecursiveDivision := proc(K) option remember; local k, l;
   if K=1 then 1
   else
      k := iquo(K, 2);
      l := K - k;
      procname(l) + 2*toomcook4(k) + 2*k * (l-k) + procname(k)
   fi
end:

# estimate number of multiplies of PolyGcd(F, G)
# with deg(F)=n and deg(G)=n-1
PolyGcd := proc(n)
   if n<=1 then 0 # gcd is G
   else HalfGcd(n,0) + PolyGcd(ceil(n/2))
   fi
end:

# deg(F)=n and deg(G)=n-1, assumes return cofactors
# and matrix if flag=1
HalfGcd := proc(n,flag) option remember; local k, l, c;
   if n<=1 then 0
   else
      k := ceil(n/2);
      l := ceil(n/4);
      c := procname(k, 1) # return 2x2 matrix with degrees n/4
      + 8*toomcook4(l) # 4 multiplies of n/2 * n/4
      + procname(k, flag) # 2nd call
      + 4*toomcook4(l); # 4 multiplies of n/4 * n/4
      if flag=1 then c:=c+7*toomcook4(l) fi; # multiply matrices
      c
   fi
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
pp1_stage1 := proc(n, A0, B1) local a0, p, a, q;
   # suggested default choice from Montgomery
   if A0=0 then a0:=2/7 mod n else a0:=A0 fi;
   if igcd(a0^2-4, n)<>1 then ERROR("igcd(a0^2-4, n)<>1") fi;
   if isprime(n) and numtheory[jacobi](a0^2-4, n)=1 then
      lprint("Warning: jacobi(a0^2-4, n) = 1")
   fi;
   p := 2;
   a := a0;
   while p <= B1 do
      q := p;
      while q*p <= B1 do q:=q*p od;
      a := Powering["P+1"](a, q, n);
      p := nextprime(p);
   od;
   p:=igcd(a-2, n);
   if p<>1 then lprint("****** Found factor in stage 1: ", p) fi;
   a
end:

Powering["P+1"] := proc(P0, p, n) local l, i, P, Q, R;
   l := convert(p-1, base, 2);
   P := P0;
   Q := 2;
   for i from nops(l) by -1 to 1 do
      if l[i] = 1 then # (i,i-1) to (2i,2i-1)
         Q := P*Q-P0 mod n;
         P := P^2-2 mod n;
      else # (i,i-1) to (2i-1,2i-2)
         P := P*Q-P0 mod n;
         Q := Q^2-2 mod n;
      fi
   od;
   P
end:
