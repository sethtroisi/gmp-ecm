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

list_mul_opt := proc(n) option remember;
   if n<=1 then n
   elif member(n, {2,5,6,7,8,17,18,23,24,29,30}) then # Karatsuba
      2*procname(ceil(n/2))+procname(floor(n/2))
   elif member(n, {3,9,10,11,12,20,21,25,26,27}) then # toomcook3
      4*procname(ceil(n/3))+procname(n-2*ceil(n/3))
   else # toomcook4
      6*procname(ceil(n/4))+procname(n-3*ceil(n/4))
   fi
end:

# number of scalar multiplies from karatsuba
karatsuba := proc(n) option remember;
   if n<=1 then n
   else 2*procname(ceil(n/2))+procname(floor(n/2))
   fi
end:

# number of scalar multiplies from karatsuba, short product
karatsuba_short := proc(n) option remember;
   if n<=1 then n
   else procname(ceil(n/2))+2*procname(floor(n/2))
   fi
end:

# number of scalar multiplies from Karatsuba
# Mulders' short product (optimal cutoff)
karatsuba_short_mulders := proc(n) option remember; local p, m;
   if n<=1 then n
   else
      m := infinity;
      for p from ceil(n/2) to n-1 do
         m := min(m, karatsuba(p)+2*procname(n-p))
      od;
      m
   fi
end:

# get m terms, with entries of n terms
karatsuba_short2 := proc(n0, m0) option remember; local n, m;
   n := n0;
   m := m0;
   if m>2*n-1 then m:=2*n-1
   elif m<n then n:=m
   fi;
   if n<=1 then n
   else
      procname(ceil(n/2), ceil(m/2)) # evaluation at t=0
      + procname(ceil(n/2), ceil((m-1)/2)) # eval. at t=1
      + procname(ceil((n-1)/2), ceil((m-1)/2)) # t = inf
   fi
end:

# number of scalar multiplies from toomcook3, as implemented in ecm-5.0.1
toomcook3 := proc(n) option remember; local l, k;
   if member (n, {0,1,2,4}) then karatsuba(n)
   else
      l := iquo(n + 2, 3);
      k := n - 2*l;
      4*procname(l) + procname(k)
   fi
end:

# number of scalar multiplies for Toom3
# Mulders' short product (optimal cutoff)
toomcook3_short_mulders := proc(n) option remember; local p, m, c, s;
   if n<=1 then n
   else
      m := infinity;
      for p from ceil(n/2) to n do
         c := toomcook3(p)+2*procname(n-p);
         if c<m then m:=c; s:={p}
         elif c=m then s:=s union {p}
         fi
      od;
#      lprint(n, s);
      m
   fi
end:

# conjectured optimal cutoff: largest of 3^k or 2*3^k that is between n/2 and n
# seems to be true: works for all n<=1000
toomcook3_short_mulders2 := proc(n) option remember; local p;
   if n<=1 then n
   else
      p := 1;
      while 3*p<=n do p:=3*p od;
      p := floor(n/p)*p;
      toomcook3(p)+2*procname(n-p);
   fi
end:

# number of scalar multiplies from toomcook3, short product
# odd-even variant: short (a0(x^3) + x*a1(x^3)+ x^2*a2(x^3), n)
toomcook3_short := proc(n) option remember;
   if n <= 2 or n = 5 then karatsuba_short(n)
   else
      if n mod 3 = 2 then # consider (x*a(x))^2
         procname(iquo(n,3)) + 4*procname(iquo(n+2,3))
      else procname(ceil(n/3)) + 4*procname(ceil((n-1)/3))
      fi
   fi
end:

# get m terms, with entries of n terms
toomcook3_short2 := proc(n0, m0) option remember; local n, m, c1, c2;
   n := n0;
   m := m0;
   if m>2*n-1 then m:=2*n-1
   elif m<n then n:=m
   fi;
   if n<=1 then n
   elif member([n,m],{[2,2],[2,3],[4,5],[4,6],[4,7],[5,5]}) then karatsuba_short2(n,m)
   else
      c1 := toomcook3_short2(ceil(n/3), ceil(m/3)) # evaluation at t=0
      + 3*toomcook3_short2(ceil(n/3), ceil((m-1)/3)) # eval. at t=1, -1, 2
      + toomcook3_short2(ceil((n-2)/3), ceil((m-1)/3)); # t = inf
      # shift entries by x: n->n+1, m->m+2
      c2 := toomcook3_short2(ceil((n+1)/3)-1, ceil((m-4)/3)) # eval. at 0
      + 3*toomcook3_short2(ceil((n+1)/3), ceil((m+1)/3)) # eval. at t=1, -1, 2
      + toomcook3_short2(ceil((n-1)/3), ceil((m+1)/3)); # t = inf
      c1 := min(c1, c2);
      if karatsuba_short2(n,m) < c1 then
         lprint("karatsuba_short2 faster for ", n, m)
      fi;
      c1
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
   if member(n,{0,1,2,3,5,6,9,17,18,25,26,27,77,78,79,80,81}) then toomcook3(n)
   else
      l := iquo(n + 3, 4);
      k := n - 3 * l;
      6*procname(l) + procname(k)
   fi
end:

# find optimal method between kara, toom3 and toom4
find_opt := proc(nmax) local n, T, kara, toom3, toom4;
   T[0]:=0;
   T[1]:=1;
   for n from 2 to nmax do
      kara := 2*T[ceil(n/2)]+T[floor(n/2)];
      toom3 := 4*T[ceil(n/3)]+T[n-2*ceil(n/3)];
      if n>=4 and n<>5 then
         toom4 := 6*T[ceil(n/4)]+T[n-3*ceil(n/4)]
      else
         toom4 := kara;
      fi;
      if kara<=min(toom3,toom4) then lprint(n, "karatsuba", kara); T[n]:=kara
      elif toom3<=toom4 then lprint(n, "toomcook3", toom3); T[n]:=toom3
      else T[n]:=toom4
      fi
   od;
end:

# number of scalar multiplies for Toom3
# Mulders' short product (optimal cutoff)
toomcook4_short_mulders := proc(n) option remember; local p, m, c, s;
   if n<=1 then n
   else
      m := infinity;
      for p from ceil(n/2) to n do
         c := toomcook4(p)+2*procname(n-p);
         if c<m then m:=c; s:={p}
         elif c=m then s:=s union {p}
         fi
      od;
      lprint(n, s);
      m
   fi
end:

# conjectured optimal cutoff: 4^k or 2*4^k or 3*4^k that is between n/2 and n
# works for almost all n: exceptions are n=32 (p=27),
# n=125-130 (108)
toomcook4_short_mulders2 := proc(n) option remember; local p;
   if n<=1 then n
   else
      p := 1;
      while 4*p<=n do p:=4*p od;
      p:=floor(n/p)*p;
      toomcook4(p)+2*toomcook4_short_mulders(n-p);
   fi
end:

# number of scalar multiplies from toomcook4, short product
toomcook4_short := proc(n) option remember;
   if n <= 3 or member (n, {6,7,10,22,30,31,90,91,94}) then toomcook3_short(n)
   else
      min(procname(ceil(n/4)) + 6*procname(ceil((n-1)/4)),
          procname(ceil((n+1)/4)-1) + 6*procname(ceil(n/4)))
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
   elif k=l+1 then  toomcook4(l) + l # special important case
   else toomcook4(l) + list_mul(max(k-l, l), min(k-l, l))
   fi
end:

# number of multiplies of PolyFromRoots
PolyFromRoots := proc(k) option remember; local l, m;
   if k<=1 then 0
   elif k=2 then 1
   else
      m := iquo(k, 2);
      l := k - m; # l = k or l = k + 1
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

Reduce := proc(n) 2*toomcook4(n-1) + 1 end:

rootsF := proc(d, S)
   2.5 * 6 * (d/6) * S
end:

rootsG := proc(dF, S)
   2.5 * 6 * dF * S
end:

step2_cost := proc(B2, d, S) option remember; local dF, k, a, b;
   dF := numtheory[phi](d)/2;
   k := ceil (B2 / d / dF);
   a := PolyFromRoots(dF) + PolyInvert(dF-1) - toomcook4(dF) -
   Reduce(dF)
   + PolyEval(dF) + rootsF(d, S);
   b := PolyFromRoots(dF) # Building G from its roots
   + toomcook4(dF)   # Computing G * H
   + Reduce(dF) # Reducing G * H mod F
   + rootsG(dF, S);
   k, a, b, a + k * b;
end:

PolyInvert := proc(n) option remember; local k, l, v;
   if n=1 then 0
   else
      k := iquo(n, 2);
      l := n - k;
      v := procname(l) + toomcook4(l) + toomcook4(k);
      if k > 1 then v := v + list_mul (l-1, k-1) fi;
      v
   fi
end:

PolyEval := proc(k) option remember; local m, l, v;
   if k=1 then 0
   else
      m := iquo(k, 2);
      l := k - m;
      v := RecursiveDivision(m);
      if k > 2*m then v := v + m fi;
      v + RecursiveDivision(l) + procname(l) + procname(m)
   fi
end:

# output list with increasing phi(d) and decreasing step2_cost(d)
gen_bestD := proc(d0, d1) local l, d, c, p, i, j;
   l := [[d0,numtheory[phi](d0),step2_cost(d0)]];
   for d from d0+6 by 6 to d1 do
      p := numtheory[phi](d);
      c := step2_cost(d);
      for i to nops(l) while p > l[i][2] and c < l[i][3] do od;
      # now i > nops(l) or phi(d) <= phi(l[i]) or c >= step2_cost(l[i])
      if i > nops(l) then l:=[op(l),[d,p,c]]
      elif p <= l[i][2] then # p <= l[j][2] for j >= i
         for j from i to nops(l) while c <= l[j][3] do od;
         l:=[op(1..i-1, l), [d,p,c], op(j..nops(l), l)]
      else # p > l[i][2] and c >= step2_cost(l[i])
      fi
   od;
   l
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

############################ ecm ###################################

# converts (x:1:1) to Weierstrass form (mod p)
# returns [X,Y,A]
montgomery_to_weierstrass := proc(x, a, p) local g;
   g := (x^3 + a*x^2 + x) mod p;
   [(3*x+a)/(3*g) mod p, 1/g mod p, (3-a^2)/(3*g^2) mod p]
end:

addW := proc(x1, y1, x2, y2, n) local u, v, p;
   u := x2-x1;
   v := 1/u mod n; # 1/(x2-x1)
   p := (y2-y1)*v mod n; # lambda=(y2-y1)/(x2-x1)
   u := p*p-x1; # lambda^2-x1
   v := u-x2 mod n; # lambda^2-x1-x2
   u := (x1-v)*p; # (2x1+x2-lambda^2)*lambda
   [v, u-y1 mod n]
end:

# (x::y) -> 2*(x::y)
dupW := proc(x, y, n, a) local u, v, p;
   v := 1/(2*y) mod n;
   u := 3*x^2+a mod n;
   p := u*v mod n;
   u := p^2;
   v := 2*x;
   u := u-v mod n;
   [u, (x-u)*p-y mod n]
end:

# (x::y) -> k*(x::y)
mulW := proc(x, y, k, n, a) local l, P, i;
   if k=1 then [x,y]
   else # k >= 3
      l := convert(k, base, 2);
      P := [x, y];
      for i from nops(l)-1 by -1 to 1 do
         P := dupW(op(P), n, a);
         if l[i]=1 then P := addW(op(P), x, y, n) fi
      od;
      P
   fi
end:

##############################################################################

# odd-even variant
kara_short_mul := proc(a, b, n)
local a0, a1, b0, b1, c0, c1, c2, p, q, r, i, res;
   if n = 0 then []
   elif n = 1 then [a[1]*b[1]]
   else
      p := ceil(n/2);
      q := ceil((n-1)/2);
      r := q;
      a0 := [seq(a[2*i-1],i=1..p)];
      b0 := [seq(b[2*i-1],i=1..p)];
      a1 := [seq(a[2*i], i=1..q)];
      b1 := [seq(b[2*i], i=1..q)];
      c0 := procname(a0, b0, p);
      c1 := procname(a0[1..q]+a1, b0[1..q]+b1, q);
      c2 := procname(a1, b1, r);
      c1 := c1 - c0[1..q] - [op(c2),0$(q-r)];
      res := [0$n];
      for i to p do res[2*i-1]:=c0[i] od;
      for i to min(r,iquo(n-1,2)) do res[2*i+1]:=res[2*i+1]+c2[i] od;
      for i to q do res[2*i]:=c1[i] od;
      res
   fi
end:

# odd-even variant
toom3_short_mul := proc(a, b, n)
local a0, a1, a2, b0, b1, b2, c0, c1, c2, c3, c4, p, q, r, i, res;
   if n = 0 then []
   elif n = 1 then [a[1]*b[1]]
   elif n = 2 then [a[1]*b[1], a[1]*b[2]+a[2]*b[1]]
   else
      p := ceil(n/3);
      q := ceil((n-1)/3);
      r := ceil((n-2)/3);
      a0 := [seq(a[3*i-2],i=1..p)];
      b0 := [seq(b[3*i-2],i=1..p)];
      a1 := [seq(a[3*i-1], i=1..q)];
      b1 := [seq(b[3*i-1], i=1..q)];
      a2 := [seq(a[3*i], i=1..r), 0$(q-r)];
      b2 := [seq(b[3*i], i=1..r), 0$(q-r)];
      c0 := procname(a0, b0, p); # 0
      c1 := procname(a0[1..q]+a1+a2, b0[1..q]+b1+b2, q); # 1
      c2 := procname(a0[1..q]-a1+a2, b0[1..q]-b1+b2, q); # -1
      c3 := procname(a0[1..q]+2*a1+4*a2, b0[1..q]+2*b1+4*b2, q); # 2
      c4 := procname(a2, b2, q);
      c1 := c1 - c0[1..q] - c4;    # d1+d2+d3
      c2 := c2 - c0[1..q] - c4;    # -d1+d2-d3
      c3 := c3 - c0[1..q] - 16*c4; # 2*d1+4*d2+8*d3
      c1 := (c1 + c2)/2; # d2
      c2 := c2 - c1;     # -d1-d3
      c3 := c3 - 4*c1;   # 2*d1+8*d3
      c3 := c3 + 2*c2;   # 6*d3
      c3 := c3/6;        # d3
      c2 := -c2-c3;      # d1
      res := [0$n];
      for i to p do res[3*i-2]:=c0[i] od;
      for i to q do res[3*i-1]:=res[3*i-1]+c2[i] od;
      for i to ceil((n-2)/3) do res[3*i]:=res[3*i]+c1[i] od;
      for i to ceil((n-3)/3) do res[3*i+1]:=res[3*i+1]+c3[i] od;
      for i to ceil((n-4)/3) do res[3*i+2]:=res[3*i+2]+c4[i] od;
      res
   fi
end:
