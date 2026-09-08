#!/usr/bin/env python3
"""Complete exact Python verification for the detailed Fujita note.

Run: python verify_fujita_further_extensions.py
Requires Python 3.10 or newer and SymPy.
Install SymPy with: python -m pip install sympy
No compiler, generated executable, NumPy, or other source file is used.
All mathematical calculations use exact integer or rational arithmetic. The finite sieves use Python integers as sets of allowed third values.
Use --export-json PATH to save coefficients and audit data for the TeX note.
"""
from math import factorial, comb, isqrt
import argparse
import json
import time
from pathlib import Path
import sympy as S
AUDIT = {}

def multiply(a, b, n):
    return [S.expand(sum(a[i] * b[j-i]
                         for i in range(max(0, j-len(b)+1), min(j+1, len(a)))))
            for j in range(n+1)]

def logarithm(a, n):
    assert a[0] == 1
    result = [S.Integer(0)] * (n+1)
    for j in range(1, n+1):
        result[j] = S.expand(a[j] - sum(S.Rational(k, j)*result[k]*a[j-k]
                                       for k in range(1, j)))
    return result

def exponential(a, n):
    assert a[0] == 0
    result = [S.Integer(1)]
    for j in range(1, n+1):
        result.append(S.expand(sum(k*a[k]*result[j-k]
                                   for k in range(1, j+1))/j))
    return result

def todd_from_powers(ps):
    n = len(ps)-1
    logarithm_td = [S.Integer(0)] * (n+1)
    logarithm_td[1] = ps[1]/2
    for j in range(2, n+1, 2):
        logarithm_td[j] = -S.bernoulli(j)*ps[j]/(j*S.factorial(j))
    return exponential(logarithm_td, n)

def chern_from_powers(ps):
    n = len(ps)-1
    cs = [S.Integer(1)]
    for k in range(1, n+1):
        cs.append(S.expand(sum((-1)**(j-1)*ps[j]*cs[k-j]
                               for j in range(1, k+1))/k))
    return cs

def exterior_euler_characteristics(ps, td, maximum):
    n = len(ps)-1
    # ch(Omega^1) = sum_i exp(-alpha_i).
    ch1 = [(-1)**j*ps[j]/S.factorial(j) for j in range(n+1)]
    exterior = [[S.Integer(1)]+[S.Integer(0)]*n]
    chis = [td[n]]
    for p in range(1, maximum+1):
        # Newton recurrence for exterior powers; Adams_a multiplies ch_j by a^j.
        cp = [S.expand(sum((-1)**(a-1)
                          *sum(exterior[p-a][j-k]*a**k*ch1[k]
                               for k in range(j+1))
                          for a in range(1, p+1))/p)
              for j in range(n+1)]
        exterior.append(cp)
        chis.append(S.expand(sum(cp[j]*td[n-j] for j in range(n+1))))
    return chis

def coindex4_todd_and_even_powers(n, m):
    r = n-3
    t = S.Symbol('t')
    z = t*(t+r)
    a = (r+2)*(r+3)*(m-r)-(r+1)
    b = r*(r+1)*(r+2)*(r+3)
    P = S.Poly(S.prod(t+j for j in range(1, r))*(z*z+a*z+b)/S.factorial(n), t)
    assert P.eval(0) == 1 and P.eval(1) == m
    assert all(P.eval(-j) == 0 for j in range(1, r))
    assert S.expand(P.as_expr().subs(t, -r-t)-(-1)**n*P.as_expr()) == 0
    td = [S.factorial(n-j)*P.coeff_monomial(t**(n-j)) for j in range(n+1)]
    ahat = multiply(td, [(-S.Rational(r, 2))**j/S.factorial(j)
                         for j in range(n+1)], n)
    assert all(ahat[j] == 0 for j in range(1, n+1, 2))
    la = logarithm(ahat, n)
    ps = [S.Integer(n)] + [S.Integer(0)]*n
    ps[1] = S.Integer(r)
    for j in range(2, n+1, 2):
        ps[j] = S.cancel(la[j]/(-S.bernoulli(j)/(j*S.factorial(j))))

    # Independent derivation using the reciprocal Hilbert numerator.
    delta = n+1-m
    B = [S.Integer(1)] + [S.Integer(0)]*n
    for j in range(2, n+1, 2):
        B[j] = 2*S.Rational(2**j-delta, factorial(j))
    lb = logarithm(B, n)
    for j in range(2, n+1, 2):
        assert ps[j] == n+1-j*S.factorial(j)*lb[j]/S.bernoulli(j)
    return td, ps

def derive_equations(n, m):
    r = n-3
    original_td, ps = coindex4_todd_and_even_powers(n, m)
    odd = list(S.symbols(' '.join('s'+str(j) for j in range(3, n+1, 2))))
    for j, v in zip(range(3, n+1, 2), odd):
        ps[j] = v
    cs = chern_from_powers(ps)
    target = S.Rational(n*(n+1)**2, 2*r)
    high1, high2 = odd[-2:]
    sol1 = S.solve(cs[n-1]-target, high1)[0]
    sol2 = S.solve((cs[n]-(n+1)).subs(high1, sol1), high2)[0]
    substitutions = {high1: sol1, high2: sol2}
    ps = [S.expand(v.subs(substitutions)) for v in ps]
    cs_check = chern_from_powers(ps)
    assert S.expand(cs_check[n-1]-target) == 0
    assert S.expand(cs_check[n]-(n+1)) == 0
    td = todd_from_powers(ps)
    assert td == original_td
    variables = odd[:-2]
    chis = exterior_euler_characteristics(ps, td, len(variables))
    equations = [S.Poly(chis[p]-(-1)**p, *variables).clear_denoms()[1].as_expr()
                 for p in range(1, len(variables)+1)]
    return variables, equations

def projective_space_controls():
    for n in [7, 9, 11]:
        ps = [S.Integer(n)]+[S.Integer(n+1)]*n
        td = todd_from_powers(ps)
        cs = chern_from_powers(ps)
        assert cs == [S.binomial(n+1, j) for j in range(n+1)]
        chis = exterior_euler_characteristics(ps, td, (n-1)//2)
        assert chis == [S.Integer((-1)**p) for p in range(len(chis))]
    print('CP7, CP9, and CP11 positive controls passed.', flush=True)

def signature_polynomial(n, r):
    """Return the signature-minus-one polynomial, derived by residues."""
    assert n % 2 == 0 and (n+1-r) % 2 == 0
    s = (n+1-r)//2
    variables = list(S.symbols(' '.join('h'+str(i) for i in range(1, s))))
    values = [S.Integer(1)] + variables
    h = []
    for i in range(s):
        h.append(S.expand(values[i] - sum(h[j]*S.binomial(n+i-j, n)
                                         for j in range(i))))
    h = h + [1-2*sum(h)] + list(reversed(h))
    t = S.Symbol('t')
    d = len(h)-1
    N = S.Poly(S.expand(sum(h[j]*(1-t)**j*(1+t)**(d-j)
                           for j in range(d+1))), t)
    D = S.Poly(S.expand((1-t*t)*sum(h[j]*(1-t)**(2*j)*(1+t)**(2*d-2*j)
                                   for j in range(d+1))), t)
    assert N.nth(0) == D.nth(0) == 1
    NN = N*N
    a = []
    for i in range(n//2+1):
        a.append(S.expand(NN.nth(2*i)-sum(D.nth(2*j)*a[i-j]
                                        for j in range(1, i+1))))
    raw = S.Poly(a[-1]-1, *variables)
    return variables, h, raw.primitive()[1], raw

def check_signature_reconstruction(n, r, variables, h, raw):
    """Independent Hilbert-polynomial -> Todd -> power sums -> L check."""
    t = S.Symbol('t')
    universal_L = S.series(t/S.tanh(t), t, 0, n+1).removeO().expand()
    logL = logarithm([universal_L.coeff(t, j) for j in range(n+1)], n)
    for vals in ([0]*len(variables), list(range(1, len(variables)+1))):
        subs = dict(zip(variables, vals))
        hs = [S.sympify(e).subs(subs) for e in h]
        P = S.Poly(S.expand(sum(hs[q]*S.prod(t+i-q for i in range(1, n+1))
                               /S.factorial(n) for q in range(len(hs)))), t)
        assert [P.eval(i+1) for i in range(len(vals))] == vals
        td = [S.factorial(n-j)*P.nth(n-j) for j in range(n+1)]
        ltd = logarithm(td, n)
        powers = [S.Integer(0)]*(n+1)
        for j in range(2, n+1, 2):
            powers[j] = -ltd[j]*j*S.factorial(j)/S.bernoulli(j)
        actual = exponential([logL[j]*powers[j] for j in range(n+1)], n)[n]
        assert actual == raw.as_expr().subs(subs)+1

def c9_eliminant():
    u,v,w,a,b,c,d,z = S.symbols('u v w a b c d z')
    ps = [S.Integer(9),S.Integer(2),u,a,v,b,w,c,z,d]
    zsol = S.solve(todd_from_powers(ps)[9]-1, z)[0]
    ps = [S.expand(e.subs(z,zsol)) for e in ps]
    cs = chern_from_powers(ps)
    csol = S.solve(cs[8]-225,c)[0]
    dsol = S.solve((cs[9]-10).subs(c,csol),d)[0]
    ps = [S.expand(e.subs({c:csol,d:dsol})) for e in ps]
    cs = chern_from_powers(ps)
    assert cs[8] == 225 and cs[9] == 10
    chi = exterior_euler_characteristics(ps,todd_from_powers(ps),2)
    equations = [S.Poly(chi[i]-(-1)**i,u,v,w,a,b).clear_denoms()[1]
                 .primitive()[1].as_expr() for i in (1,2)]
    assert all(S.degree(e,b)==1 for e in equations)
    R = S.Poly(S.resultant(*equations,b),u,v,w,a)
    F = R.primitive()[1]
    AUDIT['c9_elimination'] = {
        's8_latex': S.latex(zsol),
        'equations': [polynomial_record(S.Poly(E,u,v,w,a,b),b) for E in equations],
        'resultant_content': str(R.content()),
    }
    assert S.degree(F.as_expr(),a)==4 and len(F.terms())==60
    # Independently recover the even power sums from four Hilbert polynomials.
    t = S.Symbol('t')
    for m,k,l in [(0,0,0),(1,2,3),(9,521,19692),(3,217,9421)]:
        h = [1,m-10,k-10*m+45,l-10*k+45*m-120]
        h = h+[1-2*sum(h)]+list(reversed(h))
        P = S.Poly(S.expand(sum(h[j]*S.prod(t+i-j for i in range(1,10))
                               /S.factorial(9) for j in range(9))),t)
        assert [P.eval(i) for i in range(4)] == [1,m,k,l]
        td = [S.factorial(9-j)*P.nth(9-j) for j in range(10)]
        ltd = logarithm(td,9)
        T,Q,H = l-6*k+14*m-14,l+6*k-34*m+46,l+54*k+134*m-434
        subs = {u:10-24*T,
                v:10+240*Q-1440*T*T,
                w:10-504*H+15120*T*Q-60480*T*T*T}
        assert [subs[u],subs[v],subs[w]] == [-24*ltd[2],2880*ltd[4],-181440*ltd[6]]
        assert [e.subs(subs) for e in todd_from_powers(ps)] == td
    print('C9 index 2: quartic derived; independent Hilbert checks passed.',flush=True)
    return [u,v,w],a,F

def coindex6_even():
    # For each m=0,...,n the listed prime excludes EVERY integer k.
    certificates={
      6:[5,3,3,5,3,3,7],
      8:[11,11,7,7,11,7,7,19,11],
      10:[3,3,7,3,3,7,3,3,13,3,3],
      14:[11,5,5,11,7,7,5,5,7,23,23,5,5,23,11],
      20:[19,3,11,7,3,13,13,3,17,17,3,17,11,3,17,11,3,7,13,3,17],
      50: [
          11,5,11,47,7,7,5,17,7,61,11,5,7,11,23,
          7,5,47,7,7,29,5,7,23,11,7,5,11,17,7,
          17,5,7,7,13,11,5,19,11,7,7,5,37,7,11,
          17,5,7,19,11,7,
      ]
    }
    A,B,t=S.symbols('A B t')
    h=[1,A,B,-1-2*A-2*B,B,A,1]
    N=S.Poly(S.expand(sum(h[j]*(1-t)**j*(1+t)**(6-j) for j in range(7))),t)
    D=S.Poly(S.expand((1-t*t)*sum(h[j]*(1-t)**(2*j)*(1+t)**(12-2*j)
                                for j in range(7))),t)
    # Coefficients in u=t^2, each as constant + coefficient*A + coefficient*B.
    coeff=lambda E,count:[[int(S.Poly(E.nth(2*j),A,B).nth(*ex))
                           for ex in [(0,0),(1,0),(0,1)]] for j in range(count)]
    nc,dc=coeff(N,4),coeff(D,8)
    assert [r+5 for r in S.divisors(90) if r%2] == list(certificates)
    for n,primes in certificates.items():
        assert len(primes)==n+1
        for m,p in enumerate(primes):
            assert S.isprime(p)
            for k in range(p):
                aa=m-n-1;bb=k-(n+1)*m+n*(n+1)//2
                nn=[(c+a*aa+b*bb)%p for c,a,b in nc]
                dd=[(c+a*aa+b*bb)%p for c,a,b in dc]
                numerator=[sum(nn[i]*nn[j-i] for i in range(max(0,j-3),min(j+1,4)))%p
                           for j in range(7)]
                out=[]
                for j in range(n//2+1):
                    value=numerator[j] if j<7 else 0
                    value-=sum(dd[q]*out[j-q] for q in range(1,min(j,7)+1))
                    out.append(value%p)
                assert out[-1]!=1,(n,m,k,p)
        print(f'n={n}, r={n-5}: every h0(L)=0,...,{n} excluded.',flush=True)

def c7_coefficients(m, k):
    u = k - 4*m
    A = 18*u + 101
    B = 1728*u**2 + 16392*u + 47544 - 4320*m
    C = (528768*u**4 + 10284192*u**3 + 77812128*u**2
         + 271230552*u + 365944288
         - m*(1399680*u**2 + 13478400*u + 32447520))
    return A, B, C

def c7_finite_certificate():
    total = nonnegative = squares = 0
    per_m = []
    candidates = []
    for m in range(9):
        count = 0
        for k in range(136):
            A, B, C = c7_coefficients(m, k)
            assert A != 0
            D = B*B - 4*A*C
            total += 1
            if D >= 0:
                nonnegative += 1
                count += 1
                candidates.append((m, k, D))
                assert isqrt(D)**2 != D, (m, k, D)
        per_m.append(count)
    assert total == 1224 and nonnegative == 114
    assert per_m == [0, 0, 4, 7, 12, 16, 20, 25, 30]
    print('Index 2: 1224 pairs; 114 nonnegative discriminants; 0 squares.')
    print('Nonnegative counts for m=0,...,8:', per_m)

    # An alternative certificate using only quadratic residues.
    expected = [90, 60, 33, 23, 16, 8, 3, 0]
    for prime, count in zip([5, 7, 11, 13, 17, 19, 23, 29], expected):
        residues = {j*j % prime for j in range(prime)}
        candidates = [q for q in candidates if q[2] % prime in residues]
        assert len(candidates) == count
        print('After square-residue test modulo', prime, ':', count)
    assert not candidates

    # Index 4: the reduced equation is 4*b^2+3*b*m+6 == 0 mod 9.
    assert all((4*b*b+3*b*m+6) % 9 != 0
               for b in range(9) for m in range(9))
    print('Index 4: no solution of 4*b^2+3*b*m+6 == 0 mod 9.')

def c7_symbolic_certificate():
    import sympy as s
    r,a,b,c,d,e,f,t,m,k,u = s.symbols('r a b c d e f t m k u')
    cs = [1,r,a,b,c,d,e,f]

    # Newton identities: p_j is the j-th power sum of the seven Chern roots.
    ps = [s.Integer(7)]
    for j in range(1,8):
        ps.append(s.expand(sum((-1)**(i-1)*cs[i]*ps[j-i]
                               for i in range(1,j))
                           + (-1)**(j-1)*j*cs[j]))

    # log Todd = p1*z/2 - p2*z^2/24 + p4*z^4/2880 - p6*z^6/181440.
    logtd = {1:ps[1]/2, 2:-ps[2]/24, 4:ps[4]/2880, 6:-ps[6]/181440}
    td = [s.Integer(1)]
    for j in range(1,8):
        td.append(s.expand(sum(i*v*td[j-i]
                               for i,v in logtd.items() if i<=j)/j))
    P = s.expand(sum(t**j/s.factorial(j)*td[7-j] for j in range(8)))

    ch1 = [(-1)**j*ps[j]/s.factorial(j) for j in range(8)]
    ch2 = [s.expand((sum(ch1[i]*ch1[j-i] for i in range(j+1))
                     - 2**j*ch1[j])/2) for j in range(8)]
    chi1 = s.expand(sum(ch1[j]*td[7-j] for j in range(8)))
    chi2 = s.expand(sum(ch2[j]*td[7-j] for j in range(8)))

    E1 = r**3*c-3*r*a*c-r**2*d+3*b*c-3*a*d+7728
    E2 = (-2*r**5*a+10*r**3*a**2+2*r**4*b-10*r*a**3
          -11*r**2*a*b-2*r**3*c+r*b**2+9*r*a*c+2*r**2*d+120512)
    Ps = (12*t**7+42*t**6*r+42*t**5*(r**2+a)+105*t**4*r*a
          +14*t**3*(-r**4+4*r**2*a+3*a**2+r*b-c)
          +21*t**2*(-r**3*a+3*r*a**2+r**2*b-r*c)
          +t*(2*r**6-12*r**4*a+11*r**2*a**2+5*r**3*b+10*a**3
              +11*r*a*b-5*r**2*c-b**2-9*a*c-2*r*d+2*e)
          +60480)/60480
    assert s.expand(P-Ps-td[7]+1) == 0

    # Independently regenerate the two Chern identities from HRR.
    for index in [2,4,8]:
        ss = {r:index, e:224//index, f:8}
        assert s.expand(E2.subs(ss)+120960*(td[7].subs(ss)-1)) == 0
        assert s.expand(E1.subs(ss)+7200*(td[7].subs(ss)-1)
                        -1440*(chi1.subs(ss)+1)) == 0
        assert s.expand((chi2.subs(ss)-1)+6*(td[7].subs(ss)-1)
                        -3*(chi1.subs(ss)+1)) == 0
        G = s.groebner([td[7].subs(ss)-1, chi1.subs(ss)+1,
                       chi2.subs(ss)-1], d,c,b,a, domain=s.QQ)
        assert G.reduce(E1.subs(ss))[1] == 0
        assert G.reduce(E2.subs(ss))[1] == 0
        H = s.groebner([E1.subs(ss), E2.subs(ss)], d,c,b,a, domain=s.QQ)
        for v in [td[7]-1, chi1+1, chi2-1]:
            assert H.reduce(v.subs(ss))[1] == 0
    print('Todd, Omega^1, Omega^2 computations independently match the Chern identities.')

    standard = {cs[j]:s.binomial(8,j) for j in range(1,8)}
    assert s.expand(P.subs(standard)
                    - s.prod(t+j for j in range(1,8))/s.factorial(7)) == 0
    print('CP7 Hilbert polynomial check passed.')

    F4 = (9*b**2*m-5*b**2+1728*b*m**2-17160*b*m+36612*b
          +264384*m**4-4279824*m**3+25277544*m**2-64549008*m+60164376)
    ss = {r:4,e:56,f:8}
    eq4 = [E1.subs(ss),E2.subs(ss),Ps.subs(ss).subs(t,-1),
           Ps.subs(ss).subs(t,1)-m]
    G4 = s.groebner(eq4,d,c,b,a,m,domain=s.QQ)
    for consequence in [a-12*m+44, F4]:
        assert G4.reduce(consequence)[1] == 0
    remainder = s.Poly(F4-(4*b**2+3*b*m+6),b,m)
    assert all(int(z)%9 == 0 for z in remainder.coeffs())
    print('Index 4 elimination and reduction modulo 9 verified.')

    A,B,C = c7_coefficients(m,k)
    F2 = s.expand(A*b*b+B*b+C)
    ss = {r:2,e:112,f:8}
    eq2 = [E1.subs(ss),E2.subs(ss),Ps.subs(ss).subs(t,1)-m,
           Ps.subs(ss).subs(t,2)-k]
    G2 = s.groebner(eq2,d,c,b,a,m,k,domain=s.QQ)
    for consequence in [a-12*k+48*m-58, F2]:
        assert G2.reduce(consequence)[1] == 0
    print('Index 2 elimination verified; exactly the polynomial used by the finite certificate.')

def c8_F3(m, u):
    return (1152*u**4 + 40920*u**3 + (557147-4032*m)*u**2
            + (3443306-71232*m)*u + 1800*m**2 - 325410*m + 8148917)

def c8_F1(m, k, u):
    v = k - 3*m
    return (1152*u**4 - 23592*u**3 + 173099*u**2 - 537930*u
            + 596995 + 1800*v**2
            + (-4032*u**2 + 41664*u - 100434)*v
            + (3024*u - 16000)*m)

def c8_arithmetic_certificate():
    # No restriction at all on the integer u is needed.
    moduli3 = [11, 11, 7, 7, 11, 7, 7, 19, 11, 7]
    for m, q in enumerate(moduli3):
        assert all(c8_F3(m, u) % q != 0 for u in range(q)), (m, q)
    print('Index 3: for every m=0,...,9, a modulus excludes every integer u.')

    # m=h^0(L) is in [0,9], k=h^0(2L) is in [0,264].
    pairs = [(m, k) for m in range(10) for k in range(265)]
    moduli1 = [4, 9, 16, 25, 27, 49, 11, 13, 17, 19, 23, 29, 31, 37, 43]
    expected = [1325, 662, 332, 160, 148, 80, 64, 39, 25, 14, 8, 3, 2, 1, 0]
    print('Index 1: initial pairs:', len(pairs))
    for q, count in zip(moduli1, expected):
        pairs = [(m,k) for m,k in pairs
                 if any(c8_F1(m,k,u) % q == 0 for u in range(q))]
        assert len(pairs) == count, (q, len(pairs), count)
        print('After modulus', q, ':', len(pairs))
    assert not pairs
    print('Both nonmaximal indices are excluded.')

def c8_symbolic_certificate():
    import sympy as s
    from math import factorial
    t,z,w,A,B,C,m,k,l,u = s.symbols('t z w A B C m k l u')
    base = (z**4+A*z**3+B*z**2+C*z+40320)/40320

    # Derive the universal conversion from log A-hat to log L.
    log_ahat = s.series(s.log(w/(2*s.sinh(w/2))),w,0,10).removeO()
    log_L = s.series(s.log(w/s.tanh(w)),w,0,10).removeO()
    weights = [s.cancel(log_L.coeff(w,2*j)/log_ahat.coeff(w,2*j))
               for j in range(1,5)]
    assert weights == [-8, -224, -3968, -65024]
    print('Universal signature conversion independently derived.')

    expected_abc = {
        3: {A:56*k-280*m+492, B:-112*k+2240*m-6036,
            C:-448*k+5600*m+6128},
        1: {A:-280*k+56*l+504*m-300, B:3920*k-448*l-9072*m+5708,
            C:-6720*k+672*l+36288*m-30384},
    }

    for index in [3,1,9]:
        if index == 9:
            P = s.Poly(s.prod(t+j for j in range(1,9))/factorial(8),t)
        else:
            equations = [base.subs(z,1+index)-m, base.subs(z,4+2*index)-k]
            equations += ([base.subs(z,-2)] if index==3
                          else [base.subs(z,9+3*index)-l])
            abc = s.solve(equations,[A,B,C])
            assert all(s.expand(abc[v]-expected_abc[index][v])==0
                       for v in [A,B,C])
            P = s.Poly(s.expand(base.subs(abc).subs(z,t*(t+index))),t)
            assert s.expand(P.as_expr().subs(t,-index-t)-P.as_expr())==0
            assert P.eval(0)==1
            assert s.expand(P.eval(1)-m)==0
            assert s.expand(P.eval(2)-k)==0
            if index==3:
                assert P.eval(-1)==P.eval(-2)==0
            else:
                assert s.expand(P.eval(3)-l)==0

        # HRR: td_j=(8-j)! [t^(8-j)] P(t).
        td = [factorial(8-j)*P.coeff_monomial(t**(8-j)) for j in range(9)]
        ahat = [s.expand(sum((-s.Rational(index,2))**i/s.factorial(i)*td[j-i]
                            for i in range(j+1))) for j in range(9)]
        assert ahat[0]==1
        assert all(ahat[j]==0 for j in [1,3,5,7])
        U,V,W,Z = [ahat[j] for j in [2,4,6,8]]
        q = [U, V-U**2/2, W-U*V+U**3/3,
             Z-U*W-V**2/2+U**2*V-U**4/4]
        l1,l2,l3,l4 = [weights[j]*q[j] for j in range(4)]
        signature = s.expand(l4+l1*l3+l2**2/2+l1**2*l2/2+l1**4/24)

        if index==3:
            assert s.expand((signature-1).subs(k,u+5*m)-32*c8_F3(m,u))==0
        elif index==1:
            assert s.expand((signature-1).subs(l,u+5*k-9*m)-32*c8_F1(m,k,u))==0
        else:
            assert signature==1
            assert td[1]==s.Rational(9,2)
        print('Index',index,': Hilbert polynomial and signature identity verified.')

def high_bernoulli_certificate():
    delta = S.Symbol('delta')
    B = [S.Integer(1)] + [S.Integer(0)]*12
    for j in range(2, 13, 2):
        B[j] = 2*(S.Integer(2)**j-delta)/S.factorial(j)
    f = S.Poly(S.factorial(12)*logarithm(B, 12)[12], delta)
    assert all(c.q == 1 for c in f.all_coeffs())
    expected = -297*delta*(delta-1)*(delta-2)*(delta**3+325*delta**2-326*delta+237)
    assert all(c % 691 == 0 for c in S.Poly(f.as_expr()-expected, delta).all_coeffs())
    assert S.bernoulli(12) == -S.Rational(691, 2730)
    assert S.gcd(32760, 691) == 1
    assert [a for a in range(691) if f.eval(a) % 691 == 0] == [0, 1, 2]
    print('Bernoulli obstruction: delta is 0, 1, or 2 modulo 691.', flush=True)

def high_dimension_and_euler_certificates():
    t = S.Symbol('t')
    indices = S.divisors(24)
    assert [r+3 for r in indices] == [4, 5, 6, 7, 9, 11, 15, 27]
    # Rational generating functions for the Euler differences of the two
    # weighted complete intersections. Their coefficients are given in the note.
    g1 = (1+t)*(1+4*t)/((1-t)**2*(1+9*t))-1/(1-t)**2
    assert S.cancel(g1-S.Rational(2, 5)*(1/(1+9*t)-1/(1-t))) == 0
    g2 = (1+t)**2*(1+2*t)**2/((1-t)**2*(1+5*t)**2)-1/(1-t)**2
    g2_explicit = (S.Rational(4, 25)+S.Rational(26, 75)/(1+5*t)
                   +S.Rational(4, 25)/(1+5*t)**2-S.Rational(2, 3)/(1-t))
    assert S.cancel(g2-g2_explicit) == 0
    for n in [15, 27]:
        difference1 = S.Rational(2, 5)*((-9)**n-1)
        difference2 = S.Rational((12*n+38)*(-5)**n-50, 75)
        assert difference1.q == difference2.q == 1
        assert difference1 < 0 and difference2 < 0
    print('Dimensions 15 and 27: both weighted models have Euler number '
          'different from n+1.', flush=True)

def polynomial_record(poly, root=None):
    record = {'variables': [str(v) for v in poly.gens],
              'terms': [[list(e), int(c)] for e, c in poly.terms()]}
    if root is not None:
        p = S.Poly(poly.as_expr(), root)
        record['root'] = str(root)
        record['coefficients_latex'] = [S.latex(p.nth(i)) for i in range(p.degree()+1)]
    return record


def integer_has_root(coefficients, modulus):
    """Horner evaluation, including the zero polynomial and composite moduli."""
    for residue in range(modulus):
        value = 0
        for coefficient in reversed(coefficients):
            value = (value*residue+coefficient) % modulus
        if value == 0:
            return True
    return False


class PolynomialRootTest:
    """Specialize a polynomial in three parameters, then test its last variable.

    All coefficients are Python integers. Reducing a parameter before evaluation
    gives the same result because the original polynomial has integer coefficients.
    """
    def __init__(self, poly, parameters, root, transform=None):
        f = S.Poly(poly.as_expr(), root)
        self.terms = [[(int(c), e) for e, c in S.Poly(f.nth(j), *parameters).terms()]
                      for j in range(f.degree()+1)]
        self.max_power = max(max(e) for row in self.terms for c, e in row)
        self.transform = transform
        self.modulus = None

    def set_modulus(self, modulus):
        self.modulus = modulus
        self.reduced = [[(c % modulus, e) for c, e in row if c % modulus]
                        for row in self.terms]

    def coefficients(self, m, k, ell):
        p = self.modulus
        values = self.transform(m, k, ell, p) if self.transform else (m, k, ell)
        powers = []
        for value in values:
            row = [1]
            for _ in range(self.max_power):
                row.append(row[-1]*value % p)
            powers.append(row)
        x, y, z = powers
        return [sum(c*x[e[0]]*y[e[1]]*z[e[2]] for c, e in row) % p
                for row in self.reduced]

    def possible(self, m, k, ell):
        return integer_has_root(self.coefficients(m, k, ell), self.modulus)


def repeat_residues(residue_bits, modulus, length):
    """Bit i of the result is bit (i mod modulus) of residue_bits.

    The quotient is 1 + 2^p + ... + 2^{p(q-1)}. Multiplication copies the
    p-bit block into disjoint positions; masking truncates the last block.
    """
    blocks = (length+modulus-1)//modulus
    repeats = ((1 << (blocks*modulus))-1)//((1 << modulus)-1)
    return residue_bits*repeats & ((1 << length)-1)


def bit_positions(bits):
    """Enumerate every set bit once, in increasing order."""
    while bits:
        lowest = bits & -bits
        yield lowest.bit_length()-1
        bits ^= lowest


def sieve_three_parameters(test, bounds, primes, label, expected, switch_at=500000):
    """Exact exhaustive sieve, initially using bit sets to avoid huge lists.

    At every stage a record (m,k,bits) represents exactly those triples (m,k,l)
    for which bit l is 1. Once few triples survive, expand them into ordinary
    Python tuples. This changes the representation, not the mathematical set.
    """
    M, K, L = bounds
    rows = [(m, k, (1 << L)-1) for m in range(M) for k in range(K)]
    triples = None
    history = [['initial', M*K*L]]
    print(f'{label}: {M*K*L:,} initial triples.', flush=True)
    for p in primes:
        test.set_modulus(p)
        cache = {}
        if triples is None:
            kept = []
            for m, k, bits in rows:
                key = (m, k % p)
                if key not in cache:
                    residues = sum(1 << ell for ell in range(p)
                                   if test.possible(m, k % p, ell))
                    cache[key] = repeat_residues(residues, p, L)
                bits &= cache[key]
                if bits:
                    kept.append((m, k, bits))
            rows = kept
            count = sum(bits.bit_count() for m, k, bits in rows)
            if count <= switch_at:
                triples = [(m, k, ell) for m, k, bits in rows for ell in bit_positions(bits)]
                assert len(triples) == count
                rows = None
        else:
            kept = []
            for m, k, ell in triples:
                key = (m, k % p, ell % p)
                if key not in cache:
                    cache[key] = test.possible(*key)
                if cache[key]:
                    kept.append((m, k, ell))
            triples = kept
            count = len(triples)
        history.append([int(p), count])
        if p in expected:
            assert count == expected[p], (label, p, count, expected[p])
        print(f'{label}: modulo {p}, {count:,} triples survive.', flush=True)
        if count == 0:
            break
    assert count == 0, (label, count)
    AUDIT[label] = history


def c9_transform(m, k, ell, modulus):
    T = ell-6*k+14*m-14
    Q = ell+6*k-34*m+46
    J = ell+54*k+134*m-434
    return ((10-24*T) % modulus,
            (10+240*Q-1440*T*T) % modulus,
            (10-504*J+15120*T*Q-60480*T*T*T) % modulus)


def small_signature_sieve(poly, variables, n, last):
    m, k, ell = variables
    f = S.Poly(poly.as_expr(), ell)
    terms = [[(int(c), e) for e,c in S.Poly(f.nth(j), m, k).terms()]
             for j in range(f.degree()+1)]
    pairs = [(M,K) for M in range(n+1) for K in range(2**n+n+1)]
    history = [['initial', len(pairs)]]
    for p in S.primerange(2, last+1):
        cache = {}
        kept = []
        for M,K in pairs:
            key = (M,K % p)
            if key not in cache:
                cs = [sum(c*pow(M,e[0],p)*pow(K % p,e[1],p) for c,e in row) % p
                      for row in terms]
                cache[key] = integer_has_root(cs,p)
            if cache[key]:
                kept.append((M,K))
        pairs = kept
        history.append([int(p), len(pairs)])
        if p in (31,61,97,last):
            print(f'n={n}: modulo {p}, {len(pairs)} pairs survive.',flush=True)
        if not pairs:
            break
    assert not pairs
    AUDIT[f'n={n} signature pairs'] = history


def representation_controls():
    # Test every small periodic set directly, including truncated last periods.
    for p in range(1,8):
        for residues in range(1 << p):
            for length in range(1,20):
                actual = repeat_residues(residues,p,length)
                expected = sum(1 << i for i in range(length) if residues >> (i % p) & 1)
                assert actual == expected
                assert list(bit_positions(actual)) == [i for i in range(length) if actual >> i & 1]
    # Compare bit sets to literal triple enumeration using the actual C9 polynomial.
    parameters,root,F = c9_eliminant()
    test = PolynomialRootTest(F,parameters,root,c9_transform)
    literal = {(m,k,l) for m in range(3) for k in range(7) for l in range(19)}
    rows = [(m,k,(1 << 19)-1) for m in range(3) for k in range(7)]
    for p in [3,5,7,11,13,17,19]:
        test.set_modulus(p)
        literal = {q for q in literal if test.possible(*q)}
        rows = [(m,k,bits & repeat_residues(sum(1 << l for l in range(p)
                                                if test.possible(m,k % p,l)),p,19))
                for m,k,bits in rows]
        assert literal == {(m,k,l) for m,k,bits in rows for l in bit_positions(bits)}
    print('Python bit-set representation agrees with literal enumeration.',flush=True)
    return parameters,root,F


def high_index_arithmetic():
    certificates = {
        7: [7,5,11,7,13,7,5,7],
        9: [11,7,7,19,11,17,7,11,7,7],
        11:[13,7,23,13,7,7,11,29,7,17,11,7],
    }
    records = {}
    for n,primes in certificates.items():
        records[n] = []
        for m,p in enumerate(primes):
            variables,equations = derive_equations(n,m)
            G = S.groebner(equations,*reversed(variables),domain=S.QQ)
            last = G.polys[-1].as_expr()
            assert last.free_symbols == {variables[0]}
            F = S.Poly(last,variables[0]).clear_denoms()[1].primitive()[1]
            assert F.degree() == {7:2,9:4,11:9}[n]
            assert G.reduce(F.as_expr())[1] == 0
            residues = [int(F.eval(a) % p) for a in range(p)]
            assert all(residues)
            records[n].append({'m':m,'prime':p,'coefficients':[int(c) for c in F.all_coeffs()],
                               'residue_values':residues})
        print(f'n={n}, r={n-3}: every section count excluded.',flush=True)
    AUDIT['high_index_eliminants'] = records
    high_bernoulli_certificate()
    high_dimension_and_euler_certificates()


def compact_c10_polynomial(raw, variables):
    # A short expression in moments of the symmetric Hilbert numerator.
    m,k,l,j = variables
    T,Q,J,V = S.symbols('T Q J V')
    moments = [j-7*l+20*k-28*m+14,
               j+5*l-40*k+80*m-46,
               j+53*l+80*k-568*m+434,
               j+245*l+3800*k+7280*m-11326]
    solution = S.solve([moments[i]-v for i,v in enumerate([T,Q,J,V])],variables)
    compact = S.Poly(S.expand(raw.as_expr().subs(solution, simultaneous=True)),T,Q,J,V)
    denominator,integer_poly = compact.clear_denoms()
    content,primitive = integer_poly.primitive()
    assert S.expand(compact.as_expr().subs(dict(zip([T,Q,J,V],moments)))-raw.as_expr()) == 0
    AUDIT['c10_compact'] = {'latex':S.latex(primitive.as_expr()),
                           'factor':str(content/denominator),
                           'polynomial':polynomial_record(primitive)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--export-json',type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    projective_space_controls()
    c7_symbolic_certificate()
    c7_finite_certificate()
    c8_symbolic_certificate()
    c8_arithmetic_certificate()
    high_index_arithmetic()
    coindex6_even()
    parameters,root,F = representation_controls()
    AUDIT['c9_quartic'] = polynomial_record(F,root)
    test = PolynomialRootTest(F,parameters,root,c9_transform)
    sieve_three_parameters(test,(10,522,19693),list(S.primerange(13,224)),
                           'C9 r=2',{23:20583262,97:16021,149:136,181:3,223:0})
    variables,h,F,raw = signature_polynomial(10,1)
    assert len(F.terms()) == 126 and raw.content() == 8
    check_signature_reconstruction(10,1,variables,h,raw)
    AUDIT['c10_signature'] = polynomial_record(F,variables[-1])
    compact_c10_polynomial(raw,variables)
    test = PolynomialRootTest(F,variables[:3],variables[3])
    sieve_three_parameters(test,(11,1035,59060),list(S.primerange(3,240)),
                           'C10 r=1',{31:8661101,97:14926,149:148,181:7,239:0})
    variables,h,F,raw = signature_polynomial(14,7)
    assert len(F.terms()) == 120 and F.total_degree() == 7
    check_signature_reconstruction(14,7,variables,h,raw)
    AUDIT['n14_signature'] = polynomial_record(F,variables[-1])
    small_signature_sieve(F,variables,14,113)
    assert [r+7 for r in S.divisors(224) if r % 2] == [8,14]
    AUDIT['elapsed_seconds'] = round(time.monotonic()-started,3)
    if args.export_json:
        args.export_json.write_text(json.dumps(AUDIT,indent=2)+'\n')
    print('ALL PYTHON CHECKS PASSED.',flush=True)
    print('Elapsed seconds:',AUDIT['elapsed_seconds'],flush=True)


if __name__ == '__main__':
    main()
