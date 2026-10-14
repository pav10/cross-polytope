#!/usr/bin/env python3
"""Exact checks for cross_polytope_audit.tex. Python 3, standard library only.

This verifies finite arithmetic, not the global packing conjecture and not Lean.
Run: python3 verify.py [--output certificate.json]
All mathematical acceptance checks use integers or fractions.Fraction.
"""
from fractions import Fraction as F
from itertools import product, combinations, combinations_with_replacement, permutations
from collections import Counter, defaultdict
from math import floor
from decimal import Decimal, localcontext
import argparse
import json

A = [[20,53,-53,-42],[-62,-22,42,-42],[-22,-22,-62,62],[20,53,42,53]]
COUNTEREXAMPLE = [[238,631,-631,-500],[-738,-262,500,-500],
                 [-262,-262,-738,738],[238,631,500,631]]
M0 = [[2,2,-2,-2],[-4,0,2,-2],[0,0,-4,4],[2,2,2,2]]
S = [[-2,1,-1,0],[2,-2,0,0],[-2,-2,2,-2],[-2,1,0,1]]
U = [[0,-1,-1,0],[0,0,0,0],[0,0,0,0],[0,1,0,-1]]
V = [[1,0,0,1],[0,0,-1,-1],[0,0,0,0],[-1,0,1,0]]
TYPES = [(0,0,0,1),(0,1,1,2),(1,1,1,2),(1,1,1,3),(1,1,2,4),(1,2,2,5)]
PATTERNS = [(0,0,1,1),(0,1,1,1),(1,1,1,1),(0,1,1,2),
            (1,1,1,2),(1,1,2,2),(1,1,1,3),(1,1,2,3),
            (1,2,2,3),(1,1,2,4),(1,2,3,4),(1,2,2,5),
            (1,2,3,5),(2,3,4,5)]
CAPS = [F(19,10),F(11,4),F(17,10),F(19,10)]

def det(A):
    if not A:
        return 1
    return sum((-1)**j*A[0][j]*det([r[:j]+r[j+1:] for r in A[1:]])
               for j in range(len(A)))

def inv(A):
    n=len(A); d=det(A)
    assert d
    return [[F((-1)**(i+j)*det([r[:i]+r[i+1:] for k,r in enumerate(A) if k!=j]),d)
             for j in range(n)] for i in range(n)]

def mul(q,A):
    return tuple(sum(q[i]*A[i][j] for i in range(len(A))) for j in range(len(A)))

def l1(v): return sum(map(abs,v))
def positive(v): return next((x for x in v if x),0)>0
def flatten(A): return [x for r in A for x in r]

def short_pairs(A,radius):
    Ai=inv(A)
    bounds=[floor(radius*max(abs(Ai[i][j]) for i in range(4))) for j in range(4)]
    pairs=[]
    for q in product(*[range(-b,b+1) for b in bounds]):
        if positive(q) and l1(mul(q,A))<=radius:
            pairs.append(q)
            assert l1(mul(q,A))==radius
    return bounds,pairs

def rref(A):
    R=[[F(x) for x in row] for row in A]; m=len(R); n=len(R[0]); i=0; piv=[]
    for j in range(n):
        k=next((k for k in range(i,m) if R[k][j]),None)
        if k is None: continue
        R[i],R[k]=R[k],R[i]; d=R[i][j]; R[i]=[x/d for x in R[i]]
        for k in range(m):
            if k!=i and R[k][j]:
                d=R[k][j]; R[k]=[x-d*y for x,y in zip(R[k],R[i])]
        piv.append(j); i+=1
        if i==m: break
    return R,piv

def diagonals(p,n):
    if n==1:
        yield (p,); return
    for d in range(1,p+1):
        if p%d==0:
            for ds in diagonals(p//d,n-1): yield (d,)+ds

def hnfs(p,n=4):
    ij=[(i,j) for i in range(n) for j in range(i+1,n)]
    for ds in diagonals(p,n):
        for off in product(*[range(ds[j]) for i,j in ij]):
            H=[[0]*n for _ in range(n)]
            for i in range(n): H[i][i]=ds[i]
            for (i,j),a in zip(ij,off): H[i][j]=a
            yield H

def residue_group(gens,p,require_clean=False):
    """Generate a subgroup modulo p; cleanliness means Lee weight > p."""
    n=len(gens[0]); zero=(0,)*n; seen={zero}; pending=[zero]
    while pending:
        v=pending.pop()
        for g in gens:
            w=tuple((x+y)%p for x,y in zip(v,g))
            if w in seen: continue
            if require_clean and sum(min(x,p-x) for x in w)<=p: return None
            seen.add(w); pending.append(w)
    return tuple(sorted(seen))

def dual_residues(H,p):
    n=len(H); J=[[0]*n for _ in range(n)]
    for i in reversed(range(n)):
        J[i][i]=p//H[i][i]
        for j in range(i+1,n):
            a=-sum(H[i][k]*J[k][j] for k in range(i+1,j+1))
            assert a%H[i][i]==0
            J[i][j]=a//H[i][i]
    gens=[tuple(J[j][i]%p for j in range(n)) for i in range(n)]
    if any(any(g) and sum(min(x,p-x) for x in g)<=p for g in gens): return None
    G=residue_group(gens,p,True)
    if G is not None: assert len(G)==p
    return G

def classification():
    counts={}; found=set()
    for p in range(1,25):
        total=good=0
        for H in hnfs(p):
            total+=1; G=dual_residues(H,p)
            if G is not None: found.add((p,G)); good+=1
        counts[p]=(total,good)
    expected=set()
    for a in TYPES:
        p=a[3]; g=(-a[0],-a[1],-a[2],1)
        for perm in permutations(range(4)):
            for signs in product((-1,1),repeat=4):
                h=tuple(signs[i]*g[perm[i]]%p for i in range(4))
                G=residue_group([h],p)
                assert len(G)==p
                expected.add((p,G))
    assert found==expected
    assert sum(x[0] for x in counts.values())==201677
    assert len(found)==54
    return counts

def clean(V):
    p=abs(det(V))
    assert p
    if p==1: return True
    # This cutoff is justified by classification(), which is run first.
    if p>5: return False
    Vi=inv(V)
    gens=[tuple(int(p*x)%p for x in row) for row in Vi]
    G=residue_group(gens,p,True)
    if G is not None: assert len(G)==p
    return G is not None

def finite_tests():
    out=[]
    for a in TYPES:
        p=a[3]; P=[[1,0,0,0],[0,1,0,0],[0,0,1,0],list(a)]
        vectors=[]
        for r in product(range(-5,6),repeat=4):
            if not positive(r): continue
            if any((r[i]+a[i]*r[3])%p for i in range(3)): continue
            q=tuple((r[i]+a[i]*r[3])//p for i in range(3))+(r[3],)
            if list(q) in P or [-x for x in q] in P: continue
            good=True
            for i in range(4):
                if not r[i]: continue
                W=[row[:] for row in P]; W[i]=list(q)
                if not clean(W): good=False; break
            if good: vectors.append(q)
        out.append(vectors)
    assert [len(T) for T in out]==[1156,182,196,65,29,16]
    # Independent type-1 generation, using only the universal index bound 24.
    patterns=[]
    for c in combinations_with_replacement(range(25),4):
        if not any(c): continue
        good=True
        for i,k in enumerate(c):
            for t in range(1,k):
                mass=t+sum(min((t*c[j])%k,(-t*c[j])%k) for j in range(4) if j!=i)
                if mass<=k: good=False; break
            if not good: break
        if good and sum(c)>1: patterns.append(c)
    assert set(patterns)==set(PATTERNS)
    assert {tuple(sorted(map(abs,q))) for q in out[0]}==set(PATTERNS)
    return out

def polynomial_det(matrices):
    n=len(matrices[0]); m=len(matrices); out=defaultdict(int)
    for choices in product(range(m),repeat=n):
        ex=tuple(choices.count(j) for j in range(1,m))
        out[ex]-=det([matrices[choices[i]][i] for i in range(n)])
    return {ex:c for ex,c in out.items() if c}

def interval(poly,lo,hi):
    assert 0<=lo<=hi
    return (sum(c*(lo**k if c>=0 else hi**k) for k,c in enumerate(poly)),
            sum(c*(hi**k if c>=0 else lo**k) for k,c in enumerate(poly)))

def quadratic_min(p,lo,hi):
    a,b,c=(list(p)+[0,0,0])[:3]
    pts=[lo,hi]
    if c>0 and lo<=-F(b,2*c)<=hi: pts.append(-F(b,2*c))
    return min(a+b*t+c*t*t for t in pts)

def candidate_checks(T1):
    bounds,P=short_pairs(A,168)
    b2,P2=short_pairs(COUNTEREXAMPLE,2000)
    assert bounds==b2==[1,2,1,1] and P==P2 and len(P)==15
    old=F(abs(det(A)),84**4); new=F(abs(det(COUNTEREXAMPLE)),1000**4)
    assert old-new==F(214636001,12155062500000000)>0
    assert new==F(50513750079,62500000000)
    sig=[tuple(1 if x>0 else -1 for x in mul(q,A)) for q in P]
    R=[[q[i]*s[j] for i in range(4) for j in range(4)] for q,s in zip(P,sig)]
    assert len(rref(R)[1])==13
    for row in R:
        assert sum(x*y for x,y in zip(row,flatten(M0)))==8
        for H in [S,U,V]: assert sum(x*y for x,y in zip(row,flatten(H)))==0
    assert len(rref([flatten(H) for H in [S,U,V]])[1])==3
    # All 1156 inequalities hold throughout the full parameter interval.
    for q in T1:
        s=tuple(1 if x>=0 else -1 for x in mul(q,A))
        a,b=mul(q,M0),mul(q,S)
        for z in [F(1,2),F(2,3)]:
            assert sum(sj*(aj+z*bj) for sj,aj,bj in zip(s,a,b))>=8
    for q,s in zip(P,sig):
        a,b=mul(q,M0),mul(q,S)
        for z in [F(1,2),F(2,3)]:
            assert all(sj*(aj+z*bj)>=0 for sj,aj,bj in zip(s,a,b))
            assert sum(sj*(aj+z*bj) for sj,aj,bj in zip(s,a,b))==8
    expected={(0,0,0):256,(1,0,0):-192,(2,0,0):192,(4,0,0):-16,
              (0,2,0):32,(1,2,0):-32,(2,2,0):16,(2,1,1):16,
              (0,0,2):32,(1,0,2):-16,(2,0,2):8}
    assert polynomial_det([M0,S,U,V])==expected
    # Multipliers for the objective -det M, modulo z^3-6z+3.
    L=[[29,-56,26],[-8,12,22],[10,-8,2],[4,4,10],[24,-28,10],
       [14,8,-22],[3,4,-4],[29,-56,26],[-9,40,-22],[10,-8,2],
       [24,-28,10],[-11,40,-8],[3,4,-4],[3,0,0],[3,0,0]]
    for i in range(4):
        for j in range(4):
            C0=[row[:j]+row[j+1:] for k,row in enumerate(M0) if k!=i]
            C1=[row[:j]+row[j+1:] for k,row in enumerate(S) if k!=i]
            g=[-(-1)**(i+j)*sum(det([C1[k] if k in I else C0[k] for k in range(3)])
                for I in combinations(range(3),d)) for d in range(4)]
            g[0]-=3*g[3]; g[1]+=6*g[3]
            assert g[:3]==[sum(L[k][d]*R[k][i*4+j] for k in range(15)) for d in range(3)]
    assert all(quadratic_min(l,F(1,2),F(2,3))>0 for l in L)
    # Uniform inverse bounds near the algebraic root, proving exactly 15 pairs.
    lo,hi=F(5239,10000),F(131,250)
    assert lo**3-6*lo+3>0 and hi**3-6*hi+3<0
    dlo=interval([256,-192,192,0,-16],lo,hi)[0]
    inverse_bounds=[]
    for j in range(4):
        maxco=F(0)
        for i in range(4):
            C0=[row[:i]+row[i+1:] for k,row in enumerate(M0) if k!=j]
            C1=[row[:i]+row[i+1:] for k,row in enumerate(S) if k!=j]
            p=[sum(det([C1[k] if k in I else C0[k] for k in range(3)])
                 for I in combinations(range(3),d)) for d in range(4)]
            a,b=interval(p,lo,hi); maxco=max(maxco,abs(a),abs(b))
        inverse_bounds.append(8*maxco/dlo)
    assert all(x<cap for x,cap in zip(inverse_bounds,CAPS))
    exact=[]; gap=F(999)
    for q in product(range(-1,2),range(-2,3),range(-1,2),range(-1,2)):
        if not positive(q): continue
        s=tuple(1 if x>=0 else -1 for x in mul(q,A)); a,b=mul(q,M0),mul(q,S)
        endpoints=[sum(sj*(aj+z*bj) for sj,aj,bj in zip(s,a,b)) for z in [lo,hi]]
        assert min(endpoints)>=8
        if endpoints==[8,8]: exact.append(q)
        else: gap=min(gap,min(endpoints)-8)
    assert exact==P and gap==F(5239,5000)
    # The six types are not mutually exclusive lattice classes.
    indices=Counter(); clean_indices=Counter()
    for t in combinations(P,4):
        W=[list(q) for q in t]; p=abs(det(W)); indices[p]+=1
        if p and clean(W): clean_indices[p]+=1
    assert dict(indices)=={0:420,1:831,2:104,3:9,4:1}
    assert dict(clean_indices)=={1:831,2:104,3:9,4:1}
    with localcontext() as ctx:
        ctx.prec=60; a,b=Decimal('0.5239'),Decimal('0.5240')
        for _ in range(210):
            z=(a+b)/2
            if z**3-6*z+3>0: a=z
            else: b=z
        z=(a+b)/2; d=(16-9*z+6*z*z)/16; density=Decimal(2)/3/d
        decimals={'z':str(z),'determinant':str(d),'density':str(density)}
    return {'old_determinant':str(old),'rational_determinant':str(new),
            'rational_density':str(F(2,3)/new),'determinant_improvement':str(old-new),
            'short_coefficient_bounds':bounds,'minimal_pairs':P,'facet_signs':sig,
            'contact_rank':13,'local_multipliers':L,'uniform_inverse_bounds':list(map(str,inverse_bounds)),
            'minimal_quadruple_indices':dict(indices),'clean_quadruple_indices':dict(clean_indices),
            'algebraic_candidate_decimal_approximations':decimals}

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output',default='certificate.json')
    args=parser.parse_args()
    counts=classification()
    tests=finite_tests()
    candidate=candidate_checks(tests[0])
    certificate={'status':'Exact finite checks passed; global optimality remains unproved.',
                 'classification_counts':counts,'test_counts':[len(T) for T in tests],
                 'type1_patterns':PATTERNS,'test_sets':tests,'candidate':candidate}
    with open(args.output,'w') as f: json.dump(certificate,f,indent=2)
    print('PASS: 201677 HNFs, 54 clean superlattices, six symmetry orbits.')
    print('PASS: six finite test sets with pair counts',certificate['test_counts'])
    print('PASS: rational counterexample, interval admissibility, exact contacts, rank, polynomial, KKT identity.')
    print('Algebraic candidate density:',candidate['algebraic_candidate_decimal_approximations']['density'])
    print('Global optimality and a Lean kernel check are not claimed.')

if __name__=='__main__': main()
