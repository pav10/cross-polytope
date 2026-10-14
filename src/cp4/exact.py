"""Exact rational linear algebra (fractions.Fraction). No floats in any accept/reject decision."""
from fractions import Fraction as F


def det(M):
    """Exact determinant by Gaussian elimination over Q (int or Fraction entries)."""
    n = len(M)
    if n == 0:
        return F(1)
    A = [[F(x) for x in row] for row in M]
    sign = 1
    for k in range(n):
        piv = next((i for i in range(k, n) if A[i][k] != 0), None)
        if piv is None:
            return F(0)
        if piv != k:
            A[k], A[piv] = A[piv], A[k]
            sign = -sign
        for i in range(k + 1, n):
            f = A[i][k] / A[k][k]
            if f:
                A[i] = [a - f * b for a, b in zip(A[i], A[k])]
    p = F(sign)
    for k in range(n):
        p *= A[k][k]
    return p


def inv(M):
    n = len(M)
    A = [[F(x) for x in row] + [F(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for k in range(n):
        piv = next((i for i in range(k, n) if A[i][k] != 0), None)
        if piv is None:
            raise ZeroDivisionError("singular matrix")
        A[k], A[piv] = A[piv], A[k]
        d = A[k][k]
        A[k] = [x / d for x in A[k]]
        for i in range(n):
            if i != k and A[i][k]:
                f = A[i][k]
                A[i] = [a - f * b for a, b in zip(A[i], A[k])]
    return [row[n:] for row in A]


def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))]
            for i in range(len(A))]


def vecmat(q, A):
    return tuple(sum(q[i] * A[i][j] for i in range(len(A))) for j in range(len(A[0])))


def rref(A):
    R = [[F(x) for x in row] for row in A]
    m = len(R)
    n = len(R[0]) if m else 0
    i = 0
    piv = []
    for j in range(n):
        k = next((k for k in range(i, m) if R[k][j]), None)
        if k is None:
            continue
        R[i], R[k] = R[k], R[i]
        d = R[i][j]
        R[i] = [x / d for x in R[i]]
        for k in range(m):
            if k != i and R[k][j]:
                d = R[k][j]
                R[k] = [x - d * y for x, y in zip(R[k], R[i])]
        piv.append(j)
        i += 1
        if i == m:
            break
    return R[:i], piv


def rank(A):
    return len(rref(A)[1]) if A else 0


def nullspace(A, n=None):
    """Basis of {x : A x = 0} (A = list of rows), exact."""
    if not A:
        return [[F(int(i == j)) for i in range(n)] for j in range(n)]
    n = len(A[0])
    R, piv = rref(A)
    basis = []
    for f in (j for j in range(n) if j not in piv):
        v = [F(0)] * n
        v[f] = F(1)
        for r, p in zip(R, piv):
            v[p] = -r[f]
        basis.append(v)
    return basis


def solve_affine(A, b):
    """(x0, nullspace basis) for A x = b, or None if inconsistent. Exact."""
    n = len(A[0])
    aug = [[F(x) for x in row] + [F(bi)] for row, bi in zip(A, b)]
    R, piv = rref(aug)
    if n in piv:
        return None
    x0 = [F(0)] * n
    for r, p in zip(R, piv):
        x0[p] = r[n]
    return x0, nullspace(A)
