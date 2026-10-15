"""LP backend for the float prunes of Stage B1 (faces.py) and the merged search (merged.py).

solve(c, A_ub, b_ub, A_eq, b_eq, bounds) -> (ok, x, fun): minimise c.x subject to A_ub x <= b_ub,
A_eq x = b_eq and bounds, with scipy.optimize.linprog semantics (ok iff an optimum was found).
Backend: highspy (direct HiGHS API, much less per-call overhead) when importable, else scipy's HiGHS.
Set CP4_LP=scipy to force scipy.  Both backends run HiGHS; scripts/check_lp_backend.py checks that they
give the same verdicts on LPs recorded from B1 probes.  Floats prune only (CLAUDE.md hard rules).
"""
import os
import numpy as np
from scipy.optimize import linprog

try:
    import highspy
    _h = highspy.Highs()
    _h.setOptionValue('output_flag', False)
except ImportError:
    highspy = None
BACKEND = 'highspy' if highspy is not None and os.environ.get('CP4_LP') != 'scipy' else 'scipy'


def solve_scipy(c, A_ub=None, b_ub=None, A_eq=None, b_eq=None, bounds=None):
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method='highs')
    return (True, res.x, res.fun) if res.status == 0 else (False, None, None)


def solve_highs(c, A_ub=None, b_ub=None, A_eq=None, b_eq=None, bounds=None):
    inf = highspy.kHighsInf
    rows = []                                   # (indices, values, lower, upper)
    for A, lo, hi in ((A_ub, None, b_ub), (A_eq, b_eq, b_eq)):
        if A is None:
            continue
        for k, r in enumerate(np.asarray(A, float)):
            nz = np.nonzero(r)[0]
            rows.append((nz, r[nz], -inf if lo is None else lo[k], hi[k]))
    lp = highspy.HighsLp()
    lp.num_col_, lp.num_row_ = len(c), len(rows)
    lp.col_cost_ = np.asarray(c, float)
    lp.col_lower_ = np.array([-inf if b[0] is None else b[0] for b in bounds], float)
    lp.col_upper_ = np.array([inf if b[1] is None else b[1] for b in bounds], float)
    lp.row_lower_ = np.array([r[2] for r in rows], float)
    lp.row_upper_ = np.array([r[3] for r in rows], float)
    lp.a_matrix_.format_ = highspy.MatrixFormat.kRowwise
    lp.a_matrix_.start_ = np.cumsum([0] + [len(r[0]) for r in rows]).astype(np.int32)
    lp.a_matrix_.index_ = np.concatenate([r[0] for r in rows] or [[]]).astype(np.int32)
    lp.a_matrix_.value_ = np.concatenate([r[1] for r in rows] or [[]]).astype(float)
    _h.clearModel()
    _h.passModel(lp)
    _h.run()
    if _h.getModelStatus() != highspy.HighsModelStatus.kOptimal:
        return False, None, None
    return True, np.array(_h.getSolution().col_value), _h.getInfo().objective_function_value


solve = solve_highs if BACKEND == 'highspy' else solve_scipy
