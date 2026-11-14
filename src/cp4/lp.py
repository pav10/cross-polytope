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
    mats, lows, ups = [], [], []                # rows: A_ub (lower -inf) first, then A_eq (lower = upper)
    if A_ub is not None:
        A_ub = np.asarray(A_ub, float)
        mats.append(A_ub); lows.append(np.full(len(A_ub), -inf)); ups.append(np.asarray(b_ub, float))
    if A_eq is not None:
        A_eq = np.asarray(A_eq, float)
        mats.append(A_eq); lows.append(np.asarray(b_eq, float)); ups.append(np.asarray(b_eq, float))
    if mats:
        A = np.vstack(mats); mask = A != 0
        start = np.concatenate([[0], np.cumsum(mask.sum(axis=1))]).astype(np.int32)
        rr, cc = np.nonzero(mask)
        index, value = cc.astype(np.int32), A[rr, cc]
        lo, up = np.concatenate(lows), np.concatenate(ups)
    else:
        start, index, value = np.zeros(1, np.int32), np.zeros(0, np.int32), np.zeros(0)
        lo = up = np.zeros(0)
    lp = highspy.HighsLp()
    lp.num_col_, lp.num_row_ = len(c), len(lo)
    lp.col_cost_ = np.asarray(c, float)
    lp.col_lower_ = np.array([-inf if b[0] is None else b[0] for b in bounds], float)
    lp.col_upper_ = np.array([inf if b[1] is None else b[1] for b in bounds], float)
    lp.row_lower_ = lo
    lp.row_upper_ = up
    lp.a_matrix_.format_ = highspy.MatrixFormat.kRowwise
    lp.a_matrix_.start_ = start
    lp.a_matrix_.index_ = index
    lp.a_matrix_.value_ = value
    _h.clearModel()
    _h.passModel(lp)
    _h.run()
    if _h.getModelStatus() != highspy.HighsModelStatus.kOptimal:
        return False, None, None
    return True, np.array(_h.getSolution().col_value), _h.getInfo().objective_function_value


solve = solve_highs if BACKEND == 'highspy' else solve_scipy
