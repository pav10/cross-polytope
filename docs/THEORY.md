# The densest lattice packing of the four-dimensional cross-polytope
## Target theorem, complete proofs of everything provable now, and the exact finite remainder

*Status date: 23 September 2026. Every statement is tagged:*
**[P]** *proved here in full;* **[P+C]** *proved, with a finite exact computation as one step (script named);*
**[Cite]** *classical or published result used as a black box;* **[Conj]** *conjectural;*
**[Num]** *numerical evidence only, not a proof.*

---

# 0. Normalization and the target statement

Let $\|x\| = \|x\|_1 = \sum_k |x_k|$ on $\mathbb R^4$, $C_4=\{\|x\|\le 1\}$, $K=2C_4$. A lattice $\Lambda\subset\mathbb R^4$ is **admissible** if $\|v\|\ge 2$ for all $v\in\Lambda\setminus\{0\}$. Lattices are written $\Lambda=\mathbb Z^4B$ (rows of $B$ = a basis). A **minimal vector** is $v\in\Lambda$ with $\|v\|=2$; $M(\Lambda)$ is the set of minimal vectors.

$$D_{\mathrm{pack}}=\inf\{\det\Lambda:\Lambda\text{ admissible}\},\qquad \delta_L(C_4)=\frac{\operatorname{vol}C_4}{D_{\mathrm{pack}}}=\frac{2/3}{D_{\mathrm{pack}}}.$$

In classical terminology $D_{\mathrm{pack}}=\Delta(2C_4)=16\,\Delta(C_4)$. *(Our old drafts called $D_{\mathrm{pack}}$ "the critical determinant of $C_4$"; that is off by a factor $16$.)*

Let $\zeta\approx 0.5239763970818660$ be the unique root of $z^3-6z+3$ in $(\tfrac12,\tfrac23)$, and put
$$M(z)=\begin{pmatrix}2-2z&2+z&-2-z&-2\\2z-4&-2z&2&-2\\-2z&-2z&2z-4&4-2z\\2-2z&2+z&2&2+z\end{pmatrix},\quad B(z)=\tfrac14M(z),\quad \Lambda_\zeta=\mathbb Z^4B(\zeta),$$
$$D_*=\det\Lambda_\zeta=\frac{16-9\zeta+6\zeta^2}{16}=0.80822000090353537\ldots,\qquad \frac{2/3}{D_*}=0.82485791730145054\ldots$$

> **Target Theorem Z.** $D_{\mathrm{pack}}=D_*$. A lattice attains it iff it is the image of $\Lambda_\zeta$ under a signed coordinate permutation.

**What this document establishes.**
- $D_{\mathrm{pack}}\le D_*$ **[P+C]** (§2).
- $\Lambda_\zeta$ is a strict local minimum **[P+C]** (§3).
- Critical lattices satisfy strong necessary conditions **[P]** (§4–§6).
- Critical lattices fall into exactly **39 strata** **[P+C]** (§6).
- The smallest stratum is completely classified **[P+C]** (§7).
- Theorem Z is **equivalent to a finite, exactly specified computation** **[P]** (§9).

The global inequality itself remains **[Conj]**.

---

# 1. Basic tools

**Lemma 1.1 (packing ⟺ admissibility) [P].** The translates $x+C_4$ ($x\in\Lambda$) have disjoint interiors iff $\Lambda$ is admissible.

*Proof.* Suppose interiors of $x+C_4$ and $y+C_4$ meet at $z$. Then $\|x-y\|\le\|x-z\|+\|z-y\|<2$.

Conversely, if $\|x-y\|<2$, put $z=(x+y)/2$. Then $\|z-x\|=\|z-y\|<1$, so $z$ is interior to both. $\square$

**Theorem 1.2 (Minkowski, strict form) [Cite].** If $S\subset\mathbb R^n$ is convex and symmetric and $\operatorname{vol}S> 2^n\det\Gamma$, then $\operatorname{int}S$ contains a nonzero point of the lattice $\Gamma$. *(All uses below are in the contrapositive: if $\operatorname{int}S$ contains no nonzero point of $\Gamma$, then $\operatorname{vol}S\le 2^n\det\Gamma$.)*

**Theorem 1.3 (van der Corput) [Cite].** If moreover $\operatorname{vol}S>k\,2^n\det\Gamma$, then $S$ contains at least $k$ pairs $\pm w$ of nonzero points of $\Gamma$.

**Corollary 1.4 [P].** Every admissible $\Lambda$ has $\det\Lambda\ge 2/3$.

*Proof.* $\operatorname{int}K$ contains no nonzero point of $\Lambda$, and $\operatorname{vol}K=16\cdot\tfrac23$. By Theorem 1.2, $16\cdot\tfrac23\le 16\det\Lambda$. $\square$

**Lemma 1.5 (primitivity) [P].** If $m\in M(\Lambda)$ and $m=kz$ with $z\in\Lambda$, $k\in\mathbb Z$, then $|k|=1$.

*Proof.* Otherwise $0<\|z\|=2/|k|<2$. $\square$

**Lemma 1.6 (parity–support) [P].** Let $m,n\in M(\Lambda)$ with $m\ne\pm n$ and $m\equiv n\pmod{2\Lambda}$. Then their supports are disjoint: $m_kn_k=0$ for every $k$.

*Proof.* Both $(m+n)/2$ and $(m-n)/2$ lie in $\Lambda\setminus\{0\}$, so each has norm $\ge 2$. The triangle inequality gives $\|m\pm n\|\le 4$, so both norms equal $4$.
- Equality in $\|m+n\|=\|m\|+\|n\|$ forces $m_kn_k\ge0$ for all $k$.
- Equality in $\|m-n\|=\|m\|+\|n\|$ forces $m_kn_k\le 0$ for all $k$.

Hence $m_kn_k=0$. $\square$

### 1.7 The finite admissibility test for an arbitrary frame

This generalizes the auditor's descent theorem to frames that are not clean, which is what the max-det frame of §6 requires.

**Setup.** Let $\Gamma\supseteq\mathbb Z^n$ be a lattice with $[\Gamma:\mathbb Z^n]=p$, so $p\Gamma\subseteq\mathbb Z^n$. Let $F$ be a seminorm on $\mathbb R^n$ with $F(e_i)=2$ for all $i$. Put $\mu(c)=\|c\|_1$.

For $c\in\mathbb R^n$ with $c_i\ne 0$, the vectors $\{e_j\}_{j\ne i}\cup\{c\}$ form a basis. So every $w$ has a unique expansion $w=\sum_{j\ne i}\alpha_je_j+\beta c$, with $\beta=w_i/c_i$ the **$c$-coefficient**.

Let $V_i(c)=\operatorname{conv}(\pm e_j\,(j\ne i),\pm c)$ and define
$$T(\Gamma)=\{c\in\Gamma\setminus0:\mu(c)\le1\}\ \cup\ \Bigl\{c\in\Gamma:\mu(c)>1,\ \forall i\text{ with }c_i\neq0:\ \text{the only points of }\Gamma\cap V_i(c)\text{ with nonzero }c\text{-coefficient are }\pm c\Bigr\}.$$

**Theorem 1.7 [P].**
- (a) If $F\ge 2$ on $T(\Gamma)$, then $F\ge2$ on $\Gamma\setminus\{0\}$.
- (b) $T(\Gamma)$ is finite. Explicitly, if $n_i$ is the number of pairs $\pm w$ of nonzero points of $\Gamma$ in $O_i=\operatorname{conv}(\pm e_j:j\ne i)$, then every $c\in T(\Gamma)$ with $\mu(c)>1$ satisfies $|c_i|\le (n_i+2)\,n!/p$ whenever $c_i\ne0$.

*Proof of (a).* Suppose some nonzero $c\in\Gamma$ has $F(c)<2$; call such $c$ *violating*. The values $\mu(\Gamma)\subseteq\frac1p\mathbb Z_{\ge0}$ form a discrete set, so we may choose a violating $c$ with $\mu(c)$ minimal.

If $\mu(c)\le1$, then $c\in T(\Gamma)$ and $F(c)\ge2$, a contradiction. So $\mu(c)>1$.

Fix $i$ with $c_i\ne0$. Let $w\in\Gamma\cap V_i(c)$ with $w\ne\pm c$ and nonzero $c$-coefficient $\beta$. Write $w=\sum_{j\ne i}\alpha_je_j+\beta c$ with $\sum|\alpha_j|+|\beta|\le 1$.

- **$|\beta|<1$:** if $|\beta|=1$ then all $\alpha_j=0$ and $w=\pm c$.
- **$w\ne0$:** its $i$-th coordinate is $\beta c_i\ne0$.
- **$w$ is violating:** since $|\beta|>0$ and $F(c)<2$, $F(w)\le 2\sum|\alpha_j|+|\beta|F(c)<2(\sum|\alpha_j|+|\beta|)\le 2$.
- **$w$ has smaller $\mu$:** $\mu(w)\le\sum|\alpha_j|+|\beta|\mu(c)\le(1-|\beta|)+|\beta|\mu(c)$. This is $<\mu(c)$, because the inequality is equivalent to $(1-|\beta|)<(1-|\beta|)\mu(c)$, which holds as $\mu(c)>1$ and $|\beta|<1$.

This contradicts the minimality of $\mu(c)$. Hence no such $w$ exists for any $i$ with $c_i\ne0$. That places $c$ in $T(\Gamma)$, while $F(c)<2$, a contradiction. $\square$

*Proof of (b).* Fix $i$ and $c\in T(\Gamma)$ with $\mu(c)>1$, $c_i\ne0$. The body $V_i(c)$ is symmetric and convex, with
$$\operatorname{vol}V_i(c)=\frac{2^n}{n!}\,|\det(e_1,\dots,c,\dots,e_n)|=\frac{2^n}{n!}|c_i|.$$

Count the pairs of $\Gamma$-points in $V_i(c)$:
- Points with $c$-coefficient $0$ lie in $O_i$: at most $n_i$ pairs.
- Points with nonzero $c$-coefficient: only $\pm c$ (one pair), by $c\in T(\Gamma)$.

So $V_i(c)$ contains at most $n_i+1$ pairs. If $\operatorname{vol}V_i(c)>(n_i+2)\,2^n\det\Gamma$, Theorem 1.3 would give at least $n_i+2$ pairs. Hence, with $\det\Gamma=1/p$, $\operatorname{vol}V_i(c)\le (n_i+2)2^n/p$. That is $|c_i|\le(n_i+2)n!/p$. So $c$ lies in a bounded region, and $T(\Gamma)$ is finite. $\square$

**Remark.** For $\Gamma=\mathbb Z^n$ (a basis frame), $\{c:\mu(c)\le1\}\cap\Gamma=\{\pm e_i\}$. The condition on points with zero $c$-coefficient is then automatic, and $T(\mathbb Z^4)$ is exactly the classical Wolff test set $T_1$ ($1156$ pairs, 14 absolute patterns). This was regenerated independently three times (our `wolff_lists.py`, the auditor's `verify.py`, and the mod-$k$ criterion).

**Corollary 1.8 [P].** If $F(x)=\|xX\|$ for a real $4\times4$ matrix $X$, and $F$ satisfies the hypotheses of Theorem 1.7(a), then $X$ is nonsingular and $|\det X|\ge\frac23\cdot p$. In particular, no "singular feasible point" can pass a correctly encoded test.

*Proof.* Theorem 1.7 gives $F\ge2$ on $\Gamma\setminus0$.

If $X$ were singular, the set $\{F<2\}$ would contain an infinite symmetric cylinder around $\ker X^{\!\top}$. That cylinder contains convex symmetric bodies of arbitrarily large volume, which by Theorem 1.2 must contain nonzero $\Gamma$-points. This is impossible, so $X$ is nonsingular.

Then $\Gamma X$ is admissible, and Corollary 1.4 gives $\det(\Gamma X)=|\det X|/p\ge 2/3$. $\square$

---

# 2. The upper bound $D_{\mathrm{pack}}\le D_*$

**Theorem 2.1 [P+C].** For every $z\in[\tfrac12,\tfrac23]$ the lattice $\mathbb Z^4B(z)$ is admissible, and
$$\det B(z)=-\frac{16-12z+12z^2-z^4}{16}.$$
Hence $D_{\mathrm{pack}}\le\min_{[1/2,2/3]}|\det B|=D_*$, the minimum being attained exactly at $z=\zeta$.

*Proof.*

*(i) The rows are minimal on the whole interval.* For $z\in[\frac12,\frac23]$ the signs of all row entries are constant, and the four row norms of $M(z)$ are:
- row 1: $(2-2z)+(2+z)+(2+z)+2=8$;
- row 2: $(4-2z)+2z+2+2=8$;
- row 3: $2z+2z+(4-2z)+(4-2z)=8$;
- row 4: $(2-2z)+(2+z)+2+(2+z)=8$.

So $F_z(x)=\|xM(z)\|$ satisfies $F_z(e_i)=8$.

*(ii) The finite test holds on the whole interval.* For each of the 1156 pairs $q\in T(\mathbb Z^4)$, let $\tau_q\in\{\pm1\}^4$ be the coordinate signs of $qM(11/21)$, with $+1$ for zeros. The scalar $\tau_q\cdot qM(z)$ is affine in $z$. It is verified in exact rational arithmetic to be $\ge8$ at $z=\frac12$ and at $z=\frac23$, hence on the whole interval. Then $F_z(q)\ge\tau_q\cdot qM(z)\ge 8$. *(Computation: `verify.py`, audited; re-run and passing.)*

*(iii) Admissibility and the determinant.* Theorem 1.7(a), with $\Gamma=\mathbb Z^4$ and the scaling $2\leftrightarrow 8$, gives admissibility of $\mathbb Z^4M(z)$ at minimum $8$, hence of $\mathbb Z^4B(z)$ at minimum $2$.

The determinant formula is a polynomial identity, checked by expansion.

*(iv) Minimizing along the family.* Put $d(z)=(16-12z+12z^2-z^4)/16$. Then $d'(z)=-\frac14(z^3-6z+3)$ and $d''(z)=\frac34(2-z^2)>0$ on the interval. The cubic takes the value $\frac18>0$ at $\frac12$ and $-\frac{19}{27}<0$ at $\frac23$, and it is strictly decreasing between them. So it has the unique root $\zeta$, which is the unique minimum of $d$.

Using $\zeta^4=6\zeta^2-3\zeta$ gives $d(\zeta)=(16-9\zeta+6\zeta^2)/16$. $\square$

**Remark 2.2 (why $L_4$ failed) [P].** $L_4=\mathbb Z^4B(11/21)$ and $d'(11/21)=-2/9261<0$, so $L_4$ is not even stationary. The rational point $z=131/250$ gives an explicit admissible lattice with $\det=50513750079/62500000000<\det L_4$. We verified this independently in exact arithmetic: $\lambda_1$ was checked by exhaustive integer search, and it has 30 minimal vectors.

---

# 3. $\Lambda_\zeta$ is a strict local minimum

We work at minimum $8$ with $Y_*=M(\zeta)$.

**Contact data** (auditor's table; every entry re-verified). For each of the 15 pairs $q$ below, $\|qM(z)\|=8$ on the whole interval, with supporting sign vector $\sigma_q$, and $\lambda_q(z)>0$ on $[\frac12,\frac23]$.

| $q$ | $\sigma_q$ | $\lambda_q(z)$ | $q$ | $\sigma_q$ | $\lambda_q(z)$ |
|---|---|---|---|---|---|
| 0001 | ++++ | $29-56z+26z^2$ | 1001 | ++−+ | $-9+40z-22z^2$ |
| 0010 | −−−+ | $-8+12z+22z^2$ | 1010 | −+−+ | $10-8z+2z^2$ |
| 0011 | −+−+ | $10-8z+2z^2$ | 1100 | −+−− | $24-28z+10z^2$ |
| 0100 | −−+− | $4+4z+10z^2$ | 1101 | −++− | $-11+40z-8z^2$ |
| 0101 | −+++ | $24-28z+10z^2$ | 1110 | −+−− | $3+4z-4z^2$ |
| 0110 | −−−+ | $14+8z-22z^2$ | 1111 | −+−+ | $3$ |
| 0111 | −+++ | $3+4z-4z^2$ | 1211 | −++− | $3$ |
| 1000 | ++−− | $29-56z+26z^2$ | | | |

**Theorem 3.1 [P+C].** $\Lambda_\zeta$ is a strict local minimum of $\det$ among admissible lattices. That is: there is $\varepsilon>0$ such that every admissible $\mathbb Z^4Y$ with $0<\|Y-Y_*\|<\varepsilon$ has $|\det Y|>|\det Y_*|$.

*Proof.*

*(1) Exactly these 15 pairs are contacts near $\zeta$, and all their coordinates are nonzero.* Interval bounds on $M(z)^{-1}$ over $[5239/10000,\,131/250]\ni\zeta$ confine every vector of norm $\le 8$ to the coefficient box $(1,2,1,1)$. In that box the 15 listed pairs have norm exactly $8$, and all other pairs have norm $\ge 8+5239/5000$ (exact, `verify.py`).

So the non-contacts have a uniform gap. Moreover every coordinate of every $qY_*$ is nonzero, so each contact lies in the relative interior of a single facet of $8C_4$.

*(2) Linearization of admissibility.* For $Y=Y_*+H$ with $H$ small, the gap in (1) means only the 15 contacts matter. For each contact, the facet in (1) is unique, so $\|q(Y_*+H)\|=\sigma_q\cdot q(Y_*+H)=8+\ell_q(H)$, where $\ell_q(H)=\sigma_q\cdot qH$.

Hence, near $Y_*$: $\mathbb Z^4Y$ is admissible $\iff$ $\ell_q(H)\ge0$ for all 15 $q$.

*(3) First-order condition (KKT) with strictly positive multipliers.* Let $f=-\det$, so $f>0$ near $Y_*$. Using cofactor expansion and reducing modulo $z^3-6z+3$ gives the exact identity
$$Df(Y_*)[H]=\sum_q\lambda_q(\zeta)\,\ell_q(H),\qquad \lambda_q(\zeta)>0.$$
*(Exact, `verify.py`.)*

*(4) The critical cone.* The 15 functionals $\ell_q$ have rank 13 (exact row reduction). Their common kernel is $E=\operatorname{span}\{S,U,V\}$, where $S=M'(z)$ is the constant $z$-derivative and
$$U=\begin{pmatrix}0&-1&-1&0\\0&0&0&0\\0&0&0&0\\0&1&0&-1\end{pmatrix},\qquad V=\begin{pmatrix}1&0&0&1\\0&0&-1&-1\\0&0&0&0\\-1&0&1&0\end{pmatrix}.$$

*(5) Second-order condition on the critical cone.* The following is a polynomial identity (exact, `verify.py`):
$$-\det(M(z)+uU+vV)=16(16-12z+12z^2-z^4)+16(z^2-2z+2)u^2+16z^2uv+8(z^2-2z+4)v^2.$$
On $E$ at $(z,u,v)=(\zeta,0,0)$:
- the $zz$-second derivative is $16\cdot(24-12\zeta^2)=192(2-\zeta^2)>0$;
- the $(u,v)$ block has diagonal entries $\ge 16\cdot\frac{10}9$ and $\ge 8\cdot\frac{28}9$, and off-diagonal $|16\zeta^2|<8$, so it is positive definite;
- there are no mixed $z$–$u$ or $z$–$v$ terms at $u=v=0$.

So $D^2f(Y_*)$ is positive definite on $E$.

*(6) Second-order sufficiency.* Suppose not. Then there are admissible $Y_*+H_k$ with $0<\|H_k\|\to0$ and $f(Y_*+H_k)\le f(Y_*)$. Pass to a subsequence with $H_k/\|H_k\|\to H$, $\|H\|=1$.

- By (2), $\ell_q(H_k)\ge0$ for all $q$.
- By (3) and Taylor's theorem, $0\ge f(Y_*+H_k)-f(Y_*)=\sum_q\lambda_q\ell_q(H_k)+O(\|H_k\|^2)$.
- Every $\lambda_q>0$ and every $\ell_q(H_k)\ge0$, so $\ell_q(H_k)=O(\|H_k\|^2)$, and in the limit $\ell_q(H)=0$ for all $q$. Thus $H\in E$.
- Expanding to second order: $0\ge\frac{f(Y_*+H_k)-f(Y_*)}{\|H_k\|^2}\ge\frac12D^2f(Y_*)[H,H]+o(1)$, because the first-order term is $\ge0$.
- The right side is $>0$ for large $k$ by (5).

This is a contradiction. Dividing by 4 transfers the statement to minimum 2. $\square$

---

# 4. Existence and spanning of a critical lattice

**Theorem 4.1 (Mahler) [Cite].** $D_{\mathrm{pack}}$ is attained. A lattice attaining it is called **critical**.

**Lemma 4.2 (spanning) [P].** The minimal vectors of a critical lattice span $\mathbb R^4$.

*Proof.* Suppose $W=\operatorname{span}M(\Lambda)\ne\mathbb R^4$. Choose a projection $\pi$ onto a complement of $W$ (so $\pi|_W=0$), and put $A_\varepsilon=I-\varepsilon\pi$.

- $A_\varepsilon$ fixes $M(\Lambda)$.
- Only finitely many $v\in\Lambda$ have $\|v\|\le3$, and those that are not minimal have norm $\ge2+\eta$ for some $\eta>0$. So for $\varepsilon$ small they keep norm $>2$ after applying $A_\varepsilon$.
- Vectors with $\|v\|>3$ satisfy $\|A_\varepsilon v\|\ge(1-\varepsilon\|\pi\|_{1\to1})\|v\|>2$ for $\varepsilon$ small.

So $\Lambda A_\varepsilon$ is admissible. But $\det(\Lambda A_\varepsilon)=(1-\varepsilon)^{4-\dim W}\det\Lambda<\det\Lambda$, a contradiction. $\square$

---

# 5. Necessary conditions at a local minimum: the ℓ¹ Voronoi conditions

Let $\Lambda_*=\mathbb Z^4B_*$ be a **local** minimum of $\det$ among admissible lattices, oriented so that $\det B_*>0$. Let $S$ be a set of representatives of the minimal pairs, in basis coordinates, and put $v_q=qB_*$.

For each $q\in S$ fix a **supporting facet**: a sign vector $\sigma_q\in\{\pm1\}^4$ with $\sigma_q\cdot v_q=2$. One always exists; it is unique exactly when $v_q$ has no zero coordinate. Put $\ell_q(H)=\sigma_q\cdot qH$.

**Lemma 5.0 (uniform gap) [P].** There is a neighbourhood $N$ of $B_*$ such that $\|qB\|>2$ for all $B\in N$ and all $q\in\mathbb Z^4\setminus(\pm S\cup\{0\})$.

*Proof.* Let $c_0=\min_{\|x\|=1}\|xB_*\|>0$. Only finitely many $q$ have $\|qB_*\|\le 3$; among these, the non-minimal ones have $\|qB_*\|\ge 2+\eta$.

For $B\in N$ we have $\|q(B-B_*)\|\le\|q\|\,\|B-B_*\|\le\|qB_*\|\,\|B-B_*\|/c_0$. Hence $\|qB\|\ge\|qB_*\|(1-\|B-B_*\|/c_0)$, which exceeds $2$ for all relevant $q$ once $N$ is small. $\square$

**Theorem 5.1 (perfection-type rank condition) [P].** For every choice of supporting facets, let $E=\bigcap_{q\in S}\ker\ell_q$. Then
$$E\cap B_*\mathrm{Sym}_4=\{0\},$$
so $\operatorname{rank}\{\ell_q\}\ge 16-6=10$. In particular $|S|\ge 10$: every local minimum has at least ten minimal pairs.

*Proof.* Let $H\in E$ and $|t|$ small. For $q\in S$: $\|q(B_*+tH)\|\ge\sigma_q\cdot q(B_*+tH)=2+t\,\ell_q(H)=2$. The other $q$ are handled by Lemma 5.0. So $B_*+tH$ is admissible for all small $t$ of **both** signs.

Suppose $0\ne H=B_*S'\in E$ with $S'$ symmetric. Then $\det(B_*+tH)=\det B_*\det(I+tS')$, with
- first derivative at $0$: $\operatorname{tr}S'$,
- second derivative at $0$: $(\operatorname{tr}S')^2-\operatorname{tr}(S'^2)$.

A two-sided local minimum forces $\operatorname{tr}S'=0$. The second derivative is then $-\operatorname{tr}(S'^2)<0$, so $\det$ strictly decreases for small $t\ne0$, a contradiction.

Hence $E\cap B_*\mathrm{Sym}_4=\{0\}$. Since $H\mapsto B_*^{-1}H$ is a linear isomorphism, $\dim E\le 16-\dim\mathrm{Sym}_4=6$. $\square$

**Theorem 5.2 (ℓ¹-eutaxy) [P].** For **every** choice of supporting facets there are multipliers $\lambda_q\ge0$ with
$$\sum_{q\in S}\lambda_q\,\sigma_q^{\!\top}v_q=I_4\qquad\text{(a $4\times4$ matrix identity)},\qquad\text{and then}\qquad\sum_q\lambda_q=2.$$

*Proof.* Let $P_\sigma=\{B:\sigma_q\cdot qB\ge2\ \forall q\in S\}$, a polyhedron containing $B_*$.

For $B\in P_\sigma\cap N$ we have $\|qB\|\ge\sigma_q\cdot qB\ge2$ for $q\in S$, and Lemma 5.0 covers the other $q$. So $P_\sigma\cap N$ consists of admissible matrices, and $B_*$ is a local minimum of the smooth function $\det$ on the polyhedron $P_\sigma$.

Linear constraints satisfy Abadie's constraint qualification, so the KKT conditions hold. Since $\nabla\det(B)=\det(B)\,B^{-\top}$ and $\nabla\ell_q=q^{\!\top}\sigma_q$, they read
$$\det(B_*)B_*^{-\top}=\sum_q\mu_q\,q^{\!\top}\sigma_q,\qquad \mu_q\ge0.$$

Transpose, multiply on the right by $B_*$, and divide by $\det B_*>0$. This gives $I=\sum\lambda_q\sigma_q^{\!\top}v_q$ with $\lambda_q=\mu_q/\det B_*\ge0$. Taking traces: $4=\sum\lambda_q\sigma_q\cdot v_q=2\sum\lambda_q$. $\square$

**Check at $\Lambda_\zeta$ [P+C].** The identity of Theorem 3.1(3) is exactly Theorem 5.2 with 15 strictly positive multipliers. Numerically, $\sum\lambda_q\sigma_q^{\!\top}v_q/|\det M(\zeta)|=I_4$ to $10^{-10}$, and $\sum\lambda_q=128-72\zeta+48\zeta^2=|\det M(\zeta)|/2$ (script `maxdet_theory.py`, part 4).

Theorems 5.1–5.2 are the ℓ¹ analogues of Voronoi's perfection and eutaxy. They are the necessary conditions that Stage B of the computation (§10) will solve.

**Theorem 5.3 (Miller) [Cite].** An admissible lattice has at most 20 minimal pairs. The maximum is attained only by lattices $C_4$-equivalent to $D_4^+$ (arXiv:2501.09245, Thm. 2).

**Corollary 5.4 [P].** Every critical lattice has between 10 and 20 minimal pairs. $\Lambda_\zeta$ has 15.

---

# 6. The maximal-determinant frame

This section replaces two things from our earlier drafts: the unproved "good basis" (Descent) theorem, and the saturation hypothesis. Both are superseded by an elementary argument.

Let $\Lambda$ be admissible with spanning minimal vectors (e.g. critical, by Lemma 4.2). Choose minimal $X_1,\dots,X_4$ with $|\det X|$ **maximal** among all quadruples of minimal vectors. Such a maximum exists and is $>0$, since $M(\Lambda)$ is finite and spanning.

Define:
- $p=|\det X|/\det\Lambda$, a positive integer (the index of $\mathbb Z^4X$ in $\Lambda$);
- $\Gamma=\Lambda X^{-1}$, the lattice in **frame coordinates**, so that $\mathbb Z^4\subseteq\Gamma$ and $[\Gamma:\mathbb Z^4]=p$;
- $F(c)=\|cX\|$, a norm with $F(e_i)=2$.

**Theorem 6.1 [P].**
- (i) **(Cramer box)** Every minimal vector has frame coordinates $c\in\Gamma\cap[-1,1]^4$.
- (ii) **(determinant bound)** For any four minimal vectors, $|\det(c^1;c^2;c^3;c^4)|\le1$ in frame coordinates.
- (iii) **(admissibility of the frame lattice)** $\Gamma\cap\{\|c\|_1<1\}=\{0\}$. Every $c\in\Gamma$ with $\|c\|_1=1$ is a minimal vector. And $p\le24$.
- (iv) **($p$-modularity)** In any basis of $\Lambda$, the minimal vectors have integer coordinates, and every $4\times4$ minor lies in $[-p,p]$.

*Proof.*
- (i) Write $v=\sum c_jX_j$. Replacing $X_i$ by $v$ multiplies the determinant by $c_i$ (Cramer). Maximality gives $|c_i|\le1$.
- (ii) This is maximality applied directly.
- (iii) If $\|c\|_1<1$, then $F(c)\le\sum|c_i|F(e_i)<2$, which contradicts admissibility. If $\|c\|_1=1$, then $F(c)\le2$, and admissibility forces equality. Finally, Theorem 1.2 applied to $\Gamma$ and the unit ℓ¹ ball gives $1/p=\det\Gamma\ge\frac{2/3}{16}$.
- (iv) Take a basis $Y$ of $\Lambda$. A minor equals $\det(v^1..v^4)/\det Y=\det(c)\cdot\det X/\det Y$, which has absolute value $|\det c|\,p\le p$. $\square$

**Theorem 6.2 (the 39 strata) [P+C].** Up to signed permutations of the frame, the superlattices $\Gamma\supseteq\mathbb Z^4$ with $[\Gamma:\mathbb Z^4]\le 24$ and $\Gamma\cap\{\|c\|_1<1\}=\{0\}$ form exactly **39 orbits**, with index
$$p\in\{1,2,\dots,10,12,16,18\}\qquad\text{(no admissible $\Gamma$ of index }11,13,14,15,17,19\text{–}24).$$

Per index, the number of orbits and the number of candidate minimal pairs $\#(\Gamma\cap[-1,1]^4\cap\{\|c\|_1\ge1\})/\pm$ are:

| $p$ | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 12 | 16 | 18 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| orbits | 1 | 3 | 2 | 7 | 2 | 4 | 1 | 9 | 2 | 1 | 2 | 4 | 1 |
| candidate pairs | 40 | 48–58 | 56–64 | 64–94 | 72 | 88–98 | 88 | 96–156 | 112–136 | 122 | 146–158 | 160–220 | 208 |

*Proof.* We enumerate $\Gamma$ through its dual in Hermite normal form (201,677 HNFs for $p\le24$) and test the Lee-weight criterion. The criterion is: $\Gamma\cap\{\|c\|_1<1\}=\{0\}$ iff every nonzero residue $r\in p\Gamma/p\mathbb Z^4$ has $\sum_j\min(r_j,p-r_j)\ge p$. This holds because the minimal ℓ¹ norm of the coset $r/p+\mathbb Z^4$ is $p^{-1}\sum_j\min(r_j,p-r_j)$. The orbits were then reduced under the order-384 group. *(Script `strata.py`, exact integer arithmetic.)*

*Cross-check:* the maximum index 18 reproduces Bantegnie's bound (Canad. J. Math. 17 (1965)) for admissible superlattices of $\mathbb Z^4$. $\square$

**Remark.** For $\Lambda_\zeta$ the max-det frame has index $p=4$. All 15 minimal vectors lie in $\Gamma\cap[-1,1]^4\cap\frac14\mathbb Z^4$, and the minors in a basis take exactly the values $\{0,1,2,3,4\}$ (`maxdet_theory.py`, part 2). So $\Lambda_\zeta$ lives in one of the seven $p=4$ strata.

---

# 7. The index-one stratum is completely classified

**Theorem 7.1 [P+C].** Let $\Lambda$ be a local minimum of $\det$ whose max-det frame has $p=1$. Then $\Lambda$ has **exactly ten** minimal pairs, and there is a basis $X_1,\dots,X_4$ of $\Lambda$ with
$$M(\Lambda)=\{\pm X_i\}\cup\{\pm(X_i-X_j)\}.$$
Equivalently, with $X_0=0$, the five points $X_0,\dots,X_4$ are pairwise at ℓ¹-distance exactly 2 (an ℓ¹-equilateral simplex). The minimal pairs form the positive roots of $A_4$, whose matroid is the graphic matroid $M(K_5)$.

*Proof.* With $p=1$ we have $\Gamma=\mathbb Z^4$. By Theorem 6.1(i)–(ii), the minimal vectors have ternary frame coordinates, and every quadruple of them has determinant in $\{0,\pm1\}$. So the minimal pairs form a family of pairwise non-parallel ternary vectors (non-parallel by Lemma 1.5). The family contains $e_1,\dots,e_4$, and all its $4\times4$ minors lie in $\{0,\pm1\}$.

We enumerated all such families exhaustively, by backtracking over the 36 remaining ternary pairs with the minor test applied at every insertion (`maxdet_theory.py`, parts 1, 3, 5):
- there are 2721 maximal families, of sizes 5–10;
- exactly 200 have size 10, and they form 3 orbits under signed permutations (one for each spanning-tree type of $K_5$ used as the frame: path, star, fork);
- each of the 200 is the set of positive roots of $A_4$ in a suitable basis (verified directly).

So $|S|\le10$. By Theorem 5.1, $|S|\ge10$. Hence $|S|=10$ and $S$ is the $A_4$ configuration.

In the basis of simple roots the positive roots are the consecutive sums $\alpha_i+\dots+\alpha_j$. Putting $X_k=\alpha_1+\dots+\alpha_k$ converts them to $X_j-X_i$ ($0\le i<j\le4$). $\square$

**Remark (literature).** The bound $|S|\le10$ is Heller's theorem (Pacific J. Math. 7 (1957)): a simple unimodular system of rank $n$ has at most $n(n+1)/2$ elements. Uniqueness of the extremal system is the classical characterization of $M(K_{n+1})$ as the unique maximum-size simple regular matroid. Our exhaustive check makes the case $n=4$ self-contained.

**Consequence.** The whole $p=1$ stratum is **a single rigid combinatorial problem**: minimize $|\det X|$ over ℓ¹-equilateral 5-point sets $\{0,X_1,..,X_4\}$ of side 2 that generate an admissible lattice.

**[Num]** A search with exact shortest-vector verification (263 restarts, `p1probe.py`) finds a minimum of about $0.8363$, which is **3.5% above $D_*$**.

---

# 8. Slicing (a certified pruning tool) [P]

**Theorem 8.1.** For admissible $\Lambda$ and primitive $w\in\Lambda^*$:
$$\det\Lambda\ \ge\ \frac{\Delta_3\bigl(K\cap\hat w^\perp\bigr)}{\|w\|_2}\ \ge\ \frac{\operatorname{vol}_3(C_4\cap\hat w^\perp)}{\|w\|_2}.$$

*Proof.* Primitivity makes $\langle w,\cdot\rangle:\Lambda\to\mathbb Z$ surjective. So $\Lambda$ is a union of translates of $\Lambda_0=\Lambda\cap w^\perp$, spaced $1/\|w\|_2$ apart, and $\det\Lambda=\det\Lambda_0/\|w\|_2$.

$\Lambda_0$ inherits admissibility for the section body. Its critical determinant is $\ge\operatorname{vol}_3/8$ by Theorem 1.2 in dimension 3, and $\operatorname{vol}_3(K\cap\hat w^\perp)=8\operatorname{vol}_3(C_4\cap\hat w^\perp)$. $\square$

**Exact section volumes (Hermite–Genocchi).** Write $C_4$ as the union of its 16 orthant simplices. For a **unit** vector $u$,
$$\operatorname{vol}_3(C_4\cap u^\perp)=\sum_{\sigma}\tfrac{1}{4!}\cdot 4\cdot\bigl[(y)_+^{3}\bigr][0,\sigma_1u_1,..,\sigma_4u_4]$$
(a divided difference). Verified exactly: $s(e_4)=\frac43$, $s(1100)=\frac{2\sqrt2}3$, $s(1110)=\frac{\sqrt3}2$, $s(1111)=\frac56$. The formula replaces floating-point hull volumes, so the prune becomes an exact certificate.

**[Num] Tightness.** At $\Lambda_\zeta$ with $w=(-1,0,0,1)B^*$, the slice lattice (6 contacts, a 3-parameter contact family) is optimal for its section in every test run. Hence $\Delta_3/\|w\|=0.808220001=D_*$: the inequality is sharp at the conjectured optimum.

Generic admissible lattices give ratio $\approx0.80$. So slicing prunes $\det\gtrsim1.0$ (trivial $\Delta_3$), or $\gtrsim0.94$ (sharp $\Delta_3$). It is a filter, not a proof.

---

# 9. The Reduction Theorem

A **datum** $\delta$ consists of:
1. a stratum $\Gamma$ from Theorem 6.2;
2. a set $S\subseteq\Gamma\cap[-1,1]^4$ of representatives of minimal pairs, with $10\le|S|\le20$, containing $e_1,..,e_4$ and every $c\in\Gamma$ with $\|c\|_1=1$, and with all quadruple determinants in $[-1,1]$;
3. **face data**: for each column index $j$, a face $\phi_j$ of the central hyperplane arrangement $\{c^\perp:c\in S\}\subset\mathbb R^4$.

Face data exists because $(cX)_j=\langle c,x^{(j)}\rangle$, where $x^{(j)}$ is the $j$-th column of $X$. So the signs of **all** contact vectors are determined by which face of the arrangement each of the four columns lies in. Given the faces, every $\|cX\|$ with $c\in S$ is a linear function of $X$.

For a datum $\delta$, let $V_\delta$ be the set of real $X$ satisfying all of the following:
- **contacts**: $\|cX\|=2$ for $c\in S$ (linear, given the face data);
- **face data**: column $j$ lies in the relative interior of $\phi_j$;
- **ℓ¹-eutaxy**: $\sum_{c\in S}\lambda_c\sigma_c^{\!\top}(cX)=I$ for some $\lambda\ge0$ (Theorem 5.2 in frame form);
- **perfection-type rank**: the contact functionals have rank $\ge10$ on $X\cdot\mathrm{Sym}$ (Theorem 5.1).

**Theorem 9.1 (Reduction) [P].**
1. The set $\Phi$ of data is finite and explicitly enumerable.
2. For each $\delta$, the function $|\det X|/p$ takes **finitely many values** on $V_\delta$. These values are real algebraic numbers.
3. $D_{\mathrm{pack}}$ is the minimum, over $\delta\in\Phi$, of those values attained at points $X\in V_\delta$ for which $\Gamma X$ is admissible. Admissibility is decidable by the finite test of Theorem 1.7 applied to $\Gamma$ with the frame $X$.

Consequently $D_{\mathrm{pack}}$ is a real algebraic number, and **Theorem Z is equivalent to the finite Computational Claim C** (§10).

*Proof.*
- **Finiteness of $\Phi$.** There are 39 strata; $\Gamma\cap[-1,1]^4$ is finite; and an arrangement has finitely many faces.
- **The critical lattice appears.** Let $\Lambda$ be critical (Theorem 4.1). Take its max-det frame; by Theorem 6.1 it lands in a stratum $\Gamma$. Its minimal pairs form an admissible $S$: the count is $10\le|S|\le20$ by Corollary 5.4, and the remaining properties of item 2 hold by Theorem 6.1. Its columns determine face data. By Theorems 5.1–5.2, its frame matrix lies in $V_\delta$. So $D_{\mathrm{pack}}$ is among the values in (3).
- **Every value is an upper bound.** Each admissible value in (3) is the determinant of an admissible lattice, so it is $\ge D_{\mathrm{pack}}$.
- **Finitely many values.** $V_\delta$ is semialgebraic. On the affine space cut out by the contacts, every point of $V_\delta$ is a KKT point of the polynomial $\det$ for a linearly constrained problem. By the semialgebraic Morse–Sard theorem (Bochnak–Coste–Roy, *Real Algebraic Geometry*, 9.6), the set of such critical values is finite. $\square$

---

# 10. The computational part: exact specification

**Computational Claim C [Conj].** For every datum $\delta$ and every $X\in V_\delta$ with $\Gamma X$ admissible, $|\det X|/p\ge D_*$. Equality holds only for the datum of $\Lambda_\zeta$ (in a $p=4$ stratum), at the signed-permutation images of $\Lambda_\zeta$.

**Already done (exact, reproducible):**
- C0.1: admissibility of the $\zeta$-family and the determinant identity (`verify.py`);
- C0.2: the local-optimality identities, rank 13, and positive multipliers (`verify.py`);
- C0.3: the 39 strata (`strata.py`);
- C0.4: the $p=1$ classification (`maxdet_theory.py`);
- C0.5: ℓ¹-eutaxy at $\zeta$ (numerical to $10^{-10}$; exact through C0.2).

**Expected outcomes [Num]:**
- $p=1$: minimum $\approx0.8363$ (+3.5%);
- the $\zeta$ datum is recovered in its $p=4$ stratum at exactly $D_*$;
- every other stratum is expected to give values $>D_*$.

The strength of this expectation comes from:
- the tightness of slicing at $\zeta$;
- the 15-contact maximality of $\zeta$ (no contact can be added);
- the absence of competitors in every search so far.

**Stages (in order, each with a validation gate):**

**Stage A — combinatorial enumeration of $(\Gamma,S)$.**
- *Task:* for each of the 39 strata, backtrack over the 40–220 candidate pairs. Prune with the quadruple-determinant bound (§6(ii)), the rank-$\ge10$ requirement, and the upper bound 20.
- *Gates:* the $p=1$ output must be exactly the $A_4$ family (known answer); $\zeta$'s $S$ must appear in its $p=4$ stratum.
- *Size:* unknown. It must be measured before Stage B is designed.
- *Cost:* 2–3 days to implement; hours of cluster time.

**Stage B — geometry for each $(\Gamma,S)$.**
- *Task:* enumerate face data, then solve the Voronoi system. On the contact-defined affine space (dimension $\le6$ by Theorem 5.1), the eutaxy condition says the gradient of $\det$ lies in the cone of the contact functionals. This is a polynomial system of degree $\le3$ in $\le6$ unknowns.
- *Method:* exact Gröbner bases or rational univariate representations (e.g. msolve); or certified homotopy continuation plus a posteriori certification.
- *Then:* keep the real solutions with $\lambda\ge0$ that lie in the prescribed faces; test admissibility by Theorem 1.7; compare $|\det|/p$ with $D_*$ **exactly**, as algebraic numbers.
- *Gates:* recover $\Lambda_\zeta$ exactly; reproduce $\approx0.8363$ for $p=1$.
- *Cost:* 1–2 weeks to implement; runtime depends on the Stage A count.

**Stage C — certification and write-up.** Exact certificates (a rational representation of each critical value, and an exact comparison with $D_*$), plus the paper. Formalization is a separate project; the natural first Lean targets are the rational counterexample and Theorem 3.1.

**Recommended first milestone: the $p=1$ pilot.** One combinatorial type, three frame embeddings, and at most a 6-dimensional problem. It exercises Stages A–C end to end on a stratum with a known expected value, and it yields the first rigorous exclusion theorem: *no critical lattice lies in the $p=1$ stratum*.
- Estimate: **1–2 weeks**.

---

# 11. What remains, and how long it should take

| item | nature | estimate |
|---|---|---|
| $p=1$ pilot (Stages A–C on one stratum) | computation + certificate | 1–2 weeks |
| Stage A on all 39 strata | computation | ~1 week incl. implementation |
| Stage B on all data | computation, the dominant risk | 2–6 weeks, depending on the Stage A count and the number of face data |
| Stage C and the paper | certification + writing | 3–6 weeks |
| **Total to a complete computer-assisted proof of Theorem Z** | | **≈ 2–4 months**, conditional on a moderate Stage A count |
| Lean formalization | separate project | months |

**Main risk.** The number of face data per $S$ is polynomial in $|S|$ (faces of an arrangement of $\le20$ hyperplanes in $\mathbb R^4$), not the $16^{|S|}$ that killed `critdet.py`. It could still be large for $p\ge8$, where the candidate pools reach 96–220 pairs. If so, the fallback is to solve the Voronoi system first and read the face data off the solutions, with certified path tracking guaranteeing completeness.

**Theory that would shrink the computation, not attempted here:**
- sharper per-$p$ upper bounds on $|S|$ from the theory of $\Delta$-modular matrices (Heller gives 10 for $p=1$; bounds for $\Delta=2,3,\dots$ exist in the literature and must be checked for applicability);
- an analytic lower bound for the $p=1$ equilateral problem.

---

# 12. Mistakes this plan is built to avoid

1. **No type exclusion.** Types are configurations, not disjoint lattice classes ($L_4$ itself contains clean frames of index 1–4). The strata here are defined by the *maximal*-determinant frame of the lattice under study, and nothing is "excluded by type".
2. **No unproved descent and no saturation.** Theorem 6.1 gives the coordinate box outright. Theorem 7.1 shows that for $p=1$ the configuration is forced, with no parity-class assumption.
3. **Complete admissibility tests only.** Use true shortest-vector enumeration, or Theorem 1.7 including its $\mu\le1$ points and the contact condition on the frame. Never use a fixed coefficient box, and never use pinned-sign relaxations to bound $\det$ from below.
4. **Never minimize $|\det|$ over sets that contain singular matrices.** Theorem 1.7 and Corollary 1.8 exclude them by construction.
5. **Minima are not at vertices.** $\zeta$ is a cubic irrational. Stage B solves critical-point systems exactly instead of enumerating vertices.
6. **Validate on a known answer before any cluster run.** Every stage has a gate ($A_4$ for $p=1$; exact recovery of $\zeta$).
7. **Do not read silence as a stall.** Long cluster jobs must emit heartbeats.

---

## References to confirm before submission
- Mahler, selection theorem; Minkowski; van der Corput (standard: Cassels; Gruber–Lekkerkerker).
- Heller, *On linear systems with integral valued solutions*, Pacific J. Math. **7** (1957) 1351–1364.
- Uniqueness of $M(K_{n+1})$ among maximum-size simple regular matroids (e.g. Oxley, *Matroid Theory*).
- Miller, *On the kissing number of the cross-polytope*, arXiv:2501.09245.
- Bantegnie, Canad. J. Math. **17** (1965) 725–730.
- Swinnerton-Dyer, Proc. Cambridge Philos. Soc. **49** (1953) 161–162 (classical form of Theorem 5.1).
- Bochnak–Coste–Roy, *Real Algebraic Geometry*, §9.6.
- The external audit note (source of Theorem 3.1's data and of the finite-test theorem in the clean case).
