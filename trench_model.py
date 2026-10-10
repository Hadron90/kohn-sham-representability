"""
Strip-model numerics for the manuscript (Section 4; Tables 1 and 2).  Authoritative script.

Model (Section 4.3): N(z) = z^2 + eps^2;  angle theta with N theta' = J (vertical flux),
flux law J' = -d N sin(2 theta);  theta odd, theta(0) = 0, theta(+inf) = pi/2.
Profiles G = sqrt(N) cos(theta), H = sqrt(N) sin(theta);  W := G''/(2G).
Formula (Prop. 4.1(b), exact):   W - d = (eps^2 - J^2)/(2 N^2) - d cos^2(theta).

Method.  The system is integrated in the complementary angle phi = pi/2 - theta, i.e.
    phi' = -J/N,   J' = -d N sin(2 phi),   phi(0) = pi/2,   J(0) = s eps^2,
so that the solver's relative tolerance applies to phi itself where phi is small (integrating theta
directly limits the accuracy of pi/2 - theta to the absolute error of theta).  Shooting on
s = theta'(0): integrate on [0, Z_bc], Z_bc = 12/sqrt(2d), with DOP853 (rtol = 1e-12, atol = 1e-30).
The global branch is selected by the far-field decay condition.  The linearisation of the angle
equation for z >> eps is (z^2 phi')' = 2 d z^2 phi, whose solutions are exp(-+kappa z)/z,
kappa = sqrt(2d); the decaying one satisfies J = -N phi' = N phi (kappa + 1/z).  We bisect on s until
    r_bc(s) := J(Z_bc) - N(Z_bc) phi(Z_bc) (kappa + 1/Z_bc)
changes sign (r_bc > 0 <-> growing mode with negative coefficient <-> overshoot).  All reported
quantities are evaluated on [-Z, Z] with Z = 8/sqrt(2d).  The relative residual r of the decay
condition at Z (Table 2) is a diagnostic: in the linearised two-mode description
phi = (A e^{-kappa z} + B e^{kappa z})/z with N ~ z^2, the growing-to-decaying amplitude ratio at
kappa Z = 8 is -9R/(16 - 7R) for the signed residual R, i.e. of magnitude 9r/16 + O(r^2).  The computed solution is NOT a validated
enclosure of the infinite-domain solution; every quantity below is an a posteriori estimate.

Quantities reported for each (d, eps):
  Table 1 (paper/table_strip.tex):
    s*eps              : theta'(0) eps = J0/eps
    W0*eps/sqrt(2d)    : compared with -2/pi, the value predicted by the inner profile (S2)
    int(W-d), int|W-d| : integrals over [-Z, Z] (trapezoid on a graded mesh)
    -sqrt(2d)          : the conjectured limit of int(W-d)
    int(W-d)_+         : integral of the positive part over [-Z, Z]  (int|W-d| = -int(W-d) + 2 int(W-d)_+)
    tail               : a posteriori estimate of int_{|z|>Z} |W-d| dz, obtained by inserting the
                         computed (theta(Z), J(Z)) into the inequalities of Prop. 4.1(c),(d):
                           |W-d| <= (eps^2 + J(Z)^2)/(2 z^4) + d cos^2(theta)   (z >= Z),
                           int_Z^inf d cos^2(theta) <= J(Z)/(sqrt(2) N(Z))      (theta >= pi/4 beyond Z),
                         so  tail = 2 [ (eps^2 + J(Z)^2)/(6 Z^3) + J(Z)/(sqrt(2) N(Z)) ],
                         which requires J(Z) >= 0 and theta(Z) >= pi/4 (both asserted).
    conv               : max(|int|W-d| - same with rtol = 1e-10|, |... - same with Z_bc = 16/sqrt(2d)|)
  Table 2 (paper/table_diag.tex):
    J0/eps - 1         : the cancellation of the eps^-2 ridge (relative accuracy of J0 = eps)
    ratio              : (J0/eps - 1) / ((2/pi) sqrt(2d) eps), compared with 1
    K_eps              : max_{|z|<=1} |eps^2 - J^2| / (eps N)   (normalised flux ratio, Lemma 5.8)
    surf               : max_{z in [1/kappa, 3/kappa]} | G(z) e^{kappa z}/G(0) - 1 |
    prof               : max_{s in [0,3]} | eps (W-d)(eps s)/(-(2/pi) sqrt(2d)) - (1+s^2)^-2 |
    J(Z)/eps, pi/2-theta(Z) : end values at Z (both must be positive for the global branch)
    r                  : |J(Z) - N(Z)(pi/2-theta(Z))(kappa+1/Z)| / J(Z)
  Implementation checks (printed only): fd_err, the max relative error of W against a
  finite-difference G''/(2G) on the inner mesh; flux_err, the max relative difference between the
  two algebraically equivalent formulas for W.

Usage:  python trench_model.py            # human-readable report
        python trench_model.py --table    # also writes paper/table_strip.tex and paper/table_diag.tex
"""
import sys, os
import numpy as np
from scipy.integrate import solve_ivp, trapezoid

def solve(d, eps, Zfac=12.0, rtol=1e-12):
    kap = np.sqrt(2*d)
    Zbc = Zfac/kap
    def rhs(z, y):
        ph, J = y
        N = z*z + eps*eps
        return [-J/N, -d*N*np.sin(2*ph)]
    def shoot(s):
        return solve_ivp(rhs, [0, Zbc], [np.pi/2, s*eps*eps], method='DOP853',
                         rtol=rtol, atol=1e-30, dense_output=True)
    def resid(sol):
        phZ, JZ = sol.y[0, -1], sol.y[1, -1]
        return JZ - (Zbc*Zbc + eps*eps)*phZ*(kap + 1.0/Zbc)
    lo, hi = 0.5/eps, 2.0/eps
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if resid(shoot(mid)) > 0:
            hi = mid
        else:
            lo = mid
        if hi - lo < 4e-16*hi:
            break
    s = 0.5*(lo+hi)
    return s, shoot(s), Zbc

def W_of(sol, z, d, eps):
    """Returns cos(theta), sin(theta), J, N, W (eq. (4.7) form), W (flux-identity form)."""
    ph, J = sol.sol(z)
    N = z*z + eps*eps
    costh, sinth = np.sin(ph), np.cos(ph)
    W = 0.5*eps*eps/N**2 - 0.5*(J/N)**2 + d*sinth**2
    Wflux = d + (eps*eps - J*J)/(2*N**2) - d*costh**2
    return costh, sinth, J, N, W, Wflux

def integrals(sol, Z, d, eps):
    z = np.concatenate([np.linspace(0, 20*eps, 4001), np.geomspace(20*eps, Z, 6000)[1:]])
    costh, sinth, J, N, W, Wflux = W_of(sol, z, d, eps)
    return dict(strength=2*trapezoid(W - d, z), absint=2*trapezoid(np.abs(W - d), z),
                pos=2*trapezoid(np.maximum(W - d, 0.0), z),
                flux_err=np.max(np.abs(W - Wflux))/np.max(np.abs(W)),
                z=z, W=W, costh=costh, J=J, N=N)

def analyze(d, eps):
    s, sol, Zbc = solve(d, eps)
    kap = np.sqrt(2*d)
    Z = 8.0/kap                      # reporting domain (see docstring)
    I = integrals(sol, Z, d, eps)
    z, W, costh, J, N = I['z'], I['W'], I['costh'], I['J'], I['N']
    W0 = W[0]
    # finite-difference check of W = G''/(2G) on the inner mesh
    zf = z[:4001]; Gf = np.sqrt(zf*zf+eps*eps)*costh[:4001]
    Gpp = np.gradient(np.gradient(Gf, zf), zf)
    fd_err = np.max(np.abs(Gpp/(2*Gf) - W[:4001])[50:-50])/np.max(np.abs(W[:4001][50:-50]))
    # end values and a posteriori tail estimate (hypotheses of Prop. 4.1(c),(d) asserted)
    phZ, JZ = sol.sol(Z)
    assert JZ > 0 and phZ > 0 and phZ <= np.pi/4, (d, eps, JZ, phZ)
    NZ = Z*Z + eps*eps
    tail = 2*((eps*eps + JZ*JZ)/(6*Z**3) + JZ/(np.sqrt(2)*NZ))
    r = abs(JZ - NZ*phZ*(kap + 1.0/Z))/JZ
    # normalised flux ratio on |z| <= 1 (Lemma 5.8)
    m1 = z <= 1.0
    Keps = np.max(np.abs(eps*eps - J[m1]**2)/(eps*N[m1]))
    # surface-state check on [1/kappa, 3/kappa]
    zs = np.linspace(1/kap, 3/kap, 401)
    Gs = np.sqrt(zs*zs+eps*eps)*np.sin(sol.sol(zs)[0])
    surf = np.max(np.abs(Gs*np.exp(kap*zs)/eps - 1))
    # inner-profile check on s in [0, 3]
    ss = np.linspace(0, 3, 301)
    Wp = W_of(sol, eps*ss, d, eps)[4]
    prof = np.max(np.abs(eps*(Wp - d)/(-(2/np.pi)*kap) - (1+ss*ss)**-2))
    # stability checks
    s2, sol2, _ = solve(d, eps, rtol=1e-10)
    absint_rtol = integrals(sol2, Z, d, eps)['absint']
    s3, sol3, _ = solve(d, eps, Zfac=16.0)
    absint_Z = integrals(sol3, Z, d, eps)['absint']
    return dict(d=d, eps=eps, s_eps=s*eps, W0=W0, ratio=W0*eps/kap, strength=I['strength'],
                absint=I['absint'], pos=I['pos'], pred=-kap, tail=tail, JZ_eps=JZ/eps, phiZ=phZ, r=r,
                conv_rtol=abs(I['absint']-absint_rtol), conv_Z=abs(I['absint']-absint_Z),
                fd_err=fd_err, flux_err=I['flux_err'], Z=Z,
                canc=s*eps-1, canc_ratio=(s*eps-1)/((2/np.pi)*kap*eps), Keps=Keps, surf=surf, prof=prof)

CASES = [(0.1, 1e-1), (0.1, 1e-2), (0.1, 1e-3), (0.01, 1e-3), (0.001, 1e-3)]

def fmt_e(x):
    m, e = f"{x:.1e}".split("e")
    return f"{m}\\cdot10^{{{int(e)}}}"

def write_table(path, spec, rows):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("% generated by trench_model.py --table ; do not edit by hand\n")
        f.write("\\begin{tabular}{@{}" + "c"*len(spec.split("&")) + "@{}}\n\\toprule\n" + spec + "\\\\\n\\midrule\n")
        f.write("\n".join(rows) + "\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    print("wrote", path)

if __name__ == "__main__":
    rows1, rows2 = [], []
    for d, eps in CASES:
        r = analyze(d, eps)
        lead = f"${d}$ & $10^{{{int(np.log10(eps))}}}$ & "
        rows1.append(lead + f"${r['s_eps']:.5f}$ & ${r['ratio']:+.4f}$ & "
                     f"${r['strength']:+.6f}$ & ${r['absint']:.6f}$ & ${r['pred']:+.6f}$ & "
                     f"${fmt_e(r['pos'])}$ & ${fmt_e(r['tail'])}$ & ${fmt_e(max(r['conv_rtol'], r['conv_Z']))}$ \\\\")
        rows2.append(lead + f"${fmt_e(r['canc'])}$ & ${r['canc_ratio']:.4f}$ & ${r['Keps']:.4f}$ & "
                     f"${fmt_e(r['surf'])}$ & ${fmt_e(r['prof'])}$ & ${fmt_e(r['JZ_eps'])}$ & "
                     f"${fmt_e(r['phiZ'])}$ & ${fmt_e(r['r'])}$ \\\\")
        print(f"d={d:<6} eps={eps:<6} Z={r['Z']:.2f} theta'(0)eps={r['s_eps']:.6f} "
              f"W0*eps/sqrt(2d)={r['ratio']:+.5f} (-2/pi={-2/np.pi:+.5f}) int(W-d)={r['strength']:+.6f} "
              f"int|W-d|={r['absint']:.6f} -sqrt(2d)={r['pred']:+.6f} int(W-d)+={r['pos']:.2e} "
              f"tail~{r['tail']:.1e} J(Z)/eps={r['JZ_eps']:.2e} pi/2-theta(Z)={r['phiZ']:.2e} r={r['r']:.1e} "
              f"conv(rtol)={r['conv_rtol']:.1e} conv(Z)={r['conv_Z']:.1e} "
              f"fd_err={r['fd_err']:.1e} flux_err={r['flux_err']:.1e} | "
              f"J0/eps-1={r['canc']:.6e} ratio={r['canc_ratio']:.4f} K_eps={r['Keps']:.4f} "
              f"surf={r['surf']:.6e} surf/(J0/eps-1)={r['surf']/r['canc']:.4f} prof={r['prof']:.2e}")
    if "--table" in sys.argv:
        here = os.path.join(os.path.dirname(os.path.abspath(__file__)), "paper")
        write_table(os.path.join(here, "table_strip.tex"),
                    r"$d$ & $\varepsilon$ & $\theta'(0)\,\varepsilon$ & $\dfrac{W(0)\varepsilon}{\sqrt{2d}}$ & "
                    r"$\int_{-Z}^{Z}(W-d)$ & $\int_{-Z}^{Z}|W-d|$ & $-\sqrt{2d}$ & $\int_{-Z}^{Z}(W-d)_+$ & tail (est.) & conv", rows1)
        write_table(os.path.join(here, "table_diag.tex"),
                    r"$d$ & $\varepsilon$ & $J_0/\varepsilon-1$ & $\dfrac{J_0/\varepsilon-1}{(2/\pi)\sqrt{2d}\,\varepsilon}$ & "
                    r"$K_\varepsilon$ & surf & prof & $J(Z)/\varepsilon$ & $\tfrac{\pi}{2}-\theta(Z)$ & $r$", rows2)
