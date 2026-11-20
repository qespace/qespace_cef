"""Fit CEF parameters to several data sets at once: a mock of the v2 workflow.

Goal: show how a fit could look with the broadcasting API of AUDIT A15
(rule (b) of Q28, parameters as a labelled array): one model function that
works for one parameter set (the fit) and for a whole grid of them (the scan
for a start value), with one batched diagonalisation per data set.

Status: mock. The v2 API does not exist yet; the section "MOCK of the v2
API" below stands in for it, written on today's public StevensBase and J_x,
J_y, J_z. Its observables use the correct formulas of the tests in
tests/test_observables.py, not today's package functions (AUDIT P1, P3, P10).
The data are synthetic, made from known parameters plus noise.

Packages: qespace_cef, numpy, lmfit (not a project dependency):
    uv run --no-project --isolated --with-editable . --with lmfit python examples/fit_cef_mock.py

Units (OPEN, AUDIT C1, C11): energies in meV, T in K, fields in T, moments in
muB per ion, susceptibility in muB/T per ion.
"""

import time
from functools import lru_cache
from typing import NamedTuple

import lmfit
import numpy as np

import qespace_cef as cef # always import as cef
from qespace_cef.cefmatrices import J_x, J_y, J_z, StevensBase

# CODATA 2018: k_B = 8.617333262e-5 eV/K, mu_B = 5.7883818060e-5 eV/T.
KB = 0.08617333262  # meV/K
MUB = 0.057883818060  # meV/T


# =============================================================================
# MOCK of the v2 API (A15).
# Shapes: S = loop shape (broadcast of all inputs), n = 2J+1, K = parameters.
# =============================================================================


class Pars(NamedTuple):
    """Stands for CEFpars in v2: names for the user, a labelled array inside."""

    names: tuple  # K names, e.g. ("B20", "B40", "B44")
    values: np.ndarray  # S_p + (K,), meV

    @classmethod
    def from_dict(cls, Bdict):
        """Each value a number or an array; they broadcast together."""
        arrays = np.broadcast_arrays(*(np.asarray(v, float) for v in Bdict.values()))
        return cls(tuple(Bdict), np.stack(arrays, axis=-1))

    def new_axis(self):
        """The same parameters with one more loop axis at the end: S_p + (1, K)."""
        return Pars(self.names, self.values[..., None, :])


@lru_cache
def _operators(J, names):
    """Stevens operators in the order of ``names``, (K, n, n), and J, (3, n, n)."""
    ops = StevensBase(J)
    stevens = np.array([ops[name[1:]] for name in names]) if names else None
    return stevens, np.array([J_x(J), J_y(J), J_z(J)])

# this should go into hamiltonian after it is solved. 
class Eigensystem(NamedTuple):
    energies: np.ndarray  # S + (n,), meV, ascending, lowest at 0
    vectors: np.ndarray  # S + (n, n), columns are eigenvectors
    J: np.ndarray  # (3, n, n), free-ion J; observables rotate only what they need
    gJ: float


# I want this to be a class
class Hamiltonian():
    def __init__(self, ion, pars, field=(0.0, 0.0, 0.0)):
        """H = sum_k B_k O_k - muB gJ field.J.

        pars.values S_p + (K,), field S_f + (3,) -> H broadcast(S_p, S_f) + (n, n).
        """
        stevens, J = _operators(ion.J, pars.names)
        H_cef = np.einsum("...k,kij->...ij", pars.values, stevens)
        H_zeeman = -MUB * ion.gJ * np.einsum("...a,aij->...ij", np.asarray(field, float), J)
        self.H = H_cef + H_zeeman


    def solve(self, ion, H):
        """Batched diagonalisation: H S + (n, n) -> Eigensystem with loop shape S."""
        w, V = np.linalg.eigh(H)
        return Eigensystem(w - w[..., :1], V, _operators(ion.J, ())[1], ion.gJ)





def _rotate(sol, op):
    """<i|op|j> in the eigenbasis: op (n, n) of the free ion -> S + (n, n)."""
    return np.einsum("...ki,kl,...lj->...ij", sol.vectors.conj(), op, sol.vectors)


def _rotate_J(sol):
    """<i|J_a|j> in the eigenbasis: S + (3, n, n)."""
    return np.einsum("...ki,akl,...lj->...aij", sol.vectors.conj(), sol.J, sol.vectors)


def _populations(energies, T):
    """Boltzmann populations: energies S + (n,), T broadcasting against S -> broadcast(S, T) + (n,)."""
    p = np.exp(-energies / (KB * np.asarray(T, float)[..., None]))
    return p / p.sum(-1, keepdims=True)


def magnetization(sol, T):
    """M = gJ sum_n p_n Re<n|J|n>: T against S -> broadcast(S, T) + (3,), muB per ion."""
    J_nn = np.einsum("...in,aij,...jn->...na", sol.vectors.conj(), sol.J, sol.vectors).real
    return sol.gJ * np.einsum("...n,...na->...a", _populations(sol.energies, T), J_nn)


def susceptibility(sol, T, direction):
    """Zero-field chi along ``direction`` (3,): T against S -> broadcast(S, T), muB/T per ion.

    Van Vleck and Curie terms; degenerate levels within 1e-9 of the level
    scale, a basis-independent mask (CLAUDE.md, "Degenerate levels").
    """
    d = np.asarray(direction, float) / np.linalg.norm(direction)
    Jn2 = np.abs(_rotate(sol, np.einsum("a,aij->ij", d, sol.J))) ** 2  # S + (n, n)
    E = sol.energies[..., None, :] - sol.energies[..., :, None]  # S + (n, n), E_m - E_n
    scale = np.maximum(np.abs(sol.energies).max(-1), 1.0)[..., None, None]
    degenerate = np.abs(E) <= 1e-9 * scale
    p = _populations(sol.energies, T)  # broadcast(S, T) + (n,)
    kT = KB * np.asarray(T, float)[..., None, None]
    curie = np.where(degenerate, p[..., :, None] / kT, 0.0)
    vleck = np.where(
        degenerate, 0.0, (p[..., :, None] - p[..., None, :]) / np.where(degenerate, 1.0, E)
    )
    return sol.gJ**2 * MUB * np.sum((curie + vleck) * Jn2, axis=(-2, -1))


def kramers_peaks(sol, T):
    """Powder neutron peaks out of the ground doublet of a Kramers ion in zero field.

    Levels come in pairs (2k, 2k+1). Peak k has the energy of pair k and the weight
    sum over i in pair 0 and j in pair k of p_i 2/3 sum_a |<i|J_a|j>|^2 (the form
    factor is the same for all peaks at one Q and is left out).
    Returns energies S + (n/2 - 1,) in meV and weights S + (n/2 - 1,), arbitrary units.
    """
    W = 2 / 3 * np.sum(np.abs(_rotate_J(sol)) ** 2, axis=-3)  # S + (n, n)
    W = _populations(sol.energies, T)[..., :, None] * W
    m = sol.energies.shape[-1] // 2
    from_ground = W[..., :2, :].sum(-2).reshape(*W.shape[:-2], m, 2).sum(-1)  # S + (m,)
    return sol.energies[..., 2::2], from_ground[..., 1:]


# =============================================================================
# Synthetic data: Ce3+ in a tetragonal field.
# =============================================================================

ion = qc.Ion("Ce")
TRUE = dict(B20=0.5, B40=0.01, B44=0.1)  # meV, mock values for the synthetic data

T_ins = 5.0  # K
T_chi = np.geomspace(2, 300, 30)  # (30,)
H_mag = np.zeros((20, 3))  # (20, 3), along c
H_mag[:, 2] = np.linspace(0.5, 10, 20)
T_mag = 2.0  # K


def model(pars):
    """All observables for parameters with loop shape S_p (a point, or a grid).

    Each data set gets one more loop axis at the end, so the results have the
    shape S_p + (data points,): (2,) peaks, (30,) chi, (20,) M(H). Two batched
    diagonalisations: zero field, S_p + (1,), and the 20 fields, S_p + (20,).
    """
    p = pars.new_axis()  # S_p + (1, K)
    sol0 = solve(ion, hamiltonian(ion, p))  # S_p + (1,)
    E, W = kramers_peaks(sol0, T_ins)  # S_p + (1, 2)
    solH = solve(ion, hamiltonian(ion, p, H_mag))  # S_p + (20,)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = W[..., 0, 1:] / W[..., 0, :1]  # S_p + (1,)
    return dict(
        E=E[..., 0, :],
        ratio=ratio,
        chi_c=susceptibility(sol0, T_chi, [0, 0, 1]),  # S_p + (30,)
        chi_a=susceptibility(sol0, T_chi, [1, 0, 0]),  # S_p + (30,)
        M_c=magnetization(solH, T_mag)[..., 2],  # S_p + (20,)
    )


rng = np.random.default_rng(0)
truth = model(Pars.from_dict(TRUE))
data = {  # (values, standard deviations), each with the shape of its data axis
    "E": (truth["E"] + rng.normal(0, 0.05, 2), np.full(2, 0.05)),  # meV
    "ratio": (truth["ratio"] * (1 + rng.normal(0, 0.1)), 0.1 * truth["ratio"]),
    "chi_c": (truth["chi_c"] * (1 + rng.normal(0, 0.02, 30)), 0.02 * truth["chi_c"]),
    "chi_a": (truth["chi_a"] * (1 + rng.normal(0, 0.02, 30)), 0.02 * truth["chi_a"]),
    "M_c": (truth["M_c"] + rng.normal(0, 0.01, 20), np.full(20, 0.01)),  # muB
}


def residuals(pars):
    """Weighted residuals, S_p + (all data points,)."""
    m = model(pars)
    return np.concatenate(
        [(m[key] - value) / sigma for key, (value, sigma) in data.items()], axis=-1
    )


# =============================================================================
# 1. Vectorised scan of the parameter space for a start value.
# =============================================================================
# B20 (21, 1, 1), B40 (1, 21, 1), B44 (1, 1, 21) broadcast to a (21, 21, 21)
# grid: the same model() as for the fit, 9261 parameter sets at once.

# MS: I do not want to worry about broadcasting here
grid = Pars.from_dict(
    dict(
        B20=np.linspace(-1, 1, 21),
        B40=np.linspace(-0.05, 0.05, 21),
        B44=np.linspace(-0.3, 0.3, 21),
    )
)
t0 = time.perf_counter()
chi2 = np.sum(residuals(grid) ** 2, axis=-1)  # (21, 21, 21)
chi2[~np.isfinite(chi2)] = np.inf
best = np.unravel_index(np.argmin(chi2), chi2.shape)
start = {name: float(v) for name, v in zip(grid.names, grid.values[best])}
print(f"scan: {chi2.size} parameter sets in {time.perf_counter() - t0:.2f} s, start {start}")


# =============================================================================
# 2. Local fit with lmfit, from the scan's best point.
# =============================================================================
# lmfit parameters carry the CEFpars names, so the conversion is valuesdict().
# Constraints are lmfit expressions, e.g. for a cubic site
# params["B44"].expr = "5 * B40".

params = lmfit.Parameters()
for name, value in start.items():
    params.add(name, value=value)

t0 = time.perf_counter()
result = lmfit.minimize(lambda prm: residuals(Pars.from_dict(prm.valuesdict())), params)
print(
    f"fit: {result.nfev} evaluations in {time.perf_counter() - t0:.2f} s, reduced chi2 {result.redchi:.2f}"
)
for name in TRUE:
    prm = result.params[name]
    print(f"  {name}: {prm.value:+.5f} +- {prm.stderr:.5f} meV (true {TRUE[name]:+.5f})")

# WISH: an analytic Jacobian for lmfit (Dfun): dE_n/dB_k = <n|O_k|n> comes from the
# same eigensystem (Hellmann-Feynman); the derivatives of chi and M need more work.
# WISH: fit whole spectra with a resolution function instead of peak lists.
# WISH: several sites or ions in one fit: one Pars per site, residuals concatenated.
