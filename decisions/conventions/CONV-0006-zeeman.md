---
id: CONV-0006
title: Zeeman term and magnetic moment
status: draft
owner: qespace_cef
supersedes:
superseded_by:
---

## Statement

The magnetic moment of an ion and its Zeeman term in a magnetic field are

    μ = −gJ μB J,        ℋ_Z = −μ·B = +gJ μB J·B,

with J the angular momentum of the ground J multiplet (CONV-0005, ħ = 1),
gJ its Landé factor and B the magnetic induction acting on the ion, in T
(CONV-0003), in the local frame of CONV-0004. The moment is antiparallel to
J: in a field the ground state has ⟨J⟩ against B and the moment along it.
No demagnetisation or Lorentz correction is applied; exchange and molecular
fields enter as their own terms.

The constants are positive magnitudes: e, μB = eħ/2m_e, gJ, and the
neutron's γ = |μ_n/μ_N| (CONV-0007). Every sign is written in the formula.

## Rationale

The minus sign is the electron's charge, −e. A charge q with angular
momentum L carries the moment (q/2m) L, so the electron's orbital moment is
−μB L, its spin moment −g_s μB S, and inside the ground multiplet
`L + g_s S ≈ gJ J` (projection theorem), which gives μ = −gJ μB J.

Two conventions are in use, and the difference comes down to where the
signs of the charge and of g are put:

1. **Electron sign** (this CONV): e > 0, g > 0, and the charge's minus
   written in μ. Wikipedia "Zeeman effect", PyCrystalField (paper II.C and
   II.E), SpinW and Sunny. CODATA's bookkeeping is the same physics with the
   sign inside g: g_e = −2.0023… and μ_e = g_e μB S.
2. **Moment parallel to J**: μ = +gJ μB J, ℋ_Z = −gJ μB J·H, with H given
   as μ0H in T. Jensen and Mackintosh (2.1.1c), McPhase (manual, Zeeman
   energy), Mantid. This absorbs the charge's sign into J.

They are related by time reversal Θ, equivalently by reversing the field:
the crystal field is even under Θ, so the energies, the magnetisation along
the field, the susceptibility and the unpolarised cross sections are the
same. They differ in what is odd under Θ: the sign of ⟨J⟩, the
eigenvectors in a field (Θ-partners), and the antisymmetric part
S^{αβ} − S^{βα} that polarised neutrons see. Mapping J → −J is no
substitute: it breaks [J_x, J_y] = iJ_z.

The electron sign is the physical one: J is the electron's angular
momentum, and a field B is the induction the ion feels.

## Consequences

- One sign for the Zeeman term, the moment, the magnetisation and the mean
  field.
- Results odd under time reversal are compared with codes of the other
  convention through Θ.
- Tests: without a crystal field, in a field B along z, the ground state
  has ⟨J_z⟩ = −J and the moment +gJ μB J along B; in a field much larger
  than the crystal field the moment points along B and approaches
  gJ μB J.
