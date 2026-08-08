# Nineteen vectors and their phase-retrieval injectivity in ℂ⁶

Companion repository for the note [`nineteen_c6.pdf`](nineteen_c6.pdf).

**Claim.** The nineteen explicit vectors listed in Table 1 of the note —
every entry a Gaussian integer — do phase retrieval in ℂ⁶: the intensities
|⟨a_j, x⟩|², j = 1, …, 19, determine every x ∈ ℂ⁶ up to a unimodular
scalar. Hence m_ℂ(6) ≤ 19 = 4d−5, one below the generic count 4d−4;
combined with the Wang–Xu lower bound m_ℂ(6) ≥ 18, the minimal
measurement number in dimension six is 18 or 19.

## Note

This frame and its certificate were extracted from a larger project on
the 4d−5 measurement problem. The full development — additional
certified frames in d = 4, 7, 8, 10, a Lean 4 formalization, and the
classification preprint — lives at
**[github.com/rkka02/phase-retrieval-4d-5](https://github.com/rkka02/phase-retrieval-4d-5)**
(archived at [doi:10.5281/zenodo.21764847](https://doi.org/10.5281/zenodo.21764847)).

## The certificate

After the elementary reduction proved in the note (§2), phase retrieval
for this frame is equivalent to a system of five quadratic equations in
five unknowns having **no real solutions**. [`certificate/`](certificate/)
contains everything needed to check this claim independently:

| File | Contents |
|---|---|
| `frame_gaussian_integer.json` | the frame, exact integers |
| `standard_system_d6.ms` | the polynomial system, msolve input format |
| `certificate.json` | the degree-14 squarefree integer eliminant, its Sturm sequence, sign counts, SHA-256 |
| `verify_exact.py` | standalone exact verifier (Python 3.10+, SymPy): rebuilds the system from scratch, recomputes the Gröbner basis, eliminant, and Sturm chain, and checks rank_ℝ L_A = 19 — all in exact rational arithmetic |
| `rerun_2026-08-05.txt` | log of the latest rerun |

Rerun it yourself (a few minutes on a laptop):

```bash
python3 certificate/verify_exact.py
```

Expected last lines: `eliminant degree=14`, `squarefree=True`,
`real roots=0`, `rank_R L_A=19`.

## AI assistance

The author received assistance from an AI assistant in the search for
candidate frames and in the preparation of the note and this repository.

## License

MIT ([`LICENSE`](LICENSE)); the note text (`nineteen_c6.tex`/`.pdf`) is CC BY 4.0.
Cite via [`CITATION.cff`](CITATION.cff). This repository is archived at
[doi:10.5281/zenodo.21798931](https://doi.org/10.5281/zenodo.21798931).
