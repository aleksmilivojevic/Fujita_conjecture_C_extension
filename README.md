# Fano manifolds with the integral cohomology ring of projective space

This repository contains the article, its TeX source, and the Python program used for its computational arguments.

A Fano manifold is a smooth projective complex manifold whose anticanonical line bundle is ample. Theorem 1.2 states that a Fano manifold of dimension `n <= 10` with the integral cohomology ring of projective space is biholomorphic to projective space.

The Fano index is the largest positive integer `r` such that `-K_X = rA` for a holomorphic line bundle `A`. Under the cohomology-ring hypothesis, `A` is the ample generator `L`. Theorem 1.3 proves the same conclusion in arbitrary dimension when `r >= n - 3`, or when `n` is even and `r >= n - 7`.

The proof starts with the Kobayashi–Ochiai characterization at index `n + 1`. Parity and Chern-number divisibility give a finite list of other possible indices. Riemann–Roch and bounds on section dimensions then give polynomial equations with bounded integer parameters. The program derives these equations and excludes the permitted parameters by integer-square or congruence tests. The cases `(n,r) = (15,12), (27,24)` also use the geometric argument about genus-two curve sections proved in the article.

## Files

| File | Contents |
| --- | --- |
| `fujita_further_extensions.pdf` | The article, including the full arguments in dimensions seven and eight and the polynomial and sieve appendices. |
| `fujita_further_extensions.tex` | Standalone TeX source. The bibliography is included; no Python source is embedded. |
| `verify_fujita_further_extensions.py` | The verification program. |
| `requirements.txt` | The SymPy version used for the recorded run. |
| `verification.log` | Output of the recorded full run. |
| `audit.json` | Derived polynomials, eliminants, residue values, and sieve counts from that run. |

## Run the verification

Python 3.10 or later is required. The recorded run used Python 3.12.13 and SymPy 1.14.0.

```sh
python -m pip install -r requirements.txt
python verify_fujita_further_extensions.py --export-json audit.json
```

The program derives its equations and executes every finite check. It ends with:

```text
ALL PYTHON CHECKS PASSED.
```

Run Python without `-O` or `-OO`, since the verification uses assertions. There is no C++ dependency. The script needs only SymPy and Python's standard library; `audit.json` is an output, not an input.

Section “Computational verification” of the article maps the calculations to the corresponding functions and explains how the large candidate sets are represented. All polynomial and residue calculations use exact integer or rational arithmetic. The geometric lemmas are proved in the article; they are not encoded in a formal proof assistant.

In `audit.json`, the high-index eliminants store coefficients in descending powers. General polynomial records store variables and pairs of exponent vectors and integer coefficients. The optional `coefficients_latex` lists run in increasing powers of the named root variable.

## Compile the article

Use pdfLaTeX with a standard TeX Live or MiKTeX installation:

```sh
pdflatex -interaction=nonstopmode -halt-on-error fujita_further_extensions.tex
pdflatex -interaction=nonstopmode -halt-on-error fujita_further_extensions.tex
```

The second pass resolves the equation and bibliography references. Compilation does not execute or include the verification program.
