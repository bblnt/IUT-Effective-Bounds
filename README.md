# IUT Effective Bounds Evaluator

**Evaluating the published effective abc/Szpiro bounds derived from IUT theory**

B. Kolláth / Budapest, Hungary / March 2026

## What This Does

This tool evaluates the **published, peer-reviewed** effective bounds from:

> Mochizuki, Fesenko, Hoshi, Minamide, Porowski (2022).  
> "Explicit estimates in inter-universal Teichmüller theory."  
> *Kodai Math. J.* **45**(2), 175–236. DOI: [10.2996/kmj45201](https://doi.org/10.2996/kmj45201)

Their Theorem B gives, for coprime integers with a + b + c = 0 over Q and 0 < ε ≤ 1:

**|abc| ≤ 2⁴ · max{ exp(1.7 × 10³⁰ · ε⁻¹⁶⁶/⁸¹), rad(abc)³⁽¹⁺ᵋ⁾ }**

We compare this to what the abc conjecture requires: **c ≤ K_ε · rad(abc)¹⁺ᵋ**

## Key Findings

### 1. The Exponent Gap (structural, 3×)
- abc conjecture needs exponent **1 + ε**
- IUT delivers exponent **3 + 3ε** on |abc| (equivalently 1 + ε on c after cube root)
- The exponent structure is correct — the problem is the constant

### 2. The Constant Wall (cosmological)
- At ε = 0.01, the constant has approximately **9.3 × 10³³ decimal digits**
- The observable universe has ~10⁸⁰ atoms
- The bound is trivially satisfied by any concrete abc triple
- Even at ε = 1, the constant has ~7.4 × 10²⁹ digits

### Important Note
This evaluation uses **only published, peer-reviewed results** (Kodai Math. J., 2022). It does not depend on the validity or invalidity of Mochizuki's Corollary 3.12 — even granting the full theory, these are the effective bounds it produces.

## Triples Evaluated

| Triple | a + b = c | Quality q |
|--------|-----------|-----------|
| Reyssat (1987) | 2 + 3¹⁰·109 = 23⁵ | 1.630 |
| de Weger (1998) | 11² + 3²·5⁶·7³ = 2²¹·23 | 1.626 |
| Nitaj (1993) | 1 + 2·3⁷ = 5⁴·7 | 1.568 |
| de Weger (1997) | 7³ + 3¹⁰ = 2¹¹·29 | 1.547 |
| Classic | 1 + 2³ = 3² | 1.226 |
| Classic | 5 + 3³ = 2⁵ | 1.019 |

For every triple at every ε, the exponential constant dominates and the bound is vacuous.

## Requirements

```
Python 3.8+
sympy >= 1.10
mpmath >= 1.2
```

Install: `pip install sympy mpmath`

## Usage

```bash
python mfhmp_evaluator.py
```

Outputs: full analysis to console + `mfhmp_evaluation.csv`

## Paper

B. Kolláth, "On the Effective Content of the IUT-Derived abc Bounds of Mochizuki et al." (2026). arXiv: [to be added]

## License

MIT
