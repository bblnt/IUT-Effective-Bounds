#!/usr/bin/env python3
"""
IUT Effective Bounds Evaluator v1.0
====================================
Evaluates the published effective abc/Szpiro bounds derived from
Inter-Universal Teichmüller Theory against known abc-triples.

This engine takes the PUBLISHED formula from:
  Mochizuki-Fesenko-Hoshi-Minamide-Porowski (2022),
  "Explicit estimates in inter-universal Teichmüller theory",
  Kodai Math. J. 45(2), 175-236. DOI: 10.2996/kmj45201

Their Theorem A states (for Q, i.e. mono-complex with d=1):
  For coprime a + b + c = 0, and 0 < ε ≤ 1:
    |abc| ≤ 2^4 · max{ exp(1.7 × 10^30 · ε^(-166/81)), rad(abc)^(3+3ε) }

We compare this to the abc conjecture's requirement:
    c ≤ K_ε · rad(abc)^(1+ε)

TWO independent problems with the MFHMP bound:
  1. EXPONENT GAP: 3+3ε vs 1+ε → ratio approaches 3 as ε→0
  2. CONSTANT WALL: exp(1.7×10^30 · ε^(-166/81)) is cosmologically large

This engine computes both gaps for known high-quality abc triples.

NO disputed mathematics is used. All formulas come from published,
peer-reviewed results in Kodai Mathematical Journal.

Author: B. Kolláth / Algorion Technologies
Date: March 2026
"""

from sympy import factorint, log as sym_log, isprime, Integer, Rational
from mpmath import mp, mpf, log as mp_log, exp as mp_exp, power as mp_power, fabs, floor
import csv
import json

mp.dps = 50  # High precision

# ============================================================
# KNOWN HIGH-QUALITY ABC TRIPLES
# Source: Bart de Smit's ABC@Home project, Abderrahmane Nitaj's tables
# Quality q(a,b,c) = log(c) / log(rad(abc))
# ============================================================

ABC_TRIPLES = [
    {
        "name": "Reyssat (1987)",
        "a": 2,
        "b": 3**10 * 109,  # = 6436341
        "c": 23**5,         # = 6436343
        "quality_published": 1.6299,
        "reference": "Reyssat, 1987"
    },
    {
        "name": "de Weger (1998)",
        "a": 11**2,          # = 121
        "b": 3**2 * 5**6 * 7**3,  # = 5375^... let me compute
        "c": 2**21 * 23,     # = 48234496... 
        "quality_published": 1.6259,
        "reference": "de Weger, 1998"
    },
    {
        "name": "Nitaj (1993) #1",
        "a": 1,
        "b": 2 * 3**7,      # = 4374
        "c": 5**4 * 7,       # = 4375
        "quality_published": 1.5679,
        "reference": "Nitaj, 1993"
    },
    {
        "name": "de Weger (1997)",
        "a": 7**3,           # = 343
        "b": 3**10,          # = 59049
        "c": 2**11 * 29,     # = 59392
        "quality_published": 1.4553,
        "reference": "de Weger, 1997"
    },
    {
        "name": "Classic: 1+8=9",
        "a": 1,
        "b": 8,              # = 2^3
        "c": 9,              # = 3^2
        "quality_published": 1.2263,
        "reference": "classical"
    },
    {
        "name": "Classic: 5+27=32",
        "a": 5,
        "b": 27,             # = 3^3
        "c": 32,             # = 2^5
        "quality_published": 1.0932,
        "reference": "classical"
    },
    {
        "name": "Classic: 2+6436341=6436343",
        "a": 2,
        "b": 3**10 * 109,
        "c": 23**5,
        "quality_published": 1.6299,
        "reference": "Reyssat 1987 (duplicate verification)"
    },
]


def compute_radical(a, b, c):
    """Compute rad(abc) = product of distinct prime factors of abc."""
    all_factors = set()
    for n in [abs(a), abs(b), abs(c)]:
        if n > 0:
            all_factors.update(factorint(n).keys())
    result = 1
    for p in all_factors:
        result *= p
    return result


def compute_quality(a, b, c):
    """Compute q(a,b,c) = log(c) / log(rad(abc))."""
    rad = compute_radical(a, b, c)
    if rad <= 1:
        return float('inf')
    return float(mp_log(mpf(c)) / mp_log(mpf(rad)))


def mfhmp_bound(a, b, c, eps):
    """
    Compute the MFHMP effective bound from Kodai 2022, Theorem A.

    For coprime a+b+c=0 over Q (d=1, mono-complex):
      |abc| ≤ 2^4 · max{ exp(C · ε^(-166/81)), rad(abc)^(3+3ε) }
    where C = 1.7 × 10^30.

    Returns dict with both terms and which dominates.
    """
    rad = compute_radical(a, b, c)
    abc_product = abs(a * b * c)

    # Term 1: the exponential constant
    C = mpf('1.7e30')
    exponent1 = C * mp_power(mpf(eps), mpf('-166') / mpf('81'))
    # We can't compute exp(10^30) directly — it has ~10^30 digits
    # Instead we work in log space
    log_term1 = float(exponent1)  # this is log of term1 (since term1 = exp(exponent1))

    # Term 2: rad(abc)^(3+3ε)
    log_term2 = float((3 + 3*eps) * mp_log(mpf(rad)))

    # The bound is 2^4 * max(term1, term2)
    log_16 = float(mp_log(mpf(16)))
    log_bound = log_16 + max(log_term1, log_term2)

    # The actual |abc|
    log_abc = float(mp_log(mpf(abc_product)))

    # For abc conjecture comparison: c ≤ K_ε · rad^(1+ε)
    log_abc_bound = float((1 + eps) * mp_log(mpf(rad)))

    return {
        'eps': eps,
        'log_term1_exp': log_term1,      # log of the exponential term
        'log_term2_rad': log_term2,       # log of the radical term
        'dominant': 'exponential' if log_term1 > log_term2 else 'radical',
        'log_mfhmp_bound': log_bound,     # log of the full MFHMP bound
        'log_abc_actual': log_abc,         # log of actual |abc|
        'mfhmp_exponent': 3 + 3*eps,      # IUT's exponent
        'abc_exponent': 1 + eps,           # abc conjecture's exponent
        'exponent_ratio': (3 + 3*eps) / (1 + eps),
        'log_abc_conj_bound': log_abc_bound,
        'digits_in_constant': int(log_term1 / float(mp_log(10))) if log_term1 > 0 else 0,
    }


def analyze_triple(triple, verbose=True):
    """Full analysis of one abc triple against MFHMP bounds."""
    a, b, c = triple['a'], triple['b'], triple['c']
    name = triple['name']

    # Verify a + b = c
    assert a + b == c, f"FAIL: {a} + {b} ≠ {c}"

    rad = compute_radical(a, b, c)
    quality = compute_quality(a, b, c)
    log_c = float(mp_log(mpf(c)))
    log_rad = float(mp_log(mpf(rad)))

    if verbose:
        print(f"\n{'='*65}")
        print(f"  {name}")
        print(f"{'='*65}")
        print(f"  a = {a}")
        print(f"  b = {b}")
        print(f"  c = {c}")
        print(f"  a + b = c: {a + b == c}")
        print(f"  rad(abc) = {rad}")
        print(f"  quality q = log(c)/log(rad) = {quality:.4f}"
              f"  (published: {triple['quality_published']})")
        print(f"  log(c) = {log_c:.4f}")
        print(f"  log(rad) = {log_rad:.4f}")

    # Analyze at multiple epsilon values
    eps_values = [0.001, 0.01, 0.1, 0.5, 1.0]
    results = []

    if verbose:
        print(f"\n  {'ε':>8} | {'IUT exp':>9} | {'abc exp':>8} | {'Ratio':>6} | "
              f"{'Dominant':>12} | {'Digits in K':>11}")
        print(f"  {'-'*70}")

    for eps in eps_values:
        r = mfhmp_bound(a, b, c, eps)
        results.append(r)
        if verbose:
            print(f"  {eps:>8.3f} | {r['mfhmp_exponent']:>9.3f} | "
                  f"{r['abc_exponent']:>8.3f} | {r['exponent_ratio']:>5.2f}× | "
                  f"{r['dominant']:>12} | {r['digits_in_constant']:>11,}")

    return {
        'name': name,
        'a': a, 'b': b, 'c': c,
        'radical': rad,
        'quality': quality,
        'log_c': log_c,
        'log_rad': log_rad,
        'bounds': results,
    }


# ============================================================
# MAIN EXECUTION
# ============================================================

if __name__ == '__main__':
    print("=" * 65)
    print("  IUT EFFECTIVE BOUNDS EVALUATOR v1.0")
    print("  B. Kolláth / Algorion Technologies / March 2026")
    print("  Source: MFHMP (2022), Kodai Math. J. 45(2), 175-236")
    print("=" * 65)

    all_results = []

    for triple in ABC_TRIPLES:
        result = analyze_triple(triple, verbose=True)
        all_results.append(result)

    # ========================================
    # SUMMARY: THE TWO STRUCTURAL PROBLEMS
    # ========================================
    print("\n" + "=" * 65)
    print("  FINDING 1: THE EXPONENT GAP")
    print("=" * 65)
    print(f"\n  The abc conjecture requires: c ≤ K_ε · rad(abc)^(1+ε)")
    print(f"  The MFHMP bound gives:      |abc| ≤ ... · rad(abc)^(3+3ε)")
    print(f"\n  Exponent comparison:")
    print(f"  {'ε':>8} | {'abc needs':>10} | {'IUT gives':>10} | {'Ratio':>7}")
    print(f"  {'-'*45}")
    for eps in [0.001, 0.01, 0.1, 0.5, 1.0]:
        abc_exp = 1 + eps
        iut_exp = 3 + 3*eps
        print(f"  {eps:>8.3f} | {abc_exp:>10.3f} | {iut_exp:>10.3f} | {iut_exp/abc_exp:>6.2f}×")
    print(f"\n  As ε → 0, the ratio → 3.000×")
    print(f"  This gap is STRUCTURAL: it comes from the passage through")
    print(f"  Szpiro's conjecture with IUT's specific constants.")
    print(f"  It cannot be reduced without changing the proof strategy.")

    print("\n" + "=" * 65)
    print("  FINDING 2: THE CONSTANT WALL")
    print("=" * 65)
    print(f"\n  The MFHMP bound contains: exp(1.7 × 10^30 × ε^(-166/81))")
    print(f"\n  At ε = 0.01:")
    C = mpf('1.7e30')
    exp_at_001 = C * mp_power(mpf('0.01'), mpf('-166') / mpf('81'))
    log10_val = float(exp_at_001 / mp_log(10))
    print(f"    The constant has approximately {log10_val:.2e} decimal digits")
    print(f"    For comparison:")
    print(f"      Atoms in the observable universe:  ~10^80")
    print(f"      This constant:                     ~10^(10^34)")
    print(f"      Ratio:                             10^(10^34 - 80)")
    print(f"\n  Even with the CORRECT exponent 1+ε, the constant renders")
    print(f"  the bound vacuous for any concretely specified abc triple.")

    print("\n" + "=" * 65)
    print("  FINDING 3: COMBINED EFFECT ON KNOWN TRIPLES")
    print("=" * 65)
    print(f"\n  For the highest-quality known abc triple (Reyssat, q=1.630):")
    print(f"    The MFHMP bound is satisfied with room to spare by a factor")
    print(f"    of approximately 10^(10^34). The bound provides zero")
    print(f"    Diophantine constraint.")
    print(f"\n  For ANY known abc triple: the exponential constant dominates")
    print(f"  the radical term at all practical values of ε. The bound")
    print(f"  reduces to |abc| ≤ (cosmologically large number), which is")
    print(f"  trivially true and informationally vacuous.")

    print("\n" + "=" * 65)
    print("  CONCLUSION")
    print("=" * 65)
    print(f"\n  The MFHMP effective bounds, derived from IUT theory, suffer")
    print(f"  from two independent structural limitations:")
    print(f"    1. An exponent 3× larger than what the abc conjecture requires")
    print(f"    2. A constant with ~10^34 digits")
    print(f"  Both are inherent to the framework as published.")
    print(f"  These limitations are independent of the validity of")
    print(f"  Corollary 3.12 — they apply even granting the full theory.")

    # ========================================
    # SAVE RESULTS
    # ========================================
    csv_path = 'mfhmp_evaluation.csv'
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['name', 'a', 'b', 'c', 'radical', 'quality',
                         'log_c', 'log_rad',
                         'eps', 'iut_exponent', 'abc_exponent',
                         'exponent_ratio', 'dominant_term',
                         'digits_in_constant'])
        for result in all_results:
            for bound in result['bounds']:
                writer.writerow([
                    result['name'], result['a'], result['b'], result['c'],
                    result['radical'], f"{result['quality']:.4f}",
                    f"{result['log_c']:.4f}", f"{result['log_rad']:.4f}",
                    bound['eps'], f"{bound['mfhmp_exponent']:.3f}",
                    f"{bound['abc_exponent']:.3f}",
                    f"{bound['exponent_ratio']:.3f}",
                    bound['dominant'],
                    bound['digits_in_constant'],
                ])
    print(f"\n  Results saved to {csv_path}")
