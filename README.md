# Shunya · शून्य

**Hunting the zeros of the Riemann zeta function, in your browser.**
Live: **https://shunya.avlokai.com**

The Riemann hypothesis says every non-trivial zero of ζ(s) has real part exactly ½. Shunya walks that
critical line: it evaluates ζ(½ + it), finds every zero in view to about ten decimal places, draws the value
of ζ spiralling through the origin each time it vanishes, and checks the number of zeros it found against
the count Riemann's formula predicts. It works anywhere from t = 10 up to t = 10⁹.

---

## What it does

- **Hardy's Z(t) graph.** A real-valued function that is zero exactly where ζ(½ + it) is. Drag to pan,
  scroll to zoom, hover to read Z(t). Zeros are marked.
- **The complex plane.** The path traced by ζ(½ + it) over the same window, a family of loops that pass
  through 0 once per zero.
- **Zero list** to 10 decimals for the visible window.
- **Riemann–von Mangoldt check.** Zeros found in the window versus θ(T)/π + 1 at its two ends.
- **Count every zero up to t = 10,000** from scratch with one button: it finds all 10,142 and reports the
  difference from the formula.
- Presets: the start, **Lehmer's pair** near t = 7005 (two zeros 0.038 apart, a famous near-miss for the
  hypothesis), t = 10⁵ and t = 10⁹. Links are shareable (`#t=…&w=…`).

## The mathematics

On the critical line ζ can be written as

```
ζ(½ + it) = e^{−iθ(t)} · Z(t)
```

where θ(t) is the Riemann–Siegel theta function and Z(t) is real. So the zeros of ζ on the line are exactly
the sign changes of Z: easy to detect and easy to pin down by bisection.

**θ(t)** is evaluated from its asymptotic expansion:

```
θ(t) = (t/2)·ln(t/2π) − t/2 − π/8 + 1/(48t) + 7/(5760t³) + 31/(80640t⁵) + …
```

**Z(t) for t ≥ 120** uses the Riemann–Siegel formula:

```
Z(t) = 2 Σ_{n=1}^{N} n^{−½} cos(θ(t) − t ln n)  +  R(t),     N = ⌊√(t/2π)⌋
```

The remainder R uses the first two correction terms, C₀ and C₁, built from
Ψ(p) = cos(2π(p² − p − 1/16)) / cos(2πp), with p the fractional part of √(t/2π). The main sum has only
√(t/2π) terms, which is why t = 10⁹ needs about 12,600 cosines instead of a billion.

**Z(t) for t < 120** is computed from ζ directly, because the asymptotic remainder is least accurate there.
It uses **Borwein's accelerated alternating series** for the Dirichlet eta function and
ζ(s) = η(s) / (1 − 2^{1−s}), with 260 terms. The Borwein weights grow like 5.83ⁿ, which overflows doubles,
so they are accumulated in log space and normalised.

**Finding zeros.** Z is sampled at intervals of 1/24 of the mean zero spacing, 2π / ln(t/2π), so close pairs
like Lehmer's are not stepped over. Each sign change is refined by bisection to a bracket of 10⁻¹¹.

**Counting.** The Riemann–von Mangoldt formula gives N(T) = θ(T)/π + 1 + S(T), where S(T) is small. If the
sign changes found on the line account for all the zeros the formula says exist, none can be off the line in
that range. This is the idea behind Turing's method, in its simplest form, and behind the large-scale
machine checks of the hypothesis.

## How it was built

`zeta.js` holds the numerics and runs in both browser and Node; `index.html` is the interface.

| Piece | Notes |
| --- | --- |
| `theta(t)` | Asymptotic series above. |
| `zetaBorwein(t)` | Log-space Borwein weights, complex accumulation, division by 1 − 2^{1−s}. |
| `Zrs(t)` | Riemann–Siegel main sum plus C₀ and C₁ (C₁ uses a finite-difference third derivative of Ψ), with a nudge off Ψ's removable singularities. |
| `Z(t)` | Borwein below 120, Riemann–Siegel above. |
| `zeros(t0, t1)` | Adaptive sampling at a fraction of the local mean gap, then bisection. |
| `count(T)` | θ(T)/π + 1. |

The page draws both plots on `<canvas>` at device-pixel resolution, fills the positive and negative lobes of
Z in two colours, and reads its palette from CSS custom properties so it follows light and dark mode. The
"count to 10,000" job runs in 250-unit slices on `setTimeout` so the page stays responsive.

### What went wrong on the way (and why it matters)

- The first version used Riemann–Siegel everywhere. The first zero came out at 14.1372 instead of 14.1347:
  the asymptotic formula is weakest at small t. That is why there is a Borwein path below 120.
- The first Borwein attempt used 560 terms in plain doubles; the weights overflowed to Infinity and silently
  lost 201 zeros. Moving the weights to log space fixed it, and the exhaustive count test now guards it.

## Verification

```bash
node test.cjs
```

- The first ten zeros match published values to better than 10⁻⁸ (in practice about 10⁻¹⁰).
- N(1000) = 649 and N(10000) = 10,142, both exact.
- Lehmer's pair near t = 7005.06 / 7005.10 is resolved as two separate zeros.

## Limits

- Double precision caps useful heights at around t = 10⁹–10¹⁰, where θ(t) itself loses digits.
- The count check uses only the main term of N(T), so it shows S(T) rather than proving bounds on it.
  This is a demonstration, not a certified verification.

## Run locally

```bash
python -m http.server 8000   # then http://localhost:8000
```

## Files

```
index.html   graphs, interaction, the count job
zeta.js      θ, Borwein ζ, Riemann–Siegel Z, zero isolation, N(T)
test.cjs     published zeros, exact counts, Lehmer's pair
```

## Name

*Śūnya* (शून्य) is Sanskrit for zero, emptiness — the word Indian mathematics gave the world for nothing.

---

Built by Sushanth Kasturi as one of a set of deliberately useless things. More at
https://sushanth.avlokai.com/projects
