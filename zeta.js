/* Riemann–Siegel Z(t) with the first two correction terms, the Riemann–Siegel
   theta function, zero isolation by sign changes + Brent-style bisection, and
   the Riemann–von Mangoldt count. Shared by the page and node tests. */
(function (root) {
  'use strict';
  const TAU = 2 * Math.PI;
  // θ(t) = arg Γ(1/4 + it/2) − (t/2) ln π, asymptotic expansion (excellent for t > 5)
  function theta(t) {
    const t2 = t * t;
    return t / 2 * Math.log(t / TAU) - t / 2 - Math.PI / 8 + 1 / (48 * t) + 7 / (5760 * t * t2) + 31 / (80640 * t2 * t2 * t);
  }
  // Ψ(p) and its derivatives for the correction terms (Gabcke's C0, C1)
  function psi(p) { return Math.cos(TAU * (p * p - p - 1 / 16)) / Math.cos(TAU * p); }
  function d3psi(p) { // third derivative of psi, by a 5-point stencil (stable for this smooth function)
    const h = 1e-3;
    return (psi(p + 2 * h) - 2 * psi(p + h) + 2 * psi(p - h) - psi(p - 2 * h)) / (2 * h * h * h);
  }
  function safePsi(p, f) {
    // psi has removable singularities at p = 1/4, 3/4; nudge off them
    if (Math.abs(Math.cos(TAU * p)) < 1e-6) p += 1e-6;
    return f(p);
  }
  /* Exact-ish ζ(1/2+it) for small t: Borwein's accelerated alternating
     series for η, then ζ = η / (1 − 2^{1−s}). Used below t = 120, where Riemann–Siegel is least accurate. */
  const BORWEIN_N = 260;
  const borweinD = (() => { // d_k/d_n, built in log space: the raw terms reach 5.83^n
    const n = BORWEIN_N, lg = new Float64Array(n + 1);
    for (let i = 1; i <= n; i++) lg[i] = lg[i - 1] + Math.log((n + i - 1) * (n - i + 1) * 4 / ((2 * i - 1) * (2 * i)));
    const mx = Math.max(...lg), d = new Float64Array(n + 1);
    let sum = 0;
    for (let i = 0; i <= n; i++) { sum += Math.exp(lg[i] - mx); d[i] = sum; }
    for (let i = 0; i <= n; i++) d[i] /= sum;
    return d;
  })();
  function zetaBorwein(t) {
    const n = BORWEIN_N; let re = 0, im = 0;
    for (let k = 0; k < n; k++) {
      // (−1)^k (d_k − d_n) / (k+1)^s,  s = 1/2 + it
      const c = (k % 2 ? -1 : 1) * (borweinD[k] - 1), l = Math.log(k + 1), m = c / Math.sqrt(k + 1);
      re += m * Math.cos(t * l); im -= m * Math.sin(t * l);
    }
    re = -re; im = -im;
    // divide by 1 − 2^{1−s} = 1 − √2·e^{−it ln2}
    const dr = 1 - Math.SQRT2 * Math.cos(t * Math.LN2), di = Math.SQRT2 * Math.sin(t * Math.LN2), dd = dr * dr + di * di;
    return [(re * dr + im * di) / dd, (im * dr - re * di) / dd];
  }
  function Z(t) {
    if (t < 120) { const [zr, zi] = zetaBorwein(t), th = theta(t); return zr * Math.cos(th) - zi * Math.sin(th); }
    return Zrs(t);
  }
  function Zrs(t) {
    const a = Math.sqrt(t / TAU), N = Math.floor(a), p = a - N, th = theta(t);
    let s = 0;
    for (let n = 1; n <= N; n++) s += Math.cos(th - t * Math.log(n)) / Math.sqrt(n);
    s *= 2;
    const c0 = safePsi(p, psi);
    const c1 = -safePsi(p, d3psi) / (96 * Math.PI * Math.PI);
    const sign = (N - 1) % 2 === 0 ? 1 : -1;
    const w = Math.pow(t / TAU, -0.25);
    return s + sign * w * (c0 + c1 * w * w);
  }
  // ζ(1/2 + it) = e^{-iθ(t)} Z(t)
  function zetaCritical(t) { const z = Z(t), th = theta(t); return [z * Math.cos(th), -z * Math.sin(th)]; }

  function refine(a, b, fa) {
    for (let i = 0; i < 60; i++) {
      const m = (a + b) / 2, fm = Z(m);
      if ((fm < 0) === (fa < 0)) { a = m; fa = fm; } else b = m;
      if (b - a < 1e-11) break;
    }
    return (a + b) / 2;
  }
  // zeros in [t0, t1]: sample finer than the mean gap, which shrinks like 2π/ln(t/2π)
  function zeros(t0, t1) {
    const out = [];
    let t = Math.max(t0, 9), z = Z(t);
    while (t < t1) {
      const gap = TAU / Math.log(Math.max(t, 20) / TAU);
      const h = gap / 24;
      const t2 = Math.min(t + h, t1), z2 = Z(t2);
      if ((z < 0) !== (z2 < 0)) out.push(refine(t, t2, z));
      t = t2; z = z2;
    }
    return out;
  }
  // Riemann–von Mangoldt main term: N(T) ≈ θ(T)/π + 1
  const count = (T) => theta(T) / Math.PI + 1;
  const api = { Z, theta, zetaCritical, zeros, count };
  if (typeof module === 'object' && module.exports) module.exports = api; else root.Zeta = api;
})(typeof globalThis !== 'undefined' ? globalThis : this);
