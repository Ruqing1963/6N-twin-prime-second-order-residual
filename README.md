# 6N Twin-Prime Second-Order Residual (Part IV)

A **second-order** conditional singular series for the twin-prime gap
distribution on the 6N ± 1 skeleton: the **spatial-compression penalty** that
resolves the long-gap (210) residual left by Part III, and a residual
**two-centre (shield-gain)** term on the short gap (42).

**Background.** Part II found the twin-gap distribution depends on ω₍>3₎(N), the
distinct-prime-factor count of the centre (42 rises with ω, 210 falls). Part III
gave the first-order conditional singular series S_{2,ω}(d) ≈ C2(d)·L_ω(d) (with
the congruence-lockdown allowance L_ω), which reproduced 30/42/60 and the
*direction* of the 210 collapse but left a stable high-ω residual on 210
(observed 0.41 where first order predicts 0.62 at ω=6).

**This part — the second-order term.** The first-order series treats a prime
q ∤ N by its unconditional Hardy–Littlewood factor. But conditioning on q ∤ N
removes the class N≡0 and crowds N into the remaining q−1 classes, where the
same k_q fatal residues are proportionally denser. This **spatial-compression
penalty**

```
    pen_q(d) = [ (q-1-k_q)/(q-1) ] / [ (q-k_q)/q ]  < 1      (for k_q=4: (q²-5q)/(q²-5q+4))
```

lowers factor-rich centres' survival on gaps locked by small primes. The
second-order weight of a gap d at a real centre N with factor set Q is

```
    w2(d,N) = C2(d) · 1{no q in Q locks d} · ∏_{q∉Q, q∤(6d±1)} pen_q(d).
```

**Result (S₁₀, 23,988,173 twin centres, each centre's real factor set).**
- **210 — resolved.** Second order tracks the observation at every stratum
  (ω=6: observed 0.413, second-order 0.426, first-order 0.357).
- **42 — corrected in direction, residual remains.** Second order turns the
  (wrongly decreasing) first-order curve upward but stalls near 1.14 while the
  observation reaches 1.55.

**Honest scope — single-centre limitation.** The second-order series is by
construction a *single-centre* condition: it conditions on N and treats N+d as
independent. The two centres are correlated mod q. The 210 residual closes
because its lock primes 11,19 rarely divide *both* centres; the 42 residual
remains because the gap's factor 7 shields *both* centres when 7|N — a
**two-centre** effect omitted here. This factor-coincidence shield gain is at
present only a qualitative hypothesis (not a verified closed form); the
two-centre conditional singular series is posed as the open problem.

> No claim is made about the infinitude of twin primes or about any prime
> k-tuple conjecture. This is a conditional, factor-resolved refinement of the
> Hardy–Littlewood gap heuristic, demonstrated empirically with its open remainder.

Part I: Zenodo doi:10.5281/zenodo.20470367 ·
Part II: doi:10.5281/zenodo.20477664 ·
Part III: doi:10.5281/zenodo.20498668

---

## Layout

```
.
├── README.md
├── LICENSE                 (MIT)
├── CITATION.cff
├── data/
│   └── so_plot_S10.csv      6dN, omega, obs, first_order, second_order  (S10)
├── code/
│   ├── second_order.py            scan + first/second-order series; emits so_plot_S{K}.csv
│   └── make_second_order_fig.py   builds the 2-panel obs/1st/2nd showdown figure
├── figures/                fig_paper4_second_order.{pdf,png}
└── paper/                  Chen_6N_Paper4.{tex,pdf} + figure
```

## Reproducing

Requirements: Python 3.8+, `numpy`, `matplotlib`.

```bash
pip install numpy matplotlib

# 1. Scan + evaluate first/second-order series. Default is S10 (~15 min).
#    Prints the obs/1st/2nd table and writes so_plot_S{K}.csv.
python code/second_order.py            # S10
MAXK=9 python code/second_order.py     # S9 (faster, for validation)

# 2. Build the showdown figure (reads ../data/so_plot_S10.csv)
cd code && python make_second_order_fig.py
```

### Definitions / conventions (same as Parts II–III)

- Twin centre: N with 6N−1, 6N+1 both prime. Gap ΔN between consecutive twin
  centres (centre-step units); physical distance 6ΔN. Attributed to ω₍>3₎ of the
  left centre.
- k_q(d): number of fatal residues of d mod q (the classes N ≡ ±6⁻¹, −d±6⁻¹),
  none equal to 0 when q does not lock d.
- Only strata with ≥ 10⁴ gaps are reported (ω = 1…6 in S₁₀).
- **Normalisation:** r(d|ω) is model-internal (each curve ÷ its own ω-merged
  mean). This differs from Part III's global C₂×envelope baseline by a constant
  scale per gap; it does not change the gradient in ω or the residual sizes.
- Engine: complete segmented-sieve factorisation + deterministic interval-sieve
  primality; S₁₀ twin count 23,988,173 matches Part I.

## License

MIT — see `LICENSE`.
