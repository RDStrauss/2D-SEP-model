# Notes from the Bochum benchmark run

Findings from benchmarking this code against an independent reimplementation
(JAX, 2-D focused transport, same grid and same scheme) on the Dröge et al.
(2014) 7 February 2010 event and two modern multi-spacecraft events.

The commits on this branch are ordered so that a **prefix can be merged**:
everything up to and including `perf: fission the k-regions` is verified
**bitwise identical** at 1, 2, 4, 8, 16, 32, 48 and 96 threads. Only the last
commit changes results.

| commit | effect | verified |
|---|---|---|
| `omp: do not hardcode the thread count` | none | bitwise |
| `bug: the phi-sweep upwind test read an undefined j` | none today | bitwise |
| `perf: interchange the two i/j/k mu-regions` | none | bitwise |
| `perf: parallelise the five f_old = f_new copies` | none | bitwise |
| `perf: fission the k-regions, collapse over (k,j)` | none | bitwise |
| `physics: evaluate G on the mu cell faces` | **changes results** | conservation test |

Performance, on one exclusive node (96 cores, gcc 12.4.0, `-O3 -fopenmp`):

| configuration | before | after | note |
|---|---|---|---|
| 121×108×49, 96 threads | 31 s | **11 s** | |
| 91×81×25 (Dröge), best | 595 s | **174 s** | at **48** threads, not 96 |

The optimised build **regresses past 48 threads** on the smaller grid (174 s →
244 s at 96). On that grid the second socket costs more than it buys.

---

## Two things we found and did NOT change

Both are real, both are in the diagnostics rather than the solver, and both are
arguably your call rather than ours — so they are documented here instead of
patched.

### 1. The observer cells are hardcoded

```fortran
indeks_r   = (1. - R(1))/delta_r + 1.                     ! line 718 etc.
indeks_phi = (33./180.*PI - phi(1) + 82./180.*PI)/delta_phi + 1.
```

Every vantage is placed at **r = 1 AU** with a fixed longitude offset (82° for
unit 801, 120° for 802, −21° for 803). Fine as a demo; wrong for reproducing a
specific event. Against the real 2010-02-07 geometry this put Earth **13.3°**
and STEREO-A **17.8°** away from where our code sampled — the two codes were
reporting intensities at different points in the heliosphere, which was most of
the residual disagreement after the μ-face fix above.

Two details for anyone changing this:

* Fortran `INTEGER` assignment **truncates**. An index expression evaluating to
  20.9 selects cell 20, and on this gradient that alone is a ~5% intensity
  error. Emit rounded integer literals, or use `NINT`.
* The `omni_time_best_connection.txt` file (unit 800) is a *different* vantage
  from `omni_time_stereo_B.txt` (unit 803) — different longitude, different
  slope (−0.0715 vs −0.1346 dex/h). We compared the wrong one for a day.

### 2. The inner cell is inflated by 1/R(1)² ≈ 400

In the r-convection block:

```fortran
f_old(1,j,k) = 1000.*exp(...)/time*exp(...)   ! line 156, injection BC, UNSCALED
f_old(i,j,k) = f_old(i,j,k)*R(i)*R(i)         ! line 176, forward,  ALL i
f_old(i,j,k) = f_new(i,j,k)/R(i)/R(i)         ! line 272, backward, ALL i
```

The back-conversion divides `f_new(1)`, which still holds the **unscaled** Reid
value the injection BC wrote, so

```
f_old(1) = Reid / R(1)^2 = Reid / 0.0025 = 400 x Reid
```

freshly, every step. Measured against our code with observers matched: peak
offset **398.8, 398.2, 378.1×** at r = 1.039, 0.989, 0.940 AU. It is an
amplitude offset with no shape component — decay rates agree to 0.02–1.1% once
observers match and the background is subtracted — but it makes an *absolute*
comparison impossible, and the 5% spread across vantages suggests the round-trip
is not the whole mechanism.

**We could not fix it.** Re-asserting the Dirichlet inner boundary after the
conversion compiles, runs clean and changes nothing, because line 156 rewrites
`f_old(1)` at the top of the next timestep regardless. An earlier attempt that
seeded `f_new` at *both* boundaries produced NaN from the first output step:
`f_new(N)` is live input to the μ operators that follow. The fix has to be in
the injection BC, and that is your design decision.

---

## Two things we tried that did not work

Recorded so nobody repeats them.

* **Pointer rotation instead of the copies.** Three separate attempts, all
  failed a bitwise check. The `DO i = 2, N` bound on the copies protects the
  injection boundary from the μ operators — it is not incidental. That 34% of
  the 96-thread runtime is still available to whoever gets it right.
* **An explicit absorbing outer boundary.** The block at lines 254–262 is
  commented out, so nothing imposes a condition at `r_max`; the outer plane is
  nonetheless updated every step by the μ operators. Adding an explicit
  absorbing BC **changed the answer**, so we left it alone. Worth knowing that
  the two treatments will diverge on any run where particles populate the outer
  domain. (Measured on the Dröge event: moving r_max from 4.5 to 12 AU with dr
  held fixed changes the 20-hour intensity by only 9%, converged by 6 AU.)

## Two configuration values worth a second look

Not changed here, since they are yours to set.

* **`h = 0.5`** in the D_μμ resonance gap at line 975,
  `F = (1-mu**2)*((ABS(MU))**(1.67-1.) + 0.5)`. Dröge et al. (2014), whose
  `lambda = 0.08` and `alpha = 0.13` this file also hardcodes, use the
  Bieber/Dröge gap **h = 0.05** — and so does your own `SEP_propagator`
  (`SEP_propagator.f90:629`) and the Primer (van den Berg, Strauss &
  Effenberger 2020, §3.3.3). On the 2010-02-07 event, 0.5 → 0.05 improves the
  shape agreement at the best-connected spacecraft by a factor 2.3
  (log-RMSD 0.156 → 0.069).
* **The sign in the focusing-length generator.** This code *reads* 1/L from
  `Parker_Focus_108.dat`, and the repository ships no `.dat` files, so we
  cannot tell which form yours used:
  `1/L = cos(psi) * [2/r ± (Omega/V)^2 r / (1 + tan^2 psi)]`.
  The **minus** sign is the physical one (the spiral makes B fall off more
  slowly outward, weakening focusing; the large-r limit then gives the textbook
  L(1 AU) ≈ 1 AU). The plus form over-focuses by **1.7×** at 1 AU (L = 0.577 vs
  0.98 AU). Everything above uses the minus form.

---

Full write-up, figures and the conservation test:
contact Frederic Effenberger, Ruhr-Universität Bochum.
