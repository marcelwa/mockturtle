# Pre-registration — Phase 3: a per-node portfolio resynthesis engine

**Frozen before the first scored job was dispatched.** The sha256 and freeze timestamp of
this file are quoted verbatim in every recipe file of this study and in `NOTES.md`, and the
file is committed and pushed on the mockturtle fork branch `portfolio-resyn` alongside the
driver it describes, so that the freeze is externally checkable.

---

## 0. Where this comes from

Phase 2 (`experiments/2026-08-26-p2-acd-resyn`, pre-registration sha256
`4183dc33cd38ccd6d5d873e6c39f10f787d0f6b3513b6a9d611b377ecac395fd`) evaluated seven
resynthesis engines on one fixed operator:

```
AIG -> lut_map(K=8, area) -> 8-LUT network -> node_resynthesis<aig>(engine) -> AIG
```

scored through a fixed, choice-free `strash; if -K 6 -a; mfs2`. It found that **no single
engine dominates**. Counting per-benchmark mapped winners over 30 paired benchmarks: ACD 10,
DSD 7, SOP-factoring 6, bi-decomposition 4, direct 3. A per-benchmark best-of-all *oracle* is
worth **-24.79 % AIG nodes / -7.05 % mapped 6-LUTs** against DSD alone — but choosing the arm
by AIG size and then scoring mapped captures only **-1.59 %** of it, because the AIG-best
engine is the mapped-best on just **10 of 30** benchmarks.

So the engines are complementary and the *selector* is the entire problem. Phase 3 moves the
selection from the benchmark to the **node**, and — this is the point — uses only criteria
that are computable at resynthesis time, so every arm is a **deployable policy**, not an
oracle. A win here is directly usable.

## 1. The operator, and what is held fixed

Unchanged from Phase 2 in every respect except the engine slot:

- `lut_map(K=8, area_oriented)`, `cut_limit=12`. K = 8 and ACD block size 6 are Phase 2's
  frozen operating points and are **not** re-tuned here.
- The shared fallback is `shannon_resynthesis` down to 4 variables then a **complete**
  4-input NPN database (`xag_complete`).
- Downstream is the identical fixed, choice-free `strash; if -K 6 -a; mfs2`. It is
  deliberately **not** `&dch -f`: Phase 2 section 4 showed `&dch -f` erases the distinction
  between four structurally different front ends (median ratio exactly 1.0000), so a study of
  front ends must not be run behind it.
- Every arm is `category: size`, `start_from: original`, `reference: null`, and every scored
  row must carry a **direct `cec -n`** verdict against the benchmark's own original.
- 31 benchmarks: EPFL 20 + ISCAS'85 11.

**One binary for everything.** `portfolio_resyn` (mockturtle fork branch `portfolio-resyn`).
The `dsd` and `acd` reference arms are *re-run through this binary*, not cited from Phase 2,
and they go through the **identical** scratch-and-replay machinery as the portfolio arms — a
single-engine arm simply evaluates one candidate instead of six. This removes two confounds
at once: the baseline's runtime is the honest 1x reference for the portfolio's overhead, and
no arm differs from another in its code path.

*Verified before freezing* (disclosed as a pilot in section 7): on `cavlc`, `portfolio_resyn
--select=dsd` reproduces Phase 2's `acd_resyn --engine=dsd` **exactly** — 1168 AIG nodes,
depth 22, both — so the new machinery is output-identical to direct emission.

## 2. The six candidates and the three rulers

Candidate engines, in the fixed order that also defines the tie-break priority:

| # | id | what |
|---|---|---|
| 0 | `dsd` | `dsd_resynthesis` over the shared fallback |
| 1 | `acd` | `acd_resynthesis` (block size 6, `use_fallback`) over the shared fallback |
| 2 | `shannon` | the shared fallback alone |
| 3 | `sopf` | `sop_factoring` |
| 4 | `bidec` | `bidecomposition_resynthesis` |
| 5 | `direct` | `dsd_resynthesis` over Shannon(4) + the `aig_complete` NPN database — exactly what `convert_klut_to_graph` does per node, i.e. Phase 2's `direct` arm |

At every node of the 8-LUT network **every** candidate is built into its own throwaway
scratch AIG (`num_vars` PIs, one PO), `cleanup_dangling`-ed (which strashes it exactly as the
replay into the destination will), and measured three ways:

- **size** — AIG nodes emitted for that node;
- **depth** — AIG depth of the emitted block, from its leaves;
- **luts** — 6-LUT count of an area-oriented `lut_map` of that block **in isolation**.

`luts` is the *derived, engine-uniform* equivalent of the cost ACD reports natively (its
cascade block count; also recorded per node in the CSV as `acd_blocks` so the two rulers can
be compared post hoc). A per-engine structural LUT count is not definable at all for `sopf`
or `bidec`, so a uniform ruler is the only way six arms can be compared, and mapping the
block in isolation is the closest computable thing in the endpoint's own units.

**Stated in advance rather than discovered afterwards:** a block whose function has <= 6
variables costs exactly one LUT for *every* engine. So `p-lut` can differ from `p-area` only
on nodes with **more than 6 inputs**, and everywhere else it degenerates to its tie-break.
The driver reports `calls_wide` (nodes with > 6 inputs) and the per-node CSV records
`num_vars`, so this fraction is a reported quantity, not an assumption. On the `cavlc` pilot
it was 28 of 38 nodes; on `log2` 1823 of 3789; on `hyp` 9469 of 42833.

## 3. The arms

| arm | selects, per node, the candidate minimising |
|---|---|
| `p-area` | (size, depth, luts, engine-index) |
| `p-depth` | (depth, size, luts, engine-index) |
| `p-lut` | (luts, size, depth, engine-index) |
| `dsd` | — single-engine reference |
| `acd` | — single-engine reference |

**Tie-break rule.** Each arm ranks candidates by the 4-tuple above, compared
lexicographically. The last component is the engine index of section 2, so **no two
candidates ever compare equal** and the selection is deterministic and reproducible. Note the
consequence, declared here: on the <= 6-variable nodes where all `luts` tie, `p-lut` falls
through to AIG size and is therefore *identical to* `p-area` on those nodes. The per-node CSV
records every candidate's three costs, so any other tie-break (engine-index only, for
instance) can be evaluated post hoc from the same run without re-running anything — and any
such re-analysis is exploratory by construction.

5 arms x 31 benchmarks = **155 scored cells**, dispatched `--jobs 8` on the local 16-core
box. A separate **serial** (`--jobs 1`) sub-run on a 10-benchmark subset measures runtime,
because wall times from a parallel sweep are not comparable.

## 4. Hypotheses, and which selector I predict wins

**H1 (headline).** `p-lut` beats the `dsd` reference on **mapped 6-LUT count**.

**H2 (the reason `p-lut` is in the list).** `p-lut` beats `p-area` on **mapped 6-LUT count**.
This is the substantive claim: that a local ruler in the *endpoint's own units* predicts the
mapped result better than AIG size does. Phase 2 supplies the motivating evidence — `sopf`
produces by far the smallest AIGs of any engine and still loses to ACD on mapped LUTs by
4.7 % geomean — so AIG size is already known to be a poor proxy at benchmark granularity.

**H3.** `p-depth` beats the `dsd` reference on **mapped levels**.

**Predicted winner on mapped 6-LUT count: `p-lut`.** Reasoning, so it can be held against
me: mapped LUT count is what we are scored on; of the three computable local rulers, `luts`
is the only one measured in that unit; and Phase 2 established that the AIG-size ruler is
actively misleading at the coarser granularity. I predict `p-area` will produce much the
smallest AIGs and will *not* convert that into mapped LUTs — quite possibly it will be
*worse* than `dsd` on LUTs. I predict `p-depth` loses on LUTs and is the only arm with a
chance on levels.

**The honest counter-argument, stated in advance.** All three rulers are *local and
additive*, and the endpoint is neither: `if -K 6 -a` re-maps globally over the whole AIG and
can merge across the node boundaries these blocks were measured inside, and it exploits
sharing between nodes that no per-node measurement can see. So it is entirely possible that
none of the three rulers moves the endpoint. That outcome is F1 below and is a publishable
negative, not a failed session.

## 5. Statistics and decision thresholds

**Primary statistic:** geometric mean of per-benchmark ratios (arm / reference), over the
benchmarks paired in both arms, reported with a **leave-one-out range** and the mover.

**Secondary statistic:** total-sum ratio, always printed with **`hyp`'s share of the
denominator** beside it (`hyp` is 47.1 % of the suite LUT total, so a total-sum statement is
close to a statement about `hyp` alone), also with leave-one-out.

**Inference:** two-sided Wilcoxon signed-rank on the paired per-benchmark counts, alpha = 0.05.

**A hypothesis is SUPPORTED iff geomean <= 0.99 AND p < 0.05**, on >= 24 paired benchmarks.
Anything else is not supported. Both statistics are reported for every comparison whatever
the verdict.

Metrics reported for **every** arm, always: AIG nodes, AIG depth (both from the driver's own
telemetry line, since the harness scores only mapped BLIFs), mapped 6-LUT count, mapped
levels (both from the harness), and runtime.

## 6. Falsification conditions

Written down now so that no post-hoc reading of the data can escape them.

- **F1.** If all three portfolio arms land inside a geomean of [0.99, 1.01] against `dsd` on
  mapped LUTs, the claim "a per-node local ruler can steer engine selection toward the mapped
  endpoint" is **falsified**, and the finding is that local rulers do not predict the global
  mapped result. It goes in `knowledge/negative-results.md` with that mechanism.
- **F2.** If `p-area` <= `p-lut` on mapped LUTs, the specific premise behind `p-lut` — that
  matching the ruler's units to the endpoint's is what matters — is **falsified**, whatever
  the absolute numbers do.
- **F3.** If any portfolio arm's selection counts are >= 95 % a single engine, that arm is
  **degenerate**: whatever it shows is a statement about that engine, not about a portfolio,
  and must be reported as such. Selection counts are therefore reported for every arm before
  any QoR claim is interpreted.
- **F4.** If the correctness gate (section 8) does not pass with **zero** honest mismatches
  **and** a sabotage control that fires on every injected break, **no QoR number from this
  study is reported at all**.
- **F5.** Runtime is a reported metric, not an afterthought. If the portfolio's resynthesis
  cost exceeds 100x the single-engine cost, the policy is not deployable at suite scale
  regardless of what it does to LUT count, and that is reported as the finding.

## 7. Disclosed pilots — exploratory, never to be quoted as results

Run before this file was frozen, to check the driver works at all:

1. `cavlc`, all five arms, unmapped: AIG sizes 1168 (`dsd`), 1298 (`acd`), 788 (`p-area`),
   839 (`p-depth`), 1194 (`p-lut`); selection counts observed; all five CEC-equivalent to the
   original by `cec -n`.
2. `log2` and `hyp`, `dsd` vs `p-lut`, unmapped, for timing: resynthesis 1.50 s -> 16.47 s
   (log2) and 6.74 s -> 74.73 s (hyp), i.e. ~11x, not the ~6x a naive count of engines
   suggests; and on `hyp` the portfolio's resynthesis (74.7 s) **exceeds** the shared mapping
   step (16.1 s), so Phase 2's "resynthesis is negligible behind mapping" does not survive to
   this operator. AIG sizes seen: `log2` 75300 (`dsd`) -> 80787 (`p-lut`); `hyp` 338657
   (`dsd`) -> 321444 (`p-lut`).

**No mapped 6-LUT number has been observed for any arm of this study.** The primary endpoint
is blind at the moment of freezing.

## 8. The correctness gate, run before any QoR job

Follows the protocol that already caught a vacuous check in this project.

- **Part 1 — honest.** With `--check`, **every** candidate of **every** engine at **every**
  node of all 31 benchmarks is simulated **exhaustively** over all 2^num_vars input patterns
  and compared against the function `node_resynthesis` asked for. Any mismatch is fatal and
  aborts the run.
- **Part 2 — sabotage, which must fire.** 20 seeds x 2 modes.
  - `minterm` (the strong one): the candidate is resynthesised from a copy of the target
    truth table with **one bit flipped**, so the emitted block is wrong on exactly one of the
    2^n patterns. That is the weakest possible error and is precisely what a sampling checker,
    or one comparing against the wrong reference, would miss — it tests that the check is
    genuinely exhaustive, not merely that it compares.
  - `invert`: the emitted signal is complemented inside the scratch network before it is
    measured, simulated and compared.
  - A run passes only if it injected **at least one** break and caught **all** of them; a
    control that never fires proves nothing and is failed explicitly.
  - A sabotage run really does build wrong circuits, so the driver **writes no output** in
    sabotage mode and it can never be scored by accident.

## 9. What is exploratory in this study

Everything not listed in sections 3 and 4. In particular: any re-analysis of the per-node
CSVs (alternative tie-breaks, per-width selection statistics, oracle bounds computed from the
recorded candidate costs), and any arm added after the dispatch. Exploratory results are
labelled as such everywhere they appear and require their own pre-registration before they
are quoted as results.

## 10. Explicitly out of scope

Emitting all six candidates as *structural choices* for the mapper to select from, rather
than picking one per node. That is the natural follow-up and is deliberately **not** started
here.
