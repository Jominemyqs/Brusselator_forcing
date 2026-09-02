# Frozen C--B--A challenge and post-budget guarded discovery

## Prespecified eight-run result

The multiclass held-out protocol was frozen before any new fine-grid outcome was
computed. It used 21 labels from the old `0.025`-spaced scan, 121 previously
unsimulated candidates on `0.275 <= lambda <= 0.300`, and four sequential
two-endpoint acquisition rounds. All previously simulated interior refinement
coordinates were excluded without reading their outcomes.

The predefined primary endpoint failed:

```text
new outcomes = 7 A + 1 C + 0 B
primary B-detection success = false
```

The first three rounds moved progressively from the initially inferred C--A
transition toward smaller `lambda`, but every queried point was A. The final
round selected

```text
lambda = 0.2753906250 -> C
lambda = 0.2763671875 -> A
width  = 0.0009765625
```

Only after closing the eight-run budget was the historical pilot opened for
comparison. Its known B point at `lambda = 0.276171875` lies inside the final
ML-selected interval. This is a useful localization result, but it does not
retroactively turn the failed B-detection endpoint into a success.

## Post-budget multiclass refinement

The final ML pair was then passed to a separate, explicitly post-budget guarded
refinement. C and A midpoint outcomes updated their own endpoints; B generated
two brackets; any other or unresolved final outcome would have stopped the
procedure without coercion.

The four additional outcomes were

| Step | Lambda | Outcome | Decision time |
|---:|---:|---:|---:|
| 1 | `0.2758789062` | C | 112 |
| 2 | `0.2761230469` | C | 96 |
| 3 | `0.2762451172` | A | 13 |
| 4 | `0.2761840820` | B | 130 |

This produces

```text
C--B: [0.2761230469, 0.2761840820]
B--A: [0.2761840820, 0.2762451172]
width of each bracket = 6.103515625e-5
```

Thus the full honest result is:

> The frozen eight-run multiclass acquisition did not directly sample the
> narrow B basin. It nevertheless localized a sub-`1e-3` C/A interval
> containing B, after which four guarded PDE simulations recovered the C--B--A
> structure and produced two exceptionally tight edge-tracking brackets.

Twelve new PDE runs were required from the start of the held-out experiment
through B discovery. Their recorded solver runtime was about `288.2 s`.

## Scientific implications

The negative primary endpoint is informative. Pair disagreement and novelty can
efficiently contract a nominal transition interval while still missing a basin
whose intersection with the sampled state family is narrower than the queried
pair spacing. A future active method should therefore include an explicit
multibasin-interior test or adaptive interval subdivision rather than equating
opposite endpoint labels with a direct boundary.

The two recovered brackets now have distinct dynamical roles:

1. The C--B bracket should be edge-tracked first and compared with the already
   verified periodic `E_BC`. Recovery of the same orbit would show that the
   narrow B tongue intersects the stable manifold of the known edge state.
2. The B--A bracket should be edge-tracked without assuming an organizing
   invariant object. Its boundary dynamics may reveal a second edge state or a
   qualitatively different transient mechanism.

Machine-readable outputs are stored in:

- `experiment_outputs/ml_heldout_CBA_geometry_v1`
- `experiment_outputs/ml_heldout_CBA_adaptive_v1`
- `experiment_outputs/ml_heldout_CBA_evaluation_v1`
- `experiment_outputs/ml_heldout_CBA_postbudget_refinement_v1`
- `experiment_outputs/ml_heldout_CBA_stage_summary_v1`
