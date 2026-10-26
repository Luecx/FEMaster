# FEMaster examples

Generated example suite with 77 input decks. The files are intentionally small and organized by model/analysis family.

Run one deck, for example:

```bash
FEMaster 01_beams/010_beam_cantilever_b33_linear.inp --output 010_beam_cantilever_b33_linear.res
```

These decks were generated from the documented FEMaster DSL and should be treated as first-pass examples. Run them through your current build and adjust element-specific details where your implementation has stricter requirements.

Bolt pretension: [095_bolt_pretension.inp](09_constraints/095_bolt_pretension.inp) demonstrates FORCE, HOLD and a subsequent external load. The [STEP variant](09_constraints/096_bolt_pretension_step.inp) uses the general step syntax. See the [pretension reference](../documentation/pretension.md).
