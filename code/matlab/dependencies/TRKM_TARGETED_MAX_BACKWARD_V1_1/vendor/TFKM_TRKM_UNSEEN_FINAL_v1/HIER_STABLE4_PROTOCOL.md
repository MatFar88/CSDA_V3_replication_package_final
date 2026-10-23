# HIER_STABLE4 full-pipeline development protocol

## Frozen challenger before this run

Primary TFKM challenger: `HIER_STABLE4`.
Sensitivity comparators from the fixed-triplet audit are not eligible to replace
it after seeing the full-pipeline result.

## Estimator change

At fixed `(r,g,alpha)`, all 100 frozen initializations are fully fitted by the
same corrected alternating TFKM algorithm. For each initialization family, the
terminal fit maximizing within-family

`ARI(mutually retained labels) + 0.2 Jaccard(trim sets)`

is retained (objective only breaks exact ties). The at-most-four family winners
are then compared by the same similarity; the most stable representative is the
candidate fit. Its final four-winner consensus score becomes TFKM native evidence
`N` for the unchanged `RANK_EN` selector.

The unequal 50/35/10/5 start allocation is therefore an exploration budget, not
a family voting weight.

## Downstream freeze

No changes are allowed to:

1. `RANK_EN + GMAX_R_E`
2. `G_STRICT_W`
3. `R_STRICT_NR2`
4. `GAIN_OBJ_L0800`
5. `RG_VOTE_THEN_ALPHA`
6. ridge, DGP, grids, seeds, or trimming rule.

TRKM is unchanged.

## Development gate

Run the 16 consumed v11 clean/bridge units first. Report, without retuning:

- ARI delta vs current;
- structural over-selection rate;
- over-r and over-g separately;
- any under-r/under-g rate;
- exact `(r,g)`;
- selected-family mix.

The purpose is mechanistic diagnosis: determine whether family-neutral basin
arbitration alters the structural selector through improved `E` and `N`, with
special attention to TFKM over-selection.

## Validation status

These results cannot be called fresh validation because HIER_STABLE4 was chosen
using v11 fixed-triplet audit outcomes. Adoption requires a new disjoint panel.
