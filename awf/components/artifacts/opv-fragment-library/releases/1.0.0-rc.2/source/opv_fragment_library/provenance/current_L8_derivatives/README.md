# L8_4-Cb and L8-Cb4 ztop assembly

This directory combines validated historical NFA fragments with the three new
GAFF2/Multiwfn-RESP side-chain fragments.

The target molecules contain two `C12H25` N-sidechains. Therefore `Y` uses the
historical pyrrole-capped `C46.top/.pdb` with cut site `Y3:6-3`; the retained
component is `C12H25` and is graph-isomorphic to all four target N-sidechains.
The active L8BO library's `C24` fragment retains only `C8H17` and is not used.

- `L8_4-Cb`: `DAUYAVY`
  - first arm: IC (`A`), short cyclobutyl sidechain (`U`), N-sidechain (`Y`)
  - second arm: IC (`A`), long branched sidechain (`V`), N-sidechain (`Y`)
- `L8-Cb4`: `DAWYAXY`
  - first arm: IC (`A`), three-carbon-spacer cyclobutyl sidechain (`W`), N-sidechain (`Y`)
  - second arm: IC (`A`), historical L8BO outer sidechain (`X`), N-sidechain (`Y`)

The arm-wise order follows the validated historical asymmetric build
`DAXYAZY` in `private-source://redacted/stbu.top`. Both builds use
`--total_charge 0`. The raw ztop coordinates are assembly intermediates; final
MD coordinates are taken from the separately audited target SDFs and reordered
to the generated topology after graph-isomorphism validation.

The first `L8_4-Cb` assembly made with the unsuitable `C24` N-sidechain is
preserved under `diagnostics/wrong_C24_N_sidechains/` as a negative case. It
contains 176 rather than 200 atoms and must not be used for MD.
