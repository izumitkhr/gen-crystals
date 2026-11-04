# MLIP energies of the generated structures

Paths below are relative to the repository root.

## Reading energies

`data/energies/<model>/gen_crystals_<generator>_energies.csv` holds single-point energies of the
structures in `data/gen_crystals_<generator>.csv`, computed with a machine-learning interatomic
potential (MLIP). `index` is the 0-based row number of the corresponding structure file.

| `<model>` | MLIP | Level of theory | Correction |
|---|---|---|---|
| `mace-omat-0-medium` | MACE-OMAT-0 (medium) | PBE/PBE+U (OMat24) | MP2020 |
| `mace-matpes-r2scan-0` | MACE-MATPES-r2SCAN-0 (MACE-OMAT-0 fine-tuned on MatPES r2SCAN) | r2SCAN | none (raw energies) |

Energies from different models are on different absolute scales (r2SCAN total energies are
several eV/atom lower than PBE), so never compare or mix them across `<model>` directories.

The main columns are:

- `status`: `ok` if the structure was evaluated. `too_close` marks structures with two atoms
  closer than 0.7 Å, which were skipped to keep the MLIP from being evaluated on overlapping
  atoms; their energy columns are empty. This is a numerical guard, not a validity filter, so
  apply your own validity criteria if you need them.
- `uncorrected_energy_total_eV`: the raw MLIP total energy.
- `mp2020_correction_eV`: the Materials Project MP2020 correction for this composition (anion
  and GGA/GGA+U mixing corrections, with the run type and Hubbard U that pymatgen's MPRelaxSet
  would assign). Zero for compositions that need no correction. MP2020 is calibrated for
  PBE/PBE+U only, so for `mace-matpes-r2scan-0` no correction is applied and this column, as well
  as `mp2020_adjustments`, `run_type` and `hubbards`, is empty.
- `energy_total_eV`, `energy_per_atom_eV`: the corrected energy, i.e. the sum of the two columns
  above (divided by `nsites` for the per-atom value); equal to the raw MLIP energy when no
  correction is applied. Use these for convex-hull analysis, but
  build the reference hull with the same MLIP and the same corrections. These are MLIP
  estimates, not DFT energies, so do not mix them with a hull built from DFT data.
- `composition`: JSON dict of element counts, so phase-diagram entries can be built without
  parsing the CIF.

To join the energies with the structures:

```python
import json
import pandas as pd
from pymatgen.core.composition import Composition

structures = pd.read_csv("data/gen_crystals_mattergen_base.csv")
energies = pd.read_csv("data/energies/mace-omat-0-medium/gen_crystals_mattergen_base_energies.csv")
df = structures.merge(energies[energies["status"] == "ok"], left_index=True, right_on="index")
composition = Composition(json.loads(df.iloc[0]["composition"]))
```

## Computing energies

The energies are produced by `scripts/compute_mlip_energies.py`, which needs only the packages
in `scripts/requirements-mlip.txt` (Python 3.12, CUDA GPU recommended):

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r scripts/requirements-mlip.txt
mkdir -p models && curl -L -o models/mace-omat-0-medium.model \
  https://github.com/ACEsuit/mace-foundations/releases/download/mace_omat_0/mace-omat-0-medium.model
sha256sum models/mace-omat-0-medium.model   # d4b14be9afa294eebdbe31a0280b26a0fa29715771e978cbfd3ca24e0d90307a
python scripts/compute_mlip_energies.py --model models/mace-omat-0-medium.model

curl -L -o models/mace-matpes-r2scan-0.model \
  https://github.com/ACEsuit/mace-foundations/releases/download/mace_matpes_0/MACE-matpes-r2scan-omat-ft.model
sha256sum models/mace-matpes-r2scan-0.model  # 8f147ecffa1d06a696e5648b13b81abda9c96fda9d0fb3faba3cd3fb27d9bba9
python scripts/compute_mlip_energies.py --model models/mace-matpes-r2scan-0.model --correction none
```

The model file name (without `.model`) becomes the output directory name and is used to look up
the expected SHA256, so keep the names above. The r2SCAN release does not publish a checksum; the
value above is that of the file downloaded on 2026-11-04 (79,470,738 bytes, matching the release
asset size). `--correction` was added after the `mace-omat-0-medium` data were produced, so that
run's `manifest.json` records the SHA256 of the earlier script version; with the default
`--correction mp2020` the current script gives the same results.
