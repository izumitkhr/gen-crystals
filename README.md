# Crystal structures generated using generative models
This repository contains crystal structure data generated using generative models.
The crystal structures were generated using the following three models:
- [MatterGen (base:Alex-MP-20)](https://github.com/microsoft/mattergen)
- [MatterGen (mp_20_base:MP-20)](https://github.com/microsoft/mattergen)
- [DiffCSP++](https://github.com/jiaor17/DiffCSP-PP)
- [WyckoffTransformer](https://github.com/SymmetryAdvantage/WyckoffTransformer) with DiffCSP++ relaxation

Each generative model was used to generate 5,000 crystal structures. For WyckoffTransformer, the dataset contains 4,998 structures that were successfully generated.

## Reading structures

The generated structures are stored in `data/gen_crystals_<generator>.csv`
(`mattergen`, `diffcsp-pp`, `wyformer`), one CIF string per row in the `cif` column:

```python
import pandas as pd
from pymatgen.core.structure import Structure

df = pd.read_csv("data/gen_crystals_mattergen.csv")
structure = Structure.from_str(df.iloc[0]["cif"], fmt="cif")
```

## Energies

Single-point MLIP energies of the generated structures are provided in `data/energies/`.
See [data/energies/README.md](data/energies/README.md) for the file format and how they were computed.
