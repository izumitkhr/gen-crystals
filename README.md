# Crystal structures generated using generative models
This repository contains crystal structure data generated using generative models.
The crystal structures were generated using the following three models:
- [MatterGen](https://github.com/microsoft/mattergen)
- [DiffCSP++](https://github.com/jiaor17/DiffCSP-PP)
- [WyckoffTransformer](https://github.com/SymmetryAdvantage/WyckoffTransformer) with DiffCSP++ relaxation

Each generative model was used to generate 5,000 crystal structures. For WyckoffTransformer, the dataset contains 4,998 structures that were successfully generated.

## Accessing data
To read a structure from the dataset (e.g., `gen_crystals_mattergen.csv`), where each row contains a CIF string in the `cif` column, you can use the following code:

```python
import pandas as pd
from pymatgen.core.structure import Structure

df = pd.read_csv("gen_crystals_mattergen.csv")
cif = df.iloc[0]["cif"]
structure = Structure.from_str(cif, fmt="cif")
