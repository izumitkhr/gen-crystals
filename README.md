# gen-crystals
This repository contains crystal structure data generated using generative models.
The crystal structures were generated using the following three models:
- [MatterGen](https://github.com/microsoft/mattergen)
- [DiffCSP++](https://github.com/jiaor17/DiffCSP-PP)
- [WyckoffTransformer](https://github.com/SymmetryAdvantage/WyckoffTransformer) with DiffCSP++ relaxation

Each generative model was used to generate 5,000 crystal structures. For WyckoffTransformer, the dataset contains 4,998 structures that were successfully generated.