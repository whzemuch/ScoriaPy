# ScoriaPy

ScoriaPy is a lightweight utility library built on top of Scanpy.
It organizes helper functions following the Scanpy namespace design:

- `scoriapy.pp` — preprocessing utilities
- `scoriapy.tl` — tools/analysis utilities
- `scoriapy.pl` — plotting utilities
- `scoriapy.utils` — generic helpers (logging, markers, etc.)

This package is intended for flexible, incremental development of Scanpy pipelines.

## Documentation

https://scoriapy.readthedocs.io/en/latest/

## Installation

```bash
pip install git+https://github.com/whzemuch/scoriapy
# upgrade later
pip install -U git+https://github.com/whzemuch/scoriapy
```

## History

- add run_jasmine function to calculate jasmine scores. Proj_202511_SunL_AgingSC_subset_Fibro_T/RankGenesGroups_adata_4zellkconverter_20260129.ipynb
