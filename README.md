# FDwavekit

## Setup environments

1) Create the base FDwavekit environment from the repo root:

```bash
conda env create -f environment.yml
```

2) Create the OpenSWPC environment inside `OpenSWPC-5.2.0`:

```bash
cd OpenSWPC-5.2.0
conda env create -f environment.yml
conda activate openswpc
cd src
bash make_psv.sh
```
