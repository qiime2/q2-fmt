# q2-fmt (fecal microbiota transplant)
Note: This software is in **Alpha release** and interfaces are subject to change.
![](https://github.com/qiime2/q2-fmt/actions/workflows/ci-dev.yaml/badge.svg)

## Demo
[![](https://raw.githubusercontent.com/qiime2/q2-fmt/master/demo/screenshot.png)
**Interactive Link**](https://view.qiime2.org/visualization/?type=html&src=https%3A%2F%2Fraw.githubusercontent.com%2Fqiime2%2Fq2-fmt%2Fmaster%2Fdemo%2Fraincloud-baseline0.qzv)

## Installation

Please note that this software is now an alpha release. Initial conda packages are available, but have not been tested in integration against our other plugins yet. For now, create a fresh conda environment for q2-fmt using the following command: 

Note: This is a `Development` install of q2-fmt
Mac OS instructions 
```bash
CONDA_SUBDIR=osx-64 conda env create \
 -n q2-fmt-2024.10 \
 -f https://raw.githubusercontent.com/qiime2/q2-fmt/dev/environment-files/2024.10-q2-fmt-environment.yml
```

Linux instructions
```bash
conda env create \
 -n q2-fmt-2024.10 \
 -f https://raw.githubusercontent.com/qiime2/q2-fmt/dev/environment-files/2024.10-q2-fmt-environment.yml
```

Then activate your new environment as usual.
```bash
conda activate q2-fmt-2024.10
```
