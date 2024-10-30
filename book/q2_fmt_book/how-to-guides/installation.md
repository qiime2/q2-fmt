## Installation

Please note that this software is now an *alpha release*. Initial conda packages are available, but have not been tested in integration against our other plugins yet. For now, create a fresh  `Development` conda environment for q2-fmt using the following command: 

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

