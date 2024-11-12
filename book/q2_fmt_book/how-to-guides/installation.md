# Installation

Please note that this software is now an *alpha release* and interfaces are subject to change. 
For now, create a fresh  `Development` conda environment for q2-fmt using the following command: 

These installation installs will have the classic qiime2 amplicon env and q2-fmt.


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

## Updating an existing env to have q2-fmt 

If you have an existing qiime2 env that you want to install q2-fmt in, please run the following command: 
```bash 
conda activate <env-name> # conda env you wish to install this plugin into

conda env update --file https://raw.githubusercontent.com/qiime2/q2-fmt/dev/environment-files/2024.10-q2-fmt-environment.yml
```

Awesome! Now, you are ready to start assessing engraftment extent with q2-fmt!
