# Installation

## pip (✨recommended)
STMiner has been uploaded to [PyPI](https://pypi.org/project/STMiner), You can install STMiner via pip:

 1. Create conda environment:

  ```bash
  conda create -n stminer python=3.10
  ```

 2. Activate the environment:
    
  ```bash
  conda activate stminer
  ```

 3. Install STMiner via pip:

  ```bash
  pip install stminer
  ```

  :::{dropdown} NOTE: If the download speed is slow, please try to specify the source
  For example:
  ```bash
  pip install stminer -i https://pypi.tuna.tsinghua.edu.cn/simple
  ```
  :::

## Install a matching source checkout
Download STMiner [here](https://github.com/xjtu-omics/STMiner.git) and unzip, then run:

```bash
python -m pip install -e .
```

Run this in the STMiner repository, not the documentation repository. The current
source dependencies include `gprofiler-official==1.0.0`. New pattern-enrichment
and spatial-pathway interfaces require a checkout containing those methods;
an older PyPI installation may not include them. Online enrichment requires
internet access; plotting an existing enrichment table does not make additional
enrichment requests.
