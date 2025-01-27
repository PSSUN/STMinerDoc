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
  pip3 install STMiner
  ```

  :::{dropdown} NOTE: If the download speed is slow, please try to specify the source
  For example:
  ```bash
  pip3 install STMiner -i https://pypi.tuna.tsinghua.edu.cn/simple
  ```
  :::

## STMiner can be cloned from GitHub and installed locally (NOT recommended)
Download STMiner [here](https://github.com/xjtu-omics/STMiner.git) and unzip, then run:

```bash
pip3 install -r requirements.txt
```
