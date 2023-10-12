# Release Guide

## How to test a pre-release on biowulf

Install the development version of champagne.

```sh
# activate the conda env for development
. "/data/CCBR_Pipeliner/db/PipeDB/Conda/etc/profile.d/conda.sh"
conda activate champagne

# go to the source on biowulf and update
cd /data/CCBR_Pipeliner/Pipelines/CHAMPAGNE/champagne-dev
git pull
# optionally switch to different branch if needed

# install the version to a hidden path (e.g. .dev, .v1.0.0.9000) in /data/CCBR_Pipeliner/Pipelines/CHAMPAGNE
cd ..
pip install ./champagne-dev -t ./.dev
# pip creates a binary file .dev/bin/champagne
# add it to your path with:
./add_to_path .dev
```
