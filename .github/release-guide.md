# Release Guide

## How to test a pre-release on biowulf

Install the development version of champagne.

```sh
# go to the source on biowulf
cd /data/CCBR_Pipeliner/Pipelines/CHAMPAGNE/champagne-dev
git pull
# optionally switch to different branch if needed
# make sure you have python >= 3.10 in your path
# install the version to a hidden path (e.g. .dev, .v1.0.0.9000) in /data/CCBR_Pipeliner/Pipelines/CHAMPAGNE
cd ..
pip install ./champagne-dev -t ./.dev
# if the .dev directory already exists you'll need to use the --upgrade flag
# pip creates a binary file .dev/bin/champagne
```
