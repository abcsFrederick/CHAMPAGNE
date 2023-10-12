#!/usr/bin/env bash

version=$1

repo_path=/data/CCBR_Pipeliner/Pipelines/CHAMPAGNE/champagne-dev/
install_path=/data/CCBR_Pipeliner/Pipelines/CHAMPAGNE/${version}
bin_path=${install_path}/bin/

. "/data/CCBR_Pipeliner/db/PipeDB/Conda/etc/profile.d/conda.sh"
conda activate py311

echo "Installing CHAMPAGNE to ${install_path}"
pip install ${repo_path} --target ${install_path}

if [[ ":$PATH:" != *":${bin_path}:"* ]];then
    export PATH="${PATH}:${bin_path}"
fi

if [[ ":$PYTHONPATH:" != *":${install_path}:"* ]];then
    export PYTHONPATH="${PYTHONPATH}:${install_path}"
fi
