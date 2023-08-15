# Overview

TODO

## 1. Getting Started

## 1.1 Introduction

TODO

The following are sub-commands available:

- initialize: initialize the pipeline
- dryrun: predict the binding of peptides to any MHC molecule
- cluster: execute the pipeline on the Biowulf HPC
- local: execute a local, interactive, session
- git: execute GitHub actions
- unlock: unlock directory
- DAG: create DAG report
- report: create SNAKEMAKE report
- testrun: copies test manifests and files to WORKDIR

## 1.2 Install Dependencies

TODO

## 1.3 Login to the cluster

TODO biowulf, frce, other cloud options?

```
# ssh into cluster's head node
ssh -Y $USER@biowulf.nih.gov
```

## 1.4 Load an interactive session

An interactive session should be started before performing any of the pipeline sub-commands, even if the pipeline is to be executed on the cluster.

```
# Start an interactive node
sinteractive --time=12:00:00 --mem=8gb  --cpus-per-task=4 --pty bash
```
