#!/bin/bash

nextflow run nextflow/main.nf \
    --from_reference=data/canFam3.fa \
    --to_reference=data/ros1.0.fa \
    --inputDir=data/crams \
    -resume
