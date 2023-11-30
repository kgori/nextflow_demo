#!/bin/bash

nextflow run nextflow/main.nf \
    --from_reference=data/refs/canfam3.fa \
    --to_reference=data/refs/ros1.0.fa \
    --inputDir=data/crams \
    --outputDir=output \
    -resume
