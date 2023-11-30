#!/bin/bash
REF=$1
BAM=$2

samtools depth --reference $REF -G 0x900 $BAM
