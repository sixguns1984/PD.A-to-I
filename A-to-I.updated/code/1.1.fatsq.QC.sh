#!/usr/bin/bash

fastp -i raw_data.R1.fastq.gz -o clean_data.R1.fastq.gz -I raw_data.R2.fastq.gz -O clean_data.R2.fastq.gz
