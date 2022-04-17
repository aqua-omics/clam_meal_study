#Author: David J. Bradshaw, Ph.D. & Nicholas Dickens, Ph.D.
#Emails: djbradshaw2015@fau.edu, dbradshaw3366@gmail.com

#snakemake -s 'filepath to snakefile'

import sys

prefixes = []

with open("manifest.csv") as in_fh:
    # skip the first line
    in_fh.readline()
    for line in in_fh:
        (sample_name, path, direction) = line.split(",")
        (prefix, sep, filename) = path.rpartition("_")
        prefixes.append(prefix)
        # skip the next line as they are all in pairs
        in_fh.readline()

print("Found {} samples.".format(len(prefixes)))




rule all:
    input:
        expand("{pathname}_R1_val_1.fq.gz", pathname=prefixes),
        expand("{pathname}_R2_val_2.fq.gz", pathname=prefixes)


rule fastqc:
    input:
        "{pathname}_R1.fastq.gz",
        "{pathname}_R2.fastq.gz"
    output:
        "{pathname}_R1_fastqc.zip",
        "{pathname}_R2_fastqc.zip"
    threads: 2
    shell:
        "fastqc {input[0]} {input[1]}"


rule trim_primers:
    input:
        "{pathname}_R1.fastq.gz",
        "{pathname}_R2.fastq.gz",
        "{pathname}_R1_fastqc.zip",
        "{pathname}_R2_fastqc.zip"
    output:
        "{pathname}_R1_val_1.fq.gz",
        "{pathname}_R2_val_2.fq.gz"
    shell:
        "trim_galore --fastqc --illumina -q 20 --length 20 --clip_R1 21 --clip_R2 21 -o `dirname {input[0]}` --paired {input[0]} {input[1]}"
