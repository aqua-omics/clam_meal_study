"""
Snakemake Snakefile for Running Dada2, based on:
https://github.com/LangilleLab/microbiome_helper/wiki/Amplicon-SOP-v2

source activate qiime2-2020.6
"""

#Authors: David J. Bradshaw, Ph.D. & Nicholas Dickens, Ph.D.
#Emails: dbradshaw2015@fau.edu, dbradshaw3366@gmail.com

#snakemake -s '/home/microbiology/NICK/forwards.snakefile'  --configfile '/home/microbiology/NICK/forwards.json'

import multiprocessing
import subprocess
import qiime2
import pandas as pd

cpu_count=int(multiprocessing.cpu_count())


#rule all:
#  input:
#          file=config['output_dir'] + "/reads_trimmed.qzv"

#rule all:
#  input:
#      file=config['output_dir'] + "/denoising-stats.qzv"

#rule all:
#  input:
#      biom=config['output_dir'] + "/phyloseq.biom"

#rule all:
#  input:
#      directory(config['output_dir'] + "/species/table_level_7_aa")

rule all:
  input:
      tsv=config['output_dir'] + "/species/level_7_aa.tsv"


#TODO: add a primer trimming step, primers can be added to the config.json


# 1.5 Import fastq files (with primers trimmed)
rule import_raw_sequences:
    input:
      manifest=config["manifest_file"]
    output:
      seqs=config['output_dir'] + "/reads_trimmed.qza"
    params:
      type="SampleData[PairedEndSequencesWithQuality]",
      input_format="PairedEndFastqManifestPhred33V2"
    shell:
      "qiime tools import --type {params.type} --input-path {input.manifest} --output-path {output.seqs} --input-format {params.input_format}"

# Visualization for 1.5
rule summarize_imported_seqs:
    input:
      seqs=config['output_dir'] + "/reads_trimmed.qza"
    output:
          file=config['output_dir'] + "/reads_trimmed.qzv"
    shell:
          "qiime demux summarize --i-data {input.seqs} --o-visualization {output.file}"

# 2 - Dada2

rule run_dada2:
    input:
      seqs=config['output_dir'] + "/reads_trimmed.qza",
      file=config['output_dir'] + "/reads_trimmed.qzv"
    output:
      table=config['output_dir'] + "/table.qza",
      seqs=config['output_dir'] + "/rep_seqs.qza",
      stats=config['output_dir'] +"/denoising-stats.qza"
    shell:
      """qiime dada2 denoise-paired --i-demultiplexed-seqs {input.seqs} \
                       --p-trunc-len-f 225 \
                       --p-trunc-len-r 225 \
                       --p-n-threads 4 \
                       --o-table {output.table} \
                       --o-representative-sequences {output.seqs} \
                       --o-denoising-stats {output.stats}"""

# Visualization for step 2
rule summary_table:
    input:
      table=config['output_dir'] + "/table.qza",
    output:
      table=config['output_dir'] + "/table_summary.qzv"
    shell:
      "qiime feature-table summarize --i-table {input.table} --o-visualization {output.table}"

rule summary_seqs:
    input:
      table=config['output_dir'] + "/table_summary.qzv",
      rep_seqs=config['output_dir'] + "/rep_seqs.qza",
      stats=config['output_dir'] +"/denoising-stats.qza"
    output:
      file=config['output_dir'] + "/rep_seqs_summary.qzv"
    shell:
      "qiime feature-table tabulate-seqs --i-data {input.rep_seqs} --o-visualization {output.file}"

rule summary_stats:
    input:
      file=config['output_dir'] + "/rep_seqs_summary.qzv",
      stats=config['output_dir'] +"/denoising-stats.qza"
    output:
      file=config['output_dir'] + "/denoising-stats.qzv"
    shell:
      "qiime metadata tabulate --m-input-file {input.stats} --o-visualization {output.file}"

rule annotate_rep_seqs:
    input:
      file=config['output_dir'] + "/denoising-stats.qzv",
      rep_seqs=config['output_dir'] + "/rep_seqs.qza",
      classifier=config["classifier"]
    output:
      taxonomy=config['output_dir'] + "/taxonomy.qza"
    shell:
      "qiime feature-classifier classify-sklearn --i-classifier {input.classifier} --i-reads {input.rep_seqs} --o-classification {output.taxonomy} --p-reads-per-batch 1000 --p-n-jobs 1 --verbose"

rule visualize_taxonomy:
    input:
      taxonomy=config['output_dir'] + "/taxonomy.qza"
    output:
      file=config['output_dir'] + "/taxonomy.qzv"
    shell:
      "qiime metadata tabulate --m-input-file {input.taxonomy} --o-visualization {output.file}"

rule filter_mito_from_table:
    input:
      file=config['output_dir'] + "/taxonomy.qzv",
      table=config['output_dir'] + "/table.qza",
      taxonomy=config['output_dir'] + "/taxonomy.qza"
    output:
      table=config['output_dir'] + "/table_nmc.qza"
    shell:
      "qiime taxa filter-table --i-table {input.table} --i-taxonomy {input.taxonomy} --p-exclude mitochondria --o-filtered-table {output.table}"

rule filter_mito_from_rep_seqs:
    input:
      rep_seqs=config['output_dir'] + "/rep_seqs.qza",
      table=config['output_dir'] + "/table_nmc.qza"
    output:
      rep_seqs=config['output_dir'] + "/rep_seqs_nmc.qza",
    shell:
      "qiime feature-table filter-seqs --i-data {input.rep_seqs} --i-table {input.table} --o-filtered-data {output.rep_seqs}"

rule visualize_nmc_table:
    input:
      rep_seqs=config['output_dir'] + "/rep_seqs_nmc.qza",
      table=config['output_dir'] + "/table_nmc.qza"
    output:
      file=config['output_dir'] + "/table_nmc.qzv"
    shell:
      "qiime feature-table summarize --i-table {input.table} --o-visualization {output.file}"

rule visualize_nmc_rep_seqs:
    input:
      file=config['output_dir'] + "/table_nmc.qzv",
      rep_seqs=config['output_dir'] + "/rep_seqs_nmc.qza"
    output:
      file=config['output_dir'] + "/rep_seqs_nmc.qzv"
    shell:
      "qiime feature-table tabulate-seqs --i-data {input.rep_seqs} --o-visualization {output.file}"




rule filter_chloro_from_table:
    input:
      rep_seqs=config['output_dir'] + "/rep_seqs_nmc.qzv",
      table=config['output_dir'] + "/table_nmc.qza",
      taxonomy=config['output_dir'] + "/taxonomy.qza"
    output:
      table=config['output_dir'] + "/table_ncp_nmc.qza"
    shell:
      "qiime taxa filter-table --i-table {input.table} --i-taxonomy {input.taxonomy} --p-exclude chloroplast --o-filtered-table {output.table}"

rule filter_chloro_from_rep_seqs:
    input:
      rep_seqs=config['output_dir'] + "/rep_seqs_nmc.qza",
      table=config['output_dir'] + "/table_ncp_nmc.qza"
    output:
      rep_seqs=config['output_dir'] + "/rep_seqs_ncp_nmc.qza"
    shell:
      "qiime feature-table filter-seqs --i-data {input.rep_seqs} --i-table {input.table} --o-filtered-data {output.rep_seqs}"

rule visualize_nmc_ncp_table:
    input:
      rep_seqs=config['output_dir'] + "/rep_seqs_ncp_nmc.qza",
      table=config['output_dir'] + "/table_ncp_nmc.qza"
    output:
      file=config['output_dir'] + "/table_ncp_nmc.qzv"
    shell:
      "qiime feature-table summarize --i-table {input.table} --o-visualization {output.file}"

rule visualize_nmc_ncp_rep_seqs:
    input:
      file=config['output_dir'] + "/table_ncp_nmc.qzv",
      rep_seqs=config['output_dir'] + "/rep_seqs_ncp_nmc.qza"
    output:
      file=config['output_dir'] + "/rep_seqs_ncp_nmc.qzv"
    shell:
      "qiime feature-table tabulate-seqs --i-data {input.rep_seqs} --o-visualization {output.file}"




rule filter_Unassigned_from_table:
    input:
      file=config['output_dir'] + "/rep_seqs_ncp_nmc.qzv",
      table=config['output_dir'] + "/table_ncp_nmc.qza",
      taxonomy=config['output_dir'] + "/taxonomy.qza"
    output:
      table=config['output_dir'] + "/table_final.qza"
    shell:
      "qiime taxa filter-table --i-table {input.table} --i-taxonomy {input.taxonomy} --p-exclude Unassigned --o-filtered-table {output.table}"

rule filter_Unassigned_from_rep_seqs:
    input:
      rep_seqs=config['output_dir'] + "/rep_seqs_ncp_nmc.qza",
      table=config['output_dir'] + "/table_final.qza"
    output:
      rep_seqs=config['output_dir'] + "/rep_seqs_final.qza"
    shell:
      "qiime feature-table filter-seqs --i-data {input.rep_seqs} --i-table {input.table} --o-filtered-data {output.rep_seqs}"

rule visualize_final_table:
    input:
      rep_seqs=config['output_dir'] + "/rep_seqs_final.qza",
      table=config['output_dir'] + "/table_final.qza"
    output:
      file=config['output_dir'] + "/table_final.qzv"
    shell:
      "qiime feature-table summarize --i-table {input.table} --o-visualization {output.file}"

rule visualize_final_rep_seqs:
    input:
      table=config['output_dir'] + "/table_final.qzv",
      rep_seqs=config['output_dir'] + "/rep_seqs_final.qza"
    output:
      file=config['output_dir'] + "/rep_seqs_final.qzv"
    shell:
      "qiime feature-table tabulate-seqs --i-data {input.rep_seqs} --o-visualization {output.file}"

rule visualize_taxa_plots:
    input:
      file=config['output_dir'] + "/rep_seqs_final.qzv",
      table=config['output_dir'] + "/table_final.qza",
      taxonomy=config['output_dir'] + "/taxonomy.qza",
      map_file=config["map_file"]
    output:
      file=config['output_dir'] + "/taxa-bar-plots.qzv"
    shell:
      "qiime taxa barplot --i-table {input.table} --i-taxonomy {input.taxonomy} --o-visualization {output.file} --m-metadata-file {input.map_file}"

rule align_mafft_rep_seqs:
    input:
      file=config['output_dir'] + "/taxa-bar-plots.qzv",
      rep_seqs=config['output_dir'] + "/rep_seqs_final.qza"
    output:
      aligned_seqs=config['output_dir'] + "/aligned-rep-seqs.qza"
    threads: cpu_count
    shell:
      "qiime alignment mafft --i-sequences {input.rep_seqs} --o-alignment {output.aligned_seqs} --verbose --p-parttree --p-n-threads 4"

rule align_mask_rep_seqs:
    input:
      aligned_seqs=config['output_dir'] + "/aligned-rep-seqs.qza"
    output:
      filtered_aligned_seqs=config['output_dir'] + "/masked-aligned-rep-seqs.qza"
    shell:
      "qiime alignment mask --i-alignment {input.aligned_seqs} --o-masked-alignment {output.filtered_aligned_seqs} --verbose"

rule create_phylogenic_tree:
    input:
      filtered_aligned_seqs=config['output_dir'] + "/masked-aligned-rep-seqs.qza"
    output:
      tree=config['output_dir'] + "/unrooted-tree.qza"
    threads: cpu_count
    shell:
      "qiime phylogeny fasttree --i-alignment {input.filtered_aligned_seqs} --o-tree {output.tree} --verbose --p-n-threads {threads}"

rule create_midpoint_phylogenic_tree:
    input:
      tree=config['output_dir'] + "/unrooted-tree.qza"
    output:
      tree=config['output_dir'] + "/rooted-tree.qza"
    shell:
      "qiime phylogeny midpoint-root --i-tree {input.tree} --o-rooted-tree {output.tree}"

# Have to look into the table-nc-ncp-nmc-dn-97-0.005.qzv file for the sample-frequency-detail.csv and find the smallest number in the second column in order to define --p-sampling-depth

rule core_metrics:
    input:
      tree=config['output_dir'] + "/rooted-tree.qza",
      table=config['output_dir'] + "/table_final.qza",
      map_file=config["map_file"]
    output:
      rarefied_table=config['output_dir'] + "/rarefied_table.qza",
      faith_pd_vector=config['output_dir'] + "/faith_pd_vector.qza",
      observed_otus_vector=config['output_dir'] + "/observed_otus_vector.qza",
      shannon_vector=config['output_dir'] + "/shannon_vector.qza",
      evenness_vector=config['output_dir'] + "/evenness_vector.qza",
      unweighted_unifrac_distance_matrix=config['output_dir'] + "/unweighted_unifrac_distance_matrix.qza",
      weighted_unifrac_distance_matrix=config['output_dir'] + "/weighted_unifrac_distance_matrix.qza",
      jaccard_distance_matrix=config['output_dir'] + "/jaccard_distance_matrix.qza",
      bray_curtis_distance_matrix=config['output_dir'] + "/bray_curtis_distance_matrix.qza",
      unweighted_unifrac_pcoa_results=config['output_dir'] + "/unweighted_unifrac_pcoa_results.qza",
      weighted_unifrac_pcoa_results=config['output_dir'] + "/weighted_unifrac_pcoa_results.qza",
      jaccard_pcoa_results=config['output_dir'] + "/jaccard_pcoa_results.qza",
      bray_curtis_pcoa_results=config['output_dir'] + "/bray_curtis_pcoa_results.qza",
      unweighted_unifrac_emperor=config['output_dir'] + "/unweighted_unifrac_emperor.qzv",
      weighted_unifrac_emperor=config['output_dir'] + "/weighted_unifrac_emperor.qzv",
      jaccard_emperor=config['output_dir'] + "/jaccard_emperor.qzv",
      bray_curtis_emperor=config['output_dir'] + "/bray_curtis_emperor.qzv"
    params:
      dir=config['output_dir'] + ""
    threads: cpu_count
    run:
      my_table = qiime2.Artifact.load(input.table)
      my_tree = qiime2.Artifact.load(input.tree)

      df = my_table.view(pd.DataFrame)
      threshold = int(min(df.sum(axis=1)))

      cmd = """qiime diversity core-metrics-phylogenetic --i-phylogeny {} --i-table {}  \
      --p-sampling-depth {} --m-metadata-file {}  \
      --o-rarefied-table {} \
      --o-faith-pd-vector {} \
      --o-observed-features-vector {} \
      --o-shannon-vector {} \
      --o-evenness-vector {} \
      --o-unweighted-unifrac-distance-matrix {} \
      --o-weighted-unifrac-distance-matrix {} \
      --o-jaccard-distance-matrix {} \
      --o-bray-curtis-distance-matrix {} \
      --o-unweighted-unifrac-pcoa-results {} \
      --o-weighted-unifrac-pcoa-results {} \
      --o-jaccard-pcoa-results {} \
      --o-bray-curtis-pcoa-results {} \
      --o-unweighted-unifrac-emperor {} \
      --o-weighted-unifrac-emperor {} \
      --o-jaccard-emperor {} \
      --o-bray-curtis-emperor {} \
      --p-n-jobs-or-threads {} --verbose""".format(input.tree, input.table, threshold, input.map_file,
      output.rarefied_table,  output.faith_pd_vector,  output.observed_otus_vector, output.shannon_vector,
      output.evenness_vector, output.unweighted_unifrac_distance_matrix, output.weighted_unifrac_distance_matrix,
      output.jaccard_distance_matrix, output.bray_curtis_distance_matrix, output.unweighted_unifrac_pcoa_results,
      output.weighted_unifrac_pcoa_results, output.jaccard_pcoa_results, output.bray_curtis_pcoa_results,
      output.unweighted_unifrac_emperor, output.weighted_unifrac_emperor, output.jaccard_emperor,
      output.bray_curtis_emperor,  threads)

      try:
          cmd_output = subprocess.check_output(cmd, shell=True)
      except subprocess.CalledProcessError as c:
          sys.exit("ERROR: command {} could not be run in the shell.\n Command failed with return code: {}.".format(c.cmd, c.returncode))

# Have to look into the table-nc-ncp-nmc-dn-97-0.005.qzv file for the sample-frequency-detail.csv and find the largest number in the second column in order to define --p-max-depth

rule alpha_rarefaction:
    input:
      tree=config['output_dir'] + "/rooted-tree.qza",
      table=config['output_dir'] + "/table_final.qza",
      map_file=config["map_file"],
      rarefied_table=config['output_dir'] + "/rarefied_table.qza"
    output:
      file=config['output_dir'] + "/alpha-rarefaction.qzv"
    run:
      my_table = qiime2.Artifact.load(input.table)
      df = my_table.view(pd.DataFrame)
      threshold = int(max(df.sum(axis=1)))
      cmd = """qiime diversity alpha-rarefaction --i-phylogeny {} \
      --i-table {} --p-max-depth {} \
      --m-metadata-file {} --o-visualization {}""".format(input.tree, input.table, threshold, input.map_file, output.file)

      try:
          cmd_output = subprocess.check_output(cmd, shell=True)
      except subprocess.CalledProcessError as c:
          sys.exit("ERROR: command {} could not be run in the shell.\n Command failed with return code: {}.".format(c.cmd, c.returncode))

#Prepare QIIME2 info for phyloseq

rule export_taxonomy:
    input:
      file=config['output_dir'] + "/alpha-rarefaction.qzv",
      taxonomy=config['output_dir'] + "/taxonomy.qza"
    output:
      taxonomy=config['output_dir'] + "/taxonomy/taxonomy.tsv"
    params:
      dir=config['output_dir'] + "/taxonomy"
    shell:
      "qiime tools export --input-path {input.taxonomy} --output-path {params.dir}"

rule export_table:
    input:
      table=config['output_dir'] + "/table_final.qza"
    output:
      table=config['output_dir'] + "/final_table/feature-table.biom"
    params:
      dir=config['output_dir'] + "/final_table"
    shell:
      "qiime tools export --input-path {input.table} --output-path {params.dir}"


rule convert_table_to_txt:
    input:
      table=config['output_dir'] + "/final_table/feature-table.biom"
    output:
      table=config['output_dir'] + "/feature-table.tsv"
    shell:
      "biom convert -i {input.table} -o {output.table} --to-tsv"

rule export_tree:
    input:
      table=config['output_dir'] + "/feature-table.tsv",
      tree=config['output_dir'] + "/rooted-tree.qza"
    output:
      tree=config['output_dir'] + "/tree/tree.nwk"
    params:
      dir=config['output_dir'] + "/tree"
    shell:
      "qiime tools export --input-path {input.tree} --output-path {params.dir}"

rule change_taxonomy_label_1:
    input:
      tree=config['output_dir'] + "/tree/tree.nwk",
      taxonomy=config['output_dir'] + "/taxonomy/taxonomy.tsv"
    output:
      taxonomy=temp(config['output_dir'] + "/taxonomy/taxonomy1.tsv")
    shell:
      "sed '1s/Feature ID/#OTUID/' {input.taxonomy} > {output.taxonomy}"

rule change_taxonomy_label_2:
    input:
      taxonomy=config['output_dir'] + "/taxonomy/taxonomy1.tsv"
    output:
      taxonomy=config['output_dir'] + "/taxonomy/taxonomy2.tsv"
    shell:
      "sed '1s/Taxon/taxonomy/' {input.taxonomy} > {output.taxonomy}"

rule combine_taxonomy_and_table:
    input:
      taxonomy=config['output_dir'] + "/taxonomy/taxonomy2.tsv",
      table=config['output_dir'] + "/final_table/feature-table.biom"
    output:
      biom=config['output_dir'] + "/table-with-taxonomy.biom"
    shell:
      "biom add-metadata -i {input.table} -o {output.biom} --observation-metadata-fp {input.taxonomy} --sc-separated taxonomy --observation-header OTUID,taxonomy,Confidence"


rule convert_taxatable_to_phyloseq:
    input:
      biom=config['output_dir'] + "/table-with-taxonomy.biom"
    output:
      biom=config['output_dir'] + "/phyloseq.biom"
    shell:
      "biom convert -i {input.biom} -o {output.biom} --to-json --table-type='OTU table'"

#Create phylum level tables

rule level_2_absolute_abundance:
    input:
      biom=config['output_dir'] + "/phyloseq.biom",
      table=config['output_dir'] + "/table_final.qza",
      taxonomy=config['output_dir'] + "/taxonomy.qza"
    output:
      table=config['output_dir'] + "/phylum/table_level_2_aa.qza"
    shell:
      "qiime taxa collapse --i-table {input.table} --i-taxonomy {input.taxonomy} --o-collapsed-table {output.table} --p-level 2"

rule level_2_summarize:
    input:
      table=config['output_dir'] + "/phylum/table_level_2_aa.qza"
    output:
      file=config['output_dir'] + "/phylum/table_level_2.qzv"
    shell:
      "qiime feature-table summarize --i-table {input.table} --o-visualization {output.file}"

rule level_2_relative_frequency:
    input:
      table=config['output_dir'] + "/phylum/table_level_2_aa.qza",
      file=config['output_dir'] + "/phylum/table_level_2.qzv"
    output:
      table=config['output_dir'] + "/phylum/table_level_2_rf.qza"
    shell:
      "qiime feature-table relative-frequency --i-table {input.table} --o-relative-frequency-table {output.table}"

rule level_2_rf_export:
     input:
      table=config['output_dir'] + "/phylum/table_level_2_rf.qza"
     output:
      directory(config['output_dir'] + "/phylum/table_level_2_rf")
     shell:
      "qiime tools export --input-path {input.table} --output-path {output}"

rule level_2_aa_export:
     input:
      (config['output_dir'] + "/phylum/table_level_2_rf"),
      table=config['output_dir'] + "/phylum/table_level_2_aa.qza"
     output:
      directory(config['output_dir'] + "/phylum/table_level_2_aa")
     shell:
      "qiime tools export --input-path {input.table} --output-path {output}"

#Create class level tables

rule level_3_absolute_abundance:
    input:
      (config['output_dir'] + "/phylum/table_level_2_aa"),
      table=config['output_dir'] + "/table_final.qza",
      taxonomy=config['output_dir'] + "/taxonomy.qza"
    output:
      table=config['output_dir'] + "/class/table_level_3_aa.qza"
    shell:
      "qiime taxa collapse --i-table {input.table} --i-taxonomy {input.taxonomy} --o-collapsed-table {output.table} --p-level 3"

rule level_3_summarize:
    input:
      table=config['output_dir'] + "/class/table_level_3_aa.qza"
    output:
      file=config['output_dir'] + "/class/table_level_3.qzv"
    shell:
      "qiime feature-table summarize --i-table {input.table} --o-visualization {output.file}"

rule level_3_relative_frequency:
    input:
      table=config['output_dir'] + "/class/table_level_3_aa.qza",
      file=config['output_dir'] + "/class/table_level_3.qzv"
    output:
      table=config['output_dir'] + "/class/table_level_3_rf.qza"
    shell:
      "qiime feature-table relative-frequency --i-table {input.table} --o-relative-frequency-table {output.table}"

rule level_3_rf_export:
     input:
      table=config['output_dir'] + "/class/table_level_3_rf.qza"
     output:
      directory(config['output_dir'] + "/class/table_level_3_rf")
     shell:
      "qiime tools export --input-path {input.table} --output-path {output}"

rule level_3_aa_export:
     input:
      (config['output_dir'] + "/class/table_level_3_rf"),
      table=config['output_dir'] + "/class/table_level_3_aa.qza"
     output:
      directory(config['output_dir'] + "/class/table_level_3_aa")
     shell:
      "qiime tools export --input-path {input.table} --output-path {output}"


#Create order level tables

rule level_4_absolute_abundance:
    input:
      (config['output_dir'] + "/class/table_level_3_aa"),
      table=config['output_dir'] + "/table_final.qza",
      taxonomy=config['output_dir'] + "/taxonomy.qza"
    output:
      table=config['output_dir'] + "/order/table_level_4_aa.qza"
    shell:
      "qiime taxa collapse --i-table {input.table} --i-taxonomy {input.taxonomy} --o-collapsed-table {output.table} --p-level 4"

rule level_4_summarize:
    input:
      table=config['output_dir'] + "/order/table_level_4_aa.qza"
    output:
      file=config['output_dir'] + "/order/table_level_4.qzv"
    shell:
      "qiime feature-table summarize --i-table {input.table} --o-visualization {output.file}"

rule level_4_relative_frequency:
    input:
      table=config['output_dir'] + "/order/table_level_4_aa.qza",
      file=config['output_dir'] + "/order/table_level_4.qzv"
    output:
      table=config['output_dir'] + "/order/table_level_4_rf.qza"
    shell:
      "qiime feature-table relative-frequency --i-table {input.table} --o-relative-frequency-table {output.table}"

rule level_4_rf_export:
     input:
      table=config['output_dir'] + "/order/table_level_4_rf.qza"
     output:
      directory(config['output_dir'] + "/order/table_level_4_rf")
     shell:
      "qiime tools export --input-path {input.table} --output-path {output}"

rule level_4_aa_export:
     input:
      (config['output_dir'] + "/order/table_level_4_rf"),
      table=config['output_dir'] + "/order/table_level_4_aa.qza"
     output:
      directory(config['output_dir'] + "/order/table_level_4_aa")
     shell:
      "qiime tools export --input-path {input.table} --output-path {output}"

#Create family level tables

rule level_5_absolute_abundance:
    input:
      (config['output_dir'] + "/order/table_level_4_aa"),
      table=config['output_dir'] + "/table_final.qza",
      taxonomy=config['output_dir'] + "/taxonomy.qza"
    output:
      table=config['output_dir'] + "/family/table_level_5_aa.qza"
    shell:
      "qiime taxa collapse --i-table {input.table} --i-taxonomy {input.taxonomy} --o-collapsed-table {output.table} --p-level 5"

rule level_5_summarize:
    input:
      table=config['output_dir'] + "/family/table_level_5_aa.qza"
    output:
      file=config['output_dir'] + "/family/table_level_5.qzv"
    shell:
      "qiime feature-table summarize --i-table {input.table} --o-visualization {output.file}"

rule level_5_relative_frequency:
    input:
      table=config['output_dir'] + "/family/table_level_5_aa.qza",
      file=config['output_dir'] + "/family/table_level_5.qzv"
    output:
      table=config['output_dir'] + "/family/table_level_5_rf.qza"
    shell:
      "qiime feature-table relative-frequency --i-table {input.table} --o-relative-frequency-table {output.table}"

rule level_5_rf_export:
     input:
      table=config['output_dir'] + "/family/table_level_5_rf.qza"
     output:
      directory(config['output_dir'] + "/family/table_level_5_rf")
     shell:
      "qiime tools export --input-path {input.table} --output-path {output}"

rule level_5_aa_export:
     input:
      (config['output_dir'] + "/family/table_level_5_rf"),
      table=config['output_dir'] + "/family/table_level_5_aa.qza"
     output:
      directory(config['output_dir'] + "/family/table_level_5_aa")
     shell:
      "qiime tools export --input-path {input.table} --output-path {output}"

#Create genus level tables

rule level_6_absolute_abundance:
    input:
      (config['output_dir'] + "/family/table_level_5_aa"),
      table=config['output_dir'] + "/table_final.qza",
      taxonomy=config['output_dir'] + "/taxonomy.qza"
    output:
      table=config['output_dir'] + "/genus/table_level_6_aa.qza"
    shell:
      "qiime taxa collapse --i-table {input.table} --i-taxonomy {input.taxonomy} --o-collapsed-table {output.table} --p-level 6"

rule level_6_summarize:
    input:
      table=config['output_dir'] + "/genus/table_level_6_aa.qza"
    output:
      file=config['output_dir'] + "/genus/table_level_6.qzv"
    shell:
      "qiime feature-table summarize --i-table {input.table} --o-visualization {output.file}"

rule level_6_relative_frequency:
    input:
      table=config['output_dir'] + "/genus/table_level_6_aa.qza",
      file=config['output_dir'] + "/genus/table_level_6.qzv"
    output:
      table=config['output_dir'] + "/genus/table_level_6_rf.qza"
    shell:
      "qiime feature-table relative-frequency --i-table {input.table} --o-relative-frequency-table {output.table}"

rule level_6_rf_export:
     input:
      table=config['output_dir'] + "/genus/table_level_6_rf.qza"
     output:
      directory(config['output_dir'] + "/genus/table_level_6_rf")
     shell:
      "qiime tools export --input-path {input.table} --output-path {output}"

rule level_6_aa_export:
     input:
      (config['output_dir'] + "/genus/table_level_6_rf"),
      table=config['output_dir'] + "/genus/table_level_6_aa.qza"
     output:
      directory(config['output_dir'] + "/genus/table_level_6_aa")
     shell:
      "qiime tools export --input-path {input.table} --output-path {output}"

#Create species level tables

rule level_7_absolute_abundance:
    input:
      (config['output_dir'] + "/genus/table_level_6_aa"),
      table=config['output_dir'] + "/table_final.qza",
      taxonomy=config['output_dir'] + "/taxonomy.qza"
    output:
      table=config['output_dir'] + "/species/table_level_7_aa.qza"
    shell:
      "qiime taxa collapse --i-table {input.table} --i-taxonomy {input.taxonomy} --o-collapsed-table {output.table} --p-level 7"

rule level_7_summarize:
    input:
      table=config['output_dir'] + "/species/table_level_7_aa.qza"
    output:
      file=config['output_dir'] + "/species/table_level_7.qzv"
    shell:
      "qiime feature-table summarize --i-table {input.table} --o-visualization {output.file}"

rule level_7_relative_frequency:
    input:
      table=config['output_dir'] + "/species/table_level_7_aa.qza",
      file=config['output_dir'] + "/species/table_level_7.qzv"
    output:
      table=config['output_dir'] + "/species/table_level_7_rf.qza"
    shell:
      "qiime feature-table relative-frequency --i-table {input.table} --o-relative-frequency-table {output.table}"

rule level_7_rf_export:
     input:
      table=config['output_dir'] + "/species/table_level_7_rf.qza"
     output:
      directory(config['output_dir'] + "/species/table_level_7_rf")
     shell:
      "qiime tools export --input-path {input.table} --output-path {output}"

rule level_7_aa_export:
     input:
      (config['output_dir'] + "/species/table_level_7_rf"),
      table=config['output_dir'] + "/species/table_level_7_aa.qza"
     output:
      directory(config['output_dir'] + "/species/table_level_7_aa")
     shell:
      "qiime tools export --input-path {input.table} --output-path {output}"



rule level_2_convert_rf_table_to_txt:
    input:
      (config['output_dir'] + "/species/table_level_7_aa"),
      table=config['output_dir'] + "/phylum/table_level_2_rf/feature-table.biom"
    output:
      tsv=config['output_dir'] + "/phylum/level_2_rf.tsv"
    shell:
      "biom convert -i {input.table} -o {output.tsv} --to-tsv"

rule level_2_convert_aa_table_to_txt:
    input:
      tsv=config['output_dir'] + "/phylum/level_2_rf.tsv",
      table=config['output_dir'] + "/phylum/table_level_2_aa/feature-table.biom"
    output:
      tsv=config['output_dir'] + "/phylum/level_2_aa.tsv"
    shell:
      "biom convert -i {input.table} -o {output.tsv} --to-tsv"

rule level_3_convert_rf_table_to_txt:
    input:
      tsv=config['output_dir'] + "/phylum/level_2_aa.tsv",
      table=config['output_dir'] + "/class/table_level_3_rf/feature-table.biom"
    output:
      tsv=config['output_dir'] + "/class/level_3_rf.tsv"
    shell:
      "biom convert -i {input.table} -o {output.tsv} --to-tsv"

rule level_3_convert_aa_table_to_txt:
    input:
      tsv=config['output_dir'] + "/class/level_3_rf.tsv",
      table=config['output_dir'] + "/class/table_level_3_aa/feature-table.biom"
    output:
      tsv=config['output_dir'] + "/class/level_3_aa.tsv"
    shell:
      "biom convert -i {input.table} -o {output.tsv} --to-tsv"

rule level_4_convert_rf_table_to_txt:
    input:
      tsv=config['output_dir'] + "/class/level_3_aa.tsv",
      table=config['output_dir'] + "/order/table_level_4_rf/feature-table.biom"
    output:
      tsv=config['output_dir'] + "/order/level_4_rf.tsv"
    shell:
      "biom convert -i {input.table} -o {output.tsv} --to-tsv"

rule level_4_convert_aa_table_to_txt:
    input:
      tsv=config['output_dir'] + "/order/level_4_rf.tsv",
      table=config['output_dir'] + "/order/table_level_4_aa/feature-table.biom"
    output:
      tsv=config['output_dir'] + "/order/level_4_aa.tsv"
    shell:
      "biom convert -i {input.table} -o {output.tsv} --to-tsv"

rule level_5_convert_rf_table_to_txt:
    input:
      tsv=config['output_dir'] + "/order/level_4_aa.tsv",
      table=config['output_dir'] + "/family/table_level_5_rf/feature-table.biom"
    output:
      tsv=config['output_dir'] + "/family/level_5_rf.tsv"
    shell:
      "biom convert -i {input.table} -o {output.tsv} --to-tsv"

rule level_5_convert_aa_table_to_txt:
    input:
      tsv=config['output_dir'] + "/family/level_5_rf.tsv",
      table=config['output_dir'] + "/family/table_level_5_aa/feature-table.biom"
    output:
      tsv=config['output_dir'] + "/family/level_5_aa.tsv"
    shell:
      "biom convert -i {input.table} -o {output.tsv} --to-tsv"

rule level_6_convert_rf_table_to_txt:
    input:
      tsv=config['output_dir'] + "/family/level_5_aa.tsv",
      table=config['output_dir'] + "/genus/table_level_6_rf/feature-table.biom"
    output:
      tsv=config['output_dir'] + "/genus/level_6_rf.tsv"
    shell:
      "biom convert -i {input.table} -o {output.tsv} --to-tsv"

rule level_6_convert_aa_table_to_txt:
    input:
      tsv=config['output_dir'] + "/genus/level_6_rf.tsv",
      table=config['output_dir'] + "/genus/table_level_6_aa/feature-table.biom"
    output:
      tsv=config['output_dir'] + "/genus/level_6_aa.tsv"
    shell:
      "biom convert -i {input.table} -o {output.tsv} --to-tsv"

rule level_7_convert_rf_table_to_txt:
    input:
      tsv=config['output_dir'] + "/genus/level_6_aa.tsv",
      table=config['output_dir'] + "/species/table_level_7_rf/feature-table.biom"
    output:
      tsv=config['output_dir'] + "/species/level_7_rf.tsv"
    shell:
      "biom convert -i {input.table} -o {output.tsv} --to-tsv"

rule level_7_convert_aa_table_to_txt:
    input:
      tsv=config['output_dir'] + "/species/level_7_rf.tsv",
      table=config['output_dir'] + "/species/table_level_7_aa/feature-table.biom"
    output:
      tsv=config['output_dir'] + "/species/level_7_aa.tsv"
    shell:
      "biom convert -i {input.table} -o {output.tsv} --to-tsv"
